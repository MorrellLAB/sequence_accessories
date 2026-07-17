#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=64gb
#SBATCH --tmp=100gb
#SBATCH -t 12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

set -euo pipefail

module load emboss/6.6.0-rocky8

# This script identifies simple sequence repeats with EMBOSS fuzznuc and
# creates a mask file with padded intervals around each repeat.
#
# It scans for:
# 1) Mononucleotide repeats with length >= MIN_MONO_REPEATS
# 2) Dinucleotide repeats with repeat count >= MIN_DINUC_REPEATS
# 3) Trinucleotide repeats with repeat count >= MIN_TRINUC_REPEATS
#
# Output can be BED or GFF.

# Dependencies
# EMBOSS fuzznuc (6.6+)

#------------------
# User provided input arguments

# Input FASTA
FASTA="/path/to/genome.fasta"

# Output directory and file prefix
OUT_DIR="/path/to/mask_output"
OUT_PREFIX="genome_repeat_mask"

# Output format: bed or gff
OUTPUT_FORMAT="bed"

# Flank size on each side of the repeat
FLANK_BP="20"

# Repeat thresholds
MIN_MONO_REPEATS="10"
MIN_DINUC_REPEATS="10"
MIN_TRINUC_REPEATS="10"

# Merge overlapping/adjacent padded intervals: yes or no
MERGE_INTERVALS="yes"

#------------------
mkdir -p "${OUT_DIR}"

if [[ ! -f "${FASTA}" ]]; then
    echo "ERROR: FASTA not found: ${FASTA}" >&2
    exit 1
fi

if ! command -v fuzznuc >/dev/null 2>&1; then
    echo "ERROR: fuzznuc not found on PATH." >&2
    exit 1
fi

if [[ ! "${OUTPUT_FORMAT}" =~ ^(bed|gff)$ ]]; then
    echo "ERROR: OUTPUT_FORMAT must be 'bed' or 'gff'." >&2
    exit 1
fi

if [[ ! "${MERGE_INTERVALS}" =~ ^(yes|no)$ ]]; then
    echo "ERROR: MERGE_INTERVALS must be 'yes' or 'no'." >&2
    exit 1
fi

if ! [[ "${FLANK_BP}" =~ ^[0-9]+$ && "${MIN_MONO_REPEATS}" =~ ^[0-9]+$ && "${MIN_DINUC_REPEATS}" =~ ^[0-9]+$ && "${MIN_TRINUC_REPEATS}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: FLANK_BP, MIN_MONO_REPEATS, MIN_DINUC_REPEATS, and MIN_TRINUC_REPEATS must be non-negative integers." >&2
    exit 1
fi

build_repeat_pattern() {
    local unit="$1"
    local repeat_count="$2"
    local i
    local pattern=""
    for ((i=0; i<repeat_count; i++)); do
        pattern+="${unit}"
    done
    echo "${pattern}"
}

reverse_complement_motif() {
    local motif="$1"
    local reversed=""
    local base
    local i

    for ((i=${#motif}-1; i>=0; i--)); do
        base=${motif:i:1}
        case "${base}" in
            A) reversed+="T" ;;
            C) reversed+="G" ;;
            G) reversed+="C" ;;
            T) reversed+="A" ;;
            *)
                echo "ERROR: Unsupported base '${base}' in motif '${motif}'." >&2
                exit 1
                ;;
        esac
    done

    echo "${reversed}"
}

# Treat a repeat as the same class across cyclic shifts and reverse complements.
# This collapses motifs like AG/GA and AT/TA to a single scan target.
canonical_repeat_motif() {
    local motif="$1"
    local motif_len=${#motif}
    local rotated
    local reverse_complement
    local candidate
    local best=""
    local i

    for ((i=0; i<motif_len; i++)); do
        rotated="${motif:i}${motif:0:i}"
        reverse_complement=$(reverse_complement_motif "${rotated}")

        for candidate in "${rotated}" "${reverse_complement}"; do
            if [[ -z "${best}" || "${candidate}" < "${best}" ]]; then
                best="${candidate}"
            fi
        done
    done

    echo "${best}"
}

is_canonical_motif() {
    local motif="$1"
    local canonical

    canonical=$(canonical_repeat_motif "${motif}")
    [[ "${motif}" == "${canonical}" ]]
}

run_fuzznuc_for_motif() {
    local motif="$1"
    local repeat_class="$2"
    local out_gff="$3"

    # -complement Y finds matches on both strands.
    fuzznuc \
        -sequence "${FASTA}" \
        -pattern "${motif}" \
        -complement Y \
        -rformat2 gff \
        -outfile "${out_gff}"

    awk -v cls="${repeat_class}" -v motif="${motif}" 'BEGIN { FS=OFS="\t" }
        $0 !~ /^#/ && NF >= 8 {
            print $1, $4, $5, $7, cls, motif
        }
    ' "${out_gff}" >> "${TMP_HITS}"
}

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

SEQ_LENGTHS="${TMP_DIR}/seq_lengths.tsv"
TMP_HITS="${TMP_DIR}/raw_hits.tsv"
PADDED_BED6="${TMP_DIR}/padded_hits.bed6"
MERGED_BED3="${TMP_DIR}/merged_hits.bed3"

: > "${TMP_HITS}"

# Build sequence lengths from FASTA for boundary clamping.
awk 'BEGIN { name=""; len=0 }
    /^>/ {
        if (name != "") {
            print name "\t" len
        }
        name = substr($1, 2)
        len = 0
        next
    }
    {
        gsub(/[[:space:]]/, "", $0)
        len += length($0)
    }
    END {
        if (name != "") {
            print name "\t" len
        }
    }
' "${FASTA}" > "${SEQ_LENGTHS}"

# Mononucleotide repeats.
for base in A C G T; do
    motif=$(build_repeat_pattern "${base}" "${MIN_MONO_REPEATS}")
    run_fuzznuc_for_motif "${motif}" "mono" "${TMP_DIR}/mono_${base}.gff"
done

# Dinucleotide repeats (excluding AA/CC/GG/TT since they are mono-like).
for b1 in A C G T; do
    for b2 in A C G T; do
        if [[ "${b1}" == "${b2}" ]]; then
            continue
        fi
        motif=$(build_repeat_pattern "${b1}${b2}" "${MIN_DINUC_REPEATS}")
        if ! is_canonical_motif "${motif}"; then
            continue
        fi
        run_fuzznuc_for_motif "${motif}" "dinuc" "${TMP_DIR}/dinuc_${b1}${b2}.gff"
    done
done

# Trinucleotide repeats (excluding AAA/CCC/GGG/TTT since they are mono-like).
for b1 in A C G T; do
    for b2 in A C G T; do
        for b3 in A C G T; do
            if [[ "${b1}" == "${b2}" && "${b2}" == "${b3}" ]]; then
                continue
            fi
            motif=$(build_repeat_pattern "${b1}${b2}${b3}" "${MIN_TRINUC_REPEATS}")
            if ! is_canonical_motif "${motif}"; then
                continue
            fi
            run_fuzznuc_for_motif "${motif}" "trinuc" "${TMP_DIR}/trinuc_${b1}${b2}${b3}.gff"
        done
    done
done

if [[ ! -s "${TMP_HITS}" ]]; then
    echo "No repeats found at the requested thresholds." >&2
    exit 0
fi

# Clamp padded intervals to [1, seq_len] for GFF and [0, seq_len] for BED start/end.
awk 'NR==FNR { seqlen[$1]=$2; next }
    {
        chr=$1
        start=$2
        end=$3
        strand=$4
        cls=$5
        motif=$6

        if (!(chr in seqlen)) {
            next
        }

        bed_start = start - 1 - flank
        if (bed_start < 0) {
            bed_start = 0
        }

        bed_end = end + flank
        if (bed_end > seqlen[chr]) {
            bed_end = seqlen[chr]
        }

        if (bed_start < bed_end) {
            name = cls "_" motif
            print chr, bed_start, bed_end, name, ".", strand
        }
    }
' flank="${FLANK_BP}" "${SEQ_LENGTHS}" "${TMP_HITS}" \
    | sort -k1,1 -k2,2n -k3,3n > "${PADDED_BED6}"

if [[ "${MERGE_INTERVALS}" == "yes" ]]; then
    awk 'BEGIN { OFS="\t" }
        NR==1 {
            c=$1; s=$2; e=$3; next
        }
        {
            if ($1 == c && $2 <= e + 1) {
                if ($3 > e) {
                    e = $3
                }
            } else {
                print c, s, e
                c=$1; s=$2; e=$3
            }
        }
        END {
            if (NR > 0) {
                print c, s, e
            }
        }
    ' "${PADDED_BED6}" > "${MERGED_BED3}"
fi

if [[ "${OUTPUT_FORMAT}" == "bed" ]]; then
    OUT_FILE="${OUT_DIR}/${OUT_PREFIX}.bed"
    if [[ "${MERGE_INTERVALS}" == "yes" ]]; then
        cp "${MERGED_BED3}" "${OUT_FILE}"
    else
        cp "${PADDED_BED6}" "${OUT_FILE}"
    fi
else
    OUT_FILE="${OUT_DIR}/${OUT_PREFIX}.gff"
    {
        echo "##gff-version 3"
        if [[ "${MERGE_INTERVALS}" == "yes" ]]; then
            awk 'BEGIN { OFS="\t" }
                {
                    id = "mask_" NR
                    gff_start = $2 + 1
                    gff_end = $3
                    print $1, "fuzznuc_mask", "repeat_mask", gff_start, gff_end, ".", ".", ".", "ID=" id ";Note=merged_simple_repeat_mask"
                }
            ' "${MERGED_BED3}"
        else
            awk 'BEGIN { OFS="\t" }
                {
                    id = "mask_" NR
                    gff_start = $2 + 1
                    gff_end = $3
                    print $1, "fuzznuc_mask", "repeat_mask", gff_start, gff_end, ".", $6, ".", "ID=" id ";Class=" $4
                }
            ' "${PADDED_BED6}"
        fi
    } > "${OUT_FILE}"
fi

if [[ ! -s "${OUT_FILE}" ]]; then
    echo "ERROR: Output file is missing or empty: ${OUT_FILE}" >&2
    exit 1
fi

echo "Wrote mask file: ${OUT_FILE}"
