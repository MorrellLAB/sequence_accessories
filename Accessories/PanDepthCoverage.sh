#!/bin/bash -l
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32g
#SBATCH --tmp=32g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -u
set -o pipefail

# --- Modules ---
module load parallel
# PanDepth is a standalone binary with no standard HPC module.
# Download from https://github.com/HuiyangYu/PanDepth/releases and ensure
# the pandepth binary is on your PATH before running this script.

#   Set some defaults
declare -a pandepth_deps=(pandepth parallel)
OUTPUT_DEFAULT="$(pwd -P)/Coverage"
PROJECT_DEFAULT='PanDepthCoverage'
THREADS_DEFAULT=3
FEATURE_DEFAULT='CDS'

#   A usage message
function pandepthUsage() {
    echo -e "\
PanDepthCoverage: Calculate coverage of BAM files using PanDepth.\n\
\n\
Arguments:  --sample-list=<sample_list> [--gff=gff_file] [--bed=bed_file] \n\
            [--window-size=window_size] [--feature=feature_type] \n\
            [--min-mapq=min_mapq] [--threads=threads] \n\
            [--project=project] [--outdirectory=outdirectory] \n\
Where:      <sample_list>       is a list of BAM/CRAM files \n\
            [gff_file]          is an optional GFF/GTF file for gene-level coverage \n\
            [bed_file]          is an optional BED file for region-level coverage \n\
            [window_size]       is an optional window size in bp for windowed coverage \n\
            [feature_type]      is the GFF/GTF feature type to parse: CDS or exon [default: ${FEATURE_DEFAULT}] \n\
            [min_mapq]          is the minimum mapping quality [default: 0] \n\
            [threads]           is the number of pandepth threads per sample [default: ${THREADS_DEFAULT}] \n\
            [project]           is an optional prefix for output files [default: '${PROJECT_DEFAULT}'] \n\
            [outdirectory]      is an optional output directory [default: ${OUTPUT_DEFAULT}] \n\
\n\
Note: --gff, --bed, and --window-size are mutually exclusive. If none is given,\n\
      PanDepth reports per-chromosome coverage statistics.\n\
" >&2
    exit 1
}

#   Export the function
export -f pandepthUsage

#   Validate that only one target mode is set
function validateTargetMode() {
    local gff="$1"
    local bed="$2"
    local window="$3"
    local count=0
    [[ -n "${gff}" ]] && (( count++ ))
    [[ -n "${bed}" ]] && (( count++ ))
    [[ -n "${window}" ]] && (( count++ ))
    if [[ "${count}" -gt 1 ]]; then
        echo "Error: --gff, --bed, and --window-size are mutually exclusive. Please specify only one." >&2
        exit 1
    fi
}

#   Export the function
export -f validateTargetMode

#   A function to run pandepth on a single BAM file
function runPanDepth() {
    local bamfile="$1"
    local outdirectory="$2"
    local project="$3"
    local gff="$4"
    local bed="$5"
    local window="$6"
    local feature="$7"
    local min_mapq="$8"
    local threads="$9"

    local sample_name
    sample_name="$(basename "${bamfile}" | sed 's/\.\(bam\|cram\|sam\)$//')"
    local out_prefix="${outdirectory}/${project}_${sample_name}"

    local args=(-i "${bamfile}" -o "${out_prefix}" -t "${threads}")
    [[ -n "${gff}" ]]    && args+=(-g "${gff}" -f "${feature}")
    [[ -n "${bed}" ]]    && args+=(-b "${bed}")
    [[ -n "${window}" ]] && args+=(-w "${window}")
    [[ -n "${min_mapq}" && "${min_mapq}" != "0" ]] && args+=(-q "${min_mapq}")

    pandepth "${args[@]}"
}

#   Export the function
export -f runPanDepth

#   Main wrapper
function PanDepthCoverage() {
    local bamlist="$1"
    local outdirectory="$2"
    local project="$3"
    local gff="$4"
    local bed="$5"
    local window="$6"
    local feature="$7"
    local min_mapq="$8"
    local threads="$9"

    mkdir -p "${outdirectory}"

    #   Use parallel to process each BAM; one pandepth thread set per sample
    parallel --verbose \
        runPanDepth {} "${outdirectory}" "${project}" "${gff}" "${bed}" "${window}" "${feature}" "${min_mapq}" "${threads}" \
        :::: "${bamlist}"
}

#   Export the function
export -f PanDepthCoverage
