#!/bin/bash -l
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=48g
#SBATCH --tmp=48g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pmorrell@umn.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

set -e
set -u
set -o pipefail

# --- Modules ---
module load bcftools/1.21

# --- Arguments ---
if [[ "$#" -ne 2 ]]; then
	echo "Usage: $0 <INPUT_FILE> <OUT_DIR>"
	exit 1
fi

INPUT_FILE="$1"
OUT_DIR="$2"

INPUT_BASENAME="$(basename "${INPUT_FILE}")"
INPUT_STEM="${INPUT_BASENAME%.*}"
OUT_BCF="${OUT_DIR}/${INPUT_STEM}_QualityFiltered.bcf"
THREADS="${SLURM_CPUS_PER_TASK:-1}"

echo "Starting Quality Filtering..."
echo "Input: ${INPUT_FILE}"
echo "Output: ${OUT_BCF}"

# --- Filtering Pipeline ---
# 1. Site-level QUAL filter
# 2. Genotype-level GQ filter (sets failing genotypes to missing)
# 3. Recalculate F_MISSING and AC based on the new genotypes
# 4. Remove sites with >70% missing genotypes OR no alternate alleles (AC=0)
# 5. Recalculate all INFO tags for the final filtered output
bcftools view --threads "${THREADS}" -i 'QUAL >= 20' "${INPUT_FILE}" | \
bcftools filter --threads "${THREADS}" -e 'FMT/GQ < 20' --set-GTs . | \
bcftools +fill-tags --threads "${THREADS}" -- -t F_MISSING,AC,AN | \
bcftools view --threads "${THREADS}" -e 'F_MISSING > 0.7 || AC == 0' | \
bcftools +fill-tags --threads "${THREADS}" -O b -o "${OUT_BCF}" -- -t all

# --- Indexing ---
echo "Indexing..."
bcftools index --threads "${THREADS}" "${OUT_BCF}"

echo "Filtering complete."

