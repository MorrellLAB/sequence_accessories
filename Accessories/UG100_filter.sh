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

# --- File Paths ---
INPUT_FILE="/scratch.global/pmorrell/Cowpea/Joint_Genotyping/Cowpea_UG100_Cohort.bcf"
OUTPUT_DIR="/scratch.global/pmorrell/Cowpea/Joint_Genotyping"
OUT_BCF="${OUTPUT_DIR}/Cowpea_UG100_QualityFiltered.bcf"

echo "Starting Quality Filtering..."
echo "Input: ${INPUT_FILE}"
echo "Output: ${OUT_BCF}"

# --- Filtering Pipeline ---
# 1. Site-level QUAL filter
# 2. Genotype-level GQ filter (sets failing genotypes to missing)
# 3. Recalculate F_MISSING based on the new genotypes
# 4. Remove sites that are now 100% missing
bcftools view --threads "${SLURM_CPUS_PER_TASK}" -i 'QUAL >= 20' "${INPUT_FILE}" | \
bcftools filter --threads "${SLURM_CPUS_PER_TASK}" -e 'FMT/GQ < 20' --set-GTs . | \
bcftools +fill-tags --threads "${SLURM_CPUS_PER_TASK}" -- -t F_MISSING | \
bcftools view --threads "${SLURM_CPUS_PER_TASK}" -e 'F_MISSING == 1' -O b -o "${OUT_BCF}"

# --- Indexing ---
echo "Indexing..."
bcftools index --threads "${SLURM_CPUS_PER_TASK}" "${OUT_BCF}"

echo "Filtering complete."

