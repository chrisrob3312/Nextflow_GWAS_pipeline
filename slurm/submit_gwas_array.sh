#!/bin/bash
#SBATCH --job-name=tractor_gwas
#SBATCH --array=1-22
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/gwas_%A_%a.out
#SBATCH --error=logs/gwas_%A_%a.err

# ============================================================================
# SLURM Array Job: Tractor GWAS by Chromosome
# ============================================================================
# Runs local ancestry-aware GWAS across 22 chromosomes in parallel
#
# Usage:
#   sbatch submit_gwas_array.sh <stratum>
#
# Where stratum is: EUR, AAC, LAT1, LAT2, OTHER, or POOLED
#
# Example:
#   sbatch submit_gwas_array.sh LAT1
#   sbatch submit_gwas_array.sh LAT2
#   sbatch submit_gwas_array.sh EUR
# ============================================================================

set -euo pipefail

# Load utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/bin/slurm_utils.sh"

# Print job info
print_slurm_info
set_thread_env

# Get parameters
STRATUM="${1:-POOLED}"
CHR="${SLURM_ARRAY_TASK_ID}"

echo "=============================================="
echo "Stratum: ${STRATUM}"
echo "Chromosome: ${CHR}"
echo "=============================================="

# Configuration
DATA_DIR="${SCRIPT_DIR}/data"
RESULTS_DIR="${SCRIPT_DIR}/results/gwas/${STRATUM}"
mkdir -p "${RESULTS_DIR}"

# Phenotype file with ancestry column
PHENO="${DATA_DIR}/phenotypes/phenotypes.tsv"

# Tractor files (per chromosome)
TRACTOR_PREFIX="${DATA_DIR}/tractor/chr${CHR}"

# Kinship matrix (cohort-specific)
KINSHIP="${DATA_DIR}/kinship/kinship_matrix.rds"

# Run Tractor-GENESIS
Rscript "${SCRIPT_DIR}/bin/tractor_genesis_adapter.R" \
    --tractor_prefix "${TRACTOR_PREFIX}" \
    --phenotype "${PHENO}" \
    --trait "OS_status" \
    --model survival \
    --time_col "OS_time" \
    --event_col "OS_status" \
    --ancestries "EUR,AFR,AMR" \
    --ref_ancestry "EUR" \
    --kinship "${KINSHIP}" \
    --covariates "age,sex,PC1,PC2,PC3,PC4,PC5" \
    --stratum "${STRATUM}" \
    --min_stratum_n 30 \
    --chromosome "${CHR}" \
    --output_prefix "${RESULTS_DIR}/chr${CHR}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --verbose

echo "Completed chromosome ${CHR} for stratum ${STRATUM}"
