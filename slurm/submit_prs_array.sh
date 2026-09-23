#!/bin/bash
#SBATCH --job-name=prs_admixed
#SBATCH --array=1-6
#SBATCH --partition=normal
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/prs_%A_%a.out
#SBATCH --error=logs/prs_%A_%a.err

# ============================================================================
# SLURM Array Job: Multi-Method PRS Calculation
# ============================================================================
# Runs 6 PRS methods in parallel:
#   1. PRS-CSx (primary)
#   2. GAUDI
#   3. DiscoDivas
#   4. SDPR_admix
#   5. MUSSEL
#   6. PROSPER
#
# Usage:
#   sbatch submit_prs_array.sh <trait>
#
# Example:
#   sbatch submit_prs_array.sh OS
#   sbatch submit_prs_array.sh relapse
#   sbatch submit_prs_array.sh MRD
# ============================================================================

set -euo pipefail

# Load utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/bin/slurm_utils.sh"

# Print job info
print_slurm_info
set_thread_env

# Map array task ID to method
declare -a METHODS=("prs_csx" "gaudi" "disco_divas" "sdpr_admix" "mussel" "prosper")
METHOD="${METHODS[$((SLURM_ARRAY_TASK_ID - 1))]}"

TRAIT="${1:-OS}"

echo "=============================================="
echo "Trait: ${TRAIT}"
echo "Method: ${METHOD}"
echo "=============================================="

# Configuration
DATA_DIR="${SCRIPT_DIR}/data"
RESULTS_DIR="${SCRIPT_DIR}/results/prs/${TRAIT}"
mkdir -p "${RESULTS_DIR}"

# GWAS summary statistics directory (ancestry-specific for PRS-CSx)
SUMSTATS_DIR="${SCRIPT_DIR}/results/gwas/sumstats_for_prs"

# LD reference (ancestry-matched)
LD_REF="${DATA_DIR}/ld_reference"

# Genotypes
GENO="${DATA_DIR}/genotypes/cohort"

# Local ancestry (for LA-aware methods)
LOCAL_ANC="${DATA_DIR}/local_ancestry/cohort.msp.tsv.gz"

# Phenotype for validation
PHENO="${DATA_DIR}/phenotypes/phenotypes.tsv"

# Run PRS calculation
Rscript "${SCRIPT_DIR}/bin/calculate_prs_admixed.R" \
    --method "${METHOD}" \
    --sumstats_dir "${SUMSTATS_DIR}" \
    --ld_ref "${LD_REF}" \
    --geno "${GENO}" \
    --local_ancestry "${LOCAL_ANC}" \
    --ancestries "EUR,AFR,AMR" \
    --ancestry_props "EUR:0.35,AFR:0.15,AMR:0.50" \
    --phenotype "${PHENO}" \
    --trait "${TRAIT}" \
    --ancestry_col "GRAF_ANC" \
    --validate \
    --stratify_validation \
    --output_prefix "${RESULTS_DIR}/${METHOD}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --verbose

# ----------------------------------------------------------------------------
# Local-ancestry PARTIAL scores: run in the PRS-CSx task once its shrunk
# per-ancestry posterior weights exist (preferred); fall back to Tractor-GENESIS
# betas from results/gwas/<trait>/POOLED.sumstats.gz if PRS-CSx produced none.
# ----------------------------------------------------------------------------
if [[ "${METHOD}" == "prs_csx" && -d "${DATA_DIR}/tractor" ]]; then
    if ls "${RESULTS_DIR}"/prs_csx*.weights.tsv.gz >/dev/null 2>&1; then
        LA_WEIGHTS="${RESULTS_DIR}"                  # directory of <..>.<anc>.weights.tsv.gz
        echo "LA partial scores: using PRS-CSx per-ancestry weights"
    else
        LA_WEIGHTS="${SCRIPT_DIR}/results/gwas/${TRAIT}/POOLED.sumstats.gz"
        echo "LA partial scores: PRS-CSx weights not found, falling back to Tractor-GENESIS betas (no shrinkage)"
    fi
    Rscript "${SCRIPT_DIR}/bin/calculate_prs_admixed.R" \
        --method la_partial \
        --weights "${LA_WEIGHTS}" \
        --tractor_prefix "${DATA_DIR}/tractor/all_chr" \
        --ancestries "EUR,AFR,AMR" \
        --phenotype "${PHENO}" \
        --trait "${TRAIT}" \
        --output_prefix "${RESULTS_DIR}/${TRAIT}" \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --verbose
fi

echo "Completed PRS method ${METHOD} for trait ${TRAIT}"
