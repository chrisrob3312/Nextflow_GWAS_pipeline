#!/bin/bash
#SBATCH --job-name=pcair_pcrelate
#SBATCH --partition=normal
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/pcair_%j.out
#SBATCH --error=logs/pcair_%j.err

# ============================================================================
# STEP 0: GENESIS PC-AiR + PC-Relate (once per cohort, before every model)
# ============================================================================
# Produces the ancestry PCs and the GRM that ALL downstream models use:
#   results/pcair/cohort.pcair.pcs.tsv            PC-AiR PCs
#   results/pcair/cohort.pcrelate.grm.rds         GRM (2 x kinship) -> --kinship
#   results/pcair/cohort.pcrelate.kinship.rds
#   results/pcair/cohort.unrelated.txt            PC-AiR unrelated set
#   results/pcair/cohort.phenotypes.with_pcs.tsv  phenotype with PC1..PCn replaced
# submit_gwas_array.sh and submit_gxg_array.sh pick these up automatically.
#
# Usage: sbatch submit_pcair.sh
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/bin/slurm_utils.sh"
print_slurm_info
set_thread_env

DATA_DIR="${SCRIPT_DIR}/data"
OUT_DIR="${SCRIPT_DIR}/results/pcair"
mkdir -p "${OUT_DIR}"

Rscript "${SCRIPT_DIR}/bin/genesis_pcair_pcrelate.R" \
    --geno "${GENO:-${DATA_DIR}/genotypes/cohort}" \
    --phenotype "${PHENO_IN:-${DATA_DIR}/phenotypes/phenotypes.tsv}" \
    --n_pcs "${N_PCS:-10}" \
    --n_pcs_pcrelate "${N_PCS_PCRELATE:-5}" \
    --kin_thresh "${KIN_THRESH:-0.0442}" \
    --ld_r2 "${LD_R2:-0.1}" \
    --iterations "${PCAIR_ITER:-2}" \
    --output_prefix "${OUT_DIR}/cohort" \
    --threads "${SLURM_CPUS_PER_TASK}"

echo "PC-AiR / PC-Relate complete: ${OUT_DIR}"
