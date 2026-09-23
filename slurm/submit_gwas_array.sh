#!/bin/bash
#SBATCH --job-name=tractor_genesis
#SBATCH --array=1-22
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/gwas_%A_%a.out
#SBATCH --error=logs/gwas_%A_%a.err

# ============================================================================
# SLURM Array Job: Tractor-GENESIS GWAS by Chromosome
# ============================================================================
# ALL traits go through the same Tractor-GENESIS adapter so results are
# comparable:
#   MRD      -> binary      (GENESIS logistic mixed model null)
#   relapse  -> survival    (Cox null with kinship frailty, time-to-event)
#   OS       -> survival
#
# Usage:
#   sbatch submit_gwas_array.sh <stratum> <trait> <model>
#
#   sbatch submit_gwas_array.sh LAT1 OS      survival
#   sbatch submit_gwas_array.sh LAT1 relapse survival
#   sbatch submit_gwas_array.sh LAT1 MRD     binary
#
# Survival traits expect columns <trait>_time and <trait>_status in the
# phenotype file (override with TIME_COL / EVENT_COL env vars).
# Output: results/gwas/<trait>/<stratum>/chr<N>.tractor_genesis.tsv.gz
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/bin/slurm_utils.sh"
print_slurm_info
set_thread_env

STRATUM="${1:?stratum required (EUR, AAC, LAT1, LAT2, EAS, SAS, OTHER, POOLED)}"
TRAIT="${2:?trait required (OS, relapse, MRD)}"
MODEL="${3:-survival}"
CHR="${SLURM_ARRAY_TASK_ID}"

echo "=============================================="
echo "Stratum: ${STRATUM}   Trait: ${TRAIT}   Model: ${MODEL}   Chr: ${CHR}"
echo "=============================================="

DATA_DIR="${SCRIPT_DIR}/data"
RESULTS_DIR="${SCRIPT_DIR}/results/gwas/${TRAIT}/${STRATUM}"
mkdir -p "${RESULTS_DIR}"

PHENO="${DATA_DIR}/phenotypes/phenotypes.tsv"
TRACTOR_PREFIX="${DATA_DIR}/tractor/${STRATUM}/chr${CHR}"
KINSHIP="${DATA_DIR}/kinship/kinship_matrix.rds"

# Ancestral populations for the Tractor decomposition
case "${STRATUM}" in
    AAC)            ANCESTRIES="EUR,AFR" ;;
    LAT1|LAT2|POOLED|OTHER) ANCESTRIES="EUR,AFR,AMR" ;;
    *)              ANCESTRIES="${ANCESTRIES:-EUR,AFR,AMR}" ;;
esac

MODEL_ARGS=""
if [[ "${MODEL}" == "survival" ]]; then
    MODEL_ARGS="--time_col ${TIME_COL:-${TRAIT}_time} --event_col ${EVENT_COL:-${TRAIT}_status}"
fi

Rscript "${SCRIPT_DIR}/bin/tractor_genesis_adapter.R" \
    --tractor_prefix "${TRACTOR_PREFIX}" \
    --phenotype "${PHENO}" \
    --trait "${TRAIT}" \
    --model "${MODEL}" \
    ${MODEL_ARGS} \
    --ancestries "${ANCESTRIES}" \
    --ref_ancestry "EUR" \
    --kinship "${KINSHIP}" \
    --covariates "age,sex,PC1,PC2,PC3,PC4,PC5" \
    --stratum "${STRATUM}" \
    --min_stratum_n 30 \
    --ancestry_config "${SCRIPT_DIR}/bin/ancestry_config.R" \
    --mac_min 10 \
    --chromosome "${CHR}" \
    --output_prefix "${RESULTS_DIR}/chr${CHR}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --verbose

echo "Completed chr${CHR} | ${TRAIT} (${MODEL}) | ${STRATUM}"
