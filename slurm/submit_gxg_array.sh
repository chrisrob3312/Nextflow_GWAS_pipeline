#!/bin/bash
#SBATCH --job-name=gxg_interaction
#SBATCH --array=0-7
#SBATCH --partition=normal
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/gxg_%A_%a.out
#SBATCH --error=logs/gxg_%A_%a.err

# ============================================================================
# SLURM Array Job: GxG (Epistasis) Interaction Tests on GWAS Hits
# ============================================================================
# END OF PIPELINE. One array task per analysis stratum:
#   0 = POOLED  1 = EUR  2 = AAC  3 = LAT1  4 = LAT2  5 = EAS  6 = SAS  7 = OTHER
# Strata with N < 30 are skipped automatically by the script.
# A final non-array job computes ancestry heterogeneity of the interactions.
#
# Usage:
#   sbatch submit_gxg_array.sh <trait> [model]
#   sbatch submit_gxg_array.sh OS survival
#   sbatch submit_gxg_array.sh relapse binary
#   sbatch submit_gxg_array.sh MRD binary
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${SCRIPT_DIR}/bin/slurm_utils.sh"
print_slurm_info
set_thread_env

TRAIT="${1:-OS}"
MODEL="${2:-survival}"
STRATUM=$(get_array_item "POOLED,EUR,AAC,LAT1,LAT2,EAS,SAS,OTHER")

echo "=============================================="
echo "Trait: ${TRAIT}   Model: ${MODEL}   Stratum: ${STRATUM}"
echo "=============================================="

DATA_DIR="${SCRIPT_DIR}/data"
GWAS_DIR="${SCRIPT_DIR}/results/gwas/${TRAIT}"     # Tractor-GENESIS output for this trait (all strata + POOLED)
RESULTS_DIR="${SCRIPT_DIR}/results/gxg/${TRAIT}"
# Phenotype with PC-AiR PCs (step 0) if present, else the raw phenotype file
PCAIR_PHENO="${SCRIPT_DIR}/results/pcair/cohort.phenotypes.with_pcs.tsv"
PHENO="${PHENO:-$([[ -f "${PCAIR_PHENO}" ]] && echo "${PCAIR_PHENO}" || echo "${DATA_DIR}"/phenotypes/phenotypes.tsv)}"
echo "Phenotype: ${PHENO}"
mkdir -p "${RESULTS_DIR}"

# Hits come from EVERY GWAS variant for this trait: pooled meta + each stratum
SUMSTATS=$(ls "${GWAS_DIR}"/POOLED.sumstats.gz "${GWAS_DIR}"/{EUR,AAC,LAT1,LAT2,EAS,SAS,OTHER}.sumstats.gz 2>/dev/null | paste -sd, -)

SURV_ARGS=""
if [[ "${MODEL}" == "survival" ]]; then
    SURV_ARGS="--time_col ${TRAIT}_time --event_col ${TRAIT}_status"
fi

# Optional: Tractor dosages for local-ancestry-aware epistasis (pooled task only)
TRACTOR_ARGS=""
if [[ "${STRATUM}" == "POOLED" && -d "${DATA_DIR}/tractor" ]]; then
    TRACTOR_ARGS="--tractor_prefix ${DATA_DIR}/tractor/all_chr --ancestries EUR,AFR,AMR"
fi

# ---- customisable (env vars) ----------------------------------------------
COVARIATES="${COVARIATES:-age,sex,PC1,PC2,PC3,PC4,PC5}"   # any phenotype columns
ANCESTRY_COL="${ANCESTRY_COL:-GRAF_ANC}"
KNOWN_LOCI="${KNOWN_LOCI:-${SCRIPT_DIR}/assets/known_leukemia_risk_loci.tsv}"
LOCI_ARGS="--known_loci ${KNOWN_LOCI}"
[[ "${GXG_USE_DEFAULT_LOCI:-true}" == "false" ]] && LOCI_ARGS="${LOCI_ARGS} --no_default_loci"
[[ -n "${GXG_CUSTOM_VARIANTS:-}" ]] && LOCI_ARGS="${LOCI_ARGS} --custom_variants ${GXG_CUSTOM_VARIANTS}"   # GRCh38 list
PRUNE_MODE="${GXG_PRUNE_MODE:-variant}"     # variant = keep strongest per LD cluster; pair = keep all
# ---------------------------------------------------------------------------

Rscript "${SCRIPT_DIR}/bin/run_gxg_interaction.R" \
    --sumstats "${SUMSTATS}" \
    --p_threshold "${GXG_P_THRESHOLD:-1e-5}" \
    --max_hits "${GXG_MAX_HITS:-200}" \
    ${LOCI_ARGS} \
    --prune_mode "${PRUNE_MODE}" \
    --geno "${DATA_DIR}/genotypes/cohort" \
    --phenotype "${PHENO}" \
    --trait "${TRAIT}" \
    --model "${MODEL}" \
    ${SURV_ARGS} \
    --covariates "${COVARIATES}" \
    --ancestry_col "${ANCESTRY_COL}" \
    --ancestry_config "${SCRIPT_DIR}/bin/ancestry_config.R" \
    --stratum "${STRATUM}" \
    --min_stratum_n 30 \
    ${TRACTOR_ARGS} \
    --output_prefix "${RESULTS_DIR}/${TRAIT}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --verbose

echo "Completed GxG for ${TRAIT} / ${STRATUM}"

# Last array task submits the ancestry-heterogeneity job (all strata in one run)
if [[ "${SLURM_ARRAY_TASK_ID}" == "${SLURM_ARRAY_TASK_MAX}" ]]; then
    sbatch --dependency=afterany:${SLURM_ARRAY_JOB_ID} \
           --job-name=gxg_het_${TRAIT} --partition=normal --time=12:00:00 \
           --mem=32G --cpus-per-task=16 \
           --output="${SCRIPT_DIR}/logs/gxg_het_${TRAIT}_%j.out" \
           --wrap="source ${SCRIPT_DIR}/bin/slurm_utils.sh; set_thread_env; \
                   Rscript ${SCRIPT_DIR}/bin/run_gxg_interaction.R \
                     --sumstats '${SUMSTATS}' --p_threshold ${GXG_P_THRESHOLD:-1e-5} --max_hits ${GXG_MAX_HITS:-200} \
                     ${LOCI_ARGS} --prune_mode ${PRUNE_MODE} \
                     --geno ${DATA_DIR}/genotypes/cohort \
                     --phenotype ${PHENO} \
                     --trait ${TRAIT} --model ${MODEL} ${SURV_ARGS} \
                     --covariates ${COVARIATES} --ancestry_col ${ANCESTRY_COL} \
                     --ancestry_config ${SCRIPT_DIR}/bin/ancestry_config.R --min_stratum_n 30 \
                     --output_prefix ${RESULTS_DIR}/${TRAIT}.het --threads 16"
fi
