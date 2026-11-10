#!/bin/bash
#SBATCH --job-name=gxg_interaction
#SBATCH --array=0-7
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
# --partition, --account, --output and --error come from the sbatch command line
# (slurm/submit_full_pipeline.sh, or the direct-submit form in docs/HPC_SITE.md).

# ============================================================================
# SLURM Array Job: GxG (Epistasis) Interaction Tests on GWAS Hits
# ============================================================================
# END OF PIPELINE. One array task per analysis stratum: task 0 = POOLED, then one per stratum of
# STRATA_LIST (default: the bundle's ancestry/strata.tsv, sorted). The launcher passes
# --array=0-<number of strata>; by hand, pass it too (the #SBATCH line above is only a default).
# Strata with N < MIN_STRATUM_N (default 50) are skipped automatically by the script.
# A final non-array job computes ancestry heterogeneity of the interactions.
#
# Usage (normally submitted by slurm/submit_full_pipeline.sh; by hand, from the repo root,
# with the site arguments of docs/HPC_SITE.md "Submitting one script by hand"):
#   sbatch <site args> slurm/submit_gxg_array.sh <trait> [model]
#   sbatch <site args> slurm/submit_gxg_array.sh os  survival
#   sbatch <site args> slurm/submit_gxg_array.sh efs survival
#   sbatch <site args> slurm/submit_gxg_array.sh mrd binary
# ============================================================================

set -euo pipefail

# Repo root. Slurm runs a batch script from a copy in its spool directory, so this file's own
# path does not lead back to the repo inside a job. The launcher exports GWAS_REPO; a direct
# sbatch from the repo root leaves the repo in SLURM_SUBMIT_DIR; `bash slurm/...` uses the path.
for SCRIPT_DIR in "${GWAS_REPO:-}" "${SLURM_SUBMIT_DIR:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; do
    [[ -f "${SCRIPT_DIR}/bin/slurm_utils.sh" ]] && break
done
[[ -f "${SCRIPT_DIR}/bin/slurm_utils.sh" ]] || { echo "ERROR: cannot find the repo (submit from the repo root, or export GWAS_REPO)" >&2; exit 1; }
source "${SCRIPT_DIR}/bin/slurm_utils.sh"
rr_begin gxg              # run record row first, so a failure from here on is recorded too
source "${SCRIPT_DIR}/bin/activate_env.sh" gwas-r     # conda env gwas-r, TMPDIR on disk (docs/HPC_SITE.md)
print_slurm_info
set_thread_env

TRAIT="${1:?trait required (mrd, efs, os)}"
MODEL="${2:-survival}"
GWAS_STRATA="${STRATA_LIST:-$(bundle_strata)}" || exit 1
STRATUM=$(get_array_item "POOLED,${GWAS_STRATA}")

echo "=============================================="
echo "Trait: ${TRAIT}   Model: ${MODEL}   Stratum: ${STRATUM}"
echo "=============================================="

DATA_DIR="${SCRIPT_DIR}/data"
GWAS_DIR="${GWAS_RESULTS}/gwas/${TRAIT}"     # Tractor-GENESIS output for this trait (all strata + POOLED)
RESULTS_DIR="${GWAS_RESULTS}/gxg/${TRAIT}"
# Phenotype with PC-AiR PCs (step 0) if present, else the raw phenotype file
PCAIR_PHENO="${SCRIPT_DIR}/results/pcair/cohort.phenotypes.with_pcs.tsv"
PHENO="${PHENO:-$([[ -f "${PCAIR_PHENO}" ]] && echo "${PCAIR_PHENO}" || echo "${PHENO_IN:-${HANDOFF_DIR:-${SCRIPT_DIR}/data/handoff}/phenotypes/phenotypes.tsv}")}"
echo "Phenotype: ${PHENO}"
mkdir -p "${RESULTS_DIR}"

# Hits come from EVERY GWAS variant for this trait: pooled meta + each stratum. A stratum with no
# sumstats (skipped below MIN_STRATUM_N, or not analysed) is left out and named; none at all stops
# the task. (The old `ls ... | paste` died under pipefail as soon as one file was missing.)
SUMSTATS_FILES=(); SUMSTATS_MISSING=()
for s in POOLED ${GWAS_STRATA//,/ }; do
    f="${GWAS_DIR}/${s}.sumstats.gz"; [[ "$s" == POOLED ]] && f="${GWAS_DIR}/POOLED.ivw_descriptive.sumstats.gz"
    if [[ -s "$f" ]]; then SUMSTATS_FILES+=("$f"); else SUMSTATS_MISSING+=("${s}"); fi
done
if [[ ${#SUMSTATS_FILES[@]} -eq 0 ]]; then
    echo "ERROR: no sumstats for ${TRAIT}: none of ${GWAS_DIR}/POOLED.ivw_descriptive.sumstats.gz, {${GWAS_STRATA}}.sumstats.gz exists" >&2
    exit 1
fi
if [[ ${#SUMSTATS_MISSING[@]} -gt 0 ]]; then
    echo "Sumstats not found, left out of hit selection: ${SUMSTATS_MISSING[*]} (${GWAS_DIR}/<stratum>.sumstats.gz)"
fi
SUMSTATS="$(IFS=,; echo "${SUMSTATS_FILES[*]}")"

SURV_ARGS=""
if [[ "${MODEL}" == "survival" ]]; then
    SURV_ARGS="--time_col $(trait_time_col "${TRAIT}") --event_col $(trait_event_col "${TRAIT}")"
fi

# Optional: Tractor dosages for local-ancestry-aware epistasis (pooled task only). The ancestry list
# sits next to the Tractor files it describes: data/tractor/ancestries.txt for data/tractor/all_chr.*
TRACTOR_ARGS=""
if [[ "${STRATUM}" == "POOLED" && -d "${DATA_DIR}/tractor" ]]; then
    GXG_ANCESTRIES="$(read_stratum_ancestries "${DATA_DIR}/tractor/ancestries.txt")" || exit 1
    TRACTOR_ARGS="--tractor_prefix ${DATA_DIR}/tractor/all_chr --ancestries ${GXG_ANCESTRIES}"
fi

# ---- customisable (env vars) ----------------------------------------------
# numeric phenotype columns; default the study's list (HANDOFF_DIR/phenotypes/covariates.txt) + PC1..PC5
COVARIATES="${COVARIATES:-$(default_covariates)}" || exit 1
ANCESTRY_COL="${ANCESTRY_COL:-$(rr_default ANCESTRY_COL)}"
KNOWN_LOCI="${KNOWN_LOCI:-${SCRIPT_DIR}/assets/known_leukemia_risk_loci.tsv}"
LOCI_ARGS="--known_loci ${KNOWN_LOCI}"
[[ "${GXG_USE_DEFAULT_LOCI:-$(rr_default GXG_USE_DEFAULT_LOCI)}" == "false" ]] && LOCI_ARGS="${LOCI_ARGS} --no_default_loci"
[[ -n "${GXG_CUSTOM_VARIANTS:-}" ]] && LOCI_ARGS="${LOCI_ARGS} --custom_variants ${GXG_CUSTOM_VARIANTS}"   # GRCh38 list
PRUNE_MODE="${GXG_PRUNE_MODE:-$(rr_default GXG_PRUNE_MODE)}" # conditional = test LD-correlated hits conditional on kept hits; variant = strongest only; pair = keep all
COND_P="${GXG_COND_P:-$(rr_default GXG_COND_P)}"                # conditional P below which an LD-correlated hit is an independent signal
# ---------------------------------------------------------------------------

rr_args "${TRAIT} ${MODEL} ${STRATUM}" "${PHENO}" "${SUMSTATS_FILES[@]}"
Rscript "${SCRIPT_DIR}/bin/run_gxg_interaction.R" \
    --sumstats "${SUMSTATS}" \
    --p_threshold "${GXG_P_THRESHOLD:-$(rr_default GXG_P_THRESHOLD)}" \
    --max_hits "${GXG_MAX_HITS:-$(rr_default GXG_MAX_HITS)}" \
    ${LOCI_ARGS} \
    --prune_mode "${PRUNE_MODE}" \
    --cond_p_threshold "${COND_P}" \
    --geno "${DATA_DIR}/genotypes/cohort" \
    --phenotype "${PHENO}" \
    --trait "${TRAIT}" \
    --model "${MODEL}" \
    ${SURV_ARGS} \
    --covariates "${COVARIATES}" \
    --ancestry_col "${ANCESTRY_COL}" \
    --ancestry_config "${SCRIPT_DIR}/bin/ancestry_config.R" \
    --stratum "${STRATUM}" \
    --min_stratum_n "${MIN_STRATUM_N:-$(rr_default MIN_STRATUM_N)}" \
    ${TRACTOR_ARGS} \
    --output_prefix "${RESULTS_DIR}/${TRAIT}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --verbose

echo "Completed GxG for ${TRAIT} / ${STRATUM}"

# Last array task submits the ancestry-heterogeneity job (all strata in one run). It goes to the
# partition and account this task runs under (SLURM_JOB_PARTITION / SLURM_JOB_ACCOUNT are set by
# Slurm in every job) unless SLURM_PARTITION / SLURM_ACCOUNT were exported at launch.
# Unquoted heredoc (<<EOF): every ${...} below is filled in now, when this task submits the job,
# so the new job receives plain values. A trailing \\ becomes a line continuation in that job.
if [[ "${SLURM_ARRAY_TASK_ID}" == "${SLURM_ARRAY_TASK_MAX}" ]]; then
    sbatch --export=ALL \
           --partition="${SLURM_PARTITION:-${SLURM_JOB_PARTITION}}" --account="${SLURM_ACCOUNT:-${SLURM_JOB_ACCOUNT}}" \
           --dependency=afterany:"${SLURM_ARRAY_JOB_ID}" \
           --job-name="gxg_het_${TRAIT}" --time=12:00:00 --mem=32G --cpus-per-task=16 \
           --output="${LOG_DIR:-${SCRIPT_DIR}/logs}/gxg_het_${TRAIT}_%j.out" <<EOF
#!/bin/bash
set -euo pipefail
# This job is not an array task: drop the array variables inherited from the task that submitted it.
unset SLURM_ARRAY_JOB_ID SLURM_ARRAY_TASK_ID SLURM_ARRAY_TASK_COUNT SLURM_ARRAY_TASK_MIN SLURM_ARRAY_TASK_MAX SLURM_ARRAY_TASK_STEP
source "${SCRIPT_DIR}/bin/slurm_utils.sh"
rr_begin gxg_het "${TRAIT} ${MODEL}" "${PHENO}"
source "${SCRIPT_DIR}/bin/activate_env.sh" gwas-r
set_thread_env
Rscript "${SCRIPT_DIR}/bin/run_gxg_interaction.R" \\
    --sumstats '${SUMSTATS}' --p_threshold ${GXG_P_THRESHOLD:-$(rr_default GXG_P_THRESHOLD)} --max_hits ${GXG_MAX_HITS:-$(rr_default GXG_MAX_HITS)} \\
    ${LOCI_ARGS} --prune_mode ${PRUNE_MODE} --cond_p_threshold ${COND_P} \\
    --geno "${DATA_DIR}/genotypes/cohort" \\
    --phenotype "${PHENO}" \\
    --trait ${TRAIT} --model ${MODEL} ${SURV_ARGS} \\
    --covariates ${COVARIATES} --ancestry_col ${ANCESTRY_COL} \\
    --ancestry_config "${SCRIPT_DIR}/bin/ancestry_config.R" --min_stratum_n ${MIN_STRATUM_N:-$(rr_default MIN_STRATUM_N)} \\
    --output_prefix "${RESULTS_DIR}/${TRAIT}.het" --threads 16
EOF
fi
