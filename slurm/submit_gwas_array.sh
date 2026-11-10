#!/bin/bash
#SBATCH --job-name=tractor_genesis
#SBATCH --array=1-22
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
# --partition, --account, --output and --error come from the sbatch command line
# (slurm/submit_full_pipeline.sh, or the direct-submit form in docs/HPC_SITE.md).

# ============================================================================
# SLURM Array Job: Tractor-GENESIS GWAS by Chromosome
# ============================================================================
# ALL traits go through the same Tractor-GENESIS adapter so results are
# comparable:
#   MRD      -> binary      (GENESIS logistic mixed model null)
#   relapse  -> survival    (Cox null with kinship frailty, time-to-event)
#   OS       -> survival
#
# Usage (normally submitted by slurm/submit_full_pipeline.sh; by hand, from the repo root,
# with the site arguments of docs/HPC_SITE.md "Submitting one script by hand"):
#   sbatch <site args> slurm/submit_gwas_array.sh <stratum> <trait> <model>
#
#   sbatch <site args> slurm/submit_gwas_array.sh LAT1 OS      survival
#   sbatch <site args> slurm/submit_gwas_array.sh LAT1 relapse survival
#   sbatch <site args> slurm/submit_gwas_array.sh LAT1 MRD     binary
#
# Survival traits expect columns <trait>_time and <trait>_status in the
# phenotype file (override with TIME_COL / EVENT_COL env vars).
# The ancestries come from data/tractor/<stratum>/ancestries.txt; a stratum without
# that file is refused (docs/HANDOFF_CONTRACT.md).
# Output: results/gwas/<trait>/<stratum>/chr<N>.tractor_genesis.tsv.gz
#
# DRY_RUN=1 prints the Rscript command one task would run and stops; nothing is
# activated, created or run. Outside an array it shows chromosome 1:
#   DRY_RUN=1 bash slurm/submit_gwas_array.sh EUR OS survival
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
DRY_RUN="${DRY_RUN:-0}"
if [[ "${DRY_RUN}" != 1 ]]; then
    rr_begin gwas         # run record row first, so a failure from here on is recorded too
    source "${SCRIPT_DIR}/bin/activate_env.sh" gwas-r     # conda env gwas-r, TMPDIR on disk (docs/HPC_SITE.md)
    print_slurm_info
    set_thread_env
fi

STRATUM="${1:?stratum required (EUR, AAC, LAT1, LAT2, EAS, SAS, OTHER, POOLED)}"
TRAIT="${2:?trait required (OS, relapse, MRD)}"
MODEL="${3:-survival}"
CHR="${SLURM_ARRAY_TASK_ID:-}"
if [[ -z "${CHR}" ]]; then
    if [[ "${DRY_RUN}" == 1 ]]; then
        CHR=1; echo "(dry run outside an array: showing chromosome 1 of 1-22)"
    else
        echo "ERROR: SLURM_ARRAY_TASK_ID is empty; this script runs as a 22-task array (chromosome = task id)" >&2
        exit 1
    fi
fi

echo "=============================================="
echo "Stratum: ${STRATUM}   Trait: ${TRAIT}   Model: ${MODEL}   Chr: ${CHR}"
echo "=============================================="

DATA_DIR="${SCRIPT_DIR}/data"
RESULTS_DIR="${SCRIPT_DIR}/results/gwas/${TRAIT}/${STRATUM}"
[[ "${DRY_RUN}" == 1 ]] || mkdir -p "${RESULTS_DIR}"

# Prefer PC-AiR / PC-Relate outputs (step 0) so every model adjusts for the
# same ancestry PCs and the same admixture-robust GRM; fall back to data/
PCAIR_DIR="${SCRIPT_DIR}/results/pcair"
if [[ -z "${PHENO:-}" ]]; then
    [[ -f "${PCAIR_DIR}/cohort.phenotypes.with_pcs.tsv" ]] && PHENO="${PCAIR_DIR}/cohort.phenotypes.with_pcs.tsv" \
        || PHENO="${DATA_DIR}/phenotypes/phenotypes.tsv"
fi
if [[ -z "${KINSHIP:-}" ]]; then
    [[ -f "${PCAIR_DIR}/cohort.pcrelate.grm.rds" ]] && KINSHIP="${PCAIR_DIR}/cohort.pcrelate.grm.rds" \
        || KINSHIP="${DATA_DIR}/kinship/kinship_matrix.rds"
fi
echo "Phenotype: ${PHENO}"; echo "GRM: ${KINSHIP}"
TRACTOR_PREFIX="${DATA_DIR}/tractor/${STRATUM}/chr${CHR}"

# Ancestral populations for the Tractor decomposition: this stratum's own list, written next to
# its Tractor files by the extraction step (docs/HANDOFF_CONTRACT.md). No default: a stratum
# without the file stops here, naming the path.
ANC_FILE="${DATA_DIR}/tractor/${STRATUM}/ancestries.txt"
STRATUM_ANCESTRIES="$(read_stratum_ancestries "${ANC_FILE}")" || exit 1
REF_ANCESTRY="EUR"      # the adapter's reference for the LA terms; it must be in the list
if [[ ",${STRATUM_ANCESTRIES}," != *",${REF_ANCESTRY},"* ]]; then
    echo "ERROR: ${ANC_FILE} lists ${STRATUM_ANCESTRIES}, without the reference ancestry ${REF_ANCESTRY}" >&2
    exit 1
fi

# The command as an array: one element per argument, so nothing is split or globbed by accident
CMD=(Rscript "${SCRIPT_DIR}/bin/tractor_genesis_adapter.R"
    --tractor_prefix "${TRACTOR_PREFIX}"
    --phenotype "${PHENO}"
    --trait "${TRAIT}"
    --model "${MODEL}")
if [[ "${MODEL}" == "survival" ]]; then
    CMD+=(--time_col "${TIME_COL:-${TRAIT}_time}" --event_col "${EVENT_COL:-${TRAIT}_status}")
fi
CMD+=(--ancestries "${STRATUM_ANCESTRIES}"
    --ref_ancestry "${REF_ANCESTRY}"
    --kinship "${KINSHIP}"
    --covariates "${COVARIATES:-age,sex,PC1,PC2,PC3,PC4,PC5}"
    --stratum "${STRATUM}"
    --min_stratum_n "${MIN_STRATUM_N:-50}"
    --ancestry_config "${SCRIPT_DIR}/bin/ancestry_config.R"
    --mac_min "${MAC_MIN:-10}"
    --chromosome "${CHR}"
    --output_prefix "${RESULTS_DIR}/chr${CHR}"
    --threads "${SLURM_CPUS_PER_TASK:-1}"
    --verbose)

if [[ "${DRY_RUN}" == 1 ]]; then
    echo "[DRY-RUN] $(print_cmd "${CMD[@]}")"
    exit 0
fi
# Run record: the command line and inputs (phenotype, GRM, ancestry list, this chromosome's Tractor files)
rr_args "$(print_cmd "${CMD[@]}")" "${PHENO}" "${KINSHIP}" "${ANC_FILE}" "${TRACTOR_PREFIX}".*
"${CMD[@]}"

echo "Completed chr${CHR} | ${TRAIT} (${MODEL}) | ${STRATUM}"
