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
#   mrd      -> binary      (GENESIS logistic mixed model null)
#   efs      -> survival    (Cox null with kinship frailty, time-to-event)
#   os       -> survival
#
# Usage (normally submitted by slurm/submit_full_pipeline.sh; by hand, from the repo root,
# with the site arguments of docs/HPC_SITE.md "Submitting one script by hand"):
#   sbatch <site args> slurm/submit_gwas_array.sh <stratum> <trait> <model>
#
#   sbatch <site args> slurm/submit_gwas_array.sh LAT1 os  survival
#   sbatch <site args> slurm/submit_gwas_array.sh LAT1 efs survival
#   sbatch <site args> slurm/submit_gwas_array.sh LAT1 mrd binary
#
# A survival trait reads columns <trait>_time and <trait>_event (efs_time/efs_event, os_time/os_event,
# as docs/HANDOFF_CONTRACT.md names them); TIME_COL_<trait> / EVENT_COL_<trait> override one trait.
# The adapter checks every sample in the stratum's Tractor files against STRATA_FILE
# (default $HANDOFF_DIR/ancestry/strata.tsv) when it exists, and
# skips the stratum (exit 0, no results) below MIN_STRATUM_N analysable samples; the run
# record reads its counts from results/gwas/<trait>/<stratum>/chr<N>.run_counts.tsv.
# Numbered Tractor files (anc0, ...) are read only with ANC_INDEX_MAP (default
# $HANDOFF_DIR/local_ancestry/flare_index_map.tsv).
# The ancestries come from data/tractor/<stratum>/ancestries.txt; a stratum without
# that file is refused (docs/HANDOFF_CONTRACT.md).
# Output: results/gwas/<trait>/<stratum>/chr<N>.tractor_genesis.tsv.gz
#
# DRY_RUN=1 prints the Rscript command one task would run and stops; nothing is
# activated, created or run. Outside an array it shows chromosome 1:
#   DRY_RUN=1 bash slurm/submit_gwas_array.sh EUR os survival
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
TRAIT="${2:?trait required (mrd, efs, os)}"
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
RESULTS_DIR="${GWAS_RESULTS}/gwas/${TRAIT}/${STRATUM}"
[[ "${DRY_RUN}" == 1 ]] || mkdir -p "${RESULTS_DIR}"

# Phenotype with PCs and the GRM come from step 0 (slurm/submit_pcair.sh): this stratum's
# within-stratum PCs (results/pcair/<S>.phenotypes.with_pcs.tsv), the cohort PCs for POOLED, and
# the shipped GRM, unchanged. No silent fallback: without step 0's output the task stops. PHENO /
# KINSHIP set explicitly win. With STRUCTURE_SOURCE=compute only cohort PCs exist, and every
# stratum uses them (logged).
PCAIR_DIR="${SCRIPT_DIR}/results/pcair"
HANDOFF_DIR="${HANDOFF_DIR:-${SCRIPT_DIR}/data/handoff}"
# PC_SCALE (B-05): stratum (default) = the stratum's within-stratum PCs; cohort = the cohort PCs
PC_SCALE="${PC_SCALE:-$(rr_default PC_SCALE)}"
USE_COHORT_PCS=0
[[ "${STRATUM}" == POOLED || "${PC_SCALE}" == cohort ]] && USE_COHORT_PCS=1
if [[ -z "${PHENO:-}" ]]; then
    if [[ ${USE_COHORT_PCS} == 0 && -f "${PCAIR_DIR}/${STRATUM}.phenotypes.with_pcs.tsv" ]]; then
        PHENO="${PCAIR_DIR}/${STRATUM}.phenotypes.with_pcs.tsv"
    elif [[ ( ${USE_COHORT_PCS} == 1 || "${STRUCTURE_SOURCE:-$(rr_default STRUCTURE_SOURCE)}" == compute ) && -f "${PCAIR_DIR}/cohort.phenotypes.with_pcs.tsv" ]]; then
        PHENO="${PCAIR_DIR}/cohort.phenotypes.with_pcs.tsv"
        if [[ "${STRATUM}" != POOLED ]]; then
            if [[ "${PC_SCALE}" == cohort ]]; then echo "Note: PC_SCALE=cohort, so ${STRATUM} uses cohort-scale PCs"
            else echo "Note: STRUCTURE_SOURCE=compute, so ${STRATUM} uses cohort-scale PCs"; fi
        fi
    elif [[ "${DRY_RUN}" == 1 ]]; then
        PHENO="${PCAIR_DIR}/$([[ ${USE_COHORT_PCS} == 1 ]] && echo cohort || echo "${STRATUM}").phenotypes.with_pcs.tsv"
        echo "(dry run: ${PHENO} is written by step 0, which has not run)"
    else
        echo "ERROR: no phenotype with PCs for ${STRATUM} in ${PCAIR_DIR} (run step 0, slurm/submit_pcair.sh, first)" >&2
        exit 1
    fi
fi
if [[ -z "${KINSHIP:-}" ]]; then
    KINSHIP="${PCAIR_DIR}/cohort.pcrelate.grm.rds"
    if [[ ! -e "${KINSHIP}" && "${DRY_RUN}" != 1 ]]; then
        echo "ERROR: no GRM at ${KINSHIP} (run step 0, slurm/submit_pcair.sh, first)" >&2
        exit 1
    fi
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

# Covariates: COVARIATES, else the study's list (HANDOFF_DIR/phenotypes/covariates.txt) + PC1..PC5
COVARIATES="${COVARIATES:-$(HANDOFF_DIR="${HANDOFF_DIR}" default_covariates)}" || exit 1

# The command as an array: one element per argument, so nothing is split or globbed by accident
CMD=(Rscript "${SCRIPT_DIR}/bin/tractor_genesis_adapter.R"
    --tractor_prefix "${TRACTOR_PREFIX}"
    --phenotype "${PHENO}"
    --trait "${TRAIT}"
    --model "${MODEL}")
if [[ "${MODEL}" == "survival" ]]; then
    CMD+=(--time_col "$(trait_time_col "${TRAIT}")" --event_col "$(trait_event_col "${TRAIT}")")
fi
STRATA_FILE="${STRATA_FILE:-${HANDOFF_DIR}/ancestry/strata.tsv}"
if [[ -f "${STRATA_FILE}" ]]; then
    CMD+=(--strata_file "${STRATA_FILE}")
else
    echo "Stratum membership: ${STRATA_FILE} not found; the adapter checks a 'stratum' phenotype column if there is one, else records the check as not run"
fi
ANC_INDEX_MAP="${ANC_INDEX_MAP:-${HANDOFF_DIR}/local_ancestry/flare_index_map.tsv}"
[[ -f "${ANC_INDEX_MAP}" ]] && CMD+=(--ancestry_index_map "${ANC_INDEX_MAP}")
CMD+=(--ancestries "${STRATUM_ANCESTRIES}"
    --ref_ancestry "${REF_ANCESTRY}"
    --kinship "${KINSHIP}"
    --covariates "${COVARIATES}"
    --stratum "${STRATUM}"
    --min_stratum_n "${MIN_STRATUM_N:-$(rr_default MIN_STRATUM_N)}"
    --ancestry_config "${SCRIPT_DIR}/bin/ancestry_config.R"
    --mac_min "${MAC_MIN:-$(rr_default MAC_MIN)}"
    --chromosome "${CHR}"
    --output_prefix "${RESULTS_DIR}/chr${CHR}"
    --threads "${SLURM_CPUS_PER_TASK:-1}"
    --verbose)

if [[ "${DRY_RUN}" == 1 ]]; then
    echo "[DRY-RUN] QC gate: bin/qc_gate.sh --require ${QC_REQUIRE_STRUCTURE} --dir ${GWAS_QC_DIR}"
    echo "[DRY-RUN] $(print_cmd "${CMD[@]}")"
    exit 0
fi
# QC gate (CP-11): the structure step's verdicts must be there and must not have failed
qc_gate_or_exit "${QC_REQUIRE_STRUCTURE}" "GWAS ${TRAIT} / ${STRATUM} chr${CHR}"
# Run record: the command line and inputs (phenotype, GRM, ancestry list, strata, index map, this
# chromosome's Tractor files), and the adapter's counts file (read when the task ends)
INPUTS=("${PHENO}" "${KINSHIP}" "${ANC_FILE}" "${STRATA_FILE}" "${TRACTOR_PREFIX}".*)
[[ -f "${ANC_INDEX_MAP}" ]] && INPUTS+=("${ANC_INDEX_MAP}")
rr_args "$(print_cmd "${CMD[@]}")" "${INPUTS[@]}"
rr_sidecar "${RESULTS_DIR}/chr${CHR}.run_counts.tsv"
"${CMD[@]}"

echo "Completed chr${CHR} | ${TRAIT} (${MODEL}) | ${STRATUM}"
