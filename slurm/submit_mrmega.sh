#!/bin/bash
#SBATCH --job-name=mrmega
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
# --partition, --account, --output and --error come from the sbatch command line
# (slurm/submit_full_pipeline.sh, or the direct-submit form in docs/HPC_SITE.md).

# ============================================================================
# STEP 2b: MR-MEGA across strata, one job per trait (CP-20)            RUN_MRMEGA
# ============================================================================
# K = the strata (the bundle's ancestry/strata.tsv, restricted to STRATA_LIST) whose
# results/gwas/<trait>/<S>.sumstats.gz exists; a stratum skipped below MIN_STRATUM_N has none.
# MR-MEGA regresses effects on T axes of genetic variation and needs K >= T + 2, so
#   pc_axes = min(K - 2, MRMEGA_PCS)          MRMEGA_PCS default 3
#   K >= 3   MR-MEGA with --pc pc_axes
#   K == 2   pc_axes = 0: a fixed-effect meta-analysis of the two strata, NOT MR-MEGA, and the
#            output's header says so ("## pc_axes=0: fixed-effect meta-analysis, not MR-MEGA")
#   K <  2   the trait is excluded (K < pc_axes + 2), and the log says why
#
# Usage (normally submitted by slurm/submit_full_pipeline.sh with RUN_MRMEGA=true):
#   sbatch <site args> slurm/submit_mrmega.sh <trait>
#   DRY_RUN=1 bash slurm/submit_mrmega.sh <trait>     logs K, pc_axes and the decision; runs nothing
# The per-stratum MR-MEGA input (one effect per variant and stratum: BETA_TOTAL / SE_TOTAL) is
# written by bin/sumstats_to_mrmega.py (CP-19), which also checks EA = ALT / NEA = REF against the
# bundle's variant table. Results: results/gwas -> results/mrmega/<trait>/<trait>.mrmega.result.tsv
# (K >= 3) or <trait>.fixed_effect.tsv (K = 2), each with "##" lines saying what it is.
# ============================================================================
set -euo pipefail

# Repo root (Slurm runs a batch script from a spool copy; see docs/HPC_SITE.md)
for SCRIPT_DIR in "${GWAS_REPO:-}" "${SLURM_SUBMIT_DIR:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; do
    [[ -f "${SCRIPT_DIR}/bin/slurm_utils.sh" ]] && break
done
[[ -f "${SCRIPT_DIR}/bin/slurm_utils.sh" ]] || { echo "ERROR: cannot find the repo (submit from the repo root, or export GWAS_REPO)" >&2; exit 1; }
source "${SCRIPT_DIR}/bin/slurm_utils.sh"
DRY_RUN="${DRY_RUN:-0}"
[[ "${DRY_RUN}" == 1 ]] || rr_begin mrmega

TRAIT="${1:?trait required (mrd, efs, os)}"
MRMEGA_PCS="${MRMEGA_PCS:-$(rr_default MRMEGA_PCS)}"
[[ "${MRMEGA_PCS}" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: MRMEGA_PCS=${MRMEGA_PCS}; use a whole number >= 1" >&2; exit 2; }
GWAS_DIR="${GWAS_RESULTS}/gwas/${TRAIT}"
STRATA_RUN="${STRATA_LIST:-$(bundle_strata)}" || exit 1

PRESENT=(); ABSENT=()
for s in ${STRATA_RUN//,/ }; do
    if [[ -s "${GWAS_DIR}/${s}.sumstats.gz" ]]; then PRESENT+=("$s"); else ABSENT+=("$s"); fi
done
K=${#PRESENT[@]}
list() { if [[ $# -gt 0 ]]; then local IFS=,; echo "$*"; else echo none; fi; }
echo "MR-MEGA ${TRAIT}: strata with sumstats: $(list ${PRESENT[@]+"${PRESENT[@]}"}); without: $(list ${ABSENT[@]+"${ABSENT[@]}"})"
if (( K >= 3 )); then
    PC_AXES=$(( K - 2 < MRMEGA_PCS ? K - 2 : MRMEGA_PCS )); MODE=mrmega
    echo "MR-MEGA ${TRAIT}: K=${K} pc_axes=${PC_AXES} (min(K-2, MRMEGA_PCS=${MRMEGA_PCS}))"
elif (( K == 2 )); then
    PC_AXES=0; MODE=fixed_effect
    echo "MR-MEGA ${TRAIT}: K=2 pc_axes=0: fixed-effect meta-analysis, not MR-MEGA"
else
    echo "MR-MEGA ${TRAIT}: excluded: K=${K} < pc_axes + 2 = 2 (at least two strata with sumstats are needed)"
    exit 0
fi
[[ "${DRY_RUN}" == 1 ]] && exit 0

# One effect per variant and stratum: BETA_TOTAL / SE_TOTAL (the total-dose test), as MR-MEGA input,
# EA = ALT and NEA = REF checked against the bundle's variant table (bin/sumstats_to_mrmega.py, CP-19)
CONVERT="${SCRIPT_DIR}/bin/sumstats_to_mrmega.py"
OUT_DIR="${GWAS_RESULTS}/mrmega/${TRAIT}"
mkdir -p "${OUT_DIR}"
BIM="${HANDOFF_DIR:-${SCRIPT_DIR}/data/handoff}/genotypes/cohort.hardcall.maf01.bim"
VAR_ARGS=(); [[ -s "${BIM}" ]] && VAR_ARGS=(--variants "${BIM}")
SS_FILES=(); for s in "${PRESENT[@]}"; do SS_FILES+=("${GWAS_DIR}/${s}.sumstats.gz"); done
rr_args "K=${K} pc_axes=${PC_AXES} mode=${MODE} strata=$(list "${PRESENT[@]}")" "${SS_FILES[@]}"
source "${SCRIPT_DIR}/bin/activate_env.sh" gwas-py
IN_FILES=()
for s in "${PRESENT[@]}"; do
    f="${OUT_DIR}/${s}.mrmega_in.tsv"
    python3 "${CONVERT}" --input "${GWAS_DIR}/${s}.sumstats.gz" --output "$f" ${VAR_ARGS[@]+"${VAR_ARGS[@]}"}
    IN_FILES+=("$f")
done

if [[ "${MODE}" == fixed_effect ]]; then
    python3 "${CONVERT}" --fixed_effect --inputs "$(list "${IN_FILES[@]}")" --strata "$(list "${PRESENT[@]}")" \
        --output "${OUT_DIR}/${TRAIT}.fixed_effect.tsv"
    echo "MR-MEGA ${TRAIT}: K=2: ${OUT_DIR}/${TRAIT}.fixed_effect.tsv (fixed effect, not MR-MEGA)"
    exit 0
fi

# MR-MEGA reads plain tables: it gets copies without the "##" lines (they stay in the converter's
# files), listed one per line; --qt because the inputs carry BETA / SE
LIST="${OUT_DIR}/mrmega_inputs.in"
: > "${LIST}"
for f in "${IN_FILES[@]}"; do
    grep -v '^##' "$f" > "${f%.tsv}.plain.tsv"
    echo "${f%.tsv}.plain.tsv" >> "${LIST}"
done
source "${SCRIPT_DIR}/bin/activate_env.sh" mrmega
MR-MEGA -i "${LIST}" --pc "${PC_AXES}" --qt -o "${OUT_DIR}/${TRAIT}.mrmega"
[[ -s "${OUT_DIR}/${TRAIT}.mrmega.result" ]] || { echo "ERROR: MR-MEGA wrote no ${OUT_DIR}/${TRAIT}.mrmega.result" >&2; exit 1; }
{ printf '## method=MR-MEGA\n## pc_axes=%s\n## K=%s\n## strata=%s\n## beta_source=BETA_TOTAL\n' \
         "${PC_AXES}" "${K}" "$(list "${PRESENT[@]}")"
  cat "${OUT_DIR}/${TRAIT}.mrmega.result"; } > "${OUT_DIR}/${TRAIT}.mrmega.result.tsv"
echo "MR-MEGA ${TRAIT}: ${OUT_DIR}/${TRAIT}.mrmega.result.tsv (K=${K}, pc_axes=${PC_AXES})"
