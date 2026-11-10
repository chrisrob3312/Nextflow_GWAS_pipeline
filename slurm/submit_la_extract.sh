#!/bin/bash
#SBATCH --job-name=la_extract
#SBATCH --array=1-22
#SBATCH --time=6:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
# --partition, --account, --output and --error come from the sbatch command line
# (slurm/submit_full_pipeline.sh, or the direct-submit form in docs/HPC_SITE.md).

# ============================================================================
# FLARE -> Tractor files, one chromosome per array task, in the lai env (CP-16)
# ============================================================================
# Reads the bundle's FLARE output, HANDOFF_DIR/local_ancestry/flare/chr<N>.anc.vcf.gz, and its
# index map, HANDOFF_DIR/local_ancestry/flare_index_map.tsv, and writes the cohort's Tractor files
#   data/tractor/cohort/chr<N>.hapcount.<LABEL>.txt.gz, data/tractor/cohort/chr<N>.dosage.<LABEL>.txt.gz
# with bin/la_flare_to_tractor.sh (Tractor's extract_tracts_flare.py, then labels from the index
# map, then bin/la_tractor_check.py). The extraction is single-threaded: one CPU per task.
# slurm/submit_tractor_layout.sh then splits these per stratum.
#
#   sbatch <site args> slurm/submit_la_extract.sh
#   DRY_RUN=1 bash slurm/submit_la_extract.sh        prints the command for chromosome 1; runs nothing
# ============================================================================
set -euo pipefail
for SCRIPT_DIR in "${GWAS_REPO:-}" "${SLURM_SUBMIT_DIR:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; do
    [[ -f "${SCRIPT_DIR}/bin/slurm_utils.sh" ]] && break
done
[[ -f "${SCRIPT_DIR}/bin/slurm_utils.sh" ]] || { echo "ERROR: cannot find the repo (submit from the repo root, or export GWAS_REPO)" >&2; exit 1; }
source "${SCRIPT_DIR}/bin/slurm_utils.sh"
DRY_RUN="${DRY_RUN:-0}"
if [[ "${DRY_RUN}" != 1 ]]; then
    rr_begin la_extract
    source "${SCRIPT_DIR}/bin/activate_env.sh" lai
fi

CHR="${SLURM_ARRAY_TASK_ID:-1}"
HANDOFF_DIR="${HANDOFF_DIR:-${SCRIPT_DIR}/data/handoff}"
VCF="${HANDOFF_DIR}/local_ancestry/flare/chr${CHR}.anc.vcf.gz"
MAP="${HANDOFF_DIR}/local_ancestry/flare_index_map.tsv"
OUT="${SCRIPT_DIR}/data/tractor/cohort/chr${CHR}"
CMD=(bash "${SCRIPT_DIR}/bin/la_flare_to_tractor.sh" "${VCF}" "${MAP}" "${OUT}")
if [[ "${DRY_RUN}" == 1 ]]; then echo "[DRY-RUN] $(print_cmd "${CMD[@]}")"; exit 0; fi
rr_args "$(print_cmd "${CMD[@]}")" "${VCF}" "${MAP}"
"${CMD[@]}"
echo "Completed chr${CHR}: ${OUT}.{hapcount,dosage}.<label>.txt.gz"
