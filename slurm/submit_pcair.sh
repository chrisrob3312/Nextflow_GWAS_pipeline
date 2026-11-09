#!/bin/bash
#SBATCH --job-name=pcair_pcrelate
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
# --partition, --account, --output and --error come from the sbatch command line
# (slurm/submit_full_pipeline.sh, or the direct-submit form in docs/HPC_SITE.md).

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
# Usage (normally submitted by slurm/submit_full_pipeline.sh; by hand, from the repo root,
# with the site arguments of docs/HPC_SITE.md "Submitting one script by hand"):
#   sbatch <site args> slurm/submit_pcair.sh
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
rr_begin pcair            # run record row first, so a failure from here on is recorded too
source "${SCRIPT_DIR}/bin/activate_env.sh" gwas-r     # conda env gwas-r, TMPDIR on disk (docs/HPC_SITE.md)
print_slurm_info
set_thread_env

DATA_DIR="${SCRIPT_DIR}/data"
OUT_DIR="${SCRIPT_DIR}/results/pcair"
mkdir -p "${OUT_DIR}"

GENO_PREFIX="${GENO:-${DATA_DIR}/genotypes/cohort}"
PHENO_FILE="${PHENO_IN:-${DATA_DIR}/phenotypes/phenotypes.tsv}"
CMD=(Rscript "${SCRIPT_DIR}/bin/genesis_pcair_pcrelate.R"
    --geno "${GENO_PREFIX}"
    --phenotype "${PHENO_FILE}"
    --n_pcs "${N_PCS:-10}"
    --n_pcs_pcrelate "${N_PCS_PCRELATE:-5}"
    --kin_thresh "${KIN_THRESH:-0.0442}"
    --ld_r2 "${LD_R2:-0.1}"
    --iterations "${PCAIR_ITER:-2}"
    --output_prefix "${OUT_DIR}/cohort"
    --threads "${SLURM_CPUS_PER_TASK}")
rr_args "$(print_cmd "${CMD[@]}")" "${GENO_PREFIX}.bed" "${GENO_PREFIX}.bim" "${GENO_PREFIX}.fam" "${PHENO_FILE}"
"${CMD[@]}"

echo "PC-AiR / PC-Relate complete: ${OUT_DIR}"
