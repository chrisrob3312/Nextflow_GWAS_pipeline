#!/bin/bash
#SBATCH --job-name=structure
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
# --partition, --account, --output and --error come from the sbatch command line
# (slurm/submit_full_pipeline.sh, or the direct-submit form in docs/HPC_SITE.md).

# ============================================================================
# STEP 0: population structure (PCs and GRM) for every model, once per launch
# ============================================================================
# STRUCTURE_SOURCE=bundle (default): ingest what the handoff bundle ships (CP-07;
# docs/HANDOFF_CONTRACT.md rule 4) with bin/ingest_structure.R. Nothing is computed; a missing
# or mismatched file stops the job, naming it. Writes to results/pcair/:
#   cohort.phenotypes.with_pcs.tsv        phenotype + cohort PCs            (POOLED)
#   <S>.phenotypes.with_pcs.tsv           phenotype + within-stratum PCs    (stratum S)
#   cohort.pcrelate.grm.rds / .kinship.rds   links to the shipped files, unchanged (-> --kinship)
# and results/strata/<S>.samples.txt (each stratum's members, from ancestry/strata.tsv).
#
# STRUCTURE_SOURCE=compute: run bin/genesis_pcair_pcrelate.R on the hard calls instead. It is
# kept as the comparator (CP-18); it writes cohort-scale PCs only, which every stratum then uses.
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

HANDOFF_DIR="${HANDOFF_DIR:-${SCRIPT_DIR}/data/handoff}"
OUT_DIR="${SCRIPT_DIR}/results/pcair"
mkdir -p "${OUT_DIR}"
PHENO_FILE="${PHENO_IN:-${HANDOFF_DIR}/phenotypes/phenotypes.tsv}"

case "${STRUCTURE_SOURCE:-$(rr_default STRUCTURE_SOURCE)}" in
    bundle)
        CMD=(Rscript "${SCRIPT_DIR}/bin/ingest_structure.R"
            --handoff_dir "${HANDOFF_DIR}"
            --phenotype "${PHENO_FILE}"
            --out_dir "${OUT_DIR}"
            --strata_dir "${SCRIPT_DIR}/results/strata")
        INPUTS=("${HANDOFF_DIR}/genotypes/cohort.hardcall.maf01.fam" "${HANDOFF_DIR}/ancestry/strata.tsv"
                "${HANDOFF_DIR}/structure/cohort.pcair.pcs.tsv" "${HANDOFF_DIR}/structure/cohort.pcrelate.grm.rds"
                "${HANDOFF_DIR}/structure/cohort.pcrelate.kinship.rds" "${PHENO_FILE}")
        ;;
    compute)
        GENO_PREFIX="${GENO:-${HANDOFF_DIR}/genotypes/cohort.hardcall.maf01}"
        CMD=(Rscript "${SCRIPT_DIR}/bin/genesis_pcair_pcrelate.R"
            --geno "${GENO_PREFIX}"
            --phenotype "${PHENO_FILE}"
            --n_pcs "${N_PCS:-$(rr_default N_PCS)}"
            --n_pcs_pcrelate "${N_PCS_PCRELATE:-$(rr_default N_PCS_PCRELATE)}"
            --kin_thresh "${KIN_THRESH:-$(rr_default KIN_THRESH)}"
            --ld_r2 "${LD_R2:-$(rr_default LD_R2)}"
            --iterations "${PCAIR_ITER:-$(rr_default PCAIR_ITER)}"
            --output_prefix "${OUT_DIR}/cohort"
            --threads "${SLURM_CPUS_PER_TASK:-1}")
        INPUTS=("${GENO_PREFIX}.bed" "${GENO_PREFIX}.bim" "${GENO_PREFIX}.fam" "${PHENO_FILE}")
        ;;
    *)
        echo "REFUSED: STRUCTURE_SOURCE=${STRUCTURE_SOURCE}; use bundle or compute" >&2
        exit 2
        ;;
esac
rr_args "$(print_cmd "${CMD[@]}")" "${INPUTS[@]}"
"${CMD[@]}"

# The comparator knows no strata: write each stratum's members from the bundle's strata.tsv, as the
# ingest does, so the combine step can check them for overlap (CP-22)
if [[ "${STRUCTURE_SOURCE:-$(rr_default STRUCTURE_SOURCE)}" == compute ]]; then
    mkdir -p "${SCRIPT_DIR}/results/strata"
    awk -F'\t' -v d="${SCRIPT_DIR}/results/strata" 'NR > 1 && $2 != "" {print $1 > (d "/" $2 ".samples.txt")}' \
        "${HANDOFF_DIR}/ancestry/strata.tsv"
fi

echo "Structure (${STRUCTURE_SOURCE:-$(rr_default STRUCTURE_SOURCE)}) ready in ${OUT_DIR}"
