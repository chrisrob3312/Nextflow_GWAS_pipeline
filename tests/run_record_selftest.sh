#!/bin/bash
# ============================================================================
# run_record_selftest.sh: the B-03 run record end to end, on the cluster, in minutes
# ============================================================================
# Launches like slurm/submit_full_pipeline.sh does (one run_id, params.tsv, jobs.tsv), submits
# two stub jobs (tests/stub_step.sh; the first exits 0, the second exits 3) and the finalise job
# with --dependency=afterany on both. When the finalise job has ended, trace.tsv has two rows,
# each with an exit status, env + yml sha256 and MaxRSS.
#
# From the repo root, with SLURM_PARTITION, SLURM_ACCOUNT, GWAS_ENVS_DIR (and GWAS_SCRATCH) set:
#   tests/run_record_selftest.sh               # prints the run directory, A
#   MAC_MIN=30 tests/run_record_selftest.sh    # a second launch: a new run directory, B
#   diff A/params.tsv B/params.tsv             # exactly one row differs: MAC_MIN
# ============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${REPO}"
export GWAS_REPO="${REPO}"
for v in SLURM_PARTITION SLURM_ACCOUNT; do
    [[ -n "${!v:-}" ]] || { echo "REFUSED: ${v} is not set (docs/HPC_SITE.md)" >&2; exit 2; }
done
export LOG_DIR="${LOG_DIR:-${REPO}/logs}"
mkdir -p "${LOG_DIR}"

source "${REPO}/bin/run_record.sh"
rr_refuse_unknown_switches || exit 2
rr_resolve_switches
RR_DEFERRED="selftest: stub steps only"
GWAS_RUN_ID="$(rr_claim_run_dir "${REPO}/results/run_record")"
GWAS_RUN_DIR="${REPO}/results/run_record/${GWAS_RUN_ID}"
export GWAS_RUN_ID GWAS_RUN_DIR
rr_write_launch_record "${GWAS_RUN_DIR}" "$0" "$@"

SITE=(--parsable --export=ALL --partition="${SLURM_PARTITION}" --account="${SLURM_ACCOUNT}")
j1="$(sbatch "${SITE[@]}" --output="${LOG_DIR}/rr_stub_%j.out" "${REPO}/tests/stub_step.sh" 0)"; j1="${j1%%;*}"
rr_note_job "${GWAS_RUN_DIR}" "stub exit 0" "${j1}"
j2="$(sbatch "${SITE[@]}" --output="${LOG_DIR}/rr_stub_%j.out" "${REPO}/tests/stub_step.sh" 3)"; j2="${j2%%;*}"
rr_note_job "${GWAS_RUN_DIR}" "stub exit 3" "${j2}"
jf="$(sbatch "${SITE[@]}" --dependency=afterany:"${j1}":"${j2}" --output="${LOG_DIR}/run_record_finalize_%j.out" \
      "${REPO}/bin/run_record_finalize.sh" "${GWAS_RUN_DIR}")"; jf="${jf%%;*}"

echo "run directory: ${GWAS_RUN_DIR}"
echo "jobs: stub ${j1} (exit 0), stub ${j2} (exit 3), finalise ${jf} (after both)"
echo "when job ${jf} has ended:  column -t -s \$'\\t' ${GWAS_RUN_DIR}/trace.tsv"
