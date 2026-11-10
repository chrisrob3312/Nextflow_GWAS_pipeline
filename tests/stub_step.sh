#!/bin/bash
#SBATCH --job-name=rr_stub
#SBATCH --time=00:05:00
#SBATCH --mem=500M
#SBATCH --cpus-per-task=1
# --partition, --account and --output come from the sbatch command line.

# A stub pipeline step for the run-record test (tests/run_record_selftest.sh): it goes through
# exactly what a real step does (slurm_utils.sh, activate_env.sh, rr_begin) and then, instead of
# an analysis, runs one python3 line in the gwas-py env and exits with the code it was given.
#   stub_step.sh <exit_code>
set -euo pipefail

# Repo root (Slurm runs this file from a spool copy; see docs/HPC_SITE.md)
for SCRIPT_DIR in "${GWAS_REPO:-}" "${SLURM_SUBMIT_DIR:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; do
    [[ -f "${SCRIPT_DIR}/bin/slurm_utils.sh" ]] && break
done
[[ -f "${SCRIPT_DIR}/bin/slurm_utils.sh" ]] || { echo "ERROR: cannot find the repo (submit from the repo root, or export GWAS_REPO)" >&2; exit 1; }
source "${SCRIPT_DIR}/bin/slurm_utils.sh"
CODE="${1:-0}"
rr_begin stub "exit ${CODE}" "${SCRIPT_DIR}/envs.tsv"
source "${SCRIPT_DIR}/bin/activate_env.sh" gwas-py
python3 -c 'import sys; print("stub step: python", sys.version.split()[0], "in", sys.prefix)'
exit "${CODE}"
