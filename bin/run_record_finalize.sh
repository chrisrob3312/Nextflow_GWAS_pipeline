#!/bin/bash
#SBATCH --job-name=run_record_finalize
#SBATCH --time=00:30:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
# --partition, --account, --output and --dependency=afterany:<every job> come from the launcher.

# ============================================================================
# run_record_finalize.sh: trace.tsv and summary.md for one launch (B-03)
# ============================================================================
# Collects every task's row from <run_dir>/trace.d/, adds State, MaxRSS and Elapsed from sacct,
# and writes <run_dir>/trace.tsv and <run_dir>/summary.md. The launcher submits it with
# --dependency=afterany on every job it submitted, so it runs once all of them have ended,
# however they ended. Safe to run again by hand at any time; it rewrites both files:
#   bin/run_record_finalize.sh results/run_record/<run_id>
# ============================================================================
set -euo pipefail

# Repo root (Slurm runs this file from a spool copy; see docs/HPC_SITE.md)
for REPO in "${GWAS_REPO:-}" "${SLURM_SUBMIT_DIR:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; do
    [[ -f "${REPO}/bin/run_record.sh" ]] && break
done
[[ -f "${REPO}/bin/run_record.sh" ]] || { echo "ERROR: cannot find the repo (run from the repo root, or export GWAS_REPO)" >&2; exit 1; }

source "${REPO}/bin/run_record.sh"
rr_finalize "${1:?usage: run_record_finalize.sh <run_dir>}"
