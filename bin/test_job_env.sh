#!/bin/bash
#SBATCH --job-name=test_job_env
#SBATCH --array=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
# --partition, --account and --output come from the sbatch command line (docs/HPC_SITE.md).

# ============================================================================
# test_job_env.sh: does a job see what the pipeline relies on? (M-04)
# ============================================================================
# After `source bin/activate_env.sh gwas-r`, inside a one-task array job:
#   1. SLURM_JOB_ID, SLURM_ARRAY_TASK_ID and SLURM_CPUS_PER_TASK are non-empty, in bash and in
#      R (tractor_genesis_adapter.R and genesis_pcair_pcrelate.R read them with Sys.getenv);
#   2. TMPDIR is on disk, not tmpfs, and R's tempdir() is inside it;
#   3. data.table::fread() of a ~2 GB .gz, which first decompresses the file into TMPDIR,
#      finishes without ENOSPC ("No space left on device").
# The .gz is generated in TMPDIR unless TEST_GZ names a real one (e.g. a chromosome's Tractor
# dosage file); TEST_GZ_BYTES sets the generated size (default 2e9).
#
# Usage, from the repo root (with the site settings of docs/HPC_SITE.md exported):
#   sbatch --partition="$SLURM_PARTITION" --account="$SLURM_ACCOUNT" \
#          --output="$PWD/logs/test_job_env_%A_%a.out" bin/test_job_env.sh
# Exit 0 only if every check passed; the last line of the log says PASS or FAIL.
# ============================================================================
set -uo pipefail     # no -e: every check runs and is reported

# Repo root (Slurm runs this file from a spool copy; see docs/HPC_SITE.md)
for REPO in "${GWAS_REPO:-}" "${SLURM_SUBMIT_DIR:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; do
    [[ -f "${REPO}/bin/activate_env.sh" ]] && break
done
[[ -f "${REPO}/bin/activate_env.sh" ]] || { echo "ERROR: cannot find the repo (submit from the repo root, or export GWAS_REPO)" >&2; exit 1; }

source "${REPO}/bin/activate_env.sh" gwas-r || { echo "test_job_env: FAIL (gwas-r did not activate)"; exit 1; }

bad=0
echo "== bash, after conda activate"
for v in SLURM_JOB_ID SLURM_ARRAY_TASK_ID SLURM_CPUS_PER_TASK; do
    if [[ -n "${!v:-}" ]]; then echo "ok    ${v}=${!v}"; else echo "FAIL  ${v} is empty"; bad=1; fi
done
echo "TMPDIR=${TMPDIR} is $(stat -f -c %T "${TMPDIR}"); free there: $(df -Ph "${TMPDIR}" | awk 'NR==2 {print $4}')"
echo "/tmp is $(stat -f -c %T /tmp) (for comparison)"

echo ""
echo "== R"
R_TEST=$(cat <<'EOF'
ok <- TRUE
for (v in c("SLURM_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_CPUS_PER_TASK")) {
    x <- Sys.getenv(v)
    cat(sprintf("%-5s R sees %s=%s\n", if (nzchar(x)) "ok" else "FAIL", v, x))
    if (!nzchar(x)) ok <- FALSE
}
td <- normalizePath(tempdir()); tmp <- normalizePath(Sys.getenv("TMPDIR"))
inside <- startsWith(td, paste0(tmp, "/"))
cat(sprintf("%-5s R tempdir() %s is inside TMPDIR %s\n", if (inside) "ok" else "FAIL", td, tmp))
if (!inside) ok <- FALSE

suppressPackageStartupMessages({ library(data.table); library(R.utils) })
setDTthreads(max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")), na.rm = TRUE))
gz <- Sys.getenv("TEST_GZ")
if (!nzchar(gz)) {
    target <- as.numeric(Sys.getenv("TEST_GZ_BYTES", "2e9"))
    rows <- ceiling(target / 72)          # this table compresses to about 72 bytes per row
    gz <- tempfile(fileext = ".tsv.gz")
    cat(sprintf("writing %.0f rows of random numbers to %s ...\n", rows, gz))
    x <- data.table(CHR = sample.int(22L, rows, TRUE), POS = sample.int(2.5e8, rows, TRUE),
                    matrix(runif(rows * 8), rows))
    fwrite(x, gz, sep = "\t", compress = "gzip")
    rm(x); invisible(gc())
}
cat(sprintf("fread %s (%.2f GB) ...\n", gz, file.size(gz) / 1e9))
t0 <- Sys.time()
res <- tryCatch(fread(gz, showProgress = FALSE), error = function(e) e)
if (inherits(res, "error")) {
    msg <- conditionMessage(res)
    enospc <- grepl("No space left|disk is full", msg, ignore.case = TRUE)
    cat(sprintf("FAIL  fread: %s%s\n", if (enospc) "ENOSPC: " else "", msg))
    ok <- FALSE
} else {
    cat(sprintf("ok    fread read %d rows x %d columns in %.0f s, decompressing into %s\n",
                nrow(res), ncol(res), as.numeric(difftime(Sys.time(), t0, units = "secs")), tempdir()))
}
if (!nzchar(Sys.getenv("TEST_GZ"))) unlink(gz)
quit(save = "no", status = if (ok) 0L else 1L)
EOF
)
Rscript -e "$R_TEST" || bad=1

echo ""
if [[ $bad -eq 0 ]]; then echo "test_job_env: PASS"; else echo "test_job_env: FAIL"; fi
exit "$bad"
