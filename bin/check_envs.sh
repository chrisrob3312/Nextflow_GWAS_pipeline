#!/bin/bash
#SBATCH --job-name=check_envs
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
# --partition, --account and --output come from the sbatch command line (docs/HPC_SITE.md).

# ============================================================================
# check_envs.sh: activate each of the four conda envs and show it can do its job
# ============================================================================
#   gwas-r   library() on every R package the SLURM path loads (bin/check_r_libs.R)
#   gwas-py  import pandas, numpy, scipy
#   lai      java -version (FLARE needs Java 11 or later); the imports Tractor's scripts use
#   mrmega   MR-MEGA on PATH
#   and envs.tsv still matches envs/*.yml (bin/envs_tsv.sh --check)
#
# Every failure is a hard failure: all checks run, then the script exits 1 if any failed and
# the last line says FAIL. coxme and metafor get their own message, because without them the
# adapter and the meta-analysis still run, quietly, on a weaker model.
#
# Usage, from the repo root (with the site settings of docs/HPC_SITE.md exported):
#   sbatch --partition="$SLURM_PARTITION" --account="$SLURM_ACCOUNT" \
#          --output="$PWD/logs/check_envs_%j.out" bin/check_envs.sh
#   bash bin/check_envs.sh          # the same checks, on a login node
# ============================================================================
set -uo pipefail     # no -e: every check runs, failures are collected and reported at the end

# Repo root (Slurm runs this file from a spool copy; see docs/HPC_SITE.md)
for REPO in "${GWAS_REPO:-}" "${SLURM_SUBMIT_DIR:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; do
    [[ -f "${REPO}/bin/activate_env.sh" ]] && break
done
[[ -f "${REPO}/bin/activate_env.sh" ]] || { echo "ERROR: cannot find the repo (run from the repo root, or export GWAS_REPO)" >&2; exit 1; }

# The R packages the SLURM path loads: the 14 of the plan, plus R.utils, without which
# data.table::fread cannot read a .gz at all.
R_PACKAGES=(GENESIS SeqArray SeqVarTools SNPRelate GWASTools coxme survival metafor
            data.table R.utils optparse argparse qqman ggplot2 Matrix)
PY_GWAS=(pandas numpy scipy)
PY_LAI=(numpy pandas argparse contextlib gzip logging re os)   # Tractor's extract_tracts_flare.py: stdlib only

PY_CHECK='import importlib, sys
bad = 0
for m in sys.argv[1:]:
    try:
        mod = importlib.import_module(m)
        print("ok      %-12s %s" % (m, getattr(mod, "__version__", "")))
    except Exception as e:
        print("MISSING %-12s %s" % (m, e))
        bad = 1
print("python %s at %s" % (sys.version.split()[0], sys.executable))
sys.exit(bad)'

FAILED=()

# check_env ENV COMMAND...: activate ENV in a subshell ( ... ) so the activation does not carry
# over to the next env, then run the check command there. Prints the output; returns its status.
check_env() {
    local env="$1"; shift
    echo ""
    echo "== ${env}"
    (
        source "${REPO}/bin/activate_env.sh" "$env" || exit 10
        "$@"
    )
}

out="$(check_env gwas-r Rscript "${REPO}/bin/check_r_libs.R" "${R_PACKAGES[@]}" 2>&1)"
rc=$?
echo "$out"
if [[ $rc -ne 0 ]]; then
    FAILED+=(gwas-r)
    if grep -q "^MISSING coxme " <<< "$out"; then
        echo "HARD FAIL: coxme is missing; the adapter would fit the Cox null WITHOUT the kinship frailty"
    fi
    if grep -q "^MISSING metafor " <<< "$out"; then
        echo "HARD FAIL: metafor is missing; meta_analyze_strata.R would fall back to fixed-effect IVW"
    fi
fi

check_env gwas-py python3 -c "$PY_CHECK" "${PY_GWAS[@]}" || FAILED+=(gwas-py)

lai_check() {
    local v
    v="$(java -version 2>&1 | head -1)" || { echo "MISSING java: $v"; return 1; }
    echo "ok      java         ${v}"
    python3 -c "$PY_CHECK" "${PY_LAI[@]}"
}
check_env lai lai_check || FAILED+=(lai)

mrmega_check() {
    local bin
    bin="$(command -v MR-MEGA)" || { echo "MISSING MR-MEGA: not on PATH"; return 1; }
    echo "ok      MR-MEGA      ${bin}"
}
check_env mrmega mrmega_check || FAILED+=(mrmega)

echo ""
echo "== envs.tsv"
bash "${REPO}/bin/envs_tsv.sh" --check || FAILED+=(envs.tsv)

echo ""
if [[ ${#FAILED[@]} -eq 0 ]]; then
    echo "check_envs: PASS (gwas-r, gwas-py, lai, mrmega; envs.tsv matches the ymls)"
    exit 0
fi
echo "check_envs: FAIL (${FAILED[*]})"
exit 1
