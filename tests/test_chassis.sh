#!/usr/bin/env bash
# ============================================================================
# test_chassis.sh: local checks for the SLURM / bash chassis
# ============================================================================
# Covers every done-when of CP-01, M-04, CP-02, CP-03, CP-28 and B-03 that can be checked
# without a cluster. No Slurm and no conda needed: fake sbatch, sacct and conda (tests/fakes/)
# go first on PATH, and the conda envs are empty directories with stub programs in them.
# The behavioural checks run on a copy of the repo in a temporary directory, so nothing is
# written into this checkout. Light and sequential: no R analysis runs.
#
#   tests/test_chassis.sh [work_dir]        exit 0 only if every check passed
# ============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/.." && pwd)"
W="${1:-$(mktemp -d)}"; mkdir -p "$W"; W="$(cd "$W" && pwd)"
pass=0; failn=0
# Under pipefail, `producer | grep -q` can fail at random: grep -q exits at the first match and the
# producer, if still writing, dies of SIGPIPE. Multi-line checks below match against "$(producer)".
ok()      { echo "  ok   $1"; pass=$((pass + 1)); }
bad()     { echo "  FAIL $1"; failn=$((failn + 1)); }
check()   { local name="$1"; shift; if "$@" >/dev/null 2>&1; then ok "$name"; else bad "$name"; fi; }
section() { echo ""; echo "== $*"; }
has()     { grep -qF -- "$2" <<< "$1"; }     # has "$output" "literal text"

# --- a clean environment: nothing from the caller's shell may change a result
for v in SLURM_PARTITION SLURM_ACCOUNT SLURM_TMPDIR SLURM_JOB_ID SLURM_ARRAY_JOB_ID SLURM_ARRAY_TASK_ID \
         SLURM_CPUS_PER_TASK SLURM_SUBMIT_DIR GWAS_REPO GWAS_RUN_ID GWAS_RUN_DIR GWAS_ENV LOG_DIR DRY_RUN \
         GWAS_CONDA_MODULE CONDA_PREFIX SBATCH_PARTITION SBATCH_ACCOUNT STRATUM_MIN_N PC_SCALE MRMEGA_PCS BUNDLE_RULES; do
    unset "$v"
done
while IFS='|' read -r name _; do unset "$name"; done < <(bash -c "source '${REPO}/bin/run_record.sh'; printf '%s\n' \"\${RR_SWITCHES[@]}\"")
unset -f conda 2>/dev/null

# --- the copy of the repo, fixtures, fakes
C="${W}/repo"
rm -rf "${C:?}" "${W:?}/slurm" "${W:?}/envs" "${W:?}/scratch"; mkdir -p "$C"
tar -C "$REPO" --exclude=./.git --exclude=./results --exclude=./logs --exclude=./data --exclude=./work -cf - . | tar -C "$C" -xf -
REAL_RSCRIPT="$(command -v Rscript || true)"     # before the fakes go on PATH
export FAKE_SLURM_DIR="${W}/slurm"; mkdir -p "$FAKE_SLURM_DIR"
export PATH="${C}/tests/fakes:${PATH}"
export GWAS_ENVS_DIR="${W}/envs" GWAS_SCRATCH="${W}/scratch"
mkdir -p "${W}/envs/"{gwas-r,gwas-py,lai,mrmega}/bin "${W}/scratch"
# Stub Rscript inside the gwas-r env: emulates bin/check_r_libs.R (packages in FAKE_R_MISSING are
# missing), and stands in for every analysis script (prints its arguments, exits FAKE_RSCRIPT_EXIT).
cat > "${W}/envs/gwas-r/bin/Rscript" <<EOF
#!/usr/bin/env bash
# Scripts that need only base R + data.table run in the real R (when there is one)
case "\${1:-}" in */ingest_structure.R|*/qc_assert.R|*/meta_analyze_strata.R) [[ -n "${REAL_RSCRIPT}" ]] && exec "${REAL_RSCRIPT}" "\$@" ;; esac
EOF
cat >> "${W}/envs/gwas-r/bin/Rscript" <<'EOF'
if [[ "${1:-}" == */check_r_libs.R ]]; then
    shift; bad=0
    for p in "$@"; do
        if [[ " ${FAKE_R_MISSING:-} " == *" $p "* ]]; then echo "MISSING $(printf '%-12s' "$p") there is no package called '$p'"; bad=1
        else echo "ok      $(printf '%-12s' "$p") 0.0"; fi
    done
    exit "$bad"
fi
echo "stub Rscript: $*"
# As the adapter: write <output_prefix>.run_counts.tsv (status from FAKE_ADAPTER_STATUS)
if [[ "${1:-}" == */tractor_genesis_adapter.R ]]; then
    prev=""; for a in "$@"; do [[ "$prev" == --output_prefix ]] && out="$a"; prev="$a"; done
    st="${FAKE_ADAPTER_STATUS:-ok}"
    sp=""; for a in "$@"; do [[ "$sp" == --stratum && ",${FAKE_SKIP_STRATUM:-}," == *",$a,"* ]] && st=skipped_below_min_stratum_n; sp="$a"; done
    if [[ "$st" == ok ]]; then nt=900; else nt=NA; fi
    printf 'key\tvalue\nstatus\t%s\nstratum\tEUR\nstratum_check\tpassed\nmin_stratum_n\t50\nn_null\t%s\nn_variants\t1000\nn_pass_mac\t950\nn_tested\t%s\nn_failed\t2\n' \
        "$st" "${FAKE_N_NULL:-120}" "$nt" > "${out}.run_counts.tsv"
    # and, when it ran, a genomic-order table with its "## covariates=" line (FAKE_META_CHR<n> replaces
    # that line for chromosome n)
    nout=0; sp=""; for a in "$@"; do [[ "$sp" == --stratum && ",${FAKE_NO_OUTPUT:-}," == *",$a,"* ]] && nout=1; sp="$a"; done
    if [[ "$st" == ok && $nout == 0 ]]; then
        cov=""; chr=""; prev=""; for a in "$@"; do [[ "$prev" == --covariates ]] && cov="$a"; [[ "$prev" == --chromosome ]] && chr="$a"; prev="$a"; done
        mv="FAKE_META_CHR${chr}"
        printf '%s\nCHR\tPOS\tID\tREF\tALT\tN\tP_JOINT\tMAC_EUR\tBETA_EUR\tSE_EUR\tP_EUR\n%s\t100\t%s:100\tA\tG\t120\t0.5\t30\t0.1\t0.05\t0.04\n' \
            "${!mv:-## covariates=${cov}}" "$chr" "$chr" | gzip > "${out}.tractor_genesis.genomic_order.tsv.gz"
    fi
    # and its QC verdicts (bin/qc_assert.R): PASS, INFO for the model when skipped, FAIL for the id
    # named in FAKE_ADAPTER_QC_FAIL
    mkdir -p "${GWAS_QC_DIR}/assertions.d"
    for id in adapter_stratum_check adapter_sample_join adapter_grm_cover adapter_model_n; do
        v=PASS; [[ "$id" == adapter_model_n && "$st" != ok ]] && v=INFO; [[ "$id" == "${FAKE_ADAPTER_QC_FAIL:-}" ]] && v=FAIL
        printf '%s\t%s\t%s\tstub adapter %s\n' "$id" "$v" "${FAKE_N_NULL:-120}" "$out" > "${GWAS_QC_DIR}/assertions.d/${id}.stub.$$.tsv"
    done
fi
# As the structure step without a real R: file its verdict
if [[ "${1:-}" == */ingest_structure.R || "${1:-}" == */genesis_pcair_pcrelate.R ]]; then
    mkdir -p "${GWAS_QC_DIR}/assertions.d"
    printf 'pcair_pc_merge\t%s\t8\tstub structure step\n' "${FAKE_PC_MERGE:-PASS}" > "${GWAS_QC_DIR}/assertions.d/pcair_pc_merge.stub.$$.tsv"
    printf 'duplicates_removed\tPASS\t0\tstub structure step\n' > "${GWAS_QC_DIR}/assertions.d/duplicates_removed.stub.$$.tsv"
fi
exit "${FAKE_RSCRIPT_EXIT:-0}"
EOF
printf '#!/bin/sh\necho "openjdk version \\"21.0.10\\" (stub)" >&2\n' > "${W}/envs/lai/bin/java"
# stub MR-MEGA: with -o <root>, writes <root>.result (a header and one row) as MR-MEGA does
cat > "${W}/envs/mrmega/bin/MR-MEGA" <<'EOF'
#!/bin/sh
echo "MR-MEGA (stub) $*"
o=""; while [ $# -gt 0 ]; do [ "$1" = "-o" ] && o="$2"; shift; done
[ -n "$o" ] && printf 'MarkerName\tChromosome\tPosition\tEA\tNEA\tNcohort\tP-value_association\n1:1:A:G\t1\t1\tG\tA\t3\t0.5\n' > "$o.result"
exit 0
EOF
for e in gwas-r lai; do printf '#!/bin/sh\necho "bcftools 1.25 (stub)"\necho "Using htslib 1.25"\n' > "${W}/envs/${e}/bin/bcftools"; done
cp "${C}/tests/toy/fake_extract_tracts_flare.py" "${W}/envs/lai/bin/extract_tracts_flare.py"   # Tractor's script, stood in for
chmod +x "${W}/envs/gwas-r/bin/Rscript" "${W}/envs/lai/bin/java" "${W}/envs/mrmega/bin/MR-MEGA" "${W}"/envs/*/bin/bcftools "${W}/envs/lai/bin/extract_tracts_flare.py"
STRATA_DEFAULT=(EUR AAC LAT1 LAT2 OTHER)
for s in "${STRATA_DEFAULT[@]}"; do mkdir -p "${C}/data/tractor/${s}"; done
printf 'EUR\n' > "${C}/data/tractor/EUR/ancestries.txt"
printf 'AFR\nEUR\n' > "${C}/data/tractor/AAC/ancestries.txt"
for s in LAT1 LAT2 OTHER; do printf 'AFR\nAMR\nEUR\n' > "${C}/data/tractor/${s}/ancestries.txt"; done
printf 'AFR\nAMR\nEUR\n' > "${C}/data/tractor/ancestries.txt"       # cohort-wide list (GxG pooled, PRS)
# the stub bundle's strata (ancestry/strata.tsv: EUR AFR AMR MIX), the launcher's default STRATA_LIST
for s in AFR AMR MIX; do mkdir -p "${C}/data/tractor/${s}"; done
printf 'AFR\nEUR\n' > "${C}/data/tractor/AFR/ancestries.txt"
for s in AMR MIX; do printf 'AFR\nAMR\nEUR\n' > "${C}/data/tractor/${s}/ancestries.txt"; done
# the handoff bundle (HANDOFF_DIR default <repo>/data/handoff): the stub the checker's tests use
bash "${C}/scripts/test/make_stub_bundle.sh" "${C}/data/handoff" >/dev/null 2>&1 || echo "  (make_stub_bundle.sh failed)"
L="${C}/slurm/submit_full_pipeline.sh"
G="${C}/slurm/submit_gwas_array.sh"
SITE=(SLURM_PARTITION=testpart SLURM_ACCOUNT=testacct)

echo "work dir: ${W}"

# ============================================================================
section "syntax: bash -n on every chassis script"
for f in "$REPO"/slurm/*.sh "$REPO"/bin/*.sh "$REPO"/shell/*.sh "$REPO"/tests/*.sh "$REPO"/tests/fakes/*; do
    check "bash -n ${f#"$REPO"/}" bash -n "$f"
done
if command -v shellcheck >/dev/null 2>&1; then
    for f in "$REPO"/slurm/*.sh "$REPO"/bin/*.sh; do check "shellcheck ${f#"$REPO"/}" shellcheck -S error "$f"; done
else
    echo "  (shellcheck is not installed: skipped)"
fi

# ============================================================================
section "CP-02 / CP-28 greps (on this checkout)"
out="$(grep -rn 'partition=normal\|--output=logs/' "$REPO/slurm/")"
[[ -z "$out" ]] && ok "grep partition=normal / --output=logs/ in slurm/: nothing" || bad "grep found: $out"
out="$(grep -n ANCESTRIES "$REPO/bin/slurm_utils.sh")"
[[ -z "$out" ]] && ok "grep ANCESTRIES bin/slurm_utils.sh: nothing" || bad "grep found: $out"

section "CP-01: every Rscript/python3 call is preceded by a source of bin/activate_env.sh"
# Within each file, and within each heredoc job, the activate_env.sh line must come before every
# non-comment line calling Rscript or python3.
out="$(cd "$REPO" && grep -rln 'Rscript\|python3' slurm/ bin/slurm_utils.sh | while read -r f; do
    awk -v f="$f" '
        /<<.?EOF/ { act = 0 }                                   # a heredoc job starts afresh
        /activate_env\.sh/ && !/^[[:space:]]*#/ { act = 1 }
        /Rscript|python3/ && !/^[[:space:]]*#/ && !/activate_env/ { n++; if (!act) print f ":" FNR ": " $0 }
        END { if (n == 0) print f ": no call found?" }' "$f"; done)"
[[ -z "$out" ]] && ok "every call has an activate_env.sh source before it" || bad "not preceded: ${out}"
n_calls="$(cd "$REPO" && grep -rn 'Rscript\|python3' slurm/ bin/slurm_utils.sh | grep -cv '^[^:]*:[0-9]*:[[:space:]]*#')"
echo "       (${n_calls} lines checked)"

section "CP-01: envs.tsv and the env specs"
check "envs.tsv: four envs" test "$(tail -n +2 "$REPO/envs.tsv" | wc -l)" -eq 4
check "envs.tsv: yml_sha256 matches envs/*.yml (bin/envs_tsv.sh --check)" bash "$REPO/bin/envs_tsv.sh" --check
for e in gwas-r gwas-py lai mrmega; do check "envs/${e}.yml exists, prefix uses \${GWAS_ENVS_DIR}" grep -qP "^${e}\t\\\$\{GWAS_ENVS_DIR\}/${e}\t" "$REPO/envs.tsv"; done
for p in r-base bioconductor-genesis bioconductor-seqarray bioconductor-seqvartools bioconductor-snprelate \
         bioconductor-gwastools r-coxme r-survival r-metafor r-data.table r-optparse r-argparse r-qqman r-ggplot2 r-matrix r-r.utils \
         r-coloc r-susier bcftools; do
    grep -q "^  - ${p}[=<>]" "$REPO/envs/gwas-r.yml" || bad "gwas-r.yml lacks ${p}"
done
ok "gwas-r.yml lists the 15 planned packages + r-r.utils, r-coloc, r-susier, bcftools"
check "lai.yml has bcftools (the extraction step's site restriction)" grep -q "^  - bcftools=" "$REPO/envs/lai.yml"

# ============================================================================
section "CP-01 / M-04: bin/activate_env.sh (fake conda)"
AE="${C}/bin/activate_env.sh"
out="$(bash -c 'set -eu; source "$1" gwas-py >/dev/null; echo "opts=$- prefix=$CONDA_PREFIX env=$GWAS_ENV sha=$GWAS_ENV_YML_SHA256 tmp=$TMPDIR"' _ "$AE" 2>&1)"
has "$out" "prefix=${W}/envs/gwas-py" && ok "activates the prefix from envs.tsv" || bad "activation: $out"
[[ "$out" == *"opts="*e*u* ]] && ok "caller's set -e and set -u are back on afterwards" || bad "options after: $out"
has "$out" "sha=$(awk -F'\t' '$1 == "gwas-py" {print $3}' "$REPO/envs.tsv")" && ok "exports GWAS_ENV_YML_SHA256 from envs.tsv" || bad "sha: $out"
has "$out" "tmp=${W}/scratch/tmp" && ok "TMPDIR = \$GWAS_SCRATCH/tmp when the site sets no SLURM_TMPDIR" || bad "TMPDIR: $out"
out="$(bash -c 'source "$1" gwas-py >/dev/null; echo "opts=$-"' _ "$AE" 2>&1)"
[[ "$out" != *e* && "$out" != *u* ]] && ok "caller without set -e/-u stays without them" || bad "options leaked: $out"
mkdir -p "${W}/nodetmp"
out="$(SLURM_TMPDIR="${W}/nodetmp" bash -c 'source "$1" gwas-py >/dev/null; echo "tmp=$TMPDIR"' _ "$AE" 2>&1)"
has "$out" "tmp=${W}/nodetmp" && ok "TMPDIR = \$SLURM_TMPDIR when the site sets it" || bad "SLURM_TMPDIR: $out"
if [[ "$(stat -f -c %T /dev/shm 2>/dev/null)" == tmpfs ]]; then
    out="$(env -u GWAS_SCRATCH TMPDIR=/dev/shm bash -c 'source "$1" gwas-py; echo "rc=$?"' _ "$AE" 2>&1)"
    has "$out" "rc=1" && has "$out" "/dev/shm is tmpfs" && ok "refuses a TMPDIR on tmpfs, saying so" || bad "tmpfs: $out"
fi
mv "${W}/envs/lai" "${W}/envs/lai.away"
out="$(bash -c 'source "$1" lai; echo "rc=$? caller continues"' _ "$AE" 2>&1)"
has "$out" "env 'lai'" && has "$out" "${W}/envs/lai does not exist" && ok "missing prefix: names the env and the path" || bad "missing prefix: $out"
has "$out" "rc=1 caller continues" && ok "missing prefix: returns 1, does not exit the caller's shell" || bad "return: $out"
out="$(bash -c 'set -e; source "$1" lai; echo SHOULD-NOT-PRINT' _ "$AE" 2>&1)"; rc=$?
[[ $rc -ne 0 && "$out" != *SHOULD-NOT-PRINT* ]] && ok "under set -e the job stops at the failed activation" || bad "set -e: rc=$rc $out"
mv "${W}/envs/lai.away" "${W}/envs/lai"
out="$(bash "$AE" gwas-r 2>&1)"; rc=$?
[[ $rc -eq 2 ]] && has "$out" "must be sourced" && ok "refuses to be run instead of sourced" || bad "run directly: rc=$rc $out"
out="$(env -u GWAS_ENVS_DIR bash -c 'source "$1" gwas-r' _ "$AE" 2>&1)"
has "$out" "GWAS_ENVS_DIR is not set" && ok "names GWAS_ENVS_DIR when it is needed and unset" || bad "GWAS_ENVS_DIR: $out"

section "CP-01: bin/check_envs.sh (stub programs in the fake envs)"
out="$(bash "${C}/bin/check_envs.sh" 2>&1)"; rc=$?
[[ $rc -eq 0 ]] && has "$out" "check_envs: PASS" && ok "all present: exit 0, PASS" || bad "all present: rc=$rc $(tail -3 <<< "$out")"
out="$(FAKE_R_MISSING="coxme metafor" bash "${C}/bin/check_envs.sh" 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && ok "coxme and metafor missing: non-zero exit" || bad "coxme/metafor missing: rc=$rc"
has "$out" "HARD FAIL: coxme" && has "$out" "HARD FAIL: metafor" && ok "coxme and metafor named as hard failures" || bad "messages: $out"
out="$(FAKE_R_MISSING="qqman" bash "${C}/bin/check_envs.sh" 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "check_envs: FAIL (gwas-r)" && ok "any other missing R package also fails" || bad "qqman: rc=$rc"
out="$(FAKE_R_MISSING="hyprcoloc" bash "${C}/bin/check_envs.sh" 2>&1)"; rc=$?
[[ $rc -eq 0 ]] && has "$out" "note    hyprcoloc    not installed" && ok "hyprcoloc missing: noted, not a failure (installed by hand, RUN_COLOC only)" || bad "hyprcoloc note: rc=$rc"
out="$(FAKE_R_MISSING="coloc" bash "${C}/bin/check_envs.sh" 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "MISSING coloc" && ok "coloc missing: fails" || bad "coloc missing: rc=$rc"
mv "${W}/envs/lai/bin/bcftools" "${W}/bcftools.away"
out="$(bash "${C}/bin/check_envs.sh" 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "MISSING bcftools" && has "$out" "check_envs: FAIL (lai)" && ok "bcftools missing from lai: fails (a copy elsewhere on PATH does not count)" || bad "bcftools lai: rc=$rc $(tail -2 <<< "$out")"
mv "${W}/bcftools.away" "${W}/envs/lai/bin/bcftools"
mv "${W}/envs/mrmega/bin/MR-MEGA" "${W}/MR-MEGA.away"
out="$(bash "${C}/bin/check_envs.sh" 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "MISSING MR-MEGA" && ok "MR-MEGA not on PATH in mrmega: fails" || bad "MR-MEGA: rc=$rc"
mv "${W}/MR-MEGA.away" "${W}/envs/mrmega/bin/MR-MEGA"
if command -v Rscript >/dev/null 2>&1; then
    out="$(Rscript "${C}/bin/check_r_libs.R" utils NoSuchPackageXYZ 2>&1)"; rc=$?
    [[ $rc -eq 1 ]] && has "$out" "ok      utils" && has "$out" "MISSING NoSuchPackageXYZ" \
        && ok "bin/check_r_libs.R with real R: reports ok / MISSING, exits 1" || bad "check_r_libs.R: rc=$rc $out"
fi

section "M-04: bin/test_job_env.sh (small generated .gz, real R)"
if command -v Rscript >/dev/null 2>&1 && Rscript -e 'quit(status = !all(sapply(c("data.table", "R.utils"), requireNamespace, quietly = TRUE)))' >/dev/null 2>&1; then
    mv "${W}/envs/gwas-r/bin/Rscript" "${W}/Rscript.stub"            # use the real R for this one
    out="$(SLURM_JOB_ID=7 SLURM_ARRAY_TASK_ID=1 SLURM_CPUS_PER_TASK=2 TEST_GZ_BYTES=2e5 bash "${C}/bin/test_job_env.sh" 2>&1)"; rc=$?
    [[ $rc -eq 0 ]] && has "$out" "test_job_env: PASS" && has "$out" "R sees SLURM_ARRAY_TASK_ID=1" && has "$out" "is inside TMPDIR ${W}/scratch/tmp" \
        && ok "prints the three SLURM variables from bash and R; fread decompresses inside TMPDIR" || bad "test_job_env: rc=$rc $(tail -5 <<< "$out")"
    out="$(SLURM_JOB_ID=7 SLURM_CPUS_PER_TASK=2 TEST_GZ_BYTES=2e5 bash "${C}/bin/test_job_env.sh" 2>&1)"; rc=$?
    [[ $rc -ne 0 ]] && has "$out" "FAIL  SLURM_ARRAY_TASK_ID is empty" && ok "an empty SLURM_ARRAY_TASK_ID fails the test" || bad "empty var: rc=$rc"
    mv "${W}/Rscript.stub" "${W}/envs/gwas-r/bin/Rscript"
else
    echo "  (R with data.table and R.utils not available here: skipped)"
fi

# ============================================================================
section "CP-02: launcher refuses without SLURM_PARTITION / SLURM_ACCOUNT"
out="$(cd /tmp && env SLURM_PARTITION=testpart "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "REFUSED: SLURM_ACCOUNT is not set" && ok "SLURM_ACCOUNT unset: refused" || bad "account: rc=$rc $out"
out="$(cd /tmp && env SLURM_ACCOUNT=testacct "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "REFUSED: SLURM_PARTITION is not set" && ok "SLURM_PARTITION unset: refused" || bad "partition: rc=$rc $out"

section "CP-02 / CP-03: default dry run from /tmp"
out="$(cd /tmp && env "${SITE[@]}" "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -eq 0 ]] && ok "dry run exits 0" || bad "dry run: rc=$rc $(tail -5 <<< "$out")"
sb="$(grep '\[DRY-RUN\]' <<< "$out")"
nsb="$(wc -l <<< "$sb")"
[[ "$(grep -c -- '--partition=testpart --account=testacct' <<< "$sb")" -eq "$nsb" ]] && ok "partition and account from the environment on all ${nsb} sbatch lines" || bad "site args missing on some line"
logs="$(grep -o -- '--\(output\|error\)=[^ ]*' <<< "$sb")"
[[ -n "$logs" ]] && ! grep -qv -- "=${C}/logs/" <<< "$logs" && ok "every --output/--error is absolute, under ${C}/logs" || bad "log paths: $(grep -v -- "=${C}/logs/" <<< "$logs" | head -3)"
has "$out" "deferred: PRS COLOC VALIDATION GXG MRMEGA (RUN_*=false)" && ok "prints 'deferred: PRS COLOC VALIDATION GXG MRMEGA (RUN_*=false)'" || bad "deferred line missing"
grep -q 'submit_prs_array\|colocalization\|PRS validation\|submit_gxg_array' <<< "$sb" && bad "a deferred step was submitted" || ok "no sbatch for PRS, coloc, validation or GxG"
has "$sb" "run record finalise: sbatch" && grep -q -- '--dependency=afterany:dry1:dry2:' <<< "$(grep 'finalise' <<< "$sb")" && ok "finalise job waits on every job (afterany)" || bad "finalise line"
[[ ! -e "${C}/results" ]] && ok "a dry run creates no run record" || bad "dry run wrote ${C}/results"

section "CP-03: RUN_* gates"
out="$(cd /tmp && env "${SITE[@]}" RUN_PRS=true "$L" --dry-run 2>&1)"
grep -q 'submit_prs_array.sh os' <<< "$(grep '\[DRY-RUN\]' <<< "$out")" && ok "RUN_PRS=true prints the PRS sbatch" || bad "RUN_PRS=true: no PRS sbatch"
has "$out" "deferred: COLOC VALIDATION GXG MRMEGA (RUN_*=false)" && ok "RUN_PRS=true: deferred list shrinks" || bad "deferred with RUN_PRS"
out="$(cd /tmp && env "${SITE[@]}" RUN_GXG=true "$L" --dry-run 2>&1)"
gx="$(grep '\[DRY-RUN\] GxG' <<< "$out")"
[[ "$(grep -c -- '--dependency=afterok:dry14:dry15:dry16 ' <<< "$gx")" -eq 3 ]] && grep -q -- '--array=0-4 ' <<< "$gx" \
    && ok "RUN_GXG without PRS: GxG waits on the three combine jobs; one task for POOLED + each of the 4 bundle strata" || bad "GxG dependency: $gx"
grep -q -- 'afterok: \|afterok::\|afterok:$\|--dependency= ' <<< "$out" && bad "an empty job id in a dependency" || ok "no empty job id in any dependency"
out="$(cd /tmp && env "${SITE[@]}" RUN_PRS=true RUN_COLOC=true RUN_VALIDATION=true RUN_GXG=true RUN_MRMEGA=true "$L" --dry-run 2>&1)"
has "$out" "deferred: none" && [[ "$(grep -c '\[DRY-RUN\] \(PRS\|colocalization\|PRS validation\|GxG\)' <<< "$out")" -eq 10 ]] && ok "all RUN_*=true: every step submitted (coloc once per trait)" || bad "all on"
out="$(cd /tmp && env "${SITE[@]}" RUN_VALIDATION=true "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "RUN_VALIDATION=true needs RUN_PRS=true" && ok "RUN_VALIDATION without RUN_PRS: refused" || bad "validation alone: rc=$rc"
out="$(cd /tmp && env "${SITE[@]}" RUN_PRS=yes "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "RUN_PRS=yes; use true or false" && ok "RUN_PRS=yes: refused" || bad "RUN_PRS=yes: rc=$rc"
out="$(cd /tmp && env "${SITE[@]}" STRATUM_MIN_N=50 "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "the scripts read MIN_STRATUM_N" && ok "STRATUM_MIN_N (read by no script): refused" || bad "STRATUM_MIN_N: rc=$rc"

section "Traits and survival columns follow the handoff contract (mrd; efs, os with <trait>_time / <trait>_event)"
check "defaults: TRAITS_LIST mrd,efs,os and SURVIVAL_TRAITS efs,os" bash -c "grep -q '\"TRAITS_LIST|mrd,efs,os|' '$REPO/bin/run_record.sh' && grep -q '\"SURVIVAL_TRAITS|efs,os|' '$REPO/bin/run_record.sh'"
out="$(cd /tmp && env "${SITE[@]}" "$L" --dry-run 2>&1)"
grep -q 'GWAS mrd (binary) / EUR' <<< "$out" && grep -q 'GWAS efs (survival) / EUR' <<< "$out" && grep -q 'GWAS os (survival) / EUR' <<< "$out" \
    && ok "dry run: mrd binary, efs and os survival" || bad "trait models: $(grep 'GWAS ' <<< "$out" | head -3)"
out="$(cd /tmp && DRY_RUN=1 bash "$G" EUR efs survival 2>&1)"
grep -q -- '--time_col efs_time --event_col efs_event ' <<< "$out" && ok "efs reads efs_time / efs_event" || bad "efs columns: $(tail -1 <<< "$out")"
out="$(cd /tmp && DRY_RUN=1 TIME_COL_efs=efs_days bash "$G" EUR efs survival 2>&1)"
grep -q -- '--time_col efs_days --event_col efs_event ' <<< "$out" && ok "TIME_COL_efs overrides efs only" || bad "override: $(tail -1 <<< "$out")"
out="$(cd /tmp && env "${SITE[@]}" TIME_COL=x "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "set TIME_COL_<trait>" && ok "a global TIME_COL is refused (each trait has its own)" || bad "TIME_COL: rc=$rc"
out="$(cd /tmp && env "${SITE[@]}" TRAITS_LIST=mrd,dfs SURVIVAL_TRAITS=dfs "$L" 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "dfs: column dfs_time is not in" && ok "a launch refuses a trait whose columns the phenotype file lacks" || bad "missing columns: rc=$rc $(tail -2 <<< "$out")"
rm -rf "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR "$L" 2>&1)"
RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"
grep -qP '^SURVIVAL_COLUMNS\tefs:efs_time/efs_event,os:os_time/os_event\tderived$' "${RD%/}/params.tsv" && ok "params.tsv records each survival trait's columns" || bad "params: $(tail -2 "${RD%/}/params.tsv")"
rm -rf "${C:?}/results" "${C:?}/logs"

# ============================================================================
section "CP-28: submit_gwas_array.sh reads data/tractor/<stratum>/ancestries.txt"
out="$(cd /tmp && DRY_RUN=1 bash "$G" EUR os survival 2>&1)"; rc=$?
[[ $rc -eq 0 ]] && grep -q -- '--ancestries EUR ' <<< "$(grep '^\[DRY-RUN\]' <<< "$out")" && ok "ancestries.txt 'EUR': prints --ancestries EUR" || bad "EUR: rc=$rc $out"
out="$(cd /tmp && DRY_RUN=1 SLURM_ARRAY_TASK_ID=5 bash "$G" AAC mrd binary 2>&1)"
grep -q -- '--ancestries AFR,EUR .*--chromosome 5 ' <<< "$(grep '^\[DRY-RUN\]' <<< "$out")" && ok "two labels, one per line: --ancestries AFR,EUR (chr from the task id)" || bad "AAC: $out"
mv "${C}/data/tractor/LAT1/ancestries.txt" "${W}/anc.away"
out="$(cd /tmp && DRY_RUN=1 bash "$G" LAT1 os survival 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "${C}/data/tractor/LAT1/ancestries.txt does not exist" && ok "without the file: exits non-zero naming it" || bad "absent: rc=$rc $out"
mv "${W}/anc.away" "${C}/data/tractor/LAT1/ancestries.txt"
mv "${C}/data/tractor/AMR/ancestries.txt" "${W}/anc.away"
out="$(cd /tmp && env "${SITE[@]}" "$L" --dry-run 2>&1)"
has "$out" "WARNING: missing ${C}/data/tractor/AMR/ancestries.txt" && ok "launcher dry run warns about the missing file" || bad "launcher warning"
out="$(cd /tmp && env "${SITE[@]}" "$L" 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "REFUSED: missing ${C}/data/tractor/AMR/ancestries.txt" && ok "launcher refuses a real launch without it" || bad "launcher refusal: rc=$rc"
mv "${W}/anc.away" "${C}/data/tractor/AMR/ancestries.txt"
rm -rf "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && DRY_RUN=1 MAC_MIN=30 MIN_STRATUM_N=50 bash "$G" EUR os survival 2>&1)"
grep -q -- '--min_stratum_n 50 .*--mac_min 30 ' <<< "$(grep '^\[DRY-RUN\]' <<< "$out")" && ok "MAC_MIN and MIN_STRATUM_N reach the adapter" || bad "switches: $out"
check "shell/README.md says 'debug-only, not the execution path'" grep -q "debug-only, not the execution path" "$REPO/shell/README.md"
for f in "$REPO"/shell/*.sh; do check "$(basename "$f"): one-line DEBUG-ONLY header" bash -c "sed -n 2p '$f' | grep -q '^# DEBUG-ONLY, not the execution path'"; done

# ============================================================================
section "B-03 / B-05: one default per switch, in RR_SWITCHES; scripts read it with rr_default"
out="$(bash -c 'source "$1/bin/run_record.sh"; printf "%s\n" "${RR_SWITCHES[@]}"' _ "$C" | while IFS='|' read -r name def _; do
    [[ "$def" == "@job" ]] && continue
    grep -rhoE "\\\$\{${name}:?-[^}]*\}" "$C/slurm/" "$C/bin/slurm_utils.sh" | sort -u | while read -r m; do
        got="${m#*-}"; got="${got%\}}"
        [[ "$got" == "\$(rr_default ${name})" || ( -z "$got" && -z "$def" ) ]] || echo "${name}: script default '${got}' instead of \$(rr_default ${name})"
    done
done)"
[[ -z "$out" ]] && ok "no script carries its own default for a switch: every \${NAME:-...} is \$(rr_default NAME)" || bad "literal defaults: $out"
out="$(bash -c 'source "$1/bin/run_record.sh"; printf "%s\n" "${RR_SWITCHES[@]}"' _ "$C" | awk -F'|' 'NF != 4 || $4 == "" {print $1}')"
[[ -z "$out" ]] && ok "every RR_SWITCHES entry has a domain (NAME|DEFAULT|what|DOMAIN)" || bad "no domain: $out"

section "MIN_STRATUM_N = 50 (decided 10 Nov 2026); MAC_MIN stays 10, marked open"
check "RR_SWITCHES: MIN_STRATUM_N default 50" grep -q '"MIN_STRATUM_N|50|' "$REPO/bin/run_record.sh"
check "RR_SWITCHES: MAC_MIN default 10, noted as open against the contract's 20" grep -q '"MAC_MIN|10|.*open: the contract says 20' "$REPO/bin/run_record.sh"
check "adapter --min_stratum_n default 50" grep -q 'min_stratum_n"), type = "integer", default = 50' "$REPO/bin/tractor_genesis_adapter.R"
check "GxG --min_stratum_n default 50" grep -q '"--min_stratum_n", type = "integer", default = 50' "$REPO/bin/run_gxg_interaction.R"
check "ancestry_config.R MIN_STRATUM_N <- 50" grep -q '^MIN_STRATUM_N <- 50' "$REPO/bin/ancestry_config.R"
check "nextflow.config min_stratum_n = 50" grep -qE '^ +min_stratum_n += 50 ' "$REPO/nextflow.config"
check "contract and ENV_MRMEGA name the same floor (50)" bash -c "grep -q 'floor: 50 per stratum' '$REPO/docs/HANDOFF_CONTRACT.md' && grep -q 'floor (50, the same as .MIN_STRATUM_N' '$REPO/docs/ENV_MRMEGA.md'"
out="$(awk '/source\(opt\$ancestry_config\)/ {s = NR} /^MIN_STRATUM_N <- opt\$min_stratum_n/ {a = NR} END {print (a > s && s > 0) ? "after" : "before"}' "$REPO/bin/run_gxg_interaction.R")"
[[ "$out" == after ]] && ok "GxG: --min_stratum_n is applied after ancestry_config.R, so the option wins" || bad "GxG: ancestry_config.R overrides --min_stratum_n"

section "B-03: a real launch (fake sbatch, jobs not run)"
rm -rf "${FAKE_SLURM_DIR:?}"/*
out="$(cd /tmp && env "${SITE[@]}" "$L" 2>&1)"; rc=$?
[[ $rc -eq 0 ]] && ok "launch exits 0" || bad "launch: rc=$rc $(tail -5 <<< "$out")"
RUN1="$(ls -d "${C}"/results/run_record/*/ 2>/dev/null | head -1)"; RUN1="${RUN1%/}"
[[ -n "$RUN1" && "$(basename "$RUN1")" =~ ^[0-9]{8}-[0-9]{4}_[0-9a-z]+_1$ ]] && ok "run directory $(basename "$RUN1") (date-minute, git sha, suffix)" || bad "run dir: '$RUN1'"
for f in params.tsv run_info.tsv jobs.tsv; do check "${f} written at launch" test -s "${RUN1}/${f}"; done
check "params.tsv: MIN_STRATUM_N 50 (default)" grep -qP '^MIN_STRATUM_N\t50\tdefault$' "${RUN1}/params.tsv"
n_sb="$(wc -l < "${FAKE_SLURM_DIR}/sbatch.log")"
[[ "$(tail -n +2 "${RUN1}/jobs.tsv" | wc -l)" -eq $((n_sb - 1)) ]] && ok "jobs.tsv lists the $((n_sb - 1)) jobs finalise reports on" || bad "jobs.tsv rows"
grep -qP "^finalise_job\\t$(tail -1 "${FAKE_SLURM_DIR}/sbatch.log" | sed 's/^id=\([0-9]*\).*/\1/')\$" "${RUN1}/run_info.tsv" && ok "run_info.tsv names the finalise job" || bad "finalise_job in run_info"
[[ "$(grep -c -- '--partition=testpart --account=testacct' "${FAKE_SLURM_DIR}/sbatch.log")" -eq "$n_sb" ]] && ok "every real sbatch call carries partition and account" || bad "sbatch.log site args"
last="$(tail -1 "${FAKE_SLURM_DIR}/sbatch.log")"
ids="$(awk '{sub(/^id=/, "", $1); print $1}' "${FAKE_SLURM_DIR}/sbatch.log" | head -n -1 | paste -sd: -)"
has "$last" "--dependency=afterany:${ids} " && has "$last" "run_record_finalize.sh ${RUN1}" && ok "finalise submitted last, afterany on all ${n_sb} - 1 others" || bad "finalise: $last"
heredocs=0; badh=""
for j in "${FAKE_SLURM_DIR}"/jobs/*.sh; do
    head -1 "$j" | grep -q '^#!/bin/bash' || continue
    grep -q '^SCRIPT_DIR="${GWAS_REPO}"$' "$j" || continue      # the launcher's heredoc jobs
    heredocs=$((heredocs + 1))
    bash -n "$j" 2>/dev/null || badh+=" $(basename "$j"):syntax"
    a="$(grep -n 'activate_env.sh' "$j" | head -1 | cut -d: -f1)"; r="$(grep -n '^Rscript' "$j" | head -1 | cut -d: -f1)"
    [[ -z "$r" || ( -n "$a" && "$a" -lt "$r" ) ]] || badh+=" $(basename "$j"):activate-after-Rscript"
    grep -q 'partition=normal\|output=logs/' "$j" && badh+=" $(basename "$j"):site-path"
done
[[ $heredocs -gt 0 && -z "$badh" ]] && ok "the ${heredocs} heredoc job scripts sbatch received: valid bash, activate before Rscript" || bad "heredoc jobs (${heredocs}):${badh}"

section "B-03: a second launch with MAC_MIN=30"
before="$(cd "$RUN1" && find . -type f -exec md5sum {} + | sort)"
out="$(cd /tmp && env "${SITE[@]}" MAC_MIN=30 "$L" 2>&1)"; rc=$?
RUN2="$(ls -d "${C}"/results/run_record/*/ | sort | tail -1)"; RUN2="${RUN2%/}"
[[ $rc -eq 0 && "$RUN2" != "$RUN1" ]] && ok "new run directory $(basename "$RUN2")" || bad "second launch: rc=$rc $RUN2"
d="$(diff "${RUN1}/params.tsv" "${RUN2}/params.tsv" | grep '^[<>]')"
[[ "$(wc -l <<< "$d")" -eq 2 ]] && has "$d" $'< MAC_MIN\t10\tdefault' && has "$d" $'> MAC_MIN\t30\tenvironment' \
    && ok "params.tsv differs in exactly one row: MAC_MIN 10 (default) -> 30 (environment)" || bad "params diff: $d"
after="$(cd "$RUN1" && find . -type f -exec md5sum {} + | sort)"
[[ "$before" == "$after" ]] && ok "the first run directory is untouched" || bad "first run directory changed"
has "$out" "Switches you set: MAC_MIN=30" && ok "the launch log names the switch you set" || bad "launch log switches"

section "B-03: two launches in the same minute on the same commit do not collide"
out="$(bash -c 'source "$1/bin/run_record.sh"; rr_claim_run_dir "$2" 20990101-0000; rr_claim_run_dir "$2" 20990101-0000' _ "$C" "${W}/claim" 2>&1)"
a="$(sed -n 1p <<< "$out")"; b="$(sed -n 2p <<< "$out")"
[[ "$a" == 20990101-0000_*_1 && "$b" == 20990101-0000_*_2 && -d "${W}/claim/$a" && -d "${W}/claim/$b" ]] && ok "same stamp twice: ${a} then ${b}" || bad "claim: $out"

section "B-03: two stub jobs through the run record (tests/run_record_selftest.sh, jobs run by the fake sbatch)"
rm -rf "${FAKE_SLURM_DIR:?}"/*
out="$(cd "$C" && env "${SITE[@]}" FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 bash tests/run_record_selftest.sh 2>&1)"; rc=$?
SELF="$(sed -n 's/^run directory: //p' <<< "$out")"
[[ $rc -eq 0 && -s "${SELF}/trace.tsv" ]] && ok "self-test ran; trace.tsv written by the finalise job" || bad "self-test: rc=$rc $out"
T="${SELF}/trace.tsv"
[[ "$(tail -n +2 "$T" | wc -l)" -eq 2 ]] && ok "trace.tsv has two rows" || bad "rows: $(cat "$T")"
col() { awk -F'\t' -v c="$1" 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == c) k = i; next} {print $k}' "$T" | paste -sd' ' -; }
[[ "$(col exit_status)" == "0 3" || "$(col exit_status)" == "3 0" ]] && ok "exit status recorded: $(col exit_status)" || bad "exit_status: $(col exit_status)"
[[ "$(col env)" == "gwas-py gwas-py" ]] && ok "env recorded: gwas-py" || bad "env: $(col env)"
sha="$(awk -F'\t' '$1 == "gwas-py" {print $3}' "$REPO/envs.tsv")"
[[ "$(col yml_sha256)" == "${sha} ${sha}" ]] && ok "yml sha256 from envs.tsv on both rows" || bad "yml_sha256: $(col yml_sha256)"
[[ "$(col max_rss)" == "2048K 2048K" ]] && ok "MaxRSS from sacct on both rows" || bad "max_rss: $(col max_rss)"
[[ "$(col script_sha256 | awk '{print length($1)}')" -eq 64 ]] && ok "script sha256 recorded" || bad "script_sha256: $(col script_sha256)"
has "$(col inputs_md5)" "envs.tsv=$(md5sum "${C}/envs.tsv" | cut -d' ' -f1)" && ok "input md5 recorded" || bad "inputs_md5: $(col inputs_md5)"
[[ "$(col n_tested)" == "NA NA" && "$(col task_status)" == "NA NA" ]] && ok "no counts file for a stub step: n_* and task_status NA" || bad "n_tested: $(col n_tested)"
grep -q '| stub | 2 | 1 | 1 | 0 | 0 |' "${SELF}/summary.md" && ok "summary.md: stub step, 2 tasks, 1 ok, 1 non-zero" || bad "summary: $(cat "${SELF}/summary.md")"
out2="$(cd "$C" && env "${SITE[@]}" MAC_MIN=30 FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 bash tests/run_record_selftest.sh 2>&1)"
SELF2="$(sed -n 's/^run directory: //p' <<< "$out2")"
d="$(diff "${SELF}/params.tsv" "${SELF2}/params.tsv" | grep -c '^[<>]')"
[[ "$SELF2" != "$SELF" && "$d" -eq 2 ]] && ok "self-test again with MAC_MIN=30: new directory, params.tsv differs in one row" || bad "self-test 2: $SELF2 d=$d"

section "B-03: failures that happen before the analysis still leave a row"
rm -rf "${FAKE_SLURM_DIR:?}"/*; mkdir -p "${W}/rr_fail/trace.d"
mv "${W}/envs/gwas-py" "${W}/gwas-py.away"
out="$(cd "$C" && env "${SITE[@]}" GWAS_REPO="$C" GWAS_RUN_DIR="${W}/rr_fail" FAKE_SBATCH_RUN=1 \
        sbatch --parsable --output="${W}/rr_fail.out" tests/stub_step.sh 0 2>&1)"
mv "${W}/gwas-py.away" "${W}/envs/gwas-py"
row="$(cat "${W}"/rr_fail/trace.d/stub.*.tsv 2>/dev/null)"
[[ "$(cut -f14 <<< "$row")" == 1 && "$(cut -f7 <<< "$row")" == none ]] && ok "env missing: the task's row says exit 1, env none" || bad "activation failure row: $row"
mkdir -p "${W}/rr_term/trace.d"
GWAS_RUN_DIR="${W}/rr_term" SLURM_JOB_ID=4242 bash -c 'source "$1/bin/slurm_utils.sh"; rr_begin term_test; sleep 5 & wait $!' _ "$C" >/dev/null 2>&1 &
tp=$!; for _ in 1 2 3 4 5 6 7 8 9 10; do [[ -s "${W}/rr_term/trace.d/term_test.4242.tsv" ]] && break; sleep 0.2; done
kill -TERM "$tp"; wait "$tp" 2>/dev/null
[[ "$(cut -f14 "${W}/rr_term/trace.d/term_test.4242.tsv")" == 143 ]] && ok "SIGTERM (time limit, scancel): the row says exit 143" || bad "SIGTERM row: $(cat "${W}/rr_term/trace.d/term_test.4242.tsv")"
# A job the launcher submitted: task 1 ran and wrote a row, task 2 died before rr_begin, tasks 3-22
# were cancelled while pending. Finalise lists all three kinds.
mkdir -p "${W}/rr_sacct/trace.d"
printf 'label\tjob_id\nGWAS X\t9001\n' > "${W}/rr_sacct/jobs.tsv"
printf 'gwas\t9001_1\t9100\t9001\t1\thost\tgwas-r\tabc\tdef\targs\tNA\tt0\tt1\t0\tNA\tNA\tNA\tNA\n' > "${W}/rr_sacct/trace.d/gwas.9001_1.tsv"
printf '9001_1|COMPLETED|0:0||00:01:00\n9001_1.batch|COMPLETED|0:0|1024K|00:01:00\n9001_2|OUT_OF_MEMORY|0:9||00:00:30\n9001_[3-22]|CANCELLED by 0|0:0||00:00:00\n' > "${FAKE_SLURM_DIR}/sacct.db"
RR_SACCT_WAIT=0 bash "${C}/bin/run_record_finalize.sh" "${W}/rr_sacct" >/dev/null 2>&1
T="${W}/rr_sacct/trace.tsv"
awk -F'\t' 'NR > 1 && $2 == "9001_1" && $14 == "0" && $20 == "1024K"' "$T" | grep -q . && ok "the task that ran: exit 0 from the task, MaxRSS from sacct" || bad "ran: $(cat "$T")"
awk -F'\t' 'NR > 1 && $2 == "9001_2" && $14 == "137" && $19 == "OUT_OF_MEMORY" && $22 == "sacct"' "$T" | grep -q . && ok "the task killed before rr_begin: exit 137 (signal 9) from sacct" || bad "oom: $(cat "$T")"
awk -F'\t' 'NR > 1 && $2 == "9001_[3-22]" && $14 == "" && $19 ~ /^CANCELLED/' "$T" | grep -q . && ok "the cancelled tasks: listed, no exit status, state CANCELLED" || bad "cancelled: $(cat "$T")"

section "B-03: end to end, one stratum x one trait, every job run by the fake sbatch"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=os FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 "$L" 2>&1)"; rc=$?
E2E="$(ls -d "${C}"/results/run_record/*/ 2>/dev/null | head -1)"; E2E="${E2E%/}"
T="${E2E}/trace.tsv"
[[ $rc -eq 0 && -s "$T" ]] && ok "launch ran PC-AiR, 22 GWAS tasks, combine and finalise" || bad "e2e: rc=$rc $(tail -5 <<< "$out")"
steps="$(awk -F'\t' 'NR > 1 {n[$1]++} END {for (s in n) printf "%s=%d ", s, n[s]}' "$T" | tr ' ' '\n' | sort | paste -sd' ' -)"
[[ "$steps" == "combine=1 gwas=22 pcair=1" ]] && ok "trace rows per step: ${steps}" || bad "steps: ${steps}"
[[ "$(awk -F'\t' 'NR > 1 && $14 != "0"' "$T" | wc -l)" -eq 0 ]] && ok "every task exit 0, recorded by the task itself" || bad "exits: $(awk -F'\t' 'NR>1{print $1, $14}' "$T" | sort | uniq -c)"
grep -q -- '--ancestries EUR .*--mac_min 10 ' <<< "$(awk -F'\t' 'NR > 1 && $1 == "gwas" {print $10}' "$T")" && ok "gwas rows carry the adapter command line" || bad "gwas args"
grep -q 'ERROR: cannot find the repo' "${C}"/logs/*.out 2>/dev/null && bad "a job could not find the repo from the spool copy" || ok "every job found the repo while running from a spool copy"
[[ "$(grep -l 'activate_env: gwas-r' "${C}"/logs/gwas_os_EUR_*.out 2>/dev/null | wc -l)" -eq 22 ]] && ok "all 22 GWAS tasks activated gwas-r" || bad "activation in GWAS logs"

section "CP-02: one script submitted by hand (no launcher, no GWAS_REPO)"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"; mkdir -p "${C}/logs"
out="$(cd "$C" && env "${SITE[@]}" FAKE_SBATCH_RUN=1 sbatch --parsable --partition=testpart --account=testacct \
        --output="${C}/logs/pcair_%j.out" slurm/submit_pcair.sh 2>&1)"
jid="$(tail -1 <<< "$out")"
grep -q 'structure ingested from bundle\|stub Rscript: .*ingest_structure.R' "${C}/logs/pcair_${jid}.out" 2>/dev/null \
    && ok "from the repo root: the job finds the repo through SLURM_SUBMIT_DIR and runs" || bad "direct submit: $(cat "${C}/logs/pcair_${jid}.out" 2>&1 | tail -3)"
ls "${C}/results/run_record/direct_${jid}/trace.d/"pcair.*.tsv >/dev/null 2>&1 && ok "its trace row goes to results/run_record/direct_${jid}/" || bad "direct trace row"
out="$(cd /tmp && env "${SITE[@]}" FAKE_SBATCH_RUN=1 sbatch --parsable --output="${W}/pcair_elsewhere.out" "${C}/slurm/submit_pcair.sh" 2>&1)"
grep -q 'cannot find the repo' "${W}/pcair_elsewhere.out" && ok "from another directory: stops with 'cannot find the repo'" || bad "elsewhere: $(cat "${W}/pcair_elsewhere.out")"

# ============================================================================
section "Adapter: --stratum, --min_stratum_n, no positional Tractor files (bin/tractor_io.R)"
if command -v Rscript >/dev/null 2>&1 && Rscript -e 'quit(status = !requireNamespace("data.table", quietly = TRUE))' >/dev/null 2>&1; then
    out="$(Rscript "$REPO/tests/test_tractor_io.R" "$REPO" 2>&1)"; rc=$?
    while IFS= read -r line; do
        case "$line" in "ok   "*) ok "tractor_io.R: ${line#ok   }" ;; "FAIL "*) bad "tractor_io.R: ${line#FAIL }" ;; esac
    done <<< "$out"
    [[ $rc -eq 0 ]] || bad "tests/test_tractor_io.R exited $rc: $(tail -3 <<< "$out")"
    out="$(Rscript -e 'for (f in commandArgs(TRUE)) invisible(parse(f)); cat("parsed\n")' \
           "$REPO"/bin/tractor_genesis_adapter.R "$REPO"/bin/calculate_prs_admixed.R "$REPO"/bin/run_gxg_interaction.R "$REPO"/bin/tractor_io.R 2>&1)"
    has "$out" parsed && ok "the changed R scripts parse" || bad "R parse: $out"
else
    echo "  (R with data.table not available here: tractor_io.R checks skipped)"
fi
A="$REPO/bin/tractor_genesis_adapter.R"
grep -q 'anc_code' "$A" && bad "adapter still derives a file number from the position in --ancestries" || ok "adapter: no file number taken from the position in --ancestries"
grep -q 'source(file.path(.adapter_dir, "tractor_io.R"))' "$A" && ok "adapter sources bin/tractor_io.R" || bad "adapter does not source tractor_io.R"
grep -q 'check_stratum_membership(opt$stratum, tractor_ids, opt$strata_file, pheno)' "$A" && ok "adapter: --stratum checked against the Tractor files' samples" || bad "adapter: no stratum check"
grep -q 'if (n_null < opt$min_stratum_n)' "$A" && grep -q 'quit(save = "no", status = 0)' "$A" && ok "adapter: below --min_stratum_n, skip with exit 0" || bad "adapter: no min_stratum_n gate"
grep -q 'run_chunk <- function(ch) test_chunk(ch, null_ops, ancestries, ref_anc, opt$mac_min, opt$stratum)' "$A" && grep -q 'if (!is.null(stratum)) row$STRATUM <- stratum' "${A%/*}/tractor_model.R" \
    && ok "adapter: outputs labelled with a STRATUM column (the first; result_row in bin/tractor_model.R)" || bad "adapter: no STRATUM label"
grep -q 'write_run_counts(opt$output_prefix, run_counts)' "$A" && ok "adapter: writes <output_prefix>.run_counts.tsv" || bad "adapter: no counts file"
grep -q '".ancdose.", k' "$REPO/bin/calculate_prs_admixed.R" && bad "PRS still reads Tractor files by position" || ok "PRS LA partial: no file number from the position either"
for m in modules/local/tractor_genesis_adapter/main.nf modules/local/prs_admixed/main.nf; do
    grep -q -- '--ancestry_index_map ancestry_index_map.tsv' "$REPO/$m" && ok "${m}: passes its declared order as an explicit index map" || bad "${m}: no index map"
done
grep -q 'EUR,AFR,AMR' "$REPO/slurm/submit_gxg_array.sh" "$REPO/slurm/submit_prs_array.sh" && bad "GxG/PRS still hard-code EUR,AFR,AMR" || ok "GxG and PRS no longer hard-code EUR,AFR,AMR"
grep -q '^TRAITS=' "$REPO/bin/slurm_utils.sh" && bad "slurm_utils.sh still sets TRAITS" || ok "slurm_utils.sh no longer sets or exports TRAITS"

section "GWAS array: strata file and index map reach the adapter when present"
out="$(cd /tmp && DRY_RUN=1 STRATA_FILE=/nonexistent/strata.tsv bash "$G" EUR os survival 2>&1)"
has "$out" "/nonexistent/strata.tsv not found" && ! grep -q -- '--strata_file' <<< "$(grep '^\[DRY-RUN\]' <<< "$out")" && ok "no strata file: not passed, and the log says the check falls back" || bad "no strata file: $out"
out="$(cd /tmp && DRY_RUN=1 bash "$G" EUR os survival 2>&1 | grep '^\[DRY-RUN\]')"
has "$out" "--strata_file ${C}/data/handoff/ancestry/strata.tsv" && has "$out" "--ancestry_index_map ${C}/data/handoff/local_ancestry/flare_index_map.tsv" && has "$out" "--min_stratum_n 50 " \
    && ok "the bundle's strata.tsv and flare_index_map.tsv (HANDOFF_DIR) and MIN_STRATUM_N 50 passed" || bad "args: $out"

section "CP-07: step 0 ingests the bundle's structure (no PC-AiR); refuses rather than recomputes"
if [[ -n "${REAL_RSCRIPT}" ]] && "${REAL_RSCRIPT}" -e 'quit(status = !requireNamespace("data.table", quietly = TRUE))' >/dev/null 2>&1 \
   && ! grep -q "placeholder rds" "${C}/data/handoff/structure/cohort.pcrelate.grm.rds"; then
    rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"; mkdir -p "${C}/logs"
    pcair() { (cd "$C" && env "${SITE[@]}" GWAS_REPO="$C" FAKE_SBATCH_RUN=1 "$@" sbatch --parsable --output="${C}/logs/pcair_%j.out" slurm/submit_pcair.sh 2>&1 | head -1); }
    jid="$(pcair)"; lg="${C}/logs/pcair_${jid}.out"
    grep -q "structure ingested from bundle" "$lg" && ! grep -q "genesis_pcair_pcrelate" "$lg" && grep -q "^${jid}|COMPLETED|0:0" "${FAKE_SLURM_DIR}/sacct.db" \
        && ok "exit 0, logs 'structure ingested from bundle', PC-AiR not run" || bad "ingest: $(tail -5 "$lg")"
    grep -q "dimnames == .fam IID" "$lg" && ok "GRM and kinship dimnames checked against the .fam IIDs" || bad "dimnames check: $(grep -i dimnames "$lg")"
    for S in AFR AMR EUR MIX; do [[ -s "${C}/results/pcair/${S}.phenotypes.with_pcs.tsv" && -s "${C}/results/strata/${S}.samples.txt" ]] || bad "no ${S} outputs"; done
    eur_h="$(head -1 "${C}/results/pcair/EUR.phenotypes.with_pcs.tsv")"; coh_h="$(head -1 "${C}/results/pcair/cohort.phenotypes.with_pcs.tsv")"
    [[ "$eur_h" == *$'\tPC1\tPC2\tPC3\tPC4\tPC5' && "$coh_h" == *$'\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6' && "$(wc -l < "${C}/results/pcair/EUR.phenotypes.with_pcs.tsv")" -eq 4 ]] \
        && ok "per stratum: its members (EUR: 3) with its own within-stratum PCs; cohort PCs for POOLED" || bad "with_pcs: $eur_h / $coh_h"
    [[ "$(cat "${C}/results/strata/EUR.samples.txt" | paste -sd, -)" == "S001,S002,S003" ]] && ok "results/strata/<S>.samples.txt from strata.tsv" || bad "samples.txt"
    [[ "$(readlink "${C}/results/pcair/cohort.pcrelate.grm.rds")" == "${C}/data/handoff/structure/cohort.pcrelate.grm.rds" ]] \
        && ok "the GRM is the shipped file, linked unchanged (the adapter's --kinship)" || bad "grm link: $(ls -l "${C}/results/pcair/")"
    out="$(cd /tmp && DRY_RUN=1 bash "$G" EUR os survival 2>&1 | grep '^\[DRY-RUN\]')"
    has "$out" "--phenotype ${C}/results/pcair/EUR.phenotypes.with_pcs.tsv" && has "$out" "--kinship ${C}/results/pcair/cohort.pcrelate.grm.rds" \
        && ok "GWAS EUR task: within-stratum phenotype and the shipped GRM" || bad "gwas inputs: $out"
    out="$(cd /tmp && DRY_RUN=1 bash "$G" POOLED os survival 2>&1)"   # POOLED has no ancestries.txt in the fixture; only check the phenotype choice
    has "$out" "Phenotype: ${C}/results/pcair/cohort.phenotypes.with_pcs.tsv" && ok "POOLED task: cohort PCs" || bad "pooled pheno: $(grep Phenotype <<< "$out")"
    mv "${C}/data/handoff/structure/strata/AMR.pcair.pcs.tsv" "${W}/AMR.away"
    jid="$(pcair)"; lg="${C}/logs/pcair_${jid}.out"
    grep -q "REFUSED: structure ingest: missing or empty .*structure/strata/AMR.pcair.pcs.tsv" "$lg" && grep -q "^${jid}|FAILED|2:0" "${FAKE_SLURM_DIR}/sacct.db" \
        && ok "a stratum without its PC file: refused, naming the file" || bad "missing stratum PCs: $(tail -3 "$lg")"
    mv "${W}/AMR.away" "${C}/data/handoff/structure/strata/AMR.pcair.pcs.tsv"
    mv "${C}/data/handoff/structure" "${W}/structure.away"
    jid="$(pcair)"; lg="${C}/logs/pcair_${jid}.out"
    grep -q "REFUSED: structure ingest: missing or empty .*structure/cohort.pcair.pcs.tsv" "$lg" && ! grep -q "genesis_pcair_pcrelate" "$lg" \
        && ok "structure/ removed: refused, naming the file, nothing recomputed" || bad "no structure: $(tail -3 "$lg")"
    out="$(cd /tmp && env "${SITE[@]}" "$L" 2>&1)"; rc=$?
    [[ $rc -ne 0 ]] && has "$out" "REFUSED: missing ${C}/data/handoff/structure/cohort.pcair.pcs.tsv" && ok "the launcher refuses before queuing, too" || bad "launcher: rc=$rc"
    mv "${W}/structure.away" "${C}/data/handoff/structure"
    cp "${C}/data/handoff/structure/cohort.pcrelate.grm.rds" "${W}/grm.keep"
    "${REAL_RSCRIPT}" -e 'm <- readRDS(commandArgs(TRUE)[1]); dimnames(m) <- list(paste0("X", rownames(m)), paste0("X", colnames(m))); saveRDS(m, commandArgs(TRUE)[1])' "${C}/data/handoff/structure/cohort.pcrelate.grm.rds"
    jid="$(pcair)"; lg="${C}/logs/pcair_${jid}.out"
    grep -q "dimnames are not the .fam IIDs" "$lg" && ok "a GRM keyed by other IDs: refused" || bad "grm dimnames: $(tail -3 "$lg")"
    cp "${W}/grm.keep" "${C}/data/handoff/structure/cohort.pcrelate.grm.rds"
    out="$(cd /tmp && env "${SITE[@]}" STRUCTURE_SOURCE=compute "$L" --dry-run 2>&1)"
    has "$out" "STEP 0: Submitting structure (compute: PC-AiR / PC-Relate)" && ok "STRUCTURE_SOURCE=compute runs the comparator instead, only when asked" || bad "compute: $(grep 'STEP 0' <<< "$out")"
    out="$(cd /tmp && env "${SITE[@]}" STRUCTURE_SOURCE=recompute "$L" --dry-run 2>&1)"; rc=$?
    [[ $rc -ne 0 ]] && has "$out" "STRUCTURE_SOURCE=recompute; use bundle or compute" && ok "STRUCTURE_SOURCE outside {bundle, compute}: refused" || bad "bogus source: rc=$rc"
    rm -rf "${C:?}/results" "${C:?}/logs"
else
    echo "  (real R with data.table, or real .rds files in the stub bundle, not available: CP-07 checks skipped)"
fi

section "Run record: the adapter's counts and a stratum skipped below MIN_STRATUM_N"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=os FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 "$L" 2>&1)"
RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"; T="${RD%/}/trace.tsv"
row="$(awk -F'\t' 'NR > 1 && $1 == "gwas" {print $15, $16, $17, $18, $23; exit}' "$T")"
[[ "$row" == "120 900 950 2 ok" ]] && ok "gwas rows: n_null 120, n_tested 900, n_pass_mac 950, n_failed 2, task_status ok (from run_counts.tsv)" || bad "counts row: '$row'"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=os FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 FAKE_ADAPTER_STATUS=skipped_below_min_stratum_n FAKE_N_NULL=31 "$L" 2>&1)"
RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"
grep -q '| gwas | 22 | 22 | 0 | 0 | 22 |' "${RD%/}/summary.md" && ok "summary.md: 22 gwas tasks exit 0, all 22 skipped below MIN_STRATUM_N" || bad "summary: $(grep '| gwas' "${RD%/}/summary.md")"
grep -q 'Stratum EUR: skipped on every chromosome (below MIN_STRATUM_N); not pooled' "${C}"/logs/combine_os_*.out && ok "combine says the stratum was skipped below MIN_STRATUM_N and is not pooled" || bad "combine message: $(cat "${C}"/logs/combine_os_*.out)"

section "CP-11: qc_assert / qc_gate, sample accounting at every join, gates between steps"
QG="${C}/bin/qc_gate.sh"; QD="${W}/qc_unit"; rm -rf "${QD:?}"; mkdir -p "${QD}/assertions.d"
out="$("$QG" --require pcair_pc_merge --dir "$QD" 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "missing assertion pcair_pc_merge" && ok "gate: no pcair_pc_merge row -> non-zero, 'missing assertion'" || bad "missing: rc=$rc $out"
printf 'pcair_pc_merge\tPASS\t8\tevery row has PCs\n' > "${QD}/assertions.d/a.tsv"
out="$("$QG" --require pcair_pc_merge --dir "$QD" 2>&1)"; rc=$?
[[ $rc -eq 0 ]] && ok "gate: the row says PASS -> exit 0" || bad "pass: rc=$rc $out"
printf 'pcair_pc_merge\tFAIL\t8\t2 of 8 rows have no PCs\n' > "${QD}/assertions.d/b.tsv"
out="$("$QG" --require pcair_pc_merge --dir "$QD" 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "assertion pcair_pc_merge failed: 2 of 8 rows have no PCs" && ok "gate: a FAIL row -> non-zero, naming the id and the reason" || bad "fail: rc=$rc $out"
out="$(QC_GATE_OVERRIDE=other,pcair_pc_merge "$QG" --require pcair_pc_merge --dir "$QD" 2>&1)"; rc=$?
[[ $rc -eq 0 ]] && has "$out" "overridden: pcair_pc_merge" && grep -q $'^pcair_pc_merge\tOVERRIDE\tNA\tQC_GATE_OVERRIDE by ' "${QD}"/assertions.d/*.override.*.tsv \
    && ok "QC_GATE_OVERRIDE=<id> proceeds and records an OVERRIDE row" || bad "override: rc=$rc $out"
out="$("$QG" --dir "$QD" 2>&1)"; rc=$?
[[ $rc -eq 2 ]] && has "$out" "--require is empty" && ok "a gate that requires nothing is refused (no skip-all)" || bad "empty require: rc=$rc"
"$QG" --collect --dir "$QD" >/dev/null && [[ "$(head -1 "${QD}/qc_assertions.tsv")" == $'id\tverdict\tn\tmessage' && "$(wc -l < "${QD}/qc_assertions.tsv")" -eq 4 ]] \
    && [[ "$(head -1 "${QD}/sample_accounting.tsv")" == $'step\tstratum\ttrait\tchromosome\tn_in\tn_out\tn_removed\treason' ]] \
    && ok "--collect: qc_assertions.tsv (id verdict n message) and sample_accounting.tsv headers and rows" || bad "collect: $(cat "${QD}/qc_assertions.tsv")"
# Every id the adapter files is on the combine gate's list, and the other way round
filed="$(grep -o 'qc_assert("adapter_[a-z_]*"' "${C}/bin/tractor_genesis_adapter.R" | sed 's/qc_assert("//; s/"//' | sort -u | paste -sd, -)"
req="$(bash -c "source '${C}/bin/slurm_utils.sh'; echo \"\${QC_REQUIRE_GWAS}\"" | tr ',' '\n' | sort -u | paste -sd, -)"
[[ -n "$filed" && "$filed" == "$req" ]] && ok "the adapter files exactly the ids the combine gate requires ($req)" || bad "adapter ids '$filed' vs gate '$req'"
grep -q 'stop("sample ", off_grm\[1\], " missing from GRM "' "${C}/bin/tractor_genesis_adapter.R" \
    && ! grep -q 'Warning: Kinship matrix missing' "${C}/bin/tractor_genesis_adapter.R" && ok "adapter: a joined sample missing from the GRM is a FAIL and a stop, not a warning" || bad "adapter GRM stop"
! grep -q 'all.x = TRUE' "${C}/bin/genesis_pcair_pcrelate.R" && grep -q 'merge_pcs_checked(ph, pcs' "${C}/bin/genesis_pcair_pcrelate.R" \
    && ok "comparator: the all.x merge with a stdout-only count is gone; merge_pcs_checked files pcair_pc_merge" || bad "comparator merge"
if [[ -n "${REAL_RSCRIPT}" ]] && "${REAL_RSCRIPT}" -e 'quit(status = !requireNamespace("data.table", quietly = TRUE))' >/dev/null 2>&1; then
    rm -rf "${QD:?}"
    out="$(GWAS_QC_DIR="$QD" "${REAL_RSCRIPT}" -e '
        suppressPackageStartupMessages(library(data.table)); source(commandArgs(TRUE)[1])
        ph  <- data.table(FID = c("A","B","C"), IID = c("A","B","C"), sample_id = c("A","B","C"), os_time = 1:3)
        pcs <- data.table(sample_id = c("A","B","D"), PC1 = c(.1,.2,.3), PC2 = c(.4,.5,.6))
        m <- merge_pcs_checked(ph, pcs, stratum = "EUR", what = "toy")
        cat("COLS", names(m), "\n")' "${C}/bin/qc_assert.R" 2>&1)"
    "$QG" --collect --dir "$QD" >/dev/null
    has "$out" "COLS FID IID sample_id os_time PC1 PC2" && grep -q $'^pcair_pc_merge\tFAIL\t3\ttoy (EUR): 1 of 3 rows have no PCs (not in the genotype set?), e.g. C$' "${QD}/qc_assertions.tsv" \
        && grep -q $'^pc_merge\tEUR\tNA\tNA\t3\t2\t1\t' "${QD}/sample_accounting.tsv" \
        && ok "merge_pcs_checked: a phenotype row without PCs is a FAIL naming it, with a pc_merge accounting row 3 -> 2" || bad "merge_pcs_checked: $out $(cat "${QD}"/*.tsv)"
    # The toy PLINK set the cluster checks use: a valid SNP-major .bed, T_DUP equal to T001
    "${REAL_RSCRIPT}" "${C}/tests/toy/make_toy_plink.R" "${W}/toy" 10 40 duplicate >/dev/null 2>&1
    out="$("${REAL_RSCRIPT}" -e '
        f <- commandArgs(TRUE)[1]; n <- 11; m <- 40; nb <- ceiling(n / 4)
        b <- as.integer(readBin(paste0(f, ".bed"), "raw", 1e5)); stopifnot(identical(b[1:3], c(108L, 27L, 1L))); b <- b[-(1:3)]
        g <- function(s) { x <- b[(s - 1) * nb + 1:nb]; c(2, NA, 1, 0)[unlist(lapply(x, function(y) (y %/% 4^(0:3)) %% 4))[1:n] + 1] }
        G <- sapply(1:m, g); ph <- read.delim(paste0(f, ".phenotypes.tsv"))
        cat(length(b) == m * nb, !anyNA(G), identical(G[1, ], G[11, ]), nrow(ph), "\n")' "${W}/toy" 2>&1)"
    [[ "$out" == "TRUE TRUE TRUE 14 " ]] && [[ "$(wc -l < "${W}/toy.bim")" -eq 40 && "$(tail -1 "${W}/toy.fam" | cut -d' ' -f2)" == T_DUP ]] \
        && ok "tests/toy/make_toy_plink.R: a valid SNP-major .bed, T_DUP copies T001, 3 phenotype rows not genotyped" || bad "toy plink: $out"
    # End to end through the fake Slurm: step 0 (real ingest) files pcair_pc_merge, the GWAS tasks
    # pass their gate, the stub adapter files its four ids, combine passes its gate
    rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
    out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=os FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 "$L" 2>&1)"
    RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"; RD="${RD%/}"
    ex="$(awk -F'\t' 'NR > 1 {print $1, $14}' "${RD}/trace.tsv" | sort | uniq -c | tr -s ' ' | paste -sd'|' -)"
    [[ "$ex" == " 1 combine 0| 22 gwas 0| 1 pcair 0" ]] && grep -q "QC GATE: passed (pcair_pc_merge,duplicates_removed)" "${C}"/logs/gwas_os_EUR_*_1.out \
        && grep -q "QC GATE: passed (adapter_stratum_check,adapter_sample_join,adapter_grm_cover,adapter_model_n)" "${C}"/logs/combine_os_*.out \
        && ok "launch: GWAS tasks gate on pcair_pc_merge and duplicates_removed, combine on the adapter's four ids; all pass" || bad "e2e gates: $ex"
    QA="${RD}/qc/qc_assertions.tsv"; SA="${RD}/qc/sample_accounting.tsv"
    [[ "$(awk -F'\t' '$1 == "pcair_pc_merge" && $2 == "PASS"' "$QA" | wc -l)" -eq 5 ]] && grep -q $'^structure_ingest\tPASS\t8\t' "$QA" && grep -q $'^duplicates_removed\tPASS\t0\t' "$QA" \
        && ok "run record qc/qc_assertions.tsv: pcair_pc_merge PASS for cohort + 4 strata, duplicates_removed and structure_ingest PASS" || bad "assertions: $(cut -f1,2 "$QA" | sort | uniq -c)"
    grep -q $'^phenotype_to_fam\tcohort\tNA\tNA\t8\t8\t0\t' "$SA" && grep -q $'^stratum\tEUR\tNA\tNA\t8\t3\t5\t' "$SA" && grep -q $'^pc_merge\tEUR\tNA\tNA\t3\t3\t0\t' "$SA" \
        && ok "qc/sample_accounting.tsv: phenotype -> .fam 8 -> 8, EUR stratum 8 -> 3, PC merge 3 -> 3" || bad "accounting: $(cat "$SA")"
    grep -q '^## QC assertions' "${RD}/summary.md" && grep -q '| pcair_pc_merge | 5 | 5 | 0 | 0 | 0 |' "${RD}/summary.md" \
        && ok "summary.md: a QC assertions table" || bad "summary qc: $(sed -n '/## QC/,$p' "${RD}/summary.md")"
    cp "${C}/results/pcair/qc_assertions.tsv" "${W}/struct_qc.keep"
    # An adapter join that fails stops combine; QC_GATE_OVERRIDE=<id> lets it through, recorded
    rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results/run_record" "${C:?}/logs"
    out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=os FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 FAKE_ADAPTER_QC_FAIL=adapter_grm_cover "$L" 2>&1)"
    RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"; RD="${RD%/}"
    [[ "$(awk -F'\t' 'NR > 1 && $1 == "combine" {print $14}' "${RD}/trace.tsv")" == 3 ]] && grep -q "QC GATE: assertion adapter_grm_cover failed: stub adapter" "${C}"/logs/combine_os_*.out \
        && grep -q -- '- FAIL adapter_grm_cover: stub adapter' "${RD}/summary.md" \
        && ok "a FAIL from the adapter: combine stops (exit 3) naming the id; summary.md lists it" || bad "combine gate: $(tail -3 "${C}"/logs/combine_os_*.out)"
    rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results/run_record" "${C:?}/logs"
    out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=os FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 FAKE_ADAPTER_QC_FAIL=adapter_grm_cover QC_GATE_OVERRIDE=adapter_grm_cover "$L" 2>&1)"
    RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"; RD="${RD%/}"
    [[ "$(awk -F'\t' 'NR > 1 && $1 == "combine" {print $14}' "${RD}/trace.tsv")" == 0 ]] && grep -q $'^adapter_grm_cover\tOVERRIDE\t' "${RD}/qc/qc_assertions.tsv" \
        && grep -q $'^QC_GATE_OVERRIDE\tadapter_grm_cover\tenvironment$' "${RD}/params.tsv" \
        && ok "QC_GATE_OVERRIDE=adapter_grm_cover: combine runs; OVERRIDE row in qc_assertions.tsv, the switch in params.tsv" || bad "override e2e"
    # SKIP_PCAIR=true: the reused structure's verdicts are carried in, marked; without them the GWAS gate stops
    rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results/run_record" "${C:?}/logs"
    out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=os FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 SKIP_PCAIR=true "$L" 2>&1)"
    RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"; RD="${RD%/}"
    grep -q $'^pcair_pc_merge\tPASS\t.*\[carried from .*results/pcair/qc_assertions.tsv, SKIP_PCAIR=true\]$' "${RD}/qc/qc_assertions.tsv" \
        && [[ "$(awk -F'\t' 'NR > 1 && $1 == "gwas" && $14 == 0' "${RD}/trace.tsv" | wc -l)" -eq 22 ]] \
        && ok "SKIP_PCAIR=true: the structure verdicts are carried into the run, marked as carried; GWAS runs" || bad "carry: $(head -3 "${RD}/qc/qc_assertions.tsv")"
    rm -f "${C}/results/pcair/qc_assertions.tsv"; rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results/run_record" "${C:?}/logs"
    out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=os FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 SKIP_PCAIR=true "$L" 2>&1)"
    RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"; RD="${RD%/}"
    has "$out" "qc_assertions.tsv not found; the GWAS jobs' QC gate will stop on a missing pcair_pc_merge" \
        && [[ "$(awk -F'\t' 'NR > 1 && $1 == "gwas" && $14 == 3' "${RD}/trace.tsv" | wc -l)" -eq 22 ]] \
        && grep -q "QC GATE: missing assertion pcair_pc_merge" "${C}"/logs/gwas_os_EUR_*_1.err \
        && ok "SKIP_PCAIR=true with no verdicts to carry: warned at launch; every GWAS task stops at its gate (exit 3)" || bad "no carry: $(awk -F'\t' 'NR > 1 {print $1, $14}' "${RD}/trace.tsv" | sort | uniq -c)"
    cp "${W}/struct_qc.keep" "${C}/results/pcair/qc_assertions.tsv"
    rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
else
    echo "  (real R with data.table not available: CP-11 R-side checks skipped)"
fi
out="$(cd /tmp && DRY_RUN=1 bash "$G" EUR os survival 2>&1)"
has "$out" "[DRY-RUN] QC gate: bin/qc_gate.sh --require pcair_pc_merge" && ok "DRY_RUN=1 GWAS task shows its gate" || bad "dry gate: $out"

section "CP-08: covariates exist and are numeric; the list used is written above every results table"
A="${C}/bin/tractor_genesis_adapter.R"
l_chk="$(grep -n 'check_covariates(pheno, covariates_requested' "$A" | cut -d: -f1)"; l_rd="$(grep -n 'tractor_read_full(tractor_files' "$A" | head -1 | cut -d: -f1)"
l_fit="$(grep -n 'nullmod <- fitNullModel(' "$A" | cut -d: -f1)"; l_cox="$(grep -n 'cox_null <- coxph(' "$A" | cut -d: -f1)"
[[ -n "$l_chk" && "$l_chk" -lt "$l_rd" && "$l_chk" -lt "$l_fit" && "$l_chk" -lt "$l_cox" ]] && grep -q 'if (length(covariates) > 0) check_covariates(pheno, covariates, opt\$phenotype)' "$A" \
    && ok "adapter: covariates checked right after the phenotype is read, before Tractor files or any fit (again with global-ancestry columns)" || bad "adapter check order: $l_chk $l_rd $l_fit $l_cox"
[[ "$(grep -c 'fwrite_sumstats(' "$A")" -eq 6 && "$(grep -c '^[^#]*[^_]fwrite(' "$A")" -eq 0 ]] && ok "adapter: all six results tables written through fwrite_sumstats (## covariates=)" || bad "adapter writers: $(grep -n 'fwrite' "$A")"
if [[ -n "${REAL_RSCRIPT}" ]] && "${REAL_RSCRIPT}" -e 'quit(status = !all(sapply(c("data.table", "optparse"), requireNamespace, quietly = TRUE)))' >/dev/null 2>&1; then
    out="$("${REAL_RSCRIPT}" -e '
        suppressPackageStartupMessages(library(data.table)); source(file.path(commandArgs(TRUE)[1], "qc_assert.R"))
        ph <- data.table(sample_id = c("A","B"), age = c(5, 6), sex = c(1L, 2L), PC1 = c(.1, .2), sexc = c("M","F"))
        r <- function(e) tryCatch({ eval(e); "passed" }, error = function(err) conditionMessage(err))
        cat(r(quote(check_covariates(ph, c("age","sex","PC1","PCX"), "pheno.tsv"))), "\n")
        setnames(ph, c("sex", "sexc"), c("sexn", "sex"))
        cat(r(quote(check_covariates(ph, c("age","sex"), "pheno.tsv"))), "\n")
        cat(r(quote(check_covariates(ph, c("age","sexn","PC1"), "pheno.tsv"))), "\n")' "${C}/bin" 2>&1)"
    has "$out" "covariate PCX not found in pheno.tsv" && has "$out" "covariate sex is not numeric (character in pheno.tsv, e.g. 'M')" && has "$out" "passed" \
        && ok "check_covariates: 'covariate PCX not found', 'sex is not numeric', numeric columns pass" || bad "check_covariates: $out"
    TD="${W}/cp08"; rm -rf "${TD:?}"; mkdir -p "$TD"
    out="$("${REAL_RSCRIPT}" -e '
        suppressPackageStartupMessages(library(data.table)); a <- commandArgs(TRUE); source(file.path(a[1], "sumstats_io.R"))
        x <- data.table(CHR = c(1L, 1L), POS = c(100L, 200L), ID = c("1:100", "1:200"), REF = "A", ALT = "G", N = 10L, P_JOINT = c(.5, .01), BETA_EUR = c(.1, -.2), SE_EUR = c(.05, .05))
        for (f in file.path(a[2], c("t.tsv.gz", "t.tsv"))) {
            fwrite_sumstats(x, f, c(covariates = "sex,age_at_diagnosis,PC1"))
            m <- sumstats_meta(f); y <- fread(f)
            cat(basename(f), m, "|", sumstats_meta_value(m, "covariates"), "|", identical(names(y), names(x)), nrow(y), "\n") }' "${C}/bin" "$TD" 2>&1)"
    has "$out" "t.tsv.gz ## covariates=sex,age_at_diagnosis,PC1 | sex,age_at_diagnosis,PC1 | TRUE 2" && has "$out" "t.tsv ## covariates=sex,age_at_diagnosis,PC1 | sex,age_at_diagnosis,PC1 | TRUE 2" \
        && ok "fwrite_sumstats: '## covariates=' line, then the table (gz and plain); fread reads past it" || bad "fwrite_sumstats: $out"
    # meta_analyze_strata.R carries the strata's covariates line: once when they agree, per stratum otherwise
    cp "${TD}/t.tsv.gz" "${TD}/EUR.sumstats.gz"; cp "${TD}/t.tsv.gz" "${TD}/AFR.sumstats.gz"
    "${REAL_RSCRIPT}" "${C}/bin/meta_analyze_strata.R" --input_dir "$TD" --strata EUR,AFR --output "${TD}/P1.gz" --threads 1 >/dev/null 2>&1
    m1="$(zcat "${TD}/P1.gz" 2>/dev/null | grep '^## covariates')"
    (cd "$TD" && zcat t.tsv.gz | sed 's/^## covariates=.*/## covariates=sex,PC1/' | gzip > AFR.sumstats.gz)
    "${REAL_RSCRIPT}" "${C}/bin/meta_analyze_strata.R" --input_dir "$TD" --strata EUR,AFR --output "${TD}/P2.gz" --threads 1 >/dev/null 2>&1
    m2="$(zcat "${TD}/P2.gz" 2>/dev/null | grep '^## covariates' | paste -sd'|' -)"
    [[ "$m1" == "## covariates=sex,age_at_diagnosis,PC1" && "$m2" == "## covariates_EUR=sex,age_at_diagnosis,PC1|## covariates_AFR=sex,PC1" ]] \
        && ok "meta_analyze_strata.R: POOLED carries '## covariates=' (one line if the strata agree, one per stratum if not)" || bad "meta header: '$m1' / '$m2'"
else
    echo "  (real R with data.table and optparse not available: CP-08 R-side checks skipped)"
fi
if python3 -c 'import pandas' >/dev/null 2>&1; then
    TD="${W}/cp08py"; rm -rf "${TD:?}"; mkdir -p "$TD"
    printf '## covariates=sex,PC1\nCHR\tPOS\tID\tP\tBETA\n1\t100\t1:100\t1e-9\t0.1\n1\t200\t1:200\t0.5\t0.2\n' | gzip > "${TD}/in.tsv.gz"
    python3 -I "${C}/bin/filter_gwas_results.py" --input "${TD}/in.tsv.gz" --output "${TD}/out.tsv" --significant "${TD}/sig.tsv" >/dev/null 2>&1
    [[ "$(head -1 "${TD}/sig.tsv" 2>/dev/null)" == $'CHR\tPOS\tID\tP\tBETA' && "$(wc -l < "${TD}/sig.tsv")" -eq 2 ]] \
        && ok "filter_gwas_results.py reads past the '## covariates=' line" || bad "python reader: $(cat "${TD}/sig.tsv" 2>&1)"
fi
grep -q "/^##/ && !h {next} !h {h = 1; next}" "${C}/modules/local/enhanced_annotation/main.nf" && ok "ANNOVAR input (Nextflow) skips '##' lines and the header" || bad "annovar awk"
# The launch: covariates from the bundle unless set; every name checked before queuing
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=os FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 "$L" 2>&1)"
RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"; RD="${RD%/}"
COVS="sex,age_at_diagnosis,cd19_pos_pct,storage_days,PC1,PC2,PC3,PC4,PC5"
grep -q $'^COVARIATES\t'"${COVS}"$'\tbundle$' "${RD}/params.tsv" && grep -q -- "--covariates ${COVS} " "${C}"/logs/gwas_os_EUR_*_1.out \
    && ok "COVARIATES unset: the bundle's covariates.txt + PC1..PC5, recorded in params.tsv (set_by bundle) and passed to the adapter" || bad "default covariates: $(grep COVARIATES "${RD}/params.tsv")"
S="${C}/results/gwas/os/EUR.sumstats.gz"
[[ "$(zcat "$S" | head -1)" == "## covariates=${COVS}" && "$(zcat "$S" | grep -c '^## ')" -eq 1 && "$(zcat "$S" | grep -c '^CHR')" -eq 1 && "$(zcat "$S" | wc -l)" -eq 24 ]] \
    && ok "combine: '## covariates=' kept once, then one header and the 22 chromosomes' rows" || bad "combined: $(zcat "$S" 2>&1 | head -4)"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=os FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 'FAKE_META_CHR7=## covariates=sex' "$L" 2>&1)"
RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"; RD="${RD%/}"
[[ "$(awk -F'\t' 'NR > 1 && $1 == "combine" {print $14}' "${RD}/trace.tsv")" == 1 ]] && grep -q "chr7.tractor_genesis.genomic_order.tsv.gz has other '##' lines than" "${C}"/logs/combine_os_*.out \
    && ok "combine: a chromosome with another covariate list stops the job, naming the file" || bad "meta mismatch: $(tail -3 "${C}"/logs/combine_os_*.out)"
out="$(cd /tmp && env "${SITE[@]}" COVARIATES=age,sex,PC1,PCX "$L" --dry-run 2>&1)"
has "$out" "WARNING: covariate age is not a column of ${C}/data/handoff/phenotypes/phenotypes.tsv" && has "$out" "WARNING: covariate PCX is not a column of ${C}/data/handoff/phenotypes/phenotypes.tsv" \
    && ok "COVARIATES=age,sex,PC1,PCX --dry-run: names age and PCX (not in the phenotype / PC files)" || bad "dry covariates: $(grep -i covariate <<< "$out")"
out="$(cd /tmp && env "${SITE[@]}" COVARIATES=age,sex,PC1,PCX "$L" 2>&1)"; rc=$?
[[ $rc -eq 2 ]] && has "$out" "REFUSED: covariate PCX is not a column of" && ok "... and a launch refuses before queuing (exit 2)" || bad "launch covariates: rc=$rc"
out="$(cd /tmp && env "${SITE[@]}" COVARIATES=sex,PC6 "$L" --dry-run 2>&1)"
has "$out" "covariate PC6 is not a column of ${C}/data/handoff/structure/strata/EUR.pcair.pcs.tsv" && ! has "$out" "PC6 is not a column of ${C}/data/handoff/structure/cohort.pcair.pcs.tsv" \
    && ok "a PC the stratum files lack is named per file (the cohort file has PC6)" || bad "PC6: $(grep PC6 <<< "$out")"
mv "${C}/data/handoff/phenotypes/covariates.txt" "${W}/cov.away"
out="$(cd /tmp && env "${SITE[@]}" "$L" 2>&1)"; rc=$?
[[ $rc -eq 2 ]] && has "$out" "REFUSED: no covariates: COVARIATES is not set and ${C}/data/handoff/phenotypes/covariates.txt is missing" \
    && ok "no COVARIATES and no covariates.txt: refused" || bad "no covariates.txt: rc=$rc $(grep -i covariat <<< "$out")"
mv "${W}/cov.away" "${C}/data/handoff/phenotypes/covariates.txt"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"

section "CP-13: one of each duplicate / MZ pair (kinship > 0.3536) leaves the GWAS phenotype"
check "comparator: duplicates selected from the PC-Relate kinship, call rate from the pruned SNPs, removed before the PC merge" \
    grep -q 'dups <- select_duplicates(kin_mat, call_rate)' "${C}/bin/genesis_pcair_pcrelate.R"
check "comparator: writes <prefix>.removed_duplicates.tsv and files duplicates_removed" \
    grep -q 'qc_assert("duplicates_removed", TRUE' "${C}/bin/genesis_pcair_pcrelate.R"
check "the GWAS gate requires duplicates_removed" grep -q '^QC_REQUIRE_STRUCTURE="pcair_pc_merge,duplicates_removed"$' "${C}/bin/slurm_utils.sh"
if [[ -n "${REAL_RSCRIPT}" ]] && "${REAL_RSCRIPT}" -e 'quit(status = !requireNamespace("data.table", quietly = TRUE))' >/dev/null 2>&1 \
   && ! grep -q "placeholder rds" "${C}/data/handoff/structure/cohort.pcrelate.grm.rds"; then
    out="$("${REAL_RSCRIPT}" -e '
        source(file.path(commandArgs(TRUE)[1], "duplicates.R")); f <- file.path(commandArgs(TRUE)[2], "cr")
        n <- 6; m <- 8; G <- matrix(0L, m, n); G[1:2, 2] <- 1L; G[1:4, 5] <- 1L; pad <- cbind(G, matrix(0L, m, 2))   # 01 = missing
        bytes <- as.vector(sapply(1:m, function(i) c(sum(pad[i, 1:4] * 4^(0:3)), sum(pad[i, 5:8] * 4^(0:3)))))
        con <- file(paste0(f, ".bed"), "wb"); writeBin(as.raw(c(0x6c, 0x1b, 0x01, bytes)), con); close(con)
        ids <- paste0("S", 1:6); write.table(data.frame(ids, ids, 0, 0, 0, -9), paste0(f, ".fam"), quote = FALSE, row.names = FALSE, col.names = FALSE)
        cr <- bed_call_rate(f, ids); cat("CR", cr, "\n")
        k <- diag(0.5, 6); dimnames(k) <- list(ids, ids); cat("CLEAN", nrow(select_duplicates(k)), "\n")
        k["S1", "S2"] <- k["S2", "S1"] <- 0.49; k["S2", "S3"] <- k["S3", "S2"] <- 0.48; k["S5", "S6"] <- k["S6", "S5"] <- 0.5; k["S4", "S6"] <- k["S6", "S4"] <- 0.3
        d <- select_duplicates(k, cr); cat("DUP", paste(d$kept_iid, d$removed_iid, sep = ">"), "\n"); cat("WHY", d$reason[d$removed_iid == "S3"], "\n")' \
        "${C}/bin" "$W" 2>&1)"
    has "$out" "CR 1 0.75 1 1 0.5 1" && has "$out" "CLEAN 0" && has "$out" "DUP S1>S3 S1>S2 S6>S5" && has "$out" "grouped with S1 through S2" \
        && ok "select_duplicates / bed_call_rate: groups joined by pairs > 0.3536 keep the highest call rate (ties: first ID); 0.3 is not a duplicate" || bad "duplicates unit: $out"
    rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"; mkdir -p "${C}/logs"
    jid="$(pcair)"
    D="${C}/results/pcair/cohort.removed_duplicates.tsv"
    [[ "$(cat "$D")" == $'kept_iid\tremoved_iid\tkinship\treason' ]] && grep -q $'^duplicates_removed\tPASS\t0\tshipped kinship checked for pairs > 0.3536: 0 sample(s) removed .*upstream: duplicates_removed=no, 0 listed' "${C}/results/pcair/qc_assertions.tsv" \
        && ok "clean bundle: header-only cohort.removed_duplicates.tsv, duplicates_removed PASS naming the upstream record" || bad "clean: $(cat "$D" "${C}/results/pcair/qc_assertions.tsv" 2>&1)"
    KIN="${C}/data/handoff/structure/cohort.pcrelate.kinship.rds"; cp "$KIN" "${W}/kin.keep"
    "${REAL_RSCRIPT}" -e 'f <- commandArgs(TRUE)[1]; m <- readRDS(f); m["S001", "S002"] <- m["S002", "S001"] <- 0.49; saveRDS(m, f)' "$KIN"
    jid="$(pcair)"
    [[ "$(tail -n +2 "$D" | cut -f1,2)" == $'S001\tS002' ]] && ! grep -q '^S002\b' <(cut -f3 "${C}/results/pcair/EUR.phenotypes.with_pcs.tsv") \
        && ! cut -f3 "${C}/results/pcair/cohort.phenotypes.with_pcs.tsv" | grep -qx S002 && cut -f3 "${C}/results/pcair/EUR.phenotypes.with_pcs.tsv" | grep -qx S001 \
        && grep -q $'^duplicates_removed\tPASS\t1\t' "${C}/results/pcair/qc_assertions.tsv" && grep -q $'^duplicates\tcohort\tNA\tNA\t8\t7\t1\t' "${C}"/results/qc/accounting.d/duplicates.*.tsv \
        && ok "S001-S002 at 0.49: S002 removed (equal call rate, first ID kept), absent from every with_pcs phenotype; PASS; accounting 8 -> 7" || bad "dup: $(cat "$D"; cut -f3 "${C}/results/pcair/EUR.phenotypes.with_pcs.tsv" | paste -sd' ')"
    # call rate decides: S001 missing every call in the stub .bed -> S002 stays
    BED="${C}/data/handoff/genotypes/cohort.hardcall.maf01.bed"; cp "$BED" "${W}/bed.keep"
    { printf '\x6c\x1b\x01'; for v in 1 2 3 4; do printf '\x01\x00'; done; } > "$BED"     # S001 = code 01 (missing) at all 4 sites
    jid="$(pcair)"
    [[ "$(tail -n +2 "$D" | cut -f1,2)" == $'S002\tS001' ]] && grep -q "kept the higher call rate (1.0000 vs 0.0000)" "$D" \
        && ok "the higher call rate stays (S001 has no calls in the .bed: S002 kept)" || bad "call rate tiebreak: $(cat "$D")"
    cp "${W}/bed.keep" "$BED"
    # the bundle's own record contradicted: FAIL (the GWAS gate then stops)
    MAN="${C}/data/handoff/manifest.tsv"; cp "$MAN" "${W}/man.keep"; sed -i 's/^# duplicates_removed\tno$/# duplicates_removed\tyes/' "$MAN"
    jid="$(pcair)"
    grep -q $'^duplicates_removed\tFAIL\t1\t.*the manifest says duplicates_removed=yes, yet 1 pair(s) > 0.3536 remain in the shipped kinship' "${C}/results/pcair/qc_assertions.tsv" \
        && ok "manifest says duplicates_removed=yes but a pair remains: duplicates_removed FAIL" || bad "manifest yes: $(grep duplicates_removed "${C}/results/pcair/qc_assertions.tsv")"
    cp "${W}/man.keep" "$MAN"; cp "${W}/kin.keep" "$KIN"
    UP="${C}/data/handoff/qc/removed_duplicates.tsv"; cp "$UP" "${W}/up.keep"; printf 'S004\tS003\t0.5\tduplicate\n' >> "$UP"
    jid="$(pcair)"
    grep -q $'^duplicates_removed\tFAIL\t0\t.*qc/removed_duplicates.tsv lists 1 removed sample(s) still in the .fam (e.g. S003)' "${C}/results/pcair/qc_assertions.tsv" \
        && ok "upstream lists a removed sample that is still in the .fam: FAIL" || bad "upstream list: $(grep duplicates_removed "${C}/results/pcair/qc_assertions.tsv")"
    cp "${W}/up.keep" "$UP"
    rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
else
    echo "  (real R with data.table, or real .rds files in the stub bundle, not available: CP-13 R-side checks skipped)"
fi

section "CP-14: --p-col with P_JOINT first in the plot and filter scripts; the QQ lambda is the plotted column's"
QQ="${C}/bin/create_qq_plot.R"; MH="${C}/bin/create_manhattan_plot.R"
! grep -q 'lambda_gc <- as.numeric(readLines(args$lambda_gc)' "$QQ" "$MH" && grep -q 'writeLines(format(lambda_gc, digits = 6), paste0(args$output_prefix, ".lambda.txt"))' "$QQ" \
    && ok "QQ: the computed lambda is never overwritten by --lambda-gc; it is written to <prefix>.lambda.txt (Manhattan: same lambda in the title)" || bad "lambda overwrite still there"
[[ "$(grep -c 'fread(cmd = paste("zcat", shQuote(args$sumstats)))' "$QQ" "$MH" | awk -F: '{s += $2} END {print s}')" -eq 2 ]] && ok "plot scripts read .gz through zcat (path quoted)" || bad "zcat"
if [[ -n "${REAL_RSCRIPT}" ]]; then
    out="$("${REAL_RSCRIPT}" -e '
        for (f in commandArgs(TRUE)) {
            for (x in parse(f)) if (is.call(x) && identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name("resolve_p_col"))) eval(x)
            nm <- c("CHROM", "POS", "ID", "P_EUR", "P_JOINT")
            r <- function(...) tryCatch(resolve_p_col(...), error = function(e) paste("ERR", conditionMessage(e)))
            cat(basename(f), r(nm), r(nm, "P_EUR"), r(nm, "p_eur"), r(c("CHR", "P"), NULL), "|", r(nm, "P_AFR"), "|", r(c("CHR", "BETA")), "\n") }' "$QQ" "$MH" 2>&1)"
    want=" P_JOINT P_EUR P_EUR P | ERR --p-col P_AFR is not a column (P columns: P_EUR, P_JOINT) | ERR Could not detect the P-value column"
    has "$out" "create_qq_plot.R${want}" && has "$out" "create_manhattan_plot.R${want}" \
        && ok "resolve_p_col (both plot scripts): P_JOINT by default, --p-col P_EUR, unknown --p-col named with the P columns there are" || bad "resolver: $out"
fi
if python3 -c 'import pandas' >/dev/null 2>&1 && [[ -n "${REAL_RSCRIPT}" ]]; then
    TD="${W}/cp14"; rm -rf "${TD:?}"; mkdir -p "$TD"
    "${REAL_RSCRIPT}" -e 'set.seed(1); n <- 100; d <- data.frame(CHROM = rep(1:2, each = 50), POS = rep(100000L + 10000L * (0:49), 2), ID = paste0("v", 1:n), P_JOINT = runif(n), P_EUR = runif(n))
        d$P_JOINT[3] <- 1e-9; d$P_EUR[7] <- 2e-9; d$P_EUR[8] <- 3e-6
        con <- gzfile(commandArgs(TRUE)[1], "w"); writeLines("## covariates=sex,PC1", con); write.table(d, con, sep = "\t", quote = FALSE, row.names = FALSE); close(con)' "${TD}/s.tsv.gz"
    F="${C}/bin/filter_gwas_results.py"
    a="$(python3 -I "$F" --input "${TD}/s.tsv.gz" --output "${TD}/a.tsv.gz" --significant "${TD}/a.sig.tsv" 2>&1)"
    b="$(python3 -I "$F" --input "${TD}/s.tsv.gz" --output "${TD}/b.tsv.gz" --significant "${TD}/b.sig.tsv" --p-col P_EUR 2>&1)"
    c="$(python3 -I "$F" --input "${TD}/s.tsv.gz" --output "${TD}/c.tsv.gz" --significant "${TD}/c.sig.tsv" --p-col P_AFR 2>&1)"; rc=$?
    has "$a" "P column: P_JOINT" && [[ "$(cut -f3 "${TD}/a.sig.tsv" | paste -sd, -)" == "ID,v3" ]] \
        && has "$b" "P column: P_EUR" && [[ "$(cut -f3 "${TD}/b.sig.tsv" | paste -sd, -)" == "ID,v7" && "$(zcat "${TD}/b.tsv.gz" | cut -f3 | paste -sd, -)" == "ID,v7,v8" ]] \
        && [[ $rc -ne 0 ]] && has "$c" "--p-col P_AFR is not a column (P columns: P_JOINT, P_EUR)" \
        && ok "filter_gwas_results.py on 100 rows (CHROM POS ID P_JOINT P_EUR): P_JOINT by default, --p-col P_EUR, an unknown --p-col refused" || bad "filter: $a / $b / $c"
fi

section "M-05: distance-based loci (500 kb, not LD-based) in filter_gwas_results.py"
if python3 -c 'import pandas' >/dev/null 2>&1; then
    TD="${W}/m05"; rm -rf "${TD:?}"; mkdir -p "$TD"
    printf 'CHR\tPOS\tID\tP_JOINT\tSTRATUM\n1\t1000000\ta\t1e-7\tEUR\n1\t1100000\tb\t1e-6\tEUR\n1\t1200000\tc\t5e-6\tEUR\n1\t3200000\td\t2e-6\tEUR\n1\t5000000\te\t0.5\tEUR\n2\t1000000\tf\t0.3\tEUR\n1\t1050000\tg\t3e-6\tAFR\n' > "${TD}/s.tsv"
    out="$(python3 -I "${C}/bin/filter_gwas_results.py" --input "${TD}/s.tsv" --output "${TD}/o.tsv.gz" --significant "${TD}/o.significant.tsv" 2>&1)"
    Lc="${TD}/o.loci.tsv"
    head -1 "$Lc" | grep -q "^## method=distance-based, not LD-based: a locus is every suggestive variant within 500 kb" && grep -qx '## window_kb=500' "$Lc" \
        && [[ "$(grep -v '^##' "$Lc" | cut -f1-4,8,9 | paste -sd'|' -)" == $'STRATUM\tLOCUS\tCHR\tINDEX_ID\tN_VARIANTS\tP_JOINT|AFR\tAFR:1\t1\tg\t1\t3e-06|EUR\tEUR:1\t1\ta\t3\t1e-07|EUR\tEUR:2\t1\td\t1\t2e-06' ]] \
        && ok "three EUR variants within 200 kb + one 2 Mb away: two EUR loci (index a, d), the AFR variant its own locus; header states the window and 'distance-based'" || bad "loci: $(cat "$Lc")"
    [[ "$(zcat "${TD}/o.tsv.gz" | awk -F'\t' 'NR > 1 {print $3 ":" $6 ":" $8}' | sort | paste -sd' ' -)" == "a:EUR:1:True b:EUR:1:False c:EUR:1:False d:EUR:2:True g:AFR:1:True" ]] \
        && ok "the variant list carries LOCUS, LOCUS_INDEX and IS_INDEX" || bad "variant list: $(zcat "${TD}/o.tsv.gz")"
    out="$(python3 -I "${C}/bin/filter_gwas_results.py" --input "${TD}/s.tsv" --output "${TD}/w.tsv.gz" --significant "${TD}/w.significant.tsv" --window-kb 50 2>&1)"
    [[ "$(grep -vc '^##' "${TD}/w.loci.tsv")" -eq 6 ]] && grep -qx '## window_kb=50' "${TD}/w.loci.tsv" && ok "--window-kb 50: each EUR variant its own locus" || bad "window 50: $(cat "${TD}/w.loci.tsv")"
fi
grep -q 'path("${prefix}.loci.tsv"),        emit: loci' "${C}/modules/local/gwas_filter/main.nf" && ok "Nextflow GWAS_FILTER emits the locus table" || bad "module loci"

section "CP-22: STEP 2 reads the bundle's strata, checks them disjoint and complete, writes POOLED.ivw_descriptive"
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR,AAC "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -eq 2 ]] && has "$out" "REFUSED: stratum AAC (STRATA_LIST) is not in ${C}/data/handoff/ancestry/strata.tsv (it has: AFR,AMR,EUR,MIX)" \
    && ok "a stratum the bundle's strata.tsv does not have: refused" || bad "unknown stratum: rc=$rc $(tail -2 <<< "$out")"
if [[ -n "${REAL_RSCRIPT}" ]] && "${REAL_RSCRIPT}" -e 'quit(status = !all(sapply(c("data.table", "optparse"), requireNamespace, quietly = TRUE)))' >/dev/null 2>&1 \
   && ! grep -q "placeholder rds" "${C}/data/handoff/structure/cohort.pcrelate.grm.rds"; then
    e2e() { rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results/run_record" "${C:?}/results/gwas" "${C:?}/logs"
            out="$(cd /tmp && env "${SITE[@]}" TRAITS_LIST=os FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 "$@" "$L" 2>&1)"
            RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"; RD="${RD%/}"
            CX="$(awk -F'\t' 'NR > 1 && $1 == "combine" {print $14}' "${RD}/trace.tsv")"; CL="$(cat "${C}"/logs/combine_os_*.out 2>/dev/null)"; }
    rm -rf "${C:?}/results"
    e2e
    grep -q $'^STRATA_LIST\tAFR,AMR,EUR,MIX\tbundle$' "${RD}/params.tsv" && ok "STRATA_LIST unset: the bundle's strata (AFR,AMR,EUR,MIX), set_by bundle in params.tsv" || bad "default strata: $(grep STRATA_LIST "${RD}/params.tsv")"
    P="${C}/results/gwas/os/POOLED.ivw_descriptive.sumstats.gz"
    [[ "$CX" == 0 ]] && has "$CL" "Strata disjoint: AFR,AMR,EUR,MIX (8 samples)" && [[ -s "$P" && ! -e "${C}/results/gwas/os/POOLED.sumstats.gz" ]] \
        && [[ "$(zcat "$P" | grep '^## \(strata\|strata_skipped\|backgrounds\|covariates\)=' | paste -sd'|' -)" == "## strata=AFR,AMR,EUR,MIX|## strata_skipped=none|## backgrounds=EUR|## covariates=sex,age_at_diagnosis,cd19_pos_pct,storage_days,PC1,PC2,PC3,PC4,PC5" ]] \
        && ok "four strata: disjoint, all 22 chromosomes each, pooled into POOLED.ivw_descriptive.sumstats.gz with ## strata= / strata_skipped= / backgrounds= / covariates=" || bad "pooled: CX=$CX $(zcat "$P" 2>&1 | head -5) $CL"
    e2e STRATA_LIST=EUR,AFR,AMR FAKE_SKIP_STRATUM=AFR
    [[ "$CX" == 0 ]] && has "$CL" "Stratum AFR: skipped on every chromosome (below MIN_STRATUM_N); not pooled" && zcat "$P" | grep -qx '## strata=EUR,AMR' && zcat "$P" | grep -qx '## strata_skipped=AFR' \
        && ok "a stratum skipped below MIN_STRATUM_N: left out, and named in ## strata_skipped=" || bad "skipped stratum: CX=$CX $(zcat "$P" 2>&1 | head -3)"
    e2e STRATA_LIST=EUR,AFR,AMR FAKE_NO_OUTPUT=AMR
    [[ "$CX" == 1 ]] && has "$CL" "ERROR: stratum AMR: no results for 22 chromosome(s): chr1 chr2" \
        && ok "one stratum's chromosome files gone: STEP 2 exits non-zero naming the stratum (no silent continue)" || bad "missing stratum: CX=$CX $CL"
    e2e STRATA_LIST=EUR,AFR FAKE_SKIP_STRATUM=AFR
    [[ "$CX" == 0 ]] && has "$CL" "POOLED not written: 1 stratum(s) analysed (EUR); pooling needs two" && [[ ! -e "$P" ]] \
        && ok "one stratum left: no POOLED file, and the log says why" || bad "single stratum: CX=$CX $CL"
    printf 'S001\n' >> "${C}/results/strata/AFR.samples.txt"
    e2e STRATA_LIST=EUR,AFR SKIP_PCAIR=true
    [[ "$CX" == 1 ]] && has "$CL" "ERROR: strata not disjoint: 1 sample(s); e.g. EUR and AFR share S001" \
        && ok "an overlapping samples.txt pair halts STEP 2 with 'strata not disjoint'" || bad "overlap: CX=$CX $CL"
    rm -rf "${C:?}/results/strata"; mkdir -p "${C}/logs"
    jid="$(pcair STRUCTURE_SOURCE=compute)"
    [[ "$(cat "${C}/results/strata/EUR.samples.txt" 2>/dev/null | paste -sd, -)" == "S001,S002,S003" && -s "${C}/results/strata/MIX.samples.txt" ]] \
        && ok "STRUCTURE_SOURCE=compute: results/strata/<S>.samples.txt written from strata.tsv too" || bad "compute strata: $(ls "${C}/results/strata" 2>&1)"
    rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
else
    echo "  (real R with data.table and optparse, or real .rds files in the stub bundle, not available: CP-22 end-to-end checks skipped)"
fi

section "CP-20: MR-MEGA axes = min(K-2, MRMEGA_PCS) from the strata with sumstats"
M="${C}/slurm/submit_mrmega.sh"; GD="${C}/results/gwas/os"
rm -rf "${C:?}/results"; mkdir -p "$GD"
for s in EUR AFR AMR; do echo x | gzip > "${GD}/${s}.sumstats.gz"; done
out="$(cd /tmp && DRY_RUN=1 bash "$M" os 2>&1)"; rc=$?
[[ $rc -eq 0 ]] && has "$out" "MR-MEGA os: K=3 pc_axes=1 (min(K-2, MRMEGA_PCS=3))" && has "$out" "strata with sumstats: AFR,AMR,EUR; without: MIX" \
    && ok "three stratum files: the dry run logs K=3 pc_axes=1" || bad "K=3: $out"
echo x | gzip > "${GD}/MIX.sumstats.gz"; out4="$(cd /tmp && DRY_RUN=1 bash "$M" os 2>&1)"; out4b="$(cd /tmp && DRY_RUN=1 MRMEGA_PCS=1 bash "$M" os 2>&1)"
has "$out4" "K=4 pc_axes=2" && has "$out4b" "K=4 pc_axes=1 (min(K-2, MRMEGA_PCS=1))" && ok "K=4: pc_axes 2, capped at MRMEGA_PCS=1 when set" || bad "K=4: $out4 / $out4b"
rm -f "${GD}/AMR.sumstats.gz" "${GD}/MIX.sumstats.gz"
out="$(cd /tmp && DRY_RUN=1 bash "$M" os 2>&1)"
has "$out" "MR-MEGA os: K=2 pc_axes=0: fixed-effect meta-analysis, not MR-MEGA" && ok "two stratum files: pc_axes=0 and the fixed-effect caveat" || bad "K=2: $out"
rm -f "${GD}/AFR.sumstats.gz"
out="$(cd /tmp && DRY_RUN=1 bash "$M" os 2>&1)"; rc=$?
[[ $rc -eq 0 ]] && has "$out" "MR-MEGA os: excluded: K=1 < pc_axes + 2 = 2" && ok "one stratum file: the trait is logged as excluded" || bad "K=1: rc=$rc $out"
out="$(cd /tmp && DRY_RUN=1 STRATA_LIST=EUR,AFR bash "$M" os 2>&1)"
has "$out" "strata with sumstats: EUR; without: AFR" && ok "STRATA_LIST limits the strata counted" || bad "strata list: $out"
rm -rf "${C:?}/results"
out="$(cd /tmp && env "${SITE[@]}" RUN_MRMEGA=true "$L" --dry-run 2>&1)"
mm="$(grep '\[DRY-RUN\] MR-MEGA' <<< "$out")"
[[ "$(wc -l <<< "$mm")" -eq 3 ]] && grep -q 'MR-MEGA mrd: sbatch .*--dependency=afterok:dry14 .*submit_mrmega.sh mrd' <<< "$mm" \
    && ok "RUN_MRMEGA=true: one MR-MEGA job per trait, each after its own combine job" || bad "mrmega dry: $mm"
out="$(cd /tmp && env "${SITE[@]}" MRMEGA_PCS=0 "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -eq 2 ]] && has "$out" "REFUSED: MRMEGA_PCS=0; use a whole number >= 1" && ok "MRMEGA_PCS outside 1, 2, ...: refused" || bad "MRMEGA_PCS=0: rc=$rc"

section "B-05: switches checked against their domains at launch; RUN_LABEL; PC_SCALE; MIX_RETAIN"
out="$(cd /tmp && env "${SITE[@]}" MAC_MIN=abc "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -eq 2 ]] && has "$out" "REFUSED: MAC_MIN=abc; use a whole number >= 1" && ! grep -q '\[DRY-RUN\]' <<< "$out" \
    && ok "MAC_MIN=abc --dry-run: refused before any sbatch" || bad "MAC_MIN=abc: rc=$rc $(tail -2 <<< "$out")"
out="$(cd /tmp && env "${SITE[@]}" PC_SCALE=global MIX_RETAIN=no KIN_THRESH=0 GXG_PRUNE_MODE=all "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -eq 2 ]] && has "$out" "REFUSED: PC_SCALE=global; use stratum or cohort" && has "$out" "REFUSED: MIX_RETAIN=no; use true or false" \
    && has "$out" "REFUSED: KIN_THRESH=0; use a number > 0" && has "$out" "REFUSED: GXG_PRUNE_MODE=all; use conditional, variant or pair" \
    && ok "every switch outside its domain is named, all at once" || bad "domains: rc=$rc $out"
pa="$(cd /tmp && env "${SITE[@]}" MAC_MIN=20 RUN_LABEL=a "$L" --dry-run 2>&1)"; pb="$(cd /tmp && env "${SITE[@]}" MAC_MIN=30 RUN_LABEL=b "$L" --dry-run 2>&1)"
ta="$(sed -n '/^params.tsv (dry run: not written):$/,/^end of params.tsv$/p' <<< "$pa")"; tb="$(sed -n '/^params.tsv (dry run: not written):$/,/^end of params.tsv$/p' <<< "$pb")"
d="$(diff <(echo "$ta") <(echo "$tb") | grep '^[<>]')"
has "$pa" "Results:   ${C}/results/a (RUN_LABEL=a)" && has "$pb" "Results:   ${C}/results/b (RUN_LABEL=b)" && [[ "$(wc -l <<< "$d")" -eq 2 ]] \
    && grep -q $'^<   MAC_MIN\t20\tenvironment$' <<< "$d" && grep -q $'^>   MAC_MIN\t30\tenvironment$' <<< "$d" \
    && ok "MAC_MIN=20 RUN_LABEL=a vs MAC_MIN=30 RUN_LABEL=b: results/a and results/b; params.tsv differs in exactly one row" || bad "labels: $d"
out="$(cd /tmp && env "${SITE[@]}" RUN_LABEL=gwas "$L" --dry-run 2>&1)"; rc=$?; out2="$(cd /tmp && env "${SITE[@]}" "RUN_LABEL=a b" "$L" --dry-run 2>&1)"; rc2=$?
[[ $rc -eq 2 && $rc2 -eq 2 ]] && has "$out" "RUN_LABEL=gwas is the name of a results folder" && has "$out2" "REFUSED: RUN_LABEL=a b; use letters, digits" \
    && ok "RUN_LABEL: one plain folder name, not a results folder's name" || bad "bad labels: $rc $rc2"
out="$(cd /tmp && env "${SITE[@]}" MIX_RETAIN=false "$L" --dry-run 2>&1)"
has "$out" "Strata:    AFR AMR EUR" && ! grep -q 'GWAS .* / MIX' <<< "$out" && ok "MIX_RETAIN=false: MIX left out of the default strata" || bad "MIX out: $(grep Strata <<< "$out")"
out="$(cd /tmp && env "${SITE[@]}" MIX_RETAIN=false STRATA_LIST=EUR,MIX "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -eq 2 ]] && has "$out" "REFUSED: STRATA_LIST names MIX, but MIX_RETAIN=false" && ok "... and refused when STRATA_LIST names it" || bad "MIX listed: rc=$rc"
if [[ -n "${REAL_RSCRIPT}" ]] && ! grep -q "placeholder rds" "${C}/data/handoff/structure/cohort.pcrelate.grm.rds"; then
    rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
    out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR,AFR TRAITS_LIST=os RUN_LABEL=mac30 MAC_MIN=30 FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 "$L" 2>&1)"
    RD="$(ls -d "${C}"/results/mac30/run_record/*/ 2>/dev/null | head -1)"; RD="${RD%/}"
    [[ -s "${C}/results/mac30/gwas/os/POOLED.ivw_descriptive.sumstats.gz" && -s "${C}/results/mac30/gwas/os/EUR/chr1.run_counts.tsv" && -s "${C}/results/pcair/EUR.phenotypes.with_pcs.tsv" ]] \
        && [[ ! -e "${C}/results/gwas" && ! -e "${C}/results/run_record" ]] && grep -q $'^run_label\tmac30$' "${RD}/run_info.tsv" && grep -q -- '--mac_min 30 ' "${C}"/logs/gwas_os_EUR_*_1.out \
        && ok "RUN_LABEL=mac30: GWAS, pooled file and run record under results/mac30/; the structure stays shared in results/pcair" || bad "label run: $(ls "${C}/results" "${C}/results/mac30" 2>&1)"
    out="$(cd /tmp && DRY_RUN=1 PC_SCALE=cohort bash "$G" EUR os survival 2>&1)"
    has "$out" "--phenotype ${C}/results/pcair/cohort.phenotypes.with_pcs.tsv" && ok "PC_SCALE=cohort: a stratum's GWAS reads the cohort PCs" || bad "PC_SCALE=cohort: $(grep phenotype <<< "$out")"
    out="$(cd /tmp && DRY_RUN=1 bash "$G" EUR os survival 2>&1)"
    has "$out" "--phenotype ${C}/results/pcair/EUR.phenotypes.with_pcs.tsv" && ok "PC_SCALE=stratum (default): its within-stratum PCs" || bad "PC_SCALE default: $(grep phenotype <<< "$out")"
    rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
fi

section "CP-17 correctness: the adapter itself on toy Tractor files (survival; runs here without GENESIS)"
if [[ -n "${REAL_RSCRIPT}" ]] && "${REAL_RSCRIPT}" -e 'quit(status = !all(sapply(c("data.table", "optparse", "survival", "Matrix"), requireNamespace, quietly = TRUE)))' >/dev/null 2>&1; then
    out="$("${REAL_RSCRIPT}" "${C}/tests/test_tractor_model.R" 2>&1)"; rc=$?
    [[ $rc -eq 0 ]] && ok "tests/test_tractor_model.R: $(tail -1 <<< "$out") (minor-allele MAC, stratum gate, Cox score test == coxph's)" || bad "test_tractor_model.R: $out"
    TD="${W}/cp17"; rm -rf "${TD:?}"; mkdir -p "$TD"
    "${REAL_RSCRIPT}" "${C}/tests/toy/make_toy_tractor.R" "${TD}/toy" 120 60 AFR,EUR na_rows=5 >/dev/null
    ad() { (cd "$TD" && GWAS_QC_DIR="${TD}/qc" "${REAL_RSCRIPT}" "${C}/bin/tractor_genesis_adapter.R" --tractor_prefix toy --phenotype toy.phenotypes.tsv \
            --trait os --model survival --time_col os_time --event_col os_event --ancestries AFR,EUR --ref_ancestry EUR \
            --covariates sex,age_at_diagnosis --stratum TOY --min_stratum_n 50 --chromosome 22 "$@" 2>&1); }
    out="$(ad --sample_include toy.samples.txt --output_prefix out/chr22)"; rc=$?
    O="${TD}/out/chr22.tractor_genesis.genomic_order.tsv.gz"
    [[ $rc -eq 0 ]] && has "$out" "Complete cases: 115 of 120 joined samples; 5 dropped for a missing os_time / os_event / sex / age_at_diagnosis value" \
        && [[ "$(zcat "$O" | grep -v '^##' | head -1 | cut -f1-8)" == $'STRATUM\tCHR\tPOS\tID\tREF\tALT\tN\tMAC_STRATUM' ]] \
        && [[ "$(zcat "$O" | grep -v '^##' | awk -F'\t' 'NR > 1 && $7 != 115' | wc -l)" -eq 0 ]] \
        && grep -qx $'n_null\t115' "${TD}/out/chr22.run_counts.tsv" && grep -qx $'n_failed\t0' "${TD}/out/chr22.run_counts.tsv" \
        && ok "5 rows with NA os_time: reported and dropped before the fit; N 115 on every row, no failed test; STRATUM is the first column" || bad "toy run: rc=$rc $(tail -5 <<< "$out")"
    [[ "$(zcat "$O" | grep -v '^##' | awk -F'\t' 'NR > 1 && $11 == "NA"' | wc -l)" -eq 0 ]] && ok "no all-NA rows: every toy variant has a P_JOINT" || bad "NA rows: $(zcat "$O" | awk -F'\t' '$11 == "NA"' | head -2)"
    head -100 "${TD}/toy.samples.txt" > "${TD}/keep100.txt"
    out="$(ad --sample_include keep100.txt --output_prefix out100/chr22)"
    has "$out" "samples sample_include: 120 -> 100" && grep -qx $'n_null\t95' "${TD}/out100/chr22.run_counts.tsv" \
        && ok "--sample_include (the stratum's samples.txt): 120 -> 100 joined, 95 complete cases" || bad "sample_include: $(grep -i 'sample_include\|n_null' <<< "$out")"
    "${REAL_RSCRIPT}" -e 'ids <- readLines(commandArgs(TRUE)[1]); ids <- ids[ids != "T0007"]; m <- diag(0.5, length(ids)); dimnames(m) <- list(ids, ids); saveRDS(m, commandArgs(TRUE)[2])' "${TD}/toy.samples.txt" "${TD}/kin.rds"
    out="$(ad --kinship kin.rds --output_prefix outk/chr22)"; rc=$?
    [[ $rc -ne 0 ]] && has "$out" "sample T0007 missing from GRM kin.rds (1 joined sample(s) are not in it: T0007)" && ! has "$out" "subscript out of bounds" && ! has "$out" "Fitting null model" \
        && ok "a sample absent from the GRM: 'sample T0007 missing from GRM', before any fit (no subscript error)" || bad "grm missing: rc=$rc $(tail -3 <<< "$out")"
    grep -q 'CMD+=(--sample_include "${SAMPLE_INCLUDE}")' "${C}/slurm/submit_gwas_array.sh" && ok "submit_gwas_array.sh passes results/strata/<S>.samples.txt as --sample_include" || bad "gwas array sample_include"
else
    echo "  (real R with data.table, optparse, survival and Matrix not available: CP-17 adapter checks skipped)"
fi

section "CP-17 performance: chunks in parallel, one null per stratum x trait, no Cox null without the GRM"
if [[ -n "${REAL_RSCRIPT}" ]] && "${REAL_RSCRIPT}" -e 'quit(status = !all(sapply(c("data.table", "optparse", "survival", "Matrix"), requireNamespace, quietly = TRUE)))' >/dev/null 2>&1; then
    TD="${W}/cp17p"; rm -rf "${TD:?}"; mkdir -p "$TD"
    "${REAL_RSCRIPT}" "${C}/tests/toy/make_toy_tractor.R" "${TD}/toy" 120 60 AFR,EUR >/dev/null
    ad() { (cd "$TD" && GWAS_QC_DIR="${TD}/qc" "${REAL_RSCRIPT}" "${C}/bin/tractor_genesis_adapter.R" --tractor_prefix toy --phenotype toy.phenotypes.tsv \
            --trait os --model survival --time_col os_time --event_col os_event --ancestries AFR,EUR --ref_ancestry EUR \
            --covariates sex,age_at_diagnosis --stratum TOY --min_stratum_n 50 --chromosome 22 --chunk_size 7 "$@" 2>&1); }
    o1="$(ad --threads 1 --output_prefix t1/chr22)"; o4="$(ad --threads 4 --output_prefix t4/chr22)"
    cmp="$("${REAL_RSCRIPT}" -e 'suppressPackageStartupMessages(library(data.table)); a <- fread(commandArgs(TRUE)[1]); b <- fread(commandArgs(TRUE)[2])
            num <- names(a)[sapply(a, is.numeric)]; d <- max(abs(as.matrix(a[, ..num]) - as.matrix(b[, ..num])), na.rm = TRUE)
            cat(identical(dim(a), dim(b)) && identical(is.na(a), is.na(b)) && identical(a$ID, b$ID) && d <= 1e-10, "\n")' \
            "${TD}/t1/chr22.tractor_genesis.genomic_order.tsv.gz" "${TD}/t4/chr22.tractor_genesis.genomic_order.tsv.gz" 2>&1)"
    has "$o4" "Chunks: 9 | threads: 4" && [[ "$cmp" == "TRUE " ]] && ok "--threads 4 (9 chunks of 7, forked) gives the --threads 1 results to 1e-10" || bad "threads: $cmp $(grep -i 'chunk\|error' <<< "$o4")"
    oa="$(ad --null_model_rds cache/null_model.rds --output_prefix c/chr22)"; ob="$(ad --null_model_rds cache/null_model.rds --output_prefix c/chr21)"
    has "$oa" "null saved to cache/null_model.rds" && has "$ob" "null loaded from cache/null_model.rds" && ! has "$ob" "Fitted Cox" \
        && ok "--null_model_rds: the first chromosome fits and saves the null, the second logs 'null loaded'" || bad "null cache: $(grep null <<< "$ob")"
    oc="$(ad --null_model_rds cache/null_model.rds --covariates sex --output_prefix c/chr20)"
    has "$oc" "was fitted on other samples or settings: refitting" && ok "... a cache fitted with other covariates is not reused" || bad "stale cache: $(grep null <<< "$oc")"
    "${REAL_RSCRIPT}" -e 'ids <- readLines(commandArgs(TRUE)[1]); m <- diag(0.5, length(ids)); dimnames(m) <- list(ids, ids); saveRDS(m, commandArgs(TRUE)[2])' "${TD}/toy.samples.txt" "${TD}/kin.rds"
    if ! "${REAL_RSCRIPT}" -e 'quit(status = !requireNamespace("coxme", quietly = TRUE))' >/dev/null 2>&1; then
        out="$(ad --kinship kin.rds --output_prefix k/chr22)"; rc=$?
        [[ $rc -ne 0 ]] && has "$out" "coxme required for a Cox null with a GRM" && ok "no coxme, with a GRM: exits non-zero, 'coxme required'" || bad "coxme stop: rc=$rc"
        out="$(ad --kinship kin.rds --allow_unrelated_cox --output_prefix k2/chr22)"; rc=$?
        [[ $rc -eq 0 ]] && has "$out" "WARNING: --allow_unrelated_cox and no coxme: the Cox null is fitted WITHOUT the GRM" && ok "... and runs, saying so, with --allow_unrelated_cox" || bad "allow_unrelated_cox: rc=$rc"
    fi
    out="$(cd /tmp && DRY_RUN=1 bash "$G" EUR os survival 2>&1)"
    has "$out" "--null_model_rds ${C}/results/gwas/os/EUR/null_model.rds" && ok "submit_gwas_array.sh: one null cache per stratum x trait (results/gwas/<trait>/<S>/null_model.rds)" || bad "gwas null cache arg"
else
    echo "  (real R with data.table, optparse, survival and Matrix not available: CP-17 performance checks skipped)"
fi

section "B-01: a one-ancestry stratum (k = 1) through the adapter, no hapcount file"
if [[ -n "${REAL_RSCRIPT}" ]] && "${REAL_RSCRIPT}" -e 'quit(status = !all(sapply(c("data.table", "optparse", "survival", "Matrix"), requireNamespace, quietly = TRUE)))' >/dev/null 2>&1; then
    TD="${W}/b01"; rm -rf "${TD:?}"; mkdir -p "$TD"
    "${REAL_RSCRIPT}" "${C}/tests/toy/make_toy_tractor.R" "${TD}/eur" 120 40 EUR >/dev/null; rm -f "${TD}/eur.hapcount.EUR.txt.gz"
    out="$(cd "$TD" && GWAS_QC_DIR="${TD}/qc" "${REAL_RSCRIPT}" "${C}/bin/tractor_genesis_adapter.R" --tractor_prefix eur --phenotype eur.phenotypes.tsv \
            --trait os --model survival --time_col os_time --event_col os_event --ancestries EUR --ref_ancestry EUR \
            --covariates sex,age_at_diagnosis --stratum EUR --min_stratum_n 50 --chromosome 22 --output_prefix out/chr22 2>&1)"; rc=$?
    chk="$("${REAL_RSCRIPT}" -e 'suppressPackageStartupMessages(library(data.table)); a <- fread(commandArgs(TRUE)[1])
            cat(identical(unique(a$STRATUM), "EUR"), names(a)[1], length(setdiff(grep("^BETA_", names(a), value = TRUE), "BETA_TOTAL")), identical(a$P_JOINT, a$P_EUR), unique(a$DF_JOINT), "\n")' \
            "${TD}/out/chr22.tractor_genesis.genomic_order.tsv.gz" 2>&1)"
    [[ $rc -eq 0 && "$chk" == "TRUE STRATUM 1 TRUE 1 " ]] && has "$out" "EUR: (one ancestry: hapcount 2 for every sample), eur.dosage.EUR.txt.gz" && ! grep -qi "hapcount file\|no hapcount" <<< "$out" \
        && ok "EUR-only toy stratum: STRATUM=EUR first, one BETA_EUR/SE_EUR/P_EUR set, P_JOINT identical to P_EUR (1 df), no hapcount file needed" || bad "k=1: rc=$rc chk=$chk $(tail -3 <<< "$out")"
fi

section "M-02: the Tractor matrices streamed chunk by chunk (default), the full read kept for comparison"
if [[ -n "${REAL_RSCRIPT}" ]] && "${REAL_RSCRIPT}" -e 'quit(status = !all(sapply(c("data.table", "optparse", "survival", "Matrix"), requireNamespace, quietly = TRUE)))' >/dev/null 2>&1; then
    TD="${W}/m02"; rm -rf "${TD:?}"; mkdir -p "$TD"
    "${REAL_RSCRIPT}" "${C}/tests/toy/make_toy_tractor.R" "${TD}/toy" 120 300 AFR,EUR >/dev/null
    ad() { (cd "$TD" && GWAS_QC_DIR="${TD}/qc" "${REAL_RSCRIPT}" "${C}/bin/tractor_genesis_adapter.R" --tractor_prefix toy --phenotype toy.phenotypes.tsv \
            --trait os --model survival --time_col os_time --event_col os_event --ancestries AFR,EUR --ref_ancestry EUR \
            --covariates sex,age_at_diagnosis --stratum TOY --min_stratum_n 50 --chromosome 22 --chunk_size 64 "$@" 2>&1); }
    os="$(ad --read_mode stream --threads 2 --output_prefix s/chr22)"; of="$(ad --read_mode full --output_prefix f/chr22)"
    same="$("${REAL_RSCRIPT}" -e 'suppressPackageStartupMessages(library(data.table)); a <- fread(commandArgs(TRUE)[1]); b <- fread(commandArgs(TRUE)[2]); cat(nrow(a), isTRUE(all.equal(a, b, tolerance = 0)), "\n")' \
            "${TD}/s/chr22.tractor_genesis.genomic_order.tsv.gz" "${TD}/f/chr22.tractor_genesis.genomic_order.tsv.gz" 2>&1)"
    [[ "$same" == "300 TRUE " ]] && has "$os" "read mode: stream" && has "$os" "Chunks: 5 | threads: 2" && has "$of" "Reading the Tractor files whole (--read_mode full)" \
        && ok "300 variants in chunks of 64: --read_mode stream (default, 2 threads) gives exactly the --read_mode full results" || bad "stream vs full: $same"
    has "$os" "Peak memory of this process:" && ok "the adapter reports its peak memory (VmHWM) at the end" || bad "peak memory line"
    zcat "${TD}/toy.dosage.EUR.txt.gz" | sed '101d' | gzip > "${TD}/d.tmp" && mv "${TD}/d.tmp" "${TD}/toy.dosage.EUR.txt.gz"
    out="$(ad --output_prefix x/chr22)"; rc=$?
    [[ $rc -ne 0 ]] && has "$out" "Tractor files out of step in variants 65-128: toy.dosage.EUR.txt.gz (row 36: 22:201000:22:201000:A:G vs 22:200000:22:200000:A:G)" \
        && ok "a dose file missing a row: stops, naming the file and the first variant out of step (no silent misalignment)" || bad "out of step: rc=$rc $(tail -2 <<< "$out")"
else
    echo "  (real R with data.table, optparse, survival and Matrix not available: M-02 checks skipped)"
fi

section "CP-19: the total-dose 1-df test and EAF in every table; the MR-MEGA input; STEP 2b runs"
if [[ -n "${REAL_RSCRIPT}" ]] && "${REAL_RSCRIPT}" -e 'quit(status = !all(sapply(c("data.table", "optparse", "survival", "Matrix"), requireNamespace, quietly = TRUE)))' >/dev/null 2>&1 \
   && python3 -c 'import pandas' >/dev/null 2>&1; then
    TD="${W}/cp19"; rm -rf "${TD:?}"; mkdir -p "$TD"
    "${REAL_RSCRIPT}" "${C}/tests/toy/make_toy_tractor.R" "${TD}/a" 120 100 AFR,EUR >/dev/null
    "${REAL_RSCRIPT}" "${C}/tests/toy/make_toy_tractor.R" "${TD}/e" 120 100 EUR id_prefix=T >/dev/null; rm -f "${TD}/e.hapcount.EUR.txt.gz"
    ad() { local pre=$1 anc=$2 st=$3; shift 3; (cd "$TD" && GWAS_QC_DIR="${TD}/qc" "${REAL_RSCRIPT}" "${C}/bin/tractor_genesis_adapter.R" --tractor_prefix "$pre" \
            --phenotype "${pre}.phenotypes.tsv" --trait os --model survival --time_col os_time --event_col os_event --ancestries "$anc" --ref_ancestry EUR \
            --covariates sex,age_at_diagnosis --stratum "$st" --min_stratum_n 50 --chromosome 22 --output_prefix "$st/chr22" "$@" 2>&1); }
    ad a AFR,EUR AAA >/dev/null; ad e EUR EUR >/dev/null
    chk="$("${REAL_RSCRIPT}" -e 'suppressPackageStartupMessages(library(data.table)); a <- fread(commandArgs(TRUE)[1]); e <- fread(commandArgs(TRUE)[2])
            cat(all(c("EAF", "BETA_TOTAL", "SE_TOTAL", "P_TOTAL", "EAF_AFR", "EAF_EUR") %in% names(a)), sum(!is.na(a$BETA_TOTAL)),
                identical(e$BETA_TOTAL, e$BETA_EUR) && identical(e$SE_TOTAL, e$SE_EUR), "\n")' \
            "${TD}/AAA/chr22.tractor_genesis.genomic_order.tsv.gz" "${TD}/EUR/chr22.tractor_genesis.genomic_order.tsv.gz" 2>&1)"
    [[ "$chk" == "TRUE 100 TRUE " ]] && ok "toy sumstats carry EAF, BETA_TOTAL / SE_TOTAL / P_TOTAL and EAF_<anc>; EUR-only stratum: BETA_TOTAL identical to BETA_EUR" || bad "total columns: $chk"
    hd="$(zcat "${TD}/AAA/chr22.tractor_genesis.genomic_order.tsv.gz" | grep '^##' | cut -d= -f1 | paste -sd' ' -)"
    [[ "$hd" == "## stratum ## trait ## model ## n ## method ## covariates" ]] && grep -q '^## method=.*null: Cox proportional hazards (Breslow ties), no GRM$' <<< "$(zcat "${TD}/AAA/chr22.tractor_genesis.genomic_order.tsv.gz" | head -8)" \
        && ok "'##' lines: stratum, trait, model, n, method (with the null actually fitted), covariates" || bad "header lines: $hd"
    zcat "${TD}/a.dosage.AFR.txt.gz" | awk -F'\t' 'NR > 1 {print $1 "\t" $3 "\t0\t" $2 "\t" $5 "\t" $4}' > "${TD}/a.bim"
    out="$(python3 -I "${C}/bin/sumstats_to_mrmega.py" --input "${TD}/AAA/chr22.tractor_genesis.genomic_order.tsv.gz" --output "${TD}/AAA.mrmega_in.tsv" --variants "${TD}/a.bim" 2>&1)"; rc=$?
    [[ $rc -eq 0 ]] && has "$out" "orientation: 100 of 100 rows checked against" && grep -qx '## beta_source=BETA_TOTAL' "${TD}/AAA.mrmega_in.tsv" \
        && [[ "$(grep -v '^##' "${TD}/AAA.mrmega_in.tsv" | head -1)" == $'MARKERNAME\tCHROMOSOME\tPOSITION\tEA\tNEA\tEAF\tN\tBETA\tSE' ]] \
        && [[ "$(grep -v '^##' "${TD}/AAA.mrmega_in.tsv" | awk -F'\t' 'NR > 1 && !($4 == "G" && $5 == "A")' | wc -l)" -eq 0 ]] \
        && ok "sumstats_to_mrmega.py: 100 rows, EA/NEA = ALT/REF on every row against the variant table; header names BETA_TOTAL as the source" || bad "converter: rc=$rc $out"
    awk -F'\t' 'BEGIN {OFS = "\t"} {split($2, p, ":"); $2 = p[1] ":" p[2] ":" p[4] ":" p[3]; print}' "${TD}/a.bim" > "${TD}/swapped.bim"
    out="$(python3 -I "${C}/bin/sumstats_to_mrmega.py" --input "${TD}/AAA/chr22.tractor_genesis.genomic_order.tsv.gz" --output "${TD}/x.tsv" --variants "${TD}/swapped.bim" 2>&1)"; rc=$?
    [[ $rc -ne 0 ]] && has "$out" "ERROR: orientation: 100 row(s) do not match" && ok "... a variant table with ref and alt the other way round: refused" || bad "swapped: rc=$rc $out"
    # STEP 2b for real (fake Slurm, stub MR-MEGA): K=3 -> MR-MEGA with pc_axes 1; K=2 -> fixed effect
    rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"; mkdir -p "${C}/results/gwas/os" "${C}/logs"
    for s in AFR AMR EUR; do zcat "${TD}/AAA/chr22.tractor_genesis.genomic_order.tsv.gz" | gzip > "${C}/results/gwas/os/${s}.sumstats.gz"; done
    mm() { (cd "$C" && env "${SITE[@]}" GWAS_REPO="$C" LOG_DIR="${C}/logs" FAKE_SBATCH_RUN=1 sbatch --parsable --output="${C}/logs/mm_%j.out" slurm/submit_mrmega.sh os 2>&1 | head -1); }
    jid="$(mm)"; MO="${C}/results/mrmega/os/os.mrmega.result.tsv"
    [[ "$(grep '^##' "$MO" 2>/dev/null | paste -sd'|' -)" == "## method=MR-MEGA|## pc_axes=1|## K=3|## strata=AFR,AMR,EUR|## beta_source=BETA_TOTAL" ]] \
        && grep -q "MR-MEGA (stub) -i ${C}/results/mrmega/os/mrmega_inputs.in --pc 1 --qt -o ${C}/results/mrmega/os/os.mrmega" "${C}/logs/mm_${jid}.out" \
        && ! grep -q '^##' "${C}/results/mrmega/os/EUR.mrmega_in.plain.tsv" \
        && ok "K=3: MR-MEGA --pc 1 --qt on '##'-free copies; result with method / pc_axes / K / strata / beta_source lines" || bad "mrmega K=3: $(tail -3 "${C}/logs/mm_${jid}.out")"
    rm -f "${C}/results/gwas/os/AMR.sumstats.gz"
    jid="$(mm)"; FE="${C}/results/mrmega/os/os.fixed_effect.tsv"
    head -1 "$FE" 2>/dev/null | grep -qx '## pc_axes=0: fixed-effect meta-analysis, not MR-MEGA' && grep -qx '## strata=AFR,EUR' "$FE" \
        && ok "K=2: os.fixed_effect.tsv, its first line 'pc_axes=0: fixed-effect meta-analysis, not MR-MEGA'" || bad "mrmega K=2: $(tail -3 "${C}/logs/mm_${jid}.out")"
    rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
else
    echo "  (real R with data.table, optparse, survival and Matrix, or pandas, not available: CP-19 checks skipped)"
fi

section "CP-23: meta_analyze_strata.R: no single-stratum 'pooled' values, MAC over contributors, METHOD, sorted on the pooled P"
check "STEP 2 passes --method explicitly" grep -q -- '--method random' "${C}/slurm/submit_full_pipeline.sh"
if [[ -n "${REAL_RSCRIPT}" ]] && "${REAL_RSCRIPT}" -e 'quit(status = !all(sapply(c("data.table", "optparse"), requireNamespace, quietly = TRUE)))' >/dev/null 2>&1; then
    TD="${W}/cp23"; rm -rf "${TD:?}"; mkdir -p "$TD"
    "${REAL_RSCRIPT}" -e '
        suppressPackageStartupMessages(library(data.table)); source(file.path(commandArgs(TRUE)[1], "sumstats_io.R")); set.seed(3)
        mk <- function(s, ids) { n <- length(ids); data.table(STRATUM = s, CHR = 1L, POS = 100L * seq_along(ids), ID = ids, REF = "A", ALT = "G", N = 100L,
            MAC_STRATUM = 50, EAF = 0.2, BETA_TOTAL = rnorm(n, 0.1, 0.05), SE_TOTAL = 0.1, P_TOTAL = runif(n), STAT_JOINT = 1, DF_JOINT = 2L, P_JOINT = runif(n),
            MAC_AFR = c(30, NA, rep(20, n - 2)), EAF_AFR = 0.2, BETA_AFR = c(0.2, NA, rnorm(n - 2, 0, 0.1)), SE_AFR = c(0.1, NA, rep(0.1, n - 2)), Z_AFR = 1, P_AFR = 0.3,
            MAC_EUR = 25, EAF_EUR = 0.2, BETA_EUR = rnorm(n, 0, 0.1), SE_EUR = 0.1, Z_EUR = 1, P_EUR = 0.3, Q_HET = 1, DF_HET = 1L, P_HET = 0.3, I2 = 0) }
        a <- mk("S1", c("v1", "v2", "v3", "solo")); b <- mk("S2", c("v1", "v2", "v3"))
        a[ID == "v3", `:=`(BETA_TOTAL = 0.9, SE_TOTAL = 0.05)]; b[ID == "v3", `:=`(BETA_TOTAL = 0.9, SE_TOTAL = 0.05)]
        d <- commandArgs(TRUE)[2]; fwrite_sumstats(a, file.path(d, "S1.sumstats.gz"), c(covariates = "sex")); fwrite_sumstats(b, file.path(d, "S2.sumstats.gz"), c(covariates = "sex"))
        b[, STRATUM := "S9"]; fwrite_sumstats(b, file.path(d, "S3.sumstats.gz"), c(covariates = "sex"))' "${C}/bin" "$TD"
    "${REAL_RSCRIPT}" "${C}/bin/meta_analyze_strata.R" --input_dir "$TD" --strata S1,S2 --method fixed --output "${TD}/P.gz" --threads 1 >/dev/null 2>&1
    out="$("${REAL_RSCRIPT}" -e 'suppressPackageStartupMessages(library(data.table)); p <- fread(commandArgs(TRUE)[1])
        solo <- p[ID == "solo"]; cat(solo$N_STRATA, is.na(solo$BETA_TOTAL) && is.na(solo$P_TOTAL) && is.na(solo$BETA_AFR) && is.na(solo$BETA_EUR), "|",
        p[ID == "v1"]$MAC_AFR, p[ID == "v2"]$MAC_AFR, p[ID == "v3"]$MAC_AFR, "|", p$ID[1], p$P_TOTAL[1] == min(p$P_TOTAL, na.rm = TRUE), "|",
        unique(p[N_STRATA == 2]$METHOD), "|", all(c("Q_HET_NAIVE", "P_HET_NAIVE", "I2_NAIVE") %in% names(p)), any(c("Q_HET", "P_HET", "I2") %in% names(p)), "\n")' "${TD}/P.gz" 2>&1)"
    [[ "$out" == "1 TRUE | 60 NA 40 | v3 TRUE | IVW | TRUE FALSE " ]] && zcat "${TD}/P.gz" | head -6 | grep -qx '## sorted_on=P_TOTAL (pooled total-dose effect)' \
        && ok "a variant in one stratum: N_STRATA 1, every pooled value NA; MAC over contributing strata (60, NA, 40); METHOD set; top row the smallest P_TOTAL; heterogeneity *_NAIVE" || bad "meta CP-23: $out"
    out="$("${REAL_RSCRIPT}" "${C}/bin/meta_analyze_strata.R" --input_dir "$TD" --strata S1,S3 --method fixed --output "${TD}/Q.gz" --threads 1 2>&1)"; rc=$?
    [[ $rc -ne 0 ]] && has "$out" "is pooled as stratum S3 but its STRATUM column says S9" && ok "the input's STRATUM column is kept and checked, not re-derived" || bad "STRATUM check: rc=$rc"
fi

section "CP-16 (extraction): FLARE -> Tractor files, labelled through the index map, checked"
TD="${W}/cp16a"; rm -rf "${TD:?}"; mkdir -p "${TD}/handoff/local_ancestry/flare"
python3 -I "${C}/tests/toy/make_toy_flare.py" "${TD}/handoff/local_ancestry/flare/chr22.anc.vcf.gz" 120 80 EUR,AFR,AMR chr22 3 >/dev/null
printf 'anc_index\tlabel\n0\tEUR\n1\tAFR\n2\tAMR\n' > "${TD}/handoff/local_ancestry/flare_index_map.tsv"
X="${C}/tests/toy/fake_extract_tracts_flare.py"
out="$(bash "${C}/bin/la_flare_to_tractor.sh" "${TD}/handoff/local_ancestry/flare/chr22.anc.vcf.gz" "${TD}/handoff/local_ancestry/flare_index_map.tsv" "${TD}/out/chr22" "$X" 2>&1)"; rc=$?
n=0; for lab in AFR AMR EUR; do for kind in hapcount dosage; do [[ -s "${TD}/out/chr22.${kind}.${lab}.txt.gz" ]] && n=$((n + 1)); done; done
[[ $rc -eq 0 && $n -eq 6 ]] && has "$out" "hapcounts sum to 2, 0 <= dosage <= hapcount everywhere; first 80 records agree with" && has "$out" "the index map matches its ##ANCESTRY line" \
    && [[ "$(zcat "${TD}/out/chr22.hapcount.EUR.txt.gz" | head -1 | cut -f1-6)" == $'CHROM\tPOS\tID\tREF\tALT\tS0001' ]] \
    && ok "chr22: a hapcount AND a dosage file for each of AFR, AMR, EUR (CHROM POS ID REF ALT then IIDs); sums 2, dosage <= hapcount, records agree with GT + AN1/AN2" || bad "extraction: rc=$rc n=$n $(tail -2 <<< "$out")"
printf 'anc_index\tlabel\n0\tAFR\n1\tEUR\n2\tAMR\n' > "${TD}/alphabetical_map.tsv"
out="$(bash "${C}/bin/la_flare_to_tractor.sh" "${TD}/handoff/local_ancestry/flare/chr22.anc.vcf.gz" "${TD}/alphabetical_map.tsv" "${TD}/bad/chr22" "$X" 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "but the FLARE VCF's own ##ANCESTRY line says" && ok "an index map in another order than FLARE's own (e.g. alphabetical): refused" || bad "map order: rc=$rc"
out="$(FAKE_EXTRACT_DOSAGE_ONLY=1 bash "${C}/bin/la_flare_to_tractor.sh" "${TD}/handoff/local_ancestry/flare/chr22.anc.vcf.gz" "${TD}/handoff/local_ancestry/flare_index_map.tsv" "${TD}/d/chr22" "$X" 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "the extraction wrote no hapcount file for ancestry index 0 (EUR)" && ok "an extraction that writes dosage files only: refused (the regression needs both)" || bad "dosage only: rc=$rc"
if [[ -n "${REAL_RSCRIPT}" ]] && "${REAL_RSCRIPT}" -e 'quit(status = !all(sapply(c("data.table", "optparse", "survival", "Matrix"), requireNamespace, quietly = TRUE)))' >/dev/null 2>&1; then
    "${REAL_RSCRIPT}" "${C}/tests/toy/make_toy_tractor.R" "${TD}/ph" 120 20 EUR id_prefix=S >/dev/null    # its phenotype file (IDs S0001..)
    out="$(cd "$TD" && GWAS_QC_DIR="${TD}/qc" "${REAL_RSCRIPT}" "${C}/bin/tractor_genesis_adapter.R" --tractor_prefix out/chr22 --phenotype ph.phenotypes.tsv \
            --trait os --model survival --time_col os_time --event_col os_event --ancestries AFR,AMR,EUR --ref_ancestry EUR \
            --covariates sex --stratum ALL --min_stratum_n 50 --chromosome 22 --mac_min 1 --output_prefix res/chr22 2>&1)"; rc=$?
    [[ $rc -eq 0 ]] && [[ "$(grep -c '(named by label)' <<< "$out")" -eq 3 ]] && ok "the adapter reads them by label, no naming fallback (3 ancestries)" || bad "adapter on extraction: rc=$rc $(grep -i 'named\|error' <<< "$out" | head -3)"
fi
mkdir -p "${C}/logs"
out="$(cd "$C" && env "${SITE[@]}" GWAS_REPO="$C" HANDOFF_DIR="${TD}/handoff" FAKE_SBATCH_RUN=1 sbatch --parsable --array=22 --output="${C}/logs/lae_%A_%a.out" slurm/submit_la_extract.sh 2>&1)"
jid="$(head -1 <<< "$out")"
grep -q "Completed chr22" "${C}/logs/lae_${jid}_22.out" 2>/dev/null && [[ -s "${C}/data/tractor/cohort/chr22.dosage.AMR.txt.gz" ]] && grep -q "activated lai\|lai" "${C}/logs/lae_${jid}_22.out" \
    && ok "slurm/submit_la_extract.sh: chromosome 22's task, in the lai env, writes data/tractor/cohort/chr22.*" || bad "la_extract array: $(tail -3 "${C}/logs/lae_${jid}_22.out" 2>&1)"
rm -rf "${C:?}/data/tractor/cohort" "${C:?}/logs" "${C:?}/results"

section "GxG: missing sumstats are left out and named; none at all stops the task"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"; mkdir -p "${C}/results/gwas/os" "${C}/logs"
echo x | gzip > "${C}/results/gwas/os/POOLED.ivw_descriptive.sumstats.gz"; echo x | gzip > "${C}/results/gwas/os/EUR.sumstats.gz"
out="$(cd "$C" && env "${SITE[@]}" GWAS_REPO="$C" LOG_DIR="${C}/logs" FAKE_SBATCH_RUN=1 sbatch --parsable --array=0-4 \
        --output="${C}/logs/gxg_%A_%a.out" slurm/submit_gxg_array.sh os survival 2>&1)"
jid="$(head -1 <<< "$out")"
lg="${C}/logs/gxg_${jid}_0.out"
has "$(cat "$lg" 2>/dev/null)" "left out of hit selection: AFR AMR MIX" && grep -q "stub Rscript: .*--sumstats ${C}/results/gwas/os/POOLED.ivw_descriptive.sumstats.gz,${C}/results/gwas/os/EUR.sumstats.gz " "$lg" \
    && ok "two of five sumstats present (POOLED + the bundle's 4 strata): the task runs on those two and names the three left out" || bad "gxg partial: $(tail -5 "$lg" 2>&1)"
grep -q -- "--tractor_prefix ${C}/data/tractor/all_chr --ancestries AFR,AMR,EUR" "$lg" && ok "pooled task: --ancestries from data/tractor/ancestries.txt" || bad "gxg ancestries: $(grep stub "$lg")"
grep -q -- "--time_col os_time --event_col os_event " "$lg" && ok "GxG survival columns: os_time / os_event" || bad "gxg columns: $(grep stub "$lg")"
[[ "$(grep -l 'Completed GxG' "${C}"/logs/gxg_"${jid}"_*.out | wc -l)" -eq 5 ]] && ok "all 5 GxG tasks completed (pipefail no longer kills them)" || bad "gxg tasks completed: $(grep -L 'Completed GxG' "${C}"/logs/gxg_"${jid}"_*.out)"
rm -f "${C}"/results/gwas/os/*.sumstats.gz
out="$(cd "$C" && env "${SITE[@]}" GWAS_REPO="$C" LOG_DIR="${C}/logs" FAKE_SBATCH_RUN=1 sbatch --parsable \
        --output="${C}/logs/gxgnone_%A_%a.out" slurm/submit_gxg_array.sh os survival 2>&1)"
jid="$(head -1 <<< "$out")"
grep -q "ERROR: no sumstats for os" "${C}/logs/gxgnone_${jid}_0.out" && grep -q "^${jid}_0|FAILED|1:0" "${FAKE_SLURM_DIR}/sacct.db" \
    && ok "no sumstats at all: the task fails, saying so" || bad "gxg none: $(tail -3 "${C}/logs/gxgnone_${jid}_0.out")"

section "PRS: ancestries from data/tractor/ancestries.txt; proportions must name the same ones"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/logs"; mkdir -p "${C}/logs"
out="$(cd "$C" && env "${SITE[@]}" GWAS_REPO="$C" FAKE_SBATCH_RUN=1 sbatch --parsable --output="${C}/logs/prs_%A_%a.out" slurm/submit_prs_array.sh os 2>&1)"
jid="$(head -1 <<< "$out")"
grep -q -- "--ancestries AFR,AMR,EUR --ancestry_props EUR:0.35,AFR:0.15,AMR:0.50" "${C}/logs/prs_${jid}_1.out" && ok "PRS gets --ancestries AFR,AMR,EUR from the file" || bad "prs ancestries: $(grep stub "${C}/logs/prs_${jid}_1.out")"
grep -q -- "--trait os --trait_type survival --time_col os_time --event_col os_event " "${C}/logs/prs_${jid}_1.out" && ok "PRS validation gets os_time / os_event" || bad "prs columns: $(grep stub "${C}/logs/prs_${jid}_1.out")"
printf 'AFR\nEUR\n' > "${C}/data/tractor/ancestries.txt"
out="$(cd "$C" && env "${SITE[@]}" GWAS_REPO="$C" FAKE_SBATCH_RUN=1 sbatch --parsable --output="${C}/logs/prsbad_%A_%a.out" slurm/submit_prs_array.sh os 2>&1)"
jid="$(head -1 <<< "$out")"
grep -q "PRS_ANCESTRY_PROPS (EUR:0.35,AFR:0.15,AMR:0.50) names AFR,AMR,EUR, but .*lists AFR,EUR" "${C}/logs/prsbad_${jid}_1.out" && ok "proportions naming other ancestries than the file: refused" || bad "prs mismatch: $(tail -2 "${C}/logs/prsbad_${jid}_1.out")"
printf 'AFR\nAMR\nEUR\n' > "${C}/data/tractor/ancestries.txt"

section "Colocalization: per trait, refuses without coloc / susieR / hyprcoloc"
out="$(cd /tmp && env "${SITE[@]}" RUN_COLOC=true "$L" --dry-run 2>&1)"
cl="$(grep '\[DRY-RUN\] colocalization' <<< "$out")"
[[ "$(wc -l <<< "$cl")" -eq 3 ]] && grep -q 'colocalization mrd .*--dependency=afterok:dry14 .*coloc_mrd_%A_%a.out' <<< "$cl" \
    && grep -q 'colocalization os .*--dependency=afterok:dry16 ' <<< "$cl" && ok "one coloc array per trait, each waiting on its own combine job" || bad "coloc dry run: $cl"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=os RUN_COLOC=true FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 "$L" 2>&1)"
grep -q "stub Rscript: .*run_colocalization.R --gwas ${C}/results/gwas/os/POOLED.ivw_descriptive.sumstats.gz .*--output_prefix ${C}/results/coloc/os/eqtl " "${C}"/logs/coloc_os_*_1.out \
    && ok "coloc reads results/gwas/os/POOLED.ivw_descriptive.sumstats.gz and writes results/coloc/os/" || bad "coloc paths: $(grep stub "${C}"/logs/coloc_os_*_1.out)"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=os RUN_COLOC=true FAKE_R_MISSING=hyprcoloc FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 "$L" 2>&1)"
RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"
grep -q "REFUSED: colocalization needs the R packages coloc, susieR and hyprcoloc" "${C}"/logs/coloc_os_*_1.out && grep -q "MISSING hyprcoloc" "${C}"/logs/coloc_os_*_1.out \
    && [[ "$(awk -F'\t' 'NR > 1 && $1 == "coloc" && $14 == 3' "${RD%/}/trace.tsv" | wc -l)" -eq 6 ]] \
    && ok "hyprcoloc missing: all 6 coloc tasks refuse (exit 3), naming it" || bad "coloc refusal: $(tail -3 "${C}"/logs/coloc_os_*_1.out)"

# ============================================================================
echo ""
echo "passed ${pass}, failed ${failn}"
[[ $failn -eq 0 ]]
