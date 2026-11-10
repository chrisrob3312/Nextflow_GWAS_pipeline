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
export FAKE_SLURM_DIR="${W}/slurm"; mkdir -p "$FAKE_SLURM_DIR"
export PATH="${C}/tests/fakes:${PATH}"
export GWAS_ENVS_DIR="${W}/envs" GWAS_SCRATCH="${W}/scratch"
mkdir -p "${W}/envs/"{gwas-r,gwas-py,lai,mrmega}/bin "${W}/scratch"
# Stub Rscript inside the gwas-r env: emulates bin/check_r_libs.R (packages in FAKE_R_MISSING are
# missing), and stands in for every analysis script (prints its arguments, exits FAKE_RSCRIPT_EXIT).
cat > "${W}/envs/gwas-r/bin/Rscript" <<'EOF'
#!/usr/bin/env bash
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
    if [[ "$st" == ok ]]; then nt=900; else nt=NA; fi
    printf 'key\tvalue\nstatus\t%s\nstratum\tEUR\nstratum_check\tpassed\nmin_stratum_n\t50\nn_null\t%s\nn_variants\t1000\nn_pass_mac\t950\nn_tested\t%s\nn_failed\t2\n' \
        "$st" "${FAKE_N_NULL:-120}" "$nt" > "${out}.run_counts.tsv"
fi
exit "${FAKE_RSCRIPT_EXIT:-0}"
EOF
printf '#!/bin/sh\necho "openjdk version \\"21.0.10\\" (stub)" >&2\n' > "${W}/envs/lai/bin/java"
printf '#!/bin/sh\necho "MR-MEGA (stub)"\n' > "${W}/envs/mrmega/bin/MR-MEGA"
for e in gwas-r lai; do printf '#!/bin/sh\necho "bcftools 1.25 (stub)"\necho "Using htslib 1.25"\n' > "${W}/envs/${e}/bin/bcftools"; done
chmod +x "${W}/envs/gwas-r/bin/Rscript" "${W}/envs/lai/bin/java" "${W}/envs/mrmega/bin/MR-MEGA" "${W}"/envs/*/bin/bcftools
STRATA_DEFAULT=(EUR AAC LAT1 LAT2 OTHER)
for s in "${STRATA_DEFAULT[@]}"; do mkdir -p "${C}/data/tractor/${s}"; done
printf 'EUR\n' > "${C}/data/tractor/EUR/ancestries.txt"
printf 'AFR\nEUR\n' > "${C}/data/tractor/AAC/ancestries.txt"
for s in LAT1 LAT2 OTHER; do printf 'AFR\nAMR\nEUR\n' > "${C}/data/tractor/${s}/ancestries.txt"; done
printf 'AFR\nAMR\nEUR\n' > "${C}/data/tractor/ancestries.txt"       # cohort-wide list (GxG pooled, PRS)
L="${C}/slurm/submit_full_pipeline.sh"
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
has "$out" "deferred: PRS COLOC VALIDATION GXG (RUN_*=false)" && ok "prints 'deferred: PRS COLOC VALIDATION GXG (RUN_*=false)'" || bad "deferred line missing"
grep -q 'submit_prs_array\|colocalization\|PRS validation\|submit_gxg_array' <<< "$sb" && bad "a deferred step was submitted" || ok "no sbatch for PRS, coloc, validation or GxG"
has "$sb" "run record finalise: sbatch" && grep -q -- '--dependency=afterany:dry1:dry2:' <<< "$(grep 'finalise' <<< "$sb")" && ok "finalise job waits on every job (afterany)" || bad "finalise line"
[[ ! -e "${C}/results" ]] && ok "a dry run creates no run record" || bad "dry run wrote ${C}/results"

section "CP-03: RUN_* gates"
out="$(cd /tmp && env "${SITE[@]}" RUN_PRS=true "$L" --dry-run 2>&1)"
grep -q 'submit_prs_array.sh OS' <<< "$(grep '\[DRY-RUN\]' <<< "$out")" && ok "RUN_PRS=true prints the PRS sbatch" || bad "RUN_PRS=true: no PRS sbatch"
has "$out" "deferred: COLOC VALIDATION GXG (RUN_*=false)" && ok "RUN_PRS=true: deferred list shrinks" || bad "deferred with RUN_PRS"
out="$(cd /tmp && env "${SITE[@]}" RUN_GXG=true "$L" --dry-run 2>&1)"
gx="$(grep '\[DRY-RUN\] GxG' <<< "$out")"
[[ "$(grep -c -- '--dependency=afterok:dry17:dry18:dry19 ' <<< "$gx")" -eq 3 ]] && ok "RUN_GXG without PRS: GxG waits on the three combine jobs" || bad "GxG dependency: $gx"
grep -q -- 'afterok: \|afterok::\|afterok:$\|--dependency= ' <<< "$out" && bad "an empty job id in a dependency" || ok "no empty job id in any dependency"
out="$(cd /tmp && env "${SITE[@]}" RUN_PRS=true RUN_COLOC=true RUN_VALIDATION=true RUN_GXG=true "$L" --dry-run 2>&1)"
has "$out" "deferred: none" && [[ "$(grep -c '\[DRY-RUN\] \(PRS\|colocalization\|PRS validation\|GxG\)' <<< "$out")" -eq 10 ]] && ok "all RUN_*=true: every step submitted (coloc once per trait)" || bad "all on"
out="$(cd /tmp && env "${SITE[@]}" RUN_VALIDATION=true "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "RUN_VALIDATION=true needs RUN_PRS=true" && ok "RUN_VALIDATION without RUN_PRS: refused" || bad "validation alone: rc=$rc"
out="$(cd /tmp && env "${SITE[@]}" RUN_PRS=yes "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "RUN_PRS=yes; use true or false" && ok "RUN_PRS=yes: refused" || bad "RUN_PRS=yes: rc=$rc"
out="$(cd /tmp && env "${SITE[@]}" STRATUM_MIN_N=50 "$L" --dry-run 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "the scripts read MIN_STRATUM_N" && ok "STRATUM_MIN_N (read by no script): refused" || bad "STRATUM_MIN_N: rc=$rc"

# ============================================================================
section "CP-28: submit_gwas_array.sh reads data/tractor/<stratum>/ancestries.txt"
G="${C}/slurm/submit_gwas_array.sh"
out="$(cd /tmp && DRY_RUN=1 bash "$G" EUR OS survival 2>&1)"; rc=$?
[[ $rc -eq 0 ]] && grep -q -- '--ancestries EUR ' <<< "$(grep '^\[DRY-RUN\]' <<< "$out")" && ok "ancestries.txt 'EUR': prints --ancestries EUR" || bad "EUR: rc=$rc $out"
out="$(cd /tmp && DRY_RUN=1 SLURM_ARRAY_TASK_ID=5 bash "$G" AAC MRD binary 2>&1)"
grep -q -- '--ancestries AFR,EUR .*--chromosome 5 ' <<< "$(grep '^\[DRY-RUN\]' <<< "$out")" && ok "two labels, one per line: --ancestries AFR,EUR (chr from the task id)" || bad "AAC: $out"
mv "${C}/data/tractor/LAT1/ancestries.txt" "${W}/anc.away"
out="$(cd /tmp && DRY_RUN=1 bash "$G" LAT1 OS survival 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "${C}/data/tractor/LAT1/ancestries.txt does not exist" && ok "without the file: exits non-zero naming it" || bad "absent: rc=$rc $out"
out="$(cd /tmp && env "${SITE[@]}" "$L" --dry-run 2>&1)"
has "$out" "WARNING: missing ${C}/data/tractor/LAT1/ancestries.txt" && ok "launcher dry run warns about the missing file" || bad "launcher warning"
out="$(cd /tmp && env "${SITE[@]}" "$L" 2>&1)"; rc=$?
[[ $rc -ne 0 ]] && has "$out" "REFUSED: missing ${C}/data/tractor/LAT1/ancestries.txt" && ok "launcher refuses a real launch without it" || bad "launcher refusal: rc=$rc"
mv "${W}/anc.away" "${C}/data/tractor/LAT1/ancestries.txt"
rm -rf "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && DRY_RUN=1 MAC_MIN=30 MIN_STRATUM_N=50 bash "$G" EUR OS survival 2>&1)"
grep -q -- '--min_stratum_n 50 .*--mac_min 30 ' <<< "$(grep '^\[DRY-RUN\]' <<< "$out")" && ok "MAC_MIN and MIN_STRATUM_N reach the adapter" || bad "switches: $out"
check "shell/README.md says 'debug-only, not the execution path'" grep -q "debug-only, not the execution path" "$REPO/shell/README.md"
for f in "$REPO"/shell/*.sh; do check "$(basename "$f"): one-line DEBUG-ONLY header" bash -c "sed -n 2p '$f' | grep -q '^# DEBUG-ONLY, not the execution path'"; done

# ============================================================================
section "B-03: the switch table matches the scripts' own defaults"
out="$(bash -c 'source "$1/bin/run_record.sh"; printf "%s\n" "${RR_SWITCHES[@]}"' _ "$C" | while IFS='|' read -r name def _; do
    [[ "$def" == "@job" ]] && continue
    grep -rhoE "\\\$\{${name}:-[^}]*\}" "$C/slurm/" | sort -u | while read -r m; do
        got="${m#*:-}"; got="${got%\}}"
        [[ "$got" == "$def" ]] || echo "${name}: table '${def}', script '${got}'"
    done
done)"
[[ -z "$out" ]] && ok "every \${NAME:-default} in slurm/ equals RR_SWITCHES" || bad "defaults differ: $out"

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
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=OS FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 "$L" 2>&1)"; rc=$?
E2E="$(ls -d "${C}"/results/run_record/*/ 2>/dev/null | head -1)"; E2E="${E2E%/}"
T="${E2E}/trace.tsv"
[[ $rc -eq 0 && -s "$T" ]] && ok "launch ran PC-AiR, 22 GWAS tasks, combine and finalise" || bad "e2e: rc=$rc $(tail -5 <<< "$out")"
steps="$(awk -F'\t' 'NR > 1 {n[$1]++} END {for (s in n) printf "%s=%d ", s, n[s]}' "$T" | tr ' ' '\n' | sort | paste -sd' ' -)"
[[ "$steps" == "combine=1 gwas=22 pcair=1" ]] && ok "trace rows per step: ${steps}" || bad "steps: ${steps}"
[[ "$(awk -F'\t' 'NR > 1 && $14 != "0"' "$T" | wc -l)" -eq 0 ]] && ok "every task exit 0, recorded by the task itself" || bad "exits: $(awk -F'\t' 'NR>1{print $1, $14}' "$T" | sort | uniq -c)"
grep -q -- '--ancestries EUR .*--mac_min 10 ' <<< "$(awk -F'\t' 'NR > 1 && $1 == "gwas" {print $10}' "$T")" && ok "gwas rows carry the adapter command line" || bad "gwas args"
grep -q 'ERROR: cannot find the repo' "${C}"/logs/*.out 2>/dev/null && bad "a job could not find the repo from the spool copy" || ok "every job found the repo while running from a spool copy"
[[ "$(grep -l 'activate_env: gwas-r' "${C}"/logs/gwas_OS_EUR_*.out 2>/dev/null | wc -l)" -eq 22 ]] && ok "all 22 GWAS tasks activated gwas-r" || bad "activation in GWAS logs"

section "CP-02: one script submitted by hand (no launcher, no GWAS_REPO)"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"; mkdir -p "${C}/logs"
out="$(cd "$C" && env "${SITE[@]}" FAKE_SBATCH_RUN=1 sbatch --parsable --partition=testpart --account=testacct \
        --output="${C}/logs/pcair_%j.out" slurm/submit_pcair.sh 2>&1)"
jid="$(tail -1 <<< "$out")"
grep -q 'stub Rscript: .*genesis_pcair_pcrelate.R' "${C}/logs/pcair_${jid}.out" 2>/dev/null \
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
grep -q 'results\[, STRATUM := opt$stratum\]' "$A" && ok "adapter: outputs labelled with a STRATUM column" || bad "adapter: no STRATUM label"
grep -q 'write_run_counts(opt$output_prefix, run_counts)' "$A" && ok "adapter: writes <output_prefix>.run_counts.tsv" || bad "adapter: no counts file"
grep -q '".ancdose.", k' "$REPO/bin/calculate_prs_admixed.R" && bad "PRS still reads Tractor files by position" || ok "PRS LA partial: no file number from the position either"
for m in modules/local/tractor_genesis_adapter/main.nf modules/local/prs_admixed/main.nf; do
    grep -q -- '--ancestry_index_map ancestry_index_map.tsv' "$REPO/$m" && ok "${m}: passes its declared order as an explicit index map" || bad "${m}: no index map"
done
grep -q 'EUR,AFR,AMR' "$REPO/slurm/submit_gxg_array.sh" "$REPO/slurm/submit_prs_array.sh" && bad "GxG/PRS still hard-code EUR,AFR,AMR" || ok "GxG and PRS no longer hard-code EUR,AFR,AMR"
grep -q '^TRAITS=' "$REPO/bin/slurm_utils.sh" && bad "slurm_utils.sh still sets TRAITS" || ok "slurm_utils.sh no longer sets or exports TRAITS"

section "GWAS array: strata file and index map reach the adapter when present"
out="$(cd /tmp && DRY_RUN=1 bash "$G" EUR OS survival 2>&1)"
has "$out" "data/ancestry/strata.tsv not found" && ! grep -q -- '--strata_file' <<< "$(grep '^\[DRY-RUN\]' <<< "$out")" && ok "no strata file: not passed, and the log says the check falls back" || bad "no strata file: $out"
mkdir -p "${C}/data/ancestry"; printf 'sample_id\tstratum\nS1\tEUR\n' > "${C}/data/ancestry/strata.tsv"
printf 'anc_index\tlabel\n0\tAFR\n1\tAMR\n2\tEUR\n' > "${C}/data/tractor/flare_index_map.tsv"
out="$(cd /tmp && DRY_RUN=1 bash "$G" EUR OS survival 2>&1 | grep '^\[DRY-RUN\]')"
has "$out" "--strata_file ${C}/data/ancestry/strata.tsv" && has "$out" "--ancestry_index_map ${C}/data/tractor/flare_index_map.tsv" && has "$out" "--min_stratum_n 50 " \
    && ok "strata file, index map and MIN_STRATUM_N 50 passed" || bad "args: $out"
rm -f "${C}/data/ancestry/strata.tsv" "${C}/data/tractor/flare_index_map.tsv"

section "Run record: the adapter's counts and a stratum skipped below MIN_STRATUM_N"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=OS FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 "$L" 2>&1)"
RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"; T="${RD%/}/trace.tsv"
row="$(awk -F'\t' 'NR > 1 && $1 == "gwas" {print $15, $16, $17, $18, $23; exit}' "$T")"
[[ "$row" == "120 900 950 2 ok" ]] && ok "gwas rows: n_null 120, n_tested 900, n_pass_mac 950, n_failed 2, task_status ok (from run_counts.tsv)" || bad "counts row: '$row'"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=OS FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 FAKE_ADAPTER_STATUS=skipped_below_min_stratum_n FAKE_N_NULL=31 "$L" 2>&1)"
RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"
grep -q '| gwas | 22 | 22 | 0 | 0 | 22 |' "${RD%/}/summary.md" && ok "summary.md: 22 gwas tasks exit 0, all 22 skipped below MIN_STRATUM_N" || bad "summary: $(grep '| gwas' "${RD%/}/summary.md")"
grep -q 'skipped: N below MIN_STRATUM_N' "${C}"/logs/combine_OS_*.out && ok "combine names the likely reason for a stratum without results" || bad "combine message"

section "GxG: missing sumstats are left out and named; none at all stops the task"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"; mkdir -p "${C}/results/gwas/OS" "${C}/logs"
echo x | gzip > "${C}/results/gwas/OS/POOLED.sumstats.gz"; echo x | gzip > "${C}/results/gwas/OS/EUR.sumstats.gz"
out="$(cd "$C" && env "${SITE[@]}" GWAS_REPO="$C" LOG_DIR="${C}/logs" FAKE_SBATCH_RUN=1 sbatch --parsable \
        --output="${C}/logs/gxg_%A_%a.out" slurm/submit_gxg_array.sh OS survival 2>&1)"
jid="$(head -1 <<< "$out")"
lg="${C}/logs/gxg_${jid}_0.out"
has "$(cat "$lg" 2>/dev/null)" "left out of hit selection: AAC LAT1 LAT2 EAS SAS OTHER" && grep -q "stub Rscript: .*--sumstats ${C}/results/gwas/OS/POOLED.sumstats.gz,${C}/results/gwas/OS/EUR.sumstats.gz " "$lg" \
    && ok "two of eight sumstats present: the task runs on those two and names the six left out" || bad "gxg partial: $(tail -5 "$lg" 2>&1)"
grep -q -- "--tractor_prefix ${C}/data/tractor/all_chr --ancestries AFR,AMR,EUR" "$lg" && ok "pooled task: --ancestries from data/tractor/ancestries.txt" || bad "gxg ancestries: $(grep stub "$lg")"
[[ "$(grep -l 'Completed GxG' "${C}"/logs/gxg_"${jid}"_*.out | wc -l)" -eq 8 ]] && ok "all 8 GxG tasks completed (pipefail no longer kills them)" || bad "gxg tasks completed: $(grep -L 'Completed GxG' "${C}"/logs/gxg_"${jid}"_*.out)"
rm -f "${C}"/results/gwas/OS/*.sumstats.gz
out="$(cd "$C" && env "${SITE[@]}" GWAS_REPO="$C" LOG_DIR="${C}/logs" FAKE_SBATCH_RUN=1 sbatch --parsable \
        --output="${C}/logs/gxgnone_%A_%a.out" slurm/submit_gxg_array.sh OS survival 2>&1)"
jid="$(head -1 <<< "$out")"
grep -q "ERROR: no sumstats for OS" "${C}/logs/gxgnone_${jid}_0.out" && grep -q "^${jid}_0|FAILED|1:0" "${FAKE_SLURM_DIR}/sacct.db" \
    && ok "no sumstats at all: the task fails, saying so" || bad "gxg none: $(tail -3 "${C}/logs/gxgnone_${jid}_0.out")"

section "PRS: ancestries from data/tractor/ancestries.txt; proportions must name the same ones"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/logs"; mkdir -p "${C}/logs"
out="$(cd "$C" && env "${SITE[@]}" GWAS_REPO="$C" FAKE_SBATCH_RUN=1 sbatch --parsable --output="${C}/logs/prs_%A_%a.out" slurm/submit_prs_array.sh OS 2>&1)"
jid="$(head -1 <<< "$out")"
grep -q -- "--ancestries AFR,AMR,EUR --ancestry_props EUR:0.35,AFR:0.15,AMR:0.50" "${C}/logs/prs_${jid}_1.out" && ok "PRS gets --ancestries AFR,AMR,EUR from the file" || bad "prs ancestries: $(grep stub "${C}/logs/prs_${jid}_1.out")"
printf 'AFR\nEUR\n' > "${C}/data/tractor/ancestries.txt"
out="$(cd "$C" && env "${SITE[@]}" GWAS_REPO="$C" FAKE_SBATCH_RUN=1 sbatch --parsable --output="${C}/logs/prsbad_%A_%a.out" slurm/submit_prs_array.sh OS 2>&1)"
jid="$(head -1 <<< "$out")"
grep -q "PRS_ANCESTRY_PROPS (EUR:0.35,AFR:0.15,AMR:0.50) names AFR,AMR,EUR, but .*lists AFR,EUR" "${C}/logs/prsbad_${jid}_1.out" && ok "proportions naming other ancestries than the file: refused" || bad "prs mismatch: $(tail -2 "${C}/logs/prsbad_${jid}_1.out")"
printf 'AFR\nAMR\nEUR\n' > "${C}/data/tractor/ancestries.txt"

section "Colocalization: per trait, refuses without coloc / susieR / hyprcoloc"
out="$(cd /tmp && env "${SITE[@]}" RUN_COLOC=true "$L" --dry-run 2>&1)"
cl="$(grep '\[DRY-RUN\] colocalization' <<< "$out")"
[[ "$(wc -l <<< "$cl")" -eq 3 ]] && grep -q 'colocalization OS .*--dependency=afterok:dry17 .*coloc_OS_%A_%a.out' <<< "$cl" \
    && grep -q 'colocalization MRD .*--dependency=afterok:dry19 ' <<< "$cl" && ok "one coloc array per trait, each waiting on its own combine job" || bad "coloc dry run: $cl"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=OS RUN_COLOC=true FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 "$L" 2>&1)"
grep -q "stub Rscript: .*run_colocalization.R --gwas ${C}/results/gwas/OS/POOLED.sumstats.gz .*--output_prefix ${C}/results/coloc/OS/eqtl " "${C}"/logs/coloc_OS_*_1.out \
    && ok "coloc reads results/gwas/OS/POOLED.sumstats.gz and writes results/coloc/OS/" || bad "coloc paths: $(grep stub "${C}"/logs/coloc_OS_*_1.out)"
rm -rf "${FAKE_SLURM_DIR:?}"/* "${C:?}/results" "${C:?}/logs"
out="$(cd /tmp && env "${SITE[@]}" STRATA_LIST=EUR TRAITS_LIST=OS RUN_COLOC=true FAKE_R_MISSING=hyprcoloc FAKE_SBATCH_RUN=1 RR_SACCT_WAIT=0 "$L" 2>&1)"
RD="$(ls -d "${C}"/results/run_record/*/ | head -1)"
grep -q "REFUSED: colocalization needs the R packages coloc, susieR and hyprcoloc" "${C}"/logs/coloc_OS_*_1.out && grep -q "MISSING hyprcoloc" "${C}"/logs/coloc_OS_*_1.out \
    && [[ "$(awk -F'\t' 'NR > 1 && $1 == "coloc" && $14 == 3' "${RD%/}/trace.tsv" | wc -l)" -eq 6 ]] \
    && ok "hyprcoloc missing: all 6 coloc tasks refuse (exit 3), naming it" || bad "coloc refusal: $(tail -3 "${C}"/logs/coloc_OS_*_1.out)"

# ============================================================================
echo ""
echo "passed ${pass}, failed ${failn}"
[[ $failn -eq 0 ]]
