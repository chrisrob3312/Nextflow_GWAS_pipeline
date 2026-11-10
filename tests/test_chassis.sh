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
has "$sb" "run record finalise: sbatch" && grep 'finalise' <<< "$sb" | grep -q -- '--dependency=afterany:dry1:dry2:' && ok "finalise job waits on every job (afterany)" || bad "finalise line"
[[ ! -e "${C}/results" ]] && ok "a dry run creates no run record" || bad "dry run wrote ${C}/results"

section "CP-03: RUN_* gates"
out="$(cd /tmp && env "${SITE[@]}" RUN_PRS=true "$L" --dry-run 2>&1)"
grep '\[DRY-RUN\]' <<< "$out" | grep -q 'submit_prs_array.sh OS' && ok "RUN_PRS=true prints the PRS sbatch" || bad "RUN_PRS=true: no PRS sbatch"
has "$out" "deferred: COLOC VALIDATION GXG (RUN_*=false)" && ok "RUN_PRS=true: deferred list shrinks" || bad "deferred with RUN_PRS"
out="$(cd /tmp && env "${SITE[@]}" RUN_GXG=true "$L" --dry-run 2>&1)"
gx="$(grep '\[DRY-RUN\] GxG' <<< "$out")"
[[ "$(grep -c -- '--dependency=afterok:dry17:dry18:dry19 ' <<< "$gx")" -eq 3 ]] && ok "RUN_GXG without PRS: GxG waits on the three combine jobs" || bad "GxG dependency: $gx"
grep -q -- 'afterok: \|afterok::\|afterok:$\|--dependency= ' <<< "$out" && bad "an empty job id in a dependency" || ok "no empty job id in any dependency"
out="$(cd /tmp && env "${SITE[@]}" RUN_PRS=true RUN_COLOC=true RUN_VALIDATION=true RUN_GXG=true "$L" --dry-run 2>&1)"
has "$out" "deferred: none" && [[ "$(grep -c '\[DRY-RUN\] \(PRS\|colocalization\|PRS validation\|GxG\)' <<< "$out")" -eq 8 ]] && ok "all RUN_*=true: every step submitted" || bad "all on"
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
[[ $rc -eq 0 ]] && grep '^\[DRY-RUN\]' <<< "$out" | grep -q -- '--ancestries EUR ' && ok "ancestries.txt 'EUR': prints --ancestries EUR" || bad "EUR: rc=$rc $out"
out="$(cd /tmp && DRY_RUN=1 SLURM_ARRAY_TASK_ID=5 bash "$G" AAC MRD binary 2>&1)"
grep '^\[DRY-RUN\]' <<< "$out" | grep -q -- '--ancestries AFR,EUR .*--chromosome 5 ' && ok "two labels, one per line: --ancestries AFR,EUR (chr from the task id)" || bad "AAC: $out"
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
grep '^\[DRY-RUN\]' <<< "$out" | grep -q -- '--min_stratum_n 50 .*--mac_min 30 ' && ok "MAC_MIN and MIN_STRATUM_N reach the adapter" || bad "switches: $out"
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
[[ "$(col n_tested)" == "NA NA" ]] && ok "adapter counts NA (the adapter writes no sidecar)" || bad "n_tested: $(col n_tested)"
grep -q '| stub | 2 | 1 | 1 | 0 |' "${SELF}/summary.md" && ok "summary.md: stub step, 2 tasks, 1 ok, 1 non-zero" || bad "summary: $(cat "${SELF}/summary.md")"
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
awk -F'\t' 'NR > 1 && $1 == "gwas" {print $10}' "$T" | grep -q -- '--ancestries EUR .*--mac_min 10 ' && ok "gwas rows carry the adapter command line" || bad "gwas args"
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
echo ""
echo "passed ${pass}, failed ${failn}"
[[ $failn -eq 0 ]]
