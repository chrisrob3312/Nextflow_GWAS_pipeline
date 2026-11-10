# shellcheck shell=bash
# ============================================================================
# run_record.sh: what ran, where, in which env, on which inputs, and how it ended (B-03)
# ============================================================================
# Sourced, never run: by bin/slurm_utils.sh, so every job has it, and by the launcher.
#
# One directory per launch, results/run_record/<run_id>/:
#   params.tsv    every number-changing switch: value, and whether you set it or it is the default
#   run_info.tsv  run id, git commit and whether the tree had uncommitted changes, launch time,
#                 host, partition, account, deferred steps
#   jobs.tsv      every job the launcher submitted (label, job id)
#   trace.d/      one file per task, written only by that task (see "Why no locking" below)
#   trace.tsv     every task's row, plus State, MaxRSS and Elapsed from sacct, and the task's
#                 own status (e.g. skipped_below_min_stratum_n)                (finalise)
#   summary.md    tasks per step, failures, deferred steps                       (finalise)
#
# run_id is <yyyymmdd-HHMM>_<git commit, 7 characters>_<n>. The launcher computes it once and
# exports GWAS_RUN_ID and GWAS_RUN_DIR; jobs inherit both (sbatch --export=ALL), so every job of
# one launch writes into the same directory. <n> keeps two launches in the same minute on the same
# commit apart: the directory is claimed with `mkdir`, which fails if the directory already exists,
# so a second launch gets the next n and never writes into the first one's directory.
#
# Why no locking: concurrent array tasks never append to a shared file. Each task writes its own
# trace.d/<step>.<slurm id>.tsv (to a temporary name, then `mv`, which replaces the file in one
# step), and finalise concatenates them. `flock` would also work, but on Lustre/GPFS/NFS it depends
# on mount options; one file per task needs nothing from the filesystem.
#
# In a job, first thing after sourcing slurm_utils.sh, before activating the env, so that a task
# that fails while activating still leaves a row:
#   rr_begin STEP ["ARGS" [INPUT_FILE...]]
#   rr_args "ARGS" [INPUT_FILE...]          later, once the command line and inputs are known
#   rr_sidecar FILE                         the counts file the program writes (the adapter's
#                                           <output_prefix>.run_counts.tsv), read when the task ends
# rr_begin writes the task's row with an empty end time and exit status, then sets an EXIT trap
# (a command bash runs when the script ends, however it ends) that rewrites the row with the end
# time, the exit status and the env that was active. A job killed by Slurm (time limit, scancel) gets SIGTERM first; the TERM trap turns
# that into exit status 143 so the row says so. A job killed outright (SIGKILL, e.g. out of memory)
# keeps the empty exit status, and finalise fills it in from sacct.
# ============================================================================

# Sourced once per shell: a second source (slurm_utils.sh sources this file too) would reset
# RR_SET_BY, the record of which switches you set.
if [[ -n "${_RR_LOADED:-}" ]]; then return 0; fi
_RR_LOADED=1

_RR_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Number-changing switches, with the defaults the scripts use. One line each:
#   NAME|DEFAULT|what it changes
# DEFAULT "@job" means the job chooses at run time (e.g. PHENO is the PC-AiR output if that exists
# when the job starts), so nothing is exported and params.tsv records it only if you set it.
# tests/test_chassis.sh checks that each default here equals the ${NAME:-default} in the scripts.
# Asked for in the plan but read by no script, so not listed: STRATUM_MIN_N (the scripts use
# MIN_STRATUM_N), PC_SCALE, MRMEGA_PCS, BUNDLE_RULES; rr_refuse_unknown_switches stops a launch
# that sets one of them, rather than let it be silently ignored.
RR_SWITCHES=(
    "STRATA_LIST|EUR,AAC,LAT1,LAT2,OTHER|strata analysed"
    "TRAITS_LIST|mrd,efs,os|traits analysed (phenotype columns, as docs/HANDOFF_CONTRACT.md names them)"
    "SURVIVAL_TRAITS|efs,os|traits run as time-to-event, with columns <trait>_time and <trait>_event; the rest are binary"
    "COVARIATES|age,sex,PC1,PC2,PC3,PC4,PC5|model covariates (GWAS, GxG)"
    "MIN_STRATUM_N|50|smallest stratum analysed (--min_stratum_n); 50 decided by the user, 10 Nov 2026"
    "MAC_MIN|10|minor allele count needed to test an ancestry's dose term (--mac_min); open: the contract says 20"
    "RUN_PRS|false|step 3, PRS"
    "RUN_COLOC|false|step 4, colocalization"
    "RUN_VALIDATION|false|step 5, PRS validation"
    "RUN_GXG|false|step 6, GxG"
    "SKIP_PCAIR|false|skip step 0 and use data/kinship and the phenotype file's PCs"
    "N_PCS|10|PC-AiR PCs written and used as covariates"
    "N_PCS_PCRELATE|5|PCs used inside PC-Relate"
    "KIN_THRESH|0.0442|kinship threshold for the PC-AiR unrelated set"
    "LD_R2|0.1|LD pruning r2 for PC-AiR / PC-Relate"
    "PCAIR_ITER|2|PC-AiR / PC-Relate iterations"
    "GXG_P_THRESHOLD|1e-5|GWAS P below which a hit enters GxG"
    "GXG_MAX_HITS|200|cap on GxG hits"
    "GXG_PRUNE_MODE|conditional|how LD-correlated GxG hits are selected"
    "GXG_COND_P|1e-4|conditional P that keeps an LD-correlated GxG hit"
    "GXG_USE_DEFAULT_LOCI|true|include the known-leukemia-loci table in GxG"
    "GXG_CUSTOM_VARIANTS||extra GRCh38 variant list for GxG"
    "ANCESTRY_COL|GRAF_ANC|phenotype column holding the ancestry label (GxG)"
    "PHENO|@job|phenotype file (default: PC-AiR output if present, else data/phenotypes/phenotypes.tsv)"
    "KINSHIP|@job|GRM (default: PC-Relate output if present, else data/kinship/kinship_matrix.rds)"
    "GENO|@job|PLINK prefix for PC-AiR (default data/genotypes/cohort)"
    "PHENO_IN|@job|phenotype file PC-AiR reads (default data/phenotypes/phenotypes.tsv)"
    "KNOWN_LOCI|@job|known-loci table for GxG (default assets/known_leukemia_risk_loci.tsv)"
    "STRATA_FILE|@job|sample_id/stratum table the adapter checks --stratum against (default data/ancestry/strata.tsv)"
    "ANC_INDEX_MAP|@job|anc_index/label map for numbered Tractor files (default data/tractor/flare_index_map.tsv)"
    "PRS_ANCESTRY_PROPS|EUR:0.35,AFR:0.15,AMR:0.50|global ancestry proportions weighting the combined PRS; labels must match data/tractor/ancestries.txt"
)
declare -A RR_SET_BY=()

# Refuse a launch that sets a switch name no script reads.
rr_refuse_unknown_switches() {
    local name
    for name in STRATUM_MIN_N PC_SCALE MRMEGA_PCS BUNDLE_RULES TIME_COL EVENT_COL; do
        if [[ -n "${!name+x}" ]]; then
            if [[ "$name" == STRATUM_MIN_N ]]; then
                echo "REFUSED: STRATUM_MIN_N is set, but the scripts read MIN_STRATUM_N; set that instead" >&2
            elif [[ "$name" == TIME_COL || "$name" == EVENT_COL ]]; then
                echo "REFUSED: ${name} is set, but each survival trait has its own columns now: set ${name}_<trait> (e.g. ${name}_efs) instead" >&2
            else
                echo "REFUSED: ${name} is set, but no script reads it (it would be ignored); unset it" >&2
            fi
            return 1
        fi
    done
}

# Fill in every unset switch with its default and export it, so every job of the launch sees the
# same value; remember which ones you set (RR_SET_BY). Call once, early, in the launcher.
rr_resolve_switches() {
    local entry name default
    for entry in "${RR_SWITCHES[@]}"; do
        IFS='|' read -r name default _ <<< "$entry"
        if [[ -n "${!name+x}" ]]; then
            RR_SET_BY[$name]=environment
            export "${name}"
        elif [[ "$default" == "@job" ]]; then
            RR_SET_BY[$name]=job
        else
            RR_SET_BY[$name]=default
            export "${name}=${default}"
        fi
    done
}

# params.tsv to stdout: param, value, set_by (environment | default | job).
rr_params_tsv() {
    local entry name value
    printf 'param\tvalue\tset_by\n'
    for entry in "${RR_SWITCHES[@]}"; do
        IFS='|' read -r name _ <<< "$entry"
        value="${!name-}"
        printf '%s\t%s\t%s\n' "$name" "${value//[$'\t\n']/ }" "${RR_SET_BY[$name]:-unresolved}"
    done
    # Derived rows the launcher adds (e.g. each survival trait's columns): "name<TAB>value<TAB>set_by"
    if [[ -n "${RR_EXTRA_PARAMS:-}" ]]; then printf '%s\n' "${RR_EXTRA_PARAMS}"; fi
}

_rr_git_sha()   { git -C "${_RR_REPO}" rev-parse HEAD 2>/dev/null || echo nogit; }
_rr_git_short() { git -C "${_RR_REPO}" rev-parse --short=7 HEAD 2>/dev/null || echo nogit; }

# rr_claim_run_dir ROOT [STAMP]: create ROOT/<run_id> and print the run_id. STAMP defaults to the
# current minute (yyyymmdd-HHMM).
rr_claim_run_dir() {
    local root="$1" stamp="${2:-$(date +%Y%m%d-%H%M)}" sha n
    sha="$(_rr_git_short)"
    mkdir -p "$root" || return 1
    for ((n = 1; n <= 999; n++)); do
        if mkdir "${root}/${stamp}_${sha}_${n}" 2>/dev/null; then
            echo "${stamp}_${sha}_${n}"
            return 0
        fi
    done
    echo "rr_claim_run_dir: no free run directory under ${root} for ${stamp}_${sha}" >&2
    return 1
}

# rr_preview_run_id ROOT: the run_id a launch now would get, without creating anything (dry run).
rr_preview_run_id() {
    local root="$1" stamp sha n
    stamp="$(date +%Y%m%d-%H%M)"
    sha="$(_rr_git_short)"
    for ((n = 1; n <= 999; n++)); do
        [[ -e "${root}/${stamp}_${sha}_${n}" ]] || { echo "${stamp}_${sha}_${n}"; return 0; }
    done
}

# rr_write_launch_record RUN_DIR LAUNCHER_COMMAND...: params.tsv, run_info.tsv, an empty jobs.tsv.
rr_write_launch_record() {
    local dir="$1"; shift
    local dirty=unknown
    if git -C "${_RR_REPO}" rev-parse HEAD >/dev/null 2>&1; then
        if [[ -n "$(git -C "${_RR_REPO}" status --porcelain --untracked-files=no 2>/dev/null)" ]]; then dirty=yes; else dirty=no; fi
    fi
    mkdir -p "${dir}/trace.d"
    rr_params_tsv > "${dir}/params.tsv"
    {
        printf 'key\tvalue\n'
        printf 'run_id\t%s\n' "$(basename "$dir")"
        printf 'git_commit\t%s\n' "$(_rr_git_sha HEAD)"
        printf 'uncommitted_changes\t%s\n' "$dirty"
        printf 'launched_at\t%s\n' "$(date -Is)"
        printf 'launch_host\t%s\n' "${HOSTNAME:-unknown}"
        printf 'launcher\t%s\n' "$*"
        printf 'partition\t%s\n' "${SLURM_PARTITION:-}"
        printf 'account\t%s\n' "${SLURM_ACCOUNT:-}"
        printf 'log_dir\t%s\n' "${LOG_DIR:-}"
        printf 'deferred\t%s\n' "${RR_DEFERRED:-}"
    } > "${dir}/run_info.tsv"
    printf 'label\tjob_id\n' > "${dir}/jobs.tsv"
}

# rr_note_job RUN_DIR LABEL JOB_ID: one line of jobs.tsv per submitted job (launcher only).
rr_note_job() { printf '%s\t%s\n' "$2" "$3" >> "$1/jobs.tsv"; }

# ---------------------------------------------------------------------------
# In a job
# ---------------------------------------------------------------------------
RR_COLUMNS="step	slurm_id	job_id	array_job_id	array_task_id	host	env	yml_sha256	script_sha256	args	inputs_md5	start	end	exit_status	n_null	n_tested	n_pass_mac	n_failed"

# The program's counts file, set with rr_sidecar: bin/tractor_genesis_adapter.R writes
# <output_prefix>.run_counts.tsv (key<TAB>value; keys in bin/tractor_io.R, write_run_counts).
# Prints n_null, n_tested, n_pass_mac, n_failed and the task status, NA for anything absent.
rr_sidecar() { RR_SIDECAR="$1"; rm -f -- "$1"; }   # removed now, so a stale file from an earlier run is never read
_rr_sidecar_counts() {
    local f="${RR_SIDECAR:-}"
    if [[ -n "$f" && -s "$f" ]]; then
        awk -F'\t' '{v[$1] = $2} END {
            n = split("n_null n_tested n_pass_mac n_failed status", k, " ")
            for (i = 1; i <= n; i++) printf "%s%s", (i > 1 ? "\t" : ""), ((k[i] in v) && v[k[i]] != "" ? v[k[i]] : "NA") }' "$f"
    else
        printf 'NA\tNA\tNA\tNA\tNA'
    fi
}

_rr_clean() { local v="$1"; v="${v//[$'\t\n\r']/ }"; printf '%s' "${v:-NA}"; }

# inputs_md5: "path=md5;path=md5", paths relative to the repo when inside it.
_rr_md5_list() {
    local f out="" name m
    for f in "$@"; do
        [[ -n "$f" ]] || continue
        name="${f#"${_RR_REPO}"/}"
        if [[ -f "$f" ]]; then
            m="$(md5sum "$f" 2>/dev/null)" || m="unreadable"
            out+="${name}=${m%% *};"
        else
            out+="${name}=missing;"
        fi
    done
    printf '%s' "${out%;}"
}

_rr_write_row() {   # END EXIT_STATUS
    local tmp="${RR_PART}.tmp.$$"
    {
        printf '%s\t' "$(_rr_clean "$RR_STEP")" "$RR_SLURM_ID" "$(_rr_clean "${SLURM_JOB_ID:-}")" \
            "$(_rr_clean "${SLURM_ARRAY_JOB_ID:-}")" "$(_rr_clean "${SLURM_ARRAY_TASK_ID:-}")" \
            "$(_rr_clean "${HOSTNAME:-}")" "$(_rr_clean "${GWAS_ENV:-none}")" "$(_rr_clean "${GWAS_ENV_YML_SHA256:-}")" \
            "$RR_SCRIPT_SHA" "$(_rr_clean "$RR_ARGS")" "$(_rr_clean "$RR_INPUTS")" "$RR_START" "$1" "$2"
        _rr_sidecar_counts
        printf '\n'
    } > "$tmp" 2>/dev/null && mv -f "$tmp" "$RR_PART" 2>/dev/null
    return 0
}

_rr_on_exit() {
    local rc=$?
    _rr_write_row "$(date -Is)" "$rc"
}

# rr_begin STEP "ARGS" [INPUT_FILE...]
rr_begin() {
    RR_STEP="${1:?rr_begin: step name}"; RR_ARGS="${2:-}"
    shift
    if [[ $# -gt 0 ]]; then shift; fi
    local dir
    if [[ -n "${GWAS_RUN_DIR:-}" ]]; then
        dir="$GWAS_RUN_DIR"
    else
        # Submitted by hand, not by the launcher: one directory per Slurm job (all tasks of an array share it)
        GWAS_RUN_ID="direct_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-local$$}}"
        dir="${_RR_REPO}/results/run_record/${GWAS_RUN_ID}"
        export GWAS_RUN_ID GWAS_RUN_DIR="$dir"
    fi
    mkdir -p "${dir}/trace.d"
    if [[ -n "${SLURM_ARRAY_JOB_ID:-}" && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
        RR_SLURM_ID="${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    else
        RR_SLURM_ID="${SLURM_JOB_ID:-local$$}"
    fi
    RR_PART="${dir}/trace.d/${RR_STEP//[^A-Za-z0-9_.-]/_}.${RR_SLURM_ID}.tsv"
    RR_SCRIPT_SHA="$(sha256sum "$0" 2>/dev/null)" || RR_SCRIPT_SHA="NA"   # $0: the script this job runs
    RR_SCRIPT_SHA="${RR_SCRIPT_SHA%% *}"
    RR_INPUTS="$(_rr_md5_list "$@")"
    RR_START="$(date -Is)"
    _rr_write_row "" ""
    trap '_rr_on_exit' EXIT
    trap 'exit 143' TERM      # Slurm's time limit and scancel send SIGTERM first
    trap 'exit 130' INT
    trap 'exit 129' HUP
    echo "run record: ${RR_PART}"
}

# rr_args "ARGS" [INPUT_FILE...]: fill in the args and inputs_md5 columns once they are known.
rr_args() {
    RR_ARGS="${1:-}"
    if [[ $# -gt 0 ]]; then shift; fi
    RR_INPUTS="$(_rr_md5_list "$@")"
    _rr_write_row "" ""
}

# ---------------------------------------------------------------------------
# Finalise: trace.tsv and summary.md for one run directory (bin/run_record_finalize.sh)
# ---------------------------------------------------------------------------
# sacct output to TSV keyed on JobID: State, ExitCode, MaxRSS, Elapsed (one line per job or step).
_rr_sacct() {   # IDS(comma-separated) OUTFILE
    local ids="$1" out="$2" try
    : > "$out"
    [[ -n "$ids" ]] || return 0
    if ! command -v sacct >/dev/null 2>&1; then
        echo "rr_finalize: sacct is not on PATH; State, MaxRSS and Elapsed are NA" >&2
        return 0
    fi
    for ((try = 1; try <= ${RR_SACCT_TRIES:-3}; try++)); do
        sacct -j "$ids" --noheader --parsable2 --format=JobID,State,ExitCode,MaxRSS,Elapsed > "$out" 2>/dev/null || true
        # Accounting can lag a job's end by a few seconds: retry while a finished job has no MaxRSS yet
        if ! awk -F'|' '$1 ~ /\.batch$/ && $4 == "" {missing = 1} END {exit !missing}' "$out"; then break; fi
        sleep "${RR_SACCT_WAIT:-30}"
    done
}

rr_finalize() {
    local dir="$1" ids parts=()
    [[ -d "${dir}/trace.d" ]] || { echo "rr_finalize: ${dir} is not a run record directory" >&2; return 1; }
    shopt -s nullglob
    parts=("${dir}"/trace.d/*.tsv)
    shopt -u nullglob
    if [[ ${#parts[@]} -gt 0 ]]; then cat "${parts[@]}" > "${dir}/.rows.tsv"; else : > "${dir}/.rows.tsv"; fi

    # every Slurm id: the tasks that wrote a row, and every job the launcher submitted
    ids="$( { awk -F'\t' '$2 !~ /^local/ {print $2}' "${dir}/.rows.tsv"
              [[ -f "${dir}/jobs.tsv" ]] && awk -F'\t' 'NR > 1 && $2 != "" {print $2}' "${dir}/jobs.tsv"; } | sort -u | paste -sd, - )"
    _rr_sacct "$ids" "${dir}/.sacct.txt"

    # trace.tsv: each task row + State, MaxRSS, Elapsed from sacct; exit status from sacct when the
    # task never wrote one (killed outright); a row for each submitted job that never started a task.
    # sacct lines: <id>|State|ExitCode|MaxRSS|Elapsed for the job, then <id>.batch, <id>.extern, ...
    # The job line carries State, ExitCode and Elapsed; MaxRSS is on the .batch step line.
    local files=("${dir}/.sacct.txt" "${dir}/.rows.tsv")
    [[ -f "${dir}/jobs.tsv" ]] && files+=("${dir}/jobs.tsv")
    awk -F'\t' -v OFS='\t' -v header="${RR_COLUMNS}	state	max_rss	elapsed	exit_from	task_status" '
        FILENAME ~ /\.sacct\.txt$/ {
            split($0, f, "|"); id = f[1]; base = id; sub(/\..*$/, "", base)
            if (id == base)          { state[id] = f[2]; code[id] = f[3]; el[id] = f[5] }
            else if (id ~ /\.batch$/) { if (f[4] != "") rss[base] = f[4] }
            else if (f[4] != "" && !(base in other)) other[base] = f[4]
            next
        }
        FILENAME ~ /\.rows\.tsv$/ { row[++n] = $0; seen[$2] = 1; seen[$4] = 1; next }
        FILENAME ~ /jobs\.tsv$/ && FNR > 1 { label[$2] = $1; next }
        function norow(lbl, id, st, ex, from, elv) {     # a task or job with no row of its own
            print lbl, id, id, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", ex, \
                  "NA", "NA", "NA", "NA", st, "NA", elv, from, "NA"
        }
        END {
            print header
            for (i = 1; i <= n; i++) {
                split(row[i], c, "\t"); id = c[2]; from = "trap"
                if (c[14] == "") {
                    from = "none"
                    if (id in code) { split(code[id], e, ":"); c[14] = (e[2] > 0) ? 128 + e[2] : e[1]; from = "sacct" }
                }
                line = c[1]; for (k = 2; k <= 18; k++) line = line OFS c[k]
                r = (id in rss) ? rss[id] : ((id in other) ? other[id] : "NA")
                print line, (id in state ? state[id] : "NA"), r, (id in el ? el[id] : "NA"), from, (c[19] == "" ? "NA" : c[19])
            }
            # Slurm knows of a task (or a whole job) that wrote no row: it never started (pending when
            # cancelled, e.g. after a failed dependency), or it died before rr_begin.
            for (k in state) {
                base = k; sub(/_.*$/, "", base)
                if (!(base in label) || (k in seen)) continue
                if (state[k] ~ /^(CANCELLED|PENDING)/) norow(label[base], k, state[k], "", "none", el[k])
                else { split(code[k], e, ":"); norow(label[base], k, state[k], (e[2] > 0) ? 128 + e[2] : e[1], "sacct", el[k]) }
                listed[base] = 1
            }
            for (j in label) if (!(j in seen) && !(j in listed)) norow(label[j], j, "NA", "", "none", "NA")
        }' "${files[@]}" > "${dir}/trace.tsv.tmp"
    mv -f "${dir}/trace.tsv.tmp" "${dir}/trace.tsv"
    rm -f "${dir}/.rows.tsv"
    _rr_summary "$dir" > "${dir}/summary.md.tmp" && mv -f "${dir}/summary.md.tmp" "${dir}/summary.md"
    echo "rr_finalize: ${dir}/trace.tsv ($(( $(wc -l < "${dir}/trace.tsv") - 1 )) rows), ${dir}/summary.md"
}

_rr_info() { [[ -f "$1/run_info.tsv" ]] && awk -F'\t' -v k="$2" '$1 == k {print $2; exit}' "$1/run_info.tsv"; }

_rr_summary() {
    local dir="$1"
    echo "# Run record $(basename "$dir")"
    echo ""
    echo "- launched: $(_rr_info "$dir" launched_at) on $(_rr_info "$dir" launch_host); git $(_rr_info "$dir" git_commit) (uncommitted changes: $(_rr_info "$dir" uncommitted_changes))"
    echo "- partition / account: $(_rr_info "$dir" partition) / $(_rr_info "$dir" account)"
    echo "- deferred: $(_rr_info "$dir" deferred)"
    echo "- finalised: $(date -Is)"
    if [[ -f "${dir}/params.tsv" ]]; then
        echo "- switches you set (every switch is in params.tsv): $(awk -F'\t' 'NR > 1 && $3 == "environment" {printf "%s%s=%s", sep, $1, $2; sep = ", "}' "${dir}/params.tsv")"
    fi
    echo ""
    echo "## Tasks by step"
    echo ""
    echo "| step | tasks | exit 0 | non-zero | no exit status | skipped, stratum below MIN_STRATUM_N | largest MaxRSS | longest Elapsed |"
    echo "|---|---|---|---|---|---|---|---|"
    awk -F'\t' '
        function mib(x,   v, u) { if (x == "" || x == "NA") return -1; v = x + 0; u = substr(x, length(x))
            if (u == "K") return v / 1024; if (u == "M") return v; if (u == "G") return v * 1024; if (u == "T") return v * 1048576; return v / 1048576 }
        function secs(x,   d, p, n) { if (x == "" || x == "NA") return -1; d = 0
            if (x ~ /-/) { split(x, p, "-"); d = p[1]; x = p[2] }
            n = split(x, p, ":"); if (n == 3) return d * 86400 + p[1] * 3600 + p[2] * 60 + p[3]; if (n == 2) return d * 86400 + p[1] * 60 + p[2]; return x + 0 }
        NR == 1 { next }
        { s = $1; if (!(s in t)) order[++k] = s; t[s]++
          if ($14 == "") none[s]++; else if ($14 == 0) ok[s]++; else bad[s]++
          if ($23 == "skipped_below_min_stratum_n") skip[s]++
          r = mib($20);  if (r >= 0 && (!(s in maxr) || r > maxr[s])) { maxr[s] = r; maxrs[s] = $20 }
          e = secs($21); if (e >= 0 && (!(s in maxe) || e > maxe[s])) { maxe[s] = e; maxes[s] = $21 } }
        END { for (i = 1; i <= k; i++) { s = order[i]
            printf "| %s | %d | %d | %d | %d | %d | %s | %s |\n", s, t[s], ok[s], bad[s], none[s], skip[s], (s in maxrs ? maxrs[s] : "NA"), (s in maxes ? maxes[s] : "NA") } }' "${dir}/trace.tsv"
    echo ""
    echo "## Tasks that did not exit 0"
    echo ""
    awk -F'\t' 'NR > 1 && $14 != "0" { if (!h) { print "| step | slurm id | exit status | state | host |"; print "|---|---|---|---|---|"; h = 1 }
                printf "| %s | %s | %s | %s | %s |\n", $1, $2, ($14 == "" ? "none" : $14), $19, $6 }
                END { if (!h) print "None." }' "${dir}/trace.tsv"
    echo ""
    echo "## Adapter counts"
    echo ""
    echo "For gwas rows, n_null (samples in the null model), n_tested (variants with a joint P),"
    echo "n_pass_mac (variants where an ancestry passed the MAC gate), n_failed (tests that raised an"
    echo "error) and task_status come from the adapter's <output_prefix>.run_counts.tsv; NA elsewhere."
    awk -F'\t' 'NR > 1 && $1 == "gwas" { n++; if ($23 == "skipped_below_min_stratum_n") k++; if ($15 ~ /^[0-9]+$/) nn += $15; if ($16 ~ /^[0-9]+$/) nt += $16; if ($18 ~ /^[0-9]+$/) nf += $18 }
                 END { if (n) printf "\ngwas: %d tasks, %d skipped below MIN_STRATUM_N; %d variants tested, %d failed tests.\n", n, k, nt, nf }' "${dir}/trace.tsv"
    echo ""
    echo "A job submitted later by a running job (the GxG heterogeneity job) is not waited for; run"
    echo "\`bin/run_record_finalize.sh ${dir}\` again once it has ended to add its row."
}
