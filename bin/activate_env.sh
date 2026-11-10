# shellcheck shell=bash
# ============================================================================
# activate_env.sh: put one of the pipeline's conda envs on PATH, and TMPDIR on disk
# ============================================================================
# Usage, always with `source` (or its short form `.`), never run directly:
#   source bin/activate_env.sh gwas-r
#
# `source` runs this file inside the calling shell, so the activation and TMPDIR stay in effect
# for the rest of the job. Running it as `bash bin/activate_env.sh` would activate inside a child
# shell that ends at once, which achieves nothing, so that is refused.
#
# What it does, in order:
#   1. looks <env> up in envs.tsv (columns env, prefix, yml_sha256); a prefix written as
#      ${GWAS_ENVS_DIR}/<env> takes the site's env directory from the environment;
#   2. stops, naming the env and the path, if the prefix directory does not exist;
#   3. exports TMPDIR on disk, not on a tmpfs (M-04, docs/HPC_SITE.md "TMPDIR"): the first of
#      $SLURM_TMPDIR, $GWAS_SCRATCH/tmp, the current $TMPDIR (or /tmp) that is a writable
#      directory and not tmpfs/ramfs; it stops if none is;
#   4. runs `conda activate <prefix>` (loading conda's shell function first if the job shell
#      does not have it: `conda shell.bash hook`, or `module load $GWAS_CONDA_MODULE`);
#   5. exports GWAS_ENV, GWAS_ENV_PREFIX and GWAS_ENV_YML_SHA256 for the run record, and warns
#      if envs/<env>.yml no longer matches the sha256 recorded in envs.tsv.
#
# On failure it prints why and returns 1. A sourced file returns rather than exits, so the
# caller decides; a job script running under `set -e` stops right there.
# The caller's shell options are left as they were: set -e and set -u are off while this runs
# (conda's activation scripts are not written for them) and are switched back on afterwards if
# they were on before.
# ============================================================================

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "activate_env.sh must be sourced, not run: source bin/activate_env.sh <env>" >&2
    exit 2
fi

_gwas_activate_env() {
    local env_name="${1:-}"
    local repo tsv known row prefix yml_sha yml now token cand dir fstype tmp="" tmp_from="" rejected=""
    repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    tsv="${repo}/envs.tsv"

    # 1. look the env up
    if [[ ! -s "$tsv" ]]; then
        echo "activate_env: ${tsv} is missing or empty" >&2
        return 1
    fi
    known="$(awk -F'\t' 'NR > 1 {printf "%s%s", sep, $1; sep = " "}' "$tsv")"
    if [[ -z "$env_name" ]]; then
        echo "activate_env: no env named. Usage: source bin/activate_env.sh <env>   (envs: ${known})" >&2
        return 1
    fi
    row="$(awk -F'\t' -v e="$env_name" 'NR > 1 && $1 == e {print $2 "\t" $3; exit}' "$tsv")"
    if [[ -z "$row" ]]; then
        echo "activate_env: env '${env_name}' is not in ${tsv} (envs: ${known})" >&2
        return 1
    fi
    prefix="${row%%$'\t'*}"
    yml_sha="${row#*$'\t'}"
    # shellcheck disable=SC2016   # the literal text ${GWAS_ENVS_DIR}, as written in envs.tsv
    token='${GWAS_ENVS_DIR}'
    if [[ "$prefix" == *"$token"* ]]; then
        if [[ -z "${GWAS_ENVS_DIR:-}" ]]; then
            echo "activate_env: env '${env_name}': GWAS_ENVS_DIR is not set, and envs.tsv gives its prefix as ${prefix} (docs/HPC_SITE.md)" >&2
            return 1
        fi
        prefix="${prefix//"$token"/${GWAS_ENVS_DIR%/}}"
    fi

    # 2. the prefix must exist
    if [[ ! -d "$prefix" ]]; then
        echo "activate_env: env '${env_name}': prefix ${prefix} does not exist (build it: docs/HPC_SITE.md, \"Building the envs\")" >&2
        return 1
    fi

    # 3. TMPDIR on disk. data.table::fread() decompresses a .gz into TMPDIR before reading it;
    #    on a tmpfs that is RAM, and a chromosome-sized file fails with "No space left on device".
    for cand in SLURM_TMPDIR GWAS_SCRATCH TMPDIR; do
        case "$cand" in
            SLURM_TMPDIR) dir="${SLURM_TMPDIR:-}" ;;
            GWAS_SCRATCH) dir="${GWAS_SCRATCH:+${GWAS_SCRATCH%/}/tmp}"
                          [[ -n "$dir" ]] && mkdir -p "$dir" 2>/dev/null ;;
            TMPDIR)       dir="${TMPDIR:-/tmp}" ;;
        esac
        [[ -n "$dir" ]] || continue
        if [[ ! -d "$dir" || ! -w "$dir" ]]; then
            rejected="${rejected} ${cand}: ${dir} is not a writable directory;"
            continue
        fi
        fstype="$(stat -f -c %T "$dir" 2>/dev/null || echo unknown)"
        if [[ "$fstype" == tmpfs || "$fstype" == ramfs ]]; then
            rejected="${rejected} ${cand}: ${dir} is ${fstype};"
            continue
        fi
        tmp="$dir"; tmp_from="$cand"
        break
    done
    if [[ -z "$tmp" ]]; then
        echo "activate_env: no TMPDIR on disk:${rejected} set GWAS_SCRATCH to a scratch directory (docs/HPC_SITE.md, \"TMPDIR\")" >&2
        return 1
    fi
    if [[ "$tmp_from" == TMPDIR ]]; then
        echo "activate_env: note: neither SLURM_TMPDIR nor GWAS_SCRATCH is set; TMPDIR stays ${tmp} (${fstype})" >&2
    fi
    export TMPDIR="$tmp"

    # 4. conda activate. A batch job's shell usually lacks the `conda` shell function that
    #    `conda init` puts in an interactive login; `conda shell.bash hook` prints it.
    if ! declare -F conda >/dev/null 2>&1; then
        if ! command -v conda >/dev/null 2>&1 && [[ -n "${GWAS_CONDA_MODULE:-}" ]] && command -v module >/dev/null 2>&1; then
            module load "${GWAS_CONDA_MODULE}"
        fi
        if ! command -v conda >/dev/null 2>&1; then
            echo "activate_env: env '${env_name}': conda is not on PATH (load it before the job, or set GWAS_CONDA_MODULE; docs/HPC_SITE.md)" >&2
            return 1
        fi
        eval "$(conda shell.bash hook)" || { echo "activate_env: 'conda shell.bash hook' failed" >&2; return 1; }
    fi
    if ! conda activate "$prefix"; then
        echo "activate_env: env '${env_name}': conda activate ${prefix} failed" >&2
        return 1
    fi
    if [[ ! "${CONDA_PREFIX:-}" -ef "$prefix" ]]; then
        echo "activate_env: env '${env_name}': after conda activate, CONDA_PREFIX is '${CONDA_PREFIX:-}', not ${prefix}" >&2
        return 1
    fi

    # 5. what the run record reads
    export GWAS_ENV="$env_name" GWAS_ENV_PREFIX="$prefix" GWAS_ENV_YML_SHA256="$yml_sha"
    yml="${repo}/envs/${env_name}.yml"
    if [[ -f "$yml" ]]; then
        now="$(sha256sum "$yml" | cut -d' ' -f1)"
        if [[ "$now" != "$yml_sha" ]]; then
            echo "activate_env: WARNING envs/${env_name}.yml changed since envs.tsv was written (recorded ${yml_sha:0:12}, now ${now:0:12}); rebuild the env, then run bin/envs_tsv.sh" >&2
        fi
    fi
    echo "activate_env: ${env_name} (yml sha256 ${yml_sha:0:12}) at ${prefix}; TMPDIR=${TMPDIR} (${fstype}, from ${tmp_from})"
    return 0
}

_gwas_ae_opts="$-"            # the caller's one-letter shell options, e.g. "ehuB" under set -eu
set +eu
_gwas_activate_env "$@"
_gwas_ae_rc=$?
if [[ "$_gwas_ae_opts" == *e* ]]; then set -e; fi
if [[ "$_gwas_ae_opts" == *u* ]]; then set -u; fi
unset -f _gwas_activate_env
unset _gwas_ae_opts
return "$_gwas_ae_rc"
