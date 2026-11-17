# shellcheck shell=bash
# ============================================================================
# gwas_settings.sh: the user-editable settings file (config/gwas_settings.txt)
# ============================================================================
# Sourced, never run: by bin/run_record.sh (so the launcher, every job and bin/check_handoff_bundle.sh
# read the same file the same way).
#
# Which file: SETTINGS_FILE=<path> (relative to the directory you launch from), else
# <repo>/config/gwas_settings.txt when it exists, else none. SETTINGS_FILE=none reads no file
# (the built-in defaults of RR_SWITCHES only).
#
# Format (the harmonization pipeline's pipeline_settings.txt style):
#   one  name = value  per line; '#' starts a comment; blank lines are ignored
#   name is a switch of RR_SWITCHES (bin/run_record.sh), or one of
#     trait.<name> = binary | quantitative | survival[:<time column>:<event column>]
#                    declares a trait and its model; TRAITS_LIST defaults to every declared trait,
#                    in file order; a survival trait without columns reads <name>_time / <name>_event
#     TIME_COL_<trait> / EVENT_COL_<trait> = <column>   one survival trait's columns (the older form)
#     arm.<label> = SWITCH=value [SWITCH=value ...]
#                    a named set of switches, applied by ARM=<label> (e.g. an MRD-positive arm)
#   an empty value sets the switch to empty (COX_STRATA = means no stratification); leave a line
#   out, or comment it, for the built-in default
# Every value is checked against its switch's domain here, so a refusal names the line. An unknown
# name, a name given twice and a malformed line are refused, naming the line: a setting is never
# silently ignored.
#
# Precedence, per switch: built-in default < settings file < arm.<label> (when ARM=<label>) <
# environment. params.tsv records which one won (set_by default / settings / arm:<label> / environment).
# ============================================================================
if [[ -n "${_GS_LIB_LOADED:-}" ]]; then return 0; fi
_GS_LIB_LOADED=1

declare -A GS_VALUE=()     # name -> value (settings lines, plus what trait.<name> lines derive)
declare -A GS_LINE=()      # name -> line number in the file
declare -A GS_ARM=()       # arm label -> "SWITCH=value SWITCH=value ..."
declare -A GS_ARM_LINE=()  # arm label -> line number
GS_ARMS=()                 # arm labels in file order
GS_TRAITS=()               # declared traits (trait.<name> lines) in file order
GS_FILE=""                 # the file read; empty: none
GS_STATE=unread            # unread | loaded | error

# gs_path: the settings file this shell reads, or nothing
gs_path() {
    local repo; repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    if [[ "${SETTINGS_FILE:-}" == none ]]; then return 0
    elif [[ -n "${SETTINGS_FILE:-}" ]]; then printf '%s' "${SETTINGS_FILE}"
    elif [[ -f "${repo}/config/gwas_settings.txt" ]]; then printf '%s' "${repo}/config/gwas_settings.txt"; fi
}

_gs_is_switch() {   # NAME: is it a switch of RR_SWITCHES?
    local entry
    for entry in "${RR_SWITCHES[@]}"; do [[ "${entry%%|*}" == "$1" ]] && return 0; done
    return 1
}
_gs_trim() { local v="$1"; v="${v#"${v%%[![:space:]]*}"}"; v="${v%"${v##*[![:space:]]}"}"; printf '%s' "$v"; }

# gs_load: read and check the settings file once per shell. Returns 1 (REFUSED lines on stderr)
# when anything in it is wrong; then nothing from it is used.
gs_load() {
    case "$GS_STATE" in loaded) return 0 ;; error) return 1 ;; esac
    local f; f="$(gs_path)"
    if [[ -z "$f" ]]; then GS_STATE=loaded; return 0; fi
    if [[ ! -f "$f" ]]; then echo "REFUSED: SETTINGS_FILE=${f} does not exist" >&2; GS_STATE=error; return 1; fi
    GS_FILE="$f"
    local n=0 raw line key value bad=0 problem tok tk tv name type rest declared_types=()
    local -a parts
    _gs_bad() { echo "REFUSED: ${GS_FILE} line ${n}: $*" >&2; bad=1; }
    while IFS= read -r raw || [[ -n "$raw" ]]; do
        n=$((n + 1))
        line="${raw%%#*}"; line="${line//$'\r'/}"
        [[ "$line" =~ ^[[:space:]]*$ ]] && continue
        if [[ ! "$line" =~ ^[[:space:]]*([A-Za-z_][A-Za-z0-9_.]*)[[:space:]]*=(.*)$ ]]; then
            _gs_bad "expected 'name = value', found '$(_gs_trim "$line")'"; continue
        fi
        key="${BASH_REMATCH[1]}"; value="$(_gs_trim "${BASH_REMATCH[2]}")"
        if [[ -n "${GS_LINE[$key]+x}" ]]; then _gs_bad "${key} is already set on line ${GS_LINE[$key]}"; continue; fi
        case "$key" in
            trait.*)
                name="${key#trait.}"
                [[ "$name" =~ ^[A-Za-z][A-Za-z0-9_]*$ ]] || { _gs_bad "trait name '${name}' (letters, digits and _ only)"; continue; }
                IFS=: read -ra parts <<< "$value"; type="${parts[0]:-}"
                case "$type" in
                    binary|quantitative) [[ ${#parts[@]} -eq 1 ]] || { _gs_bad "${key} = ${value}: only survival takes columns"; continue; } ;;
                    survival)
                        if [[ ${#parts[@]} -eq 3 ]]; then
                            for tok in "${parts[1]}" "${parts[2]}"; do
                                [[ "$tok" =~ ^[A-Za-z_][A-Za-z0-9_.]*$ ]] || { _gs_bad "${key} = ${value}: '${tok}' is not a column name"; continue 2; }
                            done
                            GS_VALUE[TIME_COL_${name}]="${parts[1]}"; GS_LINE[TIME_COL_${name}]=$n
                            GS_VALUE[EVENT_COL_${name}]="${parts[2]}"; GS_LINE[EVENT_COL_${name}]=$n
                        elif [[ ${#parts[@]} -ne 1 ]]; then
                            _gs_bad "${key} = ${value}: use survival, or survival:<time column>:<event column>"; continue
                        fi ;;
                    *) _gs_bad "${key} = ${value}: use binary, quantitative or survival[:<time column>:<event column>]"; continue ;;
                esac
                GS_TRAITS+=("$name"); declared_types+=("$type") ;;
            arm.*)
                name="${key#arm.}"
                [[ "$name" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ && "$name" != all ]] \
                    || { _gs_bad "arm label '${name}' (letters, digits, _ . -, one folder name; 'all' is reserved)"; continue; }
                [[ -n "$value" ]] || { _gs_bad "${key} sets no switch (use SWITCH=value ...)"; continue; }
                for tok in $value; do
                    [[ "$tok" == *=* ]] || { _gs_bad "${key}: '${tok}' is not SWITCH=value"; continue 2; }
                    tk="${tok%%=*}"; tv="${tok#*=}"
                    case "$tk" in ARM) _gs_bad "${key}: an arm cannot set ARM"; continue 2 ;; esac
                    _gs_is_switch "$tk" || { _gs_bad "${key}: unknown switch '${tk}' (the switches are listed in bin/run_record.sh, RR_SWITCHES)"; continue 2; }
                    problem="$(rr_domain_problem "$tk" "$tv")"
                    [[ -z "$problem" ]] || { _gs_bad "${key}: ${tk}=${tv}; ${problem}"; continue 2; }
                done
                GS_ARM[$name]="$value"; GS_ARM_LINE[$name]=$n; GS_ARMS+=("$name") ;;
            TIME_COL_*|EVENT_COL_*)
                [[ "$value" =~ ^[A-Za-z_][A-Za-z0-9_.]*$ ]] || { _gs_bad "${key} = ${value}: not a column name"; continue; } ;;
            *)
                _gs_is_switch "$key" || { _gs_bad "unknown setting '${key}' (a switch of bin/run_record.sh RR_SWITCHES, trait.<name>, arm.<label>, TIME_COL_<trait> or EVENT_COL_<trait>)"; continue; }
                problem="$(rr_domain_problem "$key" "$value")"
                [[ -z "$problem" ]] || { _gs_bad "${key} = ${value}; ${problem}"; continue; } ;;
        esac
        GS_VALUE[$key]="$value"; GS_LINE[$key]=$n
    done < "$f"
    # trait.<name> lines are the trait list and the model lists; the list switches may not also be set
    if [[ ${#GS_TRAITS[@]} -gt 0 && $bad -eq 0 ]]; then
        local i surv="" quant="" first_line="${GS_LINE[trait.${GS_TRAITS[0]}]}" t
        for key in SURVIVAL_TRAITS QUANTITATIVE_TRAITS; do
            [[ -z "${GS_LINE[$key]+x}" ]] || { n="${GS_LINE[$key]}"; _gs_bad "${key} and trait.<name> lines both give the trait models; keep the trait.<name> lines"; }
        done
        for i in "${!GS_TRAITS[@]}"; do
            case "${declared_types[$i]}" in
                survival) surv+="${surv:+,}${GS_TRAITS[$i]}" ;;
                quantitative) quant+="${quant:+,}${GS_TRAITS[$i]}" ;;
            esac
        done
        GS_VALUE[SURVIVAL_TRAITS]="$surv"; GS_LINE[SURVIVAL_TRAITS]=$first_line
        GS_VALUE[QUANTITATIVE_TRAITS]="$quant"; GS_LINE[QUANTITATIVE_TRAITS]=$first_line
        if [[ -n "${GS_LINE[TRAITS_LIST]+x}" ]]; then
            IFS=, read -ra parts <<< "${GS_VALUE[TRAITS_LIST]}"
            for t in "${parts[@]}"; do
                [[ " ${GS_TRAITS[*]} " == *" ${t} "* ]] || { n="${GS_LINE[TRAITS_LIST]}"; _gs_bad "TRAITS_LIST names ${t}, which no trait.${t} line declares"; }
            done
        else
            GS_VALUE[TRAITS_LIST]="$(IFS=,; echo "${GS_TRAITS[*]}")"; GS_LINE[TRAITS_LIST]=$first_line
        fi
    fi
    unset -f _gs_bad
    if [[ $bad -ne 0 ]]; then GS_STATE=error; return 1; fi
    GS_STATE=loaded
}

# gs_value NAME: the settings file's value for NAME (return 1 when the file does not set it)
gs_value() {
    gs_load 2>/dev/null || return 1
    [[ -n "${GS_LINE[$1]+x}" ]] || return 1
    printf '%s' "${GS_VALUE[$1]}"
}

# gs_arm_value ARM NAME: the value arm.<ARM> gives NAME (return 1 when it gives none)
gs_arm_value() {
    local tok
    gs_load 2>/dev/null || return 1
    [[ -n "${GS_ARM[$1]+x}" ]] || return 1
    for tok in ${GS_ARM[$1]}; do
        [[ "${tok%%=*}" == "$2" ]] && { printf '%s' "${tok#*=}"; return 0; }
    done
    return 1
}

# gs_effective NAME: the value a launch would use: environment, else the arm (ARM), else the
# settings file, else the built-in default (for scripts outside the launcher, e.g. the bundle checker)
gs_effective() {
    if [[ -n "${!1+x}" ]]; then printf '%s' "${!1}"; else rr_default "$1"; fi
}

# gs_where NAME: " (<file> line N)" when the value of NAME came from the settings file, else nothing
gs_where() {
    gs_load 2>/dev/null || return 0
    if [[ -n "${GS_LINE[$1]+x}" ]]; then printf ' (%s line %s)' "${GS_FILE}" "${GS_LINE[$1]}"; fi
}
