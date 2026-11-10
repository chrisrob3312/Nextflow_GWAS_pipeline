#!/usr/bin/env bash
# ============================================================================
# qc_gate.sh: stop unless every required QC assertion is there and did not fail (CP-11)
# ============================================================================
#   bin/qc_gate.sh --require id1[,id2,...] [--dir DIR]
#   bin/qc_gate.sh --collect [--dir DIR]     rebuild DIR/qc_assertions.tsv and DIR/sample_accounting.tsv
#                                            from the per-row files (the run record's finalise job)
# Reads DIR/qc_assertions.tsv (header id verdict n message) and the per-row files under
# DIR/assertions.d/ that bin/qc_assert.R writes. DIR defaults
# to $GWAS_QC_DIR, else results/qc.
#
# For each required id:
#   no row at all            -> "missing assertion <id>": the check did not run, which is not a pass
#   any row with FAIL        -> "assertion <id> failed: <message>"
#   rows PASS or INFO only   -> fine
# Exit 0 only if every required id is fine. QC_GATE_OVERRIDE=<id>[,<id>...] lets a named failed or
# missing assertion through, and records an OVERRIDE row (who, when, what it overrode) next to
# it. There is no switch that skips every check.
# ============================================================================
set -uo pipefail
REQ=""; DIR="${GWAS_QC_DIR:-results/qc}"; COLLECT=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --require) REQ="${2:-}"; shift 2 ;;
        --dir)     DIR="${2:-}"; shift 2 ;;
        --collect) COLLECT=1; shift ;;
        *) echo "usage: qc_gate.sh --require id1[,id2,...] [--dir DIR] | --collect [--dir DIR]" >&2; exit 2 ;;
    esac
done

# rows DIR SUB: every per-row file under DIR/SUB, sorted, one row per line. find -exec cat {} +
# concatenates the files; sort -u drops a row seen twice.
rows() { [[ -d "$1/$2" ]] && find "$1/$2" -name '*.tsv' -type f -exec cat {} + | sort -u; }

if [[ $COLLECT == 1 ]]; then
    [[ -d "$DIR" ]] || { echo "qc_gate.sh --collect: no QC directory ${DIR}"; exit 0; }
    { printf 'id\tverdict\tn\tmessage\n'; rows "$DIR" assertions.d; } > "${DIR}/qc_assertions.tsv"
    { printf 'step\tstratum\ttrait\tchromosome\tn_in\tn_out\tn_removed\treason\n'; rows "$DIR" accounting.d; } \
        > "${DIR}/sample_accounting.tsv"
    echo "QC: $(($(wc -l < "${DIR}/qc_assertions.tsv") - 1)) assertion row(s), $(($(wc -l < "${DIR}/sample_accounting.tsv") - 1)) sample-accounting row(s) in ${DIR}"
    exit 0
fi
[[ -n "$REQ" ]] || { echo "qc_gate.sh: --require is empty (a gate must name what it requires)" >&2; exit 2; }

# qc_assertions.tsv (if collected already) plus the per-row files, each row once
ROWS="$( { [[ -f "${DIR}/qc_assertions.tsv" ]] && tail -n +2 "${DIR}/qc_assertions.tsv"; rows "$DIR" assertions.d; } 2>/dev/null | sort -u)"

bad=0; overridden=""; IFS=',' read -ra IDS <<< "$REQ"
for id in "${IDS[@]}"; do
    rows="$(awk -F'\t' -v id="$id" '$1 == id' <<< "$ROWS")"
    if [[ -z "$rows" ]]; then
        problem="missing assertion ${id} (the check did not run)"
    elif awk -F'\t' '$2 == "FAIL" {f = 1} END {exit !f}' <<< "$rows"; then
        problem="assertion ${id} failed: $(awk -F'\t' '$2 == "FAIL" {print $4; exit}' <<< "$rows")"
    else
        echo "QC GATE ok    ${id} ($(wc -l <<< "$rows") row(s))"
        continue
    fi
    if [[ ",${QC_GATE_OVERRIDE:-}," == *",${id},"* ]]; then
        mkdir -p "${DIR}/assertions.d"
        printf '%s\tOVERRIDE\tNA\tQC_GATE_OVERRIDE by %s at %s; was: %s\n' "$id" "${USER:-unknown}" "$(date -Is)" "$problem" \
            > "${DIR}/assertions.d/${id}.override.$(date +%s).$$.tsv"
        echo "QC GATE OVERRIDE ${problem} (recorded)"
        overridden="${overridden:+${overridden},}${id}"
        continue
    fi
    echo "QC GATE: ${problem}" >&2
    bad=1
done
[[ $bad -eq 0 ]] && echo "QC GATE: passed (${REQ})${overridden:+; overridden: ${overridden}}"
exit "$bad"
