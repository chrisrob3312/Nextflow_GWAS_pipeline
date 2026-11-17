#!/usr/bin/env bash
# ============================================================================
# mrmega_lead_snps.sh: the lead variants of one or more MR-MEGA results, for the study-cohort
# sensitivity arm (VARIANT_INCLUDE; docs/WORKFLOW.md, "Study-cohort sensitivity")
# ============================================================================
# usage: bin/mrmega_lead_snps.sh --output FILE [--p P] [--window_kb KB] [--all] RESULT [RESULT ...]
#   RESULT       an MR-MEGA output (<trait>.mrmega.result: MarkerName, Chromosome, Position,
#                P-value_association, ...; slurm/submit_mrmega.sh writes results/mrmega/<trait>/)
#   --p P        association P below which a variant counts (default LEAD_SNP_P: the environment,
#                else the settings file, else 5e-8; bin/run_record.sh)
#   --window_kb  one lead per locus: walking from the smallest P, a variant within KB kb of a lead
#                already kept on its chromosome is passed over (default 500)
#   --all        every variant below P, no per-locus thinning
# Writes FILE: one MarkerName (the Tractor ID chr:pos:ref:alt) per line, the union over the RESULTs,
# with a '#' line saying where it came from (the adapter's --variant_include skips '#' lines). Says
# how many it kept; exits 1 (and writes nothing) when no variant is below P, so a sensitivity arm is
# never queued on an empty list.
# ============================================================================
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/run_record.sh"   # rr_default, gs_effective
OUT=""; P=""; WIN=500; ALL=0; INPUTS=()
while [[ $# -gt 0 ]]; do case "$1" in
    --output) OUT="$2"; shift 2 ;;
    --p) P="$2"; shift 2 ;;
    --window_kb) WIN="$2"; shift 2 ;;
    --all) ALL=1; shift ;;
    -h|--help) sed -n '2,19p' "$0"; exit 0 ;;
    -*) echo "unknown option $1" >&2; exit 2 ;;
    *) INPUTS+=("$1"); shift ;;
esac; done
[[ -n "$OUT" && ${#INPUTS[@]} -gt 0 ]] || { echo "usage: $0 --output FILE [--p P] [--window_kb KB] [--all] RESULT [RESULT ...]" >&2; exit 2; }
[[ -n "$P" ]] || P="$(gs_effective LEAD_SNP_P)"
problem="$(rr_domain_problem LEAD_SNP_P "$P")"; [[ -z "$problem" ]] || { echo "ERROR: --p ${P}: ${problem}" >&2; exit 2; }
[[ "$WIN" =~ ^[0-9]+$ ]] || { echo "ERROR: --window_kb ${WIN}: a whole number of kb" >&2; exit 2; }
for f in "${INPUTS[@]}"; do [[ -s "$f" ]] || { echo "ERROR: ${f} is missing or empty" >&2; exit 1; }; done

# below P, from every input: MarkerName, chromosome, position, P (columns found by name)
TMP="$(mktemp)"; trap 'rm -f "$TMP"' EXIT
for f in "${INPUTS[@]}"; do
    awk -F'\t' -v p="$P" -v f="$f" '
        NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "MarkerName") m = i; if ($i == "Chromosome") c = i; if ($i == "Position") q = i; if ($i == "P-value_association") pc = i }
                  if (!m || !pc || !c || !q) { print "ERROR: " f " lacks MarkerName / Chromosome / Position / P-value_association" > "/dev/stderr"; exit 3 }; next }
        $pc != "NA" && $pc != "" && $pc + 0 < p + 0 { print $m "\t" $c "\t" $q "\t" $pc }' "$f" >> "$TMP"
done
N_BELOW="$(cut -f1 "$TMP" | sort -u | wc -l)"
if [[ "$N_BELOW" -eq 0 ]]; then
    echo "No variant below P = ${P} in: ${INPUTS[*]}; ${OUT} not written" >&2; exit 1
fi
mkdir -p "$(dirname "$OUT")"
{
    printf '# lead variants: MR-MEGA P-value_association < %s%s; from %s\n' "$P" "$([[ $ALL == 1 ]] && echo ", every variant" || echo ", one per ${WIN} kb locus")" "${INPUTS[*]}"
    sort -t$'\t' -k4,4g "$TMP" | awk -F'\t' -v w="$((WIN * 1000))" -v all="$ALL" '
        seen[$1]++ { next }
        all == 1 { print $1; next }
        { for (j = 1; j <= n[$2]; j++) { d = pos[$2, j] - $3; if (d < 0) d = -d; if (d <= w) next }
          n[$2]++; pos[$2, n[$2]] = $3; print $1 }'
} > "$OUT"
echo "Lead variants: $(grep -vc '^#' "$OUT") written to ${OUT} (${N_BELOW} variant(s) below P = ${P})"
