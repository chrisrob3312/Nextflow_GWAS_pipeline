#!/usr/bin/env bash
# ============================================================================
# envs_tsv.sh: keep the yml_sha256 column of envs.tsv equal to envs/<env>.yml
# ============================================================================
#   bin/envs_tsv.sh           rewrite envs.tsv with the sha256 of each envs/<env>.yml
#   bin/envs_tsv.sh --check   change nothing; exit 1, naming the env, if any recorded sha256
#                             differs from its yml (bin/check_envs.sh runs this)
#
# Run it right after (re)building an env from its yml, so envs.tsv names the exact yml each env
# was built from. The run record copies that sha256 into every trace row.
# The equivalent by hand, for one env:  sha256sum envs/gwas-r.yml
# ============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TSV="${REPO}/envs.tsv"
MODE="${1:-update}"
[[ "$MODE" == update || "$MODE" == --check ]] || { echo "usage: $0 [--check]" >&2; exit 2; }
[[ "$(head -1 "$TSV")" == $'env\tprefix\tyml_sha256' ]] || { echo "envs.tsv: header must be: env prefix yml_sha256 (tab-separated)" >&2; exit 2; }

TMP="${TSV}.tmp.$$"
trap 'rm -f "$TMP"' EXIT
head -1 "$TSV" > "$TMP"
bad=0
while IFS=$'\t' read -r env prefix sha; do
    yml="${REPO}/envs/${env}.yml"
    [[ -f "$yml" ]] || { echo "envs.tsv: env '${env}' has no envs/${env}.yml" >&2; exit 2; }
    now="$(sha256sum "$yml" | cut -d' ' -f1)"
    if [[ "$now" != "$sha" ]]; then
        bad=1
        if [[ "$MODE" == --check ]]; then
            echo "envs.tsv: ${env} records yml_sha256 ${sha:-<empty>}, envs/${env}.yml is ${now}"
        else
            echo "envs.tsv: ${env} yml_sha256 ${sha:-<empty>} -> ${now}"
        fi
    fi
    printf '%s\t%s\t%s\n' "$env" "$prefix" "$now" >> "$TMP"
done < <(tail -n +2 "$TSV")

if [[ "$MODE" == --check ]]; then
    [[ $bad -eq 0 ]] && echo "envs.tsv matches envs/*.yml"
    exit "$bad"
fi
mv "$TMP" "$TSV"
echo "envs.tsv written ($(( $(wc -l < "$TSV") - 1 )) envs)"
