#!/usr/bin/env bash
# ============================================================================
# la_flare_to_tractor.sh: a FLARE local-ancestry VCF -> the adapter's Tractor files (CP-16)
# ============================================================================
#   bin/la_flare_to_tractor.sh <chrN.anc.vcf.gz> <flare_index_map.tsv> <out_prefix> [extractor]
#
# 1. Extraction: Tractor's own extract_tracts_flare.py (Atkinson lab; on PATH in the lai env, placed
#    there as envs/lai.yml says), run with --num-ancs K, K = rows of the index map. It writes, per
#    FLARE ancestry index k, <vcf name>.anc<k>.hapcount.txt and .anc<k>.dosage.txt: the haplotypes of
#    that ancestry per sample (0-2) and the ALT alleles on them (0-2), from the phased GT that FLARE
#    called on. Both are needed (the regression has a local-ancestry term and a dose term).
# 2. Packaging (ours): index k -> label through the index map (the panel file's order of first
#    appearance: how FLARE numbers anc0..ancK-1, never alphabetical), written to the adapter's
#    pattern <out_prefix>.hapcount.<LABEL>.txt.gz / <out_prefix>.dosage.<LABEL>.txt.gz, header
#    CHROM POS ID REF ALT <IID...>.
# 3. Check (bin/la_tractor_check.py): every file there, hapcounts sum to 2, 0 <= dosage <= hapcount,
#    and the first 200 records recomputed from the VCF's GT and AN1 / AN2 agree.
# The colleague's own FLARE-to-Tractor transform is not in this repository, so the extraction is
# Tractor's script (the plan's fallback); only renaming and packaging are done here.
# ============================================================================
set -euo pipefail
VCF="${1:?usage: la_flare_to_tractor.sh <anc.vcf.gz> <flare_index_map.tsv> <out_prefix> [extractor]}"
MAP="${2:?index map required}"
OUT="${3:?output prefix required}"
EXTRACTOR="${4:-$(command -v extract_tracts_flare.py || true)}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[[ -s "$VCF" ]] || { echo "ERROR: ${VCF} is missing or empty" >&2; exit 1; }
[[ -s "$MAP" ]] || { echo "ERROR: ${MAP} is missing or empty" >&2; exit 1; }
[[ -n "$EXTRACTOR" && -f "$EXTRACTOR" ]] || { echo "ERROR: Tractor's extract_tracts_flare.py is not on PATH (the lai env; envs/lai.yml says where it comes from)" >&2; exit 1; }
head -1 "$MAP" | grep -qP '^anc_index\tlabel$' || { echo "ERROR: ${MAP}: header must be anc_index<TAB>label" >&2; exit 1; }
K=$(( $(wc -l < "$MAP") - 1 ))
LABELS=(); while IFS=$'\t' read -r k lab; do LABELS[$k]="$lab"; done < <(tail -n +2 "$MAP")
[[ ${#LABELS[@]} -eq $K ]] || { echo "ERROR: ${MAP}: anc_index must be 0..$((K - 1))" >&2; exit 1; }

WORK="$(mktemp -d "${TMPDIR:-/tmp}/la_extract.XXXXXX")"
trap 'rm -rf "${WORK:?}"' EXIT
echo "Extracting ${VCF} (K=${K}: ${LABELS[*]}) with ${EXTRACTOR}"
python3 "$EXTRACTOR" --vcf "$VCF" --num-ancs "$K" --output-dir "$WORK"

base="$(basename "$VCF")"; base="${base%.gz}"; base="${base%.vcf}"
mkdir -p "$(dirname "$OUT")"
for k in $(seq 0 $((K - 1))); do
    for kind in hapcount dosage; do
        src=""
        for c in "${WORK}/${base}.anc${k}.${kind}.txt" "${WORK}/${base}.anc${k}.${kind}.txt.gz"; do [[ -s "$c" ]] && src="$c" && break; done
        [[ -n "$src" ]] || { echo "ERROR: the extraction wrote no ${kind} file for ancestry index ${k} (${LABELS[$k]}) in ${WORK}" >&2; ls "$WORK" >&2; exit 1; }
        dst="${OUT}.${kind}.${LABELS[$k]}.txt.gz"
        if [[ "$src" == *.gz ]]; then cp "$src" "$dst"; else gzip -c "$src" > "$dst"; fi
        echo "  anc${k} ${kind} -> $(basename "$dst")"
    done
done
python3 "${HERE}/la_tractor_check.py" --prefix "$OUT" --labels "$(IFS=,; echo "${LABELS[*]}")" --vcf "$VCF" --index_map "$MAP" --rows 200
