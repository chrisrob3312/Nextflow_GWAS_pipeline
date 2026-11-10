#!/bin/bash
#SBATCH --job-name=tractor_layout
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
# --partition, --account, --output and --error come from the sbatch command line
# (slurm/submit_full_pipeline.sh, or the direct-submit form in docs/HPC_SITE.md).

# ============================================================================
# The cohort's Tractor files -> each stratum's, between step 0 and the GWAS array (CP-16)
# ============================================================================
#   sbatch <site args> slurm/submit_tractor_layout.sh ancestries
#       one job: for each stratum, the share of each ancestry among its samples' haplotypes
#       (bin/la_tractor_layout.py shares, every 50th variant of every chromosome). The ancestries at
#       LA_MIN_SHARE or more (default 0.01) go into data/tractor/<S>/ancestries.txt, one per line. An
#       ancestries.txt that is already there (curated by hand) is kept, and a different suggestion
#       is printed next to it.
#   sbatch <site args> --array=0-<22 x strata - 1> slurm/submit_tractor_layout.sh split
#       one task per stratum x chromosome (task = stratum index x 22 + chromosome - 1): the
#       stratum's columns (results/strata/<S>.samples.txt) of data/tractor/cohort/chr<N>.* into
#       data/tractor/<S>/chr<N>.* for the ancestries in its ancestries.txt; with one ancestry, one
#       ancdose file (the summed hard-call dose) and no hapcount (bin/la_tractor_layout.py split).
# Strata: STRATA_LIST, else the bundle's strata.tsv. Labels: the bundle's flare_index_map.tsv.
# DRY_RUN=1 prints what the first task would run and stops.
# ============================================================================
set -euo pipefail
for SCRIPT_DIR in "${GWAS_REPO:-}" "${SLURM_SUBMIT_DIR:-}" "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; do
    [[ -f "${SCRIPT_DIR}/bin/slurm_utils.sh" ]] && break
done
[[ -f "${SCRIPT_DIR}/bin/slurm_utils.sh" ]] || { echo "ERROR: cannot find the repo (submit from the repo root, or export GWAS_REPO)" >&2; exit 1; }
source "${SCRIPT_DIR}/bin/slurm_utils.sh"
MODE="${1:?usage: submit_tractor_layout.sh ancestries|split}"
DRY_RUN="${DRY_RUN:-0}"
if [[ "${DRY_RUN}" != 1 ]]; then
    rr_begin "tractor_layout_${MODE}"
    source "${SCRIPT_DIR}/bin/activate_env.sh" lai
fi

HANDOFF_DIR="${HANDOFF_DIR:-${SCRIPT_DIR}/data/handoff}"
MAP="${HANDOFF_DIR}/local_ancestry/flare_index_map.tsv"
[[ -s "${MAP}" ]] || { echo "ERROR: ${MAP} is missing" >&2; exit 1; }
LABELS="$(tail -n +2 "${MAP}" | cut -f2 | paste -sd, -)"
STRATA_RUN="${STRATA_LIST:-$(bundle_strata)}" || exit 1
IFS=',' read -ra STRATA <<< "${STRATA_RUN}"
COH="${SCRIPT_DIR}/data/tractor/cohort"
LAYOUT="${SCRIPT_DIR}/bin/la_tractor_layout.py"
LA_MIN_SHARE="${LA_MIN_SHARE:-$(rr_default LA_MIN_SHARE)}"

case "${MODE}" in
ancestries)
    PREFIXES="$(ls "${COH}"/chr*.hapcount."${LABELS%%,*}".txt.gz 2>/dev/null | sed 's/\.hapcount\..*//' | paste -sd, -)"
    [[ -n "${PREFIXES}" ]] || { echo "ERROR: no ${COH}/chr<N>.hapcount.*.txt.gz (run slurm/submit_la_extract.sh first)" >&2; exit 1; }
    if [[ "${DRY_RUN}" == 1 ]]; then echo "[DRY-RUN] python3 ${LAYOUT} shares --cohort_prefixes ${PREFIXES} --labels ${LABELS} --samples <S>.samples.txt"; exit 0; fi
    for S in "${STRATA[@]}"; do
        SAMPLES="${SCRIPT_DIR}/results/strata/${S}.samples.txt"
        [[ -s "${SAMPLES}" ]] || { echo "ERROR: ${SAMPLES} is missing (step 0 writes it)" >&2; exit 1; }
        SH="$(python3 "${LAYOUT}" shares --cohort_prefixes "${PREFIXES}" --labels "${LABELS}" --samples "${SAMPLES}")"
        DERIVED="$(awk -F'\t' -v m="${LA_MIN_SHARE}" '$2 + 0 >= m + 0 {print $1}' <<< "${SH}" | sort)"
        mkdir -p "${SCRIPT_DIR}/data/tractor/${S}"
        F="${SCRIPT_DIR}/data/tractor/${S}/ancestries.txt"
        echo "Stratum ${S}: haplotype shares $(tr '\t\n' '= ' <<< "${SH}")-> $(paste -sd, <<< "${DERIVED}") (LA_MIN_SHARE ${LA_MIN_SHARE})"
        if [[ -s "$F" ]]; then
            [[ "$(sort "$F")" == "${DERIVED}" ]] || echo "  kept ${F} ($(paste -sd, "$F")), curated: the shares would give $(paste -sd, <<< "${DERIVED}")"
        else
            printf '%s\n' "${DERIVED}" > "$F"
        fi
    done
    ;;
split)
    IDX="${SLURM_ARRAY_TASK_ID:-0}"
    (( IDX < 22 * ${#STRATA[@]} )) || { echo "ERROR: task ${IDX} is past ${#STRATA[@]} strata x 22 chromosomes" >&2; exit 1; }
    S="${STRATA[$((IDX / 22))]}"; CHR=$(( IDX % 22 + 1 ))
    ANC="$(read_stratum_ancestries "${SCRIPT_DIR}/data/tractor/${S}/ancestries.txt")" || exit 1
    FOLD="$([[ ",${ANC}," == *",EUR,"* ]] && echo EUR || echo "${ANC%%,*}")"
    CMD=(python3 "${LAYOUT}" split --cohort_prefix "${COH}/chr${CHR}" --labels "${LABELS}" --keep "${ANC}"
         --samples "${SCRIPT_DIR}/results/strata/${S}.samples.txt" --out_prefix "${SCRIPT_DIR}/data/tractor/${S}/chr${CHR}" --fold_into "${FOLD}")
    if [[ "${DRY_RUN}" == 1 ]]; then echo "[DRY-RUN] $(print_cmd "${CMD[@]}")"; exit 0; fi
    rr_args "$(print_cmd "${CMD[@]}")" "${COH}/chr${CHR}".* "${SCRIPT_DIR}/results/strata/${S}.samples.txt"
    "${CMD[@]}"
    ;;
*) echo "ERROR: mode ${MODE}: ancestries or split" >&2; exit 2 ;;
esac
