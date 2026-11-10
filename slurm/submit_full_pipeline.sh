#!/bin/bash
# ============================================================================
# Master SLURM Submission Script
# ============================================================================
# Submits complete multi-ancestry GWAS pipeline with proper dependencies
#
# Workflow:
#   0. Structure: ingest the bundle's PCs and GRM (STRUCTURE_SOURCE=bundle, default), or
#      compute them with PC-AiR / PC-Relate (=compute); skip with SKIP_PCAIR=true
#   1. GWAS: Run Tractor-GENESIS for each ancestry stratum (22 chr × 5 strata)
#   2. Combine: Merge chromosomes, run meta-analysis for OTHER
#   3. PRS: Run 6 methods × 3 traits                         only with RUN_PRS=true
#   4. Colocalization per trait, all QTL types                 only with RUN_COLOC=true
#   5. Validation: Compare methods overall and by ancestry    only with RUN_VALIDATION=true (needs RUN_PRS)
#   6. GxG interaction on the GWAS hits                       only with RUN_GXG=true
# Steps 3-6 are off by default; every launch prints which were deferred.
#
# Usage (from any directory):
#   export SLURM_PARTITION=<partition> SLURM_ACCOUNT=<account>   # required, docs/HPC_SITE.md
#   slurm/submit_full_pipeline.sh [--dry-run]
#
# Partition and account are put on every sbatch command line by submit() below (an #SBATCH
# line cannot read a variable), and every job log is an absolute path under LOG_DIR
# (default <repo>/logs). To submit one child script by hand, see docs/HPC_SITE.md.
#
# Run record (bin/run_record.sh): every launch gets results/run_record/<run_id>/ with
# params.tsv (every number-changing switch) written now, one trace row per task written by
# the task, and trace.tsv + summary.md written by a finalise job that waits on everything.
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Work from the repo root whatever directory this was launched from. sbatch records the
# directory it was called from as each job's SLURM_SUBMIT_DIR and starting directory.
cd "${SCRIPT_DIR}/.."
REPO="$PWD"
export GWAS_REPO="$REPO"     # how every job finds the repo (Slurm runs a job script from a spool copy)

case "${1:-}" in
    "")        DRY_RUN=0 ;;
    --dry-run) DRY_RUN=1 ;;
    *)         echo "usage: $0 [--dry-run]" >&2; exit 2 ;;
esac

# Site settings come from the environment, never from this file (docs/HPC_SITE.md)
for v in SLURM_PARTITION SLURM_ACCOUNT; do
    if [[ -z "${!v:-}" ]]; then
        echo "REFUSED: ${v} is not set. Export SLURM_PARTITION and SLURM_ACCOUNT before launching (docs/HPC_SITE.md)." >&2
        exit 2
    fi
done
LOG_DIR="${LOG_DIR:-${REPO}/logs}"
[[ "$LOG_DIR" == /* ]] || LOG_DIR="${REPO}/${LOG_DIR}"     # a relative LOG_DIR is taken from the repo root
export LOG_DIR
[[ $DRY_RUN == 1 ]] || mkdir -p "${LOG_DIR}"

# ---------------------------------------------------------------------------
# Customisable via environment variables (defaults for the leukemia cohort). The defaults of
# every number-changing switch live in one table, RR_SWITCHES in bin/run_record.sh;
# rr_resolve_switches fills in the unset ones and exports them, so every job sees the value
# that params.tsv records.
#   STRATA_LIST="EUR,AFR"                   strata analysed; default: every stratum in the bundle's
#                                           ancestry/strata.tsv (strata are assigned upstream)
#   TRAITS_LIST="mrd,efs,os"               phenotype columns, as the handoff contract names them
#   SURVIVAL_TRAITS="efs,os"                traits run as time-to-event, each with its own columns
#                                           <trait>_time and <trait>_event (TIME_COL_<trait> /
#                                           EVENT_COL_<trait> override one trait); the rest are binary
#   COVARIATES="sex,age_at_diagnosis,PC1"  numeric phenotype columns, used by GWAS and GxG; default:
#                                           the bundle's phenotypes/covariates.txt, then PC1..PC5
#   GXG_CUSTOM_VARIANTS=/path/list.txt      optional GRCh38 variant list for GxG
#   GXG_USE_DEFAULT_LOCI=true|false         include the default known-leukemia-loci table
#   LOG_DIR=/abs/path                       where every job log goes (default <repo>/logs)
#   RUN_PRS / RUN_COLOC / RUN_VALIDATION / RUN_GXG = true|false   steps 3-6 (default false)
#   MIN_STRATUM_N=50                        smallest stratum analysed (decided 10 Nov 2026)
#   MAC_MIN=10                              minor allele count per ancestry term; open question: the
#                                           handoff contract says 20, the code default is still 10
# ---------------------------------------------------------------------------
source "${REPO}/bin/slurm_utils.sh"    # helpers, and bin/run_record.sh (the switch table, rr_*)
rr_refuse_unknown_switches || exit 2
rr_resolve_switches
# Every switch inside its domain (B-05; RR_SWITCHES in bin/run_record.sh), before anything else
rr_check_domains || exit 2
# RUN_LABEL (B-05): a labelled run (a sensitivity arm) writes its results under results/<RUN_LABEL>/
if [[ -n "${RUN_LABEL:-}" ]]; then
    [[ "${RUN_LABEL}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || { echo "REFUSED: RUN_LABEL=${RUN_LABEL}; use letters, digits, _ . - (one folder name)" >&2; exit 2; }
    case "${RUN_LABEL}" in
        gwas|prs|coloc|gxg|pcair|strata|run_record|qc|mrmega)
            echo "REFUSED: RUN_LABEL=${RUN_LABEL} is the name of a results folder; choose another" >&2; exit 2 ;;
    esac
fi
export GWAS_RESULTS="${REPO}/results${RUN_LABEL:+/${RUN_LABEL}}"
export HANDOFF_DIR="${HANDOFF_DIR:-${REPO}/data/handoff}"   # the handoff bundle (docs/HPC_SITE.md)

# Strata (CP-22): assigned upstream, listed in the bundle's ancestry/strata.tsv. STRATA_LIST unset
# means all of them (set_by bundle); a stratum strata.tsv does not have is refused.
STRATA_TSV="${HANDOFF_DIR}/ancestry/strata.tsv"
[[ -s "${STRATA_TSV}" ]] || { echo "REFUSED: ${STRATA_TSV} is missing or empty (the strata come from the handoff bundle)" >&2; exit 2; }
BUNDLE_STRATA="$(bundle_strata)"
if [[ "${RR_SET_BY[STRATA_LIST]}" == job ]]; then
    # MIX_RETAIN=false leaves the bundle's MIX stratum out (B-05)
    STRATA_LIST="$(tr ',' '\n' <<< "${BUNDLE_STRATA}" | { if [[ "${MIX_RETAIN}" == false ]]; then grep -vx MIX; else cat; fi; } | paste -sd, -)"
    export STRATA_LIST; RR_SET_BY[STRATA_LIST]=bundle
fi
rr_check_domains STRATA_LIST || exit 2
IFS=',' read -ra STRATA <<< "${STRATA_LIST}"
for stratum in "${STRATA[@]}"; do
    [[ ",${BUNDLE_STRATA}," == *",${stratum},"* ]] \
        || { echo "REFUSED: stratum ${stratum} (STRATA_LIST) is not in ${STRATA_TSV} (it has: ${BUNDLE_STRATA})" >&2; exit 2; }
    [[ "${MIX_RETAIN}" == false && "${stratum}" == MIX ]] \
        && { echo "REFUSED: STRATA_LIST names MIX, but MIX_RETAIN=false" >&2; exit 2; }
done
IFS=',' read -ra TRAITS <<< "${TRAITS_LIST}"

# Steps 2b and 3-6 run only when asked for (their RUN_* are true / false: rr_check_domains above)
if [[ "${RUN_VALIDATION}" == true && "${RUN_PRS}" != true ]]; then
    echo "REFUSED: RUN_VALIDATION=true needs RUN_PRS=true (validation summarises the PRS jobs of this launch)" >&2
    exit 2
fi
# Each stratum's Tractor ancestry list (CP-28) must exist before its 22-task arrays are queued:
# submit_gwas_array.sh reads it in every task and refuses without it. A dry run lists what is
# missing; a launch refuses.
MISSING_ANC=()
for stratum in "${STRATA[@]}"; do
    [[ -s "${REPO}/data/tractor/${stratum}/ancestries.txt" ]] || MISSING_ANC+=("${REPO}/data/tractor/${stratum}/ancestries.txt")
done
# PRS, and GxG's pooled task when there are Tractor files, read the cohort-wide list instead
if [[ "${RUN_PRS}" == true || ( "${RUN_GXG}" == true && -d "${REPO}/data/tractor" ) ]]; then
    [[ -s "${REPO}/data/tractor/ancestries.txt" ]] || MISSING_ANC+=("${REPO}/data/tractor/ancestries.txt")
fi
if [[ ${#MISSING_ANC[@]} -gt 0 ]]; then
    if [[ $DRY_RUN == 1 ]]; then
        printf 'WARNING: missing %s (a launch would refuse)\n' "${MISSING_ANC[@]}"
    else
        printf 'REFUSED: missing %s (one per stratum; docs/HANDOFF_CONTRACT.md)\n' "${MISSING_ANC[@]}" >&2
        exit 2
    fi
fi

# Every trait's columns must be in the phenotype file: a binary trait its own column, a survival
# trait <trait>_time and <trait>_event (or its TIME_COL_<trait> / EVENT_COL_<trait>). Checked
# here, on the header, before anything is queued. A dry run warns; a launch refuses.
PHENO_SOURCE="${PHENO_IN:-${HANDOFF_DIR}/phenotypes/phenotypes.tsv}"
PHENO_HEADER=""
[[ -s "${PHENO_SOURCE}" ]] && PHENO_HEADER="$(head -1 "${PHENO_SOURCE}" | tr '\t' '\n')"
SURV_COLS=()
PHENO_PROBLEMS=()
for trait in "${TRAITS[@]}"; do
    [[ "$trait" =~ ^[A-Za-z][A-Za-z0-9_]*$ ]] || { echo "REFUSED: trait name '${trait}' (letters, digits and _ only)" >&2; exit 2; }
    if trait_is_survival "$trait"; then
        need=("$(trait_time_col "$trait")" "$(trait_event_col "$trait")")
        SURV_COLS+=("${trait}:${need[0]}/${need[1]}")
    else
        need=("$trait")
    fi
    if [[ -s "${PHENO_SOURCE}" ]]; then
        for c in "${need[@]}"; do
            grep -qxF -- "$c" <<< "${PHENO_HEADER}" || PHENO_PROBLEMS+=("${trait}: column ${c} is not in ${PHENO_SOURCE}")
        done
    fi
done
[[ -s "${PHENO_SOURCE}" ]] || PHENO_PROBLEMS+=("phenotype file ${PHENO_SOURCE} is missing or empty")
if [[ ${#PHENO_PROBLEMS[@]} -gt 0 ]]; then
    if [[ $DRY_RUN == 1 ]]; then
        printf 'WARNING: %s (a launch would refuse)\n' "${PHENO_PROBLEMS[@]}"
    else
        printf 'REFUSED: %s\n' "${PHENO_PROBLEMS[@]}" >&2
        exit 2
    fi
fi
# Step 0 ingests the bundle's structure (CP-07) unless STRUCTURE_SOURCE=compute: its files must be there
if [[ "${STRUCTURE_SOURCE}" == bundle && "${SKIP_PCAIR}" != true ]]; then
    MISSING_STRUCT=()
    for f in genotypes/cohort.hardcall.maf01.fam ancestry/strata.tsv structure/cohort.pcair.pcs.tsv \
             structure/cohort.pcrelate.grm.rds structure/cohort.pcrelate.kinship.rds; do
        [[ -s "${HANDOFF_DIR}/${f}" ]] || MISSING_STRUCT+=("${HANDOFF_DIR}/${f}")
    done
    if [[ ${#MISSING_STRUCT[@]} -gt 0 ]]; then
        if [[ $DRY_RUN == 1 ]]; then printf 'WARNING: missing %s (a launch would refuse)\n' "${MISSING_STRUCT[@]}"
        else printf 'REFUSED: missing %s (STRUCTURE_SOURCE=bundle reads the shipped structure; it never recomputes it)\n' "${MISSING_STRUCT[@]}" >&2; exit 2; fi
    fi
fi
# Covariates (CP-08): unset means the study's list plus PC1..PC5, resolved once here so every job
# uses the same list and params.tsv records it (set_by bundle). Every name must be a column the
# models will see: a PC in the PC files step 0 lays out (or one of the N_PCS PCs it computes),
# anything else in the phenotype file. Whether the values are numeric is checked in each task.
if [[ "${RR_SET_BY[COVARIATES]}" == job ]]; then
    if COVARIATES="$(default_covariates)"; then RR_SET_BY[COVARIATES]=bundle; else COVARIATES=""; fi
    export COVARIATES
fi
[[ -z "${COVARIATES}" ]] || rr_check_domains COVARIATES || exit 2
COV_PROBLEMS=()
[[ -n "${COVARIATES}" ]] || COV_PROBLEMS+=("no covariates: COVARIATES is not set and ${HANDOFF_DIR}/phenotypes/covariates.txt is missing")
PC_FILES=()        # the files whose PC columns the tasks will use, where they exist already
PC_SETS=(cohort); [[ "${PC_SCALE}" == stratum ]] && PC_SETS+=("${STRATA[@]}")   # PC_SCALE=cohort: cohort PCs only
for S in "${PC_SETS[@]}"; do
    if [[ "${SKIP_PCAIR}" == true ]]; then f="${REPO}/results/pcair/${S}.phenotypes.with_pcs.tsv"
    elif [[ "${STRUCTURE_SOURCE}" == bundle && "$S" == cohort ]]; then f="${HANDOFF_DIR}/structure/cohort.pcair.pcs.tsv"
    elif [[ "${STRUCTURE_SOURCE}" == bundle ]]; then f="${HANDOFF_DIR}/structure/strata/${S}.pcair.pcs.tsv"
    else continue; fi
    [[ -s "$f" ]] && PC_FILES+=("$f")
done
IFS=',' read -ra COV_NAMES <<< "${COVARIATES}"
for c in "${COV_NAMES[@]}"; do
    if [[ "$c" =~ ^PC([0-9]+)$ ]]; then
        if [[ "${STRUCTURE_SOURCE}" == compute && "${SKIP_PCAIR}" != true ]]; then
            (( BASH_REMATCH[1] <= N_PCS )) || COV_PROBLEMS+=("covariate ${c}: STRUCTURE_SOURCE=compute writes N_PCS=${N_PCS} PCs")
        fi
        for f in ${PC_FILES[@]+"${PC_FILES[@]}"}; do
            grep -qxF -- "$c" <<< "$(head -1 "$f" | tr '\t' '\n')" || COV_PROBLEMS+=("covariate ${c} is not a column of ${f}")
        done
    elif [[ -s "${PHENO_SOURCE}" ]]; then
        grep -qxF -- "$c" <<< "${PHENO_HEADER}" || COV_PROBLEMS+=("covariate ${c} is not a column of ${PHENO_SOURCE}")
    fi
done
if [[ ${#COV_PROBLEMS[@]} -gt 0 ]]; then
    if [[ $DRY_RUN == 1 ]]; then printf 'WARNING: %s (a launch would refuse)\n' "${COV_PROBLEMS[@]}"
    else printf 'REFUSED: %s\n' "${COV_PROBLEMS[@]}" >&2; exit 2; fi
fi
RR_EXTRA_PARAMS="$(printf 'SURVIVAL_COLUMNS\t%s\tderived\nHANDOFF_DIR\t%s\tderived' "$(IFS=,; echo "${SURV_COLS[*]:-none}")" "${HANDOFF_DIR}")"

DEFERRED=()
for step in PRS COLOC VALIDATION GXG MRMEGA; do
    v="RUN_${step}"
    [[ "${!v}" == true ]] || DEFERRED+=("${step}")
done
RR_DEFERRED="${DEFERRED[*]:-none}"

# Run record: the run_id is computed once, here, and exported with the directory, so every job
# of this launch writes into the same results/run_record/<run_id>/.
RR_ROOT="${GWAS_RESULTS}/run_record"
if [[ $DRY_RUN == 1 ]]; then
    GWAS_RUN_ID="$(rr_preview_run_id "${RR_ROOT}")"
else
    GWAS_RUN_ID="$(rr_claim_run_dir "${RR_ROOT}")"
fi
GWAS_RUN_DIR="${RR_ROOT}/${GWAS_RUN_ID}"
export GWAS_RUN_ID GWAS_RUN_DIR
# Every job of this launch files its QC verdicts and sample counts here (CP-11; bin/qc_assert.R)
export GWAS_QC_DIR="${GWAS_RUN_DIR}/qc"
[[ $DRY_RUN == 1 ]] || rr_write_launch_record "${GWAS_RUN_DIR}" "$0" "$@"

trait_model() { if trait_is_survival "$1"; then echo survival; else echo binary; fi; }

# QTL types for colocalization
QTL_TYPES=("eqtl" "sqtl" "pqtl" "mqtl" "caqtl" "hqtl")

echo "=============================================="
echo "Multi-Ancestry GWAS Pipeline Submission"
echo "=============================================="
echo "Strata:    ${STRATA[*]}"
echo "Traits:    ${TRAITS[*]}"
echo "Partition: ${SLURM_PARTITION}"
echo "Account:   ${SLURM_ACCOUNT}"
echo "Logs:      ${LOG_DIR}"
echo "Results:   ${GWAS_RESULTS}${RUN_LABEL:+ (RUN_LABEL=${RUN_LABEL})}"
if [[ $DRY_RUN == 1 ]]; then
    echo "Record:    ${GWAS_RUN_DIR} (dry run: not created)"
    # what params.tsv would hold, so two dry runs can be compared row by row
    echo "params.tsv (dry run: not written):"
    rr_params_tsv | sed 's/^/  /'
    echo "end of params.tsv"
else
    echo "Record:    ${GWAS_RUN_DIR}"
fi
echo "Switches you set: $(awk -F'\t' 'NR > 1 && $3 == "environment" {printf "%s%s=%s", sep, $1, $2; sep = " "} END {if (!sep) printf "none"}' <(rr_params_tsv))"
echo "=============================================="
# Printed on every launch, so the log shows these steps were left out on purpose
if [[ ${#DEFERRED[@]} -gt 0 ]]; then
    echo "deferred: ${DEFERRED[*]} (RUN_*=false)"
else
    echo "deferred: none"
fi

# ---------------------------------------------------------------------------
# submit LABEL SBATCH_ARGS... : every sbatch goes through here, so the site arguments are
# added in one place. Sets JOB to the new job id (a placeholder "dryN" under --dry-run, so the
# printed commands still show which job waits for which).
#   --export=ALL               the job sees this launcher's environment (GWAS_REPO, LOG_DIR, switches)
#   --kill-on-invalid-dep=yes  a job whose afterok dependency failed is cancelled, not left pending forever
# A heredoc written after a submit call is the job script: sbatch reads its script from
# standard input when no script file is named.
# ---------------------------------------------------------------------------
SBATCH_SITE=(--parsable --export=ALL --partition="${SLURM_PARTITION}" --account="${SLURM_ACCOUNT}"
             --kill-on-invalid-dep=yes)
JOB=""
DRY_N=0
ALL_JOBS=()          # every job id, for the finalise job's afterany dependency
submit() {
    local label="$1"; shift
    if [[ $DRY_RUN == 1 ]]; then
        DRY_N=$((DRY_N + 1)); JOB="dry${DRY_N}"
        printf '  [DRY-RUN] %s: sbatch' "$label"; printf ' %q' "${SBATCH_SITE[@]}" "$@"; printf '\n'
    else
        JOB=$(sbatch "${SBATCH_SITE[@]}" "$@")
        JOB="${JOB%%;*}"                       # --parsable prints "jobid" or "jobid;cluster"
        echo "  ${label}: job ${JOB}"
        [[ "${NOTE_JOB:-1}" == 0 ]] || rr_note_job "${GWAS_RUN_DIR}" "${label}" "${JOB}"
    fi
    ALL_JOBS+=("${JOB}")
}

# ============================================================================
# STEP 0: GENESIS PC-AiR + PC-Relate (ancestry PCs + GRM for EVERY model)
# ============================================================================
echo ""
echo "STEP 0: Submitting structure (${STRUCTURE_SOURCE}: $([[ ${STRUCTURE_SOURCE} == bundle ]] && echo "ingest from ${HANDOFF_DIR}" || echo "PC-AiR / PC-Relate"))..."
DEPEND_PCAIR=""
if [[ "${SKIP_PCAIR}" != "true" ]]; then
    submit "structure (${STRUCTURE_SOURCE})" \
        --output="${LOG_DIR}/pcair_%j.out" --error="${LOG_DIR}/pcair_%j.err" \
        "${SCRIPT_DIR}/submit_pcair.sh"
    DEPEND_PCAIR="--dependency=afterok:${JOB}"
else
    echo "  Skipped (SKIP_PCAIR=true): reusing the structure already in results/pcair"
    # Carry that structure run's QC verdicts into this run, marked as carried, so the GWAS gate
    # sees them. Without them the gate stops on "missing assertion": nothing checked this structure.
    STRUCT_QC="${REPO}/results/pcair/qc_assertions.tsv"
    if [[ ! -f "${STRUCT_QC}" ]]; then
        echo "  WARNING: ${STRUCT_QC} not found; the GWAS jobs' QC gate will stop on a missing ${QC_REQUIRE_STRUCTURE}"
    elif [[ $DRY_RUN == 1 ]]; then
        echo "  [DRY-RUN] would carry $(( $(wc -l < "${STRUCT_QC}") - 1 )) QC verdict(s) from ${STRUCT_QC}"
    else
        mkdir -p "${GWAS_QC_DIR}/assertions.d"
        awk -F'\t' -v OFS='\t' -v src="${STRUCT_QC}" 'NR > 1 { $4 = $4 " [carried from " src ", SKIP_PCAIR=true]"; print }' \
            "${STRUCT_QC}" > "${GWAS_QC_DIR}/assertions.d/structure.carried.tsv"
        echo "  QC verdicts carried from ${STRUCT_QC}"
    fi
fi

# ============================================================================
# STEP 1: Tractor-GENESIS GWAS for every trait x stratum (array over chromosomes)
# ============================================================================
# ALL traits use the same adapter so they are comparable:
#   mrd -> binary (logistic null), efs/os -> survival (Cox null + PC-Relate GRM)
echo ""
echo "STEP 1: Submitting Tractor-GENESIS GWAS jobs..."
declare -A GWAS_JOBS

for trait in "${TRAITS[@]}"; do
    model=$(trait_model "$trait")
    for stratum in "${STRATA[@]}"; do
        submit "GWAS ${trait} (${model}) / ${stratum}, 22 chromosomes" ${DEPEND_PCAIR:+"$DEPEND_PCAIR"} \
            --output="${LOG_DIR}/gwas_${trait}_${stratum}_%A_%a.out" \
            --error="${LOG_DIR}/gwas_${trait}_${stratum}_%A_%a.err" \
            "${SCRIPT_DIR}/submit_gwas_array.sh" "${stratum}" "${trait}" "${model}"
        GWAS_JOBS["${trait}.${stratum}"]="$JOB"
    done
done

# ============================================================================
# STEP 2: per trait: strata disjoint and complete, chromosomes combined, strata pooled -> POOLED
# ============================================================================
echo ""
echo "STEP 2: Submitting combine + meta-analysis jobs (one per trait)..."

declare -A COMBINE_JOBS

for trait in "${TRAITS[@]}"; do
    dep=""
    for stratum in "${STRATA[@]}"; do
        j="${GWAS_JOBS["${trait}.${stratum}"]:-}"
        [[ -n "$j" ]] && dep="${dep}:${j}"
    done
    DEPEND_GWAS=""; [[ -n "$dep" ]] && DEPEND_GWAS="--dependency=afterok${dep}"

    # Quoted heredoc (<<'EOF'): the body below reaches sbatch exactly as written, nothing in it
    # is expanded here. Launch-time values reach the job through its environment instead
    # (sbatch --export=ALL): GWAS_REPO, STRATA_LIST, HANDOFF_DIR, and GWAS_TRAIT set on the submit line.
    GWAS_TRAIT="$trait" submit "combine + meta-analysis ${trait}" ${DEPEND_GWAS:+"$DEPEND_GWAS"} \
        --job-name="combine_${trait}" --time=4:00:00 --mem=32G --cpus-per-task=8 \
        --output="${LOG_DIR}/combine_${trait}_%j.out" <<'EOF'
#!/bin/bash
# Check the strata, combine chromosomes per stratum, pool the strata -> POOLED.ivw_descriptive
set -eu    # no pipefail: 'zcat | awk {exit}' below stops zcat early (SIGPIPE) by design
SCRIPT_DIR="${GWAS_REPO}"
source "${SCRIPT_DIR}/bin/slurm_utils.sh"
rr_begin combine "${GWAS_TRAIT} ${STRATA_LIST}"
# QC gate (CP-11): every GWAS task's join checks must be there and must not have failed
qc_gate_or_exit "${QC_REQUIRE_GWAS}" "combine ${GWAS_TRAIT}"
source "${SCRIPT_DIR}/bin/activate_env.sh" gwas-r
GWAS_DIR="${GWAS_RESULTS}/gwas/${GWAS_TRAIT}"
STRATA_TSV="${HANDOFF_DIR}/ancestry/strata.tsv"
IFS=',' read -ra STRATA <<< "${STRATA_LIST}"
join() { local IFS=,; echo "$*"; }          # join a b c -> a,b,c

# 1. The strata are disjoint (CP-22): no sample in two results/strata/<S>.samples.txt (step 0
#    writes each stratum's members from the bundle's strata.tsv). awk names the first shared sample.
SAMPLE_FILES=()
for stratum in "${STRATA[@]}"; do
    f="${SCRIPT_DIR}/results/strata/${stratum}.samples.txt"
    [[ -s "$f" ]] || { echo "ERROR: ${f} is missing (step 0 writes each stratum's members from ${STRATA_TSV})" >&2; exit 1; }
    SAMPLE_FILES+=("$f")
done
overlap="$(awk 'FNR == 1 {s = FILENAME; sub(/.*\//, "", s); sub(/\.samples\.txt$/, "", s)}
               ($1 in seen) && seen[$1] != s {n++; if (!ex) ex = seen[$1] " and " s " share " $1; next}
               {seen[$1] = s}
               END {if (n) print n " sample(s); e.g. " ex}' "${SAMPLE_FILES[@]}")"
[[ -z "$overlap" ]] || { echo "ERROR: strata not disjoint: ${overlap} (results/strata/<S>.samples.txt)" >&2; exit 1; }
echo "Strata disjoint: $(join "${STRATA[@]}") ($(cat "${SAMPLE_FILES[@]}" | wc -l) samples)"

# 2. Each stratum has all 22 chromosomes' results, or was skipped below MIN_STRATUM_N on every
#    chromosome (each task's chr<N>.run_counts.tsv says which). Anything else stops the job,
#    naming the stratum: a missing stratum is never pooled around in silence.
ANALYSED=(); SKIPPED=(); PROBLEMS=()
for stratum in "${STRATA[@]}"; do
    d="${GWAS_DIR}/${stratum}"; n_ok=0; n_skip=0; missing=()
    for c in $(seq 1 22); do
        st="$(awk -F'\t' '$1 == "status" {print $2}' "${d}/chr${c}.run_counts.tsv" 2>/dev/null || true)"
        case "$st" in
            skipped_below_min_stratum_n) n_skip=$((n_skip + 1)) ;;
            ok) if [[ -s "${d}/chr${c}.tractor_genesis.genomic_order.tsv.gz" ]]; then n_ok=$((n_ok + 1)); else missing+=("chr${c}"); fi ;;
            *)  missing+=("chr${c}") ;;
        esac
    done
    if [[ $n_skip -eq 22 ]]; then
        SKIPPED+=("$stratum"); echo "Stratum ${stratum}: skipped on every chromosome (below MIN_STRATUM_N); not pooled"
    elif [[ ${#missing[@]} -gt 0 ]]; then
        PROBLEMS+=("stratum ${stratum}: no results for ${#missing[@]} chromosome(s): ${missing[*]}")
    elif [[ $n_skip -gt 0 ]]; then
        PROBLEMS+=("stratum ${stratum}: skipped below MIN_STRATUM_N on ${n_skip} chromosome(s) but analysed on ${n_ok}")
    else
        ANALYSED+=("$stratum")
    fi
done
if [[ ${#PROBLEMS[@]} -gt 0 ]]; then printf 'ERROR: %s\n' "${PROBLEMS[@]}" >&2; exit 1; fi

# 3. Concatenate chromosomes per analysed stratum (genomic order files -> one sumstats per
#    stratum). A file may start with "## key=value" lines (e.g. "## covariates=..."): they are kept
#    once, and every chromosome must carry the same ones. awk prints the leading "##" lines (meta)
#    or, after them, the column header (hdr) or the data rows (rows).
meta() { zcat "$1" | awk '/^##/ {print; next} {exit}'; }
hdr()  { zcat "$1" | awk '!/^##/ {print; exit}'; }
rows() { zcat "$1" | awk '/^##/ && !h {next} h {print} {h = 1}'; }
for stratum in ${ANALYSED[@]+"${ANALYSED[@]}"}; do
    files=("${GWAS_DIR}/${stratum}"/chr*.tractor_genesis.genomic_order.tsv.gz)
    echo "Combining ${GWAS_TRAIT} / ${stratum}..."
    m0="$(meta "${files[0]}")"
    for f in "${files[@]}"; do
        [[ "$(meta "$f")" == "$m0" ]] || { echo "ERROR: $f has other '##' lines than ${files[0]}:" >&2; meta "$f" >&2; exit 1; }
    done
    { [[ -z "$m0" ]] || printf '%s\n' "$m0"; hdr "${files[0]}"; for f in "${files[@]}"; do rows "$f"; done; } |
        gzip > "${GWAS_DIR}/${stratum}.sumstats.gz"
done

# 4. Pool the analysed strata (the same ancestry background across strata). Descriptive: an IVW /
#    random-effects summary of per-stratum estimates, named so (POOLED.ivw_descriptive), with the
#    strata pooled, the strata skipped and the backgrounds in its "##" lines.
POOLED_OUT="${GWAS_DIR}/POOLED.ivw_descriptive.sumstats.gz"
if [[ ${#ANALYSED[@]} -lt 2 ]]; then
    echo "POOLED not written: ${#ANALYSED[@]} stratum(s) analysed ($(join ${ANALYSED[@]+"${ANALYSED[@]}"})); pooling needs two"
    rm -f "${POOLED_OUT}"
else
    echo "Pooling ${GWAS_TRAIT} across $(join "${ANALYSED[@]}")..."
    Rscript "${SCRIPT_DIR}/bin/meta_analyze_strata.R" \
        --input_dir "${GWAS_DIR}" \
        --strata "$(join "${ANALYSED[@]}")" \
        --strata_skipped "$(join ${SKIPPED[@]+"${SKIPPED[@]}"})" \
        --method random \
        --output "${POOLED_OUT}" \
        --threads 8
fi
EOF
    COMBINE_JOBS[$trait]="$JOB"
done

# ============================================================================
# STEP 2b: MR-MEGA across strata, per trait (after its combine job)    RUN_MRMEGA
# ============================================================================
# slurm/submit_mrmega.sh counts K, the strata with sumstats, and uses min(K-2, MRMEGA_PCS) axes; K=2
# is a fixed-effect meta-analysis (and says so), K<2 excludes the trait (CP-20).
echo ""
if [[ "${RUN_MRMEGA}" == true ]]; then
    echo "STEP 2b: Submitting MR-MEGA jobs (one per trait; K and pc_axes decided in the job)..."
    for trait in "${TRAITS[@]}"; do
        submit "MR-MEGA ${trait}" ${COMBINE_JOBS[$trait]:+--dependency=afterok:"${COMBINE_JOBS[$trait]}"} \
            --output="${LOG_DIR}/mrmega_${trait}_%j.out" --error="${LOG_DIR}/mrmega_${trait}_%j.err" \
            "${SCRIPT_DIR}/submit_mrmega.sh" "${trait}"
    done
else
    echo "STEP 2b: MR-MEGA deferred (RUN_MRMEGA=false)"
fi

# Single dependency handle for downstream steps (all traits combined). Steps 3, 4 and 6 all
# wait on the combine jobs, so this is set here, outside every RUN_* gate.
COMBINE_JOB=""
for trait in "${TRAITS[@]}"; do
    [[ -n "${COMBINE_JOBS[$trait]:-}" ]] && COMBINE_JOB="${COMBINE_JOB:+${COMBINE_JOB}:}${COMBINE_JOBS[$trait]}"
done
DEPEND_COMBINE=""
[[ -n "${COMBINE_JOB}" ]] && DEPEND_COMBINE="--dependency=afterok:${COMBINE_JOB}"

# ============================================================================
# STEP 3: PRS for each trait (depends on GWAS combine)                RUN_PRS
# ============================================================================
echo ""
declare -A PRS_JOBS
if [[ "${RUN_PRS}" == true ]]; then
echo "STEP 3: Submitting PRS jobs..."

for trait in "${TRAITS[@]}"; do
    submit "PRS ${trait} (6 methods)" ${DEPEND_COMBINE:+"$DEPEND_COMBINE"} \
        --output="${LOG_DIR}/prs_${trait}_%A_%a.out" --error="${LOG_DIR}/prs_${trait}_%A_%a.err" \
        "${SCRIPT_DIR}/submit_prs_array.sh" "${trait}"
    PRS_JOBS[$trait]="$JOB"
done
else
    echo "STEP 3: PRS deferred (RUN_PRS=false)"
fi

# ============================================================================
# STEP 4: Colocalization (parallel with PRS)                          RUN_COLOC
# ============================================================================
echo ""
if [[ "${RUN_COLOC}" == true ]]; then
echo "STEP 4: Submitting colocalization jobs (one 6-task array per trait)..."

# Per trait: each array colocalizes that trait's POOLED meta-analysis (written by its combine job,
# results/gwas/<trait>/POOLED.ivw_descriptive.sumstats.gz) and waits on that combine job only. Quoted heredoc, as
# in step 2; GWAS_TRAIT is set on the submit line.
for trait in "${TRAITS[@]}"; do
GWAS_TRAIT="$trait" submit "colocalization ${trait} (6 QTL types)" \
    ${COMBINE_JOBS[$trait]:+--dependency=afterok:"${COMBINE_JOBS[$trait]}"} \
    --job-name="coloc_${trait}" --array=1-6 --time=24:00:00 --mem=64G --cpus-per-task=16 \
    --output="${LOG_DIR}/coloc_${trait}_%A_%a.out" <<'EOF'
#!/bin/bash
set -euo pipefail
SCRIPT_DIR="${GWAS_REPO}"
source "${SCRIPT_DIR}/bin/slurm_utils.sh"
rr_begin coloc
source "${SCRIPT_DIR}/bin/activate_env.sh" gwas-r

# run_colocalization.R quietly skips a method whose R package is missing; refuse instead.
# coloc and susieR come with envs/gwas-r.yml; hyprcoloc is installed by hand (command there).
if ! Rscript "${SCRIPT_DIR}/bin/check_r_libs.R" coloc susieR hyprcoloc; then
    echo "REFUSED: colocalization needs the R packages coloc, susieR and hyprcoloc in the gwas-r env;" \
         "the MISSING line above names what is absent. coloc and susieR come with envs/gwas-r.yml;" \
         "hyprcoloc is installed by hand with the command in envs/gwas-r.yml." >&2
    exit 3
fi

# Map array ID to QTL type
QTL_TYPES=("eqtl" "sqtl" "pqtl" "mqtl" "caqtl" "hqtl")
QTL_TYPE="${QTL_TYPES[$((SLURM_ARRAY_TASK_ID - 1))]}"
GWAS_FILE="${GWAS_RESULTS}/gwas/${GWAS_TRAIT}/POOLED.ivw_descriptive.sumstats.gz"
OUT_DIR="${GWAS_RESULTS}/coloc/${GWAS_TRAIT}"
mkdir -p "${OUT_DIR}"
rr_args "${GWAS_TRAIT} ${QTL_TYPE}" "${GWAS_FILE}"

echo "Running colocalization for ${GWAS_TRAIT} / ${QTL_TYPE}..."

# Run all methods in parallel (coloc.susie, HyPrColoc, OPERA; OPERA falls back to an R
# multi-QTL SMR when the opera binary is not on PATH)
Rscript "${SCRIPT_DIR}/bin/run_colocalization.R" \
    --gwas "${GWAS_FILE}" \
    --qtl_dir "${SCRIPT_DIR}/data/qtl" \
    --qtl_types "${QTL_TYPE}" \
    --methods coloc_susie,hyprcoloc,opera \
    --run_opera \
    --output_prefix "${OUT_DIR}/${QTL_TYPE}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --verbose
EOF
done
else
    echo "STEP 4: colocalization deferred (RUN_COLOC=false)"
fi

# ============================================================================
# STEP 5: Validation and comparison (depends on PRS)                  RUN_VALIDATION
# ============================================================================
echo ""
if [[ "${RUN_VALIDATION}" == true ]]; then
echo "STEP 5: Submitting validation job..."

# Build PRS dependency
PRS_DEP=""
for trait in "${TRAITS[@]}"; do
    if [[ -n "${PRS_JOBS[$trait]:-}" ]]; then
        PRS_DEP="${PRS_DEP}:${PRS_JOBS[$trait]}"
    fi
done

if [[ -n "$PRS_DEP" ]]; then
    DEPEND_VAL="--dependency=afterok${PRS_DEP}"
else
    DEPEND_VAL=""
fi

# Quoted heredoc, as in step 2; the traits come from TRAITS_LIST in the job's environment.
submit "PRS validation" ${DEPEND_VAL:+"$DEPEND_VAL"} \
    --job-name=validation --time=4:00:00 --mem=32G --cpus-per-task=8 \
    --output="${LOG_DIR}/validation_%j.out" <<'EOF'
#!/bin/bash
set -eu
SCRIPT_DIR="${GWAS_REPO}"
source "${SCRIPT_DIR}/bin/slurm_utils.sh"
rr_begin validation "${TRAITS_LIST}"
source "${SCRIPT_DIR}/bin/activate_env.sh" gwas-r

echo "Comparing PRS methods across strata..."

# Combine all validation results
IFS=',' read -ra TRAITS <<< "${TRAITS_LIST}"
for trait in "${TRAITS[@]}"; do
    cat "${GWAS_RESULTS}/prs/${trait}"/*.validation.tsv |
        awk 'NR==1 || !/^method/' > \
        "${GWAS_RESULTS}/prs/${trait}_all_validation.tsv"
done

# Generate summary report
echo "Generating summary report..."
# ... report generation code
EOF
else
    echo "STEP 5: validation deferred (RUN_VALIDATION=false)"
fi

# ============================================================================
# STEP 6: GxG interaction tests on GWAS hits (END OF PIPELINE)
# ============================================================================
# Needs only the combined GWAS sumstats (not PRS/coloc), so it runs in
# parallel with steps 3-5. One array task per stratum, per trait.    RUN_GXG
echo ""
declare -A GXG_JOBS
if [[ "${RUN_GXG}" == true ]]; then
echo "STEP 6: Submitting GxG interaction jobs..."

for trait in "${TRAITS[@]}"; do
    model=$(trait_model "$trait")
    # one task per stratum plus POOLED (task 0): --array on the command line overrides the script's
    submit "GxG ${trait} (${model}; POOLED + ${#STRATA[@]} strata)" ${DEPEND_COMBINE:+"$DEPEND_COMBINE"} \
        --array="0-${#STRATA[@]}" --output="${LOG_DIR}/gxg_${trait}_%A_%a.out" --error="${LOG_DIR}/gxg_${trait}_%A_%a.err" \
        "${SCRIPT_DIR}/submit_gxg_array.sh" "${trait}" "${model}"
    GXG_JOBS[$trait]="$JOB"
done
else
    echo "STEP 6: GxG deferred (RUN_GXG=false)"
fi

# ============================================================================
# Run record finalise: once every job above has ended, however it ended (afterany), pull
# State, MaxRSS and Elapsed from sacct into trace.tsv and write summary.md
# ============================================================================
echo ""
echo "Submitting the run record finalise job..."
DEPEND_ALL="$(IFS=:; echo "${ALL_JOBS[*]}")"
# NOTE_JOB=0: the finalise job is not one of the jobs it reports on, so it stays out of jobs.tsv
NOTE_JOB=0 submit "run record finalise" ${DEPEND_ALL:+--dependency=afterany:"${DEPEND_ALL}"} \
    --output="${LOG_DIR}/run_record_finalize_%j.out" \
    "${REPO}/bin/run_record_finalize.sh" "${GWAS_RUN_DIR}"
[[ $DRY_RUN == 1 ]] || printf 'finalise_job\t%s\n' "${JOB}" >> "${GWAS_RUN_DIR}/run_info.tsv"

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "=============================================="
echo "Pipeline Submitted"
echo "=============================================="
echo "GWAS jobs: ${#TRAITS[@]} traits × ${#STRATA[@]} strata × 22 chromosomes"
[[ "${RUN_PRS}" == true ]]   && echo "PRS jobs: ${#TRAITS[@]} traits × 6 methods"
[[ "${RUN_COLOC}" == true ]] && echo "Colocalization: ${#TRAITS[@]} traits × ${#QTL_TYPES[@]} QTL types × 3 methods"
[[ "${RUN_GXG}" == true ]]   && echo "GxG interaction: ${#TRAITS[@]} traits × (POOLED + ${#STRATA[@]} strata) (+ ancestry heterogeneity)"
if [[ ${#DEFERRED[@]} -gt 0 ]]; then echo "deferred: ${DEFERRED[*]} (RUN_*=false)"; fi
echo "Logs: ${LOG_DIR}"
echo "Run record: ${GWAS_RUN_DIR} (trace.tsv and summary.md once the finalise job has run)"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Cancel all: scancel -u \$USER"
echo "=============================================="
