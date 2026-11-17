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
#   3. PRS: Run 6 methods × each trait                       only with RUN_PRS=true
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
# The settings file (config/gwas_settings.txt unless SETTINGS_FILE names another; bin/gwas_settings.sh):
# a relative SETTINGS_FILE is read from the directory you launch from, so it is made absolute here,
# before the cd below, and exported so every job reads the same file
if [[ -n "${SETTINGS_FILE:-}" && "${SETTINGS_FILE}" != none && "${SETTINGS_FILE}" != /* ]]; then SETTINGS_FILE="${PWD}/${SETTINGS_FILE}"; fi
[[ -z "${SETTINGS_FILE+x}" ]] || export SETTINGS_FILE
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
# Settings: config/gwas_settings.txt (or SETTINGS_FILE), one 'name = value' per line, this study's
# values; an environment variable of the same name beats it, and the built-in defaults of
# RR_SWITCHES (bin/run_record.sh) fill in what neither sets. rr_resolve_switches resolves and
# exports every switch, so every job sees the value that params.tsv records (with where it came
# from: default / settings / arm:<label> / environment). Nothing study-specific is built in.
# Customisable via the settings file or environment variables:
#   ARM=mrd_pos                             one arm.<label> of the settings file (e.g. the MRD-positive
#                                           EFS/OS GWAS); RUN_LABEL defaults to the label; ARM=all
#                                           runs the whole cohort, then every arm, one launch each
#   STRATA_LIST="EUR,AFR"                   strata analysed; default: every stratum in the bundle's
#                                           ancestry/strata.tsv (strata are assigned upstream)
#   TRAITS_LIST="efs,os"                    traits run (default: every trait.<name> of the settings file)
#   SURVIVAL_TRAITS / QUANTITATIVE_TRAITS   the models (default: the trait.<name> lines); a survival
#                                           trait reads <trait>_time and <trait>_event (TIME_COL_<trait> /
#                                           EVENT_COL_<trait> override one trait); the rest are binary
#   TIME_ORIGIN=eoi|diagnosis               survival clock: eoi = from the landmark in EOI_COL,
#                                           <trait>_time - EOI_COL, leaving out (and counting) a sample
#                                           whose event or censoring is at or before it
#   COX_STRATA=<column>                     survival only: each level gets its own baseline hazard
#   COX_EXTRA_COVARIATES=mrd_level          survival only: covariates added to COVARIATES
#   FACTOR_COVARIATES=cyto_risk:neutral_unknown   categorical covariates, column:reference level
#   SAMPLE_SUBSET=study_group=REDIAL        restrict the GWAS, PRS validation and GxG to the rows whose
#                                           column has one of the values; needs a RUN_LABEL
#   COVARIATES="sex,age_at_diagnosis,PC1"  phenotype columns, used by GWAS and GxG; default: the
#                                           bundle's phenotypes/covariates.txt, then PC1..PC<N_PCS_MODEL>
#   GXG_CUSTOM_VARIANTS=/path/list.txt      optional GRCh38 variant list for GxG
#   GXG_USE_DEFAULT_LOCI=true|false         include the default known-leukemia-loci table
#   LOG_DIR=/abs/path                       where every job log goes (default <repo>/logs)
#   RUN_PRS / RUN_COLOC / RUN_VALIDATION / RUN_GXG = true|false   steps 3-6 (default false)
#   MIN_STRATUM_N=50                        smallest stratum analysed (decided 10 Nov 2026)
#   MAC_MIN=10                              minor allele count per ancestry term; open question: the
#                                           handoff contract says 20, the code default is still 10
# ---------------------------------------------------------------------------
source "${REPO}/bin/slurm_utils.sh"    # helpers, and bin/run_record.sh (the switch table, rr_*)
gs_load || exit 2                      # the settings file: every line checked, refusals name the line
# ARM=all: the whole cohort, then the arms ARM_ALL names (default: every arm.<label> of the settings
# file, in file order), one launch (and one run record, under results/<label>/) each. Every one is
# dry-run first, so a refusal in the last stops the launch before the first is queued.
_arm_req="${GS_VALUE[ARM]:-}"; [[ -z "${ARM+x}" ]] || _arm_req="${ARM}"      # the environment beats the settings file
if [[ "${_arm_req}" == all ]]; then
    [[ -z "${RUN_LABEL:-}" ]] || { echo "REFUSED: ARM=all gives each arm its label as RUN_LABEL; unset RUN_LABEL" >&2; exit 2; }
    [[ ${#GS_ARMS[@]} -gt 0 ]] || { echo "REFUSED: ARM=all, but ${GS_FILE:-no settings file is read, so there} is no arm.<label> line" >&2; exit 2; }
    _all="${GS_VALUE[ARM_ALL]:-}"; [[ -z "${ARM_ALL+x}" ]] || _all="${ARM_ALL}"
    [[ -n "${_all}" ]] || _all="$(IFS=,; echo "${GS_ARMS[*]}")"
    IFS=',' read -ra _arms <<< "${_all}"
    for a in "${_arms[@]}"; do
        [[ -n "${GS_ARM[$a]+x}" ]] || { echo "REFUSED: ARM_ALL names ${a}, which no arm.${a} line of ${GS_FILE} defines (arms: ${GS_ARMS[*]})" >&2; exit 2; }
    done
    if [[ $DRY_RUN == 0 ]]; then
        for a in "" "${_arms[@]}"; do
            ARM="$a" "${SCRIPT_DIR}/submit_full_pipeline.sh" --dry-run > /dev/null 2> "${LOG_DIR}/.arm_check.$$" \
                || { echo "REFUSED: ARM=all: the ${a:-whole-cohort} launch would be refused; nothing was queued:" >&2; cat "${LOG_DIR}/.arm_check.$$" >&2; rm -f "${LOG_DIR}/.arm_check.$$"; exit 2; }
        done
        rm -f "${LOG_DIR}/.arm_check.$$"
    fi
    for a in "" "${_arms[@]}"; do
        echo "=== ARM=${a:-(none): the whole cohort} ==="
        ARM="$a" "${SCRIPT_DIR}/submit_full_pipeline.sh" "$@" || exit $?
    done
    exit 0
fi
rr_refuse_unknown_switches || exit 2
rr_resolve_switches || exit 2
# Every switch inside its domain (B-05; RR_SWITCHES in bin/run_record.sh), before anything else
rr_check_domains || exit 2
# An arm's results go under results/<label>/ unless RUN_LABEL says otherwise
if [[ -n "${ARM}" && -z "${RUN_LABEL+x}" ]]; then RUN_LABEL="${ARM}"; fi
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
# POOLED: the whole cohort in one fit (cohort PCs, data/tractor/cohort/ Tractor files); alone, since a
# stratum inside it would be counted twice when the strata are pooled
if [[ ",${STRATA_LIST}," == *,POOLED,* && "${STRATA_LIST}" != POOLED ]]; then
    echo "REFUSED: STRATA_LIST=${STRATA_LIST}: POOLED is the whole cohort; list it alone (STRATA_LIST=POOLED)" >&2; exit 2
fi
for stratum in "${STRATA[@]}"; do
    [[ ",${BUNDLE_STRATA}," == *",${stratum},"* || "${stratum}" == POOLED ]] \
        || { echo "REFUSED: stratum ${stratum} (STRATA_LIST) is not in ${STRATA_TSV} (it has: ${BUNDLE_STRATA}; or POOLED, the whole cohort)" >&2; exit 2; }
    [[ "${MIX_RETAIN}" == false && "${stratum}" == MIX ]] \
        && { echo "REFUSED: STRATA_LIST names MIX, but MIX_RETAIN=false" >&2; exit 2; }
done
IFS=',' read -ra TRAITS <<< "${TRAITS_LIST}"
[[ ${#TRAITS[@]} -gt 0 ]] \
    || { echo "REFUSED: no traits: declare them with trait.<name> lines in the settings file (${GS_FILE:-none is read; SETTINGS_FILE}), or set TRAITS_LIST" >&2; exit 2; }
# A trait is one model: survival, quantitative, or (neither list) binary; with trait.<name> lines in
# the settings file, every trait run must be declared there
for trait in "${TRAITS[@]}"; do
    if [[ ${#GS_TRAITS[@]} -gt 0 && " ${GS_TRAITS[*]} " != *" ${trait} "* ]]; then
        echo "REFUSED: trait ${trait} (TRAITS_LIST) is not declared in ${GS_FILE} (add trait.${trait} = binary | quantitative | survival)" >&2; exit 2
    fi
    if trait_is_survival "$trait" && trait_is_quantitative "$trait"; then
        echo "REFUSED: trait ${trait} is in both SURVIVAL_TRAITS and QUANTITATIVE_TRAITS" >&2; exit 2
    fi
done
# SAMPLE_SUBSET (decided 17 Nov 2026): a sub-cohort analysis (e.g. study_group=REDIAL, or an MRD arm)
# must not overwrite the full-cohort results, so it needs its own RUN_LABEL (an ARM gives its label).
# The GWAS, PRS validation and GxG all apply it. A column that says what the subset selects on (the
# subset column, SUBSET_LINKED_COLUMNS) cannot be a trait inside it: checked here; nor a covariate:
# checked once COVARIATES is resolved, below.
SUBSET_LINKED=""
if [[ -n "${SAMPLE_SUBSET}" ]]; then
    [[ "${SAMPLE_SUBSET}" =~ ^[A-Za-z][A-Za-z0-9_.]*=[^=]+$ ]] \
        || { echo "REFUSED: SAMPLE_SUBSET=${SAMPLE_SUBSET}; expected column=value[,value...] (e.g. study_group=REDIAL)" >&2; exit 2; }
    [[ -n "${RUN_LABEL:-}" ]] \
        || { echo "REFUSED: SAMPLE_SUBSET=${SAMPLE_SUBSET} needs a RUN_LABEL (e.g. RUN_LABEL=redial_only), so the full-cohort results are not overwritten" >&2; exit 2; }
    # the subset column, and with it every SUBSET_LINKED_COLUMNS column when it is one of them
    SUBSET_LINKED=",${SAMPLE_SUBSET%%=*},"; _linked_desc="${SUBSET_LINKED_COLUMNS}"; [[ -n "${_linked_desc}" ]] || _linked_desc=none
    [[ ",${SUBSET_LINKED_COLUMNS}," == *"${SUBSET_LINKED}"* ]] && SUBSET_LINKED=",${SUBSET_LINKED_COLUMNS},"
    _bad=(); _keep=()
    for t in "${TRAITS[@]}"; do if [[ "${SUBSET_LINKED}" == *",${t},"* ]]; then _bad+=("$t"); else _keep+=("$t"); fi; done
    [[ ${#_bad[@]} -eq 0 ]] || { echo "REFUSED: SAMPLE_SUBSET=${SAMPLE_SUBSET} selects on ${SAMPLE_SUBSET%%=*} (linked columns: ${_linked_desc}), so ${_bad[*]} cannot be a trait inside it (set TRAITS_LIST=$(IFS=,; echo "${_keep[*]:-}"))" >&2; exit 2; }
fi
# COX_EXTRA_COVARIATES must not repeat the stratifying column (or, below, a covariate in COVARIATES)
if [[ -n "${COX_EXTRA_COVARIATES}" ]]; then
    IFS=',' read -ra _cx <<< "${COX_EXTRA_COVARIATES}"
    for c in "${_cx[@]}"; do
        [[ "$c" == "${COX_STRATA}" ]] && { echo "REFUSED: ${c} is both COX_STRATA and in COX_EXTRA_COVARIATES (stratify on it or adjust for it, not both)" >&2; exit 2; }
    done
fi

# Steps 2b and 3-6 run only when asked for (their RUN_* are true / false: rr_check_domains above)
if [[ "${RUN_VALIDATION}" == true && "${RUN_PRS}" != true ]]; then
    echo "REFUSED: RUN_VALIDATION=true needs RUN_PRS=true (validation summarises the PRS jobs of this launch)" >&2
    exit 2
fi
# Each stratum's Tractor ancestry list (CP-28) must exist before its 22-task arrays are queued:
# submit_gwas_array.sh reads it in every task and refuses without it. A dry run lists what is
# missing; a launch refuses.
MISSING_ANC=()
if [[ "${RUN_TRACTOR_PREP}" == true ]]; then
    # steps 0b-0d write each stratum's Tractor files and (when absent) its ancestries.txt; they
    # need the bundle's FLARE output for every chromosome and its index map
    for n in $(seq 1 22); do
        [[ -s "${HANDOFF_DIR}/local_ancestry/flare/chr${n}.anc.vcf.gz" ]] || MISSING_ANC+=("${HANDOFF_DIR}/local_ancestry/flare/chr${n}.anc.vcf.gz")
    done
    [[ -s "${HANDOFF_DIR}/local_ancestry/flare_index_map.tsv" ]] || MISSING_ANC+=("${HANDOFF_DIR}/local_ancestry/flare_index_map.tsv")
else
    for stratum in "${STRATA[@]}"; do
        anc_f="${REPO}/data/tractor/${stratum}/ancestries.txt"
        [[ "${stratum}" == POOLED ]] && anc_f="${REPO}/data/tractor/ancestries.txt"   # the cohort's list
        [[ -s "${anc_f}" ]] || MISSING_ANC+=("${anc_f}")
    done
fi
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

# Every trait's columns must be in the phenotype file: a binary or quantitative trait its own
# column, a survival trait <trait>_time and <trait>_event (or its TIME_COL_<trait> /
# EVENT_COL_<trait>), plus eoi_day under TIME_ORIGIN=eoi, the COX_STRATA column and any
# COX_EXTRA_COVARIATES; SAMPLE_SUBSET's column, with at least one row matching. Checked here,
# before anything is queued. A dry run warns; a launch refuses.
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
        [[ "${TIME_ORIGIN}" == eoi ]] && need+=("${EOI_COL}")
        [[ -n "${COX_STRATA}" ]] && need+=("${COX_STRATA}")
        [[ -n "${COX_EXTRA_COVARIATES}" ]] && { IFS=',' read -ra _cx <<< "${COX_EXTRA_COVARIATES}"; need+=("${_cx[@]}"); }
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
# VARIANT_INCLUDE (relative: from the repo root) must list something before 22 x strata tasks read it
if [[ -n "${VARIANT_INCLUDE}" ]]; then
    VI="${VARIANT_INCLUDE}"; [[ "$VI" == /* ]] || VI="${REPO}/${VI}"
    [[ -n "$(grep -v '^#' "$VI" 2>/dev/null | tr -d '[:space:]')" ]] \
        || PHENO_PROBLEMS+=("VARIANT_INCLUDE=${VARIANT_INCLUDE}: ${VI} is missing or lists no variant (bin/mrmega_lead_snps.sh writes it)")
fi
if [[ -n "${SAMPLE_SUBSET}" && -s "${PHENO_SOURCE}" ]]; then
    SUBSET_COL="${SAMPLE_SUBSET%%=*}"
    if grep -qxF -- "${SUBSET_COL}" <<< "${PHENO_HEADER}"; then
        SUBSET_N="$(awk -F'\t' -v col="${SUBSET_COL}" -v vals="${SAMPLE_SUBSET#*=}" '
            NR == 1 { for (i = 1; i <= NF; i++) if ($i == col) c = i; n = split(vals, v, ","); for (j = 1; j <= n; j++) { gsub(/^ +| +$/, "", v[j]); want[v[j]] = 1 }; next }
            ($c in want) { k++ } END { print k + 0 }' "${PHENO_SOURCE}")"
        [[ "${SUBSET_N}" -gt 0 ]] || PHENO_PROBLEMS+=("SAMPLE_SUBSET=${SAMPLE_SUBSET}: no row of ${PHENO_SOURCE} has ${SUBSET_COL} = ${SAMPLE_SUBSET#*=}")
    else
        PHENO_PROBLEMS+=("SAMPLE_SUBSET=${SAMPLE_SUBSET}: column ${SUBSET_COL} is not in ${PHENO_SOURCE}")
    fi
fi
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
# once COVARIATES is resolved (the bundle's list when unset): COX_EXTRA_COVARIATES must add, not repeat
if [[ -n "${COX_EXTRA_COVARIATES}" ]]; then
    IFS=',' read -ra _cx <<< "${COX_EXTRA_COVARIATES}"
    for c in "${_cx[@]}"; do
        [[ ",${COVARIATES}," == *",${c},"* ]] && { echo "REFUSED: ${c} is in both COVARIATES and COX_EXTRA_COVARIATES" >&2; exit 2; }
    done
fi
# FACTOR_COVARIATES' columns are in the model too (the adapter adds the ones COVARIATES does not name)
FACTOR_COLS=""
if [[ -n "${FACTOR_COVARIATES}" ]]; then FACTOR_COLS="$(tr ',' '\n' <<< "${FACTOR_COVARIATES}" | cut -d: -f1 | paste -sd, -)"; fi
# Inside a subset no covariate may be a column the subset selects on (e.g. no MRD covariate inside an
# MRD arm); the refusal says what to set instead
if [[ -n "${SUBSET_LINKED}" ]]; then
    for sw in COVARIATES COX_EXTRA_COVARIATES COX_STRATA FACTOR_COVARIATES; do
        IFS=',' read -ra _v <<< "${!sw}"; _bad=(); _keep=()
        for c in ${_v[@]+"${_v[@]}"}; do if [[ "${SUBSET_LINKED}" == *",${c%%:*},"* ]]; then _bad+=("${c%%:*}"); else _keep+=("$c"); fi; done
        [[ ${#_bad[@]} -eq 0 ]] || { echo "REFUSED: SAMPLE_SUBSET=${SAMPLE_SUBSET} selects on ${SAMPLE_SUBSET%%=*} (linked columns: ${_linked_desc}), so ${_bad[*]} cannot be a covariate inside it (set ${sw}=$(IFS=,; echo "${_keep[*]:-}"))" >&2; exit 2; }
    done
fi
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
IFS=',' read -ra COV_NAMES <<< "${COVARIATES}${FACTOR_COLS:+,${FACTOR_COLS}}"
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
RR_EXTRA_PARAMS="$(printf 'SURVIVAL_COLUMNS\t%s\tderived\nHANDOFF_DIR\t%s\tderived\nSETTINGS_FILE\t%s\tderived' \
    "$(IFS=,; echo "${SURV_COLS[*]:-none}")" "${HANDOFF_DIR}" "${GS_FILE:-none}")"
# The landmark's column (TIME_ORIGIN=eoi): the GWAS, PRS validation and GxG all count from it
if [[ ${#SURV_COLS[@]} -gt 0 && "${TIME_ORIGIN}" == eoi ]]; then
    RR_EXTRA_PARAMS+="$(printf '\nEOI_COLUMN\t%s\tderived' "${EOI_COL}")"
fi

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

trait_model() {
    if trait_is_survival "$1"; then echo survival
    elif trait_is_quantitative "$1"; then echo quantitative
    else echo binary; fi
}

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
echo "Settings:  ${GS_FILE:-none (built-in defaults; SETTINGS_FILE)}${ARM:+ (ARM=${ARM})}"
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
# STEPS 0b-0d: the Tractor files (CP-16)                                 RUN_TRACTOR_PREP
# ============================================================================
# 0b: FLARE -> the cohort's Tractor files, 22 tasks in the lai env; 0c: each stratum's ancestries
# (after 0b and step 0, which writes results/strata); 0d: each stratum x chromosome's files. The
# GWAS waits on 0d as well as on step 0. Off: data/tractor/<S>/ is prepared beforehand.
DEPEND_GWAS_IN="${DEPEND_PCAIR}"
if [[ "${RUN_TRACTOR_PREP}" == true ]]; then
    echo ""
    echo "STEPS 0b-0d: Submitting the Tractor extraction and the per-stratum layout..."
    submit "Tractor extraction (FLARE -> Tractor), 22 chromosomes" \
        --output="${LOG_DIR}/la_extract_%A_%a.out" --error="${LOG_DIR}/la_extract_%A_%a.err" \
        "${SCRIPT_DIR}/submit_la_extract.sh"
    J_EXTRACT="$JOB"
    J_STRUCT="${DEPEND_PCAIR#--dependency=afterok:}"
    submit "Tractor layout: each stratum's ancestries" --dependency=afterok:"${J_EXTRACT}${J_STRUCT:+:${J_STRUCT}}" \
        --output="${LOG_DIR}/la_ancestries_%j.out" --error="${LOG_DIR}/la_ancestries_%j.err" \
        "${SCRIPT_DIR}/submit_tractor_layout.sh" ancestries
    submit "Tractor layout: ${#STRATA[@]} strata x 22 chromosomes" --dependency=afterok:"${JOB}" --array="0-$(( 22 * ${#STRATA[@]} - 1 ))" \
        --output="${LOG_DIR}/la_layout_%A_%a.out" --error="${LOG_DIR}/la_layout_%A_%a.err" \
        "${SCRIPT_DIR}/submit_tractor_layout.sh" split
    DEPEND_GWAS_IN="--dependency=afterok:${JOB}${J_STRUCT:+:${J_STRUCT}}"
fi

# ============================================================================
# STEP 1: Tractor-GENESIS GWAS for every trait x stratum (array over chromosomes)
# ============================================================================
# ALL traits use the same adapter so they are comparable:
#   mrd -> binary (logistic null), mrd_level -> quantitative (linear null, the four MRD levels as
#   a trend), efs/os -> survival (Cox null + PC-Relate GRM; from the end of induction under
#   TIME_ORIGIN=eoi, baseline stratified on COX_STRATA)
echo ""
echo "STEP 1: Submitting Tractor-GENESIS GWAS jobs..."
declare -A GWAS_JOBS

for trait in "${TRAITS[@]}"; do
    model=$(trait_model "$trait")
    for stratum in "${STRATA[@]}"; do
        submit "GWAS ${trait} (${model}) / ${stratum}, 22 chromosomes" ${DEPEND_GWAS_IN:+"$DEPEND_GWAS_IN"} \
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
    [[ "${stratum}" == POOLED ]] && continue      # STRATA_LIST=POOLED: the whole cohort, nothing to be disjoint from
    f="${SCRIPT_DIR}/results/strata/${stratum}.samples.txt"
    [[ -s "$f" ]] || { echo "ERROR: ${f} is missing (step 0 writes each stratum's members from ${STRATA_TSV})" >&2; exit 1; }
    SAMPLE_FILES+=("$f")
done
[[ ${#SAMPLE_FILES[@]} -gt 0 ]] && overlap="$(awk 'FNR == 1 {s = FILENAME; sub(/.*\//, "", s); sub(/\.samples\.txt$/, "", s)}
               ($1 in seen) && seen[$1] != s {n++; if (!ex) ex = seen[$1] " and " s " share " $1; next}
               {seen[$1] = s}
               END {if (n) print n " sample(s); e.g. " ex}' "${SAMPLE_FILES[@]}")"
[[ -z "${overlap:-}" ]] || { echo "ERROR: strata not disjoint: ${overlap} (results/strata/<S>.samples.txt)" >&2; exit 1; }
if [[ ${#SAMPLE_FILES[@]} -gt 0 ]]; then echo "Strata disjoint: $(join "${STRATA[@]}") ($(cat "${SAMPLE_FILES[@]}" | wc -l) samples)"
else echo "Stratum POOLED: the whole cohort"; fi

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
    # chromosomes in genomic order, 1..22 (a chr* glob sorts as text: chr1, chr10, chr11, ..., chr2)
    files=(); for c in $(seq 1 22); do files+=("${GWAS_DIR}/${stratum}/chr${c}.tractor_genesis.genomic_order.tsv.gz"); done
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
