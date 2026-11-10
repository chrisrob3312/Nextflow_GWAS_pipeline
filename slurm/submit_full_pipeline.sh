#!/bin/bash
# ============================================================================
# Master SLURM Submission Script
# ============================================================================
# Submits complete multi-ancestry GWAS pipeline with proper dependencies
#
# Workflow:
#   0. PC-AiR / PC-Relate (skip with SKIP_PCAIR=true)
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
#   STRATA_LIST="EUR,AAC,LAT1,LAT2,OTHER"   ancestry strata (LAT1/LAT2 separate; OTHER = EAS+SAS+small groups)
#   TRAITS_LIST="OS,DFS,MRD"                any phenotype columns
#   SURVIVAL_TRAITS="OS,DFS"                traits run as time-to-event (<trait>_time, <trait>_status);
#                                           everything else is binary (0/1) - DFS (relapse OR death)
#                                           avoids the need for a competing-risk model
#   COVARIATES="age,sex,PC1,...,PC5"        any phenotype columns; used by GWAS and GxG
#   GXG_CUSTOM_VARIANTS=/path/list.txt      optional GRCh38 variant list for GxG
#   GXG_USE_DEFAULT_LOCI=true|false         include the default known-leukemia-loci table
#   LOG_DIR=/abs/path                       where every job log goes (default <repo>/logs)
#   RUN_PRS / RUN_COLOC / RUN_VALIDATION / RUN_GXG = true|false   steps 3-6 (default false)
#   MIN_STRATUM_N=50                        smallest stratum analysed (decided 10 Nov 2026)
#   MAC_MIN=10                              minor allele count per ancestry term; open question: the
#                                           handoff contract says 20, the code default is still 10
# ---------------------------------------------------------------------------
source "${REPO}/bin/run_record.sh"
rr_refuse_unknown_switches || exit 2
rr_resolve_switches
IFS=',' read -ra STRATA <<< "${STRATA_LIST}"
IFS=',' read -ra TRAITS <<< "${TRAITS_LIST}"

# Steps 3-6 run only when asked for. Anything but true/false is refused rather than guessed.
for step in PRS COLOC VALIDATION GXG; do
    v="RUN_${step}"
    case "${!v}" in
        true|false) ;;
        *) echo "REFUSED: ${v}=${!v}; use true or false" >&2; exit 2 ;;
    esac
done
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

DEFERRED=()
for step in PRS COLOC VALIDATION GXG; do
    v="RUN_${step}"
    [[ "${!v}" == true ]] || DEFERRED+=("${step}")
done
RR_DEFERRED="${DEFERRED[*]:-none}"

# Run record: the run_id is computed once, here, and exported with the directory, so every job
# of this launch writes into the same results/run_record/<run_id>/.
RR_ROOT="${REPO}/results/run_record"
if [[ $DRY_RUN == 1 ]]; then
    GWAS_RUN_ID="$(rr_preview_run_id "${RR_ROOT}")"
else
    GWAS_RUN_ID="$(rr_claim_run_dir "${RR_ROOT}")"
fi
GWAS_RUN_DIR="${RR_ROOT}/${GWAS_RUN_ID}"
export GWAS_RUN_ID GWAS_RUN_DIR
[[ $DRY_RUN == 1 ]] || rr_write_launch_record "${GWAS_RUN_DIR}" "$0" "$@"

trait_model() { [[ ",${SURVIVAL_TRAITS}," == *",$1,"* ]] && echo survival || echo binary; }

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
if [[ $DRY_RUN == 1 ]]; then
    echo "Record:    ${GWAS_RUN_DIR} (dry run: not created)"
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
echo "STEP 0: Submitting PC-AiR / PC-Relate..."
DEPEND_PCAIR=""
if [[ "${SKIP_PCAIR:-false}" != "true" ]]; then
    submit "PC-AiR / PC-Relate" \
        --output="${LOG_DIR}/pcair_%j.out" --error="${LOG_DIR}/pcair_%j.err" \
        "${SCRIPT_DIR}/submit_pcair.sh"
    DEPEND_PCAIR="--dependency=afterok:${JOB}"
else
    echo "  Skipped (SKIP_PCAIR=true): using data/kinship + PCs already in the phenotype file"
fi

# ============================================================================
# STEP 1: Tractor-GENESIS GWAS for every trait x stratum (array over chromosomes)
# ============================================================================
# ALL traits use the same adapter so they are comparable:
#   MRD -> binary (logistic null), DFS/OS -> survival (Cox null + PC-Relate GRM)
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
# STEP 2: Combine chromosomes per trait x stratum, meta-analyse strata -> POOLED
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
    # (sbatch --export=ALL): GWAS_REPO, STRATA_LIST, and GWAS_TRAIT set on the submit line.
    GWAS_TRAIT="$trait" submit "combine + meta-analysis ${trait}" ${DEPEND_GWAS:+"$DEPEND_GWAS"} \
        --job-name="combine_${trait}" --time=4:00:00 --mem=32G --cpus-per-task=8 \
        --output="${LOG_DIR}/combine_${trait}_%j.out" <<'EOF'
#!/bin/bash
# Combine chromosomes per stratum, then meta-analyse the strata -> POOLED
set -eu    # no pipefail: 'zcat | head -1' below stops zcat early (SIGPIPE) by design
SCRIPT_DIR="${GWAS_REPO}"
source "${SCRIPT_DIR}/bin/slurm_utils.sh"
rr_begin combine "${GWAS_TRAIT} ${STRATA_LIST}"
source "${SCRIPT_DIR}/bin/activate_env.sh" gwas-r
GWAS_DIR="${SCRIPT_DIR}/results/gwas/${GWAS_TRAIT}"
IFS=',' read -ra STRATA <<< "${STRATA_LIST}"

# Concatenate chromosomes per stratum (genomic order files -> one sumstats per stratum)
for stratum in "${STRATA[@]}"; do
    files=("${GWAS_DIR}/${stratum}"/chr*.tractor_genesis.genomic_order.tsv.gz)
    [[ -e "${files[0]}" ]] || { echo "No results for ${stratum} (skipped: N below MIN_STRATUM_N?)"; continue; }
    echo "Combining ${GWAS_TRAIT} / ${stratum}..."
    { zcat "${files[0]}" | head -1; for f in "${files[@]}"; do zcat "$f" | tail -n +2; done; } |
        gzip > "${GWAS_DIR}/${stratum}.sumstats.gz"
done

# Meta-analyse across strata (same ancestry background pooled across strata)
echo "Meta-analysing ${GWAS_TRAIT} across strata..."
Rscript "${SCRIPT_DIR}/bin/meta_analyze_strata.R" \
    --input_dir "${GWAS_DIR}" \
    --strata "${STRATA_LIST}" \
    --output "${GWAS_DIR}/POOLED.sumstats.gz" \
    --threads 8
EOF
    COMBINE_JOBS[$trait]="$JOB"
done

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
# results/gwas/<trait>/POOLED.sumstats.gz) and waits on that combine job only. Quoted heredoc, as
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
GWAS_FILE="${SCRIPT_DIR}/results/gwas/${GWAS_TRAIT}/POOLED.sumstats.gz"
OUT_DIR="${SCRIPT_DIR}/results/coloc/${GWAS_TRAIT}"
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
    cat "${SCRIPT_DIR}/results/prs/${trait}"/*.validation.tsv |
        awk 'NR==1 || !/^method/' > \
        "${SCRIPT_DIR}/results/prs/${trait}_all_validation.tsv"
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
    submit "GxG ${trait} (${model}; POOLED + 7 strata)" ${DEPEND_COMBINE:+"$DEPEND_COMBINE"} \
        --output="${LOG_DIR}/gxg_${trait}_%A_%a.out" --error="${LOG_DIR}/gxg_${trait}_%A_%a.err" \
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
[[ "${RUN_GXG}" == true ]]   && echo "GxG interaction: ${#TRAITS[@]} traits × 8 strata (+ ancestry heterogeneity)"
if [[ ${#DEFERRED[@]} -gt 0 ]]; then echo "deferred: ${DEFERRED[*]} (RUN_*=false)"; fi
echo "Logs: ${LOG_DIR}"
echo "Run record: ${GWAS_RUN_DIR} (trace.tsv and summary.md once the finalise job has run)"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Cancel all: scancel -u \$USER"
echo "=============================================="
