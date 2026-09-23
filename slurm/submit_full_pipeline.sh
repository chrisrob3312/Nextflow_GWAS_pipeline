#!/bin/bash
# ============================================================================
# Master SLURM Submission Script
# ============================================================================
# Submits complete multi-ancestry GWAS pipeline with proper dependencies
#
# Workflow:
#   1. GWAS: Run Tractor-GENESIS for each ancestry stratum (22 chr × 5 strata)
#   2. Combine: Merge chromosomes, run meta-analysis for OTHER
#   3. PRS: Run 6 methods × 3 traits
#   4. Colocalization: Run parallel methods on all QTL types
#   5. Validation: Compare methods overall and by ancestry
#
# Usage:
#   ./submit_full_pipeline.sh [--dry-run]
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DRY_RUN="${1:-}"

# Create log directory
mkdir -p "${SCRIPT_DIR}/../logs"

# ---------------------------------------------------------------------------
# Customisable via environment variables (defaults for the leukemia cohort):
#   STRATA_LIST="EUR,AAC,LAT1,LAT2,OTHER"   ancestry strata (LAT1/LAT2 separate; OTHER = EAS+SAS+small groups)
#   TRAITS_LIST="OS,DFS,MRD"                any phenotype columns
#   SURVIVAL_TRAITS="OS,DFS"                traits run as time-to-event (<trait>_time, <trait>_status);
#                                           everything else is binary (0/1) - DFS (relapse OR death)
#                                           avoids the need for a competing-risk model
#   COVARIATES="age,sex,PC1,...,PC5"        any phenotype columns; used by GWAS and GxG
#   GXG_CUSTOM_VARIANTS=/path/list.txt      optional GRCh38 variant list for GxG
#   GXG_USE_DEFAULT_LOCI=true|false         include the default known-leukemia-loci table
# ---------------------------------------------------------------------------
IFS=',' read -ra STRATA <<< "${STRATA_LIST:-EUR,AAC,LAT1,LAT2,OTHER}"
IFS=',' read -ra TRAITS <<< "${TRAITS_LIST:-OS,DFS,MRD}"
SURVIVAL_TRAITS="${SURVIVAL_TRAITS:-OS,DFS}"
export COVARIATES="${COVARIATES:-age,sex,PC1,PC2,PC3,PC4,PC5}"
export GXG_CUSTOM_VARIANTS="${GXG_CUSTOM_VARIANTS:-}"
export GXG_USE_DEFAULT_LOCI="${GXG_USE_DEFAULT_LOCI:-true}"

trait_model() { [[ ",${SURVIVAL_TRAITS}," == *",$1,"* ]] && echo survival || echo binary; }

# QTL types for colocalization
QTL_TYPES=("eqtl" "sqtl" "pqtl" "mqtl" "caqtl" "hqtl")

echo "=============================================="
echo "Multi-Ancestry GWAS Pipeline Submission"
echo "=============================================="
echo "Strata: ${STRATA[*]}"
echo "Traits: ${TRAITS[*]}"
echo "=============================================="

# Function to submit or print command
submit_job() {
    local cmd="$1"
    if [[ "$DRY_RUN" == "--dry-run" ]]; then
        echo "[DRY-RUN] $cmd"
    else
        eval "$cmd"
    fi
}

# ============================================================================
# STEP 0: GENESIS PC-AiR + PC-Relate (ancestry PCs + GRM for EVERY model)
# ============================================================================
echo ""
echo "STEP 0: Submitting PC-AiR / PC-Relate..."
DEPEND_PCAIR=""
if [[ "${SKIP_PCAIR:-false}" != "true" ]]; then
    cmd="sbatch --parsable ${SCRIPT_DIR}/submit_pcair.sh"
    if [[ "$DRY_RUN" != "--dry-run" ]]; then
        PCAIR_JOB=$(eval "$cmd")
        DEPEND_PCAIR="--dependency=afterok:${PCAIR_JOB}"
        echo "  Job ID: ${PCAIR_JOB}"
    else
        echo "  [DRY-RUN] $cmd"
    fi
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
        echo "  ${trait} (${model}) / ${stratum}: 22 chromosomes"
        cmd="sbatch --parsable ${DEPEND_PCAIR} ${SCRIPT_DIR}/submit_gwas_array.sh ${stratum} ${trait} ${model}"
        if [[ "$DRY_RUN" != "--dry-run" ]]; then
            GWAS_JOBS["${trait}.${stratum}"]=$(eval "$cmd")
            echo "    Job ID: ${GWAS_JOBS["${trait}.${stratum}"]}"
        else
            echo "    [DRY-RUN] $cmd"
        fi
    done
done

# ============================================================================
# STEP 2: Combine chromosomes per trait x stratum, meta-analyse strata -> POOLED
# ============================================================================
echo ""
echo "STEP 2: Submitting combine + meta-analysis jobs (one per trait)..."

declare -A COMBINE_JOBS
STRATA_CSV=$(IFS=,; echo "${STRATA[*]}")

for trait in "${TRAITS[@]}"; do
    dep=""
    for stratum in "${STRATA[@]}"; do
        j="${GWAS_JOBS["${trait}.${stratum}"]:-}"
        [[ -n "$j" ]] && dep="${dep}:${j}"
    done
    DEPEND_GWAS=""; [[ -n "$dep" ]] && DEPEND_GWAS="--dependency=afterok${dep}"

    COMBINE_CMD="sbatch --parsable ${DEPEND_GWAS} << 'EOF'
#!/bin/bash
#SBATCH --job-name=combine_${trait}
#SBATCH --partition=normal
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/combine_${trait}_%j.out

SCRIPT_DIR=\"${SCRIPT_DIR}/..\"
GWAS_DIR=\"\${SCRIPT_DIR}/results/gwas/${trait}\"

# Concatenate chromosomes per stratum (genomic order files -> one sumstats per stratum)
for stratum in ${STRATA[*]}; do
    files=(\${GWAS_DIR}/\${stratum}/chr*.tractor_genesis.genomic_order.tsv.gz)
    [[ -e \"\${files[0]}\" ]] || { echo \"No results for \${stratum} (skipped: N < 30?)\"; continue; }
    echo \"Combining ${trait} / \${stratum}...\"
    { zcat \"\${files[0]}\" | head -1; for f in \"\${files[@]}\"; do zcat \"\$f\" | tail -n +2; done; } | \\
        gzip > \${GWAS_DIR}/\${stratum}.sumstats.gz
done

# Meta-analyse across strata (same ancestry background pooled across strata)
echo \"Meta-analysing ${trait} across strata...\"
Rscript \${SCRIPT_DIR}/bin/meta_analyze_strata.R \\
    --input_dir \${GWAS_DIR} \\
    --strata ${STRATA_CSV} \\
    --output \${GWAS_DIR}/POOLED.sumstats.gz \\
    --threads 8
EOF"

    if [[ "$DRY_RUN" != "--dry-run" ]]; then
        COMBINE_JOBS[$trait]=$(eval "$COMBINE_CMD")
        echo "  ${trait} combine Job ID: ${COMBINE_JOBS[$trait]}"
    else
        echo "  [DRY-RUN] Submit combine job for ${trait}"
    fi
done

# Single dependency handle for downstream steps (all traits combined)
COMBINE_JOB=""
for trait in "${TRAITS[@]}"; do
    [[ -n "${COMBINE_JOBS[$trait]:-}" ]] && COMBINE_JOB="${COMBINE_JOB:+${COMBINE_JOB}:}${COMBINE_JOBS[$trait]}"
done

# ============================================================================
# STEP 3: PRS for each trait (depends on GWAS combine)
# ============================================================================
echo ""
echo "STEP 3: Submitting PRS jobs..."

if [[ -n "${COMBINE_JOB:-}" ]]; then
    DEPEND_PRS="--dependency=afterok:${COMBINE_JOB}"
else
    DEPEND_PRS=""
fi

declare -A PRS_JOBS

for trait in "${TRAITS[@]}"; do
    echo "  Submitting PRS for ${trait} (6 methods)..."

    cmd="sbatch --parsable ${DEPEND_PRS} ${SCRIPT_DIR}/submit_prs_array.sh ${trait}"
    if [[ "$DRY_RUN" != "--dry-run" ]]; then
        PRS_JOBS[$trait]=$(eval "$cmd")
        echo "    Job ID: ${PRS_JOBS[$trait]}"
    else
        echo "    [DRY-RUN] $cmd"
    fi
done

# ============================================================================
# STEP 4: Colocalization (parallel with PRS)
# ============================================================================
echo ""
echo "STEP 4: Submitting colocalization jobs..."

# Colocalization can start after GWAS combine
COLOC_CMD="sbatch --parsable ${DEPEND_PRS:-} << 'EOF'
#!/bin/bash
#SBATCH --job-name=colocalization
#SBATCH --array=1-6
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/coloc_%A_%a.out

SCRIPT_DIR=\"${SCRIPT_DIR}/..\"

# Map array ID to QTL type
declare -a QTL_TYPES=(\"eqtl\" \"sqtl\" \"pqtl\" \"mqtl\" \"caqtl\" \"hqtl\")
QTL_TYPE=\"\${QTL_TYPES[\$((SLURM_ARRAY_TASK_ID - 1))]}\"

echo \"Running colocalization for \${QTL_TYPE}...\"

# Run all methods in parallel (coloc.susie, HyPrColoc, OPERA)
Rscript \${SCRIPT_DIR}/bin/run_colocalization.R \\
    --gwas \${SCRIPT_DIR}/results/gwas/POOLED.sumstats.gz \\
    --qtl_dir \${SCRIPT_DIR}/data/qtl \\
    --qtl_types \${QTL_TYPE} \\
    --methods coloc_susie,hyprcoloc,opera \\
    --run_opera true \\
    --output_prefix \${SCRIPT_DIR}/results/coloc/\${QTL_TYPE} \\
    --threads \${SLURM_CPUS_PER_TASK} \\
    --verbose
EOF"

if [[ "$DRY_RUN" != "--dry-run" ]]; then
    COLOC_JOB=$(eval "$COLOC_CMD")
    echo "  Colocalization Job ID: ${COLOC_JOB}"
else
    echo "  [DRY-RUN] Submit colocalization array job"
fi

# ============================================================================
# STEP 5: Validation and comparison (depends on PRS)
# ============================================================================
echo ""
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

VALIDATE_CMD="sbatch --parsable ${DEPEND_VAL} << 'EOF'
#!/bin/bash
#SBATCH --job-name=validation
#SBATCH --partition=normal
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/validation_%j.out

SCRIPT_DIR=\"${SCRIPT_DIR}/..\"

echo \"Comparing PRS methods across strata...\"

# Combine all validation results
for trait in OS relapse MRD; do
    cat \${SCRIPT_DIR}/results/prs/\${trait}/*.validation.tsv | \\
        awk 'NR==1 || !/^method/' > \\
        \${SCRIPT_DIR}/results/prs/\${trait}_all_validation.tsv
done

# Generate summary report
echo \"Generating summary report...\"
# ... report generation code
EOF"

if [[ "$DRY_RUN" != "--dry-run" ]]; then
    VALIDATE_JOB=$(eval "$VALIDATE_CMD")
    echo "  Validation Job ID: ${VALIDATE_JOB}"
else
    echo "  [DRY-RUN] Submit validation job"
fi

# ============================================================================
# STEP 6: GxG interaction tests on GWAS hits (END OF PIPELINE)
# ============================================================================
# Needs only the combined GWAS sumstats (not PRS/coloc), so it runs in
# parallel with steps 3-5. One array task per stratum, per trait.
echo ""
echo "STEP 6: Submitting GxG interaction jobs..."

declare -A GXG_JOBS

for trait in "${TRAITS[@]}"; do
    model=$(trait_model "$trait")
    echo "  Submitting GxG for ${trait} (${model}; POOLED + 7 strata)..."

    cmd="sbatch --parsable ${DEPEND_PRS} ${SCRIPT_DIR}/submit_gxg_array.sh ${trait} ${model}"
    if [[ "$DRY_RUN" != "--dry-run" ]]; then
        GXG_JOBS[$trait]=$(eval "$cmd")
        echo "    Job ID: ${GXG_JOBS[$trait]}"
    else
        echo "    [DRY-RUN] $cmd"
    fi
done

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "=============================================="
echo "Pipeline Submitted"
echo "=============================================="
echo "GWAS jobs: ${#STRATA[@]} strata × 22 chromosomes"
echo "PRS jobs: ${#TRAITS[@]} traits × 6 methods"
echo "Colocalization: ${#QTL_TYPES[@]} QTL types × 3 methods"
echo "GxG interaction: ${#TRAITS[@]} traits × 8 strata (+ ancestry heterogeneity)"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Cancel all: scancel -u \$USER"
echo "=============================================="
