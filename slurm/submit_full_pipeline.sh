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

# Ancestry strata to analyze
# LAT1/LAT2 run separately (>60% of cohort)
# EUR, AAC run separately
# OTHER = EAS + SAS + small groups (<30)
STRATA=("EUR" "AAC" "LAT1" "LAT2" "OTHER")

# Traits to analyze
TRAITS=("OS" "relapse" "MRD")

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
# STEP 1: GWAS for each ancestry stratum (parallelized by chromosome)
# ============================================================================
echo ""
echo "STEP 1: Submitting GWAS jobs..."
declare -A GWAS_JOBS

for stratum in "${STRATA[@]}"; do
    echo "  Submitting GWAS for ${stratum} (22 chromosomes)..."

    cmd="sbatch --parsable ${SCRIPT_DIR}/submit_gwas_array.sh ${stratum}"
    if [[ "$DRY_RUN" != "--dry-run" ]]; then
        GWAS_JOBS[$stratum]=$(eval "$cmd")
        echo "    Job ID: ${GWAS_JOBS[$stratum]}"
    else
        echo "    [DRY-RUN] $cmd"
    fi
done

# ============================================================================
# STEP 2: Combine chromosomes and meta-analyze (depends on GWAS)
# ============================================================================
echo ""
echo "STEP 2: Submitting chromosome combine jobs..."

# Build dependency string
GWAS_DEP=""
for stratum in "${STRATA[@]}"; do
    if [[ -n "${GWAS_JOBS[$stratum]:-}" ]]; then
        GWAS_DEP="${GWAS_DEP}:${GWAS_JOBS[$stratum]}"
    fi
done

if [[ -n "$GWAS_DEP" ]]; then
    DEPEND_GWAS="--dependency=afterok${GWAS_DEP}"
else
    DEPEND_GWAS=""
fi

# Submit combine job
COMBINE_CMD="sbatch --parsable ${DEPEND_GWAS} << 'EOF'
#!/bin/bash
#SBATCH --job-name=combine_gwas
#SBATCH --partition=normal
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/combine_gwas_%j.out

SCRIPT_DIR=\"${SCRIPT_DIR}/..\"

# Combine chromosomes for each stratum
for stratum in EUR AAC LAT1 LAT2 OTHER; do
    echo \"Combining \${stratum}...\"
    cat \${SCRIPT_DIR}/results/gwas/\${stratum}/chr*.sumstats.tsv | \\
        awk 'NR==1 || !/^CHR/' | sort -k1,1 -k2,2n | \\
        gzip > \${SCRIPT_DIR}/results/gwas/\${stratum}.sumstats.gz
done

# Meta-analyze OTHER with primary strata for pooled results
echo \"Meta-analyzing for pooled results...\"
Rscript \${SCRIPT_DIR}/bin/meta_analyze_strata.R \\
    --input_dir \${SCRIPT_DIR}/results/gwas \\
    --strata EUR,AAC,LAT1,LAT2,OTHER \\
    --output \${SCRIPT_DIR}/results/gwas/POOLED.sumstats.gz
EOF"

if [[ "$DRY_RUN" != "--dry-run" ]]; then
    COMBINE_JOB=$(eval "$COMBINE_CMD")
    echo "  Combine Job ID: ${COMBINE_JOB}"
else
    echo "  [DRY-RUN] Submit combine job with dependencies"
fi

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
# Summary
# ============================================================================
echo ""
echo "=============================================="
echo "Pipeline Submitted"
echo "=============================================="
echo "GWAS jobs: ${#STRATA[@]} strata × 22 chromosomes"
echo "PRS jobs: ${#TRAITS[@]} traits × 6 methods"
echo "Colocalization: ${#QTL_TYPES[@]} QTL types × 3 methods"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Cancel all: scancel -u \$USER"
echo "=============================================="
