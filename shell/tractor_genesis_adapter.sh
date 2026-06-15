#!/bin/bash
# Shell wrapper for Tractor → GENESIS Adapter
# Maintains parity with modules/local/tractor_genesis_adapter/main.nf
#
# DUAL STRUCTURE NOTE:
# This shell script calls the same core R script (bin/tractor_genesis_adapter.R)
# that the Nextflow module uses. Changes to analysis logic should be made in
# bin/tractor_genesis_adapter.R only. This wrapper handles I/O and logging.

set -euo pipefail

# ============================================
# Configuration
# ============================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${SCRIPT_DIR}/../bin"
TRACTOR_GENESIS_SCRIPT="${BIN_DIR}/tractor_genesis_adapter.R"

# Default values
GDS_FILE=""
WEIGHTS_FILE=""
PHENOTYPE_FILE=""
KINSHIP_FILE=""
TRAIT=""
MODEL="binary"
ANCESTRIES="EUR,AFR"
COVARIATES=""
TIME_COL="time"
EVENT_COL="event"
OUTPUT_PREFIX="tractor_genesis"
OUTPUT_DIR="."
VERBOSE=""

# ============================================
# Usage
# ============================================
usage() {
    cat <<EOF
Tractor → GENESIS Adapter Shell Wrapper

Usage: $(basename "$0") [OPTIONS]

Required:
    -g, --gds FILE          Tractor-prepared GDS file
    -p, --phenotype FILE    Phenotype file (TSV)
    -t, --trait NAME        Trait column name

Optional:
    -w, --weights FILE      Ancestry weights RDS file
    -k, --kinship FILE      Kinship matrix RDS file
    -m, --model TYPE        Model type: survival, binary, quantitative (default: binary)
    -a, --ancestries LIST   Comma-separated ancestries (default: EUR,AFR)
    -c, --covariates LIST   Comma-separated covariate columns
    --time_col NAME         Time column for survival (default: time)
    --event_col NAME        Event column for survival (default: event)
    -o, --output_prefix STR Output file prefix (default: tractor_genesis)
    -d, --output_dir DIR    Output directory (default: .)
    -v, --verbose           Enable verbose output
    -h, --help              Show this help

Model Types:
    survival      Cox proportional hazards (requires --time_col, --event_col)
    binary        Logistic mixed model for case-control
    quantitative  Linear mixed model for continuous traits

Ancestry Options:
    2-way: EUR,AFR (African American admixture)
    3-way: EUR,AFR,AMR (Latino admixture)

Examples:
    # Binary trait analysis for African American cohort
    $(basename "$0") -g data.gds -p pheno.tsv -t relapse -m binary -a EUR,AFR

    # Survival analysis for Latino cohort
    $(basename "$0") -g data.gds -p pheno.tsv -t os -m survival \\
        -a EUR,AFR,AMR --time_col os_time --event_col os_event

    # Quantitative MRD analysis with kinship adjustment
    $(basename "$0") -g data.gds -p pheno.tsv -t mrd_level -m quantitative \\
        -k kinship.rds -c "age,sex,PC1,PC2,PC3,PC4,PC5"

EOF
}

# ============================================
# Parse arguments
# ============================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        -g|--gds)
            GDS_FILE="$2"
            shift 2
            ;;
        -w|--weights)
            WEIGHTS_FILE="$2"
            shift 2
            ;;
        -p|--phenotype)
            PHENOTYPE_FILE="$2"
            shift 2
            ;;
        -k|--kinship)
            KINSHIP_FILE="$2"
            shift 2
            ;;
        -t|--trait)
            TRAIT="$2"
            shift 2
            ;;
        -m|--model)
            MODEL="$2"
            shift 2
            ;;
        -a|--ancestries)
            ANCESTRIES="$2"
            shift 2
            ;;
        -c|--covariates)
            COVARIATES="$2"
            shift 2
            ;;
        --time_col)
            TIME_COL="$2"
            shift 2
            ;;
        --event_col)
            EVENT_COL="$2"
            shift 2
            ;;
        -o|--output_prefix)
            OUTPUT_PREFIX="$2"
            shift 2
            ;;
        -d|--output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE="--verbose"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            exit 1
            ;;
    esac
done

# ============================================
# Validate inputs
# ============================================
if [[ -z "$GDS_FILE" || -z "$PHENOTYPE_FILE" || -z "$TRAIT" ]]; then
    echo "Error: Required arguments missing"
    usage
    exit 1
fi

if [[ ! -f "$GDS_FILE" ]]; then
    echo "Error: GDS file not found: $GDS_FILE"
    exit 1
fi

if [[ ! -f "$PHENOTYPE_FILE" ]]; then
    echo "Error: Phenotype file not found: $PHENOTYPE_FILE"
    exit 1
fi

if [[ ! -f "$TRACTOR_GENESIS_SCRIPT" ]]; then
    echo "Error: Core R script not found: $TRACTOR_GENESIS_SCRIPT"
    exit 1
fi

# Validate model type
case "$MODEL" in
    survival|binary|quantitative)
        ;;
    *)
        echo "Error: Invalid model type: $MODEL"
        echo "Valid options: survival, binary, quantitative"
        exit 1
        ;;
esac

# Create output directory
mkdir -p "$OUTPUT_DIR"

# ============================================
# Build command
# ============================================
CMD="Rscript ${TRACTOR_GENESIS_SCRIPT}"
CMD+=" --gds ${GDS_FILE}"
CMD+=" --phenotype ${PHENOTYPE_FILE}"
CMD+=" --trait ${TRAIT}"
CMD+=" --model ${MODEL}"
CMD+=" --ancestries ${ANCESTRIES}"
CMD+=" --output_prefix ${OUTPUT_DIR}/${OUTPUT_PREFIX}"

if [[ -n "$WEIGHTS_FILE" && -f "$WEIGHTS_FILE" ]]; then
    CMD+=" --weights ${WEIGHTS_FILE}"
fi

if [[ -n "$KINSHIP_FILE" && -f "$KINSHIP_FILE" ]]; then
    CMD+=" --kinship ${KINSHIP_FILE}"
fi

if [[ -n "$COVARIATES" ]]; then
    CMD+=" --covariates ${COVARIATES}"
fi

if [[ "$MODEL" == "survival" ]]; then
    CMD+=" --time_col ${TIME_COL}"
    CMD+=" --event_col ${EVENT_COL}"
fi

if [[ -n "$VERBOSE" ]]; then
    CMD+=" ${VERBOSE}"
fi

# ============================================
# Execute
# ============================================
echo "=============================================="
echo "Tractor → GENESIS Adapter"
echo "=============================================="
echo "GDS file: $GDS_FILE"
echo "Phenotype: $PHENOTYPE_FILE"
echo "Trait: $TRAIT"
echo "Model: $MODEL"
echo "Ancestries: $ANCESTRIES"
echo "Output: ${OUTPUT_DIR}/${OUTPUT_PREFIX}"
echo "=============================================="
echo ""

# Log file
LOG_FILE="${OUTPUT_DIR}/${OUTPUT_PREFIX}.log"
echo "Running analysis (log: $LOG_FILE)..."

# Execute with logging
$CMD 2>&1 | tee "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]}

if [[ $EXIT_CODE -eq 0 ]]; then
    echo ""
    echo "=============================================="
    echo "Analysis completed successfully"
    echo "=============================================="
    echo "Output files:"
    ls -la "${OUTPUT_DIR}/${OUTPUT_PREFIX}"* 2>/dev/null || true
else
    echo ""
    echo "=============================================="
    echo "Analysis failed with exit code: $EXIT_CODE"
    echo "=============================================="
    exit $EXIT_CODE
fi
