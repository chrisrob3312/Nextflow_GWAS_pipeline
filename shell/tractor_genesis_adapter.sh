#!/bin/bash
# ============================================================================
# TractorGENESIS Shell Wrapper
# Local ancestry-aware mixed model GWAS
# ============================================================================
# Maintains parity with modules/local/tractor_genesis_adapter/main.nf
#
# SUPPORTS:
#   - 2-way admixture (e.g., AAC: AFR-EUR)
#   - 3-way admixture (e.g., Latino: EUR-AFR-AMR)
#   - N-way admixture (any ancestry backgrounds)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${SCRIPT_DIR}/../bin"
TRACTOR_GENESIS="${BIN_DIR}/tractor_genesis_adapter.R"

usage() {
    cat <<EOF
TractorGENESIS: Local Ancestry-Aware Mixed Model GWAS

Combines Tractor's local ancestry decomposition with GENESIS mixed models.
Supports 2-way, 3-way, or N-way admixture with any ancestry backgrounds.

Usage: $(basename "$0") [OPTIONS]

REQUIRED:
    --tractor_prefix PREFIX   Prefix for Tractor output files
                              (expects .hapcount.{anc}.txt.gz and .dosage.{anc}.txt.gz)
    -p, --phenotype FILE      Phenotype file (TSV with sample_id column)
    -t, --trait NAME          Trait column name
    -a, --ancestries LIST     Comma-separated ancestries matching Tractor coding
                              e.g., "EUR,AFR" (2-way) or "EUR,AFR,AMR" (3-way)

OPTIONAL INPUT:
    --msp FILE                RFMix MSP file (alternative to Tractor files)
    --gds FILE                GDS file (required with --msp)
    -k, --kinship FILE        Kinship matrix RDS (for relatedness)

MODEL:
    -m, --model TYPE          Model type: binary, quantitative, survival
                              [default: binary]
    --time_col NAME           Time column for survival [default: time]
    --event_col NAME          Event column for survival [default: event]
    --ref_ancestry NAME       Reference ancestry for LA terms [default: first]

COVARIATES:
    -c, --covariates LIST     Comma-separated covariate columns
    --global_ancestry_prefix  Prefix for global ancestry columns (e.g., "prop_")

OUTPUT:
    -o, --output_prefix STR   Output prefix [default: tractor_genesis]
    -d, --output_dir DIR      Output directory [default: .]

RUNTIME:
    --threads INT             Number of threads [default: 1]
    -v, --verbose             Verbose output
    -h, --help                Show this help

EXAMPLES:
    # 2-way admixture (African American: AFR-EUR)
    $(basename "$0") --tractor_prefix data/aac \\
        -p phenotypes.tsv -t relapse -a "EUR,AFR" \\
        -k kinship.rds -c "age,sex,PC1,PC2,PC3"

    # 3-way admixture (Latino: EUR-AFR-AMR)
    $(basename "$0") --tractor_prefix data/latino \\
        -p phenotypes.tsv -t os -m survival \\
        -a "EUR,AFR,AMR" --time_col os_time --event_col os_event

    # Custom ancestry backgrounds (e.g., Southeast Asian)
    $(basename "$0") --tractor_prefix data/sea \\
        -p phenotypes.tsv -t mrd -m quantitative \\
        -a "EAS,SAS,EUR"

OUTPUT INTERPRETATION:
    P_JOINT     : Does SNP have ANY effect (across all ancestries)?
    BETA_{anc}  : Effect size when allele is on {anc} haplotype
    P_{anc}     : Is there an effect on {anc} background specifically?
    P_HET       : Do effects DIFFER across ancestries?
    I2          : Heterogeneity magnitude (0-100%)

EOF
}

# Defaults
TRACTOR_PREFIX=""
MSP_FILE=""
GDS_FILE=""
PHENOTYPE=""
KINSHIP=""
TRAIT=""
MODEL="binary"
ANCESTRIES=""
REF_ANCESTRY=""
COVARIATES=""
GLOBAL_ANC_PREFIX=""
TIME_COL="time"
EVENT_COL="event"
OUTPUT_PREFIX="tractor_genesis"
OUTPUT_DIR="."
THREADS=1
VERBOSE=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tractor_prefix) TRACTOR_PREFIX="$2"; shift 2 ;;
        --msp) MSP_FILE="$2"; shift 2 ;;
        --gds) GDS_FILE="$2"; shift 2 ;;
        -p|--phenotype) PHENOTYPE="$2"; shift 2 ;;
        -k|--kinship) KINSHIP="$2"; shift 2 ;;
        -t|--trait) TRAIT="$2"; shift 2 ;;
        -m|--model) MODEL="$2"; shift 2 ;;
        -a|--ancestries) ANCESTRIES="$2"; shift 2 ;;
        --ref_ancestry) REF_ANCESTRY="$2"; shift 2 ;;
        -c|--covariates) COVARIATES="$2"; shift 2 ;;
        --global_ancestry_prefix) GLOBAL_ANC_PREFIX="$2"; shift 2 ;;
        --time_col) TIME_COL="$2"; shift 2 ;;
        --event_col) EVENT_COL="$2"; shift 2 ;;
        -o|--output_prefix) OUTPUT_PREFIX="$2"; shift 2 ;;
        -d|--output_dir) OUTPUT_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        -v|--verbose) VERBOSE="--verbose"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Error: Unknown option $1"; usage; exit 1 ;;
    esac
done

# Validate required
if [[ -z "$PHENOTYPE" || -z "$TRAIT" || -z "$ANCESTRIES" ]]; then
    echo "Error: --phenotype, --trait, and --ancestries are required"
    usage
    exit 1
fi

if [[ -z "$TRACTOR_PREFIX" && -z "$MSP_FILE" ]]; then
    echo "Error: --tractor_prefix or --msp is required"
    usage
    exit 1
fi

if [[ ! -f "$PHENOTYPE" ]]; then
    echo "Error: Phenotype file not found: $PHENOTYPE"
    exit 1
fi

if [[ ! -f "$TRACTOR_GENESIS" ]]; then
    echo "Error: Core R script not found: $TRACTOR_GENESIS"
    exit 1
fi

# Build command
mkdir -p "$OUTPUT_DIR"

CMD="Rscript ${TRACTOR_GENESIS}"
CMD+=" --phenotype ${PHENOTYPE}"
CMD+=" --trait ${TRAIT}"
CMD+=" --ancestries ${ANCESTRIES}"
CMD+=" --model ${MODEL}"
CMD+=" --output_prefix ${OUTPUT_DIR}/${OUTPUT_PREFIX}"

if [[ -n "$TRACTOR_PREFIX" ]]; then
    CMD+=" --tractor_prefix ${TRACTOR_PREFIX}"
fi

if [[ -n "$MSP_FILE" ]]; then
    CMD+=" --msp ${MSP_FILE}"
fi

if [[ -n "$GDS_FILE" && -f "$GDS_FILE" ]]; then
    CMD+=" --gds ${GDS_FILE}"
fi

if [[ -n "$KINSHIP" && -f "$KINSHIP" ]]; then
    CMD+=" --kinship ${KINSHIP}"
fi

if [[ -n "$REF_ANCESTRY" ]]; then
    CMD+=" --ref_ancestry ${REF_ANCESTRY}"
fi

if [[ -n "$COVARIATES" ]]; then
    CMD+=" --covariates ${COVARIATES}"
fi

if [[ -n "$GLOBAL_ANC_PREFIX" ]]; then
    CMD+=" --global_ancestry_prefix ${GLOBAL_ANC_PREFIX}"
fi

if [[ "$MODEL" == "survival" ]]; then
    CMD+=" --time_col ${TIME_COL} --event_col ${EVENT_COL}"
fi

CMD+=" --threads ${THREADS}"

if [[ -n "$VERBOSE" ]]; then
    CMD+=" ${VERBOSE}"
fi

# Execute
echo "============================================================================"
echo "TractorGENESIS: Local Ancestry-Aware Mixed Model GWAS"
echo "============================================================================"
echo "Ancestries: $ANCESTRIES"
echo "Trait: $TRAIT ($MODEL)"
echo "Output: ${OUTPUT_DIR}/${OUTPUT_PREFIX}"
echo "============================================================================"
echo ""

LOG_FILE="${OUTPUT_DIR}/${OUTPUT_PREFIX}.log"
$CMD 2>&1 | tee "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]}

if [[ $EXIT_CODE -eq 0 ]]; then
    echo ""
    echo "Analysis completed successfully"
    echo "Output files:"
    ls -la "${OUTPUT_DIR}/${OUTPUT_PREFIX}"* 2>/dev/null || true
else
    echo ""
    echo "Analysis failed with exit code: $EXIT_CODE"
    exit $EXIT_CODE
fi
