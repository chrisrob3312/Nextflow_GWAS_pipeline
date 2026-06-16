#!/bin/bash
# Shell wrapper for Cohort-Specific LD Calculation
# Maintains parity with modules/local/cohort_ld/main.nf

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${SCRIPT_DIR}/../bin"
LD_SCRIPT="${BIN_DIR}/calculate_cohort_ld.R"

# Default values
GDS_FILE=""
ANCESTRY_FILE=""
ANCESTRY_COL="ancestry"
SAMPLE_COL="sample_id"
GROUPS=""
COMBINE_GROUPS=""
REGION=""
REGIONS_FILE=""
WINDOW_KB=1000
R2_THRESHOLD=0.01
MAF_FILTER=0.01
FORMAT="ldstore"
OUTPUT_PREFIX="cohort_ld"
OUTPUT_DIR="."
THREADS=4
VERBOSE=""

usage() {
    cat <<EOF
Cohort-Specific LD Calculation

Calculates LD matrices from cohort data, stratified by user-defined groups.
Flexible grouping - not limited to specific ancestry inference tools.

Usage: $(basename "$0") [OPTIONS]

Required:
    -g, --gds FILE          GDS file with genotypes

Ancestry/Stratification:
    -a, --ancestry_file FILE  TSV with sample ancestry assignments
    --ancestry_col COL        Column name for ancestry/group (default: ancestry)
    --sample_col COL          Column name for sample ID (default: sample_id)
    --groups LIST             Comma-separated groups to process (e.g., EUR,AFR,LATINO)
                              Use "ALL" to pool all samples
    --combine_groups SPEC     Combine groups, format: "NEW=OLD1+OLD2;NEW2=OLD3+OLD4"
                              Example: "LATINO=LAT1+LAT2"

Region:
    --region REGION           Single region: "chr:start-end" or "chr"
    --regions_file FILE       BED file with regions

LD Parameters:
    --window_kb INT           LD window in kb (default: 1000)
    --r2_threshold NUM        Minimum r2 to store (default: 0.01)
    --maf_filter NUM          MAF filter (default: 0.01)

Output:
    -o, --output_prefix STR   Output prefix (default: cohort_ld)
    -d, --output_dir DIR      Output directory (default: .)
    --format FORMAT           Output format: ldstore, coloc, prs_csx, matrix (default: ldstore)

Runtime:
    --threads INT             Number of threads (default: 4)
    -v, --verbose             Verbose output
    -h, --help                Show this help

Examples:
    # Calculate LD by GRAF-ANC ancestry groups
    $(basename "$0") -g data.gds -a ancestry.tsv --ancestry_col grafanc_category \\
        --groups "EUR,AFR,LATINO,AAC"

    # Combine LAT1+LAT2 into LATINO
    $(basename "$0") -g data.gds -a ancestry.tsv --groups "EUR,AFR,LATINO" \\
        --combine_groups "LATINO=LAT1+LAT2"

    # Single region for coloc
    $(basename "$0") -g data.gds -a ancestry.tsv --groups "EUR,AFR" \\
        --region "6:28000000-34000000" --format coloc

    # Pool all samples (no stratification)
    $(basename "$0") -g data.gds --groups "ALL"

EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -g|--gds) GDS_FILE="$2"; shift 2 ;;
        -a|--ancestry_file) ANCESTRY_FILE="$2"; shift 2 ;;
        --ancestry_col) ANCESTRY_COL="$2"; shift 2 ;;
        --sample_col) SAMPLE_COL="$2"; shift 2 ;;
        --groups) GROUPS="$2"; shift 2 ;;
        --combine_groups) COMBINE_GROUPS="$2"; shift 2 ;;
        --region) REGION="$2"; shift 2 ;;
        --regions_file) REGIONS_FILE="$2"; shift 2 ;;
        --window_kb) WINDOW_KB="$2"; shift 2 ;;
        --r2_threshold) R2_THRESHOLD="$2"; shift 2 ;;
        --maf_filter) MAF_FILTER="$2"; shift 2 ;;
        --format) FORMAT="$2"; shift 2 ;;
        -o|--output_prefix) OUTPUT_PREFIX="$2"; shift 2 ;;
        -d|--output_dir) OUTPUT_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        -v|--verbose) VERBOSE="--verbose"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Error: Unknown option $1"; usage; exit 1 ;;
    esac
done

# Validate
if [[ -z "$GDS_FILE" ]]; then
    echo "Error: --gds is required"
    usage
    exit 1
fi

if [[ ! -f "$GDS_FILE" ]]; then
    echo "Error: GDS file not found: $GDS_FILE"
    exit 1
fi

if [[ ! -f "$LD_SCRIPT" ]]; then
    echo "Error: Core R script not found: $LD_SCRIPT"
    exit 1
fi

# Build command
mkdir -p "$OUTPUT_DIR"
CMD="Rscript ${LD_SCRIPT}"
CMD+=" --gds ${GDS_FILE}"

if [[ -n "$ANCESTRY_FILE" && -f "$ANCESTRY_FILE" ]]; then
    CMD+=" --ancestry_file ${ANCESTRY_FILE}"
    CMD+=" --ancestry_col ${ANCESTRY_COL}"
    CMD+=" --sample_col ${SAMPLE_COL}"
fi

if [[ -n "$GROUPS" ]]; then
    CMD+=" --groups ${GROUPS}"
fi

if [[ -n "$COMBINE_GROUPS" ]]; then
    CMD+=" --combine_groups '${COMBINE_GROUPS}'"
fi

if [[ -n "$REGION" ]]; then
    CMD+=" --region ${REGION}"
fi

if [[ -n "$REGIONS_FILE" && -f "$REGIONS_FILE" ]]; then
    CMD+=" --regions_file ${REGIONS_FILE}"
fi

CMD+=" --window_kb ${WINDOW_KB}"
CMD+=" --r2_threshold ${R2_THRESHOLD}"
CMD+=" --maf_filter ${MAF_FILTER}"
CMD+=" --format ${FORMAT}"
CMD+=" --threads ${THREADS}"
CMD+=" --output_prefix ${OUTPUT_DIR}/${OUTPUT_PREFIX}"

if [[ -n "$VERBOSE" ]]; then
    CMD+=" ${VERBOSE}"
fi

# Execute
echo "=============================================="
echo "Cohort-Specific LD Calculation"
echo "=============================================="
echo "GDS: $GDS_FILE"
echo "Groups: ${GROUPS:-ALL}"
echo "Output: ${OUTPUT_DIR}/${OUTPUT_PREFIX}"
echo "=============================================="
echo ""

LOG_FILE="${OUTPUT_DIR}/${OUTPUT_PREFIX}.log"
eval $CMD 2>&1 | tee "$LOG_FILE"

EXIT_CODE=${PIPESTATUS[0]}

if [[ $EXIT_CODE -eq 0 ]]; then
    echo ""
    echo "LD calculation completed successfully"
    ls -la "${OUTPUT_DIR}/${OUTPUT_PREFIX}"* 2>/dev/null || true
else
    echo "LD calculation failed with exit code: $EXIT_CODE"
    exit $EXIT_CODE
fi
