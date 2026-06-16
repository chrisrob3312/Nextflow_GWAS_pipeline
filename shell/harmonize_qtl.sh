#!/bin/bash
# QTL Harmonization Shell Wrapper
# Maintains parity with modules/local/qtl_harmonization/main.nf

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${SCRIPT_DIR}/../bin"
QTL_SCRIPT="${BIN_DIR}/harmonize_qtl.R"

usage() {
    cat <<EOF
QTL Harmonization Pipeline

Harmonizes QTL summary statistics from multiple sources into unified format.
Supports liftover, deduplication, filtering by tissue/ancestry/type.

Usage: $(basename "$0") [OPTIONS]

INPUT (one required):
    -i, --input FILE/DIR      Single QTL file or directory of files
    --manifest FILE           TSV manifest with columns: path,source,tissue,ancestry,qtl_type

FORMAT:
    --format FORMAT           Input format: auto, eqtl_catalogue, gtex, eqtlgen, mesa, sushie, custom
                              [default: auto]

GENOME BUILD:
    --input_build BUILD       Input genome build [default: GRCh38]
    --output_build BUILD      Output genome build [default: GRCh38]
    --chain_file FILE         Chain file for liftover

FILTERING:
    --tissues LIST            Comma-separated tissues to include
    --ancestries LIST         Comma-separated ancestries to include
    --qtl_types LIST          Comma-separated QTL types: eqtl,sqtl,pqtl,mqtl,caqtl
    --genes FILE              File with gene IDs to include
    --region REGION           Genomic region: chr:start-end
    --pvalue_threshold NUM    P-value threshold [default: 1]

OUTPUT:
    -o, --output PREFIX       Output prefix [default: harmonized_qtl]
    -d, --output_dir DIR      Output directory [default: .]
    --split_by VAR            Split output by: tissue, ancestry, qtl_type, gene, chr

DEDUPLICATION:
    --dedup_strategy STR      Deduplication: best_p, first, none [default: best_p]

CUSTOM FORMAT COLUMNS:
    --col_variant COL         Variant ID column
    --col_chr COL             Chromosome column
    --col_pos COL             Position column
    --col_gene COL            Gene/feature column
    --col_beta COL            Effect size column
    --col_pvalue COL          P-value column

EXAMPLES:
    # Harmonize GTEx eQTLs for blood tissues
    $(basename "$0") -i gtex_v8_eqtl.txt.gz --format gtex \\
        --tissues "Whole_Blood,Cells_EBV-transformed_lymphocytes"

    # Process multiple QTL files via manifest
    $(basename "$0") --manifest qtl_manifest.tsv --output_build GRCh38

    # Filter to specific region for colocalization
    $(basename "$0") -i harmonized_qtl.tsv.gz --region "6:28000000-34000000"

EOF
}

# Defaults
INPUT=""
MANIFEST=""
FORMAT="auto"
INPUT_BUILD="GRCh38"
OUTPUT_BUILD="GRCh38"
CHAIN_FILE=""
TISSUES=""
ANCESTRIES=""
QTL_TYPES=""
GENES=""
REGION=""
PVALUE_THRESHOLD=1
OUTPUT_PREFIX="harmonized_qtl"
OUTPUT_DIR="."
SPLIT_BY=""
DEDUP_STRATEGY="best_p"
VERBOSE=""

# Custom columns
COL_VARIANT=""
COL_CHR=""
COL_POS=""
COL_GENE=""
COL_BETA=""
COL_PVALUE=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input) INPUT="$2"; shift 2 ;;
        --manifest) MANIFEST="$2"; shift 2 ;;
        --format) FORMAT="$2"; shift 2 ;;
        --input_build) INPUT_BUILD="$2"; shift 2 ;;
        --output_build) OUTPUT_BUILD="$2"; shift 2 ;;
        --chain_file) CHAIN_FILE="$2"; shift 2 ;;
        --tissues) TISSUES="$2"; shift 2 ;;
        --ancestries) ANCESTRIES="$2"; shift 2 ;;
        --qtl_types) QTL_TYPES="$2"; shift 2 ;;
        --genes) GENES="$2"; shift 2 ;;
        --region) REGION="$2"; shift 2 ;;
        --pvalue_threshold) PVALUE_THRESHOLD="$2"; shift 2 ;;
        -o|--output) OUTPUT_PREFIX="$2"; shift 2 ;;
        -d|--output_dir) OUTPUT_DIR="$2"; shift 2 ;;
        --split_by) SPLIT_BY="$2"; shift 2 ;;
        --dedup_strategy) DEDUP_STRATEGY="$2"; shift 2 ;;
        --col_variant) COL_VARIANT="$2"; shift 2 ;;
        --col_chr) COL_CHR="$2"; shift 2 ;;
        --col_pos) COL_POS="$2"; shift 2 ;;
        --col_gene) COL_GENE="$2"; shift 2 ;;
        --col_beta) COL_BETA="$2"; shift 2 ;;
        --col_pvalue) COL_PVALUE="$2"; shift 2 ;;
        -v|--verbose) VERBOSE="--verbose"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Error: Unknown option $1"; usage; exit 1 ;;
    esac
done

# Validate
if [[ -z "$INPUT" && -z "$MANIFEST" ]]; then
    echo "Error: --input or --manifest required"
    usage
    exit 1
fi

if [[ ! -f "$QTL_SCRIPT" ]]; then
    echo "Error: Core R script not found: $QTL_SCRIPT"
    exit 1
fi

# Build command
mkdir -p "$OUTPUT_DIR"

CMD="Rscript ${QTL_SCRIPT}"

if [[ -n "$INPUT" ]]; then
    CMD+=" --input ${INPUT}"
fi

if [[ -n "$MANIFEST" ]]; then
    CMD+=" --manifest ${MANIFEST}"
fi

CMD+=" --format ${FORMAT}"
CMD+=" --input_build ${INPUT_BUILD}"
CMD+=" --output_build ${OUTPUT_BUILD}"
CMD+=" --dedup_strategy ${DEDUP_STRATEGY}"
CMD+=" --pvalue_threshold ${PVALUE_THRESHOLD}"
CMD+=" --output ${OUTPUT_DIR}/${OUTPUT_PREFIX}"

[[ -n "$CHAIN_FILE" ]] && CMD+=" --chain_file ${CHAIN_FILE}"
[[ -n "$TISSUES" ]] && CMD+=" --tissues ${TISSUES}"
[[ -n "$ANCESTRIES" ]] && CMD+=" --ancestries ${ANCESTRIES}"
[[ -n "$QTL_TYPES" ]] && CMD+=" --qtl_types ${QTL_TYPES}"
[[ -n "$GENES" ]] && CMD+=" --genes ${GENES}"
[[ -n "$REGION" ]] && CMD+=" --region ${REGION}"
[[ -n "$SPLIT_BY" ]] && CMD+=" --split_by ${SPLIT_BY}"

# Custom columns
[[ -n "$COL_VARIANT" ]] && CMD+=" --col_variant ${COL_VARIANT}"
[[ -n "$COL_CHR" ]] && CMD+=" --col_chr ${COL_CHR}"
[[ -n "$COL_POS" ]] && CMD+=" --col_pos ${COL_POS}"
[[ -n "$COL_GENE" ]] && CMD+=" --col_gene ${COL_GENE}"
[[ -n "$COL_BETA" ]] && CMD+=" --col_beta ${COL_BETA}"
[[ -n "$COL_PVALUE" ]] && CMD+=" --col_pvalue ${COL_PVALUE}"

[[ -n "$VERBOSE" ]] && CMD+=" ${VERBOSE}"

# Execute
echo "=============================================="
echo "QTL Harmonization Pipeline"
echo "=============================================="
echo "Output: ${OUTPUT_DIR}/${OUTPUT_PREFIX}"
echo "=============================================="
echo ""

LOG_FILE="${OUTPUT_DIR}/${OUTPUT_PREFIX}.log"
$CMD 2>&1 | tee "$LOG_FILE"

echo ""
echo "Output files:"
ls -la "${OUTPUT_DIR}/${OUTPUT_PREFIX}"* 2>/dev/null || true
