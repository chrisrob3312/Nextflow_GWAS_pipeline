#!/bin/bash
# Heritability Estimation Shell Wrapper
# Maintains parity with modules/local/heritability/main.nf

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${SCRIPT_DIR}/../bin"
H2_SCRIPT="${BIN_DIR}/estimate_heritability.R"

usage() {
    cat <<EOF
Heritability Estimation Pipeline

Estimates SNP heritability using methods appropriate for admixed populations.
For admixed cohorts, use cov_ldsc or greml with ancestry PCs.

Usage: $(basename "$0") [OPTIONS]

METHOD:
    --method METHOD           Method: ldsc, cov_ldsc, greml, bolt_reml
                              [default: cov_ldsc]

SUMMARY STAT INPUT (for LDSC methods):
    --sumstats FILE           GWAS summary statistics file
    --sumstats_format FMT     Format: auto, ldsc, gwas_ssf, regenie, saige
    --ld_scores PREFIX        LD score files prefix

INDIVIDUAL-LEVEL INPUT (for GREML/BOLT):
    --grm PREFIX              GRM prefix for GREML
    --plink PREFIX            PLINK prefix for computing GRM
    --phenotype FILE          Phenotype file
    --trait NAME              Trait column name
    --covariates LIST         Comma-separated covariate columns

ANCESTRY STRATIFICATION:
    --ancestry_file FILE      Sample ancestry assignments
    --ancestry_col COL        Ancestry column name [default: ancestry]
    --ancestries LIST         Comma-separated ancestries to analyze

GENETIC CORRELATION:
    --estimate_rg             Estimate cross-ancestry/trait rg
    --sumstats2 FILE          Second GWAS for rg estimation
    --rg_method METHOD        rg method: popcorn, s_ldxr, ldsc [default: popcorn]

OUTPUT:
    -o, --output_prefix STR   Output prefix [default: heritability]
    -d, --output_dir DIR      Output directory [default: .]

RUNTIME:
    --threads INT             Number of threads [default: 4]
    -v, --verbose             Verbose output
    -h, --help                Show this help

METHODS:
    ldsc        Standard LDSC. WARNING: Biased for admixed populations.
    cov_ldsc    Covariate-stratified LDSC. Better for admixed, stratifies by ancestry.
    greml       GCTA-GREML. Gold standard, requires individual-level data.
    bolt_reml   BOLT-LMM REML. Fast, handles relatedness.

EXAMPLES:
    # Standard LDSC (EUR reference)
    $(basename "$0") --method ldsc --sumstats gwas.sumstats.gz \\
        --ld_scores eur_w_ld_chr/

    # cov-LDSC for admixed cohort
    $(basename "$0") --method cov_ldsc --sumstats gwas.sumstats.gz \\
        --ld_scores multi_ancestry_ld/ --ancestries "EUR,AFR,LATINO"

    # GREML with ancestry PCs as covariates
    $(basename "$0") --method greml --grm cohort --phenotype pheno.tsv \\
        --trait relapse --covariates "age,sex,PC1,PC2,PC3,PC4,PC5"

    # Cross-ancestry genetic correlation
    $(basename "$0") --method ldsc --sumstats eur.sumstats.gz \\
        --sumstats2 afr.sumstats.gz --estimate_rg --rg_method popcorn

EOF
}

# Defaults
METHOD="cov_ldsc"
SUMSTATS=""
SUMSTATS_FORMAT="auto"
LD_SCORES=""
GRM=""
PLINK=""
PHENOTYPE=""
TRAIT=""
COVARIATES=""
ANCESTRY_FILE=""
ANCESTRY_COL="ancestry"
ANCESTRIES=""
ESTIMATE_RG=""
SUMSTATS2=""
RG_METHOD="popcorn"
OUTPUT_PREFIX="heritability"
OUTPUT_DIR="."
THREADS=4
VERBOSE=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --method) METHOD="$2"; shift 2 ;;
        --sumstats) SUMSTATS="$2"; shift 2 ;;
        --sumstats_format) SUMSTATS_FORMAT="$2"; shift 2 ;;
        --ld_scores) LD_SCORES="$2"; shift 2 ;;
        --grm) GRM="$2"; shift 2 ;;
        --plink) PLINK="$2"; shift 2 ;;
        --phenotype) PHENOTYPE="$2"; shift 2 ;;
        --trait) TRAIT="$2"; shift 2 ;;
        --covariates) COVARIATES="$2"; shift 2 ;;
        --ancestry_file) ANCESTRY_FILE="$2"; shift 2 ;;
        --ancestry_col) ANCESTRY_COL="$2"; shift 2 ;;
        --ancestries) ANCESTRIES="$2"; shift 2 ;;
        --estimate_rg) ESTIMATE_RG="--estimate_rg"; shift ;;
        --sumstats2) SUMSTATS2="$2"; shift 2 ;;
        --rg_method) RG_METHOD="$2"; shift 2 ;;
        -o|--output_prefix) OUTPUT_PREFIX="$2"; shift 2 ;;
        -d|--output_dir) OUTPUT_DIR="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        -v|--verbose) VERBOSE="--verbose"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Error: Unknown option $1"; usage; exit 1 ;;
    esac
done

# Validate
if [[ ! -f "$H2_SCRIPT" ]]; then
    echo "Error: Core R script not found: $H2_SCRIPT"
    exit 1
fi

# Build command
mkdir -p "$OUTPUT_DIR"

CMD="Rscript ${H2_SCRIPT}"
CMD+=" --method ${METHOD}"
CMD+=" --output_prefix ${OUTPUT_DIR}/${OUTPUT_PREFIX}"
CMD+=" --threads ${THREADS}"

[[ -n "$SUMSTATS" ]] && CMD+=" --sumstats ${SUMSTATS}"
[[ -n "$SUMSTATS_FORMAT" ]] && CMD+=" --sumstats_format ${SUMSTATS_FORMAT}"
[[ -n "$LD_SCORES" ]] && CMD+=" --ld_scores ${LD_SCORES}"
[[ -n "$GRM" ]] && CMD+=" --grm ${GRM}"
[[ -n "$PLINK" ]] && CMD+=" --plink ${PLINK}"
[[ -n "$PHENOTYPE" ]] && CMD+=" --phenotype ${PHENOTYPE}"
[[ -n "$TRAIT" ]] && CMD+=" --trait ${TRAIT}"
[[ -n "$COVARIATES" ]] && CMD+=" --covariates ${COVARIATES}"
[[ -n "$ANCESTRY_FILE" ]] && CMD+=" --ancestry_file ${ANCESTRY_FILE}"
[[ -n "$ANCESTRY_COL" ]] && CMD+=" --ancestry_col ${ANCESTRY_COL}"
[[ -n "$ANCESTRIES" ]] && CMD+=" --ancestries ${ANCESTRIES}"
[[ -n "$ESTIMATE_RG" ]] && CMD+=" ${ESTIMATE_RG}"
[[ -n "$SUMSTATS2" ]] && CMD+=" --sumstats2 ${SUMSTATS2}"
[[ -n "$RG_METHOD" ]] && CMD+=" --rg_method ${RG_METHOD}"
[[ -n "$VERBOSE" ]] && CMD+=" ${VERBOSE}"

# Execute
echo "=============================================="
echo "Heritability Estimation Pipeline"
echo "=============================================="
echo "Method: $METHOD"
echo "Output: ${OUTPUT_DIR}/${OUTPUT_PREFIX}"
echo "=============================================="
echo ""

LOG_FILE="${OUTPUT_DIR}/${OUTPUT_PREFIX}.log"
$CMD 2>&1 | tee "$LOG_FILE"

echo ""
echo "Output files:"
ls -la "${OUTPUT_DIR}/${OUTPUT_PREFIX}"* 2>/dev/null || true
