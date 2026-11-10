#!/usr/bin/env Rscript

# ============================================================================
# TractorGENESIS: Local Ancestry-Aware Mixed Model Association Testing
# ============================================================================
# Production-ready adapter integrating Tractor's local ancestry decomposition
# with GENESIS mixed model framework.
#
# SUPPORTS:
#   - 2-way admixture (e.g., AAC: AFR-EUR)
#   - 3-way admixture (e.g., Latino: EUR-AFR-AMR)
#   - N-way admixture (any number of ancestral populations)
#   - Any ancestry backgrounds (user-specified)
#
# STATISTICAL MODEL (Atkinson et al., 2021):
#   For k ancestral populations:
#   g(E[Y]) = b0 + sum_{j=1}^{k-1} b_LA_j * LA_j
#                + sum_{j=1}^{k} b_dose_j * Dose_j + covariates
#
#   Where:
#     LA_j    = Local ancestry count for ancestry j (0, 1, or 2 haplotypes)
#     Dose_j  = Allele dosage on ancestry j haplotypes (0 to LA_j)
#
# OUTPUT INTERPRETATION:
#   - P_JOINT: Tests if SNP has ANY effect (across all ancestries)
#   - BETA_{anc}: Effect size of allele on {anc} ancestry background
#   - P_{anc}: P-value for ancestry-specific effect
#   - P_HET: Tests if effects DIFFER across ancestries (heterogeneity)
#   - I2: Heterogeneity measure (0-100%, higher = more different)
#
# INPUTS:
#   From Tractor ExtractTracts:
#     - {prefix}.hapcount.{anc}.txt.gz: Haplotype counts per ancestry
#     - {prefix}.dosage.{anc}.txt.gz: Ancestry-deconvoluted allele dosages
#   (--msp / --gds, RFMix input, is not implemented: the run stops if it is given)
#
# ============================================================================

suppressPackageStartupMessages({
    library(survival)
    library(data.table)
    library(optparse)
    library(Matrix)
})
# GENESIS (fitNullModel, with GWASTools' ScanAnnotationDataFrame) is loaded below only for the
# binary / quantitative null model; the survival path needs survival (+ coxme with a GRM).

# ============================================================================
# Command-line Arguments
# ============================================================================
option_list <- list(
    # Input files
    make_option(c("--tractor_prefix"), type = "character", default = NULL,
                help = "Prefix for Tractor output files (expects .hapcount.{anc}.txt.gz and .dosage.{anc}.txt.gz)"),
    make_option(c("--msp"), type = "character", default = NULL,
                help = "RFMix MSP file (alternative to Tractor files)"),
    make_option(c("--gds"), type = "character", default = NULL,
                help = "GDS file with genotypes (required with --msp)"),
    make_option(c("-p", "--phenotype"), type = "character", default = NULL,
                help = "Phenotype file (TSV with header, must include sample_id column)"),
    make_option(c("-k", "--kinship"), type = "character", default = NULL,
                help = "Kinship matrix (RDS file, for relatedness adjustment)"),

    # Trait specification
    make_option(c("-t", "--trait"), type = "character", default = NULL,
                help = "Trait column name in phenotype file"),
    make_option(c("-m", "--model"), type = "character", default = "binary",
                help = "Model type: 'binary', 'quantitative', 'survival' [default: binary]"),
    make_option(c("--time_col"), type = "character", default = "time",
                help = "Time column for survival analysis [default: time]"),
    make_option(c("--event_col"), type = "character", default = "event",
                help = "Event column for survival (1=event, 0=censored) [default: event]"),

    # Ancestry specification (FLEXIBLE - any populations)
    make_option(c("-a", "--ancestries"), type = "character", default = NULL,
                help = "Comma-separated ancestral populations in order matching Tractor/RFMix coding (e.g., 'EUR,AFR' or 'EUR,AFR,AMR')"),
    make_option(c("--ref_ancestry"), type = "character", default = NULL,
                help = "Reference ancestry for LA terms (excluded to avoid collinearity). Default: first ancestry"),

    # Covariates
    make_option(c("-c", "--covariates"), type = "character", default = NULL,
                help = "Comma-separated covariate column names; each must be a numeric phenotype column (checked before anything is read or fitted). The list used is written as '## covariates=' above every results table"),
    make_option(c("--global_ancestry_prefix"), type = "character", default = NULL,
                help = "Prefix for global ancestry columns (e.g., 'prop_' expects prop_EUR, prop_AFR, etc.)"),

    # Output
    make_option(c("-o", "--output_prefix"), type = "character", default = "tractor_genesis",
                help = "Output file prefix [default: tractor_genesis]"),

    # Variant filters
    make_option(c("--mac_min"), type = "integer", default = 10,
                help = paste("Minimum MINOR allele count: a variant is tested only if its stratum-wide minor allele count",
                             "(min(ALT, 2n - ALT)) reaches it, and an ancestry's dose term only if that background's does",
                             "(min(ALT_a, haplotypes_a - ALT_a)) [default: 10]")),

    # Runtime
    make_option(c("--chunk_size"), type = "integer", default = 1000,
                help = "Variants per chunk [default: 1000]"),
    make_option(c("--threads"), type = "integer", default = 1,
                help = "Number of threads [default: 1]"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output"),

    # SLURM array job support
    make_option(c("--slurm_array_task_id"), type = "integer", default = NULL,
                help = "SLURM_ARRAY_TASK_ID (auto-detected if not set)"),
    make_option(c("--slurm_array_task_count"), type = "integer", default = NULL,
                help = "Total array tasks (for chunking)"),
    make_option(c("--chromosome"), type = "character", default = NULL,
                help = "Chromosome to analyze (for array parallelization)"),

    # Ancestry stratification
    make_option(c("--stratum"), type = "character", default = NULL,
                help = paste("Stratum these Tractor files belong to (EUR, AAC, LAT1, LAT2, OTHER, POOLED, ...).",
                             "Every sample in them must be in this stratum per --strata_file (or a 'stratum'",
                             "column of the phenotype file); the outputs carry it in a STRATUM column")),
    make_option(c("--strata_file"), type = "character", default = NULL,
                help = "sample_id<TAB>stratum table to check --stratum against (the handoff's ancestry/strata.tsv)"),
    make_option(c("--sample_include"), type = "character", default = NULL,
                help = paste("File of sample IDs (one per line; results/strata/<S>.samples.txt) the analysis is",
                             "restricted to, intersected with the phenotype and Tractor samples")),
    make_option(c("--min_stratum_n"), type = "integer", default = 50,
                help = paste("Fewer analysable samples than this: the stratum is skipped (exit 0, no results,",
                             "status skipped_below_min_stratum_n in <output_prefix>.run_counts.tsv) [default: 50]")),
    make_option(c("--ancestry_index_map"), type = "character", default = NULL,
                help = paste("anc_index<TAB>label table (the handoff's local_ancestry/flare_index_map.tsv), needed only",
                             "when the Tractor files are numbered (anc0, 0, ...) instead of named by ancestry")),
    make_option(c("--ancestry_config"), type = "character", default = NULL,
                help = "Path to ancestry_config.R for cohort-specific settings")
)

opt <- parse_args(OptionParser(
    option_list = option_list,
    prog = "tractor_genesis_adapter.R",
    description = "TractorGENESIS: Local ancestry-aware mixed model GWAS"
))

# ============================================================================
# Validation
# ============================================================================
if (is.null(opt$phenotype) || is.null(opt$trait) || is.null(opt$ancestries)) {
    stop("Required: --phenotype, --trait, --ancestries")
}

if (is.null(opt$tractor_prefix) && is.null(opt$msp)) {
    stop("Required: --tractor_prefix OR (--msp and --gds)")
}

if (!is.null(opt$msp) || !is.null(opt$gds)) {
    stop("--msp / --gds (RFMix input) is not implemented: pass the Tractor files with --tractor_prefix")
}
if (opt$model != "survival") suppressPackageStartupMessages({ library(GENESIS); library(GWASTools) })

# ============================================================================
# SLURM Array Job Detection
# ============================================================================
slurm_task_id <- opt$slurm_array_task_id
if (is.null(slurm_task_id)) {
    slurm_task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", ""))
    if (is.na(slurm_task_id)) slurm_task_id <- NULL
}

if (!is.null(slurm_task_id)) {
    cat("Running as SLURM array task:", slurm_task_id, "\n")

    # If chromosome specified via array, use task ID
    if (is.null(opt$chromosome) && slurm_task_id >= 1 && slurm_task_id <= 22) {
        opt$chromosome <- as.character(slurm_task_id)
        cat("  Chromosome (from array):", opt$chromosome, "\n")
    }
}

# Set thread count from SLURM if available
slurm_cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", ""))
if (!is.na(slurm_cpus) && slurm_cpus > 0) {
    opt$threads <- slurm_cpus
    cat("Using SLURM CPUs:", opt$threads, "\n")
}

# Load ancestry config if provided
if (!is.null(opt$ancestry_config) && file.exists(opt$ancestry_config)) {
    source(opt$ancestry_config)
    cat("Loaded ancestry config from:", opt$ancestry_config, "\n")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Shared helpers (file resolution, stratum check, N gate, counts sidecar), next to this script
.adapter_dir <- local({
    f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
    if (length(f) > 0) dirname(normalizePath(f[1])) else "."
})
source(file.path(.adapter_dir, "tractor_io.R"))
source(file.path(.adapter_dir, "qc_assert.R"))    # qc_assert / qc_account (CP-11), check_covariates (CP-08)
source(file.path(.adapter_dir, "sumstats_io.R"))
source(file.path(.adapter_dir, "tractor_model.R")) # the score test (genesis_ops, cox_ops, test_variant_tractor)
dir.create(dirname(opt$output_prefix), recursive = TRUE, showWarnings = FALSE)  # fwrite_sumstats: "## covariates=" above every table (CP-08)

# Parse ancestries
ancestries <- strsplit(opt$ancestries, ",")[[1]]
n_anc <- length(ancestries)
ref_anc <- if (!is.null(opt$ref_ancestry)) opt$ref_ancestry else ancestries[1]

if (!ref_anc %in% ancestries) {
    stop("Reference ancestry '", ref_anc, "' not in ancestries list")
}

# Non-reference ancestries (for LA terms)
non_ref_anc <- setdiff(ancestries, ref_anc)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║           TractorGENESIS: Local Ancestry-Aware GWAS             ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║ Model:", sprintf("%-58s", paste0(n_anc, "-way admixture (", paste(ancestries, collapse = "-"), ")")), "║\n")
cat("║ Trait:", sprintf("%-58s", opt$trait), "║\n")
cat("║ Type: ", sprintf("%-58s", opt$model), "║\n")
cat("║ Reference ancestry:", sprintf("%-45s", ref_anc), "║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# ============================================================================
# Statistical Model Display
# ============================================================================
cat("STATISTICAL MODEL:\n")
cat("  g(E[Y]) = β₀")

# Display with numbered coefficients to clarify these are SEPARATE terms
coef_num <- 1
for (anc in non_ref_anc) {
    cat(" + β", coef_num, "·LA_", anc, sep = "")
    coef_num <- coef_num + 1
}
for (anc in ancestries) {
    cat(" + β", coef_num, "·Dose_", anc, sep = "")
    coef_num <- coef_num + 1
}
cat(" + covariates\n\n")

# Show coefficient mapping
cat("COEFFICIENT MAPPING:\n")
coef_num <- 1
for (anc in non_ref_anc) {
    cat("  β", coef_num, " = effect of having ", anc, " ancestry at this locus (LA term)\n", sep = "")
    coef_num <- coef_num + 1
}
for (anc in ancestries) {
    cat("  β", coef_num, " = effect of allele WHEN carried on ", anc, " haplotype (Dose term)\n", sep = "")
    coef_num <- coef_num + 1
}
cat("\n")

cat("NOTE: Each β is an INDEPENDENT coefficient (not interaction terms).\n")
cat("      LA terms: k-1 coefficients (", paste(non_ref_anc, collapse = ", "), " vs ", ref_anc, " reference)\n", sep = "")
cat("      Dose terms: k coefficients (all ancestries)\n\n")

# ============================================================================
# Load Phenotype Data
# ============================================================================
cat("Loading phenotype data...\n")
pheno <- fread(opt$phenotype)

# Ensure sample_id column exists
if (!"sample_id" %in% names(pheno)) {
    # Try common alternatives
    id_cols <- c("IID", "SampleID", "sample", "ID")
    found <- intersect(id_cols, names(pheno))
    if (length(found) > 0) {
        setnames(pheno, found[1], "sample_id")
        cat("  Using '", found[1], "' as sample_id\n", sep = "")
    } else {
        stop("No sample_id column found. Expected: sample_id, IID, SampleID, sample, or ID")
    }
}

n_samples <- nrow(pheno)
cat("  Loaded", n_samples, "samples\n")

# Every covariate must be a numeric column of the phenotype file (CP-08): checked here, before
# anything else is read or fitted, so a misspelt or text column stops the task with its name
# instead of failing inside a model (or, as before, being coerced to NA without a word)
covariates_requested <- if (!is.null(opt$covariates)) trimws(strsplit(opt$covariates, ",")[[1]]) else character()
check_covariates(pheno, covariates_requested, opt$phenotype)

# Trait info
if (opt$model == "binary") {
    n_cases <- sum(pheno[[opt$trait]] == 1, na.rm = TRUE)
    n_controls <- sum(pheno[[opt$trait]] == 0, na.rm = TRUE)
    cat("  Cases:", n_cases, "| Controls:", n_controls, "\n")
} else if (opt$model == "survival") {
    n_events <- sum(pheno[[opt$event_col]] == 1, na.rm = TRUE)
    cat("  Events:", n_events, "| Censored:", n_samples - n_events, "\n")
} else {
    cat("  Mean:", round(mean(pheno[[opt$trait]], na.rm = TRUE), 3),
        "| SD:", round(sd(pheno[[opt$trait]], na.rm = TRUE), 3), "\n")
}

# ============================================================================
# Load Kinship Matrix
# ============================================================================
kinship <- NULL
if (!is.null(opt$kinship) && file.exists(opt$kinship)) {
    cat("\nLoading kinship matrix...\n")
    kinship <- readRDS(opt$kinship)
    cat("  Dimensions:", nrow(kinship), "x", ncol(kinship), "\n")
}

# ============================================================================
# Parse Tractor Input Files
# ============================================================================
cat("\nLoading Tractor data...\n")

# Which file holds which ancestry: named by label, or numbered and mapped through
# --ancestry_index_map; never by position in --ancestries (bin/tractor_io.R)
index_map <- if (!is.null(opt$ancestry_index_map)) read_index_map(opt$ancestry_index_map) else NULL
tractor_files <- lapply(setNames(ancestries, ancestries), function(a) list(
    hapcount = resolve_tractor_file(opt$tractor_prefix, a, "hapcount", index_map),
    dosage   = resolve_tractor_file(opt$tractor_prefix, a, "dosage",   index_map)
))
for (a in ancestries) {
    cat("  ", a, ": ", basename(tractor_files[[a]]$hapcount$path), ", ", basename(tractor_files[[a]]$dosage$path),
        " (", tractor_files[[a]]$dosage$how, ")\n", sep = "")
}

# --stratum and --min_stratum_n, from the file headers, before the matrices are read.
# Every join below writes one sample_accounting row and the checks file verdicts (CP-11); the
# combine step's gate requires adapter_stratum_check, adapter_sample_join, adapter_grm_cover and
# adapter_model_n from every task, so a skipped or failed join cannot pass unseen.
tractor_ids <- tractor_sample_ids(tractor_files[[1]]$hapcount$path)
qc_where <- sprintf("%s / %s / chr%s", opt$trait, opt$stratum %||% "unnamed", opt$chromosome %||% "all")
acct <- function(step, n_in, n_out, reason) {
    qc_account(step, n_in, n_out, reason, stratum = opt$stratum %||% NA, trait = opt$trait, chromosome = opt$chromosome %||% NA)
}
acct("phenotype_input", nrow(pheno), nrow(pheno), paste("rows in", basename(opt$phenotype)))
stratum_check <- tryCatch(check_stratum_membership(opt$stratum, tractor_ids, opt$strata_file, pheno), error = function(e) {
    qc_assert("adapter_stratum_check", FALSE, n = length(tractor_ids), message = paste(qc_where, conditionMessage(e)))
    stop(e)
})
qc_assert("adapter_stratum_check", verdict = if (stratum_check$status == "passed") "PASS" else "INFO", n = length(tractor_ids),
          message = paste0(qc_where, ": ", stratum_check$status, " (", stratum_check$source, ")"))
acct("stratum", length(tractor_ids), length(tractor_ids), paste("Tractor samples; stratum check", stratum_check$status))
cat("  Stratum check:", stratum_check$status, "(", stratum_check$source, ")\n")
joined <- intersect(pheno$sample_id, tractor_ids)
acct("tractor_join", nrow(pheno), length(joined),
     sprintf("phenotype rows also in the Tractor files (%d Tractor samples have no phenotype row)", sum(!tractor_ids %in% pheno$sample_id)))
qc_assert("adapter_sample_join", length(joined) > 0, n = length(joined),
          message = sprintf("%s: %d of %d phenotype rows and %d Tractor samples joined", qc_where, length(joined), nrow(pheno), length(tractor_ids)))
if (!is.null(opt$sample_include)) {
    if (!file.exists(opt$sample_include)) stop("--sample_include ", opt$sample_include, " does not exist")
    keep_ids <- unique(trimws(readLines(opt$sample_include)))
    n_j <- length(joined)
    joined <- intersect(joined, keep_ids)
    acct("sample_include", n_j, length(joined), paste("joined samples listed in", basename(opt$sample_include)))
}
if (!is.null(kinship)) {
    # The GRM is keyed by the .fam IID (docs/HANDOFF_CONTRACT.md rule 1): a joined sample missing
    # from it is a broken input, not something to drop quietly.
    off_grm <- setdiff(joined, rownames(kinship))
    acct("grm_intersection", length(joined), length(joined) - length(off_grm), paste("joined samples in", basename(opt$kinship)))
    qc_assert("adapter_grm_cover", length(off_grm) == 0, n = length(joined),
              message = if (length(off_grm) == 0) paste0(qc_where, ": every joined sample is in the GRM") else
                  sprintf("%s: %d joined samples are not in %s, e.g. %s", qc_where, length(off_grm), basename(opt$kinship),
                          paste(head(off_grm, 3), collapse = ", ")))
    if (length(off_grm) > 0) stop("sample ", off_grm[1], " missing from GRM ", opt$kinship, " (", length(off_grm),
                                  " joined sample(s) are not in it: ", paste(head(off_grm, 5), collapse = ", "), ")")
} else {
    qc_assert("adapter_grm_cover", verdict = "INFO", n = length(joined), message = paste0(qc_where, ": no --kinship; fitted without relatedness"))
}
# Covariates: the ones asked for, plus the global-ancestry columns when --global_ancestry_prefix
covariates <- if (length(covariates_requested) > 0) covariates_requested else NULL
if (!is.null(opt$global_ancestry_prefix)) {
    ga_cols <- intersect(paste0(opt$global_ancestry_prefix, non_ref_anc), names(pheno))
    if (length(ga_cols) > 0) {
        covariates <- c(covariates, ga_cols)
        cat("  Added global ancestry covariates:", paste(ga_cols, collapse = ", "), "\n")
    }
}
# Complete cases (CP-17): the null model and the score test use only the joined samples with the
# outcome and every covariate; the others are dropped here and counted, never left for a model
# function to drop on its own (which misaligned residuals and genotypes before).
outcome_cols <- if (opt$model == "survival") c(opt$time_col, opt$event_col) else opt$trait
missing_cols <- setdiff(outcome_cols, names(pheno))
if (length(missing_cols) > 0) stop("outcome column(s) not in ", opt$phenotype, ": ", paste(missing_cols, collapse = ", "))
cc_cols <- c(outcome_cols, covariates)
in_joined <- pheno$sample_id %in% joined
analysable <- pheno$sample_id[in_joined & complete.cases(pheno[, cc_cols, with = FALSE])]
n_null <- length(analysable)
acct("complete_cases", length(joined), n_null, paste("joined samples with", paste(cc_cols, collapse = ", ")))
if (n_null < length(joined)) cat("  Complete cases:", n_null, "of", length(joined), "joined samples;", length(joined) - n_null,
                                 "dropped for a missing", paste(cc_cols, collapse = " / "), "value\n")
acct("min_stratum_n", n_null, if (n_null < opt$min_stratum_n) 0 else n_null, paste("--min_stratum_n", opt$min_stratum_n))
run_counts <- list(status = "ok", stratum = opt$stratum %||% NA, stratum_check = stratum_check$status,
                   min_stratum_n = opt$min_stratum_n, n_null = n_null, n_variants = NA,
                   n_pass_mac = NA, n_tested = NA, n_failed = NA)
if (n_null < opt$min_stratum_n) {
    run_counts$status <- "skipped_below_min_stratum_n"
    write_run_counts(opt$output_prefix, run_counts)
    qc_assert("adapter_model_n", verdict = "INFO", n = n_null,
              message = sprintf("%s: skipped, %d analysable samples < --min_stratum_n %d", qc_where, n_null, opt$min_stratum_n))
    cat("\nSKIPPED: stratum ", opt$stratum %||% "(unnamed)", " has ", n_null,
        " analysable samples, fewer than --min_stratum_n ", opt$min_stratum_n,
        ". No association results are written; see ", opt$output_prefix, ".run_counts.tsv\n", sep = "")
    quit(save = "no", status = 0)
}

# Function to read Tractor hapcount/dosage files (paths resolved above)
read_tractor_files <- function(files, ancestries) {
    hapcount <- list()
    dosage <- list()
    variant_info <- NULL

    for (anc in ancestries) {
        hc_file  <- files[[anc]]$hapcount$path
        dos_file <- files[[anc]]$dosage$path

        cat("  Loading", anc, "...\n")
        cat("    Hapcount:", basename(hc_file), "\n")
        cat("    Dosage:", basename(dos_file), "\n")

        # Read files
        # Tractor format: first cols are variant info, rest are samples
        hc_data <- fread(hc_file)
        dos_data <- fread(dos_file)

        # Extract variant info from first file
        if (is.null(variant_info)) {
            # Variant info columns (CHR, POS, ID, REF, ALT or similar; else the first 5)
            info_cols <- tractor_info_cols(names(hc_data))
            variant_info <- hc_data[, ..info_cols]

            # Standardize column names
            names(variant_info) <- toupper(names(variant_info))
            if ("#CHROM" %in% names(variant_info)) setnames(variant_info, "#CHROM", "CHR")
            if ("CHROM" %in% names(variant_info)) setnames(variant_info, "CHROM", "CHR")
            if ("BP" %in% names(variant_info)) setnames(variant_info, "BP", "POS")
            if ("SNP" %in% names(variant_info) && !"ID" %in% names(variant_info)) {
                setnames(variant_info, "SNP", "ID")
            }
        }

        # Extract sample data (non-info columns)
        sample_cols <- setdiff(names(hc_data), info_cols)

        hapcount[[anc]] <- as.matrix(hc_data[, ..sample_cols])
        dosage[[anc]] <- as.matrix(dos_data[, ..sample_cols])
    }

    # Get sample IDs from column names
    sample_ids <- colnames(hapcount[[ancestries[1]]])

    cat("  Variants:", nrow(variant_info), "\n")
    cat("  Samples:", length(sample_ids), "\n")

    list(
        hapcount = hapcount,
        dosage = dosage,
        variant_info = variant_info,
        sample_ids = sample_ids
    )
}

# Load Tractor data
tractor_data <- read_tractor_files(tractor_files, ancestries)

# The analysed samples: complete cases (above) that are in the Tractor files, in phenotype order
common_samples <- intersect(analysable, tractor_data$sample_ids)
cat("  Samples analysed (complete cases with phenotype and Tractor data):", length(common_samples), "\n")

if (length(common_samples) == 0) {
    stop("No overlapping samples between phenotype and Tractor data!")
}

# Subset and align data
pheno <- pheno[sample_id %in% common_samples]
pheno <- pheno[match(common_samples, sample_id)]  # Align order

sample_idx <- match(common_samples, tractor_data$sample_ids)
for (anc in ancestries) {
    tractor_data$hapcount[[anc]] <- tractor_data$hapcount[[anc]][, sample_idx, drop = FALSE]
    tractor_data$dosage[[anc]] <- tractor_data$dosage[[anc]][, sample_idx, drop = FALSE]
}

n_samples <- length(common_samples)
n_variants <- nrow(tractor_data$variant_info)

# ============================================================================
# Fit Null Model (no genetic terms)
# ============================================================================
cat("\nFitting null model...\n")


if (!is.null(covariates)) {
    cat("  Covariates:", paste(covariates, collapse = ", "), "\n")
}

# Subset kinship to common samples (every one of them is in it: adapter_grm_cover above)
if (!is.null(kinship)) {
    kinship <- kinship[common_samples, common_samples]
}

# ----------------------------------------------------------------------------
# ONE score-test engine for all trait types, so MRD (binary), relapse and OS
# (time-to-event) are directly comparable:
#     U = D' P y        V = D' P D        D = [LA terms | Dose terms]
#
#   binary / quantitative : P and Py from the GENESIS mixed-model null
#                           (kinship enters through Sigma)
#   survival              : y -> martingale residuals from a Cox null
#                           (coxme frailty on 2*kinship when available);
#                           P from the Breslow risk-set information,
#                           projected on the null covariates.
#                           This is the Cox score test used by SPACox/GATE.
# ----------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a)) b else a

# Checked again with the global-ancestry columns added (CP-08); the realised list is written as
# "## covariates=" above every results table
if (length(covariates) > 0) check_covariates(pheno, covariates, opt$phenotype)
sumstats_header <- c(covariates = if (length(covariates) > 0) paste(covariates, collapse = ",") else "none")
cov_mat <- if (length(covariates) > 0) {
    m <- as.matrix(pheno[, covariates, with = FALSE]); storage.mode(m) <- "double"; m
} else NULL

# --- fit the null model ------------------------------------------------------
if (opt$model == "survival") {
    if (!all(c(opt$time_col, opt$event_col) %in% names(pheno))) {
        stop("Survival columns not found: ", opt$time_col, ", ", opt$event_col)
    }
    surv_time  <- as.numeric(pheno[[opt$time_col]])
    surv_event <- as.integer(pheno[[opt$event_col]])
    if (anyNA(surv_time) || anyNA(surv_event)) stop("missing time / event values after the complete-case filter (a bug)")

    surv_df <- data.frame(time = surv_time, event = surv_event)
    if (!is.null(cov_mat)) surv_df <- cbind(surv_df, as.data.frame(cov_mat))
    cov_str <- if (length(covariates) > 0) paste("+", paste(covariates, collapse = " + ")) else ""
    has_coxme <- requireNamespace("coxme", quietly = TRUE)

    if (!is.null(kinship) && has_coxme) {
        K <- as.matrix(kinship) * 2                     # expected relationship = 2 * kinship
        dimnames(K) <- list(common_samples, common_samples)
        surv_df$sample_id <- common_samples
        f_me <- as.formula(paste("Surv(time, event) ~ 1", cov_str, "+ (1 | sample_id)"))
        cox_null <- coxme::coxme(f_me, data = surv_df,
                                 varlist = coxme::coxmeMlist(list(K), rescale = FALSE))
        eta <- as.numeric(cox_null$linear.predictor)   # fixed effects + frailty
        cat("  Fitted Cox MIXED model (coxme; frailty on 2 x kinship)\n")
        cat("    Frailty variance:", signif(as.numeric(coxme::VarCorr(cox_null)[[1]]), 4), "\n")
    } else {
        if (!is.null(kinship) && !has_coxme) {
            cat("  WARNING: coxme not installed - Cox null fitted WITHOUT kinship. Install coxme or pass an unrelated set.\n")
        }
        f_null <- as.formula(paste("Surv(time, event) ~ 1", cov_str))
        cox_null <- coxph(f_null, data = surv_df, ties = "breslow")
        eta <- as.numeric(cox_null$linear.predictors)
        cat("  Fitted Cox proportional hazards model (Breslow ties)\n")
    }
    null_ops <- cox_ops(surv_time, surv_event, eta, cov_mat)
    nullmod  <- list(type = "cox", model = cox_null, martingale = null_ops$PY,
                     n_events = sum(surv_event == 1))
    cat("  Events:", nullmod$n_events, "| Martingale residual range:",
        paste(signif(range(null_ops$PY), 3), collapse = " to "), "\n")

} else {
    scanAnnot <- ScanAnnotationDataFrame(data.frame(
        scanID = common_samples,
        pheno[, c(opt$trait, covariates), with = FALSE]
    ))
    fam <- if (opt$model == "binary") binomial(link = "logit") else gaussian()
    nullmod <- fitNullModel(
        scanAnnot,
        outcome = opt$trait,
        covars = covariates,
        cov.mat = kinship,
        family = fam,
        verbose = opt$verbose
    )
    cat("  Fitted", if (opt$model == "binary") "logistic" else "linear",
        "mixed model (GENESIS", if (is.null(kinship)) "- no kinship)" else "with kinship)", "\n")
    null_ops <- genesis_ops(nullmod, n_samples)
}

# Samples the null model was fitted on: must be the complete cases counted above, or the score
# test would pair residuals with the wrong genotype rows
model_n <- length(null_ops$PY)
acct("model_n", n_samples, model_n, paste(opt$model, "null model"))
qc_assert("adapter_model_n", model_n == n_null, n = model_n,
          message = sprintf("%s: null model on %d samples; %d complete cases", qc_where, model_n, n_null))

# Save null model
saveRDS(nullmod, paste0(opt$output_prefix, ".null_model.rds"))

# ============================================================================
# Run Association Tests
# ============================================================================
cat("\nRunning association tests...\n")
cat("  Variants:", n_variants, "\n")
cat("  Chunk size:", opt$chunk_size, "\n")

# Process variants
n_status <- c(ok = 0L, mac = 0L, error = 0L)      # per-variant outcome, for <output_prefix>.run_counts.tsv
rows <- vector("list", n_variants)
pb_interval <- max(1, floor(n_variants / 20))
start_time <- Sys.time()

for (v in 1:n_variants) {
    # Get variant info
    var_info <- tractor_data$variant_info[v, ]

    # LA (non-reference hapcounts), dose and hapcount (every ancestry) for this variant
    dose_vec <- lapply(setNames(ancestries, ancestries), function(a) tractor_data$dosage[[a]][v, ])
    hap_vec  <- lapply(setNames(ancestries, ancestries), function(a) tractor_data$hapcount[[a]][v, ])
    la_vec   <- hap_vec[non_ref_anc]

    # Run test (same conditional score test for binary / quantitative / survival)
    res <- test_variant_tractor(
        la_vec = la_vec,
        dose_vec = dose_vec,
        hap_vec = hap_vec,
        null_ops = null_ops,
        ancestries = ancestries,
        ref_anc = ref_anc,
        mac_min = opt$mac_min
    )

    n_status[[res$status]] <- n_status[[res$status]] + 1L

    rows[[v]] <- result_row(var_info, res, ancestries, opt$stratum)

    # Progress
    if (v %% pb_interval == 0 || v == n_variants) {
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        rate <- v / elapsed
        eta <- (n_variants - v) / rate
        cat(sprintf("\r  Progress: %d/%d (%.1f%%) | %.0f var/s | ETA: %.0fs    ",
                    v, n_variants, 100 * v / n_variants, rate, eta))
    }
}
cat("\n")
results <- rbindlist(rows)       # STRATUM (with --stratum) is the first column (result_row)

# ============================================================================
# Write Results
# ============================================================================
cat("\nWriting results...\n")

# Sort results for different output files
# 1. Genomic order (for LD, fine-mapping, visualization)
# 2. P-value order (for quick interpretation)


results_genomic <- copy(results)
results_genomic[, CHR_NUM := as.integer(gsub("chr|X|Y|M", "", CHR, ignore.case = TRUE))]
results_genomic[is.na(CHR_NUM), CHR_NUM := 99]  # Handle X, Y, MT
setorder(results_genomic, CHR_NUM, POS)
results_genomic[, CHR_NUM := NULL]

results_pval <- copy(results)
setorder(results_pval, P_JOINT, na.last = TRUE)

# Main results file (SORTED BY P_JOINT - most significant first)
out_file <- paste0(opt$output_prefix, ".tractor_genesis.tsv.gz")
fwrite_sumstats(results_pval, out_file, sumstats_header)
cat("  Main results (sorted by P_JOINT):", out_file, "\n")

# Genomic order file (for downstream analysis)
genomic_file <- paste0(opt$output_prefix, ".tractor_genesis.genomic_order.tsv.gz")
fwrite_sumstats(results_genomic, genomic_file, sumstats_header)
cat("  Genomic order (for LD/fine-mapping):", genomic_file, "\n")

# Top hits file (top 1000 by P_JOINT)
n_top <- min(1000, nrow(results_pval))
top_file <- paste0(opt$output_prefix, ".tractor_genesis.top_hits.tsv")
fwrite_sumstats(results_pval[1:n_top], top_file, sumstats_header)
cat("  Top", n_top, "hits:", top_file, "\n")

# Significant hits (P_JOINT < 1e-5 OR any P_{anc} < 1e-5) - sorted by P_JOINT
sig_cols <- c("P_JOINT", paste0("P_", ancestries))
is_sig <- apply(results_pval[, ..sig_cols], 1, function(x) any(x < 1e-5, na.rm = TRUE))
if (sum(is_sig) > 0) {
    sig_results <- results_pval[is_sig]
    setorder(sig_results, P_JOINT, na.last = TRUE)
    sig_file <- paste0(opt$output_prefix, ".tractor_genesis.significant.tsv")
    fwrite_sumstats(sig_results, sig_file, sumstats_header)
    cat("  Significant hits (P < 1e-5, sorted):", sig_file, "(", sum(is_sig), "variants)\n")
}

# Heterogeneous effects (P_HET < 0.05 AND P_JOINT < 1e-5) - sorted by P_HET
is_het <- results_pval$P_HET < 0.05 & results_pval$P_JOINT < 1e-5
is_het[is.na(is_het)] <- FALSE
if (sum(is_het) > 0) {
    het_results <- results_pval[is_het]
    setorder(het_results, P_HET, na.last = TRUE)
    het_file <- paste0(opt$output_prefix, ".tractor_genesis.heterogeneous.tsv")
    fwrite_sumstats(het_results, het_file, sumstats_header)
    cat("  Heterogeneous effects (sorted by P_HET):", het_file, "(", sum(is_het), "variants)\n")
}

# Ancestry-divergent effects (significant in one ancestry but not others)
cat("\n  Checking for ancestry-divergent effects...\n")
for (anc in ancestries) {
    p_col <- paste0("P_", anc)
    other_p_cols <- paste0("P_", setdiff(ancestries, anc))

    # Significant in this ancestry (P < 5e-8) but not others (P > 0.05)
    is_anc_specific <- results_pval[[p_col]] < 5e-8
    for (other_col in other_p_cols) {
        is_anc_specific <- is_anc_specific & (results_pval[[other_col]] > 0.05 | is.na(results_pval[[other_col]]))
    }
    is_anc_specific[is.na(is_anc_specific)] <- FALSE

    if (sum(is_anc_specific) > 0) {
        anc_results <- results_pval[is_anc_specific]
        setorder(anc_results, get(p_col), na.last = TRUE)
        anc_file <- paste0(opt$output_prefix, ".tractor_genesis.", anc, "_specific.tsv")
        fwrite_sumstats(anc_results, anc_file, sumstats_header)
        cat("    ", anc, "-specific effects:", sum(is_anc_specific), "variants →", anc_file, "\n")
    }
}

# ============================================================================
# Summary Statistics
# ============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                          RESULTS SUMMARY                         ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║ Variants tested: %-47d ║\n", nrow(results)))
cat(sprintf("║ Samples analyzed: %-46d ║\n", n_samples))
cat("║                                                                  ║\n")

# Joint test summary
n_joint_sig <- sum(results$P_JOINT < 5e-8, na.rm = TRUE)
n_joint_sug <- sum(results$P_JOINT < 1e-5, na.rm = TRUE)
cat(sprintf("║ Joint test (any ancestry effect):                                ║\n"))
cat(sprintf("║   Genome-wide significant (P < 5e-8): %-26d ║\n", n_joint_sig))
cat(sprintf("║   Suggestive (P < 1e-5): %-40d ║\n", n_joint_sug))
cat("║                                                                  ║\n")

# Per-ancestry summary
cat("║ Ancestry-specific effects (P < 5e-8):                           ║\n")
for (anc in ancestries) {
    p_col <- paste0("P_", anc)
    n_sig <- sum(results[[p_col]] < 5e-8, na.rm = TRUE)
    cat(sprintf("║   %-5s: %-55d ║\n", anc, n_sig))
}
cat("║                                                                  ║\n")

# Heterogeneity summary
n_het <- sum(results$P_HET < 0.05 & results$P_JOINT < 1e-5, na.rm = TRUE)
cat(sprintf("║ Heterogeneous effects (P_HET < 0.05, P_JOINT < 1e-5): %-10d ║\n", n_het))

cat("╚══════════════════════════════════════════════════════════════════╝\n")

# ============================================================================
# Output Interpretation Guide
# ============================================================================
cat("\n")
cat("OUTPUT COLUMN INTERPRETATION:\n")
cat("─────────────────────────────────────────────────────────────────────\n")
cat("  P_JOINT      : Does this SNP have ANY effect? (", n_anc, "-df test)\n", sep = "")
cat("                 Significant = SNP affects trait through at least one ancestry\n")
cat("\n")
for (anc in ancestries) {
    cat("  BETA_", anc, sprintf("%*s", 6 - nchar(anc), ""), ": Effect size when allele is on ", anc, " haplotype\n", sep = "")
    cat("  P_", anc, sprintf("%*s", 9 - nchar(anc), ""), ": Is there an effect on ", anc, " background specifically?\n", sep = "")
}
cat("\n")
cat("  P_HET        : Do effects DIFFER across ancestries?\n")
cat("                 Significant = Effect size varies by ancestral background\n")
cat("  I2           : Heterogeneity magnitude (0-100%)\n")
cat("                 >50% = substantial heterogeneity, >75% = considerable\n")
cat("─────────────────────────────────────────────────────────────────────\n")

# Counts for the run record (bin/run_record.sh reads <output_prefix>.run_counts.tsv)
run_counts$n_variants <- n_variants
run_counts$n_pass_mac <- n_status[["ok"]] + n_status[["error"]]
run_counts$n_tested   <- sum(!is.na(results$P_JOINT))
run_counts$n_failed   <- n_status[["error"]]
cat("  Counts:", write_run_counts(opt$output_prefix, run_counts), "\n")

cat("\nDone!\n")
