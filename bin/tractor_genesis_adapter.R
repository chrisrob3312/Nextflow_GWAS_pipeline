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
#   Or from RFMix directly:
#     - {prefix}.msp.tsv.gz: Local ancestry calls (MSP format)
#     - VCF/GDS for genotypes
#
# ============================================================================

suppressPackageStartupMessages({
    library(GENESIS)
    library(SeqArray)
    library(SeqVarTools)
    library(survival)
    library(data.table)
    library(optparse)
    library(Matrix)
})

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
                help = "Comma-separated covariate column names"),
    make_option(c("--global_ancestry_prefix"), type = "character", default = NULL,
                help = "Prefix for global ancestry columns (e.g., 'prop_' expects prop_EUR, prop_AFR, etc.)"),

    # Output
    make_option(c("-o", "--output_prefix"), type = "character", default = "tractor_genesis",
                help = "Output file prefix [default: tractor_genesis]"),

    # Variant filters
    make_option(c("--mac_min"), type = "integer", default = 10,
                help = "Minimum minor allele count on an ancestry background to test that ancestry's dose term [default: 10]"),

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
                help = "Ancestry stratum to analyze (EUR, AAC, LAT1, LAT2, OTHER, or POOLED)"),
    make_option(c("--min_stratum_n"), type = "integer", default = 30,
                help = "Minimum samples per stratum [default: 30]"),
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

if (!is.null(opt$msp) && is.null(opt$gds)) {
    stop("--gds required when using --msp input")
}

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

# Function to read Tractor hapcount/dosage files
read_tractor_files <- function(prefix, ancestries) {
    hapcount <- list()
    dosage <- list()
    variant_info <- NULL

    for (i in seq_along(ancestries)) {
        anc <- ancestries[i]
        anc_code <- i - 1  # Tractor uses 0-indexed ancestry codes

        # Try different file naming conventions
        hapcount_files <- c(
            paste0(prefix, ".hapcount.", anc, ".txt.gz"),
            paste0(prefix, ".hapcount.", anc_code, ".txt.gz"),
            paste0(prefix, ".", anc, ".hapcount.txt.gz"),
            paste0(prefix, ".hapcount.", anc, ".tsv.gz"),
            paste0(prefix, ".hapcount.", anc_code, ".tsv.gz")   # TRACTOR_EXTRACT_TRACTS naming
        )

        dosage_files <- c(
            paste0(prefix, ".dosage.", anc, ".txt.gz"),
            paste0(prefix, ".dosage.", anc_code, ".txt.gz"),
            paste0(prefix, ".", anc, ".dosage.txt.gz"),
            paste0(prefix, ".ancdose.", anc, ".tsv.gz"),
            paste0(prefix, ".ancdose.", anc_code, ".tsv.gz")    # TRACTOR_EXTRACT_TRACTS naming
        )

        # Find hapcount file
        hc_file <- NULL
        for (f in hapcount_files) {
            if (file.exists(f)) {
                hc_file <- f
                break
            }
        }

        # Find dosage file
        dos_file <- NULL
        for (f in dosage_files) {
            if (file.exists(f)) {
                dos_file <- f
                break
            }
        }

        if (is.null(hc_file) || is.null(dos_file)) {
            stop("Cannot find Tractor files for ancestry '", anc, "'\n",
                 "  Tried: ", paste(hapcount_files, collapse = ", "))
        }

        cat("  Loading", anc, "...\n")
        cat("    Hapcount:", basename(hc_file), "\n")
        cat("    Dosage:", basename(dos_file), "\n")

        # Read files
        # Tractor format: first cols are variant info, rest are samples
        hc_data <- fread(hc_file)
        dos_data <- fread(dos_file)

        # Extract variant info from first file
        if (is.null(variant_info)) {
            # Identify variant info columns (typically CHR, POS, ID, REF, ALT or similar)
            info_cols <- intersect(names(hc_data),
                                   c("CHR", "CHROM", "POS", "ID", "SNP", "REF", "ALT",
                                     "chrom", "pos", "id", "snp", "ref", "alt",
                                     "#CHROM", "BP"))
            if (length(info_cols) == 0) {
                # Assume first 5 columns are variant info
                info_cols <- names(hc_data)[1:min(5, ncol(hc_data))]
            }
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
tractor_data <- read_tractor_files(opt$tractor_prefix, ancestries)

# Match samples between phenotype and Tractor data
common_samples <- intersect(pheno$sample_id, tractor_data$sample_ids)
cat("  Samples with both phenotype and genotype:", length(common_samples), "\n")

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

# Parse covariates
covariates <- NULL
if (!is.null(opt$covariates)) {
    covariates <- strsplit(opt$covariates, ",")[[1]]
}

# Add global ancestry proportions as covariates if specified
if (!is.null(opt$global_ancestry_prefix)) {
    ga_cols <- paste0(opt$global_ancestry_prefix, non_ref_anc)
    ga_cols <- intersect(ga_cols, names(pheno))
    if (length(ga_cols) > 0) {
        covariates <- c(covariates, ga_cols)
        cat("  Added global ancestry covariates:", paste(ga_cols, collapse = ", "), "\n")
    }
}

if (!is.null(covariates)) {
    cat("  Covariates:", paste(covariates, collapse = ", "), "\n")
}

# Subset kinship to common samples
if (!is.null(kinship)) {
    kin_samples <- intersect(rownames(kinship), common_samples)
    if (length(kin_samples) < length(common_samples)) {
        cat("  Warning: Kinship matrix missing", length(common_samples) - length(kin_samples), "samples\n")
    }
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

cov_mat <- if (length(covariates) > 0) {
    m <- as.matrix(pheno[, covariates, with = FALSE]); storage.mode(m) <- "double"; m
} else NULL

# --- helpers for GENESIS null-model geometry ---------------------------------
# GENESIS stores cholSigmaInv = L with SigmaInv = L L'; CX = L'X; CXCXI = CX (CX'CX)^-1
.as_L <- function(nm, n) {
    L <- nm$cholSigmaInv
    if (is.null(L)) stop("GENESIS null model has no cholSigmaInv; update GENESIS (>= 2.16)")
    if (is.matrix(L) || inherits(L, "Matrix")) return(L)
    if (length(L) == n) return(Matrix::Diagonal(n, x = as.numeric(L)))
    Matrix::Diagonal(n, x = rep(as.numeric(L)[1], n))
}
.crossL <- function(L, M) as.matrix(Matrix::crossprod(L, M))   # L'M
.multL  <- function(L, v) as.numeric(L %*% v)                    # L v

genesis_ops <- function(nm, n) {
    L <- .as_L(nm, n)
    CX <- as.matrix(nm$CX); CXCXI <- as.matrix(nm$CXCXI)
    PY <- nm$fit$resid.PY %||% nm$resid
    if (is.null(PY)) {
        Y  <- nm$fit$workingY %||% nm$workingY
        Yt <- .crossL(L, Y)
        PY <- .multL(L, Yt - CXCXI %*% crossprod(CX, Yt))
    }
    list(
        type = "genesis",
        PY = as.numeric(PY),
        # projected cross-product  D'PD  for a design matrix D (n x p)
        proj_crossprod = function(D) {
            Dt <- .crossL(L, D)
            crossprod(Dt) - crossprod(Dt, CXCXI) %*% crossprod(CX, Dt)
        }
    )
}

cox_ops <- function(time, event, eta, Xc) {
    w <- exp(eta)
    n <- length(time)
    # Breslow baseline hazard and martingale residuals under the null
    ord_asc <- order(time, -event)
    t_a <- time[ord_asc]; w_a <- w[ord_asc]; e_a <- event[ord_asc]
    S0_a <- rev(cumsum(rev(w_a)))                      # sum_{j: t_j >= t_i} w_j
    S0_a <- S0_a[match(t_a, t_a)]                      # ties share one risk set
    dL   <- ifelse(e_a == 1, 1 / S0_a, 0)
    Lam  <- ave(cumsum(dL), t_a, FUN = max)            # Lambda0(t_i) incl. all events at t_i
    mart <- numeric(n); mart[ord_asc] <- e_a - Lam * w_a

    # Risk-set machinery in DESCENDING time: risk set of t_k = positions 1..last(tie group)
    ord_d <- order(-time, event)
    t_d <- time[ord_d]; w_d <- w[ord_d]; e_d <- event[ord_d]
    last_idx <- ave(seq_along(t_d), t_d, FUN = max)
    ev_pos <- which(e_d == 1)
    Xc_d <- if (!is.null(Xc)) Xc[ord_d, , drop = FALSE] else NULL

    # Observed information (Breslow) for full design A = [D | Xc], then project on Xc
    info <- function(A_d) {
        p <- ncol(A_d)
        S0 <- cumsum(w_d)[last_idx]
        S1 <- apply(w_d * A_d, 2, cumsum)[last_idx, , drop = FALSE]
        I  <- matrix(0, p, p)
        for (a in 1:p) for (b in a:p) {
            S2 <- cumsum(w_d * A_d[, a] * A_d[, b])[last_idx]
            v  <- sum((S2 / S0 - (S1[, a] / S0) * (S1[, b] / S0))[ev_pos])
            I[a, b] <- v; I[b, a] <- v
        }
        I
    }
    q <- if (is.null(Xc_d)) 0 else ncol(Xc_d)
    list(
        type = "cox",
        PY = mart,
        proj_crossprod = function(D) {
            A <- if (q > 0) cbind(D[ord_d, , drop = FALSE], Xc_d) else D[ord_d, , drop = FALSE]
            I <- info(A)
            p <- ncol(D)
            if (q == 0) return(I)
            iD <- 1:p; iX <- (p + 1):(p + q)
            I[iD, iD] - I[iD, iX] %*% solve(I[iX, iX] + diag(1e-8, q)) %*% I[iX, iD]
        }
    )
}

# --- fit the null model ------------------------------------------------------
if (opt$model == "survival") {
    if (!all(c(opt$time_col, opt$event_col) %in% names(pheno))) {
        stop("Survival columns not found: ", opt$time_col, ", ", opt$event_col)
    }
    surv_time  <- as.numeric(pheno[[opt$time_col]])
    surv_event <- as.integer(pheno[[opt$event_col]])
    if (anyNA(surv_time) || anyNA(surv_event)) stop("Missing time/event values - filter the phenotype file first")

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

# Save null model
saveRDS(nullmod, paste0(opt$output_prefix, ".null_model.rds"))

# ============================================================================
# Association Testing Functions
# ============================================================================

# Main test function for a single variant
#
# Conditional score test (identical for binary / quantitative / survival):
#   D = [LA_1..LA_{k-1} | Dose_1..Dose_k]
#   U = D' P y ,  M = D' P D      (P, Py supplied by null_ops)
#   Dose block conditioned on the LA block (Schur complement):
#     U_c = U_G - M_GZ M_ZZ^-1 U_Z ,   V_c = M_GG - M_GZ M_ZZ^-1 M_ZG
#   Joint (k-df):    U_c' V_c^-1 U_c  ~ chi2_k       (H0: all dose effects 0)
#   Per-ancestry:    beta = V_c^-1 U_c (one-step estimate), se = sqrt(diag V_c^-1)
#   Heterogeneity:   contrasts C beta (beta_j - beta_1),  Q = (Cb)'(C V_c^-1 C')^-1 (Cb) ~ chi2_{k-1}
#                    - uses the full covariance, so correlated dose terms are handled
empty_result <- function(ancestries, n_dose, mac) list(
    joint = list(stat = NA, df = n_dose, p = NA),
    marginal = setNames(lapply(ancestries, function(a) list(ancestry = a, beta = NA, se = NA, z = NA, p = NA)), ancestries),
    het = list(Q = NA, df = n_dose - 1, p = NA, I2 = NA),
    n_eff = NA, mac = mac
)

test_variant_tractor <- function(
    la_vec,      # Named list: LA counts per ancestry (k-1 non-ref ancestries)
    dose_vec,    # Named list: Dosages per ancestry (all k ancestries)
    null_ops,    # list(PY, proj_crossprod) from the null model
    ancestries,  # All ancestry names
    ref_anc,     # Reference ancestry
    mac_min = 10
) {
    n <- length(dose_vec[[1]])
    non_ref <- setdiff(ancestries, ref_anc)
    n_dose <- length(ancestries)

    # Build design: LA terms then Dose terms; mean-impute sporadic missingness
    D <- matrix(0, nrow = n, ncol = length(non_ref) + n_dose)
    colnames(D) <- c(paste0("LA_", non_ref), paste0("Dose_", ancestries))
    for (anc in non_ref)    D[, paste0("LA_", anc)]   <- as.numeric(la_vec[[anc]])
    for (anc in ancestries) D[, paste0("Dose_", anc)] <- as.numeric(dose_vec[[anc]])
    n_missing <- sum(!complete.cases(D))
    if (n_missing > 0) {
        for (j in seq_len(ncol(D))) { miss <- is.na(D[, j]); if (any(miss)) D[miss, j] <- mean(D[, j], na.rm = TRUE) }
    }

    # Minor allele count per ancestry background
    mac <- sapply(ancestries, function(a) sum(D[, paste0("Dose_", a)]))
    names(mac) <- ancestries
    testable <- ancestries[mac >= mac_min & apply(D[, paste0("Dose_", ancestries), drop = FALSE], 2, var) > 1e-10]
    if (length(testable) == 0) return(empty_result(ancestries, n_dose, mac))

    # Drop LA columns with no variation (e.g. no AFR tracts in this chunk)
    la_cols <- paste0("LA_", non_ref)
    la_cols <- la_cols[apply(D[, la_cols, drop = FALSE], 2, var) > 1e-10]
    g_cols  <- paste0("Dose_", testable)
    Dk <- D[, c(la_cols, g_cols), drop = FALSE]

    res <- tryCatch({
        U <- as.numeric(crossprod(Dk, null_ops$PY))
        M <- as.matrix(null_ops$proj_crossprod(Dk))
        iZ <- seq_along(la_cols); iG <- length(la_cols) + seq_along(g_cols)

        if (length(iZ) > 0) {
            MZZi <- solve(M[iZ, iZ, drop = FALSE] + diag(1e-8, length(iZ)))
            U_c <- U[iG] - M[iG, iZ, drop = FALSE] %*% MZZi %*% U[iZ]
            V_c <- M[iG, iG, drop = FALSE] - M[iG, iZ, drop = FALSE] %*% MZZi %*% M[iZ, iG, drop = FALSE]
        } else { U_c <- U[iG]; V_c <- M[iG, iG, drop = FALSE] }
        V_c <- (V_c + t(V_c)) / 2
        V_inv <- solve(V_c + diag(1e-8, nrow(V_c)))

        # Joint test
        stat_joint <- as.numeric(t(U_c) %*% V_inv %*% U_c)
        df_joint <- length(iG)
        p_joint <- pchisq(stat_joint, df_joint, lower.tail = FALSE)

        # Per-ancestry one-step estimates
        beta <- as.numeric(V_inv %*% U_c); se <- sqrt(diag(V_inv))
        z <- beta / se; p <- 2 * pnorm(-abs(z))
        marginal <- setNames(lapply(ancestries, function(a) list(ancestry = a, beta = NA, se = NA, z = NA, p = NA)), ancestries)
        for (i in seq_along(testable)) {
            marginal[[testable[i]]] <- list(ancestry = testable[i], beta = beta[i], se = se[i], z = z[i], p = p[i])
        }

        # Heterogeneity across ancestry backgrounds (covariance-aware)
        if (length(testable) >= 2) {
            k <- length(testable)
            C <- cbind(-1, diag(k - 1))               # beta_j - beta_1
            Cb <- C %*% beta
            Q <- as.numeric(t(Cb) %*% solve(C %*% V_inv %*% t(C) + diag(1e-10, k - 1)) %*% Cb)
            df_het <- k - 1
            p_het <- pchisq(Q, df_het, lower.tail = FALSE)
            I2 <- max(0, (Q - df_het) / Q * 100)
        } else { Q <- NA; df_het <- n_dose - 1; p_het <- NA; I2 <- NA }

        list(joint = list(stat = stat_joint, df = df_joint, p = p_joint),
             marginal = marginal,
             het = list(Q = Q, df = df_het, p = p_het, I2 = I2),
             n_eff = n - n_missing, mac = mac)
    }, error = function(e) empty_result(ancestries, n_dose, mac))

    res
}

# ============================================================================
# Run Association Tests
# ============================================================================
cat("\nRunning association tests...\n")
cat("  Variants:", n_variants, "\n")
cat("  Chunk size:", opt$chunk_size, "\n")

# Initialize results storage
results <- data.table(
    CHR = character(),
    POS = integer(),
    ID = character(),
    REF = character(),
    ALT = character(),
    N = integer(),

    # Joint test
    STAT_JOINT = numeric(),
    DF_JOINT = integer(),
    P_JOINT = numeric()
)

# Add ancestry-specific columns dynamically
for (anc in ancestries) {
    results[[paste0("MAC_", anc)]] <- numeric()     # minor allele count on this background
    results[[paste0("BETA_", anc)]] <- numeric()
    results[[paste0("SE_", anc)]] <- numeric()
    results[[paste0("Z_", anc)]] <- numeric()
    results[[paste0("P_", anc)]] <- numeric()
}

# Heterogeneity columns
results$Q_HET <- numeric()
results$DF_HET <- integer()
results$P_HET <- numeric()
results$I2 <- numeric()

# Process variants
pb_interval <- max(1, floor(n_variants / 20))
start_time <- Sys.time()

for (v in 1:n_variants) {
    # Get variant info
    var_info <- tractor_data$variant_info[v, ]

    # Extract LA and Dose for this variant
    la_vec <- list()
    dose_vec <- list()

    for (anc in ancestries) {
        dose_vec[[anc]] <- tractor_data$dosage[[anc]][v, ]

        if (anc != ref_anc) {
            la_vec[[anc]] <- tractor_data$hapcount[[anc]][v, ]
        }
    }

    # Run test (same conditional score test for binary / quantitative / survival)
    res <- test_variant_tractor(
        la_vec = la_vec,
        dose_vec = dose_vec,
        null_ops = null_ops,
        ancestries = ancestries,
        ref_anc = ref_anc,
        mac_min = opt$mac_min
    )

    # Build result row
    row <- data.table(
        CHR = as.character(var_info$CHR),
        POS = as.integer(var_info$POS),
        ID = if ("ID" %in% names(var_info)) as.character(var_info$ID) else paste0(var_info$CHR, ":", var_info$POS),
        REF = if ("REF" %in% names(var_info)) as.character(var_info$REF) else NA_character_,
        ALT = if ("ALT" %in% names(var_info)) as.character(var_info$ALT) else NA_character_,
        N = res$n_eff,
        STAT_JOINT = res$joint$stat,
        DF_JOINT = res$joint$df,
        P_JOINT = res$joint$p
    )

    # Add ancestry-specific results
    for (anc in ancestries) {
        row[[paste0("MAC_", anc)]] <- as.numeric(res$mac[[anc]])
        row[[paste0("BETA_", anc)]] <- res$marginal[[anc]]$beta
        row[[paste0("SE_", anc)]] <- res$marginal[[anc]]$se
        row[[paste0("Z_", anc)]] <- res$marginal[[anc]]$z
        row[[paste0("P_", anc)]] <- res$marginal[[anc]]$p
    }

    # Add heterogeneity results
    row$Q_HET <- res$het$Q
    row$DF_HET <- res$het$df
    row$P_HET <- res$het$p
    row$I2 <- res$het$I2

    results <- rbind(results, row)

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
fwrite(results_pval, out_file, sep = "\t", compress = "gzip")
cat("  Main results (sorted by P_JOINT):", out_file, "\n")

# Genomic order file (for downstream analysis)
genomic_file <- paste0(opt$output_prefix, ".tractor_genesis.genomic_order.tsv.gz")
fwrite(results_genomic, genomic_file, sep = "\t", compress = "gzip")
cat("  Genomic order (for LD/fine-mapping):", genomic_file, "\n")

# Top hits file (top 1000 by P_JOINT)
n_top <- min(1000, nrow(results_pval))
top_file <- paste0(opt$output_prefix, ".tractor_genesis.top_hits.tsv")
fwrite(results_pval[1:n_top], top_file, sep = "\t")
cat("  Top", n_top, "hits:", top_file, "\n")

# Significant hits (P_JOINT < 1e-5 OR any P_{anc} < 1e-5) - sorted by P_JOINT
sig_cols <- c("P_JOINT", paste0("P_", ancestries))
is_sig <- apply(results_pval[, ..sig_cols], 1, function(x) any(x < 1e-5, na.rm = TRUE))
if (sum(is_sig) > 0) {
    sig_results <- results_pval[is_sig]
    setorder(sig_results, P_JOINT, na.last = TRUE)
    sig_file <- paste0(opt$output_prefix, ".tractor_genesis.significant.tsv")
    fwrite(sig_results, sig_file, sep = "\t")
    cat("  Significant hits (P < 1e-5, sorted):", sig_file, "(", sum(is_sig), "variants)\n")
}

# Heterogeneous effects (P_HET < 0.05 AND P_JOINT < 1e-5) - sorted by P_HET
is_het <- results_pval$P_HET < 0.05 & results_pval$P_JOINT < 1e-5
is_het[is.na(is_het)] <- FALSE
if (sum(is_het) > 0) {
    het_results <- results_pval[is_het]
    setorder(het_results, P_HET, na.last = TRUE)
    het_file <- paste0(opt$output_prefix, ".tractor_genesis.heterogeneous.tsv")
    fwrite(het_results, het_file, sep = "\t")
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
        fwrite(anc_results, anc_file, sep = "\t")
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

cat("\nDone!\n")
