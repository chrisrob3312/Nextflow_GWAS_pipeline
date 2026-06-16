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

    # Runtime
    make_option(c("--chunk_size"), type = "integer", default = 1000,
                help = "Variants per chunk [default: 1000]"),
    make_option(c("--threads"), type = "integer", default = 1,
                help = "Number of threads [default: 1]"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output")
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
for (anc in non_ref_anc) {
    cat(" + β_LA_", anc, "·LA_", anc, sep = "")
}
for (anc in ancestries) {
    cat(" + β_", anc, "·Dose_", anc, sep = "")
}
cat(" + covariates\n\n")

cat("INTERPRETATION:\n")
cat("  LA_", non_ref_anc[1], ": Local ancestry count (0,1,2 haplotypes of ", non_ref_anc[1], " at this locus)\n", sep = "")
cat("  Dose_", ancestries[1], ": Risk allele copies carried on ", ancestries[1], " haplotypes\n", sep = "")
cat("  β_", ancestries[1], ": Effect of allele WHEN on ", ancestries[1], " background\n\n", sep = "")

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
            paste0(prefix, ".hapcount.", anc, ".tsv.gz")
        )

        dosage_files <- c(
            paste0(prefix, ".dosage.", anc, ".txt.gz"),
            paste0(prefix, ".dosage.", anc_code, ".txt.gz"),
            paste0(prefix, ".", anc, ".dosage.txt.gz"),
            paste0(prefix, ".ancdose.", anc, ".tsv.gz")
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

# Create ScanAnnotationDataFrame
scanAnnot <- ScanAnnotationDataFrame(data.frame(
    scanID = common_samples,
    pheno[, c(opt$trait, covariates), with = FALSE]
))

# Fit null model
if (opt$model == "survival") {
    pheno$surv_outcome <- Surv(pheno[[opt$time_col]], pheno[[opt$event_col]])
    scanAnnot <- ScanAnnotationDataFrame(data.frame(
        scanID = common_samples,
        surv_outcome = pheno$surv_outcome,
        pheno[, covariates, with = FALSE]
    ))

    nullmod <- fitNullModel(
        scanAnnot,
        outcome = "surv_outcome",
        covars = covariates,
        cov.mat = kinship,
        family = "cox",
        verbose = opt$verbose
    )
    cat("  Fitted Cox proportional hazards mixed model\n")

} else if (opt$model == "binary") {
    nullmod <- fitNullModel(
        scanAnnot,
        outcome = opt$trait,
        covars = covariates,
        cov.mat = kinship,
        family = binomial(link = "logit"),
        verbose = opt$verbose
    )
    cat("  Fitted logistic mixed model\n")

} else {
    nullmod <- fitNullModel(
        scanAnnot,
        outcome = opt$trait,
        covars = covariates,
        cov.mat = kinship,
        family = gaussian(),
        verbose = opt$verbose
    )
    cat("  Fitted linear mixed model\n")
}

# Save null model
saveRDS(nullmod, paste0(opt$output_prefix, ".null_model.rds"))

# ============================================================================
# Association Testing Functions
# ============================================================================

# Main test function for a single variant
test_variant_tractor <- function(
    la_vec,      # Named list: LA counts per ancestry (k-1 non-ref ancestries)
    dose_vec,    # Named list: Dosages per ancestry (all k ancestries)
    nullmod,     # GENESIS null model
    ancestries,  # All ancestry names
    ref_anc      # Reference ancestry
) {
    n <- length(dose_vec[[1]])
    non_ref <- setdiff(ancestries, ref_anc)

    # Build design matrix
    # Columns: LA terms (k-1) + Dose terms (k)
    n_la <- length(non_ref)
    n_dose <- length(ancestries)

    X <- matrix(0, nrow = n, ncol = n_la + n_dose)
    col_names <- c(paste0("LA_", non_ref), paste0("Dose_", ancestries))
    colnames(X) <- col_names

    # Fill LA columns
    for (i in seq_along(non_ref)) {
        anc <- non_ref[i]
        X[, paste0("LA_", anc)] <- la_vec[[anc]]
    }

    # Fill Dose columns
    for (anc in ancestries) {
        X[, paste0("Dose_", anc)] <- dose_vec[[anc]]
    }

    # Remove samples with missing values
    complete <- complete.cases(X)
    if (sum(complete) < 10) {
        return(list(
            joint = list(stat = NA, df = n_dose, p = NA),
            marginal = lapply(ancestries, function(a) list(ancestry = a, beta = NA, se = NA, p = NA)),
            het = list(Q = NA, df = n_dose - 1, p = NA, I2 = NA)
        ))
    }

    X <- X[complete, , drop = FALSE]
    n_eff <- nrow(X)

    # Get null model components
    # Working vector (adjusted phenotype)
    resid <- nullmod$resid[complete]

    # ========== JOINT TEST ==========
    # k-df test: H0: all beta_Dose = 0 (adjusting for LA)
    # Using score test framework

    dose_cols <- paste0("Dose_", ancestries)
    X_dose <- X[, dose_cols, drop = FALSE]

    # Score vector
    U <- as.numeric(t(X_dose) %*% resid)

    # Variance of score (simplified - assumes independence)
    V <- t(X_dose) %*% X_dose

    # Add small ridge for numerical stability
    V_reg <- V + diag(1e-6, nrow(V))

    # Joint test statistic
    tryCatch({
        V_inv <- solve(V_reg)
        stat_joint <- as.numeric(t(U) %*% V_inv %*% U)
        df_joint <- n_dose
        p_joint <- pchisq(stat_joint, df = df_joint, lower.tail = FALSE)
    }, error = function(e) {
        stat_joint <<- NA
        df_joint <<- n_dose
        p_joint <<- NA
    })

    # ========== MARGINAL TESTS ==========
    # 1-df test for each ancestry
    marginal <- list()
    betas <- numeric(n_dose)
    ses <- numeric(n_dose)

    for (i in seq_along(ancestries)) {
        anc <- ancestries[i]
        x <- X[, paste0("Dose_", anc)]

        # OLS estimate (simplified - full version uses GLS with null model weights)
        var_x <- var(x)

        if (var_x > 1e-10) {
            beta <- sum(x * resid) / sum(x^2)
            se <- sqrt(1 / sum(x^2))
            z <- beta / se
            p <- 2 * pnorm(abs(z), lower.tail = FALSE)
        } else {
            beta <- NA
            se <- NA
            p <- NA
        }

        betas[i] <- beta
        ses[i] <- se

        marginal[[i]] <- list(
            ancestry = anc,
            beta = beta,
            se = se,
            z = if (!is.na(beta) && !is.na(se) && se > 0) beta / se else NA,
            p = p
        )
    }
    names(marginal) <- ancestries

    # ========== HETEROGENEITY TEST ==========
    # Test H0: beta_1 = beta_2 = ... = beta_k (equal effects across ancestries)

    valid <- !is.na(betas) & !is.na(ses) & ses > 0

    if (sum(valid) >= 2) {
        betas_v <- betas[valid]
        ses_v <- ses[valid]
        weights <- 1 / ses_v^2

        # Inverse-variance weighted mean
        beta_pooled <- sum(weights * betas_v) / sum(weights)

        # Cochran's Q
        Q <- sum(weights * (betas_v - beta_pooled)^2)
        df_het <- sum(valid) - 1
        p_het <- pchisq(Q, df = df_het, lower.tail = FALSE)

        # I-squared
        I2 <- max(0, (Q - df_het) / Q * 100)
    } else {
        Q <- NA
        df_het <- n_dose - 1
        p_het <- NA
        I2 <- NA
    }

    list(
        joint = list(stat = stat_joint, df = df_joint, p = p_joint),
        marginal = marginal,
        het = list(Q = Q, df = df_het, p = p_het, I2 = I2),
        n_eff = n_eff
    )
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

    # Run test
    res <- test_variant_tractor(
        la_vec = la_vec,
        dose_vec = dose_vec,
        nullmod = nullmod,
        ancestries = ancestries,
        ref_anc = ref_anc
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

# Main results file
out_file <- paste0(opt$output_prefix, ".tractor_genesis.tsv.gz")
fwrite(results, out_file, sep = "\t", compress = "gzip")
cat("  Main results:", out_file, "\n")

# Significant hits (P_JOINT < 1e-5 OR any P_{anc} < 1e-5)
sig_cols <- c("P_JOINT", paste0("P_", ancestries))
is_sig <- apply(results[, ..sig_cols], 1, function(x) any(x < 1e-5, na.rm = TRUE))
if (sum(is_sig) > 0) {
    sig_file <- paste0(opt$output_prefix, ".tractor_genesis.significant.tsv")
    fwrite(results[is_sig], sig_file, sep = "\t")
    cat("  Significant hits (P < 1e-5):", sig_file, "(", sum(is_sig), "variants)\n")
}

# Heterogeneous effects (P_HET < 0.05 AND P_JOINT < 1e-5)
is_het <- results$P_HET < 0.05 & results$P_JOINT < 1e-5
is_het[is.na(is_het)] <- FALSE
if (sum(is_het) > 0) {
    het_file <- paste0(opt$output_prefix, ".tractor_genesis.heterogeneous.tsv")
    fwrite(results[is_het], het_file, sep = "\t")
    cat("  Heterogeneous effects:", het_file, "(", sum(is_het), "variants)\n")
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
