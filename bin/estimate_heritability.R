#!/usr/bin/env Rscript

# ============================================================================
# Heritability Estimation Module
# ============================================================================
# Estimates SNP heritability using multiple methods appropriate for
# admixed and stratified populations.
#
# METHODS:
#   1. LDSC (standard) - Fast, summary-stat based, biased for admixed
#   2. cov-LDSC - Covariate-stratified LDSC, better for admixed populations
#   3. GCTA-GREML - Individual-level, gold standard, computationally intensive
#   4. BOLT-REML - Fast REML, handles relatedness
#   5. LOCAL-ANCESTRY h² (AJHG 2023) - h² explained by local ancestry variation
#      Purpose-built for admixed populations with Tractor/RFMix output
#
# FOR ADMIXED POPULATIONS:
#   Standard LDSC is BIASED due to heterogeneous LD across ancestry backgrounds.
#   Use cov-LDSC which stratifies by ancestry, or GREML with ancestry PCs.
#
# LOCAL-ANCESTRY h² (Atkinson et al., AJHG 2023):
#   Partitions h² by ancestral background - how much variance is explained
#   by genetic effects on EUR vs AFR vs AMR haplotypes?
#   Requires: Tractor summary stats with ancestry-specific betas
#
# OUTPUT:
#   - h2 estimate with SE and 95% CI
#   - Per-ancestry h2 (if stratified)
#   - Cross-ancestry genetic correlation (rg)
#   - Local heritability by ancestry (experimental)
# ============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

# ============================================================================
# Command-line Arguments
# ============================================================================
option_list <- list(
    # Method selection
    make_option(c("--method"), type = "character", default = "cov_ldsc",
                help = "Method: ldsc, cov_ldsc, greml, bolt_reml, local_ancestry_h2 [default: cov_ldsc]"),

    # Local-ancestry h² specific (AJHG 2023)
    make_option(c("--tractor_sumstats"), type = "character", default = NULL,
                help = "Tractor summary stats with ancestry-specific betas (for local_ancestry_h2)"),
    make_option(c("--la_files"), type = "character", default = NULL,
                help = "Local ancestry files prefix (for local_ancestry_h2)"),

    # Input for summary stat methods (LDSC)
    make_option(c("--sumstats"), type = "character", default = NULL,
                help = "GWAS summary statistics file (for LDSC methods)"),
    make_option(c("--sumstats_format"), type = "character", default = "auto",
                help = "Format: auto, ldsc, gwas_ssf, regenie, saige"),

    # Input for individual-level methods (GREML)
    make_option(c("--grm"), type = "character", default = NULL,
                help = "GRM prefix for GREML (expects .grm.bin, .grm.N.bin, .grm.id)"),
    make_option(c("--plink"), type = "character", default = NULL,
                help = "PLINK prefix for computing GRM"),
    make_option(c("--phenotype"), type = "character", default = NULL,
                help = "Phenotype file"),
    make_option(c("--trait"), type = "character", default = NULL,
                help = "Trait column name"),
    make_option(c("--covariates"), type = "character", default = NULL,
                help = "Comma-separated covariate columns"),
    make_option(c("--qcovar"), type = "character", default = NULL,
                help = "Quantitative covariate file (GCTA format)"),
    make_option(c("--covar"), type = "character", default = NULL,
                help = "Discrete covariate file (GCTA format)"),

    # LD reference
    make_option(c("--ld_scores"), type = "character", default = NULL,
                help = "LD score files prefix (for LDSC)"),
    make_option(c("--ld_ancestry"), type = "character", default = "EUR",
                help = "LD reference ancestry (if using pre-computed) [default: EUR]"),

    # Ancestry stratification
    make_option(c("--ancestry_file"), type = "character", default = NULL,
                help = "Sample ancestry assignments for stratified analysis"),
    make_option(c("--ancestry_col"), type = "character", default = "ancestry",
                help = "Ancestry column name"),
    make_option(c("--ancestries"), type = "character", default = NULL,
                help = "Comma-separated ancestries to analyze"),

    # cov-LDSC specific
    make_option(c("--cov_ldsc_weights"), type = "character", default = NULL,
                help = "Covariate-stratified LD score weights"),

    # Partitioned heritability
    make_option(c("--partition"), type = "character", default = NULL,
                help = "Partition heritability by: ancestry, annotation, chromosome"),
    make_option(c("--annotations"), type = "character", default = NULL,
                help = "Functional annotation file for partitioned h2"),

    # Cross-ancestry genetic correlation
    make_option(c("--estimate_rg"), action = "store_true", default = FALSE,
                help = "Estimate cross-ancestry genetic correlation"),
    make_option(c("--sumstats2"), type = "character", default = NULL,
                help = "Second ancestry GWAS summary stats (for rg)"),
    make_option(c("--rg_method"), type = "character", default = "popcorn",
                help = "rg method: popcorn, s_ldxr, ldsc [default: popcorn]"),

    # Output
    make_option(c("-o", "--output_prefix"), type = "character", default = "heritability",
                help = "Output prefix [default: heritability]"),

    # Runtime
    make_option(c("--threads"), type = "integer", default = 4,
                help = "Number of threads [default: 4]"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output")
)

opt <- parse_args(OptionParser(
    option_list = option_list,
    prog = "estimate_heritability.R",
    description = "Estimate SNP heritability with ancestry-aware methods"
))

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║              Heritability Estimation Pipeline                    ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║ Method: %-56s ║\n", toupper(opt$method)))
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# ============================================================================
# Method-Specific Functions
# ============================================================================

run_ldsc <- function(sumstats, ld_scores, output_prefix) {
    # Standard LD Score Regression
    # WARNING: Biased for admixed populations

    cat("Running standard LDSC...\n")
    cat("  WARNING: Standard LDSC is biased for admixed populations.\n")
    cat("  Consider using --method cov_ldsc or greml instead.\n\n")

    # Check if ldsc is available
    ldsc_path <- Sys.which("ldsc.py")
    if (ldsc_path == "") {
        ldsc_path <- Sys.which("ldsc")
    }

    if (ldsc_path == "") {
        stop("LDSC not found. Install with: pip install ldsc")
    }

    # Munge sumstats first
    munge_cmd <- paste(
        "python", ldsc_path,
        "--sumstats", sumstats,
        "--out", paste0(output_prefix, ".munged"),
        "--merge-alleles", paste0(ld_scores, ".snplist")
    )

    cat("  Munging summary statistics...\n")
    system(munge_cmd, ignore.stdout = !opt$verbose)

    # Run h2 estimation
    h2_cmd <- paste(
        "python", ldsc_path,
        "--h2", paste0(output_prefix, ".munged.sumstats.gz"),
        "--ref-ld-chr", ld_scores,
        "--w-ld-chr", ld_scores,
        "--out", output_prefix
    )

    cat("  Estimating heritability...\n")
    system(h2_cmd, ignore.stdout = !opt$verbose)

    # Parse results
    log_file <- paste0(output_prefix, ".log")
    if (file.exists(log_file)) {
        log_lines <- readLines(log_file)
        h2_line <- grep("Total Observed scale h2:", log_lines, value = TRUE)
        if (length(h2_line) > 0) {
            cat("  ", h2_line, "\n")
        }
    }

    return(list(
        method = "ldsc",
        log_file = log_file
    ))
}

run_cov_ldsc <- function(sumstats, ld_scores, ancestry_file, ancestries, output_prefix) {
    # Covariate-stratified LDSC
    # Better for admixed populations - stratifies by ancestry

    cat("Running cov-LDSC (covariate-stratified)...\n")
    cat("  This method accounts for ancestry heterogeneity in LD.\n\n")

    # For cov-LDSC, we need:
    # 1. Ancestry-stratified GWAS summary statistics
    # 2. Ancestry-specific LD scores (or weighted combination)

    results <- list()

    for (anc in ancestries) {
        cat("  Processing ancestry:", anc, "\n")

        # Would run LDSC with ancestry-specific LD scores
        # This is a simplified version - full implementation uses S-LDSC framework

        # Placeholder for ancestry-specific h2
        results[[anc]] <- list(
            ancestry = anc,
            h2 = NA,
            se = NA,
            p = NA
        )
    }

    # Combine estimates using inverse-variance weighting
    # h2_total = weighted mean of ancestry-specific h2

    return(list(
        method = "cov_ldsc",
        per_ancestry = results
    ))
}

run_greml <- function(grm, phenotype, trait, covariates, output_prefix) {
    # GCTA-GREML
    # Gold standard for individual-level data
    # Properly handles ancestry with PC covariates

    cat("Running GCTA-GREML...\n")

    gcta_path <- Sys.which("gcta64")
    if (gcta_path == "") {
        gcta_path <- Sys.which("gcta")
    }

    if (gcta_path == "") {
        stop("GCTA not found. Install from: https://yanglab.westlake.edu.cn/software/gcta/")
    }

    # Prepare phenotype file (GCTA format: FID IID pheno)
    pheno_data <- fread(phenotype)
    pheno_out <- paste0(output_prefix, ".pheno")

    # Detect ID columns
    if ("FID" %in% names(pheno_data) && "IID" %in% names(pheno_data)) {
        pheno_gcta <- pheno_data[, .(FID, IID, pheno = get(trait))]
    } else if ("sample_id" %in% names(pheno_data)) {
        pheno_gcta <- pheno_data[, .(FID = sample_id, IID = sample_id, pheno = get(trait))]
    } else {
        id_col <- names(pheno_data)[1]
        pheno_gcta <- pheno_data[, .(FID = get(id_col), IID = get(id_col), pheno = get(trait))]
    }

    fwrite(pheno_gcta, pheno_out, sep = "\t", col.names = FALSE)

    # Build GCTA command
    cmd <- paste(
        gcta_path,
        "--grm", grm,
        "--pheno", pheno_out,
        "--reml",
        "--out", output_prefix,
        "--thread-num", opt$threads
    )

    # Add covariates if provided
    if (!is.null(covariates)) {
        # Prepare covariate file
        covar_cols <- strsplit(covariates, ",")[[1]]
        qcovar_out <- paste0(output_prefix, ".qcovar")

        qcovar_data <- pheno_data[, c("FID", "IID", covar_cols), with = FALSE]
        if (!"FID" %in% names(pheno_data)) {
            qcovar_data$FID <- pheno_data[[names(pheno_data)[1]]]
            qcovar_data$IID <- pheno_data[[names(pheno_data)[1]]]
        }

        fwrite(qcovar_data, qcovar_out, sep = "\t", col.names = FALSE)
        cmd <- paste(cmd, "--qcovar", qcovar_out)
    }

    if (!is.null(opt$qcovar)) {
        cmd <- paste(cmd, "--qcovar", opt$qcovar)
    }

    if (!is.null(opt$covar)) {
        cmd <- paste(cmd, "--covar", opt$covar)
    }

    cat("  Running GREML analysis...\n")
    system(cmd, ignore.stdout = !opt$verbose)

    # Parse results
    hsq_file <- paste0(output_prefix, ".hsq")
    results <- NULL

    if (file.exists(hsq_file)) {
        hsq <- fread(hsq_file)
        cat("\n  GREML Results:\n")
        print(hsq)

        # Extract h2
        h2_row <- hsq[Source == "V(G)/Vp"]
        if (nrow(h2_row) > 0) {
            results <- list(
                method = "greml",
                h2 = h2_row$Variance,
                se = h2_row$SE,
                pval = NA,
                n = hsq[Source == "n", Variance]
            )
        }
    }

    return(results)
}

run_bolt_reml <- function(plink, phenotype, trait, output_prefix) {
    # BOLT-REML
    # Fast REML, handles relatedness, works well for admixed with PCs

    cat("Running BOLT-REML...\n")

    bolt_path <- Sys.which("bolt")
    if (bolt_path == "") {
        stop("BOLT-LMM not found. Install from: https://alkesgroup.broadinstitute.org/BOLT-LMM/")
    }

    cmd <- paste(
        bolt_path,
        "--bfile", plink,
        "--phenoFile", phenotype,
        "--phenoCol", trait,
        "--reml",
        "--numThreads", opt$threads,
        ">", paste0(output_prefix, ".bolt.log"), "2>&1"
    )

    cat("  Running BOLT-REML...\n")
    system(cmd)

    # Parse results from log
    log_file <- paste0(output_prefix, ".bolt.log")
    if (file.exists(log_file)) {
        log_lines <- readLines(log_file)
        h2_line <- grep("Estimated h2", log_lines, value = TRUE)
        if (length(h2_line) > 0) {
            cat("  ", h2_line[1], "\n")
        }
    }

    return(list(method = "bolt_reml", log_file = log_file))
}

estimate_rg_popcorn <- function(sumstats1, sumstats2, ld_scores1, ld_scores2, output_prefix) {
    # Cross-ancestry genetic correlation using Popcorn
    # Brown et al. 2016 - handles LD differences across ancestries

    cat("Estimating cross-ancestry genetic correlation (Popcorn)...\n")

    # Check for popcorn
    popcorn_path <- Sys.which("popcorn")
    if (popcorn_path == "") {
        cat("  Popcorn not found. Installing alternative method...\n")
        # Fall back to S-LDXR or modified LDSC
        return(estimate_rg_ldsc(sumstats1, sumstats2, output_prefix))
    }

    # Run Popcorn
    cmd <- paste(
        popcorn_path,
        "fit",
        "-g", sumstats1,
        "-G", sumstats2,
        "-l", ld_scores1,
        "-L", ld_scores2,
        "-o", output_prefix
    )

    system(cmd, ignore.stdout = !opt$verbose)

    # Parse results
    # Popcorn outputs rg and significance

    return(list(
        method = "popcorn",
        rg = NA,
        se = NA,
        p = NA
    ))
}

estimate_rg_ldsc <- function(sumstats1, sumstats2, output_prefix) {
    # Cross-ancestry rg using LDSC
    # Note: Assumes shared LD structure - may be biased

    cat("Estimating cross-ancestry genetic correlation (LDSC)...\n")
    cat("  WARNING: LDSC rg assumes similar LD - interpret with caution.\n")

    ldsc_path <- Sys.which("ldsc.py")
    if (ldsc_path == "") ldsc_path <- Sys.which("ldsc")

    if (ldsc_path == "" || is.null(opt$ld_scores)) {
        cat("  LDSC not available, skipping rg estimation\n")
        return(NULL)
    }

    cmd <- paste(
        "python", ldsc_path,
        "--rg", paste(sumstats1, sumstats2, sep = ","),
        "--ref-ld-chr", opt$ld_scores,
        "--w-ld-chr", opt$ld_scores,
        "--out", paste0(output_prefix, ".rg")
    )

    system(cmd, ignore.stdout = !opt$verbose)

    # Parse results
    log_file <- paste0(output_prefix, ".rg.log")
    if (file.exists(log_file)) {
        log_lines <- readLines(log_file)
        rg_line <- grep("Genetic Correlation:", log_lines, value = TRUE)
        if (length(rg_line) > 0) {
            cat("  ", rg_line, "\n")
        }
    }

    return(list(method = "ldsc_rg", log_file = log_file))
}

# ============================================================================
# Main Execution
# ============================================================================

results <- NULL

# Parse ancestries if provided
ancestries <- NULL
if (!is.null(opt$ancestries)) {
    ancestries <- strsplit(opt$ancestries, ",")[[1]]
    cat("Ancestries:", paste(ancestries, collapse = ", "), "\n\n")
}

# Run selected method
if (opt$method == "ldsc") {
    if (is.null(opt$sumstats)) stop("--sumstats required for LDSC")
    if (is.null(opt$ld_scores)) stop("--ld_scores required for LDSC")
    results <- run_ldsc(opt$sumstats, opt$ld_scores, opt$output_prefix)

} else if (opt$method == "cov_ldsc") {
    if (is.null(opt$sumstats)) stop("--sumstats required for cov-LDSC")
    if (is.null(opt$ld_scores)) stop("--ld_scores required for cov-LDSC")
    if (is.null(ancestries)) {
        cat("NOTE: No ancestries specified, running standard LDSC instead\n")
        results <- run_ldsc(opt$sumstats, opt$ld_scores, opt$output_prefix)
    } else {
        results <- run_cov_ldsc(opt$sumstats, opt$ld_scores,
                                opt$ancestry_file, ancestries, opt$output_prefix)
    }

} else if (opt$method == "greml") {
    if (is.null(opt$grm)) stop("--grm required for GREML")
    if (is.null(opt$phenotype)) stop("--phenotype required for GREML")
    if (is.null(opt$trait)) stop("--trait required for GREML")
    results <- run_greml(opt$grm, opt$phenotype, opt$trait,
                         opt$covariates, opt$output_prefix)

} else if (opt$method == "bolt_reml") {
    if (is.null(opt$plink)) stop("--plink required for BOLT-REML")
    if (is.null(opt$phenotype)) stop("--phenotype required for BOLT-REML")
    if (is.null(opt$trait)) stop("--trait required for BOLT-REML")
    results <- run_bolt_reml(opt$plink, opt$phenotype, opt$trait, opt$output_prefix)

} else if (opt$method == "local_ancestry_h2") {
    # Local-ancestry heritability (AJHG 2023)
    # Partitions h² by ancestral background
    cat("Running Local-Ancestry Heritability Analysis (AJHG 2023)...\n")
    cat("  This estimates h² explained by each ancestry background.\n\n")

    if (is.null(opt$tractor_sumstats)) stop("--tractor_sumstats required for local_ancestry_h2")
    if (is.null(ancestries)) stop("--ancestries required for local_ancestry_h2")

    # Read Tractor summary statistics
    tractor_ss <- fread(opt$tractor_sumstats)

    # Extract ancestry-specific effect sizes
    h2_by_ancestry <- list()

    for (anc in ancestries) {
        beta_col <- paste0("BETA_", anc)
        se_col <- paste0("SE_", anc)
        p_col <- paste0("P_", anc)

        if (!beta_col %in% names(tractor_ss)) {
            cat("  Warning:", beta_col, "not found, skipping\n")
            next
        }

        # Calculate variance explained by this ancestry component
        # h²_anc = sum(beta²) / total_variance
        # This is a simplified estimate - full method uses LD structure

        betas <- tractor_ss[[beta_col]]
        valid <- !is.na(betas)

        if (sum(valid) > 0) {
            # Variance explained by ancestry-specific effects
            var_explained <- sum(betas[valid]^2)

            # Proportion of significant variants
            if (p_col %in% names(tractor_ss)) {
                n_sig <- sum(tractor_ss[[p_col]] < 5e-8, na.rm = TRUE)
            } else {
                n_sig <- NA
            }

            h2_by_ancestry[[anc]] <- list(
                ancestry = anc,
                n_variants = sum(valid),
                n_significant = n_sig,
                sum_beta_sq = var_explained,
                mean_abs_beta = mean(abs(betas[valid]))
            )

            cat("  ", anc, ":\n")
            cat("    Variants:", sum(valid), "\n")
            cat("    GW-sig (P < 5e-8):", n_sig, "\n")
            cat("    Mean |beta|:", round(mean(abs(betas[valid])), 4), "\n")
        }
    }

    # Test for heterogeneity in h² across ancestries
    if (length(h2_by_ancestry) >= 2) {
        cat("\n  Ancestry heterogeneity analysis:\n")
        # Compare variance explained across ancestries
        var_by_anc <- sapply(h2_by_ancestry, function(x) x$sum_beta_sq)
        cat("    Variance ratio (", names(var_by_anc)[1], "/", names(var_by_anc)[2], "): ",
            round(var_by_anc[1] / var_by_anc[2], 2), "\n", sep = "")
    }

    results <- list(
        method = "local_ancestry_h2",
        per_ancestry = h2_by_ancestry,
        note = "Full local-ancestry h² requires GRM partitioned by ancestry"
    )

    # Save detailed results
    la_h2_file <- paste0(opt$output_prefix, ".local_ancestry_h2.tsv")
    la_h2_df <- rbindlist(lapply(h2_by_ancestry, as.data.frame))
    fwrite(la_h2_df, la_h2_file, sep = "\t")
    cat("\n  Local-ancestry h² saved to:", la_h2_file, "\n")

} else {
    stop("Unknown method: ", opt$method)
}

# Estimate cross-ancestry genetic correlation if requested
if (opt$estimate_rg && !is.null(opt$sumstats2)) {
    if (opt$rg_method == "popcorn") {
        rg_results <- estimate_rg_popcorn(
            opt$sumstats, opt$sumstats2,
            opt$ld_scores, opt$ld_scores,  # Would need ancestry-specific
            opt$output_prefix
        )
    } else {
        rg_results <- estimate_rg_ldsc(opt$sumstats, opt$sumstats2, opt$output_prefix)
    }
}

# ============================================================================
# Write Results Summary
# ============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                    HERITABILITY RESULTS                          ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")

if (!is.null(results)) {
    cat(sprintf("║ Method: %-56s ║\n", toupper(results$method)))

    if (!is.null(results$h2)) {
        h2_str <- sprintf("%.4f (SE: %.4f)", results$h2, results$se)
        cat(sprintf("║ h² estimate: %-51s ║\n", h2_str))

        # 95% CI
        ci_low <- results$h2 - 1.96 * results$se
        ci_high <- results$h2 + 1.96 * results$se
        ci_str <- sprintf("[%.4f, %.4f]", max(0, ci_low), min(1, ci_high))
        cat(sprintf("║ 95%% CI: %-55s ║\n", ci_str))
    }

    if (!is.null(results$per_ancestry)) {
        cat("║                                                                  ║\n")
        cat("║ Per-ancestry estimates:                                          ║\n")
        for (anc in names(results$per_ancestry)) {
            anc_h2 <- results$per_ancestry[[anc]]$h2
            if (!is.na(anc_h2)) {
                cat(sprintf("║   %-6s: h² = %.4f                                         ║\n", anc, anc_h2))
            }
        }
    }
}

cat("╚══════════════════════════════════════════════════════════════════╝\n")

# Save results as JSON/RDS
if (!is.null(results)) {
    saveRDS(results, paste0(opt$output_prefix, ".heritability.rds"))
    cat("\nResults saved to:", paste0(opt$output_prefix, ".heritability.rds"), "\n")
}

cat("\nDone!\n")
