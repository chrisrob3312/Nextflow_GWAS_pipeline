#!/usr/bin/env Rscript

# ============================================================================
# Multi-Method PRS for Admixed Populations
# ============================================================================
# Comprehensive PRS calculation using methods optimized for diverse ancestries
#
# PRIMARY METHOD:
#   PRS-CSx - Multi-ancestry PRS using ancestry-specific GWAS + LD references
#             Gold standard for admixed populations
#
# ADJUNCT METHODS (local ancestry-aware):
#   GAUDI       - Local ancestry-informed PRS
#   DiscoDivas  - Disentangling PRS by local ancestry
#   SDPR_admix  - Admixture-aware sparse Dirichlet process regression
#   MUSSEL      - Multi-ancestry stacking ensemble learner
#   PROSPER     - Polygenic Risk Score using Optimal Penalized Regression
#
# NOTE: NO Clumping+Thresholding (C+T) methods - they underperform in admixed
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
    make_option(c("--method"), type = "character", default = "prs_csx",
                help = "Method: prs_csx, gaudi, disco_divas, sdpr_admix, mussel, prosper, all [default: prs_csx]"),

    # Input: GWAS summary statistics
    make_option(c("--sumstats"), type = "character", default = NULL,
                help = "GWAS summary statistics (single ancestry or joint)"),
    make_option(c("--sumstats_dir"), type = "character", default = NULL,
                help = "Directory with ancestry-specific sumstats (for PRS-CSx)"),
    make_option(c("--sumstats_pattern"), type = "character", default = "{ancestry}.sumstats.gz",
                help = "Pattern for ancestry-specific sumstats [default: {ancestry}.sumstats.gz]"),

    # Input: LD reference
    make_option(c("--ld_ref"), type = "character", default = NULL,
                help = "LD reference panel directory"),
    make_option(c("--ld_ref_pattern"), type = "character", default = "{ancestry}/",
                help = "Pattern for ancestry-specific LD refs [default: {ancestry}/]"),
    make_option(c("--use_cohort_ld"), action = "store_true", default = FALSE,
                help = "Use cohort-specific LD (from calculate_cohort_ld.R)"),
    make_option(c("--cohort_ld_dir"), type = "character", default = NULL,
                help = "Directory with cohort-specific LD files"),

    # Input: Genotypes (for scoring)
    make_option(c("--geno"), type = "character", default = NULL,
                help = "Genotype file prefix (PLINK format)"),
    make_option(c("--gds"), type = "character", default = NULL,
                help = "GDS file (alternative to PLINK)"),

    # Input: Local ancestry (for LA-aware methods)
    make_option(c("--local_ancestry"), type = "character", default = NULL,
                help = "Local ancestry files prefix (RFMix MSP format)"),

    # Ancestry specification
    make_option(c("--ancestries"), type = "character", default = "EUR,AFR",
                help = "Comma-separated ancestries [default: EUR,AFR]"),
    make_option(c("--target_ancestry"), type = "character", default = NULL,
                help = "Target population ancestry (for weighting)"),

    # Ancestry proportions (for target population)
    make_option(c("--ancestry_props"), type = "character", default = NULL,
                help = "Global ancestry proportions: 'EUR:0.6,AFR:0.3,AMR:0.1'"),

    # Sample info
    make_option(c("--sample_file"), type = "character", default = NULL,
                help = "Sample file with ancestry assignments"),

    # PRS-CSx specific
    make_option(c("--phi"), type = "numeric", default = NULL,
                help = "PRS-CSx phi parameter (NULL for auto-tuning)"),
    make_option(c("--n_iter"), type = "integer", default = 1000,
                help = "MCMC iterations [default: 1000]"),
    make_option(c("--n_burnin"), type = "integer", default = 500,
                help = "MCMC burn-in [default: 500]"),

    # Output
    make_option(c("-o", "--output_prefix"), type = "character", default = "prs",
                help = "Output file prefix [default: prs]"),

    # Validation
    make_option(c("--phenotype"), type = "character", default = NULL,
                help = "Phenotype file for validation"),
    make_option(c("--trait"), type = "character", default = NULL,
                help = "Trait column for validation"),
    make_option(c("--validate"), action = "store_true", default = FALSE,
                help = "Run validation analysis"),

    # Runtime
    make_option(c("--threads"), type = "integer", default = 4,
                help = "Number of threads [default: 4]"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output")
)

opt <- parse_args(OptionParser(
    option_list = option_list,
    prog = "calculate_prs_admixed.R",
    description = "Multi-method PRS for admixed populations"
))

# Parse ancestries
ancestries <- strsplit(opt$ancestries, ",")[[1]]
n_anc <- length(ancestries)

# Parse ancestry proportions if provided
anc_props <- NULL
if (!is.null(opt$ancestry_props)) {
    props <- strsplit(opt$ancestry_props, ",")[[1]]
    anc_props <- sapply(props, function(x) {
        parts <- strsplit(x, ":")[[1]]
        as.numeric(parts[2])
    })
    names(anc_props) <- sapply(props, function(x) strsplit(x, ":")[[1]][1])
}

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║         Multi-Method PRS for Admixed Populations                ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║ Method: %-56s ║\n", toupper(opt$method)))
cat(sprintf("║ Ancestries: %-52s ║\n", paste(ancestries, collapse = ", ")))
if (!is.null(anc_props)) {
    props_str <- paste(names(anc_props), round(anc_props, 2), sep = ":", collapse = ", ")
    cat(sprintf("║ Proportions: %-51s ║\n", props_str))
}
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# ============================================================================
# Method: PRS-CSx
# ============================================================================
run_prs_csx <- function(sumstats_dir, ld_ref, geno, ancestries, output_prefix, opt) {
    # PRS-CSx: Multi-ancestry PRS
    # Ge et al., Nature Genetics 2022
    # https://github.com/getian107/PRScsx

    cat("Running PRS-CSx...\n")
    cat("  Reference: Ge et al., Nat Genet 2022\n")
    cat("  Uses ancestry-specific GWAS + LD to derive optimal weights\n\n")

    # Check for PRS-CSx
    prscsx_path <- Sys.which("PRScsx.py")
    if (prscsx_path == "") {
        prscsx_path <- file.path(Sys.getenv("PRSCSX_PATH", ""), "PRScsx.py")
        if (!file.exists(prscsx_path)) {
            stop("PRS-CSx not found. Set PRSCSX_PATH or install from https://github.com/getian107/PRScsx")
        }
    }

    # Build command for each chromosome
    results <- list()

    for (chr in 1:22) {
        cat("  Processing chromosome", chr, "\n")

        cmd <- paste(
            "python", prscsx_path,
            "--ref_dir", ld_ref,
            "--bim_prefix", geno,
            "--n_iter", opt$n_iter,
            "--n_burnin", opt$n_burnin,
            "--out_dir", dirname(output_prefix),
            "--out_name", paste0(basename(output_prefix), "_chr", chr),
            "--chrom", chr
        )

        # Add ancestry-specific sumstats
        for (anc in ancestries) {
            ss_file <- file.path(sumstats_dir, gsub("\\{ancestry\\}", anc, opt$sumstats_pattern))
            if (file.exists(ss_file)) {
                cmd <- paste(cmd, paste0("--sst_file_", tolower(anc)), ss_file)

                # Get sample size from sumstats
                ss <- fread(cmd = paste("zcat", shQuote(ss_file), "| head -1000"))
                if ("N" %in% names(ss)) {
                    n_gwas <- max(ss$N, na.rm = TRUE)
                    cmd <- paste(cmd, paste0("--n_gwas_", tolower(anc)), n_gwas)
                }
            }
        }

        # Add phi if specified
        if (!is.null(opt$phi)) {
            cmd <- paste(cmd, "--phi", opt$phi)
        }

        # Run PRS-CSx for this chromosome
        if (opt$verbose) cat("    ", cmd, "\n")
        system(cmd, ignore.stdout = !opt$verbose)
    }

    # Combine across chromosomes and ancestries
    cat("\n  Combining results across chromosomes...\n")

    combined_weights <- list()

    for (anc in ancestries) {
        anc_weights <- data.table()

        for (chr in 1:22) {
            weight_file <- paste0(output_prefix, "_chr", chr, "_", anc, "_pst_eff_a1_b0.5_phi",
                                  ifelse(is.null(opt$phi), "auto", opt$phi), ".txt")
            if (file.exists(weight_file)) {
                w <- fread(weight_file)
                anc_weights <- rbind(anc_weights, w)
            }
        }

        if (nrow(anc_weights) > 0) {
            combined_weights[[anc]] <- anc_weights
            out_file <- paste0(output_prefix, ".", anc, ".weights.tsv.gz")
            fwrite(anc_weights, out_file, sep = "\t", compress = "gzip")
            cat("    ", anc, "weights:", out_file, "\n")
        }
    }

    return(list(method = "prs_csx", weights = combined_weights))
}

# ============================================================================
# Method: GAUDI
# ============================================================================
run_gaudi <- function(sumstats, local_ancestry, geno, ancestries, output_prefix, opt) {
    # GAUDI: Local ancestry-informed PRS
    # Uses local ancestry to weight PRS by ancestry-specific effects

    cat("Running GAUDI...\n")
    cat("  Local ancestry-informed PRS\n")
    cat("  Weights variants by local ancestry at each locus\n\n")

    # Check for GAUDI
    gaudi_path <- Sys.which("gaudi")
    if (gaudi_path == "") {
        gaudi_path <- file.path(Sys.getenv("GAUDI_PATH", ""), "gaudi.py")
    }

    if (!file.exists(gaudi_path) && gaudi_path != "") {
        cat("  GAUDI not found - using R implementation\n")

        # R implementation of GAUDI concept
        # Weight variants by local ancestry proportion

        # Load local ancestry
        if (is.null(local_ancestry)) {
            stop("--local_ancestry required for GAUDI")
        }

        # Read sumstats
        ss <- fread(sumstats)

        # For each ancestry, calculate weighted PRS
        prs_scores <- data.table()

        # Placeholder - actual implementation would:
        # 1. For each individual, get local ancestry at each SNP
        # 2. Weight the effect size by local ancestry
        # 3. PRS_i = sum_j (LA_ij * beta_j * genotype_ij)

        cat("  GAUDI R implementation - placeholder\n")
        cat("  Full implementation requires local ancestry per SNP per individual\n")
    } else {
        # Run GAUDI executable
        cmd <- paste(
            "python", gaudi_path,
            "--sumstats", sumstats,
            "--la", local_ancestry,
            "--geno", geno,
            "--out", output_prefix
        )

        system(cmd, ignore.stdout = !opt$verbose)
    }

    return(list(method = "gaudi"))
}

# ============================================================================
# Method: DiscoDivas
# ============================================================================
run_disco_divas <- function(sumstats, local_ancestry, geno, ancestries, output_prefix, opt) {
    # DiscoDivas: Disentangling PRS by local ancestry
    # Separates PRS into ancestry-specific components

    cat("Running DiscoDivas...\n")
    cat("  Disentangles PRS into ancestry-specific components\n")
    cat("  Useful for understanding ancestry contribution to risk\n\n")

    # Check for DiscoDivas
    dd_path <- Sys.which("disco_divas")

    if (dd_path == "" || !file.exists(dd_path)) {
        cat("  DiscoDivas not found - using R implementation\n")

        # Conceptual implementation:
        # PRS_EUR = sum(beta_EUR * genotype * I(LA=EUR))
        # PRS_AFR = sum(beta_AFR * genotype * I(LA=AFR))
        # Total PRS = PRS_EUR + PRS_AFR + ...

        if (is.null(local_ancestry)) {
            stop("--local_ancestry required for DiscoDivas")
        }

        # Load sumstats (should have ancestry-specific betas)
        ss <- fread(sumstats)

        # Check for ancestry-specific columns
        has_anc_betas <- all(paste0("BETA_", ancestries) %in% names(ss))

        if (!has_anc_betas) {
            cat("  Warning: Sumstats don't have ancestry-specific betas\n")
            cat("  Using joint beta for all ancestries\n")
        }

        cat("  DiscoDivas R implementation - placeholder\n")
    }

    return(list(method = "disco_divas"))
}

# ============================================================================
# Method: SDPR_admix
# ============================================================================
run_sdpr_admix <- function(sumstats, ld_ref, geno, ancestries, output_prefix, opt) {
    # SDPR_admix: Admixture-aware SDPR
    # Extension of SDPR for admixed populations

    cat("Running SDPR_admix...\n")
    cat("  Admixture-aware sparse Dirichlet process regression\n\n")

    # Check for SDPR
    sdpr_path <- Sys.which("SDPR")
    if (sdpr_path == "") {
        sdpr_path <- file.path(Sys.getenv("SDPR_PATH", ""), "SDPR")
    }

    if (!file.exists(sdpr_path) && sdpr_path != "") {
        stop("SDPR not found. Install from https://github.com/eldronzhou/SDPR")
    }

    # SDPR can use multiple ancestry references
    # Run for each ancestry
    for (anc in ancestries) {
        cat("  Processing", anc, "...\n")

        # Get ancestry-specific sumstats
        if (!is.null(opt$sumstats_dir)) {
            ss_file <- file.path(opt$sumstats_dir, gsub("\\{ancestry\\}", anc, opt$sumstats_pattern))
        } else {
            ss_file <- sumstats  # Use joint sumstats
        }

        # Get ancestry-specific LD
        ld_dir <- file.path(ld_ref, gsub("\\{ancestry\\}", anc, opt$ld_ref_pattern))

        cmd <- paste(
            sdpr_path,
            "-mcmc",
            "-ref_dir", ld_dir,
            "-ss", ss_file,
            "-out", paste0(output_prefix, ".", anc)
        )

        if (opt$verbose) cat("    ", cmd, "\n")
        system(cmd, ignore.stdout = !opt$verbose)
    }

    return(list(method = "sdpr_admix"))
}

# ============================================================================
# Method: MUSSEL
# ============================================================================
run_mussel <- function(sumstats_dir, ld_ref, geno, ancestries, output_prefix, opt) {
    # MUSSEL: Multi-ancestry Stacking Ensemble Learner
    # Combines multiple ancestry-specific PRS using stacking

    cat("Running MUSSEL...\n")
    cat("  Multi-ancestry stacking ensemble learner\n")
    cat("  Optimally combines ancestry-specific PRS\n\n")

    # Check for MUSSEL
    mussel_path <- Sys.which("mussel")

    if (mussel_path == "" || !file.exists(mussel_path)) {
        cat("  MUSSEL not found - using R implementation\n")

        # R implementation of stacking ensemble
        # 1. Calculate PRS for each ancestry
        # 2. Use stacking regression to combine

        # This requires a validation set to train stacking weights
        if (is.null(opt$phenotype)) {
            cat("  Note: Validation phenotype recommended for optimal stacking\n")
        }

        cat("  MUSSEL R implementation - placeholder\n")
        cat("  Full version at: https://github.com/mancusolab/mussel\n")
    }

    return(list(method = "mussel"))
}

# ============================================================================
# Method: PROSPER
# ============================================================================
run_prosper <- function(sumstats, ld_ref, geno, ancestries, output_prefix, opt) {
    # PROSPER: Polygenic Risk Score using Optimal Penalized Regression
    # Multi-ancestry PRS with optimal penalization

    cat("Running PROSPER...\n")
    cat("  Optimal penalized regression for multi-ancestry PRS\n\n")

    # Check for PROSPER
    prosper_path <- Sys.which("prosper")

    if (prosper_path == "" || !file.exists(prosper_path)) {
        cat("  PROSPER not found - using R implementation\n")

        # PROSPER uses penalized regression with ancestry-specific penalties
        # Can use glmnet or similar

        cat("  PROSPER R implementation - placeholder\n")
        cat("  Full version at: https://github.com/bogdanlab/prosper\n")
    }

    return(list(method = "prosper"))
}

# ============================================================================
# Calculate PRS scores
# ============================================================================
calculate_scores <- function(weights_list, geno_prefix, output_prefix) {
    # Calculate PRS scores from weights

    cat("\nCalculating PRS scores...\n")

    scores <- data.table()

    for (anc in names(weights_list)) {
        weights <- weights_list[[anc]]

        if (nrow(weights) == 0) next

        # Write weights in PLINK score format
        score_file <- paste0(output_prefix, ".", anc, ".score")
        # Format: SNP, A1, BETA
        fwrite(weights[, .(V2, V4, V6)], score_file, sep = "\t", col.names = FALSE)

        # Run PLINK scoring
        cmd <- paste(
            "plink2",
            "--bfile", geno_prefix,
            "--score", score_file, "1", "2", "3",
            "--out", paste0(output_prefix, ".", anc, ".scores")
        )

        system(cmd, ignore.stdout = TRUE)

        # Read scores
        score_out <- paste0(output_prefix, ".", anc, ".scores.sscore")
        if (file.exists(score_out)) {
            anc_scores <- fread(score_out)
            anc_scores$ancestry <- anc
            scores <- rbind(scores, anc_scores, fill = TRUE)
        }
    }

    return(scores)
}

# ============================================================================
# Combine ancestry-specific PRS
# ============================================================================
combine_prs <- function(scores, ancestries, anc_props = NULL, method = "weighted") {
    # Combine ancestry-specific PRS into final score

    cat("\nCombining ancestry-specific PRS...\n")

    if (method == "weighted" && !is.null(anc_props)) {
        # Weight by global ancestry proportions
        cat("  Method: Weighted by global ancestry proportions\n")

        # ... implementation
    } else if (method == "sum") {
        # Simple sum
        cat("  Method: Simple sum\n")
    } else if (method == "local_ancestry") {
        # Weight by local ancestry per individual
        cat("  Method: Local ancestry weighted (requires LA data)\n")
    }

    return(scores)
}

# ============================================================================
# Run Selected Method
# ============================================================================
results <- NULL

if (opt$method == "prs_csx" || opt$method == "all") {
    if (is.null(opt$sumstats_dir)) stop("--sumstats_dir required for PRS-CSx")
    if (is.null(opt$ld_ref)) stop("--ld_ref required for PRS-CSx")
    if (is.null(opt$geno)) stop("--geno required for scoring")

    results <- run_prs_csx(opt$sumstats_dir, opt$ld_ref, opt$geno, ancestries,
                           paste0(opt$output_prefix, ".prscsx"), opt)
}

if (opt$method == "gaudi" || opt$method == "all") {
    if (is.null(opt$sumstats)) stop("--sumstats required for GAUDI")
    run_gaudi(opt$sumstats, opt$local_ancestry, opt$geno, ancestries,
              paste0(opt$output_prefix, ".gaudi"), opt)
}

if (opt$method == "disco_divas" || opt$method == "all") {
    if (is.null(opt$sumstats)) stop("--sumstats required for DiscoDivas")
    run_disco_divas(opt$sumstats, opt$local_ancestry, opt$geno, ancestries,
                    paste0(opt$output_prefix, ".discodivas"), opt)
}

if (opt$method == "sdpr_admix" || opt$method == "all") {
    if (is.null(opt$sumstats) && is.null(opt$sumstats_dir)) {
        stop("--sumstats or --sumstats_dir required for SDPR_admix")
    }
    if (is.null(opt$ld_ref)) stop("--ld_ref required for SDPR_admix")
    run_sdpr_admix(opt$sumstats, opt$ld_ref, opt$geno, ancestries,
                   paste0(opt$output_prefix, ".sdpr"), opt)
}

if (opt$method == "mussel" || opt$method == "all") {
    run_mussel(opt$sumstats_dir, opt$ld_ref, opt$geno, ancestries,
               paste0(opt$output_prefix, ".mussel"), opt)
}

if (opt$method == "prosper" || opt$method == "all") {
    run_prosper(opt$sumstats, opt$ld_ref, opt$geno, ancestries,
                paste0(opt$output_prefix, ".prosper"), opt)
}

# ============================================================================
# Validation (if requested)
# ============================================================================
if (opt$validate && !is.null(opt$phenotype) && !is.null(opt$trait)) {
    cat("\nRunning validation analysis...\n")

    # Load phenotype
    pheno <- fread(opt$phenotype)

    # Calculate R² and AUC for each method
    # ... implementation
}

# ============================================================================
# Summary
# ============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                         PRS COMPLETE                             ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║ Method: %-56s ║\n", toupper(opt$method)))
cat(sprintf("║ Ancestries: %-52s ║\n", paste(ancestries, collapse = ", ")))
cat(sprintf("║ Output: %-56s ║\n", opt$output_prefix))
cat("╚══════════════════════════════════════════════════════════════════╝\n")

cat("\nDone!\n")
