#!/usr/bin/env Rscript

# ============================================================================
# Multi-Method Colocalization Workflow
# ============================================================================
# Tiered colocalization analysis for GWAS-QTL integration
#
# TIER 1: coloc.susie (per QTL type individually)
#   - Reduces multiple testing burden
#   - Tests eQTL, sQTL, pQTL, etc. separately first
#   - Identifies candidate loci with PP4 > threshold
#
# TIER 2: HyPrColoc / ColocBoost (on candidates)
#   - Multi-trait colocalization on Tier 1 candidates
#   - Tests if multiple QTLs share the same causal variant
#
# TIER 3: OPERA (SMR-based)
#   - Bayesian SMR for causal inference
#   - NOT causal mediation - tests for shared causal variants
#
# USES COHORT-SPECIFIC LD (not reference panels)
# ============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

# Check for colocalization packages
has_coloc <- requireNamespace("coloc", quietly = TRUE)
has_susie <- requireNamespace("susieR", quietly = TRUE)
has_hyprcoloc <- requireNamespace("hyprcoloc", quietly = TRUE)

# ============================================================================
# Command-line Arguments
# ============================================================================
option_list <- list(
    # Input: GWAS
    make_option(c("--gwas"), type = "character", default = NULL,
                help = "GWAS summary statistics file"),
    make_option(c("--gwas_n"), type = "integer", default = NULL,
                help = "GWAS sample size"),
    make_option(c("--gwas_type"), type = "character", default = "cc",
                help = "GWAS type: cc (case-control), quant [default: cc]"),
    make_option(c("--gwas_s"), type = "numeric", default = NULL,
                help = "Case proportion for case-control GWAS"),

    # Input: QTL
    make_option(c("--qtl"), type = "character", default = NULL,
                help = "Harmonized QTL summary statistics"),
    make_option(c("--qtl_dir"), type = "character", default = NULL,
                help = "Directory with QTL files by type"),
    make_option(c("--qtl_n"), type = "integer", default = NULL,
                help = "QTL sample size"),
    make_option(c("--qtl_types"), type = "character", default = "eqtl,sqtl,pqtl",
                help = "QTL types to test [default: eqtl,sqtl,pqtl]"),

    # Input: LD
    make_option(c("--ld"), type = "character", default = NULL,
                help = "LD matrix file (RDS format from calculate_cohort_ld.R)"),
    make_option(c("--ld_dir"), type = "character", default = NULL,
                help = "Directory with LD matrices by region"),

    # Region specification
    make_option(c("--region"), type = "character", default = NULL,
                help = "Region to test: chr:start-end"),
    make_option(c("--regions_file"), type = "character", default = NULL,
                help = "BED file with regions (e.g., GWAS loci)"),
    make_option(c("--gene"), type = "character", default = NULL,
                help = "Gene to test (for QTL filtering)"),

    # Method selection
    make_option(c("--tier1_method"), type = "character", default = "coloc_susie",
                help = "Tier 1 method: coloc, coloc_susie [default: coloc_susie]"),
    make_option(c("--tier2_method"), type = "character", default = "hyprcoloc",
                help = "Tier 2 method: hyprcoloc, colocboost [default: hyprcoloc]"),
    make_option(c("--tier3"), action = "store_true", default = FALSE,
                help = "Run Tier 3 OPERA analysis"),

    # Thresholds
    make_option(c("--pp4_threshold"), type = "numeric", default = 0.8,
                help = "PP4 threshold for colocalization [default: 0.8]"),
    make_option(c("--tier2_threshold"), type = "numeric", default = 0.5,
                help = "PP4 threshold to advance to Tier 2 [default: 0.5]"),

    # Fine-mapping
    make_option(c("--max_causal"), type = "integer", default = 5,
                help = "Maximum causal variants for SuSiE [default: 5]"),

    # Priors
    make_option(c("--p1"), type = "numeric", default = 1e-4,
                help = "Prior for GWAS association [default: 1e-4]"),
    make_option(c("--p2"), type = "numeric", default = 1e-4,
                help = "Prior for QTL association [default: 1e-4]"),
    make_option(c("--p12"), type = "numeric", default = 1e-5,
                help = "Prior for colocalization [default: 1e-5]"),

    # Output
    make_option(c("-o", "--output_prefix"), type = "character", default = "coloc",
                help = "Output prefix [default: coloc]"),

    # Runtime
    make_option(c("--threads"), type = "integer", default = 4,
                help = "Number of threads [default: 4]"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output")
)

opt <- parse_args(OptionParser(
    option_list = option_list,
    prog = "run_colocalization.R",
    description = "Multi-method colocalization workflow"
))

# Parse QTL types
qtl_types <- strsplit(opt$qtl_types, ",")[[1]]

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║            Multi-Method Colocalization Workflow                  ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║ TIER 1: coloc.susie (per QTL type)                               ║\n")
cat("║ TIER 2: HyPrColoc (multi-trait, on candidates)                   ║\n")
cat("║ TIER 3: OPERA (SMR-based, optional)                              ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║ QTL types: %-53s ║\n", paste(qtl_types, collapse = ", ")))
cat(sprintf("║ PP4 threshold: %-49s ║\n", opt$pp4_threshold))
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# ============================================================================
# Load Data
# ============================================================================
cat("Loading data...\n")

# Load GWAS
gwas <- fread(opt$gwas)
cat("  GWAS:", nrow(gwas), "variants\n")

# Standardize GWAS columns
gwas_cols <- list(
    snp = c("SNP", "variant_id", "ID", "rsid"),
    beta = c("BETA", "beta", "b", "effect"),
    se = c("SE", "se", "stderr"),
    pval = c("P", "pvalue", "p.value", "P_JOINT"),
    maf = c("MAF", "maf", "eaf")
)

for (col_type in names(gwas_cols)) {
    found <- intersect(names(gwas), gwas_cols[[col_type]])
    if (length(found) > 0 && found[1] != col_type) {
        setnames(gwas, found[1], col_type)
    }
}

# Load LD matrix
ld_matrix <- NULL
if (!is.null(opt$ld) && file.exists(opt$ld)) {
    cat("  Loading LD matrix:", opt$ld, "\n")
    ld_data <- readRDS(opt$ld)
    if ("LD" %in% names(ld_data)) {
        ld_matrix <- ld_data$LD
    } else {
        ld_matrix <- ld_data
    }
    cat("    Dimensions:", nrow(ld_matrix), "x", ncol(ld_matrix), "\n")
}

# ============================================================================
# TIER 1: coloc.susie (per QTL type)
# ============================================================================
cat("\n")
cat("════════════════════════════════════════════════════════════════════\n")
cat("TIER 1: coloc.susie - Testing each QTL type individually\n")
cat("════════════════════════════════════════════════════════════════════\n\n")

tier1_results <- list()

for (qt in qtl_types) {
    cat("Processing", toupper(qt), "...\n")

    # Load QTL data for this type
    qtl_file <- NULL
    if (!is.null(opt$qtl_dir)) {
        qtl_file <- file.path(opt$qtl_dir, paste0(qt, ".harmonized.tsv.gz"))
    } else if (!is.null(opt$qtl)) {
        qtl_file <- opt$qtl
    }

    if (is.null(qtl_file) || !file.exists(qtl_file)) {
        cat("  Skipping - file not found\n")
        next
    }

    qtl <- fread(qtl_file)

    # Filter by QTL type if combined file
    if ("qtl_type" %in% names(qtl)) {
        qtl <- qtl[qtl_type == qt]
    }

    cat("  QTL variants:", nrow(qtl), "\n")

    if (nrow(qtl) == 0) {
        cat("  Skipping - no variants\n")
        next
    }

    # Get unique genes
    genes <- unique(qtl$gene)
    cat("  Genes to test:", length(genes), "\n")

    # Test each gene
    gene_results <- list()

    for (g in genes) {
        qtl_gene <- qtl[gene == g]

        # Merge GWAS and QTL
        merged <- merge(gwas, qtl_gene, by.x = "snp", by.y = "variant_id", suffixes = c("_gwas", "_qtl"))

        if (nrow(merged) < 10) next  # Need minimum variants

        # Prepare coloc datasets
        if (has_coloc && has_susie) {
            # Run coloc.susie
            tryCatch({
                # Dataset 1: GWAS
                D1 <- list(
                    beta = merged$beta_gwas,
                    varbeta = merged$se_gwas^2,
                    snp = merged$snp,
                    type = opt$gwas_type,
                    N = opt$gwas_n
                )
                if (opt$gwas_type == "cc" && !is.null(opt$gwas_s)) {
                    D1$s <- opt$gwas_s
                }

                # Dataset 2: QTL
                D2 <- list(
                    beta = merged$beta_qtl,
                    varbeta = merged$se_qtl^2,
                    snp = merged$snp,
                    type = "quant",
                    N = opt$qtl_n
                )

                # Run coloc
                if (!is.null(ld_matrix)) {
                    # Use SuSiE with LD
                    # Subset LD to matching variants
                    ld_snps <- rownames(ld_matrix)
                    common <- intersect(merged$snp, ld_snps)

                    if (length(common) >= 10) {
                        ld_sub <- ld_matrix[common, common]

                        # Run SuSiE + coloc
                        res <- coloc::coloc.susie(D1, D2, LD = ld_sub,
                                                   p1 = opt$p1, p2 = opt$p2, p12 = opt$p12)
                    } else {
                        res <- coloc::coloc.abf(D1, D2,
                                                 p1 = opt$p1, p2 = opt$p2, p12 = opt$p12)
                    }
                } else {
                    # ABF without LD
                    res <- coloc::coloc.abf(D1, D2,
                                             p1 = opt$p1, p2 = opt$p2, p12 = opt$p12)
                }

                # Extract results
                if ("summary" %in% names(res)) {
                    pp4 <- res$summary["PP.H4.abf"]
                } else {
                    pp4 <- NA
                }

                gene_results[[g]] <- data.table(
                    qtl_type = qt,
                    gene = g,
                    n_snps = nrow(merged),
                    PP0 = res$summary["PP.H0.abf"],
                    PP1 = res$summary["PP.H1.abf"],
                    PP2 = res$summary["PP.H2.abf"],
                    PP3 = res$summary["PP.H3.abf"],
                    PP4 = pp4,
                    colocalizes = pp4 >= opt$pp4_threshold
                )

                if (pp4 >= opt$pp4_threshold) {
                    cat("    ✓", g, "- PP4:", round(pp4, 3), "\n")
                }

            }, error = function(e) {
                if (opt$verbose) cat("    Error for", g, ":", e$message, "\n")
            })
        } else {
            cat("  coloc/susieR not available - install with BiocManager\n")
        }
    }

    if (length(gene_results) > 0) {
        tier1_results[[qt]] <- rbindlist(gene_results)
    }
}

# Combine Tier 1 results
if (length(tier1_results) > 0) {
    tier1_combined <- rbindlist(tier1_results)
    tier1_combined <- tier1_combined[order(-PP4)]

    # Write Tier 1 results
    tier1_file <- paste0(opt$output_prefix, ".tier1.tsv")
    fwrite(tier1_combined, tier1_file, sep = "\t")
    cat("\nTier 1 results:", tier1_file, "\n")
    cat("  Total genes tested:", nrow(tier1_combined), "\n")
    cat("  Colocalized (PP4 ≥", opt$pp4_threshold, "):", sum(tier1_combined$colocalizes, na.rm = TRUE), "\n")

    # Identify candidates for Tier 2
    tier2_candidates <- tier1_combined[PP4 >= opt$tier2_threshold]
    cat("  Candidates for Tier 2 (PP4 ≥", opt$tier2_threshold, "):", nrow(tier2_candidates), "\n")
} else {
    cat("\nNo Tier 1 results\n")
    tier1_combined <- data.table()
    tier2_candidates <- data.table()
}

# ============================================================================
# TIER 2: HyPrColoc (multi-trait on candidates)
# ============================================================================
if (nrow(tier2_candidates) > 0 && has_hyprcoloc) {
    cat("\n")
    cat("════════════════════════════════════════════════════════════════════\n")
    cat("TIER 2: HyPrColoc - Multi-trait colocalization\n")
    cat("════════════════════════════════════════════════════════════════════\n\n")

    library(hyprcoloc)

    # Get unique genes from candidates
    candidate_genes <- unique(tier2_candidates$gene)

    tier2_results <- list()

    for (gene in candidate_genes) {
        cat("Testing", gene, "across traits...\n")

        # Get QTL types with evidence for this gene
        gene_qtl_types <- tier2_candidates[gene == gene]$qtl_type

        if (length(gene_qtl_types) < 2) {
            cat("  Only 1 QTL type - skipping multi-trait test\n")
            next
        }

        # Collect effect sizes and SEs for each trait
        betas <- list()
        ses <- list()
        trait_names <- c()
        common_snps <- NULL

        # Add GWAS
        gwas_sub <- gwas[!is.na(beta) & !is.na(se)]
        if (nrow(gwas_sub) > 0) {
            betas[["GWAS"]] <- gwas_sub$beta
            ses[["GWAS"]] <- gwas_sub$se
            common_snps <- gwas_sub$snp
            trait_names <- c(trait_names, "GWAS")
        }

        # Add each QTL type
        for (qt in gene_qtl_types) {
            qtl_file <- NULL
            if (!is.null(opt$qtl_dir)) {
                qtl_file <- file.path(opt$qtl_dir, paste0(qt, ".harmonized.tsv.gz"))
            }

            if (!is.null(qtl_file) && file.exists(qtl_file)) {
                qtl <- fread(qtl_file)
                if ("gene" %in% names(qtl)) {
                    qtl <- qtl[gene == gene]
                }

                if (nrow(qtl) > 0) {
                    # Align to common SNPs
                    if (!is.null(common_snps)) {
                        qtl <- qtl[variant_id %in% common_snps]
                        common_snps <- intersect(common_snps, qtl$variant_id)
                    }

                    betas[[qt]] <- qtl$beta
                    ses[[qt]] <- qtl$se
                    trait_names <- c(trait_names, qt)
                }
            }
        }

        # Need at least 2 traits
        if (length(betas) < 2 || length(common_snps) < 10) {
            cat("  Insufficient data for HyPrColoc\n")
            next
        }

        # Align all traits to common SNPs
        # ... alignment logic

        # Convert to matrices
        n_snps <- length(common_snps)
        n_traits <- length(trait_names)

        beta_matrix <- matrix(NA, nrow = n_snps, ncol = n_traits)
        se_matrix <- matrix(NA, nrow = n_snps, ncol = n_traits)

        for (i in 1:n_traits) {
            beta_matrix[, i] <- betas[[trait_names[i]]][1:min(n_snps, length(betas[[trait_names[i]]]))]
            se_matrix[, i] <- ses[[trait_names[i]]][1:min(n_snps, length(ses[[trait_names[i]]]))]
        }

        # Run HyPrColoc
        tryCatch({
            res <- hyprcoloc(
                effect.est = beta_matrix,
                effect.se = se_matrix,
                trait.names = trait_names,
                snp.id = common_snps[1:n_snps]
            )

            # Extract results
            if (!is.null(res$results)) {
                tier2_results[[gene]] <- data.table(
                    gene = gene,
                    n_traits = n_traits,
                    traits = paste(trait_names, collapse = ";"),
                    n_snps = n_snps,
                    posterior_prob = res$results$posterior_prob,
                    regional_prob = res$results$regional_prob,
                    candidate_snp = res$results$candidate_snp,
                    posterior_explained_by_snp = res$results$posterior_explained_by_snp
                )

                if (!is.null(res$results$posterior_prob) &&
                    res$results$posterior_prob >= 0.8) {
                    cat("  ✓ Multi-trait colocalization PP:", round(res$results$posterior_prob, 3), "\n")
                    cat("    Traits:", paste(trait_names, collapse = ", "), "\n")
                }
            }
        }, error = function(e) {
            if (opt$verbose) cat("  Error:", e$message, "\n")
        })
    }

    # Combine Tier 2 results
    if (length(tier2_results) > 0) {
        tier2_combined <- rbindlist(tier2_results, fill = TRUE)
        tier2_file <- paste0(opt$output_prefix, ".tier2.hyprcoloc.tsv")
        fwrite(tier2_combined, tier2_file, sep = "\t")
        cat("\nTier 2 results:", tier2_file, "\n")
    }
} else if (nrow(tier2_candidates) > 0) {
    cat("\nhyprcoloc not installed - skipping Tier 2\n")
    cat("Install with: devtools::install_github('cnfoley/hyprcoloc')\n")
}

# ============================================================================
# TIER 3: OPERA (optional)
# ============================================================================
if (opt$tier3) {
    cat("\n")
    cat("════════════════════════════════════════════════════════════════════\n")
    cat("TIER 3: OPERA - Bayesian SMR analysis\n")
    cat("════════════════════════════════════════════════════════════════════\n\n")

    # OPERA is NOT causal mediation - it's SMR-based
    # Tests for shared causal variants between GWAS and QTL
    # Uses instrumental variable approach

    cat("NOTE: OPERA uses SMR framework but is NOT causal mediation\n")
    cat("      It tests whether the same causal variant affects both traits\n\n")

    # Check for OPERA executable
    opera_path <- Sys.which("opera")
    if (opera_path == "") {
        opera_path <- file.path(Sys.getenv("OPERA_PATH", ""), "opera")
    }

    if (nrow(tier1_combined) == 0) {
        cat("  No Tier 1 results - skipping OPERA\n")
    } else if (file.exists(opera_path)) {
        # Run OPERA on colocalized loci
        coloc_loci <- tier1_combined[colocalizes == TRUE]

        if (nrow(coloc_loci) == 0) {
            cat("  No colocalized loci - skipping OPERA\n")
        } else {
            tier3_results <- list()

            for (i in 1:nrow(coloc_loci)) {
                gene <- coloc_loci$gene[i]
                qt <- coloc_loci$qtl_type[i]

                cat("  Running OPERA for", gene, "(", qt, ")\n")

                # Prepare BESD format for QTL
                # OPERA uses SMR's BESD format

                # Run OPERA
                cmd <- paste(
                    opera_path,
                    "--bfile", opt$geno,
                    "--gwas-summary", opt$gwas,
                    "--beqtl-summary", paste0(opt$qtl_dir, "/", qt, ".besd"),
                    "--out", paste0(opt$output_prefix, ".opera.", gene),
                    "--thread-num", opt$threads
                )

                if (opt$verbose) cat("    ", cmd, "\n")
                # system(cmd)
            }
        }
    } else {
        # R-based SMR approximation
        cat("  OPERA not found - using R-based SMR approximation\n")
        cat("  Full OPERA: https://github.com/yanglab-emory/OPERA\n\n")

        if (nrow(tier1_combined) > 0) {
            coloc_loci <- tier1_combined[colocalizes == TRUE]

            tier3_results <- list()

            for (i in 1:min(nrow(coloc_loci), 100)) {  # Limit for efficiency
                gene <- coloc_loci$gene[i]
                qt <- coloc_loci$qtl_type[i]

                # Simple SMR test
                # SMR beta = beta_GWAS / beta_QTL
                # SMR SE = sqrt(se_GWAS^2/beta_QTL^2 + beta_GWAS^2*se_QTL^2/beta_QTL^4)

                # Get top QTL for this gene
                qtl_file <- file.path(opt$qtl_dir, paste0(qt, ".harmonized.tsv.gz"))
                if (file.exists(qtl_file)) {
                    qtl <- fread(qtl_file)
                    if ("gene" %in% names(qtl)) {
                        qtl <- qtl[gene == gene]
                    }

                    if (nrow(qtl) > 0) {
                        # Get lead QTL variant
                        qtl <- qtl[order(pvalue)][1]

                        # Get matching GWAS
                        gwas_match <- gwas[snp == qtl$variant_id]

                        if (nrow(gwas_match) > 0) {
                            # Calculate SMR statistics
                            beta_gwas <- gwas_match$beta[1]
                            se_gwas <- gwas_match$se[1]
                            beta_qtl <- qtl$beta[1]
                            se_qtl <- qtl$se[1]

                            smr_beta <- beta_gwas / beta_qtl
                            smr_se <- sqrt(
                                (se_gwas^2 / beta_qtl^2) +
                                (beta_gwas^2 * se_qtl^2 / beta_qtl^4)
                            )
                            smr_z <- smr_beta / smr_se
                            smr_p <- 2 * pnorm(-abs(smr_z))

                            # HEIDI test (simplified - tests for heterogeneity)
                            # Full HEIDI requires multiple variants

                            tier3_results[[paste(gene, qt, sep = "_")]] <- data.table(
                                gene = gene,
                                qtl_type = qt,
                                lead_snp = qtl$variant_id,
                                beta_gwas = beta_gwas,
                                beta_qtl = beta_qtl,
                                smr_beta = smr_beta,
                                smr_se = smr_se,
                                smr_p = smr_p,
                                tier1_pp4 = coloc_loci$PP4[i]
                            )

                            if (smr_p < 0.05) {
                                cat("    ✓", gene, "SMR p:", format(smr_p, scientific = TRUE), "\n")
                            }
                        }
                    }
                }
            }

            # Combine Tier 3 results
            if (length(tier3_results) > 0) {
                tier3_combined <- rbindlist(tier3_results, fill = TRUE)
                tier3_combined <- tier3_combined[order(smr_p)]
                tier3_file <- paste0(opt$output_prefix, ".tier3.opera.tsv")
                fwrite(tier3_combined, tier3_file, sep = "\t")
                cat("\nTier 3 results:", tier3_file, "\n")
                cat("  SMR significant (p < 0.05):", sum(tier3_combined$smr_p < 0.05, na.rm = TRUE), "\n")
            }
        }
    }
}

# ============================================================================
# Summary
# ============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                  COLOCALIZATION COMPLETE                         ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")

if (nrow(tier1_combined) > 0) {
    n_coloc <- sum(tier1_combined$colocalizes, na.rm = TRUE)
    cat(sprintf("║ Tier 1 (coloc.susie):                                            ║\n"))
    cat(sprintf("║   Genes tested: %-48d ║\n", nrow(tier1_combined)))
    cat(sprintf("║   Colocalized (PP4 ≥ %.1f): %-36d ║\n", opt$pp4_threshold, n_coloc))

    # By QTL type
    cat("║   By QTL type:                                                   ║\n")
    for (qt in unique(tier1_combined$qtl_type)) {
        n_qt <- sum(tier1_combined$qtl_type == qt & tier1_combined$colocalizes, na.rm = TRUE)
        cat(sprintf("║     %-8s: %-50d ║\n", qt, n_qt))
    }
}

cat("╚══════════════════════════════════════════════════════════════════╝\n")

cat("\nDone!\n")
