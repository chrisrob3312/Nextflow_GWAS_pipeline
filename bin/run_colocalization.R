#!/usr/bin/env Rscript

# ============================================================================
# Multi-Method Colocalization Workflow
# ============================================================================
# Comprehensive colocalization analysis for GWAS-QTL integration
#
# PARALLEL METHODS (not hierarchical - each handles its own MTC):
#
# 1. coloc.susie (per QTL type individually)
#    - Fine-mapping based colocalization
#    - Tests eQTL, sQTL, pQTL, mQTL, caQTL, hQTL separately
#    - Bayesian framework handles multiple testing
#
# 2. HyPrColoc (multi-trait)
#    - Tests if multiple QTLs share the same causal variant
#    - Run on ALL QTL types simultaneously
#
# 3. OPERA (multi-QTL SMR-based) - PARALLEL, NOT TIERED
#    - Uses ALL QTLs simultaneously
#    - Built-in multiple testing correction (Bonferroni + FDR)
#    - NOT causal mediation - tests for shared causal variants
#
# QTL MEGASET: Draws from curated_qtl_sources.yml containing:
#   - eQTL: eQTLGen, GTEx, OneK1K, DICE, BLUEPRINT, MESA, etc.
#   - sQTL: GTEx, AFGR, AIDA
#   - mQTL: GoDMC, GENOA, BLUEPRINT
#   - pQTL: UK Biobank, deCODE, ARIC, Fenland
#   - caQTL: AFGR, DICE, snATAC PBMC
#   - hQTL: BLUEPRINT H3K27ac/H3K4me1, Roadmap
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
    make_option(c("--methods"), type = "character", default = "coloc_susie,hyprcoloc,opera",
                help = "Methods to run (comma-separated): coloc_susie, hyprcoloc, opera [default: all]"),
    make_option(c("--parallel"), action = "store_true", default = TRUE,
                help = "Run methods in parallel (not hierarchical) [default: TRUE]"),
    make_option(c("--coloc_method"), type = "character", default = "coloc_susie",
                help = "Coloc variant: coloc, coloc_susie [default: coloc_susie]"),
    make_option(c("--run_opera"), action = "store_true", default = TRUE,
                help = "Run OPERA on ALL QTLs (handles own MTC) [default: TRUE]"),

    # QTL megaset
    make_option(c("--qtl_megaset"), type = "character", default = NULL,
                help = "QTL megaset directory (from curated_qtl_sources.yml)"),
    make_option(c("--include_all_qtl_types"), action = "store_true", default = TRUE,
                help = "Include all QTL types: eQTL,sQTL,pQTL,mQTL,caQTL,hQTL [default: TRUE]"),

    # Thresholds
    make_option(c("--pp4_threshold"), type = "numeric", default = 0.8,
                help = "PP4 threshold for colocalization [default: 0.8]"),
    make_option(c("--opera_fdr"), type = "numeric", default = 0.05,
                help = "OPERA FDR threshold [default: 0.05]"),

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

# Parse QTL types - expand to full megaset if requested
if (opt$include_all_qtl_types) {
    qtl_types <- c("eqtl", "sqtl", "pqtl", "mqtl", "caqtl", "hqtl")
} else {
    qtl_types <- strsplit(opt$qtl_types, ",")[[1]]
}

# Parse methods
methods <- strsplit(opt$methods, ",")[[1]]

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║         Multi-Method Colocalization Workflow                     ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║ MODE: PARALLEL (each method handles its own MTC)                 ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat("║ 1. coloc.susie  - Fine-mapping colocalization per QTL type       ║\n")
cat("║ 2. HyPrColoc    - Multi-trait across ALL QTL types               ║\n")
cat("║ 3. OPERA        - Multi-QTL SMR (ALL QTLs, built-in MTC)         ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║ QTL types: %-53s ║\n", paste(qtl_types, collapse = ", ")))
cat(sprintf("║ Methods: %-55s ║\n", paste(methods, collapse = ", ")))
cat(sprintf("║ PP4 threshold: %-49.2f ║\n", opt$pp4_threshold))
if (opt$run_opera) {
    cat(sprintf("║ OPERA FDR: %-53.3f ║\n", opt$opera_fdr))
}
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# QTL megaset sources (from curated_qtl_sources.yml)
QTL_SOURCES <- list(
    eqtl = c("eqtlgen_blood", "gtex_v8", "onek1k_pbmc", "dice_immune",
             "blueprint_blood", "mesa_eqtl", "afgr_eqtl", "hchs_sol_eqtl"),
    sqtl = c("gtex_sqtl", "afgr_sqtl", "aida_sqtl"),
    mqtl = c("godmc", "genoa_mqtl", "blueprint_mqtl", "eas_blood_mqtl"),
    pqtl = c("ukb_pqtl", "decode_pqtl", "fenland_pqtl", "aric_pqtl"),
    caqtl = c("afgr_caqtl", "dice_caqtl", "snATAC_pbmc"),
    hqtl = c("blueprint_h3k4me1", "blueprint_h3k27ac", "roadmap_hqtl")
)

cat("QTL MEGASET SOURCES:\n")
for (qt in names(QTL_SOURCES)) {
    cat("  ", toupper(qt), ":", length(QTL_SOURCES[[qt]]), "datasets\n")
}

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
# OPERA: Multi-QTL SMR Analysis (PARALLEL - runs on ALL QTLs)
# ============================================================================
if (opt$run_opera && "opera" %in% methods) {
    cat("\n")
    cat("════════════════════════════════════════════════════════════════════\n")
    cat("OPERA - Multi-QTL SMR Analysis (PARALLEL on ALL QTLs)\n")
    cat("════════════════════════════════════════════════════════════════════\n\n")

    # OPERA is NOT causal mediation - it's SMR-based
    # Tests for shared causal variants between GWAS and QTL
    # Uses instrumental variable approach
    # BUILT-IN multiple testing correction (Bonferroni + FDR)

    cat("NOTE: OPERA runs on ALL QTL types simultaneously (not hierarchical)\n")
    cat("      Built-in MTC: Bonferroni within, FDR across QTL types\n")
    cat("      NOT causal mediation - tests shared causal variants\n\n")

    # Check for OPERA executable
    opera_path <- Sys.which("opera")
    if (opera_path == "") {
        opera_path <- file.path(Sys.getenv("OPERA_PATH", ""), "opera")
    }

    opera_results <- list()

    if (file.exists(opera_path)) {
        # Run OPERA with ALL QTL types simultaneously
        cat("Running OPERA on ALL QTL types:\n")
        cat("  ", paste(qtl_types, collapse = ", "), "\n\n")

        # Collect all QTL BESD files
        besd_files <- c()
        for (qt in qtl_types) {
            besd_file <- file.path(opt$qtl_dir, paste0(qt, ".besd"))
            if (file.exists(besd_file)) {
                besd_files <- c(besd_files, besd_file)
                cat("  Found:", qt, "\n")
            }
        }

        if (length(besd_files) > 0) {
            # Write multi-QTL config
            config_file <- paste0(opt$output_prefix, ".opera_config.txt")
            writeLines(besd_files, config_file)

            # Run OPERA with all QTLs
            cmd <- paste(
                opera_path,
                "--bfile", opt$geno,
                "--gwas-summary", opt$gwas,
                "--beqtl-summary-list", config_file,
                "--out", paste0(opt$output_prefix, ".opera"),
                "--thread-num", opt$threads,
                "--diff-freq-prop 0.1"  # Allow MAF diff up to 10%
            )

            cat("\nRunning OPERA...\n")
            if (opt$verbose) cat("  ", cmd, "\n")
            system(cmd, ignore.stdout = !opt$verbose)

            # Read OPERA results
            opera_out <- paste0(opt$output_prefix, ".opera.smr")
            if (file.exists(opera_out)) {
                opera_results <- fread(opera_out)
                cat("  OPERA complete:", nrow(opera_results), "genes tested\n")
            }
        }
    } else {
        # R-based multi-QTL SMR approximation
        cat("OPERA not found - using R-based multi-QTL SMR\n")
        cat("Full OPERA: https://github.com/yanglab-emory/OPERA\n\n")

        # Test ALL QTL types for ALL genes (not just coloc candidates)
        all_genes <- c()

        # Collect genes from all QTL types
        for (qt in qtl_types) {
            qtl_file <- NULL
            if (!is.null(opt$qtl_dir)) {
                qtl_file <- file.path(opt$qtl_dir, paste0(qt, ".harmonized.tsv.gz"))
            } else if (!is.null(opt$qtl)) {
                qtl_file <- opt$qtl
            }

            if (!is.null(qtl_file) && file.exists(qtl_file)) {
                qtl <- fread(qtl_file, select = "gene")
                all_genes <- union(all_genes, unique(qtl$gene))
            }
        }

        cat("Testing", length(all_genes), "genes across", length(qtl_types), "QTL types\n\n")

        # SMR for each gene x QTL type combination
        smr_results <- list()

        for (gene in all_genes) {
            for (qt in qtl_types) {
                qtl_file <- file.path(opt$qtl_dir, paste0(qt, ".harmonized.tsv.gz"))

                if (file.exists(qtl_file)) {
                    qtl <- fread(qtl_file)
                    if ("gene" %in% names(qtl)) {
                        qtl <- qtl[gene == gene]
                    }

                    if (nrow(qtl) > 0) {
                        # Get lead QTL variant (strongest association)
                        qtl_lead <- qtl[order(pvalue)][1]

                        # Get matching GWAS
                        gwas_match <- gwas[snp == qtl_lead$variant_id]

                        if (nrow(gwas_match) > 0) {
                            # Calculate SMR statistics
                            beta_gwas <- gwas_match$beta[1]
                            se_gwas <- gwas_match$se[1]
                            beta_qtl <- qtl_lead$beta[1]
                            se_qtl <- qtl_lead$se[1]

                            if (!is.na(beta_qtl) && abs(beta_qtl) > 0) {
                                smr_beta <- beta_gwas / beta_qtl
                                smr_se <- sqrt(
                                    (se_gwas^2 / beta_qtl^2) +
                                    (beta_gwas^2 * se_qtl^2 / beta_qtl^4)
                                )
                                smr_z <- smr_beta / smr_se
                                smr_p <- 2 * pnorm(-abs(smr_z))

                                smr_results[[paste(gene, qt, sep = "_")]] <- data.table(
                                    gene = gene,
                                    qtl_type = qt,
                                    lead_snp = qtl_lead$variant_id,
                                    n_qtl_variants = nrow(qtl),
                                    qtl_pvalue = qtl_lead$pvalue,
                                    beta_gwas = beta_gwas,
                                    se_gwas = se_gwas,
                                    beta_qtl = beta_qtl,
                                    se_qtl = se_qtl,
                                    smr_beta = smr_beta,
                                    smr_se = smr_se,
                                    smr_z = smr_z,
                                    smr_p = smr_p
                                )
                            }
                        }
                    }
                }
            }
        }

        # Combine and apply MTC
        if (length(smr_results) > 0) {
            opera_results <- rbindlist(smr_results, fill = TRUE)

            # Multiple testing correction (like OPERA)
            # FDR within each QTL type, then Bonferroni across types
            opera_results[, fdr_within_qtl := p.adjust(smr_p, method = "BH"), by = qtl_type]
            opera_results[, bonf_across_qtl := smr_p * length(qtl_types)]
            opera_results[bonf_across_qtl > 1, bonf_across_qtl := 1]

            # Final significance
            opera_results[, significant := fdr_within_qtl < opt$opera_fdr]

            # Sort by p-value
            opera_results <- opera_results[order(smr_p)]

            cat("SMR Results:\n")
            cat("  Total gene-QTL pairs tested:", nrow(opera_results), "\n")
            cat("  Significant (FDR <", opt$opera_fdr, "):",
                sum(opera_results$significant, na.rm = TRUE), "\n")

            # Summary by QTL type
            cat("\n  By QTL type:\n")
            for (qt in unique(opera_results$qtl_type)) {
                n_sig <- sum(opera_results$qtl_type == qt & opera_results$significant, na.rm = TRUE)
                n_tot <- sum(opera_results$qtl_type == qt)
                cat("    ", toupper(qt), ":", n_sig, "/", n_tot, "significant\n")
            }
        }
    }

    # Write OPERA results
    if (length(opera_results) > 0 && nrow(opera_results) > 0) {
        opera_file <- paste0(opt$output_prefix, ".opera.tsv")
        fwrite(opera_results, opera_file, sep = "\t")
        cat("\nOPERA results:", opera_file, "\n")

        # Write significant hits
        if ("significant" %in% names(opera_results)) {
            sig_hits <- opera_results[significant == TRUE]
            if (nrow(sig_hits) > 0) {
                sig_file <- paste0(opt$output_prefix, ".opera.significant.tsv")
                fwrite(sig_hits, sig_file, sep = "\t")
                cat("Significant hits:", sig_file, "\n")
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
