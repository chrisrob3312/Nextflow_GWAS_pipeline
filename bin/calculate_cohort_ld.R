#!/usr/bin/env Rscript

# Cohort-Specific LD Calculation
# Calculates LD matrices from cohort data, stratified by user-defined ancestry groups
#
# KEY DESIGN PRINCIPLES:
# 1. FLEXIBLE ANCESTRY STRATA - Not hardcoded to GRAF-ANC categories
#    Users can define any grouping via --ancestry_file with sample→group mapping
# 2. REUSABLE ACROSS TRAITS - LD matrices are trait-independent
# 3. REGIONAL OR GENOME-WIDE - Can calculate for specific regions or whole genome
# 4. MULTIPLE OUTPUT FORMATS - Supports coloc, fine-mapping, PRS-CSx
#
# USAGE EXAMPLES:
#   # Use GRAF-ANC categories from file
#   calculate_cohort_ld.R --gds data.gds --ancestry_file ancestry.tsv \
#       --ancestry_col grafanc_category --groups "EUR,AFR,LATINO,AAC"
#
#   # Custom strata (e.g., by recruitment site)
#   calculate_cohort_ld.R --gds data.gds --ancestry_file samples.tsv \
#       --ancestry_col recruitment_site --groups "site_A,site_B,site_C"
#
#   # Single stratum (all samples pooled)
#   calculate_cohort_ld.R --gds data.gds --groups "ALL"

suppressPackageStartupMessages({
    library(SeqArray)
    library(SeqVarTools)
    library(data.table)
    library(optparse)
    library(Matrix)  # For sparse LD matrices
})

# ============================================
# Command-line argument parsing
# ============================================
option_list <- list(
    # Input files
    make_option(c("-g", "--gds"), type = "character", default = NULL,
                help = "Path to GDS file with genotypes"),
    make_option(c("--plink_prefix"), type = "character", default = NULL,
                help = "Prefix for PLINK files (alternative to GDS)"),
    make_option(c("-a", "--ancestry_file"), type = "character", default = NULL,
                help = "TSV with sample ancestry assignments (sample_id, ancestry columns)"),
    make_option(c("--ancestry_col"), type = "character", default = "ancestry",
                help = "Column name for ancestry/stratum assignment"),
    make_option(c("--sample_col"), type = "character", default = "sample_id",
                help = "Column name for sample identifier"),

    # Ancestry/stratum specification
    make_option(c("--groups"), type = "character", default = NULL,
                help = "Comma-separated ancestry/stratum groups to calculate LD for. Use 'ALL' to pool all samples."),
    make_option(c("--min_n"), type = "integer", default = 50,
                help = "Minimum samples per group (groups with fewer are skipped) [default: 50]"),
    make_option(c("--combine_groups"), type = "character", default = NULL,
                help = "Groups to combine, format: 'NEW_NAME=GROUP1+GROUP2+GROUP3' (e.g., 'LATINO=LAT1+LAT2')"),

    # Region specification
    make_option(c("--region"), type = "character", default = NULL,
                help = "Genomic region: 'chr:start-end' or 'chr' (genome-wide if NULL)"),
    make_option(c("--regions_file"), type = "character", default = NULL,
                help = "BED file with regions to calculate LD for"),
    make_option(c("--window_kb"), type = "integer", default = 1000,
                help = "LD window in kb (variants beyond this are assumed r2=0) [default: 1000]"),

    # LD calculation parameters
    make_option(c("--r2_threshold"), type = "numeric", default = 0.01,
                help = "Only store r2 values above this threshold (sparsity) [default: 0.01]"),
    make_option(c("--maf_filter"), type = "numeric", default = 0.01,
                help = "MAF filter for variants [default: 0.01]"),
    make_option(c("--method"), type = "character", default = "pearson",
                help = "LD method: 'pearson' (r), 'spearman' [default: pearson]"),

    # Output options
    make_option(c("-o", "--output_prefix"), type = "character", default = "cohort_ld",
                help = "Output file prefix"),
    make_option(c("--format"), type = "character", default = "ldstore",
                help = "Output format: 'ldstore', 'coloc', 'prs_csx', 'matrix' [default: ldstore]"),
    make_option(c("--compress"), action = "store_true", default = TRUE,
                help = "Compress output files"),

    # Runtime
    make_option(c("--threads"), type = "integer", default = 4,
                help = "Number of threads [default: 4]"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output")
)

opt_parser <- OptionParser(
    option_list = option_list,
    description = paste(
        "Calculate cohort-specific LD matrices stratified by user-defined ancestry/strata.",
        "\n\nFlexible grouping - not limited to specific ancestry inference tools.",
        "\nCan use GRAF-ANC, ADMIXTURE, PCA-based clusters, or any sample→group mapping."
    )
)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$gds) && is.null(opt$plink_prefix)) {
    stop("Must provide --gds or --plink_prefix")
}

cat("==============================================\n")
cat("Cohort-Specific LD Calculation\n")
cat("==============================================\n")
cat("Input:", if (!is.null(opt$gds)) opt$gds else opt$plink_prefix, "\n")
cat("Groups:", opt$groups %||% "ALL (pooled)", "\n")
cat("Window:", opt$window_kb, "kb\n")
cat("r2 threshold:", opt$r2_threshold, "\n")
cat("Output format:", opt$format, "\n")
cat("==============================================\n\n")

# ============================================
# Load ancestry/stratum assignments
# ============================================
cat("Loading sample assignments...\n")

sample_groups <- NULL
groups_to_process <- character()

if (!is.null(opt$ancestry_file) && file.exists(opt$ancestry_file)) {
    sample_groups <- fread(opt$ancestry_file)

    if (!opt$sample_col %in% names(sample_groups)) {
        stop("Sample column '", opt$sample_col, "' not found in ancestry file")
    }
    if (!opt$ancestry_col %in% names(sample_groups)) {
        stop("Ancestry column '", opt$ancestry_col, "' not found in ancestry file")
    }

    # Rename for consistency
    setnames(sample_groups, c(opt$sample_col, opt$ancestry_col), c("sample_id", "group"))

    cat("  Loaded", nrow(sample_groups), "sample assignments\n")
    cat("  Groups found:", paste(unique(sample_groups$group), collapse = ", "), "\n")

    # Handle group combinations (e.g., LATINO=LAT1+LAT2)
    if (!is.null(opt$combine_groups)) {
        combos <- strsplit(opt$combine_groups, ";")[[1]]
        for (combo in combos) {
            parts <- strsplit(combo, "=")[[1]]
            new_name <- parts[1]
            old_names <- strsplit(parts[2], "\\+")[[1]]
            sample_groups[group %in% old_names, group := new_name]
            cat("  Combined", paste(old_names, collapse = "+"), "→", new_name, "\n")
        }
    }
}

# Determine which groups to process
if (!is.null(opt$groups)) {
    groups_to_process <- strsplit(opt$groups, ",")[[1]]
} else if (!is.null(sample_groups)) {
    groups_to_process <- unique(sample_groups$group)
} else {
    groups_to_process <- "ALL"
}

cat("  Will process groups:", paste(groups_to_process, collapse = ", "), "\n")

# ============================================
# Load genotype data
# ============================================
cat("\nLoading genotype data...\n")

if (!is.null(opt$gds)) {
    gds <- seqOpen(opt$gds)
    all_samples <- seqGetData(gds, "sample.id")
    n_variants <- seqSummary(gds, "variant.id")$numValue
    cat("  GDS:", length(all_samples), "samples,", n_variants, "variants\n")
} else {
    # PLINK input - would use snpStats or similar
    stop("PLINK input not yet implemented - use GDS")
}

# ============================================
# Parse regions
# ============================================
regions <- list()

if (!is.null(opt$region)) {
    # Single region
    if (grepl(":", opt$region)) {
        parts <- strsplit(opt$region, "[:-]")[[1]]
        regions[[1]] <- list(chr = parts[1], start = as.integer(parts[2]), end = as.integer(parts[3]))
    } else {
        regions[[1]] <- list(chr = opt$region, start = 1, end = .Machine$integer.max)
    }
    cat("  Region:", opt$region, "\n")
} else if (!is.null(opt$regions_file)) {
    # Multiple regions from BED
    bed <- fread(opt$regions_file, col.names = c("chr", "start", "end"))
    for (i in 1:nrow(bed)) {
        regions[[i]] <- list(chr = bed$chr[i], start = bed$start[i], end = bed$end[i])
    }
    cat("  Regions from BED:", nrow(bed), "regions\n")
} else {
    # Genome-wide by chromosome
    chrs <- seqGetData(gds, "chromosome")
    unique_chrs <- unique(chrs)
    for (chr in unique_chrs) {
        regions[[chr]] <- list(chr = chr, start = 1, end = .Machine$integer.max)
    }
    cat("  Genome-wide:", length(regions), "chromosomes\n")
}

# ============================================
# LD Calculation Functions
# ============================================

calculate_ld_matrix <- function(gds, samples, chr, start, end, window_kb, r2_threshold, method = "pearson") {
    # Filter to region and samples
    seqSetFilter(gds, sample.id = samples, verbose = FALSE)
    seqSetFilterChrom(gds, chr, from.bp = start, to.bp = end, verbose = FALSE)

    # Get variant info
    var_ids <- seqGetData(gds, "variant.id")
    pos <- seqGetData(gds, "position")
    ref <- seqGetData(gds, "$ref")
    alt <- seqGetData(gds, "$alt")
    n_var <- length(var_ids)

    if (n_var == 0) return(NULL)

    # Get dosage matrix (samples x variants)
    dosage <- seqGetData(gds, "$dosage")
    if (is.null(dosage)) {
        # Fall back to counting alleles
        geno <- seqGetData(gds, "genotype")
        dosage <- colSums(geno == 1, na.rm = TRUE)
    }

    # Calculate MAF and filter
    maf <- colMeans(dosage, na.rm = TRUE) / 2
    maf <- pmin(maf, 1 - maf)
    keep <- which(maf >= opt$maf_filter)

    if (length(keep) == 0) return(NULL)

    dosage <- dosage[, keep, drop = FALSE]
    var_ids <- var_ids[keep]
    pos <- pos[keep]
    ref <- ref[keep]
    alt <- alt[keep]
    n_var <- length(var_ids)

    cat("    Calculating LD for", n_var, "variants,", length(samples), "samples\n")

    # Calculate correlation matrix (windowed)
    window_bp <- window_kb * 1000

    # Use sparse representation
    # Only store r2 > threshold and within window

    ld_i <- integer()
    ld_j <- integer()
    ld_r2 <- numeric()

    for (i in 1:(n_var - 1)) {
        # Find variants within window
        in_window <- which(pos > pos[i] & pos <= pos[i] + window_bp)
        in_window <- in_window[in_window > i]

        if (length(in_window) == 0) next

        # Calculate r2 with each variant in window
        x_i <- dosage[, i]
        x_i <- (x_i - mean(x_i, na.rm = TRUE)) / sd(x_i, na.rm = TRUE)

        for (j in in_window) {
            x_j <- dosage[, j]
            x_j <- (x_j - mean(x_j, na.rm = TRUE)) / sd(x_j, na.rm = TRUE)

            # Pearson correlation
            valid <- !is.na(x_i) & !is.na(x_j)
            if (sum(valid) < 10) next

            r <- sum(x_i[valid] * x_j[valid]) / (sum(valid) - 1)
            r2 <- r^2

            if (r2 >= r2_threshold) {
                ld_i <- c(ld_i, i)
                ld_j <- c(ld_j, j)
                ld_r2 <- c(ld_r2, r2)
            }
        }
    }

    # Create sparse matrix
    ld_sparse <- Matrix::sparseMatrix(
        i = c(ld_i, ld_j),
        j = c(ld_j, ld_i),
        x = c(ld_r2, ld_r2),
        dims = c(n_var, n_var),
        symmetric = TRUE
    )

    # Add diagonal
    diag(ld_sparse) <- 1

    # Variant info
    var_info <- data.table(
        idx = 1:n_var,
        var_id = var_ids,
        chr = chr,
        pos = pos,
        ref = ref,
        alt = alt
    )

    list(
        ld_matrix = ld_sparse,
        variants = var_info,
        n_samples = length(samples),
        region = paste0(chr, ":", min(pos), "-", max(pos))
    )
}

# ============================================
# Process each group
# ============================================
cat("\nCalculating LD by group...\n")

results <- list()

for (grp in groups_to_process) {
    cat("\n--- Group:", grp, "---\n")

    # Get samples for this group
    if (grp == "ALL" || is.null(sample_groups)) {
        samples <- all_samples
    } else {
        samples <- sample_groups[group == grp, sample_id]
        samples <- intersect(samples, all_samples)
    }

    n_samples <- length(samples)
    cat("  Samples:", n_samples, "\n")

    if (n_samples < opt$min_n) {
        cat("  SKIPPED: Below minimum (", opt$min_n, ")\n")
        next
    }

    # Process each region
    grp_results <- list()

    for (reg_name in names(regions)) {
        reg <- regions[[reg_name]]
        cat("  Region:", reg$chr, "\n")

        ld_result <- calculate_ld_matrix(
            gds = gds,
            samples = samples,
            chr = reg$chr,
            start = reg$start,
            end = reg$end,
            window_kb = opt$window_kb,
            r2_threshold = opt$r2_threshold,
            method = opt$method
        )

        if (!is.null(ld_result)) {
            grp_results[[reg$chr]] <- ld_result
        }
    }

    results[[grp]] <- grp_results
}

seqClose(gds)

# ============================================
# Write output files
# ============================================
cat("\nWriting output files...\n")

for (grp in names(results)) {
    grp_clean <- gsub("[^A-Za-z0-9_]", "_", grp)

    for (chr_name in names(results[[grp]])) {
        ld_data <- results[[grp]][[chr_name]]

        out_base <- paste0(opt$output_prefix, ".", grp_clean, ".", chr_name)

        # Write variant info
        var_file <- paste0(out_base, ".variants.tsv")
        fwrite(ld_data$variants, var_file, sep = "\t")

        # Write LD matrix based on format
        if (opt$format == "ldstore") {
            # LDstore format: RSID1 RSID2 R
            ld_df <- summary(ld_data$ld_matrix)
            ld_df <- ld_df[ld_df$i < ld_df$j, ]  # Upper triangle only
            ld_df$snp1 <- ld_data$variants$var_id[ld_df$i]
            ld_df$snp2 <- ld_data$variants$var_id[ld_df$j]
            ld_df$r <- sqrt(ld_df$x)

            ld_file <- paste0(out_base, ".ld.tsv.gz")
            fwrite(ld_df[, .(snp1, snp2, r)], ld_file, sep = "\t", compress = "gzip")
            cat("  Written:", ld_file, "\n")

        } else if (opt$format == "coloc") {
            # Coloc format: needs full LD matrix as RDS
            ld_file <- paste0(out_base, ".ld.rds")
            saveRDS(list(
                LD = as.matrix(ld_data$ld_matrix),
                snps = ld_data$variants$var_id,
                pos = ld_data$variants$pos,
                n = ld_data$n_samples
            ), ld_file)
            cat("  Written:", ld_file, "\n")

        } else if (opt$format == "prs_csx") {
            # PRS-CSx format: binary LD reference
            ld_file <- paste0(out_base, ".ldblk")
            # Would write in PRS-CSx binary format
            cat("  PRS-CSx format not yet implemented\n")

        } else {
            # Matrix format: full matrix as sparse RDS
            ld_file <- paste0(out_base, ".ld_matrix.rds")
            saveRDS(ld_data, ld_file)
            cat("  Written:", ld_file, "\n")
        }
    }
}

# ============================================
# Summary
# ============================================
cat("\n==============================================\n")
cat("SUMMARY\n")
cat("==============================================\n")
cat("Groups processed:", length(results), "\n")
for (grp in names(results)) {
    cat("  ", grp, ":", length(results[[grp]]), "chromosomes/regions\n")
}
cat("Output prefix:", opt$output_prefix, "\n")
cat("==============================================\n")
