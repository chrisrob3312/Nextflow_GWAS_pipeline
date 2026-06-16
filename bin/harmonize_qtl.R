#!/usr/bin/env Rscript

# ============================================================================
# QTL Harmonization Module
# ============================================================================
# Harmonizes QTL summary statistics from multiple sources into a unified format
# for downstream colocalization and functional annotation.
#
# FEATURES:
#   - Liftover between genome builds (GRCh37 ↔ GRCh38)
#   - Standardize column names and formats
#   - Deduplicate overlapping variants
#   - Filter by tissue, ancestry, QTL type
#   - Output unified TSV for coloc/HyPrColoc/OPERA
#
# SUPPORTED INPUT FORMATS:
#   - eQTL Catalogue
#   - GTEx v8
#   - eQTLGen
#   - MESA
#   - SuShiE fine-mapping
#   - Custom (with column mapping)
#
# OUTPUT FORMAT (unified):
#   variant_id, chr, pos, ref, alt, gene/feature, beta, se, pvalue,
#   maf, n, tissue, ancestry, qtl_type, source, build
# ============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

# Check for liftover capability
has_rtracklayer <- requireNamespace("rtracklayer", quietly = TRUE)
if (!has_rtracklayer) {
    cat("Note: rtracklayer not available, liftover disabled\n")
}

# ============================================================================
# Command-line Arguments
# ============================================================================
option_list <- list(
    # Input
    make_option(c("-i", "--input"), type = "character", default = NULL,
                help = "Input QTL file or directory"),
    make_option(c("--manifest"), type = "character", default = NULL,
                help = "Manifest TSV listing multiple QTL files with metadata"),
    make_option(c("--format"), type = "character", default = "auto",
                help = "Input format: auto, eqtl_catalogue, gtex, eqtlgen, mesa, sushie, custom"),

    # Column mapping (for custom format)
    make_option(c("--col_variant"), type = "character", default = NULL,
                help = "Variant ID column (for custom format)"),
    make_option(c("--col_chr"), type = "character", default = NULL,
                help = "Chromosome column"),
    make_option(c("--col_pos"), type = "character", default = NULL,
                help = "Position column"),
    make_option(c("--col_ref"), type = "character", default = NULL,
                help = "Reference allele column"),
    make_option(c("--col_alt"), type = "character", default = NULL,
                help = "Alternative allele column"),
    make_option(c("--col_gene"), type = "character", default = NULL,
                help = "Gene/feature column"),
    make_option(c("--col_beta"), type = "character", default = NULL,
                help = "Effect size column"),
    make_option(c("--col_se"), type = "character", default = NULL,
                help = "Standard error column"),
    make_option(c("--col_pvalue"), type = "character", default = NULL,
                help = "P-value column"),

    # Genome build
    make_option(c("--input_build"), type = "character", default = "GRCh38",
                help = "Input genome build [default: GRCh38]"),
    make_option(c("--output_build"), type = "character", default = "GRCh38",
                help = "Output genome build [default: GRCh38]"),
    make_option(c("--chain_file"), type = "character", default = NULL,
                help = "Chain file for liftover"),

    # Filtering
    make_option(c("--tissues"), type = "character", default = NULL,
                help = "Comma-separated tissues to include"),
    make_option(c("--ancestries"), type = "character", default = NULL,
                help = "Comma-separated ancestries to include"),
    make_option(c("--qtl_types"), type = "character", default = NULL,
                help = "Comma-separated QTL types: eqtl,sqtl,pqtl,mqtl,caqtl"),
    make_option(c("--genes"), type = "character", default = NULL,
                help = "File with gene IDs to include (one per line)"),
    make_option(c("--region"), type = "character", default = NULL,
                help = "Genomic region: chr:start-end"),

    # Quality filters
    make_option(c("--pvalue_threshold"), type = "numeric", default = 1,
                help = "P-value threshold [default: 1 (no filter)]"),
    make_option(c("--maf_threshold"), type = "numeric", default = 0,
                help = "MAF threshold [default: 0 (no filter)]"),

    # Output
    make_option(c("-o", "--output"), type = "character", default = "harmonized_qtl",
                help = "Output prefix"),
    make_option(c("--split_by"), type = "character", default = NULL,
                help = "Split output by: tissue, ancestry, qtl_type, gene, chr"),

    # Deduplication
    make_option(c("--dedup_strategy"), type = "character", default = "best_p",
                help = "Deduplication: best_p, first, none [default: best_p]"),

    # Runtime
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output")
)

opt <- parse_args(OptionParser(
    option_list = option_list,
    prog = "harmonize_qtl.R",
    description = "Harmonize QTL summary statistics into unified format"
))

# ============================================================================
# Input Format Specifications
# ============================================================================
# Define column mappings for known formats
FORMAT_SPECS <- list(
    eqtl_catalogue = list(
        variant = "variant",
        chr = "chromosome",
        pos = "position",
        ref = "ref",
        alt = "alt",
        gene = "gene_id",
        beta = "beta",
        se = "se",
        pvalue = "pvalue",
        maf = "maf",
        build = "GRCh38"
    ),
    gtex = list(
        variant = "variant_id",
        chr = "chr",
        pos = "variant_pos",
        ref = "ref",
        alt = "alt",
        gene = "gene_id",
        beta = "slope",
        se = "slope_se",
        pvalue = "pval_nominal",
        maf = "maf",
        build = "GRCh38"
    ),
    eqtlgen = list(
        variant = "SNP",
        chr = "SNPChr",
        pos = "SNPPos",
        ref = "AssessedAllele",  # Note: eQTLGen uses assessed/other
        alt = "OtherAllele",
        gene = "Gene",
        beta = "Zscore",  # Need to convert
        se = NA,
        pvalue = "Pvalue",
        maf = NA,
        build = "GRCh37",  # eQTLGen is GRCh37
        beta_is_z = TRUE
    ),
    mesa = list(
        variant = "rsid",
        chr = "chr",
        pos = "pos",
        ref = "ref",
        alt = "alt",
        gene = "gene",
        beta = "beta",
        se = "se",
        pvalue = "pval",
        maf = "maf",
        build = "GRCh38"
    ),
    sushie = list(
        variant = "snp",
        chr = "chr",
        pos = "bp",
        ref = "a1",
        alt = "a2",
        gene = "gene",
        beta = "beta",
        se = "se",
        pvalue = "pvalue",
        pip = "pip",  # SuShiE includes PIPs
        build = "GRCh38"
    )
)

# ============================================================================
# Functions
# ============================================================================

detect_format <- function(dt) {
    # Auto-detect format from column names
    cols <- tolower(names(dt))

    if ("slope" %in% cols && "variant_id" %in% cols) {
        return("gtex")
    } else if ("snpchr" %in% cols && "assessedallele" %in% cols) {
        return("eqtlgen")
    } else if ("variant" %in% cols && "chromosome" %in% cols && "gene_id" %in% cols) {
        return("eqtl_catalogue")
    } else if ("pip" %in% cols && "snp" %in% cols) {
        return("sushie")
    } else if ("rsid" %in% cols && "gene" %in% cols) {
        return("mesa")
    } else {
        return("custom")
    }
}

read_qtl_file <- function(path, format = "auto", custom_cols = NULL) {
    # Read QTL file and standardize columns
    cat("  Reading:", basename(path), "\n")

    # Handle compressed files
    if (grepl("\\.gz$", path)) {
        dt <- fread(cmd = paste("zcat", shQuote(path)), header = TRUE)
    } else {
        dt <- fread(path, header = TRUE)
    }

    # Detect format if auto
    if (format == "auto") {
        format <- detect_format(dt)
        cat("    Detected format:", format, "\n")
    }

    # Get column mapping
    if (format == "custom") {
        if (is.null(custom_cols)) {
            stop("Custom format requires column mapping")
        }
        col_map <- custom_cols
    } else {
        col_map <- FORMAT_SPECS[[format]]
    }

    # Standardize column names
    setnames_safe <- function(dt, old, new) {
        if (!is.na(old) && old %in% names(dt)) {
            setnames(dt, old, new)
        }
    }

    setnames_safe(dt, col_map$variant, "variant_id")
    setnames_safe(dt, col_map$chr, "chr")
    setnames_safe(dt, col_map$pos, "pos")
    setnames_safe(dt, col_map$ref, "ref")
    setnames_safe(dt, col_map$alt, "alt")
    setnames_safe(dt, col_map$gene, "gene")
    setnames_safe(dt, col_map$beta, "beta")
    setnames_safe(dt, col_map$se, "se")
    setnames_safe(dt, col_map$pvalue, "pvalue")

    if (!is.na(col_map$maf) && col_map$maf %in% names(dt)) {
        setnames(dt, col_map$maf, "maf")
    } else {
        dt$maf <- NA_real_
    }

    # Handle Z-score to beta conversion (eQTLGen)
    if (isTRUE(col_map$beta_is_z)) {
        # Approximate beta from Z and sample size
        # beta ≈ Z / sqrt(N) for standardized traits
        if ("NrSamples" %in% names(dt)) {
            dt$se <- 1 / sqrt(dt$NrSamples)
            dt$beta <- dt$beta * dt$se
        }
    }

    # Create variant_id if missing
    if (!"variant_id" %in% names(dt)) {
        dt$variant_id <- paste(dt$chr, dt$pos, dt$ref, dt$alt, sep = ":")
    }

    # Standardize chromosome format
    dt$chr <- gsub("^chr", "", dt$chr)

    # Set input build
    dt$input_build <- col_map$build

    # Select and order standard columns
    std_cols <- c("variant_id", "chr", "pos", "ref", "alt", "gene",
                  "beta", "se", "pvalue", "maf", "input_build")

    # Keep additional useful columns if present
    extra_cols <- intersect(names(dt), c("pip", "cs", "n", "n_samples", "NrSamples"))

    keep_cols <- intersect(c(std_cols, extra_cols), names(dt))
    dt <- dt[, ..keep_cols]

    cat("    Loaded:", nrow(dt), "associations\n")

    return(dt)
}

liftover_positions <- function(dt, from_build, to_build, chain_file) {
    if (!has_rtracklayer) {
        warning("rtracklayer not available, skipping liftover")
        return(dt)
    }

    if (from_build == to_build) {
        return(dt)
    }

    library(rtracklayer)
    library(GenomicRanges)

    cat("  Liftover:", from_build, "→", to_build, "\n")

    # Load chain file
    chain <- import.chain(chain_file)

    # Create GRanges
    gr <- GRanges(
        seqnames = paste0("chr", dt$chr),
        ranges = IRanges(start = dt$pos, width = 1),
        strand = "*"
    )

    # Liftover
    lifted <- liftOver(gr, chain)

    # Extract new positions
    n_mapped <- sum(lengths(lifted) == 1)
    cat("    Mapped:", n_mapped, "/", nrow(dt), "(", round(100 * n_mapped / nrow(dt), 1), "%)\n")

    # Filter to successfully lifted
    success <- lengths(lifted) == 1
    dt <- dt[success]
    dt$pos <- start(unlist(lifted[success]))

    return(dt)
}

deduplicate_qtl <- function(dt, strategy = "best_p") {
    # Deduplicate overlapping QTL associations
    # Key: variant_id + gene

    n_before <- nrow(dt)

    if (strategy == "none") {
        return(dt)
    }

    # Create dedup key
    dt$dedup_key <- paste(dt$variant_id, dt$gene, sep = "_")

    if (strategy == "best_p") {
        # Keep association with best p-value
        dt <- dt[order(pvalue)]
        dt <- dt[!duplicated(dedup_key)]
    } else if (strategy == "first") {
        dt <- dt[!duplicated(dedup_key)]
    }

    dt$dedup_key <- NULL

    n_after <- nrow(dt)
    if (n_before > n_after) {
        cat("  Deduplicated:", n_before, "→", n_after, "associations\n")
    }

    return(dt)
}

# ============================================================================
# Main Processing
# ============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                    QTL Harmonization Pipeline                    ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# Collect all input files
input_files <- list()
metadata <- list()

if (!is.null(opt$manifest)) {
    # Read manifest file
    cat("Reading manifest:", opt$manifest, "\n")
    manifest <- fread(opt$manifest)

    # Expected columns: path, source, tissue, ancestry, qtl_type
    for (i in 1:nrow(manifest)) {
        input_files[[i]] <- manifest$path[i]
        metadata[[i]] <- as.list(manifest[i, ])
    }
} else if (!is.null(opt$input)) {
    if (dir.exists(opt$input)) {
        # Directory of files
        files <- list.files(opt$input, pattern = "\\.(tsv|txt|gz)$", full.names = TRUE)
        for (f in files) {
            input_files[[length(input_files) + 1]] <- f
            metadata[[length(metadata) + 1]] <- list(source = basename(opt$input))
        }
    } else {
        # Single file
        input_files[[1]] <- opt$input
        metadata[[1]] <- list(source = "user_input")
    }
}

cat("Input files:", length(input_files), "\n\n")

# Process each file
all_qtl <- list()

for (i in seq_along(input_files)) {
    file_path <- input_files[[i]]
    file_meta <- metadata[[i]]

    cat("Processing file", i, "/", length(input_files), "\n")

    # Build custom column mapping if specified
    custom_cols <- NULL
    if (opt$format == "custom") {
        custom_cols <- list(
            variant = opt$col_variant,
            chr = opt$col_chr,
            pos = opt$col_pos,
            ref = opt$col_ref,
            alt = opt$col_alt,
            gene = opt$col_gene,
            beta = opt$col_beta,
            se = opt$col_se,
            pvalue = opt$col_pvalue,
            maf = NA,
            build = opt$input_build
        )
    }

    # Read file
    dt <- read_qtl_file(file_path, opt$format, custom_cols)

    # Add metadata
    dt$source <- file_meta$source %||% basename(file_path)
    dt$tissue <- file_meta$tissue %||% NA_character_
    dt$ancestry <- file_meta$ancestry %||% NA_character_
    dt$qtl_type <- file_meta$qtl_type %||% "eqtl"

    # Liftover if needed
    if (dt$input_build[1] != opt$output_build && !is.null(opt$chain_file)) {
        dt <- liftover_positions(dt, dt$input_build[1], opt$output_build, opt$chain_file)
    }
    dt$build <- opt$output_build

    # Apply filters
    if (!is.null(opt$tissues)) {
        tissues <- strsplit(opt$tissues, ",")[[1]]
        dt <- dt[tissue %in% tissues]
    }

    if (!is.null(opt$ancestries)) {
        ancestries <- strsplit(opt$ancestries, ",")[[1]]
        dt <- dt[ancestry %in% ancestries]
    }

    if (!is.null(opt$qtl_types)) {
        qtl_types <- strsplit(opt$qtl_types, ",")[[1]]
        dt <- dt[qtl_type %in% qtl_types]
    }

    if (!is.null(opt$genes)) {
        gene_list <- fread(opt$genes, header = FALSE)$V1
        dt <- dt[gene %in% gene_list]
    }

    if (!is.null(opt$region)) {
        parts <- strsplit(opt$region, "[:-]")[[1]]
        reg_chr <- parts[1]
        reg_start <- as.integer(parts[2])
        reg_end <- as.integer(parts[3])
        dt <- dt[chr == reg_chr & pos >= reg_start & pos <= reg_end]
    }

    if (opt$pvalue_threshold < 1) {
        dt <- dt[pvalue <= opt$pvalue_threshold]
    }

    if (opt$maf_threshold > 0 && "maf" %in% names(dt)) {
        dt <- dt[is.na(maf) | maf >= opt$maf_threshold]
    }

    cat("    After filters:", nrow(dt), "associations\n")

    all_qtl[[i]] <- dt
}

# Combine all QTL data
cat("\nCombining", length(all_qtl), "datasets...\n")
combined <- rbindlist(all_qtl, fill = TRUE)
cat("  Total associations:", nrow(combined), "\n")

# Deduplicate
combined <- deduplicate_qtl(combined, opt$dedup_strategy)

# ============================================================================
# Write Output
# ============================================================================
cat("\nWriting output...\n")

# Final column order
final_cols <- c("variant_id", "chr", "pos", "ref", "alt", "gene",
                "beta", "se", "pvalue", "maf", "n",
                "tissue", "ancestry", "qtl_type", "source", "build")
final_cols <- intersect(final_cols, names(combined))

# Add PIP if present (from fine-mapping)
if ("pip" %in% names(combined)) {
    final_cols <- c(final_cols, "pip")
}

combined <- combined[, ..final_cols]

# Split output if requested
if (!is.null(opt$split_by)) {
    split_var <- opt$split_by
    if (!split_var %in% names(combined)) {
        warning("Split variable '", split_var, "' not in data, writing single file")
        split_var <- NULL
    }
}

if (!is.null(opt$split_by) && opt$split_by %in% names(combined)) {
    split_vals <- unique(combined[[opt$split_by]])
    for (val in split_vals) {
        subset_dt <- combined[get(opt$split_by) == val]
        out_file <- paste0(opt$output, ".", val, ".harmonized.tsv.gz")
        fwrite(subset_dt, out_file, sep = "\t", compress = "gzip")
        cat("  Written:", out_file, "(", nrow(subset_dt), "rows)\n")
    }
} else {
    out_file <- paste0(opt$output, ".harmonized.tsv.gz")
    fwrite(combined, out_file, sep = "\t", compress = "gzip")
    cat("  Written:", out_file, "(", nrow(combined), "rows)\n")
}

# Write summary statistics
summary_stats <- combined[, .(
    n_associations = .N,
    n_genes = uniqueN(gene),
    n_variants = uniqueN(variant_id),
    min_pvalue = min(pvalue, na.rm = TRUE),
    median_pvalue = median(pvalue, na.rm = TRUE)
), by = .(source, tissue, ancestry, qtl_type)]

summary_file <- paste0(opt$output, ".summary.tsv")
fwrite(summary_stats, summary_file, sep = "\t")
cat("  Summary:", summary_file, "\n")

# ============================================================================
# Summary
# ============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                         HARMONIZATION COMPLETE                   ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║ Total associations: %-44d ║\n", nrow(combined)))
cat(sprintf("║ Unique variants: %-47d ║\n", uniqueN(combined$variant_id)))
cat(sprintf("║ Unique genes: %-50d ║\n", uniqueN(combined$gene)))
cat(sprintf("║ Genome build: %-50s ║\n", opt$output_build))
cat("║                                                                  ║\n")
cat("║ By QTL type:                                                     ║\n")
for (qt in unique(combined$qtl_type)) {
    n <- sum(combined$qtl_type == qt)
    cat(sprintf("║   %-8s: %-52d ║\n", qt, n))
}
cat("╚══════════════════════════════════════════════════════════════════╝\n")

cat("\nDone!\n")
