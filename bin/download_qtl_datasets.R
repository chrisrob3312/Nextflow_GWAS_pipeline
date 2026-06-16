#!/usr/bin/env Rscript

# =============================================================================
# QTL Dataset Download and Setup
# =============================================================================
# Auto-downloads publicly accessible QTL datasets from curated_qtl_sources.yml
# Handles: eQTL, sQTL, mQTL, pQTL, caQTL, hQTL
# Supported sources: GTEx, eQTLGen, eQTL Catalogue, GoDMC, etc.

suppressPackageStartupMessages({
    library(optparse)
    library(yaml)
    library(data.table)
    library(httr)
    library(R.utils)
})

# =============================================================================
# Command-line options
# =============================================================================
option_list <- list(
    make_option(c("--config"), type = "character", default = NULL,
                help = "QTL sources YAML config file"),
    make_option(c("--output_dir"), type = "character", default = "qtl_data",
                help = "Output directory for downloaded data"),
    make_option(c("--qtl_types"), type = "character", default = "eqtl,sqtl,pqtl,mqtl",
                help = "Comma-separated QTL types to download"),
    make_option(c("--tissues"), type = "character", default = NULL,
                help = "Comma-separated tissues to download (NULL = all)"),
    make_option(c("--ancestries"), type = "character", default = NULL,
                help = "Comma-separated ancestries to prioritize"),
    make_option(c("--priority_leukemia"), type = "logical", default = TRUE,
                help = "Prioritize leukemia-relevant datasets"),
    make_option(c("--public_only"), type = "logical", default = TRUE,
                help = "Only download publicly accessible datasets"),
    make_option(c("--max_parallel"), type = "integer", default = 4,
                help = "Max parallel downloads"),
    make_option(c("--genome_build"), type = "character", default = "GRCh38",
                help = "Target genome build (GRCh37 or GRCh38)"),
    make_option(c("--verbose"), action = "store_true", default = FALSE,
                help = "Verbose output"),
    make_option(c("--dry_run"), action = "store_true", default = FALSE,
                help = "List datasets to download without downloading"),
    make_option(c("--threads"), type = "integer", default = 4,
                help = "Number of threads")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

verbose <- opt$verbose
if (verbose) {
    cat("╔══════════════════════════════════════════════════════════════════════╗\n")
    cat("║            QTL Dataset Download and Setup                            ║\n")
    cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")
}

# =============================================================================
# Download registry - known download URLs for public datasets
# =============================================================================
DOWNLOAD_REGISTRY <- list(
    # GTEx v8
    gtex_v8 = list(
        base_url = "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data",
        eqtl_pattern = "GTEx_Analysis_v8_eQTL/{tissue}.v8.signif_variant_gene_pairs.txt.gz",
        sqtl_pattern = "GTEx_Analysis_v8_sQTL/{tissue}.v8.sqtl_signifpairs.txt.gz",
        tissues = c("Whole_Blood", "Cells_EBV-transformed_lymphocytes",
                   "Spleen", "Bone_Marrow", # Note: Bone_Marrow not in GTEx
                   "Liver", "Lung", "Heart_Left_Ventricle")
    ),

    # eQTLGen Phase II
    eqtlgen_blood = list(
        base_url = "https://www.eqtlgen.org/phase2",
        cis_eqtl = "cis-eQTL_significant_eGenes.txt.gz",
        trans_eqtl = "trans-eQTL_significant.txt.gz"
    ),

    # eQTL Catalogue (EBI)
    eqtl_catalogue = list(
        base_url = "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master",
        tabix_manifest = "data_tables/dataset_metadata.tsv",
        api_base = "https://www.ebi.ac.uk/eqtl/api/v2"
    ),

    # GoDMC (methylation QTL)
    godmc = list(
        base_url = "http://www.godmc.org.uk/data",
        cis_mqtl = "cis-mQTLs.txt.gz",
        trans_mqtl = "trans-mQTLs.txt.gz"
    ),

    # BLUEPRINT
    blueprint_blood = list(
        # EGA controlled access
        access = "controlled_ega",
        ega_datasets = c("EGAD00001005199", "EGAD00001005200")
    ),

    # OneK1K
    onek1k_pbmc = list(
        github = "https://github.com/powellgenomicslab/onek1k_phase1",
        geo = "GSE139324"
    ),

    # DICE
    dice_immune = list(
        base_url = "https://dice-database.org/download",
        api_base = "https://dice-database.org/api"
    )
)

# =============================================================================
# Helper functions
# =============================================================================

download_with_retry <- function(url, destfile, max_retries = 3,
                                 timeout = 300, verbose = FALSE) {
    for (attempt in 1:max_retries) {
        if (verbose) cat("  Attempt", attempt, "of", max_retries, "...\n")

        tryCatch({
            response <- GET(url,
                           write_disk(destfile, overwrite = TRUE),
                           timeout(timeout),
                           progress())

            if (status_code(response) == 200) {
                return(TRUE)
            } else {
                warning("HTTP status: ", status_code(response))
            }
        }, error = function(e) {
            if (attempt < max_retries) {
                Sys.sleep(2^attempt)  # Exponential backoff
            }
        })
    }
    return(FALSE)
}

create_directory_structure <- function(base_dir) {
    dirs <- c(
        file.path(base_dir, "eqtl"),
        file.path(base_dir, "sqtl"),
        file.path(base_dir, "mqtl"),
        file.path(base_dir, "pqtl"),
        file.path(base_dir, "caqtl"),
        file.path(base_dir, "hqtl"),
        file.path(base_dir, "metadata"),
        file.path(base_dir, "logs")
    )

    for (d in dirs) {
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }

    return(dirs)
}

filter_datasets <- function(config, qtl_types, tissues = NULL,
                           ancestries = NULL, priority_leukemia = TRUE,
                           public_only = TRUE) {
    selected <- list()
    qtl_types <- strsplit(qtl_types, ",")[[1]]

    for (section_name in names(config)) {
        section <- config[[section_name]]

        if (!is.list(section)) next

        for (dataset_name in names(section)) {
            dataset <- section[[dataset_name]]

            if (!is.list(dataset)) next

            # Check if enabled
            if (!isTRUE(dataset$enabled)) next

            # Check access level
            if (public_only && !is.null(dataset$access)) {
                if (!dataset$access %in% c("public", "public_restricted")) {
                    next
                }
            }

            # Check QTL types
            if (!is.null(dataset$types)) {
                if (!any(dataset$types %in% qtl_types)) next
            }

            # Check priority for leukemia
            if (priority_leukemia) {
                if (!is.null(dataset$priority_leukemia) &&
                    dataset$priority_leukemia <= 2) {
                    dataset$priority_score <- dataset$priority_leukemia
                } else {
                    dataset$priority_score <- 10
                }
            }

            # Check tissues if specified
            if (!is.null(tissues)) {
                tissue_list <- strsplit(tissues, ",")[[1]]
                if (!is.null(dataset$tissues)) {
                    if (!any(dataset$tissues %in% tissue_list)) next
                }
            }

            # Check ancestries if specified
            if (!is.null(ancestries)) {
                anc_list <- strsplit(ancestries, ",")[[1]]
                if (!is.null(dataset$ancestry)) {
                    dataset_anc <- names(dataset$ancestry)
                    if (!any(dataset_anc %in% anc_list)) {
                        dataset$priority_score <- dataset$priority_score + 5
                    }
                }
            }

            selected[[dataset_name]] <- dataset
            selected[[dataset_name]]$section <- section_name
        }
    }

    # Sort by priority
    if (length(selected) > 0) {
        priorities <- sapply(selected, function(x) x$priority_score %||% 10)
        selected <- selected[order(priorities)]
    }

    return(selected)
}

# =============================================================================
# Download functions for specific sources
# =============================================================================

download_gtex <- function(output_dir, tissues, qtl_types, verbose = FALSE) {
    results <- list()

    gtex_dir <- file.path(output_dir, "gtex_v8")
    dir.create(gtex_dir, recursive = TRUE, showWarnings = FALSE)

    base_url <- "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data"

    for (tissue in tissues) {
        if (verbose) cat("  Downloading GTEx", tissue, "...\n")

        # eQTL
        if ("eqtl" %in% qtl_types) {
            url <- paste0(base_url, "/GTEx_Analysis_v8_eQTL/",
                         tissue, ".v8.signif_variant_gene_pairs.txt.gz")
            destfile <- file.path(gtex_dir, paste0(tissue, ".eqtl.txt.gz"))

            success <- download_with_retry(url, destfile, verbose = verbose)
            results[[paste0("gtex_", tissue, "_eqtl")]] <- success
        }

        # sQTL
        if ("sqtl" %in% qtl_types) {
            url <- paste0(base_url, "/GTEx_Analysis_v8_sQTL/",
                         tissue, ".v8.sqtl_signifpairs.txt.gz")
            destfile <- file.path(gtex_dir, paste0(tissue, ".sqtl.txt.gz"))

            success <- download_with_retry(url, destfile, verbose = verbose)
            results[[paste0("gtex_", tissue, "_sqtl")]] <- success
        }
    }

    return(results)
}

download_eqtlgen <- function(output_dir, verbose = FALSE) {
    results <- list()

    eqtlgen_dir <- file.path(output_dir, "eqtlgen")
    dir.create(eqtlgen_dir, recursive = TRUE, showWarnings = FALSE)

    if (verbose) cat("  Downloading eQTLGen cis-eQTL...\n")

    # Note: eQTLGen requires registration
    # These are placeholder URLs - actual downloads require API access
    cat("    NOTE: eQTLGen requires registration at https://www.eqtlgen.org/\n")
    cat("    Please download manually and place in:", eqtlgen_dir, "\n")

    results$eqtlgen_cis <- "manual_required"
    return(results)
}

download_eqtl_catalogue <- function(output_dir, tissues = NULL, verbose = FALSE) {
    results <- list()

    eqtl_cat_dir <- file.path(output_dir, "eqtl_catalogue")
    dir.create(eqtl_cat_dir, recursive = TRUE, showWarnings = FALSE)

    if (verbose) cat("  Fetching eQTL Catalogue metadata...\n")

    # Download metadata
    metadata_url <- "https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/data_tables/dataset_metadata.tsv"
    metadata_file <- file.path(eqtl_cat_dir, "dataset_metadata.tsv")

    success <- download_with_retry(metadata_url, metadata_file, verbose = verbose)

    if (success && file.exists(metadata_file)) {
        metadata <- fread(metadata_file)

        # Filter for blood/immune tissues relevant to leukemia
        blood_tissues <- c("blood", "PBMC", "monocyte", "neutrophil",
                          "T_cell", "B_cell", "NK", "lymphocyte")

        if (!is.null(tissues)) {
            tissue_list <- strsplit(tissues, ",")[[1]]
            relevant <- metadata[tissue_label %in% tissue_list |
                                tissue_ontology_term_id %like% "CL:"]
        } else {
            relevant <- metadata[tissue_label %like% paste(blood_tissues, collapse = "|")]
        }

        if (verbose) {
            cat("    Found", nrow(relevant), "relevant datasets\n")
        }

        # Write filtered metadata
        fwrite(relevant, file.path(eqtl_cat_dir, "selected_datasets.tsv"), sep = "\t")

        results$eqtl_catalogue_metadata <- TRUE
        results$n_datasets <- nrow(relevant)
    }

    return(results)
}

download_godmc <- function(output_dir, verbose = FALSE) {
    results <- list()

    godmc_dir <- file.path(output_dir, "godmc")
    dir.create(godmc_dir, recursive = TRUE, showWarnings = FALSE)

    if (verbose) cat("  Downloading GoDMC mQTL...\n")

    # GoDMC data portal
    cat("    NOTE: GoDMC data available at http://www.godmc.org.uk/\n")
    cat("    Place downloaded files in:", godmc_dir, "\n")

    results$godmc <- "manual_required"
    return(results)
}

# =============================================================================
# Main execution
# =============================================================================

# Load config
if (is.null(opt$config)) {
    # Look for default config
    default_paths <- c(
        "assets/qtl_datasets/curated_qtl_sources.yml",
        "../assets/qtl_datasets/curated_qtl_sources.yml"
    )

    for (path in default_paths) {
        if (file.exists(path)) {
            opt$config <- path
            break
        }
    }

    if (is.null(opt$config)) {
        stop("No config file specified and default not found")
    }
}

if (verbose) cat("Loading config:", opt$config, "\n")
config <- yaml::read_yaml(opt$config)

# Create output directory structure
create_directory_structure(opt$output_dir)

# Filter datasets
selected <- filter_datasets(
    config = config,
    qtl_types = opt$qtl_types,
    tissues = opt$tissues,
    ancestries = opt$ancestries,
    priority_leukemia = opt$priority_leukemia,
    public_only = opt$public_only
)

if (verbose) {
    cat("\n")
    cat("SELECTED DATASETS (", length(selected), "):\n")
    for (name in names(selected)) {
        ds <- selected[[name]]
        cat("  •", name, "\n")
        cat("    Types:", paste(ds$types, collapse = ", "), "\n")
        if (!is.null(ds$ancestry)) {
            cat("    Ancestry:", paste(names(ds$ancestry), collapse = ", "), "\n")
        }
        if (!is.null(ds$availability)) {
            cat("    Status:", ds$availability, "\n")
        }
    }
    cat("\n")
}

if (opt$dry_run) {
    cat("DRY RUN - No downloads performed\n")
    quit(save = "no", status = 0)
}

# Download datasets
download_results <- list()
qtl_types <- strsplit(opt$qtl_types, ",")[[1]]

# GTEx
if ("gtex_v8" %in% names(selected)) {
    if (verbose) cat("\nDownloading GTEx v8...\n")
    gtex_tissues <- c("Whole_Blood", "Cells_EBV-transformed_lymphocytes")
    download_results$gtex <- download_gtex(
        opt$output_dir, gtex_tissues, qtl_types, verbose
    )
}

# eQTLGen
if ("eqtlgen_blood" %in% names(selected)) {
    if (verbose) cat("\nSetting up eQTLGen...\n")
    download_results$eqtlgen <- download_eqtlgen(opt$output_dir, verbose)
}

# eQTL Catalogue
if ("eqtl_catalogue" %in% names(selected)) {
    if (verbose) cat("\nDownloading eQTL Catalogue...\n")
    download_results$eqtl_catalogue <- download_eqtl_catalogue(
        opt$output_dir, opt$tissues, verbose
    )
}

# GoDMC
if ("godmc" %in% names(selected)) {
    if (verbose) cat("\nSetting up GoDMC...\n")
    download_results$godmc <- download_godmc(opt$output_dir, verbose)
}

# Write download summary
summary_dt <- data.table(
    dataset = names(selected),
    types = sapply(selected, function(x) paste(x$types, collapse = ";")),
    section = sapply(selected, function(x) x$section),
    access = sapply(selected, function(x) x$access %||% "public"),
    availability = sapply(selected, function(x) x$availability %||% "unknown"),
    priority_leukemia = sapply(selected, function(x) x$priority_leukemia %||% NA)
)

fwrite(summary_dt, file.path(opt$output_dir, "metadata", "download_summary.tsv"), sep = "\t")

# Write manifest
manifest <- list(
    download_date = Sys.time(),
    config_file = opt$config,
    qtl_types = qtl_types,
    genome_build = opt$genome_build,
    n_datasets_selected = length(selected),
    datasets = names(selected)
)

yaml::write_yaml(manifest, file.path(opt$output_dir, "metadata", "manifest.yml"))

if (verbose) {
    cat("\n")
    cat("══════════════════════════════════════════════════════════════════════\n")
    cat("DOWNLOAD COMPLETE\n")
    cat("  Output directory:", opt$output_dir, "\n")
    cat("  Datasets selected:", length(selected), "\n")
    cat("  Manifest:", file.path(opt$output_dir, "metadata", "manifest.yml"), "\n")
    cat("══════════════════════════════════════════════════════════════════════\n")
}
