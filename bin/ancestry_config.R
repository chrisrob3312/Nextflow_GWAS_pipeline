#!/usr/bin/env Rscript

# ============================================================================
# Ancestry Group Configuration for Admixed Latino Cohort
# ============================================================================
# GRAF-ANC Codes:
#   AAC = 107 (African American)
#   AFR = 100 (African)
#   EUR = 300 (European)
#   EAS = 500 (East Asian)
#   SAS = (South Asian)
#   AMR = 600 (Native American/Amerindigenous)
#   LAT1 = 601 (Latino type 1 - higher Indigenous)
#   LAT2 = 602 (Latino type 2 - higher European)
#
# Cohort Structure (>60% Latino):
#   - EUR: European (run separately)
#   - AAC: African American (run separately)
#   - LAT1: Latino high-Indigenous (run separately)
#   - LAT2: Latino high-European (run separately)
#   - EAS: East Asian (run separately IF N >= 30, else pool to OTHER)
#   - SAS: South Asian (run separately IF N >= 30, else pool to OTHER)
#   - OTHER: Small groups (<30) - pooled/meta-analyzed
#
# DYNAMIC: Groups with N >= MIN_STRATUM_N run separately
#          Groups with N < MIN_STRATUM_N pool to OTHER
# ============================================================================

# Minimum sample size for separate analysis
MIN_STRATUM_N <- 30

# Guaranteed primary ancestries (always attempt separate analysis)
PRIMARY_ANCESTRIES <- c("EUR", "AAC", "LAT1", "LAT2")

# Conditional ancestries (run separately if N >= MIN_STRATUM_N)
CONDITIONAL_ANCESTRIES <- c("EAS", "SAS")

# Groups that always get pooled into OTHER
ALWAYS_POOL <- c("AMR", "AFR", "UNKNOWN")

# GRAF-ANC code mapping
GRAF_ANC_CODES <- list(
    "107" = "AAC",
    "100" = "AFR",
    "300" = "EUR",
    "500" = "EAS",
    "600" = "AMR",
    "601" = "LAT1",
    "602" = "LAT2"
)

#' Assign samples to analysis strata based on ancestry and sample size
#'
#' @param ancestry_data Data frame with sample IDs and ancestry assignments
#' @param ancestry_col Column name with ancestry labels (GRAF-ANC codes or names)
#' @param min_n Minimum sample size for separate stratum (default: 30)
#' @return Data frame with analysis_stratum column added
assign_analysis_strata <- function(ancestry_data, ancestry_col = "ancestry",
                                   min_n = MIN_STRATUM_N) {

    # Convert GRAF-ANC codes to names if needed
    if (is.numeric(ancestry_data[[ancestry_col]]) ||
        all(ancestry_data[[ancestry_col]] %in% names(GRAF_ANC_CODES))) {
        ancestry_data$ancestry_name <- sapply(as.character(ancestry_data[[ancestry_col]]),
                                               function(x) GRAF_ANC_CODES[[x]] %||% "OTHER")
    } else {
        ancestry_data$ancestry_name <- ancestry_data[[ancestry_col]]
    }

    # Count samples per ancestry
    anc_counts <- table(ancestry_data$ancestry_name)

    cat("Ancestry distribution:\n")
    for (anc in names(sort(anc_counts, decreasing = TRUE))) {
        cat("  ", anc, ":", anc_counts[anc], "\n")
    }
    cat("\n")

    # Assign analysis strata
    ancestry_data$analysis_stratum <- sapply(ancestry_data$ancestry_name, function(anc) {
        n <- anc_counts[anc]
        if (is.na(n)) n <- 0

        # Primary ancestries always get their own stratum if N >= min_n
        if (anc %in% PRIMARY_ANCESTRIES) {
            if (n >= min_n) {
                return(anc)
            } else {
                cat("  Warning:", anc, "has N =", n,
                    "< min_n, pooling to OTHER\n")
                return("OTHER")
            }
        }

        # Conditional ancestries (EAS, SAS): run separately if N >= min_n
        if (anc %in% CONDITIONAL_ANCESTRIES) {
            if (n >= min_n) {
                cat("  ✓", anc, "has N =", n, ">= min_n, running separately\n")
                return(anc)
            } else {
                cat("  →", anc, "has N =", n, "< min_n, pooling to OTHER\n")
                return("OTHER")
            }
        }

        # Always-pool ancestries go to OTHER
        if (anc %in% ALWAYS_POOL) {
            return("OTHER")
        }

        # Unknown groups: pool if small
        if (n < min_n) {
            return("OTHER")
        }

        return(anc)
    })

    # Report final strata
    strata_counts <- table(ancestry_data$analysis_stratum)
    cat("\nFinal analysis strata:\n")
    for (strat in names(sort(strata_counts, decreasing = TRUE))) {
        cat("  ", strat, ":", strata_counts[strat], "\n")
    }

    # List what's in OTHER
    if ("OTHER" %in% names(strata_counts)) {
        other_breakdown <- table(ancestry_data$ancestry_name[ancestry_data$analysis_stratum == "OTHER"])
        cat("\n  OTHER contains:\n")
        for (anc in names(other_breakdown)) {
            cat("    ", anc, ":", other_breakdown[anc], "\n")
        }
    }

    return(ancestry_data)
}

#' Get strata to run for a given analysis type
#'
#' @param ancestry_data Data frame with analysis_stratum column
#' @param include_pooled Include pooled/joint analysis across all samples
#' @param include_other Include OTHER stratum (for meta-analysis)
#' @return Character vector of strata to analyze
get_analysis_strata <- function(ancestry_data, include_pooled = TRUE,
                                 include_other = TRUE) {
    strata <- unique(ancestry_data$analysis_stratum)

    # Always include primary ancestries with sufficient N
    primary_strata <- intersect(strata, PRIMARY_ANCESTRIES)

    # Optionally include OTHER
    if (include_other && "OTHER" %in% strata) {
        primary_strata <- c(primary_strata, "OTHER")
    }

    # Optionally include pooled analysis
    if (include_pooled) {
        primary_strata <- c("POOLED", primary_strata)
    }

    return(primary_strata)
}

#' Meta-analyze results across strata
#'
#' @param results_list List of results data.tables, one per stratum
#' @param method Meta-analysis method: "fixed" or "random"
#' @return Combined results with meta-analysis
meta_analyze_strata <- function(results_list, method = "random") {
    if (!requireNamespace("metafor", quietly = TRUE)) {
        cat("metafor package not installed - using simple pooling\n")
        # Simple inverse-variance weighted pooling
        return(simple_pool_results(results_list))
    }

    library(metafor)

    # Combine all results
    all_results <- rbindlist(results_list, idcol = "stratum", fill = TRUE)

    # Meta-analyze per variant/gene
    if ("variant_id" %in% names(all_results)) {
        group_col <- "variant_id"
    } else if ("gene" %in% names(all_results)) {
        group_col <- "gene"
    } else {
        group_col <- names(all_results)[1]
    }

    # Perform meta-analysis
    meta_results <- all_results[, {
        if (.N >= 2 && all(!is.na(beta)) && all(!is.na(se))) {
            # Random effects meta-analysis
            ma <- tryCatch({
                rma(yi = beta, sei = se, method = ifelse(method == "random", "REML", "FE"))
            }, error = function(e) NULL)

            if (!is.null(ma)) {
                list(
                    beta_meta = ma$beta[1],
                    se_meta = ma$se,
                    pval_meta = ma$pval,
                    i2 = ma$I2,
                    n_strata = .N,
                    strata = paste(stratum, collapse = ";")
                )
            } else {
                list(beta_meta = NA, se_meta = NA, pval_meta = NA,
                     i2 = NA, n_strata = .N, strata = paste(stratum, collapse = ";"))
            }
        } else if (.N == 1) {
            list(
                beta_meta = beta,
                se_meta = se,
                pval_meta = pvalue,
                i2 = NA,
                n_strata = 1,
                strata = stratum
            )
        } else {
            list(beta_meta = NA, se_meta = NA, pval_meta = NA,
                 i2 = NA, n_strata = .N, strata = paste(stratum, collapse = ";"))
        }
    }, by = group_col]

    return(meta_results)
}

#' Simple inverse-variance weighted pooling (fallback)
simple_pool_results <- function(results_list) {
    all_results <- rbindlist(results_list, idcol = "stratum", fill = TRUE)

    if ("variant_id" %in% names(all_results)) {
        group_col <- "variant_id"
    } else {
        group_col <- names(all_results)[1]
    }

    pooled <- all_results[, {
        weights <- 1 / se^2
        weights <- weights / sum(weights, na.rm = TRUE)

        list(
            beta_meta = sum(beta * weights, na.rm = TRUE),
            se_meta = sqrt(1 / sum(1/se^2, na.rm = TRUE)),
            n_strata = .N,
            strata = paste(stratum, collapse = ";")
        )
    }, by = group_col]

    pooled[, pval_meta := 2 * pnorm(-abs(beta_meta / se_meta))]

    return(pooled)
}

# Export for use in other scripts
if (exists("opt") && !is.null(opt$export_config)) {
    cat("\nExporting ancestry configuration...\n")
    config <- list(
        min_stratum_n = MIN_STRATUM_N,
        primary_ancestries = PRIMARY_ANCESTRIES,
        poolable_ancestries = POOLABLE_ANCESTRIES,
        graf_anc_codes = GRAF_ANC_CODES
    )
    saveRDS(config, opt$export_config)
}
