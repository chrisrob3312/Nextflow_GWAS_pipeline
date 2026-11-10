#!/usr/bin/env Rscript

# ============================================================================
# Meta-analyse Tractor-GENESIS results across ancestry strata -> POOLED
# ============================================================================
# Combines per-stratum summary statistics (EUR, AAC, LAT1, LAT2, EAS, SAS,
# OTHER) into one pooled result per variant using inverse-variance
# random-effects meta-analysis (metafor if installed, else IVW fallback).
#
# Tractor-GENESIS output has one BETA/SE per ancestry background per stratum
# (BETA_EUR, BETA_AFR, BETA_AMR ...). The same background is meta-analysed
# ACROSS strata (e.g. the AMR-background effect in LAT1 + LAT2 + OTHER), so
# the pooled file keeps the Tractor structure and stays comparable with the
# per-stratum files. BETA_TOTAL (the total-dose effect) is pooled the same way, and the output is
# sorted on its pooled P_TOTAL; P_JOINT is the strata's P_JOINT combined with Fisher's method, a
# descriptive summary. A variant (or a background) in fewer than two strata is not pooled: NA,
# never one stratum's estimate presented as pooled. METHOD says how each row was pooled.
#
# Usage:
#   meta_analyze_strata.R --input_dir results/gwas/OS --strata EUR,AAC,LAT1,LAT2,OTHER \
#                         --output results/gwas/os/POOLED.ivw_descriptive.sumstats.gz
# The output starts with "##" lines: strata= (pooled), strata_skipped= (below MIN_STRATUM_N, from
# --strata_skipped), backgrounds= (the ancestry backgrounds pooled), method=, sorted_on=, covariates=.
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(optparse) })
local({   # sumstats_meta / fwrite_sumstats ("## key=value" lines), from next to this script
    f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
    source(file.path(if (length(f) > 0) dirname(normalizePath(f[1])) else "bin", "sumstats_io.R"))
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--input_dir", type = "character"),
    make_option("--strata", type = "character", default = "EUR,AAC,LAT1,LAT2,EAS,SAS,OTHER"),
    make_option("--pattern", type = "character", default = "{stratum}.sumstats.gz",
                help = "File name pattern inside input_dir [default: {stratum}.sumstats.gz]"),
    make_option("--method", type = "character", default = "random", help = "random or fixed"),
    make_option("--min_strata", type = "integer", default = 2),
    make_option("--strata_skipped", type = "character", default = "",
                help = "Strata left out because they were skipped (below MIN_STRATUM_N); recorded in the header only"),
    make_option("--output", type = "character", default = "POOLED.ivw_descriptive.sumstats.gz"),
    make_option("--threads", type = "integer", default = 4)
)))

setDTthreads(opt$threads)
if (!opt$method %in% c("random", "fixed")) stop("--method must be random or fixed")
has_metafor <- requireNamespace("metafor", quietly = TRUE)
strata <- strsplit(opt$strata, ",")[[1]]
method_label <- if (opt$method == "fixed") "IVW fixed effect" else if (has_metafor) "random effects (metafor REML)" else
    "IVW fixed effect (metafor not installed; --method random asked for)"

cat("Meta-analysing strata:", paste(strata, collapse = ", "), "\n")
cat("Method:", method_label, "\n\n")

# Load strata. Each table carries its STRATUM column (the adapter's first column); it must name
# the stratum the file is pooled as.
dt_list <- list(); cov_by_stratum <- c()
for (s in strata) {
    f <- file.path(opt$input_dir, gsub("\\{stratum\\}", s, opt$pattern))
    if (!file.exists(f)) { cat("  ", s, ": missing (", f, ")\n"); next }
    cov_by_stratum[s] <- sumstats_meta_value(sumstats_meta(f), "covariates")
    d <- fread(f)
    if ("STRATUM" %in% names(d)) {
        other <- setdiff(unique(as.character(d$STRATUM)), s)
        if (length(other)) stop(f, " is pooled as stratum ", s, " but its STRATUM column says ", paste(other, collapse = ", "))
    } else {
        d[, STRATUM := s]
    }
    dt_list[[s]] <- d
    cat("  ", s, ":", nrow(d), "variants\n")
}
if (length(dt_list) < opt$min_strata) stop("Fewer than ", opt$min_strata, " strata found")
all <- rbindlist(dt_list, fill = TRUE)

# Ancestry backgrounds from the BETA_<anc> columns (BETA_TOTAL is the total-dose effect, not one)
bg <- setdiff(sub("^BETA_", "", grep("^BETA_", names(all), value = TRUE)), "TOTAL")
has_total <- all(c("BETA_TOTAL", "SE_TOTAL") %in% names(all))
cat("\nAncestry backgrounds:", paste(bg, collapse = ", "), "\n")
if (!has_total) cat("NOTE: no BETA_TOTAL / SE_TOTAL (tables from before CP-19): sorted on the Fisher P_JOINT instead\n")

# pool(b, se): estimate over the strata with an estimate. Fewer than two contributing strata: NA,
# never one stratum's estimate presented as pooled (CP-23). Records the method actually used.
pool <- function(b, se) {
    ok <- !is.na(b) & !is.na(se) & se > 0
    k <- sum(ok)
    if (k < 2) return(list(beta = NA_real_, se = NA_real_, p = NA_real_, i2 = NA_real_, n = k, method = NA_character_))
    if (opt$method == "random" && has_metafor) {
        m <- tryCatch(metafor::rma(yi = b[ok], sei = se[ok], method = "REML"), error = function(e) NULL)
        if (!is.null(m)) return(list(beta = as.numeric(m$beta), se = m$se, p = m$pval, i2 = m$I2, n = k, method = "REML"))
    }
    w <- 1 / se[ok]^2; bm <- sum(w * b[ok]) / sum(w); sem <- sqrt(1 / sum(w))
    Q <- sum(w * (b[ok] - bm)^2); df <- k - 1
    list(beta = bm, se = sem, p = 2 * pnorm(-abs(bm / sem)), i2 = if (Q > 0) max(0, (Q - df) / Q * 100) else 0, n = k,
         method = if (opt$method == "random" && has_metafor) "IVW (REML failed)" else "IVW")
}

fisher_p <- function(p) {
    p <- p[!is.na(p) & p > 0]
    if (length(p) == 0) return(NA_real_)
    pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE)
}

cat("Running meta-analysis on", length(unique(all$ID)), "variants...\n")
pooled <- all[, {
    out <- list(CHR = CHR[1], POS = POS[1], REF = REF[1], ALT = ALT[1],
                N = as.numeric(sum(N, na.rm = TRUE)), N_STRATA = .N, STRATA = paste(STRATUM, collapse = ";"))
    if (has_total) {
        m <- pool(BETA_TOTAL, SE_TOTAL)
        out$BETA_TOTAL <- m$beta; out$SE_TOTAL <- m$se; out$P_TOTAL <- m$p
        out$I2_STRATA_TOTAL <- m$i2; out$N_STRATA_TOTAL <- m$n
        out$METHOD <- m$method
    } else out$METHOD <- NA_character_
    out$P_JOINT <- if (.N >= 2) fisher_p(P_JOINT) else NA_real_       # Fisher's combination: descriptive
    for (a in bg) {
        b <- get(paste0("BETA_", a)); s <- get(paste0("SE_", a))
        m <- pool(b, s)
        contributing <- !is.na(b) & !is.na(s) & s > 0
        mac_col <- paste0("MAC_", a)
        # MAC over the strata that contribute to this background's pooled estimate only
        out[[mac_col]]  <- if (mac_col %in% names(.SD) && any(contributing)) as.numeric(sum(get(mac_col)[contributing], na.rm = TRUE)) else NA_real_
        out[[paste0("BETA_", a)]] <- m$beta
        out[[paste0("SE_", a)]]   <- m$se
        out[[paste0("Z_", a)]]    <- if (!is.na(m$se) && m$se > 0) m$beta / m$se else NA_real_
        out[[paste0("P_", a)]]    <- m$p
        out[[paste0("I2_STRATA_", a)]] <- m$i2
        out[[paste0("N_STRATA_", a)]]  <- m$n
    }
    out
}, by = ID]

# Heterogeneity across ancestry BACKGROUNDS of the pooled estimates: a naive Cochran's Q that
# ignores the covariance between backgrounds (the adapter's Q_HET / P_HET / I2 use it), hence the
# _NAIVE names, so the two never share a column name (CP-23)
pooled[, c("Q_HET_NAIVE", "DF_HET_NAIVE", "P_HET_NAIVE", "I2_NAIVE") := {
    b <- unlist(.SD[, paste0("BETA_", bg), with = FALSE]); s <- unlist(.SD[, paste0("SE_", bg), with = FALSE])
    ok <- !is.na(b) & !is.na(s) & s > 0
    if (sum(ok) >= 2) { w <- 1 / s[ok]^2; bm <- sum(w * b[ok]) / sum(w); Q <- sum(w * (b[ok] - bm)^2); df <- sum(ok) - 1
        list(Q, df, pchisq(Q, df, lower.tail = FALSE), if (Q > 0) max(0, (Q - df) / Q * 100) else 0) } else list(NA_real_, NA_integer_, NA_real_, NA_real_)
}, by = ID]

# Sorted, and tallied, on the pooled P (the total-dose effect pooled across strata), not on the
# Fisher P_JOINT
sort_col <- if (has_total) "P_TOTAL" else "P_JOINT"
setorderv(pooled, sort_col, na.last = TRUE)
# The strata's "## covariates=" lines (CP-08): one line when they agree, one per stratum otherwise
cov_known <- cov_by_stratum[!is.na(cov_by_stratum)]
meta <- if (length(cov_known) == 0) character() else if (length(unique(cov_known)) == 1 && length(cov_known) == length(cov_by_stratum)) {
    c(covariates = unname(cov_known[1]))
} else setNames(ifelse(is.na(cov_by_stratum), "unknown", cov_by_stratum), paste0("covariates_", names(cov_by_stratum)))
# Scope (CP-22) and method (CP-23)
meta <- c(strata = paste(names(dt_list), collapse = ","),
          strata_skipped = if (nzchar(opt$strata_skipped)) opt$strata_skipped else "none",
          backgrounds = paste(bg, collapse = ","),
          method = paste0(method_label, "; a variant or background in fewer than two strata is not pooled (NA)"),
          sorted_on = paste0(sort_col, if (has_total) " (pooled total-dose effect)" else " (Fisher; no BETA_TOTAL in the inputs)"),
          P_JOINT = "Fisher's combination of the strata's P_JOINT (descriptive)",
          meta)
fwrite_sumstats(pooled, opt$output, meta)
cat("\nPooled results:", opt$output, "\n")
cat("  Variants:", nrow(pooled), "| pooled (in two or more strata):", sum(pooled$N_STRATA >= 2), "\n")
cat(sprintf("  %s < 5e-8: %d | < 1e-5: %d\n", sort_col, sum(pooled[[sort_col]] < 5e-8, na.rm = TRUE), sum(pooled[[sort_col]] < 1e-5, na.rm = TRUE)))
