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
# per-stratum files. P_JOINT is combined with Fisher's method as a summary.
#
# Usage:
#   meta_analyze_strata.R --input_dir results/gwas/OS --strata EUR,AAC,LAT1,LAT2,OTHER \
#                         --output results/gwas/OS/POOLED.sumstats.gz
# ============================================================================

suppressPackageStartupMessages({ library(data.table); library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
    make_option("--input_dir", type = "character"),
    make_option("--strata", type = "character", default = "EUR,AAC,LAT1,LAT2,EAS,SAS,OTHER"),
    make_option("--pattern", type = "character", default = "{stratum}.sumstats.gz",
                help = "File name pattern inside input_dir [default: {stratum}.sumstats.gz]"),
    make_option("--method", type = "character", default = "random", help = "random or fixed"),
    make_option("--min_strata", type = "integer", default = 2),
    make_option("--output", type = "character", default = "POOLED.sumstats.gz"),
    make_option("--threads", type = "integer", default = 4)
)))

setDTthreads(opt$threads)
has_metafor <- requireNamespace("metafor", quietly = TRUE)
strata <- strsplit(opt$strata, ",")[[1]]

cat("Meta-analysing strata:", paste(strata, collapse = ", "), "\n")
cat("Method:", if (has_metafor) paste0(opt$method, "-effects (metafor)") else "IVW fixed-effect fallback", "\n\n")

# Load strata
dt_list <- list()
for (s in strata) {
    f <- file.path(opt$input_dir, gsub("\\{stratum\\}", s, opt$pattern))
    if (!file.exists(f)) { cat("  ", s, ": missing (", f, ")\n"); next }
    d <- fread(f)
    d[, STRATUM := s]
    dt_list[[s]] <- d
    cat("  ", s, ":", nrow(d), "variants\n")
}
if (length(dt_list) < opt$min_strata) stop("Fewer than ", opt$min_strata, " strata found")
all <- rbindlist(dt_list, fill = TRUE)

# Detect ancestry backgrounds from BETA_ columns
bg <- sub("^BETA_", "", grep("^BETA_", names(all), value = TRUE))
cat("\nAncestry backgrounds:", paste(bg, collapse = ", "), "\n")

ivw <- function(b, se) {
    ok <- !is.na(b) & !is.na(se) & se > 0
    if (sum(ok) == 0) return(list(beta = NA_real_, se = NA_real_, p = NA_real_, i2 = NA_real_, n = 0L))
    if (sum(ok) == 1) return(list(beta = b[ok], se = se[ok], p = 2 * pnorm(-abs(b[ok] / se[ok])), i2 = NA_real_, n = 1L))
    if (has_metafor) {
        m <- tryCatch(metafor::rma(yi = b[ok], sei = se[ok], method = if (opt$method == "random") "REML" else "FE"),
                      error = function(e) NULL)
        if (!is.null(m)) return(list(beta = as.numeric(m$beta), se = m$se, p = m$pval, i2 = m$I2, n = sum(ok)))
    }
    w <- 1 / se[ok]^2; bm <- sum(w * b[ok]) / sum(w); sem <- sqrt(1 / sum(w))
    Q <- sum(w * (b[ok] - bm)^2); df <- sum(ok) - 1
    list(beta = bm, se = sem, p = 2 * pnorm(-abs(bm / sem)), i2 = max(0, (Q - df) / Q * 100), n = sum(ok))
}

fisher_p <- function(p) {
    p <- p[!is.na(p) & p > 0]
    if (length(p) == 0) return(NA_real_)
    pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE)
}

cat("Running meta-analysis on", length(unique(all$ID)), "variants...\n")
pooled <- all[, {
    out <- list(CHR = CHR[1], POS = POS[1], REF = REF[1], ALT = ALT[1],
                N = sum(N, na.rm = TRUE), N_STRATA = .N, STRATA = paste(STRATUM, collapse = ";"),
                P_JOINT = fisher_p(P_JOINT))
    for (a in bg) {
        m <- ivw(get(paste0("BETA_", a)), get(paste0("SE_", a)))
        out[[paste0("MAC_", a)]]  <- if (paste0("MAC_", a) %in% names(.SD)) sum(get(paste0("MAC_", a)), na.rm = TRUE) else NA_real_
        out[[paste0("BETA_", a)]] <- m$beta
        out[[paste0("SE_", a)]]   <- m$se
        out[[paste0("Z_", a)]]    <- if (!is.na(m$se) && m$se > 0) m$beta / m$se else NA_real_
        out[[paste0("P_", a)]]    <- m$p
        out[[paste0("I2_STRATA_", a)]] <- m$i2
        out[[paste0("N_STRATA_", a)]]  <- m$n
    }
    out
}, by = ID]

# Heterogeneity across ancestry BACKGROUNDS in the pooled estimates (Cochran's Q)
pooled[, c("Q_HET", "DF_HET", "P_HET", "I2") := {
    b <- unlist(.SD[, paste0("BETA_", bg), with = FALSE]); s <- unlist(.SD[, paste0("SE_", bg), with = FALSE])
    ok <- !is.na(b) & !is.na(s) & s > 0
    if (sum(ok) >= 2) { w <- 1 / s[ok]^2; bm <- sum(w * b[ok]) / sum(w); Q <- sum(w * (b[ok] - bm)^2); df <- sum(ok) - 1
        list(Q, df, pchisq(Q, df, lower.tail = FALSE), max(0, (Q - df) / Q * 100)) } else list(NA_real_, NA_integer_, NA_real_, NA_real_)
}, by = ID]

setorder(pooled, P_JOINT, na.last = TRUE)
fwrite(pooled, opt$output, sep = "\t", compress = "gzip")
cat("\nPooled results:", opt$output, "\n")
cat("  Variants:", nrow(pooled), "\n")
cat("  P_JOINT < 5e-8:", sum(pooled$P_JOINT < 5e-8, na.rm = TRUE), "| < 1e-5:", sum(pooled$P_JOINT < 1e-5, na.rm = TRUE), "\n")
