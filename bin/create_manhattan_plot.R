#!/usr/bin/env Rscript
# Manhattan plot from ANY GWAS layout (Tractor-GENESIS, REGENIE, SAIGE, GENESIS, PLINK2, ...)
# Column handling lives in sumstats_utils.R.

suppressPackageStartupMessages({ library(data.table); library(argparse) })

.this_dir <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) ".")
.helper <- Sys.which("sumstats_utils.R"); if (!nzchar(.helper)) .helper <- file.path(.this_dir, "sumstats_utils.R")
source(.helper)

parser <- ArgumentParser(description = "Manhattan plot")
parser$add_argument("--sumstats", required = TRUE)
parser$add_argument("--lambda-gc", required = FALSE, help = "Lambda GC file from LAMBDA_GC (optional)")
parser$add_argument("--output-prefix", required = TRUE)
parser$add_argument("--title", default = "Manhattan plot")
parser$add_argument("--genome-wide", type = "double", default = 5e-8)
parser$add_argument("--suggestive", type = "double", default = 1e-5)
parser$add_argument("--p-column", default = NULL, help = "Override: e.g. P_AMR to plot an ancestry-specific Tractor test")
parser$add_argument("--per-ancestry", action = "store_true", default = FALSE,
                    help = "Tractor-GENESIS: also draw one panel per ancestry-specific P_<anc>")
parser$add_argument("--label-top", type = "integer", default = 10, help = "Label this many lead variants")
args <- parser$parse_args()

cat("Reading", args$sumstats, "\n")
ss <- read_sumstats_std(args$sumstats)
fmt <- attr(ss, "source_format")
cat("Detected format:", fmt, "|", nrow(ss), "variants\n")

lambda <- read_lambda(args$lambda_gc)
title_main <- if (!is.na(lambda)) sprintf("%s  (lambda GC = %.3f)", args$title, lambda) else args$title

# ---- plotting (base graphics, no qqman dependency) -------------------------
manhattan_base <- function(d, pcol = "P", main = "", gw = 5e-8, sug = 1e-5, label_top = 10) {
    d <- d[!is.na(get(pcol)) & get(pcol) > 0]
    d[, logp := -log10(get(pcol))]
    setorder(d, CHR, POS)
    chr_len <- d[, .(len = max(POS)), by = CHR][order(CHR)]
    chr_len[, offset := cumsum(as.numeric(shift(len, fill = 0)))]
    d <- merge(d, chr_len[, .(CHR, offset)], by = "CHR")
    d[, x := POS + offset]
    centers <- d[, .(x = mean(range(x))), by = CHR]
    cols <- ifelse(d$CHR %% 2 == 1, "#1565C0", "#90CAF9")
    ymax <- max(c(d$logp, -log10(gw)), na.rm = TRUE) * 1.08
    plot(d$x, d$logp, pch = 20, cex = 0.5, col = cols, xaxt = "n", xlab = "Chromosome",
         ylab = expression(-log[10](P)), main = main, ylim = c(0, ymax), bty = "l")
    lab <- ifelse(centers$CHR == 23, "X", ifelse(centers$CHR == 24, "Y", ifelse(centers$CHR == 25, "MT", centers$CHR)))
    axis(1, at = centers$x, labels = lab, cex.axis = 0.8, las = 2)
    abline(h = -log10(gw), col = "red", lty = 2); abline(h = -log10(sug), col = "blue", lty = 3)
    if (label_top > 0) {
        leads <- loci_from_sumstats(d[, .(CHR, POS, ID, P = get(pcol))], p_threshold = gw, window_kb = 500,
                                    fallback_threshold = sug, max_loci = label_top)
        if (nrow(leads)) {
            leads <- merge(leads, d[, .(CHR, POS, x, logp)], by = c("CHR", "POS"))
            points(leads$x, leads$logp, pch = 18, col = "#D32F2F", cex = 1.2)
            text(leads$x, leads$logp, labels = leads$ID, pos = 3, cex = 0.6, offset = 0.4)
        }
    }
}

draw_all <- function() {
    panels <- list(list(col = "P", title = title_main))
    if (!is.null(args$p_column)) {
        if (!args$p_column %in% names(ss)) stop("--p-column ", args$p_column, " not in file")
        panels <- list(list(col = args$p_column, title = paste(args$title, "-", args$p_column)))
    } else if (args$per_ancestry && fmt %in% c("tractor_genesis", "pooled_ivw")) {
        for (a in ancestries_in_sumstats(ss)) {
            pc <- paste0("P_", a)
            if (pc %in% names(ss)) panels[[length(panels) + 1]] <- list(col = pc, title = sprintf("%s - %s background (Tractor)", args$title, a))
        }
    }
    if (length(panels) > 1) par(mfrow = c(length(panels), 1), mar = c(4, 4.5, 3, 1))
    for (p in panels) manhattan_base(copy(ss), p$col, p$title, args$genome_wide, args$suggestive, args$label_top)
}

h <- 600 * (if (args$per_ancestry && fmt %in% c("tractor_genesis", "pooled_ivw")) 1 + length(ancestries_in_sumstats(ss)) else 1)
png(paste0(args$output_prefix, ".manhattan.png"), width = 1400, height = h, res = 110); draw_all(); dev.off()
pdf(paste0(args$output_prefix, ".manhattan.pdf"), width = 13, height = h / 110); draw_all(); dev.off()
cat("Manhattan plot written:", paste0(args$output_prefix, ".manhattan.{png,pdf}"), "\n")
