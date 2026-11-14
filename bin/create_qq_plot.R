#!/usr/bin/env Rscript
# QQ plot from ANY GWAS layout; column handling lives in sumstats_utils.R.
# For Tractor-GENESIS files, --per-ancestry overlays the per-ancestry P_<anc> tests.

suppressPackageStartupMessages({ library(data.table); library(argparse) })

.this_dir <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) ".")
.helper <- Sys.which("sumstats_utils.R"); if (!nzchar(.helper)) .helper <- file.path(.this_dir, "sumstats_utils.R")
source(.helper)

parser <- ArgumentParser(description = "QQ plot")
parser$add_argument("--sumstats", required = TRUE)
parser$add_argument("--lambda-gc", required = FALSE)
parser$add_argument("--output-prefix", required = TRUE)
parser$add_argument("--title", default = "QQ plot")
parser$add_argument("--per-ancestry", action = "store_true", default = FALSE)
parser$add_argument("--max-points", type = "integer", default = 200000, help = "Thin the dense low end for speed")
args <- parser$parse_args()

ss <- read_sumstats_std(args$sumstats)
fmt <- attr(ss, "source_format")
cat("Detected format:", fmt, "|", nrow(ss), "variants\n")

lambda_of <- function(p) { p <- p[!is.na(p) & p > 0 & p <= 1]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
qq_points <- function(p, max_points) {
    p <- sort(p[!is.na(p) & p > 0 & p <= 1]); n <- length(p)
    obs <- -log10(p); exp <- -log10(ppoints(n))
    lo <- -log10(qbeta(0.975, 1:n, n:1)); hi <- -log10(qbeta(0.025, 1:n, n:1))
    keep <- if (n > max_points) sort(unique(c(1:min(n, 20000), round(seq(1, n, length.out = max_points))))) else seq_len(n)
    data.table(exp = exp[keep], obs = obs[keep], lo = lo[keep], hi = hi[keep])
}

series <- list(list(name = "P", p = ss$P, col = "#1565C0"))
if (args$per_ancestry && fmt %in% c("tractor_genesis", "pooled_ivw")) {
    pal <- c("#D32F2F", "#2E7D32", "#F9A825", "#6A1B9A", "#00838F")
    for (i in seq_along(ancestries_in_sumstats(ss))) {
        a <- ancestries_in_sumstats(ss)[i]; pc <- paste0("P_", a)
        if (pc %in% names(ss)) series[[length(series) + 1]] <- list(name = paste0(a, " background"), p = ss[[pc]], col = pal[(i - 1) %% length(pal) + 1])
    }
}

lambda <- read_lambda(args$lambda_gc); if (is.na(lambda)) lambda <- lambda_of(ss$P)
main <- sprintf("%s  (lambda GC = %.3f)", args$title, lambda)

draw <- function() {
    pts <- lapply(series, function(s) qq_points(s$p, args$max_points))
    xmax <- max(sapply(pts, function(q) max(q$exp))); ymax <- max(sapply(pts, function(q) max(q$obs)))
    plot(NA, xlim = c(0, xmax), ylim = c(0, ymax * 1.05), xlab = expression(Expected ~ -log[10](P)),
         ylab = expression(Observed ~ -log[10](P)), main = main, bty = "l")
    q0 <- pts[[1]]; polygon(c(q0$exp, rev(q0$exp)), c(q0$lo, rev(q0$hi)), col = "grey88", border = NA)
    abline(0, 1, col = "red", lty = 2)
    for (i in seq_along(series)) points(pts[[i]]$exp, pts[[i]]$obs, pch = 20, cex = 0.6, col = series[[i]]$col)
    if (length(series) > 1) {
        legend("topleft", bty = "n", pch = 20, col = sapply(series, `[[`, "col"),
               legend = sprintf("%s (lambda %.3f)", sapply(series, `[[`, "name"), sapply(series, function(s) lambda_of(s$p))))
    }
}

png(paste0(args$output_prefix, ".qq.png"), width = 800, height = 800, res = 110); draw(); dev.off()
pdf(paste0(args$output_prefix, ".qq.pdf"), width = 7, height = 7); draw(); dev.off()
cat("QQ plot written; lambda GC =", round(lambda, 4), "\n")
