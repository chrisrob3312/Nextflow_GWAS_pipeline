#!/usr/bin/env Rscript
# Regional association plots (LocusZoom-style, no LD server needed) from ANY GWAS layout.
# Loci come from a GWAS_FILTER "significant" file when given, else from the sumstats
# themselves. Optional gene track from a GTF/GFF (gene features) or a BED (chr start end name).
# Tractor-GENESIS: --per-ancestry adds one row per ancestry-specific P_<anc>.

suppressPackageStartupMessages({ library(data.table); library(argparse) })

.this_dir <- tryCatch(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) ".")
.helper <- Sys.which("sumstats_utils.R"); if (!nzchar(.helper)) .helper <- file.path(.this_dir, "sumstats_utils.R")
source(.helper)

parser <- ArgumentParser(description = "Regional association plots")
parser$add_argument("--sumstats", required = TRUE, help = "Full GWAS sumstats (any layout)")
parser$add_argument("--loci", default = NULL, help = "Significant-variant file or loci table (any layout); default: derive from sumstats")
parser$add_argument("--genes", default = NULL, help = "GTF/GFF or BED for the gene track (optional)")
parser$add_argument("--output-dir", required = TRUE)
parser$add_argument("--ancestry", default = "ALL")
parser$add_argument("--window-kb", type = "integer", default = 500)
parser$add_argument("--p-threshold", type = "double", default = 5e-8)
parser$add_argument("--suggestive", type = "double", default = 1e-5)
parser$add_argument("--max-loci", type = "integer", default = 20)
parser$add_argument("--per-ancestry", action = "store_true", default = FALSE)
args <- parser$parse_args()
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

ss <- read_sumstats_std(args$sumstats)
fmt <- attr(ss, "source_format")
cat("Detected format:", fmt, "|", nrow(ss), "variants\n")

# ---- loci ------------------------------------------------------------------
loci <- NULL
if (!is.null(args$loci) && file.exists(args$loci) && file.size(args$loci) > 0) {
    sig <- tryCatch(read_sumstats_std(args$loci), error = function(e) NULL)
    if (!is.null(sig) && nrow(sig) > 0) loci <- loci_from_sumstats(sig, args$p_threshold, args$window_kb, args$suggestive, args$max_loci)
}
if (is.null(loci) || nrow(loci) == 0) loci <- loci_from_sumstats(ss, args$p_threshold, args$window_kb, args$suggestive, args$max_loci)
if (nrow(loci) == 0) { cat("No loci below", args$suggestive, "- nothing to plot\n"); quit(save = "no", status = 0) }
cat("Loci:", nrow(loci), "\n")

# ---- optional gene track ---------------------------------------------------
genes <- NULL
if (!is.null(args$genes) && file.exists(args$genes)) {
    if (grepl("\\.(gtf|gff3?)(\\.gz)?$", args$genes, ignore.case = TRUE)) {
        g <- fread(args$genes, sep = "\t", header = FALSE, quote = "", fill = TRUE, skip = "#")
        g <- g[V3 == "gene"]
        name <- sub('.*gene_name[ =]"?([^";]+)"?.*', "\\1", g$V9)
        name[name == g$V9] <- sub('.*Name=([^;]+).*', "\\1", g$V9[name == g$V9])
        genes <- data.table(CHR = chr_to_num(g$V1), start = as.integer(g$V4), end = as.integer(g$V5), name = name, strand = g$V7)
    } else {
        g <- fread(args$genes, header = FALSE)
        genes <- data.table(CHR = chr_to_num(g$V1), start = as.integer(g$V2), end = as.integer(g$V3),
                            name = if (ncol(g) >= 4) as.character(g$V4) else "", strand = if (ncol(g) >= 6) as.character(g$V6) else "+")
    }
    genes <- genes[!is.na(CHR)]
}

# ---- plotting --------------------------------------------------------------
pcols <- "P"
if (args$per_ancestry && fmt %in% c("tractor_genesis", "pooled_ivw")) {
    pcols <- c("P", intersect(paste0("P_", ancestries_in_sumstats(ss)), names(ss)))
}

plot_locus <- function(reg, pcol, title_suffix = "") {
    r <- ss[CHR == reg$CHR & POS >= reg$start & POS <= reg$end & !is.na(get(pcol)) & get(pcol) > 0]
    if (nrow(r) == 0) return(invisible(FALSE))
    r[, logp := -log10(get(pcol))]
    lead <- r[which.min(get(pcol))]
    has_genes <- !is.null(genes) && nrow(genes[CHR == reg$CHR & end >= reg$start & start <= reg$end]) > 0
    layout(matrix(1:2, ncol = 1), heights = if (has_genes) c(3, 1.2) else c(1, 0.001))
    par(mar = c(if (has_genes) 0.5 else 4, 4.5, 3, 1))
    chr_lab <- ifelse(reg$CHR == 23, "X", ifelse(reg$CHR == 24, "Y", reg$CHR))
    ymax <- max(c(r$logp, -log10(args$p_threshold))) * 1.1
    plot(r$POS / 1e6, r$logp, pch = 20, cex = 0.7, col = "#607D8B", xaxt = if (has_genes) "n" else "s",
         xlab = if (has_genes) "" else sprintf("Chromosome %s position (Mb)", chr_lab), ylab = expression(-log[10](P)),
         main = sprintf("%s  chr%s:%s  lead %s (P = %.2g)%s", args$ancestry, chr_lab, format(lead$POS, big.mark = ","), lead$ID, lead[[pcol]], title_suffix),
         ylim = c(0, ymax), bty = "l", cex.main = 0.9)
    abline(h = -log10(args$p_threshold), col = "red", lty = 2); abline(h = -log10(args$suggestive), col = "blue", lty = 3)
    points(lead$POS / 1e6, lead$logp, pch = 18, cex = 1.8, col = "#7B1FA2")
    text(lead$POS / 1e6, lead$logp, lead$ID, pos = 3, cex = 0.75)
    if (has_genes) {
        gg <- genes[CHR == reg$CHR & end >= reg$start & start <= reg$end][order(start)]
        gg[, row := (seq_len(.N) - 1) %% 4]
        par(mar = c(4, 4.5, 0.5, 1))
        plot(NA, xlim = range(r$POS / 1e6), ylim = c(-0.5, 3.8), yaxt = "n", ylab = "", xlab = sprintf("Chromosome %s position (Mb)", chr_lab), bty = "n")
        segments(gg$start / 1e6, gg$row, gg$end / 1e6, gg$row, lwd = 4, col = "#455A64")
        text((pmax(gg$start, reg$start) + pmin(gg$end, reg$end)) / 2e6, gg$row + 0.35, gg$name, cex = 0.6)
    }
    invisible(TRUE)
}

n_written <- 0
for (i in seq_len(nrow(loci))) {
    reg <- loci[i]
    for (pcol in pcols) {
        tag <- if (pcol == "P") "" else paste0(".", sub("^P_", "", pcol))
        base <- file.path(args$output_dir, sprintf("%s.chr%s_%d%s", args$ancestry, reg$CHR, reg$POS, tag))
        suffix <- if (pcol == "P") "" else sprintf("  [%s background]", sub("^P_", "", pcol))
        png(paste0(base, ".regional.png"), width = 1100, height = 750, res = 110); ok <- plot_locus(reg, pcol, suffix); dev.off()
        pdf(paste0(base, ".regional.pdf"), width = 10, height = 7); plot_locus(reg, pcol, suffix); dev.off()
        if (isTRUE(ok)) n_written <- n_written + 1
    }
}
fwrite(loci, file.path(args$output_dir, paste0(args$ancestry, ".loci.tsv")), sep = "\t")
cat("Regional plots written:", n_written, "\n")
