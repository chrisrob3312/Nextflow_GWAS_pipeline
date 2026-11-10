#!/usr/bin/env Rscript
#' Create Manhattan plot from GWAS summary statistics

suppressPackageStartupMessages({
  library(qqman)
  library(data.table)
  library(ggplot2)
  library(argparse)
})

# Parse arguments
parser <- ArgumentParser(description = "Create Manhattan plot")
parser$add_argument("--sumstats", required = TRUE, help = "Summary statistics file")
parser$add_argument("--lambda-gc", required = FALSE, help = "Lambda GC file from elsewhere: printed for comparison, never used in place of the lambda computed here")
parser$add_argument("--output-prefix", required = TRUE, help = "Output prefix")
parser$add_argument("--title", default = "Manhattan Plot", help = "Plot title")
parser$add_argument("--genome-wide", type = "double", default = 5e-8, help = "Genome-wide significance threshold")
parser$add_argument("--suggestive", type = "double", default = 1e-5, help = "Suggestive threshold")
parser$add_argument("--p-col", required = FALSE, default = NULL, help = "P column to plot [default: P_JOINT, else P / PVALUE / P_VALUE / p.value]")
args <- parser$parse_args()

# resolve_p_col(nm, requested): the P column to use. --p-col names it (exactly, else ignoring case);
# without it P_JOINT (the Tractor joint test) comes first, then the usual single-P names (CP-14)
resolve_p_col <- function(nm, requested = NULL) {
  if (!is.null(requested) && nzchar(requested)) {
    hit <- c(nm[nm == requested], nm[toupper(nm) == toupper(requested)])
    if (length(hit) == 0) stop("--p-col ", requested, " is not a column (P columns: ",
                               paste(grep("^P(_|$)|PVAL", nm, value = TRUE, ignore.case = TRUE), collapse = ", "), ")")
    return(hit[1])
  }
  for (pat in c("^P_JOINT$", "^(P|PVALUE|P_VALUE|p.value)$")) {
    hit <- grep(pat, nm, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  stop("Could not detect the P-value column (looked for P_JOINT, P, PVALUE, P_VALUE, p.value; or pass --p-col)")
}

# Read data (leading "## key=value" lines are skipped by fread)
cat("Reading summary statistics...\n")
if (grepl("\\.gz$", args$sumstats)) {
  df <- fread(cmd = paste("zcat", shQuote(args$sumstats)))
} else {
  df <- fread(args$sumstats)
}

# Detect columns
chr_col <- names(df)[grep("^(CHR|CHROM|#CHROM|chromosome)$", names(df), ignore.case = TRUE)[1]]
pos_col <- names(df)[grep("^(POS|BP|GENPOS|position)$", names(df), ignore.case = TRUE)[1]]
snp_col <- names(df)[grep("^(SNP|ID|RSID|variant_id|MarkerName)$", names(df), ignore.case = TRUE)[1]]
p_col <- resolve_p_col(names(df), args$p_col)
cat("P column:", p_col, "\n")

if (is.na(chr_col) || is.na(pos_col)) {
  stop("Could not detect required columns (CHR / CHROM, POS / BP)")
}

# Prepare data for manhattan plot
plot_df <- data.frame(
  CHR = df[[chr_col]],
  BP = df[[pos_col]],
  P = df[[p_col]],
  SNP = if (!is.na(snp_col)) df[[snp_col]] else paste0("chr", df[[chr_col]], ":", df[[pos_col]])
)

# Remove invalid entries
plot_df <- plot_df[!is.na(plot_df$P) & plot_df$P > 0 & plot_df$P <= 1, ]
plot_df <- plot_df[!is.na(plot_df$CHR) & !is.na(plot_df$BP), ]

# Convert chromosome to numeric
plot_df$CHR <- as.numeric(gsub("chr", "", plot_df$CHR, ignore.case = TRUE))
plot_df <- plot_df[!is.na(plot_df$CHR), ]

# Lambda GC of the plotted column (as in create_qq_plot.R); a file passed in is printed, not used
if (nrow(plot_df) == 0) stop("No P values in (0, 1] in column ", p_col)
lambda_gc <- median(qchisq(plot_df$P, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
if (!is.null(args$lambda_gc) && file.exists(args$lambda_gc)) {
  cat("Lambda GC passed in (", args$lambda_gc, "; not used): ", readLines(args$lambda_gc, n = 1), "\n", sep = "")
}

# Create title with the column and its lambda
plot_title <- paste0(args$title, " (", p_col, ", λGC = ", round(lambda_gc, 3), ")")

# Create PNG
cat("Creating Manhattan plot...\n")
png(paste0(args$output_prefix, ".manhattan.png"), width = 1200, height = 600, res = 100)
manhattan(
  plot_df,
  main = plot_title,
  genomewideline = -log10(args$genome_wide),
  suggestiveline = -log10(args$suggestive),
  col = c("#1565C0", "#90CAF9"),
  cex = 0.6,
  cex.axis = 0.9
)
dev.off()

# Create PDF
pdf(paste0(args$output_prefix, ".manhattan.pdf"), width = 12, height = 6)
manhattan(
  plot_df,
  main = plot_title,
  genomewideline = -log10(args$genome_wide),
  suggestiveline = -log10(args$suggestive),
  col = c("#1565C0", "#90CAF9"),
  cex = 0.6,
  cex.axis = 0.9
)
dev.off()

cat("Manhattan plot created successfully\n")
