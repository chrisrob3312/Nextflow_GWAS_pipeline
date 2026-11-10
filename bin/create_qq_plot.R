#!/usr/bin/env Rscript
#' Create QQ plot from GWAS summary statistics

suppressPackageStartupMessages({
  library(qqman)
  library(data.table)
  library(ggplot2)
  library(argparse)
})

# Parse arguments
parser <- ArgumentParser(description = "Create QQ plot")
parser$add_argument("--sumstats", required = TRUE, help = "Summary statistics file")
parser$add_argument("--lambda-gc", required = FALSE, help = "Lambda GC file from elsewhere: printed for comparison, never used in place of the lambda computed here")
parser$add_argument("--p-col", required = FALSE, default = NULL, help = "P column to plot [default: P_JOINT, else P / PVALUE / P_VALUE / p.value]")
parser$add_argument("--output-prefix", required = TRUE, help = "Output prefix")
parser$add_argument("--title", default = "QQ Plot", help = "Plot title")
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

p_col <- resolve_p_col(names(df), args$p_col)
cat("P column:", p_col, "\n")

# Extract P-values
p_values <- df[[p_col]]
p_values <- p_values[!is.na(p_values) & p_values > 0 & p_values <= 1]
if (length(p_values) == 0) stop("No P values in (0, 1] in column ", p_col)

# Lambda GC of the plotted column: the median of the P values as 1-df chi-square, over its null
# median. Always this one: a lambda passed in (--lambda-gc) is printed next to it, never put in
# its place (it used to overwrite this one, so the plot could show another column's lambda).
lambda_gc <- median(qchisq(p_values, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
if (!is.finite(lambda_gc)) stop("lambda GC is not finite for column ", p_col)
writeLines(format(lambda_gc, digits = 6), paste0(args$output_prefix, ".lambda.txt"))
cat("Lambda GC (", p_col, ", ", length(p_values), " P values): ", round(lambda_gc, 4), " -> ", args$output_prefix, ".lambda.txt\n", sep = "")
if (!is.null(args$lambda_gc) && file.exists(args$lambda_gc)) {
  cat("Lambda GC passed in (", args$lambda_gc, "; not used): ", readLines(args$lambda_gc, n = 1), "\n", sep = "")
}

# Create title with lambda
plot_title <- paste0(args$title, " (", p_col, ", λGC = ", round(lambda_gc, 3), ")")

# Create PNG
cat("Creating QQ plot...\n")
png(paste0(args$output_prefix, ".qq.png"), width = 800, height = 800, res = 100)
qq(p_values, main = plot_title)
dev.off()

# Create PDF
pdf(paste0(args$output_prefix, ".qq.pdf"), width = 8, height = 8)
qq(p_values, main = plot_title)
dev.off()

# Also create a ggplot2 version with confidence intervals
observed <- -log10(sort(p_values))
expected <- -log10(ppoints(length(p_values)))

qq_df <- data.frame(
  expected = expected,
  observed = observed
)

# Confidence intervals
n <- length(p_values)
qq_df$lower <- -log10(qbeta(0.975, 1:n, n:1))
qq_df$upper <- -log10(qbeta(0.025, 1:n, n:1))

p <- ggplot(qq_df, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey80", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_point(size = 0.8, alpha = 0.6, color = "#1565C0") +
  labs(
    title = plot_title,
    x = expression(Expected ~ -log[10](P)),
    y = expression(Observed ~ -log[10](P))
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(paste0(args$output_prefix, ".qq_ggplot.png"), p, width = 8, height = 8, dpi = 100)
ggsave(paste0(args$output_prefix, ".qq_ggplot.pdf"), p, width = 8, height = 8)

cat("QQ plots created successfully\n")
cat("Lambda GC:", round(lambda_gc, 4), "\n")
