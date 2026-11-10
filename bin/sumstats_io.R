# ============================================================================
# sumstats_io.R: "## key=value" lines at the top of a summary-statistics table (CP-08)
# ============================================================================
# Sourced by tractor_genesis_adapter.R and meta_analyze_strata.R (from the same directory).
# A sumstats table may start with lines "## key=value" (no tab in them) before its column header,
# e.g. "## covariates=sex,age_at_diagnosis,PC1,PC2". data.table::fread skips them on its own;
# the combine step keeps them once; sumstats_meta() reads them back.
# ============================================================================
suppressPackageStartupMessages(library(data.table))

# sumstats_meta(file): the leading "##" lines of a table (plain or .gz), as they are
sumstats_meta <- function(file) {
    con <- if (grepl("\\.gz$", file)) gzfile(file, "r") else file(file, "r")
    on.exit(close(con))
    out <- character()
    repeat {
        l <- readLines(con, n = 1)
        if (length(l) == 0 || !startsWith(l, "##")) break
        out <- c(out, l)
    }
    out
}

# sumstats_meta_value(meta, key): the value of "## key=value", or NA
sumstats_meta_value <- function(meta, key) {
    hit <- grep(paste0("^##\\s*", key, "="), meta, value = TRUE)
    if (length(hit) == 0) NA_character_ else sub(paste0("^##\\s*", key, "="), "", hit[1])
}

# fwrite_sumstats(x, file, meta): the "## key=value" lines named in meta (a named character
# vector, c(covariates = "sex,PC1")), then the table; gzip when file ends in .gz
fwrite_sumstats <- function(x, file, meta = character()) {
    lines <- if (length(meta)) paste0("## ", names(meta), "=", gsub("[\t\r\n]+", " ", meta)) else character()
    gz <- grepl("\\.gz$", file)
    if (length(lines)) {
        con <- if (gz) gzfile(file, "w") else file(file, "w")
        writeLines(lines, con); close(con)
    } else if (file.exists(file)) {
        file.remove(file)
    }
    # appended as a second gzip member when compressed; zcat and fread read both as one stream
    fwrite(x, file, sep = "\t", append = length(lines) > 0, col.names = TRUE,
           compress = if (gz) "gzip" else "none")
    invisible(file)
}
