#!/usr/bin/env Rscript
# Load every R package named on the command line with library(). Prints one line per package,
# "ok <package> <version>" or "MISSING <package> <error>", then the R version, library paths and
# tempdir(), and exits 1 if any package failed. bin/check_envs.sh runs it inside the gwas-r env:
#   Rscript bin/check_r_libs.R GENESIS SeqArray coxme ...
pkgs <- commandArgs(trailingOnly = TRUE)
if (length(pkgs) == 0) stop("usage: check_r_libs.R <package> [<package> ...]")

failed <- character(0)
for (p in pkgs) {
    err <- tryCatch({
        suppressPackageStartupMessages(library(p, character.only = TRUE))
        NULL
    }, error = function(e) conditionMessage(e))
    if (is.null(err)) {
        cat(sprintf("ok      %-12s %s\n", p, as.character(utils::packageVersion(p))))
    } else {
        cat(sprintf("MISSING %-12s %s\n", p, gsub("[[:space:]]+", " ", err)))
        failed <- c(failed, p)
    }
}
cat(sprintf("R %s; library paths: %s\n", getRversion(), paste(.libPaths(), collapse = ", ")))
cat(sprintf("tempdir(): %s\n", tempdir()))
quit(save = "no", status = if (length(failed) > 0) 1L else 0L)
