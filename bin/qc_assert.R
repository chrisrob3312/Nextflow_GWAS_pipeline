# ============================================================================
# qc_assert.R: file a QC verdict, and count samples at every join (CP-11)
# ============================================================================
# Sourced by the R steps:   source(file.path(<bin dir>, "qc_assert.R"))
#   qc_assert(id, pass, n, message)   one verdict row: id, verdict (PASS | FAIL | INFO), n, message
#   qc_account(step, n_in, n_out, reason, stratum, trait, chromosome)   one sample-accounting row
# or run from a shell:
#   Rscript bin/qc_assert.R --id <id> --verdict PASS|FAIL|INFO [--n N] [--message "..."]
#
# Rows go to GWAS_QC_DIR (the launcher sets <run dir>/qc; default results/qc), one small file per
# row under assertions.d/ and accounting.d/, so concurrent array tasks never write the same file.
# bin/qc_gate.sh reads them; the run record's finalise job merges them into qc_assertions.tsv
# (header id verdict n message, as in the handoff bundle) and sample_accounting.tsv.
# A gate treats a MISSING assertion exactly like a failed one: "did not run" is not "passed".
# ============================================================================

qc_dir <- function() Sys.getenv("GWAS_QC_DIR", file.path(getwd(), "results", "qc"))

.qc_tag <- local({
    k <- 0L
    function() {
        k <<- k + 1L
        job <- Sys.getenv("SLURM_ARRAY_JOB_ID")
        job <- if (nzchar(job)) paste0(job, "_", Sys.getenv("SLURM_ARRAY_TASK_ID")) else Sys.getenv("SLURM_JOB_ID", "local")
        paste0(job, ".", Sys.getpid(), ".", k)
    }
})

.qc_clean <- function(x) gsub("[\t\r\n]+", " ", if (is.null(x) || length(x) == 0 || is.na(x[1])) "NA" else as.character(x[1]))

.qc_write <- function(sub, name, fields) {
    d <- file.path(qc_dir(), sub)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    f <- file.path(d, paste0(gsub("[^A-Za-z0-9_.-]", "_", name), ".", .qc_tag(), ".tsv"))
    tmp <- paste0(f, ".tmp")
    writeLines(paste(vapply(fields, .qc_clean, ""), collapse = "\t"), tmp)
    file.rename(tmp, f)
    invisible(f)
}

# qc_copy_to(path): also keep this step's verdicts in <path> (header id verdict n message), next to
# the step's outputs. The structure step does this, so a launch that reuses its outputs
# (SKIP_PCAIR=true) can carry the verdicts into its own run record instead of losing them.
qc_copy_to <- function(path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines("id\tverdict\tn\tmessage", path)
    options(qc_copy_to = path)
}

# qc_assert(id, pass, ...): pass = TRUE / FALSE, or verdict = "INFO" for a note that is neither
qc_assert <- function(id, pass, n = NA, message = "", verdict = NULL) {
    v <- verdict %||% (if (isTRUE(pass)) "PASS" else "FAIL")
    if (!v %in% c("PASS", "FAIL", "INFO")) stop("qc_assert: verdict must be PASS, FAIL or INFO")
    fields <- list(id, v, n, message)
    .qc_write("assertions.d", id, fields)
    if (!is.null(getOption("qc_copy_to"))) {
        cat(paste(vapply(fields, .qc_clean, ""), collapse = "\t"), "\n", sep = "", file = getOption("qc_copy_to"), append = TRUE)
    }
    cat(sprintf("QC %s %s (n=%s) %s\n", v, id, .qc_clean(n), message))
    invisible(v == "PASS" || v == "INFO")
}

# qc_account(step, n_in, n_out, ...): one join's sample counts
qc_account <- function(step, n_in, n_out, reason = "", stratum = NA, trait = NA, chromosome = NA) {
    .qc_write("accounting.d", step, list(step, stratum, trait, chromosome, n_in, n_out, n_in - n_out, reason))
    cat(sprintf("samples %s: %s -> %s (%s)\n", step, n_in, n_out, reason))
    invisible(NULL)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# merge_pcs_checked(ph, pcs, ...): phenotype rows + PC columns, with the pcair_pc_merge verdict.
# Every phenotype row must get its PCs. A row without them used to be carried into the GWAS with
# NA PCs and dropped there without a word; now it is a FAIL row that the next step's gate stops
# on. ph and pcs are data.tables; pcs is keyed by sample_id. Returns the merged table, phenotype
# columns first (FID / IID stay in front, which REGENIE, SAIGE and SPA-Cox need).
merge_pcs_checked <- function(ph, pcs, id_col = "sample_id", stratum = NA, what = "phenotype",
                              assertion = "pcair_pc_merge") {
    pc_cols <- setdiff(names(pcs), "sample_id")
    out <- merge(ph, pcs, by.x = id_col, by.y = "sample_id", all.x = TRUE, sort = FALSE)
    data.table::setcolorder(out, c(names(ph), pc_cols))
    lacking <- out[[id_col]][!stats::complete.cases(as.data.frame(out)[, pc_cols, drop = FALSE])]
    n_unused <- sum(!pcs$sample_id %in% ph[[id_col]])
    qc_account("pc_merge", nrow(ph), nrow(ph) - length(lacking), stratum = stratum,
               reason = sprintf("%s rows with all of %s (%d PC rows have no %s row)", what,
                                paste(range(pc_cols), collapse = ".."), n_unused, what))
    qc_assert(assertion, length(lacking) == 0, n = nrow(ph),
              message = if (length(lacking) == 0) {
                  sprintf("%s%s: every %s row has PCs", what, if (is.na(stratum)) "" else paste0(" (", stratum, ")"), what)
              } else {
                  sprintf("%s%s: %d of %d rows have no PCs (not in the genotype set?), e.g. %s", what,
                          if (is.na(stratum)) "" else paste0(" (", stratum, ")"), length(lacking), nrow(ph),
                          paste(utils::head(lacking, 3), collapse = ", "))
              })
    out
}

# --- command line --------------------------------------------------------------------------
if (any(grepl("^--file=.*qc_assert\\.R$", commandArgs(FALSE)))) {
    a <- commandArgs(trailingOnly = TRUE)
    val <- function(flag, default = NULL) { i <- match(flag, a); if (is.na(i) || i == length(a)) default else a[i + 1] }
    id <- val("--id"); v <- val("--verdict")
    if (is.null(id) || is.null(v)) stop("usage: qc_assert.R --id <id> --verdict PASS|FAIL|INFO [--n N] [--message text]")
    qc_assert(id, pass = v == "PASS", n = val("--n", NA), message = val("--message", ""), verdict = v)
}
