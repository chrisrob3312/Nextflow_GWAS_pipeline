# ============================================================================
# tractor_io.R: shared helpers for reading Tractor files and gating a stratum
# ============================================================================
# Sourced by tractor_genesis_adapter.R and calculate_prs_admixed.R (from the same directory).
# Needs data.table. No analysis happens here; every function either returns a value or stops
# with a message that says what to fix.
# ============================================================================

suppressPackageStartupMessages(library(data.table))

# ---------------------------------------------------------------------------
# Which file holds which ancestry
# ---------------------------------------------------------------------------
# Files named by ancestry LABEL are used as they are. Files named by a NUMBER (FLARE's anc0..ancK-1,
# Tractor's extract_tracts*.py output, the Nextflow extraction module) are used only through an
# index map that says which number is which label: the handoff bundle's
# local_ancestry/flare_index_map.tsv (header anc_index<TAB>label). The position of an ancestry in
# --ancestries is never taken to be its number: that list is alphabetical, FLARE numbers by the
# reference panel's order, and the two need not agree.
TRACTOR_LABEL_PATTERNS <- list(
    hapcount = c("%s.hapcount.%s.txt.gz", "%s.%s.hapcount.txt.gz", "%s.hapcount.%s.tsv.gz"),
    dosage   = c("%s.dosage.%s.txt.gz",   "%s.%s.dosage.txt.gz",   "%s.ancdose.%s.tsv.gz")
)
TRACTOR_INDEX_PATTERNS <- list(
    hapcount = c("%s.anc%s.hapcount.txt.gz", "%s.hapcount.%s.txt.gz", "%s.hapcount.%s.tsv.gz"),
    dosage   = c("%s.anc%s.dosage.txt.gz",   "%s.dosage.%s.txt.gz",   "%s.ancdose.%s.tsv.gz")
)

# read_index_map(path): named character vector, label -> index ("EUR" = "2", ...)
read_index_map <- function(path) {
    if (!file.exists(path)) stop("ancestry index map ", path, " does not exist")
    m <- fread(path, colClasses = "character")
    if (!all(c("anc_index", "label") %in% names(m))) {
        stop(path, ": header must be anc_index<TAB>label (the handoff's local_ancestry/flare_index_map.tsv)")
    }
    if (anyDuplicated(m$label) || anyDuplicated(m$anc_index)) stop(path, ": a label or an index appears twice")
    setNames(m$anc_index, m$label)
}

# resolve_tractor_file(prefix, anc, kind, index_map): list(path, how) for kind "hapcount" or "dosage"
resolve_tractor_file <- function(prefix, anc, kind, index_map = NULL) {
    by_label <- sprintf(TRACTOR_LABEL_PATTERNS[[kind]], prefix, anc)
    hit <- by_label[file.exists(by_label)]
    if (length(hit) > 0) return(list(path = hit[1], how = "named by label"))

    if (!is.null(index_map)) {
        if (!anc %in% names(index_map)) {
            stop("ancestry '", anc, "' is not in the index map (it has: ",
                 paste(names(index_map), collapse = ", "), ")")
        }
        k <- index_map[[anc]]
        by_index <- sprintf(TRACTOR_INDEX_PATTERNS[[kind]], prefix, k)
        hit <- by_index[file.exists(by_index)]
        if (length(hit) > 0) return(list(path = hit[1], how = paste0("index ", k, " per the index map")))
        stop("no ", kind, " file for ancestry '", anc, "': tried by label ",
             paste(basename(by_label), collapse = ", "), "; by index ", k, " ",
             paste(basename(by_index), collapse = ", "))
    }

    numbered <- unlist(lapply(0:9, function(k) {
        f <- sprintf(TRACTOR_INDEX_PATTERNS[[kind]], prefix, k); f[file.exists(f)]
    }))
    msg <- paste0("no ", kind, " file named by ancestry label for '", anc, "' (tried ",
                  paste(basename(by_label), collapse = ", "), ")")
    if (length(numbered) > 0) {
        msg <- paste0(msg, ". Files numbered by position exist (", paste(basename(numbered), collapse = ", "),
                      "), but nothing says which number is which ancestry: rename them by label, or pass ",
                      "--ancestry_index_map (the handoff's local_ancestry/flare_index_map.tsv)")
    }
    stop(msg)
}

# Variant-information columns of a Tractor table; everything else is a sample
TRACTOR_INFO_COLS <- c("CHR", "CHROM", "POS", "ID", "SNP", "REF", "ALT",
                       "chrom", "pos", "id", "snp", "ref", "alt", "#CHROM", "BP")
tractor_info_cols <- function(nm) {
    info <- intersect(nm, TRACTOR_INFO_COLS)
    if (length(info) == 0) info <- nm[seq_len(min(5, length(nm)))]   # assume the first 5 columns
    info
}

# tractor_sample_ids(path): the sample columns of a Tractor table, from its header only
tractor_sample_ids <- function(path) {
    nm <- names(fread(path, nrows = 0))
    setdiff(nm, tractor_info_cols(nm))
}

# ---------------------------------------------------------------------------
# --stratum: is every sample in the Tractor files a member of that stratum?
# ---------------------------------------------------------------------------
# The membership comes from --strata_file (the handoff's ancestry/strata.tsv: sample_id, stratum)
# or, without it, from a 'stratum' column of the phenotype file. With neither, the check cannot
# run, and the result says so ("not_run"), so it is never mistaken for a pass. POOLED, and no
# --stratum at all, have nothing to check.
check_stratum_membership <- function(stratum, samples, strata_file = NULL, pheno = NULL) {
    if (is.null(stratum)) return(list(status = "not_applicable", source = "no --stratum"))
    if (toupper(stratum) == "POOLED") return(list(status = "not_applicable", source = "POOLED"))
    tab <- NULL; src <- NULL
    if (!is.null(strata_file)) {
        if (!file.exists(strata_file)) stop("--strata_file ", strata_file, " does not exist")
        tab <- fread(strata_file, colClasses = "character")
        if (!all(c("sample_id", "stratum") %in% names(tab))) stop(strata_file, ": needs columns sample_id and stratum")
        src <- strata_file
    } else if (!is.null(pheno) && "stratum" %in% names(pheno)) {
        tab <- data.table(sample_id = as.character(pheno$sample_id), stratum = as.character(pheno$stratum))
        src <- "the phenotype file's 'stratum' column"
    }
    if (is.null(tab)) {
        return(list(status = "not_run", source = "no --strata_file and no 'stratum' column in the phenotype file"))
    }
    lab <- tab$stratum[match(samples, tab$sample_id)]
    absent <- samples[is.na(lab)]
    other  <- samples[!is.na(lab) & lab != stratum]
    if (length(absent) > 0 || length(other) > 0) {
        stop(sprintf(paste0("--stratum %s: of the %d samples in the Tractor files, %d belong to another stratum in %s",
                            "%s and %d are not listed there%s"),
                     stratum, length(samples), length(other), src,
                     if (length(other)) paste0(" (e.g. ", other[1], ": ", lab[match(other[1], samples)], ")") else "",
                     length(absent), if (length(absent)) paste0(" (e.g. ", absent[1], ")") else ""))
    }
    list(status = "passed", source = src)
}

# ---------------------------------------------------------------------------
# --min_stratum_n: how many samples can enter the null model?
# ---------------------------------------------------------------------------
# A sample counts if it is in both the phenotype and the Tractor files and has no missing value in
# the outcome column(s) or in any covariate the phenotype file has.
n_analysable <- function(pheno, samples, outcome_cols, covariates = NULL) {
    keep <- pheno$sample_id %in% samples
    cols <- intersect(c(outcome_cols, covariates), names(pheno))
    if (length(cols) == 0) return(sum(keep))
    sum(complete.cases(pheno[keep, cols, with = FALSE]))
}

# ---------------------------------------------------------------------------
# The counts sidecar the run record reads: <output_prefix>.run_counts.tsv, key<TAB>value
# ---------------------------------------------------------------------------
# status          ok | skipped_below_min_stratum_n
# stratum         --stratum, or NA
# stratum_check   passed | not_run | not_applicable  (a failed check stops the run instead)
# min_stratum_n   the floor applied
# n_null          samples entering the null model (complete cases)
# n_variants      variants read
# n_pass_mac      variants where at least one ancestry passed the MAC gate
# n_tested        variants with a joint-test P value
# n_failed        variants whose test raised an error (P left NA)
write_run_counts <- function(output_prefix, counts) {
    f <- paste0(output_prefix, ".run_counts.tsv")
    vals <- vapply(counts, function(x) if (is.null(x) || length(x) == 0 || is.na(x[1])) "NA" else as.character(x[1]), "")
    # data.frame, not data.table(key = ...): 'key' is an argument of data.table() itself
    fwrite(data.frame(key = names(counts), value = unname(vals)), f, sep = "\t")
    invisible(f)
}
