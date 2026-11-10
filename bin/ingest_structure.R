#!/usr/bin/env Rscript
# ============================================================================
# ingest_structure.R: take the population structure the handoff bundle ships (CP-07)
# ============================================================================
# The harmonization pipeline computes PC-AiR and PC-Relate at both scales and ships them
# (docs/HANDOFF_CONTRACT.md, rule 4). This script checks them and lays them out for the GWAS
# step; it computes nothing. bin/genesis_pcair_pcrelate.R is kept as the comparator (CP-18).
#
# Reads, from --handoff_dir:
#   genotypes/cohort.hardcall.maf01.fam        the sample key (IID, column 2)
#   ancestry/strata.tsv                         sample_id, stratum
#   structure/cohort.pcair.pcs.tsv              sample_id, PC1..PCn   (every .fam IID)
#   structure/strata/<S>.pcair.pcs.tsv          sample_id, PC1..PCn   (exactly S's members)
#   structure/cohort.pcrelate.grm.rds / .kinship.rds   dimnames = .fam IID
#   the phenotype file (--phenotype)
# Stops, naming the file and the problem, if anything is missing or does not match. Never falls
# back to computing.
#
# Writes, to --out_dir (default results/pcair) and --strata_dir (default results/strata):
#   cohort.phenotypes.with_pcs.tsv     phenotype + cohort PCs (POOLED runs)
#   <S>.phenotypes.with_pcs.tsv        phenotype of S's members + S's within-stratum PCs
#   cohort.pcrelate.grm.rds, cohort.pcrelate.kinship.rds   symbolic links to the shipped files
#   structure_source.tsv               where each file came from
#   <strata_dir>/<S>.samples.txt       S's members, one IID per line
#   qc_assertions.tsv                  this step's QC verdicts (also filed to GWAS_QC_DIR)
# Existing PC columns (PC1, PC2, ...) of the phenotype file are replaced, never mixed.
# QC (CP-11, bin/qc_assert.R): files pcair_pc_merge once per output phenotype file (every row has
# its PCs) and structure_ingest (PASS at the end, FAIL with the reason on any refusal), and one
# sample_accounting row per join: phenotype -> .fam IIDs -> stratum -> PC merge.
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(optparse) })
`%||%` <- function(a, b) if (is.null(a)) b else a
local({
    f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
    source(file.path(if (length(f) > 0) dirname(normalizePath(f[1])) else "bin", "qc_assert.R"))
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--handoff_dir", type = "character", help = "the handoff bundle directory"),
    make_option("--phenotype", type = "character", default = NULL,
                help = "phenotype file [default: <handoff_dir>/phenotypes/phenotypes.tsv]"),
    make_option("--out_dir", type = "character", default = "results/pcair"),
    make_option("--strata_dir", type = "character", default = "results/strata")
)))
if (is.null(opt$handoff_dir)) stop("--handoff_dir is required")
H <- normalizePath(opt$handoff_dir, mustWork = TRUE)
pheno_file <- opt$phenotype %||% file.path(H, "phenotypes", "phenotypes.tsv")

refuse <- function(...) {
    qc_assert("structure_ingest", FALSE, message = paste(...))
    cat("REFUSED:", ..., "\n", file = stderr()); quit(save = "no", status = 2)
}
need <- function(rel) {
    f <- if (startsWith(rel, "/")) rel else file.path(H, rel)
    if (!file.exists(f) || file.size(f) == 0) refuse("structure ingest: missing or empty", f)
    f
}

fam_file    <- need("genotypes/cohort.hardcall.maf01.fam")
strata_file <- need("ancestry/strata.tsv")
cohort_pcs  <- need("structure/cohort.pcair.pcs.tsv")
grm_file    <- need("structure/cohort.pcrelate.grm.rds")
kin_file    <- need("structure/cohort.pcrelate.kinship.rds")
pheno_file  <- need(pheno_file)

fam_ids <- fread(fam_file, header = FALSE, colClasses = "character")[[2]]
if (anyDuplicated(fam_ids)) refuse(fam_file, "has a duplicated IID")

# --- GRM and kinship: dimnames are the .fam IIDs -----------------------------------
for (f in c(grm_file, kin_file)) {
    m <- readRDS(f)
    rn <- rownames(m); cn <- colnames(m)
    if (is.null(rn) || is.null(cn)) refuse(f, "has no dimnames (the contract keys it by .fam IID)")
    if (!identical(rn, cn)) refuse(f, "row and column names differ")
    if (length(rn) != length(fam_ids) || !setequal(rn, fam_ids)) {
        refuse(f, sprintf("dimnames are not the .fam IIDs (%d names, %d IIDs; e.g. %s)", length(rn), length(fam_ids),
                          paste(head(c(setdiff(rn, fam_ids), setdiff(fam_ids, rn)), 3), collapse = ", ")))
    }
    cat("ok  ", basename(f), ":", length(rn), "x", length(rn), "dimnames == .fam IID\n")
}

# --- strata --------------------------------------------------------------------------
strata <- fread(strata_file, colClasses = "character")
if (!all(c("sample_id", "stratum") %in% names(strata))) refuse(strata_file, "needs columns sample_id and stratum")
if (anyDuplicated(strata$sample_id)) refuse(strata_file, "lists a sample twice")
if (!setequal(strata$sample_id, fam_ids)) refuse(strata_file, "sample set differs from the .fam IIDs")
stratum_names <- sort(unique(strata$stratum))

# --- PC files ------------------------------------------------------------------------
read_pcs <- function(f, expected_ids, what) {
    p <- fread(f, colClasses = list(character = "sample_id"))
    if (names(p)[1] != "sample_id" || !all(grepl("^PC[0-9]+$", names(p)[-1])) || ncol(p) < 2) {
        refuse(f, "header must be sample_id PC1 ... PCn")
    }
    tryCatch(check_covariates(p, names(p)[-1], f), error = function(e) refuse(conditionMessage(e)))   # CP-08: numeric PCs
    if (anyDuplicated(p$sample_id) || !setequal(p$sample_id, expected_ids)) {
        refuse(f, sprintf("samples differ from %s (%d rows, %d expected)", what, nrow(p), length(expected_ids)))
    }
    p
}
pcs <- list(cohort = read_pcs(cohort_pcs, fam_ids, "the .fam IIDs"))
for (S in stratum_names) {
    f <- need(file.path("structure", "strata", paste0(S, ".pcair.pcs.tsv")))
    pcs[[S]] <- read_pcs(f, strata$sample_id[strata$stratum == S], paste0("stratum ", S, "'s members in strata.tsv"))
}

# --- outputs -------------------------------------------------------------------------
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$strata_dir, recursive = TRUE, showWarnings = FALSE)
qc_copy_to(file.path(opt$out_dir, "qc_assertions.tsv"))
pheno <- fread(pheno_file)
if (!"sample_id" %in% names(pheno)) refuse(pheno_file, "has no sample_id column")
pheno[, sample_id := as.character(sample_id)]
base <- pheno[, setdiff(names(pheno), grep("^PC[0-9]+$", names(pheno), value = TRUE)), with = FALSE]
geno_rows <- base[sample_id %in% fam_ids]
qc_account("phenotype_to_fam", nrow(base), nrow(geno_rows), stratum = "cohort",
           reason = sprintf("phenotype rows with a .fam IID (%d .fam IIDs have no phenotype row)",
                            sum(!fam_ids %in% base$sample_id)))
write_with_pcs <- function(ph, p, out, S) {
    ph <- merge_pcs_checked(ph, p, stratum = S, what = basename(out))
    fwrite(ph, out, sep = "\t")
    cat("ok  ", basename(out), ":", nrow(ph), "samples,", ncol(p) - 1, "PCs\n")
}
write_with_pcs(geno_rows, pcs$cohort, file.path(opt$out_dir, "cohort.phenotypes.with_pcs.tsv"), "cohort")
for (S in stratum_names) {
    ids <- strata$sample_id[strata$stratum == S]
    in_s <- geno_rows[sample_id %in% ids]
    qc_account("stratum", nrow(geno_rows), nrow(in_s), stratum = S,
               reason = sprintf("phenotype rows in stratum %s (%d members in strata.tsv)", S, length(ids)))
    write_with_pcs(in_s, pcs[[S]], file.path(opt$out_dir, paste0(S, ".phenotypes.with_pcs.tsv")), S)
    writeLines(ids, file.path(opt$strata_dir, paste0(S, ".samples.txt")))
}
link <- function(from, to) { unlink(to); if (!file.symlink(from, to)) refuse("could not link", to, "to", from) }
link(grm_file, file.path(opt$out_dir, "cohort.pcrelate.grm.rds"))
link(kin_file, file.path(opt$out_dir, "cohort.pcrelate.kinship.rds"))
fwrite(data.table(file = c("cohort.pcair.pcs.tsv", paste0("strata/", stratum_names, ".pcair.pcs.tsv"),
                           "cohort.pcrelate.grm.rds", "cohort.pcrelate.kinship.rds"),
                  source = c(cohort_pcs, file.path(H, "structure", "strata", paste0(stratum_names, ".pcair.pcs.tsv")),
                             grm_file, kin_file)),
       file.path(opt$out_dir, "structure_source.tsv"), sep = "\t")
qc_assert("structure_ingest", TRUE, n = length(fam_ids),
          message = sprintf("bundle structure checked and laid out: %d samples, %d strata", length(fam_ids), length(stratum_names)))
cat("structure ingested from bundle:", H, "(", length(stratum_names), "strata:", paste(stratum_names, collapse = ", "), ")\n")
