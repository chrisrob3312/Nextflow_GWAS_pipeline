#!/usr/bin/env Rscript
# Unit checks for bin/tractor_io.R (needs only data.table). Prints "ok <what>" or "FAIL <what>" per
# check and exits 1 if any failed. Run by tests/test_chassis.sh:  Rscript tests/test_tractor_io.R <repo>
args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) > 0) args[1] else "."
source(file.path(repo, "bin", "tractor_io.R"))

failed <- 0L
check <- function(what, cond) {
    if (isTRUE(cond)) cat("ok  ", what, "\n") else { cat("FAIL", what, "\n"); failed <<- failed + 1L }
}
err_of <- function(expr) tryCatch({ force(expr); "" }, error = function(e) conditionMessage(e))

d <- tempfile("tractor_io_"); dir.create(d)
tab <- data.table(CHR = "chr1", POS = 1:3, ID = paste0("v", 1:3), REF = "A", ALT = "G", S1 = 0, S2 = 1, S3 = 2)
put <- function(name) { f <- file.path(d, name); fwrite(tab, f, sep = "\t"); invisible(f) }

# --- file resolution ---------------------------------------------------------
p1 <- file.path(d, "lab")
put("lab.hapcount.EUR.txt.gz"); put("lab.EUR.dosage.txt.gz")
r <- resolve_tractor_file(p1, "EUR", "hapcount")
check("label-named hapcount file is found", basename(r$path) == "lab.hapcount.EUR.txt.gz" && r$how == "named by label")
check("label-named dosage file (<prefix>.<anc>.dosage.txt.gz) is found",
      basename(resolve_tractor_file(p1, "EUR", "dosage")$path) == "lab.EUR.dosage.txt.gz")

p2 <- file.path(d, "num")
put("num.anc0.dosage.txt.gz"); put("num.anc1.dosage.txt.gz"); put("num.anc2.dosage.txt.gz")
e <- err_of(resolve_tractor_file(p2, "EUR", "dosage"))
check("numbered files without an index map: refused", grepl("numbered by position", e) && grepl("--ancestry_index_map", e))
check("the refusal names the numbered files it saw", grepl("num.anc0.dosage.txt.gz", e))

mapf <- file.path(d, "flare_index_map.tsv")
writeLines(c("anc_index\tlabel", "0\tAFR", "1\tAMR", "2\tEUR"), mapf)
m <- read_index_map(mapf)
r <- resolve_tractor_file(p2, "EUR", "dosage", m)
check("with the index map, EUR is the file numbered 2 (not position 1 of 'AFR,EUR')",
      basename(r$path) == "num.anc2.dosage.txt.gz" && grepl("index 2", r$how))
check("with the index map, AFR is the file numbered 0", basename(resolve_tractor_file(p2, "AFR", "dosage", m)$path) == "num.anc0.dosage.txt.gz")
e <- err_of(resolve_tractor_file(p2, "SAS", "dosage", m))
check("an ancestry the map does not have: refused", grepl("not in the index map", e))
writeLines(c("anc_index\tlabel", "0\tAFR", "1\tAFR"), mapf)
check("a map with a label twice: refused", grepl("appears twice", err_of(read_index_map(mapf))))
writeLines(c("index\tname", "0\tAFR"), mapf)
check("a map without anc_index/label: refused", grepl("header must be anc_index", err_of(read_index_map(mapf))))
check("label-named files win even with a map given", basename(resolve_tractor_file(p1, "EUR", "hapcount", m)$path) == "lab.hapcount.EUR.txt.gz")

check("sample ids come from the header only", identical(tractor_sample_ids(file.path(d, "lab.hapcount.EUR.txt.gz")), c("S1", "S2", "S3")))

# --- --stratum ----------------------------------------------------------------
sf <- file.path(d, "strata.tsv")
fwrite(data.table(sample_id = c("S1", "S2", "S3", "S4"), stratum = c("EUR", "EUR", "EUR", "AAC")), sf, sep = "\t")
check("all samples in the stratum: passed", check_stratum_membership("EUR", c("S1", "S2", "S3"), sf)$status == "passed")
e <- err_of(check_stratum_membership("EUR", c("S1", "S4"), sf))
check("a sample of another stratum: refused, naming it", grepl("1 belong to another stratum", e) && grepl("S4: AAC", e))
e <- err_of(check_stratum_membership("EUR", c("S1", "S9"), sf))
check("a sample missing from the strata file: refused", grepl("1 are not listed there", e) && grepl("S9", e))
check("POOLED: not applicable", check_stratum_membership("POOLED", c("S1", "S4"), sf)$status == "not_applicable")
check("no --stratum: not applicable", check_stratum_membership(NULL, c("S1"), sf)$status == "not_applicable")
ph <- data.table(sample_id = c("S1", "S2", "S3"), OS_time = c(1, 2, NA), OS_status = c(1, 0, 1), age = c(5, NA, 7))
check("no strata file and no 'stratum' column: not_run (never a pass)", check_stratum_membership("EUR", c("S1"), NULL, ph)$status == "not_run")
ph2 <- copy(ph)[, stratum := c("EUR", "EUR", "LAT1")]
check("a 'stratum' phenotype column is used when there is no strata file",
      check_stratum_membership("EUR", c("S1", "S2"), NULL, ph2)$status == "passed" &&
      grepl("another stratum", err_of(check_stratum_membership("EUR", c("S1", "S3"), NULL, ph2))))
check("a missing --strata_file is an error", grepl("does not exist", err_of(check_stratum_membership("EUR", "S1", file.path(d, "none.tsv")))))

# --- --min_stratum_n --------------------------------------------------------------
check("complete cases over outcome + covariates", n_analysable(ph, c("S1", "S2", "S3"), c("OS_time", "OS_status"), "age") == 1)
check("only samples in the Tractor files count", n_analysable(ph, c("S1"), c("OS_time", "OS_status")) == 1)
check("a covariate the phenotype lacks is not counted against anyone", n_analysable(ph, c("S1", "S2"), "OS_status", "not_a_column") == 2)

# --- counts file ---------------------------------------------------------------
f <- write_run_counts(file.path(d, "chr1"), list(status = "skipped_below_min_stratum_n", stratum = "EUR", n_null = 12, n_tested = NA))
lines <- readLines(f)
check("run_counts.tsv: key/value, the text NA for a missing value",
      lines[1] == "key\tvalue" && "status\tskipped_below_min_stratum_n" %in% lines && "n_tested\tNA" %in% lines)

unlink(d, recursive = TRUE)
quit(save = "no", status = if (failed > 0) 1L else 0L)
