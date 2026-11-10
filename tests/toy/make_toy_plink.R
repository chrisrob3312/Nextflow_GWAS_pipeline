#!/usr/bin/env Rscript
# ============================================================================
# make_toy_plink.R: a small PLINK 1 set + phenotype file for checks on the cluster
# ============================================================================
# Base R only, so it runs anywhere. Three simulated populations (allele frequencies drawn apart)
# so PC-AiR has structure to find; 22 autosomes. Used by the cluster checks in the GWAS
# verification notes (bin/genesis_pcair_pcrelate.R on a toy set; CP-11 sample accounting, CP-13
# duplicate removal).
#   Rscript tests/toy/make_toy_plink.R <out_prefix> [n_samples=150] [n_snps=3000] [options]
# options (after the numbers, any order):
#   duplicate        add T_DUP, an exact copy of T001's genotypes (a duplicate / MZ pair)
#   extra_pheno=K    K phenotype rows for samples that are not genotyped (default 3)
# Writes <out_prefix>.bed/.bim/.fam and <out_prefix>.phenotypes.tsv (the handoff contract's
# columns: FID IID sample_id mrd efs_time efs_event os_time os_event protocol + covariates).
# ============================================================================
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 1) stop("usage: make_toy_plink.R <out_prefix> [n_samples] [n_snps] [duplicate] [extra_pheno=K]")
out <- a[1]
num <- suppressWarnings(as.integer(a[-1])); num <- num[!is.na(num)]
n <- if (length(num) >= 1) num[1] else 150L
m <- if (length(num) >= 2) num[2] else 3000L
dup <- "duplicate" %in% a
k_extra <- { x <- grep("^extra_pheno=", a, value = TRUE); if (length(x)) as.integer(sub(".*=", "", x)) else 3L }
set.seed(20261110)

ids <- sprintf("T%03d", seq_len(n))
pop <- rep(1:3, length.out = n)
p_anc <- runif(m, 0.1, 0.9)
p_pop <- sapply(1:3, function(k) pmin(pmax(rbeta(m, p_anc * 20, (1 - p_anc) * 20), 0.02), 0.98))   # m x 3
G <- sapply(seq_len(n), function(i) rbinom(m, 2, p_pop[, pop[i]]))                                 # m x n, count of A1
if (dup) { G <- cbind(G, G[, 1]); ids <- c(ids, "T_DUP"); pop <- c(pop, pop[1]) }
n_all <- length(ids)

# .bed: magic 0x6c 0x1b, 0x01 = SNP-major; per SNP, 4 samples per byte, first sample in the low
# bits; code 00 = two A1, 10 = het, 11 = two A2, 01 = missing
code <- matrix(c(0L, 2L, 3L)[3L - G], nrow = m)          # A1 count 2 -> 00, 1 -> 10, 0 -> 11
pad <- (4 - n_all %% 4) %% 4
if (pad) code <- cbind(code, matrix(0L, m, pad))
nb <- ncol(code) / 4
bytes <- integer(m * nb)
for (j in seq_len(nb)) {
    cols <- (j - 1) * 4 + 1:4
    bytes[seq(j, by = nb, length.out = m)] <- code[, cols[1]] + 4L * code[, cols[2]] + 16L * code[, cols[3]] + 64L * code[, cols[4]]
}
con <- file(paste0(out, ".bed"), "wb")
writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con); writeBin(as.raw(bytes), con); close(con)

chr <- sort(rep_len(1:22, m)); pos <- ave(seq_len(m), chr, FUN = function(x) seq_along(x) * 10000L)
a1 <- "A"; a2 <- "G"
write.table(data.frame(chr, sprintf("chr%d:%d:%s:%s", chr, pos, a2, a1), 0, pos, a1, a2),
            paste0(out, ".bim"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(data.frame(ids, ids, 0, 0, 0, -9), paste0(out, ".fam"), sep = " ", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

ph_ids <- c(ids, if (k_extra > 0) sprintf("NOGENO%02d", seq_len(k_extra)))
np <- length(ph_ids)
ph <- data.frame(FID = ph_ids, IID = ph_ids, sample_id = ph_ids,
                 mrd = rbinom(np, 1, 0.3),
                 efs_time = round(rexp(np, 0.2), 2), efs_event = rbinom(np, 1, 0.4),
                 os_time = round(rexp(np, 0.1), 2), os_event = rbinom(np, 1, 0.2),
                 protocol = sample(c("P1", "P2"), np, TRUE),
                 sex = rbinom(np, 1, 0.5), age_at_diagnosis = round(runif(np, 1, 18), 1),
                 cd19_pos_pct = round(runif(np, 50, 99), 1), storage_days = sample(100:3000, np, TRUE))
write.table(ph, paste0(out, ".phenotypes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("toy set:", out, "-", n_all, "samples,", m, "SNPs, 22 chromosomes;", np, "phenotype rows (", k_extra,
    "not genotyped )", if (dup) "; T_DUP duplicates T001" else "", "\n")
