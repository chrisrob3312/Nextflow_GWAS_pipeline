#!/usr/bin/env Rscript
# ============================================================================
# make_toy_tractor.R: small Tractor files + phenotype for running the adapter end to end
# ============================================================================
# Base R only. Each sample carries two haplotypes whose local ancestry switches in tracts along the
# chromosome; alleles are drawn from ancestry-specific frequencies. Written in the adapter's label
# layout (bin/tractor_io.R):
#   <prefix>.hapcount.<ANC>.txt.gz   CHROM POS ID REF ALT <sample...>   haplotypes of ANC (0-2)
#   <prefix>.dosage.<ANC>.txt.gz     same columns                        ALT alleles on them (0-2)
#   <prefix>.phenotypes.tsv          FID IID sample_id mrd efs_time efs_event os_time os_event protocol
#                                    sex age_at_diagnosis cd19_pos_pct storage_days PC1 PC2
#   <prefix>.samples.txt             the sample IDs, one per line
# usage: Rscript tests/toy/make_toy_tractor.R <prefix> [n_samples=120] [n_variants=200] [ancestries=AFR,EUR]
#        [chrom=22] [seed=1] [na_rows=K] [id_prefix=T]
#   na_rows=K   the first K samples get NA for mrd and os_time (complete-case checks)
# The hapcount of each sample sums to 2 across ancestries; the dosage is <= the hapcount.
# ============================================================================
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 1) stop("usage: make_toy_tractor.R <prefix> [n_samples] [n_variants] [ancestries] [key=value ...]")
prefix <- a[1]
pos_args <- a[-1][!grepl("=", a[-1])]
kv <- a[-1][grepl("=", a[-1])]
opt <- function(k, d) { x <- grep(paste0("^", k, "="), kv, value = TRUE); if (length(x)) sub("^[^=]*=", "", x[1]) else d }
n <- if (length(pos_args) >= 1) as.integer(pos_args[1]) else 120L
m <- if (length(pos_args) >= 2) as.integer(pos_args[2]) else 200L
anc <- if (length(pos_args) >= 3) strsplit(pos_args[3], ",")[[1]] else c("AFR", "EUR")
chrom <- opt("chrom", "22"); seed <- as.integer(opt("seed", "1")); na_rows <- as.integer(opt("na_rows", "0"))
id_prefix <- opt("id_prefix", "T")
set.seed(seed)
k <- length(anc)
ids <- sprintf("%s%04d", id_prefix, seq_len(n))

# allele frequencies per ancestry, apart enough that the backgrounds differ
p_base <- runif(m, 0.15, 0.6)
freq <- sapply(seq_len(k), function(j) pmin(pmax(p_base + rnorm(m, 0, 0.08), 0.08), 0.85))   # m x k
# local ancestry of each haplotype: tracts of ~25 variants, drawn from the sample's global proportions
props <- if (k == 1) matrix(1, n, 1) else { x <- matrix(rgamma(n * k, 2), n, k); x / rowSums(x) }
hap_anc <- function(i) { la <- integer(m); pos <- 1
    while (pos <= m) { len <- 1 + rpois(1, 25); la[pos:min(m, pos + len - 1)] <- sample.int(k, 1, prob = props[i, ]); pos <- pos + len }
    la }
H <- array(0L, c(m, n, k)); D <- array(0L, c(m, n, k))
for (i in seq_len(n)) for (h in 1:2) {
    la <- hap_anc(i)
    allele <- rbinom(m, 1, freq[cbind(seq_len(m), la)])
    for (j in seq_len(k)) { on <- la == j; H[on, i, j] <- H[on, i, j] + 1L; D[on, i, j] <- D[on, i, j] + allele[on] }
}

pos <- 100000L + 1000L * seq_len(m)
info <- data.frame(CHROM = chrom, POS = pos, ID = sprintf("%s:%d:A:G", chrom, pos), REF = "A", ALT = "G")
write_mat <- function(x, f) {
    con <- gzfile(f, "w")
    writeLines(paste(c(names(info), ids), collapse = "\t"), con)
    write.table(cbind(info, x), con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    close(con)
}
for (j in seq_len(k)) {
    write_mat(H[, , j], sprintf("%s.hapcount.%s.txt.gz", prefix, anc[j]))
    write_mat(D[, , j], sprintf("%s.dosage.%s.txt.gz", prefix, anc[j]))
}

# phenotype (the handoff contract's columns) with a weak genetic signal at variant 10
g <- rowSums(sapply(seq_len(k), function(j) D[10, , j]))
age <- round(runif(n, 1, 18), 1); sex <- rbinom(n, 1, 0.5)
lp <- 0.6 * g + 0.02 * age
os_time <- round(rexp(n, 0.05 * exp(lp)), 2); os_event <- rbinom(n, 1, 0.6)
ph <- data.frame(FID = ids, IID = ids, sample_id = ids,
                 mrd = rbinom(n, 1, plogis(-0.5 + 0.5 * g)),
                 efs_time = round(rexp(n, 0.08), 2), efs_event = rbinom(n, 1, 0.5),
                 os_time = os_time, os_event = os_event, protocol = sample(c("P1", "P2"), n, TRUE),
                 sex = sex, age_at_diagnosis = age, cd19_pos_pct = round(runif(n, 50, 99), 1),
                 storage_days = sample(100:3000, n, TRUE), PC1 = round(rnorm(n), 4), PC2 = round(rnorm(n), 4))
if (na_rows > 0) { ph$mrd[seq_len(na_rows)] <- NA; ph$os_time[seq_len(na_rows)] <- NA }
write.table(ph, paste0(prefix, ".phenotypes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
writeLines(ids, paste0(prefix, ".samples.txt"))
cat("toy Tractor set:", prefix, "-", n, "samples,", m, "variants on chr", chrom, ", ancestries", paste(anc, collapse = ","),
    if (na_rows > 0) paste0("; ", na_rows, " phenotype rows with NA mrd / os_time") else "", "\n")
