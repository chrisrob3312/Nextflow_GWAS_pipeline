#!/usr/bin/env Rscript

# ============================================================================
# GENESIS PC-AiR + PC-Relate: ancestry PCs and GRM for an admixed cohort
# ============================================================================
# Standard GENESIS workflow (Conomos et al. 2015, 2016):
#   1. LD-prune SNPs (SNPRelate)
#   2. KING-robust kinship (robust to admixture at close degrees)
#   3. PC-AiR : PCs from an unrelated set, projected onto relatives
#              -> ancestry PCs NOT distorted by family structure
#   4. PC-Relate : kinship conditional on the PC-AiR PCs
#              -> relatedness NOT inflated by ancestry / admixture
#   5. Iterate once: PC-AiR with PC-Relate kinship, then PC-Relate again
#
# Outputs (used by every downstream model: standard GENESIS GWAS,
# Tractor-GENESIS, GxG, PRS validation):
#   <prefix>.pcair.pcs.tsv           sample_id, PC1..PCn (PC-AiR)
#   <prefix>.pcrelate.kinship.rds    kinship matrix (scaleKin = 1)
#   <prefix>.pcrelate.grm.rds        GRM = 2 x kinship (scaleKin = 2)  <- pass as --kinship
#   <prefix>.unrelated.txt           PC-AiR unrelated set (for GLM/Cox-only steps)
#   <prefix>.phenotypes.with_pcs.tsv phenotype file with PC1..PCn replaced by PC-AiR PCs
#   <prefix>.pcair.variance.tsv      variance explained per PC
#   <prefix>.removed_duplicates.tsv  kept_iid removed_iid kinship reason: one of each pair with
#                                    kinship > 0.3536 (duplicate / MZ), left out of the phenotype
#   <dir of prefix>/qc_assertions.tsv  this run's QC verdicts (pcair_pc_merge), also filed to GWAS_QC_DIR
# ============================================================================

suppressPackageStartupMessages({
    library(optparse); library(data.table)
    library(SNPRelate); library(GWASTools); library(GENESIS)
})

opt <- parse_args(OptionParser(option_list = list(
    make_option("--geno", type = "character", default = NULL, help = "PLINK prefix (bed/bim/fam)"),
    make_option("--gds", type = "character", default = NULL, help = "SNP GDS file (alternative to --geno)"),
    make_option("--phenotype", type = "character", default = NULL, help = "Phenotype TSV (sample_id column)"),
    make_option("--n_pcs", type = "integer", default = 10, help = "PCs to output / condition on [default: 10]"),
    make_option("--n_pcs_pcrelate", type = "integer", default = 5, help = "PCs used inside PC-Relate [default: 5]"),
    make_option("--kin_thresh", type = "numeric", default = 2^(-9/2), help = "Kinship threshold for unrelated set (3rd degree) [default: 0.0442]"),
    make_option("--div_thresh", type = "numeric", default = -2^(-9/2), help = "Divergence threshold [default: -0.0442]"),
    make_option("--ld_r2", type = "numeric", default = 0.1, help = "LD pruning r2 (sqrt applied to SNPRelate threshold) [default: 0.1]"),
    make_option("--ld_window_bp", type = "integer", default = 500000, help = "LD pruning window [default: 500 kb]"),
    make_option("--maf", type = "numeric", default = 0.01, help = "MAF filter for PCA / kinship SNPs [default: 0.01]"),
    make_option("--missing_rate", type = "numeric", default = 0.05, help = "Max missing rate per SNP [default: 0.05]"),
    make_option("--iterations", type = "integer", default = 2, help = "PC-AiR / PC-Relate iterations [default: 2]"),
    make_option("--sample_include", type = "character", default = NULL, help = "Optional file of sample IDs to keep"),
    make_option(c("-o", "--output_prefix"), type = "character", default = "cohort"),
    make_option("--threads", type = "integer", default = 4),
    make_option("--seed", type = "integer", default = 1234)
)))

if (is.null(opt$geno) && is.null(opt$gds)) stop("Required: --geno or --gds")
`%||%` <- function(a, b) if (is.null(a)) b else a
local({   # QC verdicts and sample accounting (CP-11), from next to this script
    f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
    d <- if (length(f) > 0) dirname(normalizePath(f[1])) else "bin"
    source(file.path(d, "qc_assert.R")); source(file.path(d, "duplicates.R"))   # CP-13: duplicate / MZ removal
})
qc_copy_to(file.path(dirname(opt$output_prefix), "qc_assertions.tsv"))
slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
if (!is.na(slurm_cpus) && slurm_cpus > 0) opt$threads <- slurm_cpus
set.seed(opt$seed)

cat("\n╔══════════════════════════════════════════════════════════════════╗\n")
cat("║   GENESIS PC-AiR + PC-Relate: ancestry PCs and GRM               ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# ---------------------------------------------------------------------------
# 1. GDS
# ---------------------------------------------------------------------------
gds_file <- opt$gds
if (is.null(gds_file)) {
    gds_file <- paste0(opt$output_prefix, ".snp.gds")
    cat("Converting PLINK -> SNP GDS:", gds_file, "\n")
    snpgdsBED2GDS(paste0(opt$geno, ".bed"), paste0(opt$geno, ".fam"), paste0(opt$geno, ".bim"),
                  gds_file, cvt.chr = "int", verbose = FALSE)
}
gds <- snpgdsOpen(gds_file)
sample_ids <- gds_ids <- as.character(read.gdsn(index.gdsn(gds, "sample.id")))
if (!is.null(opt$sample_include)) {
    keep <- readLines(opt$sample_include); sample_ids <- intersect(sample_ids, keep)
    qc_account("sample_include", length(gds_ids), length(sample_ids), stratum = "cohort",
               reason = paste("genotyped samples listed in", basename(opt$sample_include)))
}
cat("Samples:", length(sample_ids), "\n")

# ---------------------------------------------------------------------------
# 2. LD pruning
# ---------------------------------------------------------------------------
cat("\nLD pruning (r2 <", opt$ld_r2, ", window", opt$ld_window_bp / 1000, "kb, MAF >", opt$maf, ")...\n")
pruned <- snpgdsLDpruning(gds, sample.id = sample_ids, method = "corr", slide.max.bp = opt$ld_window_bp,
                          ld.threshold = sqrt(opt$ld_r2), maf = opt$maf, missing.rate = opt$missing_rate,
                          autosome.only = TRUE, num.thread = opt$threads, verbose = FALSE)
pruned <- unlist(pruned, use.names = FALSE)
cat("  Pruned SNPs:", length(pruned), "\n")

# ---------------------------------------------------------------------------
# 3. KING-robust kinship (starting point)
# ---------------------------------------------------------------------------
cat("\nKING-robust kinship...\n")
king <- snpgdsIBDKING(gds, sample.id = sample_ids, snp.id = pruned, num.thread = opt$threads, verbose = FALSE)
KINGmat <- king$kinship; dimnames(KINGmat) <- list(king$sample.id, king$sample.id)
cat("  Pairs with kinship > 3rd degree:", sum(KINGmat[upper.tri(KINGmat)] > opt$kin_thresh), "\n")
# Call rate on the pruned SNPs: the tie-break when one of a duplicate / MZ pair is removed (CP-13).
# SNPRelate returns the samples in the GDS order, which sample_ids keeps (intersect keeps order).
call_rate <- setNames(1 - snpgdsSampMissRate(gds, sample.id = sample_ids, snp.id = pruned, with.id = FALSE), sample_ids)
snpgdsClose(gds)

# ---------------------------------------------------------------------------
# 4-5. PC-AiR -> PC-Relate, iterated
# ---------------------------------------------------------------------------
geno <- GdsGenotypeReader(gds_file)
genoData <- GenotypeData(geno)
kinobj <- KINGmat; divobj <- KINGmat

for (it in seq_len(opt$iterations)) {
    cat("\nIteration", it, "of", opt$iterations, "\n")
    cat("  PC-AiR...\n")
    pca <- pcair(genoData, kinobj = kinobj, divobj = divobj, snp.include = pruned,
                 sample.include = sample_ids, kin.thresh = opt$kin_thresh, div.thresh = opt$div_thresh,
                 num.cores = opt$threads, verbose = FALSE)
    cat("    Unrelated set:", length(pca$unrels), "| related:", length(pca$rels), "\n")
    varprop <- pca$varprop[1:opt$n_pcs]
    cat("    Variance explained PC1-", opt$n_pcs, ": ", paste(sprintf("%.1f%%", 100 * varprop), collapse = " "), "\n", sep = "")

    cat("  PC-Relate (conditioning on", opt$n_pcs_pcrelate, "PCs, training on unrelated set)...\n")
    iterator <- GenotypeBlockIterator(genoData, snpBlock = 20000, snpInclude = pruned)
    pcrel <- pcrelate(iterator, pcs = pca$vectors[, 1:opt$n_pcs_pcrelate, drop = FALSE],
                      training.set = pca$unrels, sample.include = sample_ids, verbose = FALSE)
    kinobj <- pcrelateToMatrix(pcrel, scaleKin = 1, verbose = FALSE)
    divobj <- kinobj
}
close(genoData)

kin_mat <- as.matrix(kinobj)
grm     <- kin_mat * 2

# ---------------------------------------------------------------------------
# 5b. Duplicates / MZ twins (CP-13): one member of each pair with kinship > 0.3536 leaves the
#     phenotype step 1 reads (bin/duplicates.R; the bundle path applies the same rule)
# ---------------------------------------------------------------------------
dups <- select_duplicates(kin_mat, call_rate)
dup_file <- paste0(opt$output_prefix, ".removed_duplicates.tsv")
fwrite(dups, dup_file, sep = "\t")             # header-only when nothing is removed
cat("\nDuplicates / MZ twins (kinship > ", sprintf("%.4f", DUP_KINSHIP), "): ", nrow(dups), " removed -> ", dup_file, "\n", sep = "")
qc_assert("duplicates_removed", TRUE, n = nrow(dups),
          message = sprintf("PC-Relate kinship checked for pairs > %.4f: %d sample(s) removed from the GWAS phenotype (%s)",
                            DUP_KINSHIP, nrow(dups), basename(dup_file)))

# ---------------------------------------------------------------------------
# 6. Write outputs
# ---------------------------------------------------------------------------
cat("\nWriting outputs...\n")
pcs <- as.data.table(pca$vectors[, 1:opt$n_pcs, drop = FALSE])
setnames(pcs, paste0("PC", 1:opt$n_pcs))
pcs[, sample_id := rownames(pca$vectors)]
setcolorder(pcs, "sample_id")
fwrite(pcs, paste0(opt$output_prefix, ".pcair.pcs.tsv"), sep = "\t")

fwrite(data.table(PC = paste0("PC", 1:opt$n_pcs), variance_explained = varprop),
       paste0(opt$output_prefix, ".pcair.variance.tsv"), sep = "\t")
saveRDS(kin_mat, paste0(opt$output_prefix, ".pcrelate.kinship.rds"))
saveRDS(grm,     paste0(opt$output_prefix, ".pcrelate.grm.rds"))
writeLines(pca$unrels, paste0(opt$output_prefix, ".unrelated.txt"))

# Relatedness summary
kin_up <- kin_mat[upper.tri(kin_mat)]
# Degree boundaries (Manichaikul et al. 2010): kinship 2^(-3/2) = 0.3536 between
# duplicate/MZ and 1st degree, 2^(-5/2) = 0.1768 between 1st and 2nd, 2^(-7/2) = 0.0884
# between 2nd and 3rd, 2^(-9/2) = 0.0442 below 3rd. The counts were shifted one degree
# (1st counted above 0.3536); the labels are now printed from the constants the counts use,
# as in the harmonization pipeline's Module 9.
b_dup <- 2^(-3/2); b_1st <- 2^(-5/2); b_2nd <- 2^(-7/2); b_3rd <- 2^(-9/2)
bnd <- function(x) sprintf("%.4f", x)
rel_summary <- data.table(
    degree  = c(paste0("duplicate_or_MZ (> ", bnd(b_dup), ")"),
                paste0("1st (", bnd(b_1st), "-", bnd(b_dup), "]"),
                paste0("2nd (", bnd(b_2nd), "-", bnd(b_1st), "]"),
                paste0("3rd (", bnd(b_3rd), "-", bnd(b_2nd), "]")),
    n_pairs = c(sum(kin_up > b_dup),
                sum(kin_up > b_1st & kin_up <= b_dup),
                sum(kin_up > b_2nd & kin_up <= b_1st),
                sum(kin_up > b_3rd & kin_up <= b_2nd))
)
fwrite(rel_summary, paste0(opt$output_prefix, ".pcrelate.relatedness_summary.tsv"), sep = "\t")
print(rel_summary)

# Phenotype with PC-AiR PCs (replaces any existing PC columns)
if (!is.null(opt$phenotype) && file.exists(opt$phenotype)) {
    ph <- fread(opt$phenotype)
    idc <- intersect(c("sample_id", "IID", "ID", "sample"), names(ph))[1]
    if (is.na(idc)) stop("No sample ID column in phenotype")
    ph[[idc]] <- as.character(ph[[idc]])
    old_pcs <- grep("^PC[0-9]+$", names(ph), value = TRUE)
    if (length(old_pcs) > 0) { cat("  Replacing existing", length(old_pcs), "PC columns with PC-AiR PCs\n"); ph[, (old_pcs) := NULL] }
    out_ph <- paste0(opt$output_prefix, ".phenotypes.with_pcs.tsv")
    # CP-11: every phenotype row must get its PCs. merge_pcs_checked keeps FID / IID first, writes the
    # pc_merge sample-accounting row and files pcair_pc_merge (FAIL, with examples, if a row has
    # none); the GWAS step's gate stops on a FAIL. The rows are still written, so an override
    # (QC_GATE_OVERRIDE=pcair_pc_merge) has something to proceed with.
    # Phenotype rows of samples that were never genotyped cannot enter any model: counted and
    # dropped here. A genotyped sample left without PCs (e.g. outside --sample_include) is a FAIL.
    n_ph <- nrow(ph)
    ph <- ph[ph[[idc]] %in% gds_ids]
    qc_account("phenotype_to_genotype", n_ph, nrow(ph), stratum = "cohort",
               reason = sprintf("rows of %s with a sample in the genotype file (%d genotyped samples have no phenotype row)",
                                basename(opt$phenotype), sum(!gds_ids %in% ph[[idc]])))
    n_ph <- nrow(ph)
    ph <- ph[!ph[[idc]] %in% dups$removed_iid]
    qc_account("duplicates", n_ph, nrow(ph), stratum = "cohort", reason = paste("one of each duplicate / MZ pair removed; see", basename(dup_file)))
    ph <- merge_pcs_checked(ph, pcs, id_col = idc, stratum = "cohort", what = basename(out_ph))
    # CP-08: the PC columns every model will use as covariates are there and numeric
    check_covariates(ph, paste0("PC", 1:opt$n_pcs), out_ph)
    fwrite(ph, out_ph, sep = "\t")
    cat("  Phenotype with PC-AiR PCs:", out_ph, "(", sum(!is.na(ph$PC1)), "of", nrow(ph), "samples with PCs )\n")
} else {
    qc_assert("pcair_pc_merge", FALSE, message = paste("phenotype file not given or not found (--phenotype",
                                                        opt$phenotype %||% "unset", "): no phenotype with PCs was written"))
}

cat("\nDone. Use", paste0(opt$output_prefix, ".pcrelate.grm.rds"), "as --kinship and PC1..PC",
    opt$n_pcs, "from the augmented phenotype as covariates.\n", sep = "")
