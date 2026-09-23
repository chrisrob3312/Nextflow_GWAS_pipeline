#!/usr/bin/env Rscript

# ============================================================================
# GxG (Epistasis) Interaction Testing on GWAS Hits
# ============================================================================
# END-OF-PIPELINE step. Takes the top hits from the GWAS (Tractor-GENESIS
# and/or standard GWAS) plus known leukemia risk loci, and tests all pairwise
# SNP x SNP interactions:
#
#   1. POOLED cohort   : y ~ g1 + g2 + g1:g2 + covariates (+ stratum main effect)
#   2. WITHIN STRATUM  : same model fit separately in EUR, AAC, LAT1, LAT2, ...
#   3. ANCESTRY MODIFICATION of the interaction (does the g1:g2 effect differ
#      by ancestry background?):
#        a) Cochran's Q / I2 across the stratum-specific interaction betas
#        b) Pooled LRT: model with g1:g2:stratum vs model without it
#   4. OPTIONAL local-ancestry-aware epistasis (Tractor files): interaction
#      between ancestry-deconvoluted dosages, e.g. Dose1_AMR x Dose2_AMR, to
#      ask whether the interaction only exists when both risk alleles sit on
#      the same ancestral haplotype background.
#
# Models: binary (logistic), quantitative (linear), survival (Cox)
# Multiple testing: Bonferroni over pairs + BH FDR, reported per analysis.
#
# SLURM: --stratum selects one stratum (array over strata) OR
#        --array_index/--array_total splits the pair list (array over pairs).
# ============================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
    library(survival)
    library(parallel)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

# ============================================================================
# Command-line arguments
# ============================================================================
option_list <- list(
    # Hit selection
    make_option("--sumstats", type = "character", default = NULL,
                help = "Comma-separated GWAS summary statistics files to select hits from"),
    make_option("--p_threshold", type = "numeric", default = 1e-5,
                help = "P-value threshold for selecting hits [default: 1e-5 (suggestive)]"),
    make_option("--p_col", type = "character", default = NULL,
                help = "P-value column (auto: P_JOINT, P, PVAL, LOG10P)"),
    make_option("--max_hits", type = "integer", default = 200,
                help = "Maximum number of hit SNPs to carry into pairwise testing [default: 200]"),
    make_option("--known_loci", type = "character", default = NULL,
                help = "TSV of known risk loci to always include (columns: snp, gene, ...)"),
    make_option("--snp_list", type = "character", default = NULL,
                help = "Explicit SNP list (one ID per line); skips hit selection if given"),
    make_option("--custom_variants", type = "character", default = NULL,
                help = "Custom GRCh38 variant list, one per line: rsID, chr:pos, or chr:pos:ref:alt. Always included, highest priority."),
    make_option("--no_default_loci", action = "store_true", default = FALSE,
                help = "Do not add the default known-risk-loci table given by --known_loci"),
    make_option("--prune_mode", type = "character", default = "conditional",
                help = "conditional = test each LD-correlated hit CONDITIONAL on the hits already kept; keep it if it stays significant [default]; variant = keep only the strongest per LD cluster; pair = keep all hits, skip LD pairs"),
    make_option("--cond_p_threshold", type = "numeric", default = 1e-4,
                help = "Conditional P below which an LD-correlated hit is kept as an independent signal [default: 1e-4]"),
    make_option("--min_distance_kb", type = "numeric", default = 1000,
                help = "Skip pairs on the same chromosome closer than this [default: 1000 kb]"),
    make_option("--max_pair_r2", type = "numeric", default = 0.2,
                help = "Skip pairs with genotype r2 above this (LD, not epistasis) [default: 0.2]"),

    # Genotypes
    make_option("--geno", type = "character", default = NULL,
                help = "PLINK prefix (bed/bim/fam) for the full cohort"),
    make_option("--gds", type = "character", default = NULL,
                help = "GDS file (alternative to --geno)"),
    make_option("--tractor_prefix", type = "character", default = NULL,
                help = "Tractor prefix (.hapcount.{anc}.txt.gz / .dosage.{anc}.txt.gz) for LA-aware epistasis"),
    make_option("--ancestries", type = "character", default = "EUR,AFR,AMR",
                help = "Ancestral populations in Tractor files [default: EUR,AFR,AMR]"),

    # Phenotype / model
    make_option("--phenotype", type = "character", default = NULL,
                help = "Phenotype file (TSV, sample_id column)"),
    make_option("--trait", type = "character", default = NULL,
                help = "Trait column"),
    make_option("--model", type = "character", default = "binary",
                help = "binary, quantitative, survival [default: binary]"),
    make_option("--time_col", type = "character", default = "time"),
    make_option("--event_col", type = "character", default = "event"),
    make_option("--covariates", type = "character", default = NULL,
                help = "Comma-separated covariates (age,sex,PC1..PC5)"),
    make_option("--unrelated_ids", type = "character", default = NULL,
                help = "Optional file of unrelated sample IDs (GLM/Cox assume independence)"),

    # Ancestry stratification
    make_option("--ancestry_col", type = "character", default = "GRAF_ANC",
                help = "Ancestry column in phenotype file [default: GRAF_ANC]"),
    make_option("--stratum", type = "character", default = NULL,
                help = "Run only this stratum (POOLED, EUR, AAC, LAT1, LAT2, EAS, SAS, OTHER)"),
    make_option("--min_stratum_n", type = "integer", default = 30,
                help = "Minimum N per stratum; smaller groups pool to OTHER [default: 30]"),
    make_option("--ancestry_config", type = "character", default = NULL,
                help = "Path to ancestry_config.R (cohort stratum rules)"),

    # Output / runtime
    make_option(c("-o", "--output_prefix"), type = "character", default = "gxg"),
    make_option("--threads", type = "integer", default = 4),
    make_option("--array_index", type = "integer", default = NULL,
                help = "SLURM array index to split the pair list"),
    make_option("--array_total", type = "integer", default = NULL,
                help = "Total SLURM array tasks for pair splitting"),
    make_option(c("-v", "--verbose"), action = "store_true", default = FALSE)
)

opt <- parse_args(OptionParser(option_list = option_list, prog = "run_gxg_interaction.R",
                               description = "Pairwise GxG interaction tests on GWAS hits, overall and by ancestry"))

if (is.null(opt$phenotype) || is.null(opt$trait)) stop("Required: --phenotype, --trait")
if (is.null(opt$geno) && is.null(opt$gds)) stop("Required: --geno or --gds")
if (is.null(opt$sumstats) && is.null(opt$snp_list) && is.null(opt$custom_variants) &&
    (is.null(opt$known_loci) || opt$no_default_loci)) {
    stop("Required: --sumstats and/or --custom_variants and/or --snp_list and/or --known_loci")
}
if (!opt$prune_mode %in% c("conditional", "variant", "pair")) stop("--prune_mode must be 'conditional', 'variant' or 'pair'")

# SLURM detection
slurm_cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "")))
if (!is.na(slurm_cpus) && slurm_cpus > 0) opt$threads <- slurm_cpus
if (is.null(opt$array_index)) {
    aid <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "")))
    if (!is.na(aid)) opt$array_index <- aid
}
if (is.null(opt$array_total)) {
    acount <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_COUNT", "")))
    if (!is.na(acount)) opt$array_total <- acount
}

MIN_STRATUM_N <- opt$min_stratum_n
PRIMARY_ANCESTRIES <- c("EUR", "AAC", "LAT1", "LAT2")
CONDITIONAL_ANCESTRIES <- c("EAS", "SAS")
ALWAYS_POOL <- c("AMR", "AFR", "UNKNOWN")
GRAF_ANC_CODES <- c("107" = "AAC", "100" = "AFR", "300" = "EUR", "500" = "EAS",
                    "600" = "AMR", "601" = "LAT1", "602" = "LAT2")
if (!is.null(opt$ancestry_config) && file.exists(opt$ancestry_config)) {
    source(opt$ancestry_config)
}

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║        GxG Interaction Testing on GWAS Hits (Epistasis)          ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║ Trait: %-57s ║\n", opt$trait))
cat(sprintf("║ Model: %-57s ║\n", opt$model))
cat(sprintf("║ Hit threshold: P < %-45g ║\n", opt$p_threshold))
cat(sprintf("║ Stratum: %-55s ║\n", opt$stratum %||% "ALL (pooled + each ancestry)"))
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# ============================================================================
# 1. Select hit SNPs
# ============================================================================
find_col <- function(dt, candidates) {
    hit <- intersect(candidates, names(dt))
    if (length(hit) == 0) NA_character_ else hit[1]
}

# Variant lookup (GRCh38) from the genotype file, used to resolve user-supplied
# IDs given as rsID, chr:pos or chr:pos:ref:alt to the IDs actually in the data
bim <- NULL
if (!is.null(opt$geno) && file.exists(paste0(opt$geno, ".bim"))) {
    bim <- fread(paste0(opt$geno, ".bim"), col.names = c("CHR", "ID", "CM", "POS", "A1", "A2"))
    bim[, CHR := sub("^chr", "", as.character(CHR))]
}
resolve_ids <- function(ids, chr = NULL, pos = NULL) {
    ids <- as.character(ids)
    if (is.null(bim)) return(ids)                       # no lookup possible: use as given
    out <- rep(NA_character_, length(ids))
    direct <- ids %in% bim$ID
    out[direct] <- ids[direct]
    for (k in which(!direct)) {
        parts <- strsplit(ids[k], "[:_]")[[1]]
        c_ <- if (length(parts) >= 2) sub("^chr", "", parts[1]) else if (!is.null(chr)) sub("^chr", "", as.character(chr[k])) else NA
        p_ <- if (length(parts) >= 2) suppressWarnings(as.integer(parts[2])) else if (!is.null(pos)) suppressWarnings(as.integer(pos[k])) else NA
        if (is.na(c_) || is.na(p_)) next
        cand <- bim[CHR == c_ & POS == p_]
        if (nrow(cand) == 0) next
        if (length(parts) >= 4) {                       # allele-aware match, either orientation
            a1 <- toupper(parts[3]); a2 <- toupper(parts[4])
            m <- cand[(toupper(A1) == a1 & toupper(A2) == a2) | (toupper(A1) == a2 & toupper(A2) == a1)]
            if (nrow(m) > 0) cand <- m
        }
        out[k] <- cand$ID[1]
    }
    out
}

select_hits <- function() {
    hits <- data.table()

    if (!is.null(opt$snp_list)) {
        ids <- readLines(opt$snp_list)
        ids <- ids[nzchar(ids)]
        hits <- data.table(ID = ids, CHR = NA_character_, POS = NA_integer_,
                           P = NA_real_, SOURCE = "snp_list")
    }

    # Custom GRCh38 variant list: always included, highest priority (P = -1)
    if (!is.null(opt$custom_variants) && file.exists(opt$custom_variants)) {
        raw <- readLines(opt$custom_variants)
        raw <- trimws(raw); raw <- raw[nzchar(raw) & !grepl("^#", raw)]
        raw <- sapply(strsplit(raw, "[ \t]+"), `[`, 1)       # first column only
        resolved <- resolve_ids(raw)
        n_unres <- sum(is.na(resolved))
        if (n_unres > 0) {
            cat("  WARNING:", n_unres, "custom variants not found in genotypes:",
                paste(head(raw[is.na(resolved)], 10), collapse = ", "), "\n")
        }
        ok <- !is.na(resolved)
        hits <- rbind(hits, data.table(ID = resolved[ok], CHR = NA_character_, POS = NA_integer_,
                                       P = -1, SOURCE = paste0("custom:", raw[ok])), fill = TRUE)
        cat("  Custom variants added:", sum(ok), "of", length(raw), "\n")
    }

    if (!is.null(opt$sumstats)) {
        for (f in strsplit(opt$sumstats, ",")[[1]]) {
            if (!file.exists(f)) { cat("  WARNING: missing sumstats", f, "\n"); next }
            ss <- fread(f)
            id_col  <- find_col(ss, c("ID", "SNP", "variant_id", "rsid", "MarkerName"))
            chr_col <- find_col(ss, c("CHR", "CHROM", "chr", "chromosome"))
            pos_col <- find_col(ss, c("POS", "GENPOS", "BP", "pos", "position"))
            p_col   <- opt$p_col %||% find_col(ss, c("P_JOINT", "P", "PVAL", "pvalue", "p.value", "P_META"))
            log10p  <- find_col(ss, c("LOG10P"))

            if (is.na(id_col)) { cat("  WARNING: no ID column in", f, "\n"); next }
            pvals <- if (!is.na(p_col)) ss[[p_col]] else if (!is.na(log10p)) 10^(-ss[[log10p]]) else NA_real_
            sel <- ss[!is.na(pvals) & pvals < opt$p_threshold]
            sel_p <- pvals[!is.na(pvals) & pvals < opt$p_threshold]
            if (nrow(sel) == 0) { cat("  ", basename(f), ": 0 hits\n"); next }

            hits <- rbind(hits, data.table(
                ID = as.character(sel[[id_col]]),
                CHR = if (!is.na(chr_col)) as.character(sel[[chr_col]]) else NA_character_,
                POS = if (!is.na(pos_col)) as.integer(sel[[pos_col]]) else NA_integer_,
                P = sel_p,
                SOURCE = basename(f)
            ), fill = TRUE)
            cat("  ", basename(f), ":", nrow(sel), "hits at P <", opt$p_threshold, "\n")
        }
    }

    # Default known-risk-loci table (switch off with --no_default_loci); P = 0 so it
    # ranks below custom variants but above GWAS hits when pruning
    if (!opt$no_default_loci && !is.null(opt$known_loci) && file.exists(opt$known_loci)) {
        kl <- fread(opt$known_loci)
        id_col <- find_col(kl, c("snp", "SNP", "ID", "rsid", "variant_id"))
        if (!is.na(id_col)) {
            kl_chr <- if ("chr" %in% names(kl)) as.character(kl$chr) else NULL
            kl_pos <- if ("pos_grch38" %in% names(kl)) kl$pos_grch38 else NULL
            resolved <- resolve_ids(kl[[id_col]], chr = kl_chr, pos = kl_pos)   # rsID, else GRCh38 chr:pos
            ok <- !is.na(resolved)
            if (any(!ok)) cat("  Known loci not in genotypes (skipped):",
                              paste(kl[[id_col]][!ok], collapse = ", "), "\n")
            known <- data.table(ID = resolved[ok],
                                CHR = if (!is.null(kl_chr)) kl_chr[ok] else NA_character_,
                                POS = if (!is.null(kl_pos)) as.integer(kl_pos[ok]) else NA_integer_,
                                P = 0, SOURCE = paste0("known:", if ("gene" %in% names(kl)) kl$gene[ok] else "locus"))
            hits <- rbind(hits, known, fill = TRUE)
            cat("  Known risk loci added:", nrow(known), "\n")
        }
    } else if (opt$no_default_loci) {
        cat("  Default known-loci table disabled (--no_default_loci)\n")
    }

    if (nrow(hits) == 0) stop("No hit SNPs selected")

    # Keep best P per ID (known loci get P = 0 so are always kept)
    hits <- hits[order(P)][, .SD[1], by = ID]

    # Cap
    if (nrow(hits) > opt$max_hits) {
        cat("  Capping hits from", nrow(hits), "to", opt$max_hits, "(best P first)\n")
        hits <- hits[order(P)][1:opt$max_hits]
    }
    hits
}

cat("Selecting hit SNPs...\n")
hits <- select_hits()
cat("  Total unique hits:", nrow(hits), "\n\n")

# ============================================================================
# 2. Extract genotypes for hits
# ============================================================================
extract_genotypes_plink <- function(prefix, ids, out_prefix) {
    snp_file <- paste0(out_prefix, ".hit_snps.txt")
    writeLines(ids, snp_file)
    plink2 <- Sys.which("plink2"); plink1 <- Sys.which("plink")
    if (plink2 != "") {
        cmd <- paste(plink2, "--bfile", prefix, "--extract", snp_file,
                     "--export A", "--out", paste0(out_prefix, ".hits"),
                     "--threads", opt$threads)
    } else if (plink1 != "") {
        cmd <- paste(plink1, "--bfile", prefix, "--extract", snp_file,
                     "--recode A", "--out", paste0(out_prefix, ".hits"),
                     "--threads", opt$threads)
    } else stop("plink2/plink not found for genotype extraction")
    if (opt$verbose) cat("  ", cmd, "\n")
    system(cmd, ignore.stdout = !opt$verbose, ignore.stderr = !opt$verbose)
    raw <- fread(paste0(out_prefix, ".hits.raw"))
    geno_cols <- setdiff(names(raw), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
    G <- as.matrix(raw[, ..geno_cols])
    # PLINK appends _<counted allele>; strip to ID
    colnames(G) <- sub("_[ACGTacgt0-9]+$", "", geno_cols)
    rownames(G) <- as.character(raw$IID)
    G
}

extract_genotypes_gds <- function(gds_file, ids) {
    suppressPackageStartupMessages({ library(SeqArray) })
    gds <- seqOpen(gds_file)
    on.exit(seqClose(gds))
    all_ids <- seqGetData(gds, "annotation/id")
    sel <- which(all_ids %in% ids)
    seqSetFilter(gds, variant.sel = sel, verbose = FALSE)
    G <- seqGetData(gds, "$dosage_alt")
    colnames(G) <- all_ids[sel]
    rownames(G) <- seqGetData(gds, "sample.id")
    G
}

cat("Extracting genotypes for", nrow(hits), "hits...\n")
G <- if (!is.null(opt$gds)) extract_genotypes_gds(opt$gds, hits$ID) else
     extract_genotypes_plink(opt$geno, hits$ID, opt$output_prefix)
missing_ids <- setdiff(hits$ID, colnames(G))
if (length(missing_ids) > 0) {
    cat("  WARNING:", length(missing_ids), "hits not found in genotypes (dropped)\n")
    hits <- hits[ID %in% colnames(G)]
}
G <- G[, hits$ID, drop = FALSE]
cat("  Genotype matrix:", nrow(G), "samples x", ncol(G), "SNPs\n\n")

# Fill CHR/POS from the genotype lookup if missing
if (!is.null(bim) && any(is.na(hits$CHR) | is.na(hits$POS))) {
    hits[bim, on = "ID", `:=`(CHR = ifelse(is.na(CHR), as.character(i.CHR), CHR),
                              POS = ifelse(is.na(POS), i.POS, POS))]
}
hits[, CHR := sub("^chr", "", as.character(CHR))]

# ============================================================================
# 3. Phenotype, covariates, strata
# ============================================================================
pheno <- fread(opt$phenotype)
id_col <- find_col(pheno, c("sample_id", "IID", "ID", "sample"))
if (is.na(id_col)) stop("No sample ID column in phenotype file")
pheno[[id_col]] <- as.character(pheno[[id_col]])
if (!is.null(opt$unrelated_ids)) {
    keep <- readLines(opt$unrelated_ids)
    pheno <- pheno[get(id_col) %in% keep]
    cat("Restricted to", nrow(pheno), "unrelated samples\n")
}

covs <- if (!is.null(opt$covariates)) strsplit(opt$covariates, ",")[[1]] else character(0)
missing_covs <- setdiff(covs, names(pheno))
if (length(missing_covs) > 0) stop("Covariates not in phenotype: ", paste(missing_covs, collapse = ","))

# Ancestry → stratum
if (opt$ancestry_col %in% names(pheno)) {
    anc_raw <- as.character(pheno[[opt$ancestry_col]])
    anc_name <- ifelse(anc_raw %in% names(GRAF_ANC_CODES), GRAF_ANC_CODES[anc_raw], anc_raw)
    anc_name[is.na(anc_name) | anc_name == ""] <- "UNKNOWN"
    counts <- table(anc_name)
    pheno$STRATUM <- sapply(anc_name, function(a) {
        n <- counts[a]
        if (a %in% PRIMARY_ANCESTRIES)      return(if (n >= MIN_STRATUM_N) a else "OTHER")
        if (a %in% CONDITIONAL_ANCESTRIES)  return(if (n >= MIN_STRATUM_N) a else "OTHER")
        if (a %in% ALWAYS_POOL)             return("OTHER")
        if (n < MIN_STRATUM_N)              return("OTHER")
        a
    })
    cat("Analysis strata:\n")
    print(table(pheno$STRATUM))
    cat("\n")
} else {
    cat("No ancestry column '", opt$ancestry_col, "' - pooled analysis only\n\n", sep = "")
    pheno$STRATUM <- "POOLED_ONLY"
}

# Align samples
common <- intersect(rownames(G), pheno[[id_col]])
pheno <- pheno[match(common, get(id_col))]
G <- G[common, , drop = FALSE]
cat("Samples with genotype + phenotype:", length(common), "\n\n")

# Outcome
if (opt$model == "survival") {
    if (!all(c(opt$time_col, opt$event_col) %in% names(pheno))) stop("Survival columns missing")
    keep <- !is.na(pheno[[opt$time_col]]) & !is.na(pheno[[opt$event_col]])
} else {
    if (!opt$trait %in% names(pheno)) stop("Trait '", opt$trait, "' not in phenotype")
    keep <- !is.na(pheno[[opt$trait]])
}
pheno <- pheno[keep]; G <- G[keep, , drop = FALSE]

# ============================================================================
# 3b. Independent-signal selection among hits: CONDITIONAL testing, not LD-drop
# ============================================================================
# Precedence: GWAS hits by P first, then the custom list, then the default
# known loci. Walking down that order, a candidate in LD (r2 > max_pair_r2)
# with variants already kept is NOT dropped outright. It is tested
# CONDITIONAL on those kept LD partners in the pooled cohort:
#     y ~ candidate + kept LD partners + covariates (+ stratum)
# and kept as an independent signal if its conditional P < cond_p_threshold
# (likelihood-ratio test). Otherwise it is recorded with its proxy. This is
# an individual-level stepwise conditional / joint analysis, the thing
# GCTA-COJO approximates from summary statistics and reference LD.
# Proximity alone (< min_distance_kb, low r2) never removes a variant.
# prune_mode = variant restores plain "keep the strongest per LD cluster".
fit_conditional_p <- function(g_cand, G_part, ph, covs, model) {
    d <- data.frame(cand = g_cand, G_part, ph[, c(covs, "STRATUM",
                    intersect(c(opt$trait, opt$time_col, opt$event_col), names(ph))), with = FALSE])
    part_names <- colnames(G_part)
    cov_str <- if (length(covs) > 0) paste("+", paste(covs, collapse = " + ")) else ""
    strat_str <- if (length(unique(d$STRATUM)) > 1) "+ STRATUM" else ""
    rhs0 <- paste(c(part_names, "1"), collapse = " + ")
    f0 <- as.formula(paste("Y ~", rhs0, cov_str, strat_str))
    f1 <- as.formula(paste("Y ~ cand +", rhs0, cov_str, strat_str))
    tryCatch({
        if (model == "survival") {
            d$Y <- Surv(d[[opt$time_col]], d[[opt$event_col]])
            anova(coxph(f0, data = d), coxph(f1, data = d))[2, "Pr(>|Chi|)"]
        } else if (model == "binary") {
            d$Y <- d[[opt$trait]]
            anova(glm(f0, data = d, family = binomial()), glm(f1, data = d, family = binomial()), test = "LRT")[2, "Pr(>Chi)"]
        } else {
            d$Y <- d[[opt$trait]]
            anova(lm(f0, data = d), lm(f1, data = d))[2, "Pr(>F)"]
        }
    }, error = function(e) NA_real_)
}

if (opt$prune_mode %in% c("conditional", "variant") && ncol(G) > 1) {
    r2_all <- suppressWarnings(cor(G, use = "pairwise.complete.obs")^2)
    # precedence: 0 = GWAS hit, 1 = custom list, 2 = default known loci
    pri <- ifelse(grepl("^custom", hits$SOURCE), 1L, ifelse(grepl("^known", hits$SOURCE), 2L, 0L))
    walk <- order(pri, hits$P)
    kept <- character(0)
    status <- rep(NA_character_, nrow(hits)); pruned_by <- rep(NA_character_, nrow(hits))
    r2_proxy <- rep(NA_real_, nrow(hits)); cond_p <- rep(NA_real_, nrow(hits)); cond_on <- rep(NA_character_, nrow(hits))
    n_cond <- 0L
    for (k in walk) {
        id <- hits$ID[k]
        same_chr <- kept[hits$CHR[match(kept, hits$ID)] == hits$CHR[k]]
        same_chr <- same_chr[!is.na(same_chr)]
        partners <- character(0)
        if (length(same_chr) > 0) {
            r2v <- r2_all[id, same_chr]
            partners <- same_chr[!is.na(r2v) & r2v > opt$max_pair_r2]
        }
        if (length(partners) == 0) { kept <- c(kept, id); status[k] <- "kept_independent"; next }
        j <- which.max(r2_all[id, partners]); pruned_by[k] <- partners[j]; r2_proxy[k] <- r2_all[id, partners[j]]
        if (opt$prune_mode == "variant") { status[k] <- "dropped_ld"; next }
        # conditional test on all kept LD partners
        n_cond <- n_cond + 1L
        p_c <- fit_conditional_p(G[, id], G[, partners, drop = FALSE], pheno, covs, opt$model)
        cond_p[k] <- p_c; cond_on[k] <- paste(partners, collapse = ";")
        if (!is.na(p_c) && p_c < opt$cond_p_threshold) {
            kept <- c(kept, id); status[k] <- "kept_conditional"
        } else {
            status[k] <- "dropped_conditional"
        }
    }
    hits[, `:=`(priority = pri, status = status, kept = ID %in% kept, proxy = pruned_by,
                r2_with_proxy = r2_proxy, conditional_on = cond_on, conditional_p = cond_p)]
    fwrite(hits[order(priority, P)], paste0(opt$output_prefix, ".gxg.hits_pruned.tsv"), sep = "\t")
    cat("Independent-signal selection (r2 >", opt$max_pair_r2, "; GWAS > custom > known):\n")
    cat("  kept without LD partner:", sum(status == "kept_independent", na.rm = TRUE), "\n")
    if (opt$prune_mode == "conditional") {
        cat("  conditional tests run:", n_cond, "| kept as independent (cond P <", opt$cond_p_threshold, "):",
            sum(status == "kept_conditional", na.rm = TRUE), "| dropped (explained by partner):",
            sum(status == "dropped_conditional", na.rm = TRUE), "\n")
    } else {
        cat("  dropped (LD with stronger hit):", sum(status == "dropped_ld", na.rm = TRUE), "\n")
    }
    hits <- hits[kept == TRUE]
    G <- G[, hits$ID, drop = FALSE]
} else if (opt$prune_mode == "pair") {
    cat("prune_mode = pair: all hits retained; LD pairs skipped at the pair level\n")
}

# ============================================================================
# 4. Build pair list (drop LD / proximal pairs)
# ============================================================================
n_snp <- ncol(G)
pairs <- as.data.table(t(combn(n_snp, 2)))
setnames(pairs, c("i", "j"))
pairs[, `:=`(SNP1 = colnames(G)[i], SNP2 = colnames(G)[j])]
pairs <- merge(pairs, hits[, .(SNP1 = ID, CHR1 = CHR, POS1 = POS, SRC1 = SOURCE)], by = "SNP1")
pairs <- merge(pairs, hits[, .(SNP2 = ID, CHR2 = CHR, POS2 = POS, SRC2 = SOURCE)], by = "SNP2")

# Proximity filter
pairs[, proximal := !is.na(CHR1) & !is.na(CHR2) & CHR1 == CHR2 &
                     abs(POS1 - POS2) < opt$min_distance_kb * 1000]
# LD filter (genotype r2 in the tested cohort)
r2 <- suppressWarnings(cor(G, use = "pairwise.complete.obs")^2)
pairs[, r2 := r2[cbind(i, j)]]
pairs[, in_ld := !is.na(r2) & r2 > opt$max_pair_r2]
n_all <- nrow(pairs)
pairs <- pairs[!proximal & !in_ld]
cat("Pairs:", n_all, "total,", n_all - nrow(pairs), "dropped (proximal/LD),",
    nrow(pairs), "to test\n")

# SLURM pair splitting
if (!is.null(opt$array_index) && !is.null(opt$array_total) && is.null(opt$stratum)) {
    chunk <- ceiling(nrow(pairs) / opt$array_total)
    idx0 <- (opt$array_index - 1) * chunk + 1
    idx1 <- min(opt$array_index * chunk, nrow(pairs))
    pairs <- pairs[idx0:idx1]
    cat("SLURM array", opt$array_index, "/", opt$array_total, ": testing pairs", idx0, "-", idx1, "\n")
}
cat("\n")

# ============================================================================
# 5. Interaction test functions
# ============================================================================
fit_interaction <- function(g1, g2, ph, covs, model, extra_terms = NULL) {
    d <- data.frame(g1 = g1, g2 = g2, ph[, c(covs, "STRATUM",
                                             intersect(c(opt$trait, opt$time_col, opt$event_col), names(ph))),
                                          with = FALSE])
    if (length(unique(d$g1[!is.na(d$g1)])) < 2 || length(unique(d$g2[!is.na(d$g2)])) < 2) return(NULL)
    cov_str <- if (length(covs) > 0) paste("+", paste(covs, collapse = " + ")) else ""
    extra <- if (!is.null(extra_terms)) paste("+", extra_terms) else ""

    f_null <- as.formula(paste("Y ~ g1 + g2", cov_str, extra))
    f_int  <- as.formula(paste("Y ~ g1 + g2 + g1:g2", cov_str, extra))

    tryCatch({
        if (model == "survival") {
            d$Y <- Surv(d[[opt$time_col]], d[[opt$event_col]])
            m0 <- coxph(f_null, data = d); m1 <- coxph(f_int, data = d)
            co <- summary(m1)$coefficients["g1:g2", ]
            list(beta = co["coef"], se = co["se(coef)"], p_wald = co["Pr(>|z|)"],
                 p_lrt = anova(m0, m1)[2, "Pr(>|Chi|)"], n = m1$n, n_events = m1$nevent)
        } else if (model == "binary") {
            d$Y <- d[[opt$trait]]
            m0 <- glm(f_null, data = d, family = binomial()); m1 <- glm(f_int, data = d, family = binomial())
            co <- summary(m1)$coefficients["g1:g2", ]
            list(beta = co["Estimate"], se = co["Std. Error"], p_wald = co["Pr(>|z|)"],
                 p_lrt = anova(m0, m1, test = "LRT")[2, "Pr(>Chi)"], n = nobs(m1),
                 n_events = sum(d$Y == 1, na.rm = TRUE))
        } else {
            d$Y <- d[[opt$trait]]
            m0 <- lm(f_null, data = d); m1 <- lm(f_int, data = d)
            co <- summary(m1)$coefficients["g1:g2", ]
            list(beta = co["Estimate"], se = co["Std. Error"], p_wald = co["Pr(>|t|)"],
                 p_lrt = anova(m0, m1)[2, "Pr(>F)"], n = nobs(m1), n_events = NA)
        }
    }, error = function(e) NULL)
}

# Ancestry modification of interaction (pooled 3-way LRT)
fit_interaction_by_ancestry <- function(g1, g2, ph, covs, model) {
    if (length(unique(ph$STRATUM)) < 2) return(NULL)
    d <- data.frame(g1 = g1, g2 = g2, ph[, c(covs, "STRATUM",
                                             intersect(c(opt$trait, opt$time_col, opt$event_col), names(ph))),
                                          with = FALSE])
    d$STRATUM <- factor(d$STRATUM)
    cov_str <- if (length(covs) > 0) paste("+", paste(covs, collapse = " + ")) else ""
    f0 <- as.formula(paste("Y ~ (g1 + g2 + g1:g2) + STRATUM + g1:STRATUM + g2:STRATUM", cov_str))
    f1 <- as.formula(paste("Y ~ (g1 + g2 + g1:g2) * STRATUM", cov_str))
    tryCatch({
        if (model == "survival") {
            d$Y <- Surv(d[[opt$time_col]], d[[opt$event_col]])
            m0 <- coxph(f0, data = d); m1 <- coxph(f1, data = d)
            anova(m0, m1)[2, "Pr(>|Chi|)"]
        } else if (model == "binary") {
            d$Y <- d[[opt$trait]]
            m0 <- glm(f0, data = d, family = binomial()); m1 <- glm(f1, data = d, family = binomial())
            anova(m0, m1, test = "LRT")[2, "Pr(>Chi)"]
        } else {
            d$Y <- d[[opt$trait]]
            m0 <- lm(f0, data = d); m1 <- lm(f1, data = d)
            anova(m0, m1)[2, "Pr(>F)"]
        }
    }, error = function(e) NA_real_)
}

cochran_q <- function(betas, ses) {
    ok <- !is.na(betas) & !is.na(ses) & ses > 0
    if (sum(ok) < 2) return(list(Q = NA, df = NA, p = NA, I2 = NA))
    w <- 1 / ses[ok]^2; bm <- sum(w * betas[ok]) / sum(w)
    Q <- sum(w * (betas[ok] - bm)^2); df <- sum(ok) - 1
    list(Q = Q, df = df, p = pchisq(Q, df, lower.tail = FALSE),
         I2 = max(0, (Q - df) / Q) * 100)
}

# ============================================================================
# 6. Run tests: pooled, per stratum, ancestry modification
# ============================================================================
strata_all <- sort(unique(pheno$STRATUM))
run_pooled <- is.null(opt$stratum) || opt$stratum == "POOLED"
strata_to_run <- if (is.null(opt$stratum)) strata_all else
                 if (opt$stratum == "POOLED") character(0) else opt$stratum
strata_to_run <- intersect(strata_to_run, strata_all)
if (length(strata_to_run) == 1 && strata_to_run == "POOLED_ONLY") strata_to_run <- character(0)

test_pair_set <- function(ph, G_sub, label, extra_terms = NULL) {
    cat("  [", label, "] N =", nrow(ph), "| pairs =", nrow(pairs), "\n")
    res <- mclapply(seq_len(nrow(pairs)), function(k) {
        r <- fit_interaction(G_sub[, pairs$SNP1[k]], G_sub[, pairs$SNP2[k]], ph, covs, opt$model, extra_terms)
        if (is.null(r)) return(NULL)
        data.table(analysis = label, SNP1 = pairs$SNP1[k], SNP2 = pairs$SNP2[k],
                   CHR1 = pairs$CHR1[k], POS1 = pairs$POS1[k], CHR2 = pairs$CHR2[k], POS2 = pairs$POS2[k],
                   SOURCE1 = pairs$SRC1[k], SOURCE2 = pairs$SRC2[k], r2_pair = pairs$r2[k],
                   N = r$n, N_EVENTS = r$n_events, BETA_INT = r$beta, SE_INT = r$se,
                   OR_INT = if (opt$model != "quantitative") exp(r$beta) else NA_real_,
                   P_WALD = r$p_wald, P_LRT = r$p_lrt)
    }, mc.cores = opt$threads)
    out <- rbindlist(res[!sapply(res, is.null)])
    if (nrow(out) > 0) {
        out[, P_BONF := pmin(1, P_LRT * nrow(pairs))]
        out[, P_FDR := p.adjust(P_LRT, method = "BH")]
        out <- out[order(P_LRT)]
    }
    out
}

all_results <- list()

if (run_pooled) {
    cat("POOLED cohort (stratum as covariate main effect)...\n")
    extra <- if (length(strata_all) > 1 && !("POOLED_ONLY" %in% strata_all)) "STRATUM" else NULL
    all_results[["POOLED"]] <- test_pair_set(pheno, G, "POOLED", extra)
}

for (s in strata_to_run) {
    idx <- pheno$STRATUM == s
    if (sum(idx) < MIN_STRATUM_N) { cat("  Skipping", s, "(N =", sum(idx), ")\n"); next }
    cat("Stratum", s, "...\n")
    all_results[[s]] <- test_pair_set(pheno[idx], G[idx, , drop = FALSE], s)
}

# Ancestry modification of the interaction
het_results <- data.table()
if (is.null(opt$stratum) && length(strata_to_run) >= 2) {
    cat("\nAncestry modification of interaction effects...\n")
    strat_res <- rbindlist(all_results[strata_to_run], fill = TRUE)
    if (nrow(strat_res) > 0) {
        het_results <- strat_res[, {
            q <- cochran_q(BETA_INT, SE_INT)
            best <- which.min(P_LRT)
            list(n_strata = .N,
                 strata = paste(analysis, collapse = ";"),
                 betas = paste(analysis, "=", signif(BETA_INT, 3), collapse = ";"),
                 strongest_stratum = analysis[best],
                 strongest_p = P_LRT[best],
                 Q_HET_INT = q$Q, DF_HET = q$df, P_HET_INT = q$p, I2_INT = q$I2)
        }, by = .(SNP1, SNP2)]

        # Pooled 3-way LRT (g1:g2:STRATUM)
        cat("  3-way LRT (g1:g2:ancestry) on", nrow(het_results), "pairs...\n")
        p3 <- mclapply(seq_len(nrow(het_results)), function(k) {
            fit_interaction_by_ancestry(G[, het_results$SNP1[k]], G[, het_results$SNP2[k]], pheno, covs, opt$model)
        }, mc.cores = opt$threads)
        het_results[, P_3WAY_LRT := unlist(p3)]
        het_results[, P_HET_FDR := p.adjust(P_HET_INT, method = "BH")]
        het_results[, P_3WAY_FDR := p.adjust(P_3WAY_LRT, method = "BH")]
        het_results <- het_results[order(P_3WAY_LRT)]
    }
}

# ============================================================================
# 7. Optional: local-ancestry-aware epistasis (Tractor dosages)
# ============================================================================
la_results <- data.table()
if (!is.null(opt$tractor_prefix) && run_pooled) {
    cat("\nLocal-ancestry-aware epistasis (Tractor dosage x dosage)...\n")
    tr_anc <- strsplit(opt$ancestries, ",")[[1]]
    read_tractor_rows <- function(f, ids) {
        if (!file.exists(f)) return(NULL)
        dt <- fread(f)
        idc <- find_col(dt, c("ID", "SNP", "variant_id"))
        if (is.na(idc)) { dt[, ID := paste0(dt[[1]], ":", dt[[2]])]; idc <- "ID" }
        dt <- dt[get(idc) %in% ids]
        meta_cols <- intersect(names(dt), c("CHROM", "CHR", "POS", "ID", "REF", "ALT", "SNP", "variant_id"))
        m <- as.matrix(dt[, setdiff(names(dt), meta_cols), with = FALSE])
        rownames(m) <- dt[[idc]]
        m
    }
    dose <- list()
    for (a in tr_anc) {
        f <- paste0(opt$tractor_prefix, ".dosage.", a, ".txt.gz")
        dose[[a]] <- read_tractor_rows(f, colnames(G))
        if (is.null(dose[[a]])) cat("  WARNING: missing", f, "\n")
    }
    dose <- dose[!sapply(dose, is.null)]
    if (length(dose) > 0) {
        samp <- intersect(colnames(dose[[1]]), pheno[[id_col]])
        ph_la <- pheno[match(samp, get(id_col))]
        for (a in names(dose)) {
            D <- t(dose[[a]][, samp, drop = FALSE])
            ok_pairs <- pairs[SNP1 %in% colnames(D) & SNP2 %in% colnames(D)]
            if (nrow(ok_pairs) == 0) next
            cat("  Background", a, ": Dose_", a, "(SNP1) x Dose_", a, "(SNP2) on", nrow(ok_pairs), "pairs\n", sep = "")
            res <- mclapply(seq_len(nrow(ok_pairs)), function(k) {
                r <- fit_interaction(D[, ok_pairs$SNP1[k]], D[, ok_pairs$SNP2[k]], ph_la, covs, opt$model,
                                     extra_terms = if (length(strata_all) > 1) "STRATUM" else NULL)
                if (is.null(r)) return(NULL)
                data.table(analysis = paste0("LA_", a), SNP1 = ok_pairs$SNP1[k], SNP2 = ok_pairs$SNP2[k],
                           N = r$n, BETA_INT = r$beta, SE_INT = r$se, P_WALD = r$p_wald, P_LRT = r$p_lrt)
            }, mc.cores = opt$threads)
            out <- rbindlist(res[!sapply(res, is.null)])
            if (nrow(out) > 0) {
                out[, P_FDR := p.adjust(P_LRT, method = "BH")]
                la_results <- rbind(la_results, out, fill = TRUE)
            }
        }
        if (nrow(la_results) > 0) la_results <- la_results[order(P_LRT)]
    }
}

# ============================================================================
# 8. Write outputs (sorted by significance)
# ============================================================================
suffix <- if (!is.null(opt$array_index) && is.null(opt$stratum)) paste0(".part", opt$array_index) else ""
out_files <- c()
for (nm in names(all_results)) {
    f <- paste0(opt$output_prefix, ".gxg.", nm, suffix, ".tsv")
    fwrite(all_results[[nm]], f, sep = "\t"); out_files <- c(out_files, f)
}
combined <- rbindlist(all_results, fill = TRUE)
if (nrow(combined) > 0) {
    f <- paste0(opt$output_prefix, ".gxg.all", suffix, ".tsv")
    fwrite(combined[order(P_LRT)], f, sep = "\t"); out_files <- c(out_files, f)
    sig <- combined[P_FDR < 0.05 | P_BONF < 0.05]
    if (nrow(sig) > 0) {
        f <- paste0(opt$output_prefix, ".gxg.significant", suffix, ".tsv")
        fwrite(sig[order(P_LRT)], f, sep = "\t"); out_files <- c(out_files, f)
    }
}
if (nrow(het_results) > 0) {
    f <- paste0(opt$output_prefix, ".gxg.ancestry_heterogeneity", suffix, ".tsv")
    fwrite(het_results, f, sep = "\t"); out_files <- c(out_files, f)
}
if (nrow(la_results) > 0) {
    f <- paste0(opt$output_prefix, ".gxg.local_ancestry", suffix, ".tsv")
    fwrite(la_results, f, sep = "\t"); out_files <- c(out_files, f)
}
fwrite(hits, paste0(opt$output_prefix, ".gxg.hits_tested.tsv"), sep = "\t")

# ============================================================================
# 9. Summary
# ============================================================================
cat("\n")
cat("╔══════════════════════════════════════════════════════════════════╗\n")
cat("║                     GxG ANALYSIS COMPLETE                        ║\n")
cat("╠══════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║ Hits tested: %-51d ║\n", nrow(hits)))
cat(sprintf("║ Pairs tested: %-50d ║\n", nrow(pairs)))
for (nm in names(all_results)) {
    r <- all_results[[nm]]
    n_sig <- if (nrow(r) > 0) sum(r$P_FDR < 0.05, na.rm = TRUE) else 0
    n_bonf <- if (nrow(r) > 0) sum(r$P_BONF < 0.05, na.rm = TRUE) else 0
    cat(sprintf("║ %-10s FDR<0.05: %-5d Bonferroni<0.05: %-5d              ║\n", nm, n_sig, n_bonf))
}
if (nrow(het_results) > 0) {
    cat(sprintf("║ Ancestry-modified interactions (3-way FDR<0.05): %-15d ║\n",
                sum(het_results$P_3WAY_FDR < 0.05, na.rm = TRUE)))
    top <- het_results[1]
    cat(sprintf("║ Top: %s x %s (strongest in %s)%*s║\n", top$SNP1, top$SNP2, top$strongest_stratum,
                max(1, 60 - nchar(top$SNP1) - nchar(top$SNP2) - nchar(top$strongest_stratum)), ""))
}
cat("╚══════════════════════════════════════════════════════════════════╝\n")
cat("Outputs:\n"); for (f in out_files) cat("  ", f, "\n")
