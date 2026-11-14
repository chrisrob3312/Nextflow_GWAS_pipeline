#!/usr/bin/env Rscript
# ============================================================================
# sumstats_utils.R - ONE place that knows every GWAS output layout
# ============================================================================
# Source from any R script in bin/ (bin/ is on PATH in Nextflow tasks):
#     source(Sys.which("sumstats_utils.R"))
# or, when PATH is not set up, source(file.path(dirname(this_script), "sumstats_utils.R"))
#
# read_sumstats_std(path) returns a data.table with STANDARD columns
#     CHR   integer (1-22, X=23, Y=24, MT=25)      CHR_RAW  original label
#     POS   integer
#     ID    character (falls back to CHR:POS)
#     P     numeric  (from P_JOINT | P | LOG10P | Score.pval | p.value | ...)
#     LOG10P numeric
#     BETA, SE, EA, OA, AF, N        when present
#     BETA_<anc>, SE_<anc>, P_<anc>, MAC_<anc>, P_HET, I2   when present (Tractor-GENESIS)
# attr(x, "source_format") is one of
#     tractor_genesis | regenie | saige | genesis | plink2 | spacox | bolt | mr_mega | pooled_ivw | generic
#
# Layouts handled:
#   Tractor-GENESIS : CHR POS ID REF ALT N STAT_JOINT DF_JOINT P_JOINT BETA_<anc> SE_<anc> Z_<anc> P_<anc> MAC_<anc> Q_HET P_HET I2
#   REGENIE         : CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P
#   SAIGE           : CHR POS MarkerID Allele1 Allele2 AC_Allele2 AF_Allele2 N BETA SE Tstat var p.value
#   GENESIS         : variant.id chr pos allele.index n.obs freq MAC Score Score.SE Score.Stat Score.pval Est Est.SE
#   PLINK2 glm      : #CHROM POS ID REF ALT A1 TEST OBS_CT BETA/OR SE T_STAT/Z_STAT P
#   SPA-Cox         : CHR POS MarkerID ... p.value
#   BOLT-LMM        : SNP CHR BP ALLELE1 ALLELE0 A1FREQ BETA SE P_BOLT_LMM
#   MR-MEGA / IVW   : MarkerName/ID chromosome/CHR position/POS ... P-value_association / P_JOINT
# ============================================================================

suppressPackageStartupMessages(library(data.table))

chr_to_num <- function(x) {
    x <- toupper(sub("^CHR", "", as.character(x)))
    x[x == "X"] <- "23"; x[x == "Y"] <- "24"; x[x %in% c("MT", "M")] <- "25"; x[x == "XY"] <- "26"
    suppressWarnings(as.integer(x))
}

.pick_col <- function(nm, candidates) {
    # exact (case-insensitive) match, in candidate priority order
    low <- tolower(nm)
    for (c in candidates) { i <- match(tolower(c), low); if (!is.na(i)) return(nm[i]) }
    NA_character_
}

detect_sumstats_format <- function(nm) {
    if ("P_JOINT" %in% nm && any(grepl("^BETA_[A-Za-z0-9]+$", nm))) return("tractor_genesis")
    if ("N_STRATA" %in% nm && "P_JOINT" %in% nm) return("pooled_ivw")
    if ("LOG10P" %in% nm && ("GENPOS" %in% nm || "CHROM" %in% nm)) return("regenie")
    if ("Score.pval" %in% nm) return("genesis")
    if ("p.value" %in% nm && "MarkerID" %in% nm && any(grepl("eventTime|Surv", nm, ignore.case = TRUE))) return("spacox")
    if ("p.value" %in% nm && "MarkerID" %in% nm) return("saige")
    if ("#CHROM" %in% nm && "P" %in% nm) return("plink2")
    if ("P_BOLT_LMM" %in% nm) return("bolt")
    if (any(grepl("P-value_association", nm))) return("mr_mega")
    "generic"
}

read_sumstats_std <- function(path, drop_invalid = TRUE) {
    dt <- fread(path)
    nm <- names(dt)
    fmt <- detect_sumstats_format(nm)

    chr <- .pick_col(nm, c("CHR", "CHROM", "#CHROM", "chr", "chromosome", "Chromosome"))
    pos <- .pick_col(nm, c("POS", "GENPOS", "BP", "pos", "position", "Position", "base_pair_location"))
    id  <- .pick_col(nm, c("ID", "SNP", "MarkerID", "variant_id", "variant.id", "rsid", "RSID", "MarkerName", "SNPID"))
    p   <- .pick_col(nm, c("P_JOINT", "P", "PVAL", "P_VALUE", "PVALUE", "pvalue", "p_value", "p.value", "Score.pval",
                           "P_META", "P_BOLT_LMM", "P-value_association", "Pvalue"))
    lp  <- .pick_col(nm, c("LOG10P", "LOG10_P", "neglog10_pval", "-log10P"))
    beta <- .pick_col(nm, c("BETA", "beta", "Est", "Effect", "b", "BETA_META", "beta_0"))
    se   <- .pick_col(nm, c("SE", "se", "Est.SE", "StdErr", "SE_META", "se_0"))
    ea   <- .pick_col(nm, c("ALT", "ALLELE1", "A1", "Allele2", "EA", "effect_allele", "allele.index"))
    oa   <- .pick_col(nm, c("REF", "ALLELE0", "A2", "Allele1", "OA", "other_allele"))
    af   <- .pick_col(nm, c("A1FREQ", "AF_Allele2", "freq", "EAF", "AF", "MAF", "A1_FREQ"))
    n    <- .pick_col(nm, c("N", "n", "n.obs", "OBS_CT", "N_total", "Nsample"))

    if (is.na(chr) || is.na(pos)) stop("Cannot find chromosome/position columns in ", path, " (have: ", paste(head(nm, 15), collapse = ", "), ")")
    if (is.na(p) && is.na(lp)) stop("Cannot find a P-value column in ", path, " (have: ", paste(head(nm, 15), collapse = ", "), ")")

    P <- if (!is.na(p)) suppressWarnings(as.numeric(dt[[p]])) else 10^(-suppressWarnings(as.numeric(dt[[lp]])))
    out <- data.table(CHR_RAW = as.character(dt[[chr]]), POS = suppressWarnings(as.integer(dt[[pos]])), P = P)
    out[, CHR := chr_to_num(CHR_RAW)]
    out[, ID := if (!is.na(id)) as.character(dt[[id]]) else paste0(CHR_RAW, ":", POS)]
    out[, LOG10P := -log10(P)]
    if (!is.na(beta)) out[, BETA := suppressWarnings(as.numeric(dt[[beta]]))]
    if (!is.na(se))   out[, SE := suppressWarnings(as.numeric(dt[[se]]))]
    if (!is.na(ea))   out[, EA := as.character(dt[[ea]])]
    if (!is.na(oa))   out[, OA := as.character(dt[[oa]])]
    if (!is.na(af))   out[, AF := suppressWarnings(as.numeric(dt[[af]]))]
    if (!is.na(n))    out[, N := suppressWarnings(as.numeric(dt[[n]]))]

    # Tractor-GENESIS / pooled per-ancestry columns, carried through unchanged
    anc_cols <- grep("^(BETA|SE|Z|P|MAC|I2_STRATA|N_STRATA)_[A-Za-z0-9]+$", nm, value = TRUE)
    anc_cols <- setdiff(anc_cols, c("P_JOINT", "P_META", "BETA_META", "SE_META", "P_HET", "P_BOLT_LMM", "P_VALUE"))
    for (cc in anc_cols) out[[cc]] <- dt[[cc]]
    for (cc in intersect(c("P_HET", "I2", "Q_HET", "STAT_JOINT", "DF_JOINT"), nm)) out[[cc]] <- dt[[cc]]

    setcolorder(out, c("CHR", "POS", "ID", "P", "LOG10P", "CHR_RAW"))
    if (drop_invalid) out <- out[!is.na(P) & P > 0 & P <= 1 & !is.na(CHR) & !is.na(POS)]
    setattr(out, "source_format", fmt)
    out
}

ancestries_in_sumstats <- function(x) sub("^BETA_", "", grep("^BETA_[A-Za-z0-9]+$", names(x), value = TRUE))

read_lambda <- function(path) {
    # LAMBDA_GC writes "lambda_gc<TAB>1.03"; also accept a bare number
    if (is.null(path) || !nzchar(path) || !file.exists(path)) return(NA_real_)
    lines <- readLines(path, warn = FALSE)
    hit <- grep("lambda", lines, ignore.case = TRUE, value = TRUE)
    if (length(hit)) {
        v <- suppressWarnings(as.numeric(regmatches(hit[1], regexpr("[0-9]+\\.?[0-9]*(e[-+]?[0-9]+)?$", hit[1]))))
        if (!is.na(v)) return(v)
    }
    suppressWarnings(as.numeric(trimws(lines[1])))
}

# Standardised, loci-only table from a GWAS_FILTER "significant" file or any sumstats:
# one row per independent locus (greedy, window_kb), lead variant by P
loci_from_sumstats <- function(x, p_threshold = 5e-8, window_kb = 500, fallback_threshold = 1e-5, max_loci = 50) {
    sig <- x[P < p_threshold]
    if (nrow(sig) == 0) sig <- x[P < fallback_threshold]
    sig <- sig[order(P)]
    w <- window_kb * 1000L; loci <- list()
    while (nrow(sig) > 0 && length(loci) < max_loci) {
        lead <- sig[1]
        loci[[length(loci) + 1]] <- lead[, .(CHR, POS, ID, P, start = pmax(1L, POS - w), end = POS + w)]
        sig <- sig[!(CHR == lead$CHR & abs(POS - lead$POS) <= w)]
    }
    rbindlist(loci)
}
