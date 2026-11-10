# ============================================================================
# duplicates.R: one member of every duplicate / MZ-twin pair leaves the analysis (CP-13)
# ============================================================================
# Sourced by ingest_structure.R (the bundle's kinship) and genesis_pcair_pcrelate.R (the
# comparator's PC-Relate kinship), so both structure paths apply one rule:
#   a pair with kinship > 2^(-3/2) = 0.3536 (scaleKin = 1; the duplicate / MZ boundary of
#   Manichaikul et al. 2010) is a duplicate. In each group of samples joined by such pairs, the
#   sample with the highest call rate stays; ties go to the first sample ID in sort order. The
#   others are removed from the phenotype the GWAS reads, so step 1 never sees them. The GRM and
#   PC files are left as they are (a sample without a phenotype row is not analysed).
# Writes nothing itself; the callers write <out>/cohort.removed_duplicates.tsv with the bundle's
# header (kept_iid removed_iid kinship reason; header-only when nothing is removed).
# ============================================================================

DUP_KINSHIP <- 2^(-3/2)

# select_duplicates(kin, call_rate, thresh): data.frame(kept_iid, removed_iid, kinship, reason),
# zero rows when no pair is above thresh. kin: square kinship matrix with dimnames (scaleKin = 1).
# call_rate: named numeric vector (sample -> call rate), or NULL when not available.
select_duplicates <- function(kin, call_rate = NULL, thresh = DUP_KINSHIP) {
    out <- data.frame(kept_iid = character(), removed_iid = character(), kinship = numeric(),
                      reason = character(), stringsAsFactors = FALSE)
    k <- as.matrix(kin)
    ids <- rownames(k)
    if (is.null(ids) || !identical(ids, colnames(k))) stop("select_duplicates: the kinship matrix needs matching dimnames")
    hit <- which(k > thresh & upper.tri(k), arr.ind = TRUE)
    if (nrow(hit) == 0) return(out)
    # groups of samples joined by duplicate pairs (union-find over the pairs)
    parent <- seq_along(ids)
    find <- function(i) { while (parent[i] != i) i <- parent[i]; i }
    for (r in seq_len(nrow(hit))) {
        a <- find(hit[r, 1]); b <- find(hit[r, 2])
        if (a != b) parent[max(a, b)] <- min(a, b)
    }
    involved <- sort(unique(as.vector(hit)))
    group <- vapply(involved, find, 1L)
    cr <- setNames(rep(NA_real_, length(ids)), ids)
    if (!is.null(call_rate)) cr[intersect(names(call_rate), ids)] <- call_rate[intersect(names(call_rate), ids)]
    for (g in unique(group)) {
        mem <- ids[involved[group == g]]
        crm <- ifelse(is.na(cr[mem]), -Inf, cr[mem])
        ord <- order(-crm, mem)
        keep <- mem[ord[1]]
        for (rm in mem[ord[-1]]) {
            why <- if (is.na(cr[keep]) || is.na(cr[rm])) "call rate not available; kept the first ID"
                   else if (cr[keep] > cr[rm]) sprintf("kept the higher call rate (%.4f vs %.4f)", cr[keep], cr[rm])
                   else sprintf("equal call rate (%.4f); kept the first ID", cr[keep])
            # kinship: with the kept sample, or, in a group of three or more not all paired with
            # it, the largest kinship with any group member (named)
            others <- setdiff(mem, rm); partner <- others[which.max(k[rm, others])]
            via <- if (k[keep, rm] > thresh) "" else sprintf("; grouped with %s through %s", keep, partner)
            out[nrow(out) + 1, ] <- list(keep, rm, max(k[keep, rm], k[rm, partner]),
                                         sprintf("kinship > %.4f (duplicate / MZ twin); %s%s", thresh, why, via))
        }
    }
    out
}

# bed_call_rate(prefix, ids, max_snps): call rate of each of ids from a PLINK 1 .bed (SNP-major),
# over up to max_snps variants read in 10 evenly spaced blocks (enough to break a tie between two
# samples; the whole file of an imputed cohort would be gigabytes). Base R only.
bed_call_rate <- function(prefix, ids, max_snps = 50000L, n_blocks = 10L) {
    bed <- paste0(prefix, ".bed")
    fam <- read.table(paste0(prefix, ".fam"), colClasses = "character", header = FALSE)[[2]]
    idx <- match(ids, fam)
    if (anyNA(idx)) stop("bed_call_rate: not in ", prefix, ".fam: ", paste(ids[is.na(idx)], collapse = ", "))
    nb <- ceiling(length(fam) / 4)
    m <- (file.size(bed) - 3) / nb
    if (m < 1 || m != round(m)) stop("bed_call_rate: ", bed, " does not have the size of a SNP-major .bed for this .fam")
    con <- file(bed, "rb"); on.exit(close(con))
    if (!identical(readBin(con, "raw", 3), as.raw(c(0x6c, 0x1b, 0x01)))) stop("bed_call_rate: ", bed, " is not a SNP-major PLINK 1 .bed")
    block <- max(1, ceiling(min(max_snps, m) / n_blocks))
    starts <- unique(floor(seq(0, m - block, length.out = min(n_blocks, m))))
    byte_i <- (idx - 1) %/% 4 + 1
    shift <- ((idx - 1) %% 4) * 2
    missing <- numeric(length(ids)); total <- 0
    for (s in starts) {
        seek(con, 3 + s * nb)
        r <- matrix(as.integer(readBin(con, "raw", block * nb)), nrow = nb)[byte_i, , drop = FALSE]
        codes <- matrix(bitwAnd(bitwShiftR(r, shift), 3L), nrow = nrow(r))    # 2 bits per sample; 01 = missing
        missing <- missing + rowSums(codes == 1L)
        total <- total + ncol(r)
    }
    setNames(1 - missing / total, ids)
}
