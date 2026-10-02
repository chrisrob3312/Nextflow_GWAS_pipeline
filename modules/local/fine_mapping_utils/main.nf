// Fine-mapping utilities
//   EXTRACT_REGIONS : windows around GWAS signals, one sumstats slice per region
//   SUSIE_RSS       : SuSiE from summary statistics + LD computed from a PLINK
//                     reference (cohort-specific or ancestry-matched panel)
//   CREDIBLE_SETS   : merge credible-set tables from any fine-mapping method

process EXTRACT_REGIONS {
    tag "${meta.id} - ${meta.trait}"
    label 'process_low'

    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-r:latest' :
        'your-registry/gwas-r:latest' }"

    input:
    tuple val(meta), path(sumstats)
    val window_kb
    val p_threshold

    output:
    tuple val(meta), path("${prefix}.regions.tsv"),            emit: regions
    tuple val(meta), path("${prefix}.region_*.sumstats.tsv.gz"), emit: region_sumstats, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    #!/usr/bin/env Rscript
    library(data.table)
    ss <- fread("${sumstats}")
    fc <- function(c) { h <- intersect(c, names(ss)); if (length(h)) h[1] else NA }
    chr <- fc(c("CHR","CHROM","#CHROM","chr")); pos <- fc(c("POS","GENPOS","BP","pos"))
    pc  <- fc(c("P_JOINT","P","PVAL","pvalue","p.value","P_META")); lp <- fc(c("LOG10P"))
    if (is.na(chr) || is.na(pos)) stop("No CHR/POS columns")
    p <- if (!is.na(pc)) ss[[pc]] else if (!is.na(lp)) 10^(-ss[[lp]]) else stop("No P column")
    ss[, CHR_ := sub("^chr", "", as.character(get(chr)))]; ss[, POS_ := as.integer(get(pos))]; ss[, P_ := p]
    sig <- ss[!is.na(P_) & P_ < ${p_threshold}][order(P_)]
    w <- ${window_kb} * 1000L
    regions <- list()
    while (nrow(sig) > 0) {
        lead <- sig[1]
        regions[[length(regions) + 1]] <- data.table(region = length(regions) + 1, CHR = lead\$CHR_,
            start = max(1L, lead\$POS_ - w), end = lead\$POS_ + w, lead_pos = lead\$POS_, lead_p = lead\$P_)
        sig <- sig[!(CHR_ == lead\$CHR_ & abs(POS_ - lead\$POS_) <= w)]
    }
    reg <- rbindlist(regions)
    if (nrow(reg) == 0) { reg <- data.table(region = integer(), CHR = character(), start = integer(), end = integer(), lead_pos = integer(), lead_p = numeric()) }
    fwrite(reg, "${prefix}.regions.tsv", sep = "\\t")
    for (i in seq_len(nrow(reg))) {
        r <- reg[i]
        slice <- ss[CHR_ == r\$CHR & POS_ >= r\$start & POS_ <= r\$end]
        slice[, c("CHR_", "POS_", "P_") := NULL]
        fwrite(slice, sprintf("${prefix}.region_%03d_chr%s_%d_%d.sumstats.tsv.gz", i, r\$CHR, r\$start, r\$end), sep = "\\t", compress = "gzip")
    }
    cat("Regions:", nrow(reg), "\\n")
    writeLines(c('"${task.process}":', '    extract_regions: "1.0.0"'), "versions.yml")
    """
}

process SUSIE_RSS {
    tag "${meta.id} - ${meta.trait}"
    label 'process_medium'

    conda "conda-forge::r-data.table conda-forge::r-susier bioconda::plink"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-r:latest' :
        'your-registry/gwas-r:latest' }"

    input:
    tuple val(meta), path(region_files, stageAs: 'regions/*')
    path ld_reference_dir      // directory with <ancestry>.bed/bim/fam (cohort or panel)
    val max_causal

    output:
    tuple val(meta), path("${prefix}.susie.tsv.gz"),           emit: results
    tuple val(meta), path("${prefix}.susie.credible_sets.tsv"), emit: credible_sets
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    def anc = meta.ancestry ?: 'ALL'
    """
    #!/usr/bin/env Rscript
    library(data.table); library(susieR)
    files <- list.files("regions", full.names = TRUE)
    ld_dir <- "${ld_reference_dir}"
    bfile <- c(file.path(ld_dir, "${anc}"), file.path(ld_dir, "ALL"), file.path(ld_dir, "POOLED"))
    bfile <- bfile[file.exists(paste0(bfile, ".bed"))][1]
    if (is.na(bfile)) stop("No PLINK LD reference for ${anc} in ", ld_dir, " (need <ancestry>.bed/bim/fam or ALL.bed)")
    res <- list(); cs_all <- list()
    for (f in files) {
        ss <- fread(f)
        idc <- intersect(c("ID","SNP","variant_id","rsid"), names(ss))[1]
        bc <- intersect(c("BETA","beta","BETA_EUR"), names(ss))[1]; sec <- intersect(c("SE","se","SE_EUR"), names(ss))[1]
        nc <- intersect(c("N","n"), names(ss))[1]
        if (is.na(idc) || is.na(bc) || is.na(sec)) next
        ss <- ss[!is.na(get(bc)) & !is.na(get(sec)) & get(sec) > 0]
        writeLines(ss[[idc]], "snps.txt")
        system(paste("plink --bfile", bfile, "--extract snps.txt --r square gz --out ld --silent"), ignore.stdout = TRUE)
        if (!file.exists("ld.ld.gz")) next
        R <- as.matrix(fread("ld.ld.gz", header = FALSE))
        # plink writes LD in bim order of the extracted SNPs
        bim <- fread(paste0(bfile, ".bim"), header = FALSE); keep <- bim\$V2[bim\$V2 %in% ss[[idc]]]
        ss <- ss[match(keep, get(idc))]; R[is.na(R)] <- 0; diag(R) <- 1
        z <- ss[[bc]] / ss[[sec]]; n <- if (!is.na(nc)) max(ss[[nc]], na.rm = TRUE) else 1000
        fit <- tryCatch(susie_rss(z = z, R = R, n = n, L = ${max_causal}), error = function(e) NULL)
        if (is.null(fit)) next
        region <- sub("\\\\.sumstats\\\\.tsv\\\\.gz\$", "", basename(f))
        out <- data.table(region = region, ID = ss[[idc]], z = z, PIP = fit\$pip, cs = NA_integer_)
        if (!is.null(fit\$sets\$cs)) for (k in seq_along(fit\$sets\$cs)) out[fit\$sets\$cs[[k]], cs := k]
        res[[region]] <- out
        if (!is.null(fit\$sets\$cs)) cs_all[[region]] <- out[!is.na(cs)][, .(region, cs, ID, PIP)]
    }
    results <- rbindlist(res); fwrite(results, "${prefix}.susie.tsv.gz", sep = "\\t", compress = "gzip")
    cs <- rbindlist(cs_all); if (nrow(cs) == 0) cs <- data.table(region = character(), cs = integer(), ID = character(), PIP = numeric())
    cs[, method := "susie_rss"]; cs[, ancestry := "${anc}"]
    fwrite(cs, "${prefix}.susie.credible_sets.tsv", sep = "\\t")
    cat("Regions fine-mapped:", length(res), "| credible sets:", if (nrow(cs)) length(unique(paste(cs\$region, cs\$cs))) else 0, "\\n")
    writeLines(c('"${task.process}":', paste0('    susieR: "', packageVersion("susieR"), '"')), "versions.yml")
    """
}

process CREDIBLE_SETS {
    tag "${meta.id ?: meta.trait}"
    label 'process_single'

    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-r:latest' :
        'your-registry/gwas-r:latest' }"

    input:
    tuple val(meta), path(cs_files, stageAs: 'cs/*')

    output:
    tuple val(meta), path("${prefix}.credible_sets.tsv"), emit: credible_sets
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id ?: meta.trait}"
    """
    #!/usr/bin/env Rscript
    library(data.table)
    fs <- list.files("cs", full.names = TRUE)
    d <- rbindlist(lapply(fs, function(f) { x <- fread(f); x[, source_file := basename(f)]; x }), fill = TRUE)
    if (nrow(d) == 0) d <- data.table(note = "no credible sets")
    fwrite(d, "${prefix}.credible_sets.tsv", sep = "\\t")
    writeLines(c('"${task.process}":', '    credible_sets: "1.0.0"'), "versions.yml")
    """
}
