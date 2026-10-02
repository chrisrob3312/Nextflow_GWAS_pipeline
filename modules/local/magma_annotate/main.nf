// MAGMA annotation (SNP -> gene) and functional summary

process MAGMA_ANNOTATE {
    tag "${meta.id ?: meta.trait} - ${meta.trait}"
    label 'process_low'

    conda "bioconda::magma=1.10 conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magma:1.10--h9f5acd7_0' :
        'quay.io/biocontainers/magma:1.10--h9f5acd7_0' }"

    input:
    tuple val(meta), path(sumstats)
    path gene_loc           // MAGMA gene location file (NCBI build matching genome_build)
    val window              // e.g. '35,10' (kb upstream, downstream)

    output:
    tuple val(meta), path("${prefix}.pval.txt"), path("${prefix}.genes.annot"), emit: annotated
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id ?: meta.trait}.${meta.trait}"
    """
    # SNP location (SNP CHR BP) and p-value (SNP P) files from any sumstats layout
    Rscript - <<'EOF'
    library(data.table)
    ss <- fread("${sumstats}")
    fc <- function(c) { h <- intersect(c, names(ss)); if (length(h)) h[1] else NA }
    idc <- fc(c("ID","SNP","variant_id","rsid")); chr <- fc(c("CHR","CHROM","#CHROM"))
    pos <- fc(c("POS","GENPOS","BP")); pc <- fc(c("P_JOINT","P","PVAL","pvalue","P_META")); lp <- fc(c("LOG10P"))
    if (any(is.na(c(idc, chr, pos)))) stop("Need ID, CHR, POS columns")
    p <- if (!is.na(pc)) ss[[pc]] else if (!is.na(lp)) 10^(-ss[[lp]]) else stop("No P column")
    ok <- !is.na(p)
    fwrite(data.table(SNP = ss[[idc]][ok], CHR = sub("^chr", "", as.character(ss[[chr]][ok])), BP = ss[[pos]][ok]),
           "${prefix}.snploc.txt", sep = "\\t", col.names = FALSE)
    fwrite(data.table(SNP = ss[[idc]][ok], P = p[ok]), "${prefix}.pval.txt", sep = "\\t")
    EOF

    magma --annotate window=${window} \\
          --snp-loc ${prefix}.snploc.txt \\
          --gene-loc ${gene_loc} \\
          --out ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magma: \$(magma --version 2>&1 | head -1 | sed 's/MAGMA version: //')
    END_VERSIONS
    """
}

process FUNCTIONAL_SUMMARY {
    tag "functional_summary"
    label 'process_single'

    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-r:latest' :
        'your-registry/gwas-r:latest' }"

    input:
    path result_files, stageAs: 'results/*'

    output:
    path "functional_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    fs <- list.files("results", full.names = TRUE)
    rows <- lapply(fs, function(f) {
        x <- tryCatch(fread(f), error = function(e) NULL)
        if (is.null(x) || nrow(x) == 0) return(NULL)
        pc <- intersect(c("P", "P_JOINT", "pvalue"), names(x))[1]
        data.table(file = basename(f), n_rows = nrow(x),
                   n_signif_0.05_bonf = if (!is.na(pc)) sum(x[[pc]] < 0.05 / nrow(x), na.rm = TRUE) else NA_integer_,
                   top = if (!is.na(pc) && "GENE" %in% names(x)) x[which.min(x[[pc]])][["GENE"]] else NA_character_)
    })
    out <- rbindlist(rows); if (nrow(out) == 0) out <- data.table(note = "no functional results")
    fwrite(out, "functional_summary.tsv", sep = "\\t")
    writeLines(c('"${task.process}":', '    functional_summary: "1.0.0"'), "versions.yml")
    """
}
