// Genomic inflation factor per GWAS
process LAMBDA_GC {
    tag "${meta.id} - ${meta.trait}"
    label 'process_single'

    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-r:latest' :
        'your-registry/gwas-r:latest' }"

    input:
    tuple val(meta), path(sumstats)

    output:
    tuple val(meta), path("${prefix}.lambda.txt"), emit: lambda
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    #!/usr/bin/env Rscript
    library(data.table)
    ss <- fread("${sumstats}")
    pc <- intersect(c("P_JOINT","P","PVAL","pvalue","p.value","P_META"), names(ss))[1]
    lp <- intersect(c("LOG10P"), names(ss))[1]
    p <- if (!is.na(pc)) ss[[pc]] else if (!is.na(lp)) 10^(-ss[[lp]]) else stop("No P column")
    p <- p[!is.na(p) & p > 0 & p <= 1]
    chisq <- qchisq(p, df = 1, lower.tail = FALSE)
    lambda <- median(chisq) / qchisq(0.5, df = 1)
    writeLines(c(paste0("lambda_gc\\t", signif(lambda, 4)), paste0("n_variants\\t", length(p))), "${prefix}.lambda.txt")
    cat("lambda GC:", signif(lambda, 4), "\\n")
    writeLines(c('"${task.process}":', '    lambda_gc: "1.0.0"'), "versions.yml")
    """
}
