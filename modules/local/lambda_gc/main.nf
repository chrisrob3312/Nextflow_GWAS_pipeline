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
    # Column handling for every GWAS layout lives in bin/sumstats_utils.R
    source(Sys.which("sumstats_utils.R"))
    ss <- read_sumstats_std("${sumstats}")
    p <- ss\$P
    chisq <- qchisq(p, df = 1, lower.tail = FALSE)
    lambda <- median(chisq) / qchisq(0.5, df = 1)
    writeLines(c(paste0("lambda_gc\\t", signif(lambda, 4)), paste0("n_variants\\t", length(p)),
                 paste0("source_format\\t", attr(ss, "source_format"))), "${prefix}.lambda.txt")
    cat("lambda GC:", signif(lambda, 4), "(", attr(ss, "source_format"), ")\\n")
    writeLines(c('"${task.process}":', '    lambda_gc: "1.0.0"'), "versions.yml")
    """
}
