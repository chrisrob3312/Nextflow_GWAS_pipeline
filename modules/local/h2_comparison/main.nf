// Combine heritability and genetic-correlation results across strata
process H2_COMPARISON {
    tag "h2_comparison"
    label 'process_single'

    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-r:latest' :
        'your-registry/gwas-r:latest' }"

    input:
    path h2_files, stageAs: 'h2/*'
    path rg_files, stageAs: 'rg/*'

    output:
    path "h2_comparison.tsv", emit: comparison
    path "rg_comparison.tsv", emit: rg_table, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)
    read_any <- function(dir) {
        fs <- list.files(dir, full.names = TRUE)
        rbindlist(lapply(fs, function(f) { x <- tryCatch(fread(f), error = function(e) NULL); if (is.null(x)) return(NULL); x[, source_file := basename(f)]; x }), fill = TRUE)
    }
    h2 <- read_any("h2"); if (nrow(h2) == 0) h2 <- data.table(note = "no h2 results")
    fwrite(h2, "h2_comparison.tsv", sep = "\\t")
    rg <- read_any("rg"); if (nrow(rg) > 0) fwrite(rg, "rg_comparison.tsv", sep = "\\t")
    cat("h2 rows:", nrow(h2), "| rg rows:", nrow(rg), "\\n")
    writeLines(c('"${task.process}":', '    h2_comparison: "1.0.0"'), "versions.yml")
    """
}
