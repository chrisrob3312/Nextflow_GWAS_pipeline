// Ancestry calls -> per-ancestry genotype sets
//
// ANCESTRY CALLS FILE CONTRACT (GRAF output or user-supplied):
//   tab- or comma-separated, header, a sample ID column (sample_id | IID | ID)
//   and an ancestry column (ancestry | GRAF_ANC | population | ancestry_group).
//   GRAF-ANC numeric codes are mapped to labels (107 AAC, 100 AFR, 300 EUR,
//   500 EAS, 400 SAS, 600 AMR, 601 LAT1, 602 LAT2); labels pass through.
// ANCESTRY_CALLS_PREPARE normalises to:  sample_id <TAB> ancestry

process ANCESTRY_CALLS_PREPARE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-r:latest' :
        'your-registry/gwas-r:latest' }"

    input:
    tuple val(meta), path(calls)

    output:
    tuple val(meta), path("${prefix}.ancestry.txt"), emit: ancestry_calls
    tuple val(meta), path("${prefix}.ancestry_counts.tsv"), emit: counts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    library(data.table)
    d <- fread("${calls}")
    idc <- intersect(c("sample_id", "IID", "ID", "sample", "Sample", "SampleID"), names(d))[1]
    anc <- intersect(c("ancestry", "GRAF_ANC", "Ancestry", "population", "ancestry_group", "AncGroup", "group"), names(d))[1]
    if (is.na(idc) || is.na(anc)) {
        if (ncol(d) >= 2) { idc <- names(d)[1]; anc <- names(d)[2] } else stop("Need sample ID and ancestry columns")
    }
    codes <- c("107" = "AAC", "100" = "AFR", "300" = "EUR", "500" = "EAS", "400" = "SAS",
               "600" = "AMR", "601" = "LAT1", "602" = "LAT2", "200" = "MEN", "700" = "OCN", "800" = "MIX")
    a <- as.character(d[[anc]])
    a <- ifelse(a %in% names(codes), codes[a], a)
    a <- gsub("[^A-Za-z0-9_]", "_", a)
    a[is.na(a) | a == ""] <- "UNKNOWN"
    out <- data.table(sample_id = as.character(d[[idc]]), ancestry = a)
    fwrite(out, "${prefix}.ancestry.txt", sep = "\\t")
    cnt <- out[, .N, by = ancestry][order(-N)]
    fwrite(cnt, "${prefix}.ancestry_counts.tsv", sep = "\\t")
    cat("Ancestry groups:\\n"); print(cnt)
    writeLines(c('"${task.process}":', paste0('    r-data.table: "', packageVersion("data.table"), '"')), "versions.yml")
    """
}

process ANCESTRY_STRATIFY {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'quay.io/biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam), path(ancestry_calls)
    val min_n

    output:
    // One bed/bim/fam per ancestry group named <prefix>.<ancestry>.{bed,bim,fam};
    // the subworkflow splits these into one [meta+ancestry, bed, bim, fam] per group
    tuple val(meta), path("${prefix}.*.bed"), path("${prefix}.*.bim"), path("${prefix}.*.fam"), emit: stratified_genotypes
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # ancestry_calls: sample_id <TAB> ancestry (header)
    tail -n +2 ${ancestry_calls} | cut -f2 | sort -u > ancestry_groups.txt

    while read ancestry; do
        awk -v anc="\$ancestry" 'NR > 1 && \$2 == anc {print \$1, \$1}' ${ancestry_calls} > \${ancestry}_samples.txt
        n=\$(wc -l < \${ancestry}_samples.txt)
        if [ "\$n" -lt ${min_n} ]; then
            echo "Skipping \$ancestry: N=\$n < ${min_n}"
            continue
        fi
        plink2 --bfile ${bed.baseName} --keep \${ancestry}_samples.txt --make-bed \\
               --out ${prefix}.\${ancestry} --threads ${task.cpus} ${args}
    done < ancestry_groups.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed -n '1p' | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
