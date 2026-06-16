// QTL Dataset Setup and Harmonization Module
// Downloads/prepares QTL datasets from curated_qtl_sources.yml

process QTL_DOWNLOAD {
    tag "qtl_download"
    label 'process_medium'

    conda "conda-forge::r-yaml conda-forge::r-data.table conda-forge::r-httr"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/qtl-tools:latest' :
        'your-registry/qtl-tools:latest' }"

    input:
    path config
    val qtl_types
    val tissues
    val priority_leukemia

    output:
    path "qtl_data", emit: qtl_dir
    path "qtl_data/metadata/manifest.yml", emit: manifest
    path "qtl_data/metadata/download_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def tissues_arg = tissues ? "--tissues ${tissues}" : ''
    """
    download_qtl_datasets.R \\
        --config ${config} \\
        --output_dir qtl_data \\
        --qtl_types ${qtl_types} \\
        ${tissues_arg} \\
        --priority_leukemia ${priority_leukemia} \\
        --public_only true \\
        --genome_build GRCh38 \\
        --verbose \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r_yaml: \$(Rscript -e "cat(as.character(packageVersion('yaml')))" 2>/dev/null || echo "unknown")
        r_httr: \$(Rscript -e "cat(as.character(packageVersion('httr')))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}

process QTL_HARMONIZE {
    tag "${meta.id} - ${qtl_type}"
    label 'process_medium'

    conda "conda-forge::r-data.table conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/qtl-tools:latest' :
        'your-registry/qtl-tools:latest' }"

    input:
    tuple val(meta), path(qtl_file)
    val qtl_type
    val source_format
    val target_build

    output:
    tuple val(meta), path("${prefix}.harmonized.tsv.gz"), emit: harmonized
    tuple val(meta), path("${prefix}.harmonized.tbi"), emit: index, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${qtl_type}"
    """
    harmonize_qtl.R \\
        --input ${qtl_file} \\
        --format ${source_format} \\
        --qtl_type ${qtl_type} \\
        --output ${prefix}.harmonized.tsv.gz \\
        --target_build ${target_build} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        harmonize_qtl: "1.0.0"
    END_VERSIONS
    """
}

process QTL_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::tabix"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://biocontainers/tabix:1.11' :
        'biocontainers/tabix:1.11' }"

    input:
    tuple val(meta), path(qtl_file)

    output:
    tuple val(meta), path(qtl_file), path("*.tbi"), emit: indexed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Sort and index for fast region lookups
    zcat ${qtl_file} | sort -k1,1 -k2,2n | bgzip -c > ${qtl_file.baseName}.sorted.gz
    tabix -s 1 -b 2 -e 2 ${qtl_file.baseName}.sorted.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(tabix --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """
}

process QTL_MERGE_ANCESTRY {
    tag "${meta.id} - ${qtl_type}"
    label 'process_medium'

    // Merge QTL datasets across ancestries for multi-ancestry colocalization
    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/qtl-tools:latest' :
        'your-registry/qtl-tools:latest' }"

    input:
    tuple val(meta), path(qtl_files)
    val qtl_type

    output:
    tuple val(meta), path("${prefix}.merged.tsv.gz"), emit: merged
    tuple val(meta), path("${prefix}.ancestry_summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.${qtl_type}"
    """
    #!/usr/bin/env Rscript

    library(data.table)

    cat("Merging QTL datasets across ancestries\\n")

    # Load all QTL files
    files <- strsplit("${qtl_files}", " ")[[1]]
    cat("Input files:", length(files), "\\n")

    all_qtl <- rbindlist(lapply(files, function(f) {
        dt <- fread(f)
        # Extract ancestry from filename if present
        if (grepl("EUR|AFR|AMR|EAS|SAS|HISP", f)) {
            anc <- regmatches(f, regexpr("EUR|AFR|AMR|EAS|SAS|HISP", f))
            dt[, ancestry := anc]
        } else {
            dt[, ancestry := "UNKNOWN"]
        }
        return(dt)
    }), fill = TRUE)

    cat("Total variants:", nrow(all_qtl), "\\n")
    cat("Ancestries:", paste(unique(all_qtl\$ancestry), collapse = ", "), "\\n")

    # Write merged file
    fwrite(all_qtl, "${prefix}.merged.tsv.gz", sep = "\\t", compress = "gzip")

    # Write ancestry summary
    summary_dt <- all_qtl[, .(
        n_variants = .N,
        n_genes = length(unique(gene_id))
    ), by = ancestry]

    fwrite(summary_dt, "${prefix}.ancestry_summary.tsv", sep = "\\t")

    writeLines(c(
        '"${task.process}":',
        '    qtl_merge: "1.0.0"'
    ), "versions.yml")
    """
}

process QTL_FILTER_REGION {
    tag "${meta.id} - ${region}"
    label 'process_low'

    conda "bioconda::tabix conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/qtl-tools:latest' :
        'your-registry/qtl-tools:latest' }"

    input:
    tuple val(meta), path(qtl_file), path(qtl_index)
    val region  // chr:start-end format

    output:
    tuple val(meta), path("${prefix}.region.tsv.gz"), emit: filtered
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}.${region.replaceAll(':', '_').replaceAll('-', '_')}"
    """
    # Extract region using tabix
    tabix ${qtl_file} ${region} | gzip -c > ${prefix}.region.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(tabix --version 2>&1 | head -1 | awk '{print \$2}')
    END_VERSIONS
    """
}
