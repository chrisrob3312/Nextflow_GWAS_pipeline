// QTL Harmonization Module
// Harmonizes QTL summary statistics from multiple sources into unified format
// Supports liftover, deduplication, filtering by tissue/ancestry/type

process QTL_HARMONIZE {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::r-data.table conda-forge::r-optparse bioconda::bioconductor-rtracklayer"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/qtl-harmonize:latest' :
        'your-registry/qtl-harmonize:latest' }"

    input:
    tuple val(meta), path(qtl_files)
    path manifest                    // Optional: manifest with file metadata
    path chain_file                  // Optional: for liftover
    val output_build                 // GRCh37 or GRCh38

    output:
    tuple val(meta), path("${prefix}.harmonized.tsv.gz"), emit: harmonized
    tuple val(meta), path("${prefix}.summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def manifest_arg = manifest ? "--manifest ${manifest}" : "--input ."
    def chain_arg = chain_file ? "--chain_file ${chain_file}" : ''
    def input_build = meta.input_build ?: 'GRCh38'
    """
    # Link input files if they exist
    for f in ${qtl_files}; do
        ln -sf \$f . 2>/dev/null || true
    done

    harmonize_qtl.R \\
        ${manifest_arg} \\
        --input_build ${input_build} \\
        --output_build ${output_build} \\
        ${chain_arg} \\
        --output ${prefix} \\
        --dedup_strategy best_p \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        harmonize_qtl: "1.0.0"
        R: \$(Rscript -e "cat(R.version.string)")
    END_VERSIONS
    """
}

process QTL_HARMONIZE_BY_TYPE {
    tag "${meta.id} - ${qtl_type}"
    label 'process_low'

    conda "conda-forge::r-data.table conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/qtl-harmonize:latest' :
        'your-registry/qtl-harmonize:latest' }"

    input:
    tuple val(meta), path(qtl_file)
    val qtl_type                     // eqtl, sqtl, pqtl, etc.
    val format                       // Input format

    output:
    tuple val(meta), path("${prefix}.${qtl_type}.harmonized.tsv.gz"), emit: harmonized
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    harmonize_qtl.R \\
        --input ${qtl_file} \\
        --format ${format} \\
        --qtl_types ${qtl_type} \\
        --output ${prefix}.${qtl_type} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        harmonize_qtl: "1.0.0"
    END_VERSIONS
    """
}

process QTL_FILTER_REGION {
    tag "${meta.id} - ${region}"
    label 'process_low'

    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/qtl-harmonize:latest' :
        'your-registry/qtl-harmonize:latest' }"

    input:
    tuple val(meta), path(harmonized_qtl)
    val region                       // chr:start-end

    output:
    tuple val(meta), path("${prefix}.region.tsv.gz"), emit: filtered
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    harmonize_qtl.R \\
        --input ${harmonized_qtl} \\
        --region ${region} \\
        --output ${prefix}.region \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        harmonize_qtl: "1.0.0"
    END_VERSIONS
    """
}
