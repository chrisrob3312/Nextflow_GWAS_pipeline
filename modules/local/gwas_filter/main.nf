process GWAS_FILTER {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.0' :
        'quay.io/biocontainers/pandas:2.0.0' }"

    input:
    tuple val(meta), path(sumstats)

    output:
    tuple val(meta), path("${prefix}.filtered.tsv.gz"), emit: filtered
    tuple val(meta), path("${prefix}.significant.tsv"), emit: significant
    tuple val(meta), path("${prefix}.loci.tsv"),        emit: loci      // distance-based, not LD-based (500 kb)
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    """
    filter_gwas_results.py \\
        --input ${sumstats} \\
        --output ${prefix}.filtered.tsv.gz \\
        --significant ${prefix}.significant.tsv \\
        --loci ${prefix}.loci.tsv \\
        --p-threshold 5e-8 \\
        --suggestive-threshold 1e-5 \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
