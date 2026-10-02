process MAGMA_GENE {
    tag "${meta.id ?: meta.trait} - ${meta.trait}"
    label 'process_medium'

    conda "bioconda::magma=1.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/magma:1.10--h9f5acd7_0' :
        'quay.io/biocontainers/magma:1.10--h9f5acd7_0' }"

    input:
    tuple val(meta), path(pval_file), path(genes_annot)   // SNP/P file + .genes.annot from MAGMA_ANNOTATE
    path magma_reference, stageAs: 'ref/*'                // <prefix>.bed/bim/fam LD reference (ancestry-matched)
    val n_samples

    output:
    tuple val(meta), path("${prefix}.genes.out"), emit: genes
    tuple val(meta), path("${prefix}.genes.raw"), emit: raw
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id ?: meta.trait}.${meta.trait}"
    def n = meta.n ?: n_samples
    """
    BFILE=\$(ls ref/*.bed | head -1 | sed 's/\\.bed\$//')

    magma \\
        --bfile \$BFILE \\
        --pval ${pval_file} N=${n} \\
        --gene-annot ${genes_annot} \\
        --out ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        magma: \$(magma --version 2>&1 | head -1 | sed 's/MAGMA version: //')
    END_VERSIONS
    """
}
