process MANHATTAN_PLOT {
    tag "$meta.trait - $meta.ancestry"
    label 'process_low'

    // Column handling (Tractor-GENESIS, REGENIE, SAIGE, GENESIS, PLINK2, ...) is in bin/sumstats_utils.R
    conda "conda-forge::r-base=4.3 conda-forge::r-data.table conda-forge::r-argparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-plots:latest' :
        'your-registry/gwas-plots:latest' }"

    input:
    tuple val(meta), path(sumstats), path(lambda_gc)

    output:
    tuple val(meta), path("${prefix}.manhattan.png"), emit: plots
    tuple val(meta), path("${prefix}.manhattan.pdf"), emit: pdf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.trait}.${meta.ancestry ?: 'ALL'}"
    def per_anc = (params.plot_per_ancestry == null || params.plot_per_ancestry) ? '--per-ancestry' : ''
    """
    create_manhattan_plot.R \\
        --sumstats ${sumstats} \\
        --lambda-gc ${lambda_gc} \\
        --output-prefix ${prefix} \\
        --title "${meta.trait} - ${meta.ancestry ?: 'ALL'}" \\
        --genome-wide ${params.p_threshold ?: 5e-8} \\
        --suggestive ${params.suggestive_threshold ?: 1e-5} \\
        ${per_anc} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """
}
