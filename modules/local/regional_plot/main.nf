process REGIONAL_PLOT {
    tag "$meta.trait - $meta.ancestry"
    label 'process_low'

    // LocusZoom-style regional plots at each independent locus; any sumstats layout
    // (column handling in bin/sumstats_utils.R). Gene track from GTF/GFF/BED when given.
    conda "conda-forge::r-base=4.3 conda-forge::r-data.table conda-forge::r-argparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-plots:latest' :
        'your-registry/gwas-plots:latest' }"

    input:
    tuple val(meta), path(sumstats), path(significant)   // full sumstats + GWAS_FILTER significant file (loci)
    path gene_annotation                                   // GTF / GFF / BED, or []

    output:
    tuple val(meta), path("${prefix}/*.regional.png"), emit: plots, optional: true
    tuple val(meta), path("${prefix}/*.regional.pdf"), emit: pdf,   optional: true
    tuple val(meta), path("${prefix}/*.loci.tsv"),     emit: loci,  optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.trait}.${meta.ancestry ?: 'ALL'}"
    def genes_arg = gene_annotation ? "--genes ${gene_annotation}" : ''
    def loci_arg  = significant ? "--loci ${significant}" : ''
    def per_anc   = (params.plot_per_ancestry == null || params.plot_per_ancestry) ? '--per-ancestry' : ''
    """
    mkdir -p ${prefix}

    create_regional_plots.R \\
        --sumstats ${sumstats} \\
        ${loci_arg} \\
        ${genes_arg} \\
        --output-dir ${prefix} \\
        --ancestry ${meta.ancestry ?: 'ALL'} \\
        --window-kb ${params.regional_window_kb ?: 500} \\
        --p-threshold ${params.p_threshold ?: 5e-8} \\
        --suggestive ${params.suggestive_threshold ?: 1e-5} \\
        --max-loci ${params.regional_max_loci ?: 20} \\
        ${per_anc} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -1 | cut -d' ' -f3)
    END_VERSIONS
    """
}
