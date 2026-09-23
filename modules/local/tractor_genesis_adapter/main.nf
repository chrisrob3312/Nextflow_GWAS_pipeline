// Tractor -> GENESIS Adapter Module
// Local ancestry-aware association with kinship, for ALL trait types through
// ONE conditional score-test engine so results are comparable across traits:
//   binary        (MRD)          : GENESIS logistic mixed-model null
//   quantitative                 : GENESIS linear mixed-model null
//   survival      (OS, relapse)  : Cox null (coxme frailty on kinship) +
//                                  martingale-residual score test
// Thin wrapper - all logic lives in bin/tractor_genesis_adapter.R

process TRACTOR_GENESIS {
    tag "$meta.id - $meta.ancestry - $meta.trait ($model)"
    label 'process_high'

    conda "bioconda::bioconductor-genesis=2.32.0 bioconda::bioconductor-seqarray bioconda::bioconductor-seqvartools conda-forge::r-survival conda-forge::r-coxme conda-forge::r-matrix conda-forge::r-data.table conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/tractor-genesis:latest' :
        'your-registry/tractor-genesis:latest' }"

    input:
    tuple val(meta), path(ancestry_dosages), path(haplotype_counts), path(phenotype)
    path kinship_matrix
    val covariate_cols

    output:
    tuple val(meta), path("${prefix}.tractor_genesis.tsv.gz"),               emit: sumstats
    tuple val(meta), path("${prefix}.tractor_genesis.genomic_order.tsv.gz"), emit: genomic_order
    tuple val(meta), path("${prefix}.tractor_genesis.top_hits.tsv"),         emit: top_hits
    tuple val(meta), path("${prefix}.tractor_genesis.significant.tsv"),      emit: significant, optional: true
    tuple val(meta), path("${prefix}.tractor_genesis.heterogeneous.tsv"),    emit: heterogeneous, optional: true
    tuple val(meta), path("${prefix}.tractor_genesis.*_specific.tsv"),       emit: ancestry_specific, optional: true
    tuple val(meta), path("${prefix}.null_model.rds"),                       emit: null_model
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    model = meta.model ?: (meta.survival ? 'survival' : (meta.binary ? 'binary' : 'quantitative'))
    def surv_args    = model == 'survival' ? "--time_col ${meta.time_col} --event_col ${meta.event_col}" : ''
    def covar_arg    = covariate_cols ? "--covariates ${covariate_cols}" : ''
    def kinship_arg  = kinship_matrix ? "--kinship ${kinship_matrix}" : ''
    def ref_arg      = meta.ref_ancestry ? "--ref_ancestry ${meta.ref_ancestry}" : ''
    // TRACTOR_EXTRACT_TRACTS writes <id>.ancdose.<anc>.tsv.gz / <id>.hapcount.<anc>.tsv.gz
    def tractor_prefix = ancestry_dosages[0].name.replaceAll(/\.ancdose\..*$/, '')
    """
    tractor_genesis_adapter.R \\
        --tractor_prefix ${tractor_prefix} \\
        --phenotype ${phenotype} \\
        --trait ${meta.trait} \\
        --model ${model} \\
        ${surv_args} \\
        --ancestries ${meta.tractor_pops} \\
        ${ref_arg} \\
        ${covar_arg} \\
        ${kinship_arg} \\
        --mac_min ${params.tractor_mac_min ?: 10} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tractor_genesis_adapter: "2.0.0"
        genesis: \$(Rscript -e "cat(as.character(packageVersion('GENESIS')))")
        survival: \$(Rscript -e "cat(as.character(packageVersion('survival')))")
        coxme: \$(Rscript -e "cat(tryCatch(as.character(packageVersion('coxme')), error=function(e) 'not installed'))")
        R: \$(Rscript -e "cat(R.version.string)")
    END_VERSIONS
    """
}
