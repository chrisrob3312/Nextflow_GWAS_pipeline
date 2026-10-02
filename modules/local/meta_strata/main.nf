// Inverse-variance meta-analysis of per-stratum GWAS -> POOLED
// Thin wrapper around bin/meta_analyze_strata.R (same script as SLURM step 2).
// Keeps the Tractor-GENESIS structure: each ancestry background is pooled
// across strata (BETA_<anc>, SE_<anc>, P_<anc>), P_JOINT combined by Fisher.

process META_STRATA {
    tag "${trait_meta.trait}"
    label 'process_medium'

    conda "conda-forge::r-data.table conda-forge::r-optparse conda-forge::r-metafor"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-r:latest' :
        'your-registry/gwas-r:latest' }"

    input:
    tuple val(trait_meta), val(ancestries), path(sumstats_files, stageAs: 'in/*')
    val method   // 'random' or 'fixed'

    output:
    tuple val(trait_meta), path("${prefix}.POOLED.sumstats.gz"), emit: meta_results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${trait_meta.trait}"
    def anc_list = ancestries instanceof List ? ancestries : [ancestries]
    def ss_list  = sumstats_files instanceof List ? sumstats_files : [sumstats_files]
    def stage = [ss_list, anc_list].transpose().collect { f, a -> "ln -sf \$(readlink -f ${f}) strata/${a}.sumstats.gz" }.join('\n    ')
    """
    mkdir -p strata
    ${stage}

    meta_analyze_strata.R \\
        --input_dir strata \\
        --strata ${anc_list.join(',')} \\
        --method ${method} \\
        --min_strata 2 \\
        --output ${prefix}.POOLED.sumstats.gz \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meta_analyze_strata: "1.0.0"
        metafor: \$(Rscript -e "cat(tryCatch(as.character(packageVersion('metafor')), error=function(e) 'not installed'))")
    END_VERSIONS
    """
}
