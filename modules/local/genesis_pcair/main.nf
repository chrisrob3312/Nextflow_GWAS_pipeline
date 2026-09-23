// GENESIS PC-AiR + PC-Relate
// Ancestry PCs robust to relatedness (PC-AiR) and a kinship / GRM robust to
// admixture (PC-Relate). Run ONCE on the QC'd full cohort; the outputs feed
// every downstream model (standard GENESIS GWAS, Tractor-GENESIS, GxG, PRS
// validation) so all of them adjust for the same structure.

process GENESIS_PCAIR_PCRELATE {
    tag "${meta.id}"
    label 'process_high'

    conda "bioconda::bioconductor-genesis=2.32.0 bioconda::bioconductor-snprelate bioconda::bioconductor-gwastools conda-forge::r-data.table conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/tractor-genesis:latest' :
        'your-registry/tractor-genesis:latest' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path phenotype
    val n_pcs

    output:
    tuple val(meta), path("${prefix}.pcrelate.grm.rds"),           emit: grm
    tuple val(meta), path("${prefix}.pcrelate.kinship.rds"),       emit: kinship
    tuple val(meta), path("${prefix}.pcair.pcs.tsv"),              emit: pcs
    tuple val(meta), path("${prefix}.phenotypes.with_pcs.tsv"),    emit: phenotype_with_pcs
    tuple val(meta), path("${prefix}.unrelated.txt"),              emit: unrelated
    path "${prefix}.pcair.variance.tsv",                            emit: variance
    path "${prefix}.pcrelate.relatedness_summary.tsv",              emit: relatedness
    path "versions.yml",                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    genesis_pcair_pcrelate.R \\
        --geno ${bed.baseName} \\
        --phenotype ${phenotype} \\
        --n_pcs ${n_pcs} \\
        --n_pcs_pcrelate ${params.pcair_n_pcs_pcrelate ?: 5} \\
        --kin_thresh ${params.pcair_kin_thresh ?: 0.0442} \\
        --ld_r2 ${params.pcair_ld_r2 ?: 0.1} \\
        --iterations ${params.pcair_iterations ?: 2} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genesis: \$(Rscript -e "cat(as.character(packageVersion('GENESIS')))")
        snprelate: \$(Rscript -e "cat(as.character(packageVersion('SNPRelate')))")
    END_VERSIONS
    """
}
