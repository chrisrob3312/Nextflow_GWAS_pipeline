// Multi-Method PRS for Admixed Populations
// Primary: PRS-CSx | Adjuncts: GAUDI, DiscoDivas, SDPR_admix, MUSSEL, PROSPER

process PRS_CSX {
    tag "${meta.id}"
    label 'process_high'

    conda "conda-forge::python conda-forge::numpy conda-forge::scipy"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/prscsx:latest' :
        'your-registry/prscsx:latest' }"

    input:
    tuple val(meta), path(sumstats_dir)
    path ld_ref
    tuple val(meta), path(geno_files)
    val ancestries

    output:
    tuple val(meta), path("${prefix}.*.weights.tsv.gz"), emit: weights
    tuple val(meta), path("${prefix}.*.scores.sscore"), emit: scores, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    calculate_prs_admixed.R \\
        --method prs_csx \\
        --sumstats_dir ${sumstats_dir} \\
        --ld_ref ${ld_ref} \\
        --geno ${geno_files[0].baseName} \\
        --ancestries ${ancestries} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prs_csx: "1.0.0"
    END_VERSIONS
    """
}

process PRS_GAUDI {
    tag "${meta.id}"
    label 'process_high'

    // GAUDI: Local ancestry-informed PRS
    conda "conda-forge::r-data.table conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gaudi:latest' :
        'your-registry/gaudi:latest' }"

    input:
    tuple val(meta), path(sumstats)
    tuple val(meta), path(local_ancestry)
    tuple val(meta), path(geno_files)
    val ancestries

    output:
    tuple val(meta), path("${prefix}.gaudi.scores.tsv"), emit: scores
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    calculate_prs_admixed.R \\
        --method gaudi \\
        --sumstats ${sumstats} \\
        --local_ancestry ${local_ancestry} \\
        --geno ${geno_files[0].baseName} \\
        --ancestries ${ancestries} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gaudi: "1.0.0"
    END_VERSIONS
    """
}

process PRS_DISCO_DIVAS {
    tag "${meta.id}"
    label 'process_high'

    // DiscoDivas: Disentangling PRS by local ancestry
    conda "conda-forge::r-data.table conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/discodivas:latest' :
        'your-registry/discodivas:latest' }"

    input:
    tuple val(meta), path(sumstats)
    tuple val(meta), path(local_ancestry)
    tuple val(meta), path(geno_files)
    val ancestries

    output:
    tuple val(meta), path("${prefix}.discodivas.*.scores.tsv"), emit: scores
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    calculate_prs_admixed.R \\
        --method disco_divas \\
        --sumstats ${sumstats} \\
        --local_ancestry ${local_ancestry} \\
        --geno ${geno_files[0].baseName} \\
        --ancestries ${ancestries} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        disco_divas: "1.0.0"
    END_VERSIONS
    """
}

process PRS_SDPR_ADMIX {
    tag "${meta.id}"
    label 'process_high'

    // SDPR_admix: Admixture-aware sparse Dirichlet process regression
    conda "conda-forge::python conda-forge::numpy"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/sdpr:latest' :
        'your-registry/sdpr:latest' }"

    input:
    tuple val(meta), path(sumstats)
    path ld_ref
    tuple val(meta), path(geno_files)
    val ancestries

    output:
    tuple val(meta), path("${prefix}.sdpr.*.weights.tsv"), emit: weights
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    calculate_prs_admixed.R \\
        --method sdpr_admix \\
        --sumstats ${sumstats} \\
        --ld_ref ${ld_ref} \\
        --geno ${geno_files[0].baseName} \\
        --ancestries ${ancestries} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sdpr_admix: "1.0.0"
    END_VERSIONS
    """
}

process PRS_MUSSEL {
    tag "${meta.id}"
    label 'process_high'

    // MUSSEL: Multi-ancestry stacking ensemble learner
    conda "conda-forge::r-data.table conda-forge::r-glmnet"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/mussel:latest' :
        'your-registry/mussel:latest' }"

    input:
    tuple val(meta), path(sumstats_dir)
    path ld_ref
    tuple val(meta), path(geno_files)
    tuple val(meta), path(phenotype)
    val ancestries
    val trait

    output:
    tuple val(meta), path("${prefix}.mussel.scores.tsv"), emit: scores
    tuple val(meta), path("${prefix}.mussel.weights.tsv"), emit: stacking_weights
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    calculate_prs_admixed.R \\
        --method mussel \\
        --sumstats_dir ${sumstats_dir} \\
        --ld_ref ${ld_ref} \\
        --geno ${geno_files[0].baseName} \\
        --ancestries ${ancestries} \\
        --phenotype ${phenotype} \\
        --trait ${trait} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mussel: "1.0.0"
    END_VERSIONS
    """
}

process PRS_PROSPER {
    tag "${meta.id}"
    label 'process_high'

    // PROSPER: Optimal penalized regression
    conda "conda-forge::r-data.table conda-forge::r-glmnet"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/prosper:latest' :
        'your-registry/prosper:latest' }"

    input:
    tuple val(meta), path(sumstats)
    path ld_ref
    tuple val(meta), path(geno_files)
    val ancestries

    output:
    tuple val(meta), path("${prefix}.prosper.weights.tsv"), emit: weights
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    calculate_prs_admixed.R \\
        --method prosper \\
        --sumstats ${sumstats} \\
        --ld_ref ${ld_ref} \\
        --geno ${geno_files[0].baseName} \\
        --ancestries ${ancestries} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prosper: "1.0.0"
    END_VERSIONS
    """
}

process PRS_COMBINE_METHODS {
    tag "${meta.id}"
    label 'process_low'

    // Compare and combine PRS from multiple methods
    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/prs:latest' :
        'your-registry/prs:latest' }"

    input:
    tuple val(meta), path(prs_scores)
    tuple val(meta), path(phenotype)
    val trait

    output:
    tuple val(meta), path("${prefix}.prs_comparison.tsv"), emit: comparison
    tuple val(meta), path("${prefix}.best_prs.tsv"), emit: best_prs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    #!/usr/bin/env Rscript

    library(data.table)

    # Load all PRS scores
    score_files <- list.files(".", pattern = "*.scores.*.tsv|*.sscore", full.names = TRUE)
    all_scores <- lapply(score_files, fread)

    # Load phenotype
    pheno <- fread("${phenotype}")

    # Calculate R² and AUC for each method
    # ... validation code

    cat("PRS comparison - placeholder\\n")

    # Write comparison
    fwrite(data.table(method = score_files), "${prefix}.prs_comparison.tsv", sep = "\\t")
    fwrite(data.table(best = "prs_csx"), "${prefix}.best_prs.tsv", sep = "\\t")

    writeLines(c(
        '"${task.process}":',
        '    prs_combine: "1.0.0"'
    ), "versions.yml")
    """
}
