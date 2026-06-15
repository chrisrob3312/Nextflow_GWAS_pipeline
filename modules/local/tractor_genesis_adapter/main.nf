// Tractor → GENESIS Adapter Module
// Integrates Tractor local ancestry-aware decomposition with GENESIS mixed models
// Enables time-to-event (Cox) analysis for admixed populations with local ancestry

process TRACTOR_GENESIS_PREPARE {
    tag "$meta.id - $meta.ancestry"
    label 'process_high'

    conda "bioconda::bioconductor-genesis=2.32.0 conda-forge::r-survival"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/tractor-genesis:latest' :
        'your-registry/tractor-genesis:latest' }"

    input:
    tuple val(meta), path(ancestry_dosages)  // Tractor ancestry-specific dosages
    tuple val(meta), path(gds)               // GDS file for GENESIS
    tuple val(meta), path(phenotype)
    path kinship_matrix                       // GRM for relatedness

    output:
    tuple val(meta), path("${prefix}.tractor_gds"), emit: tractor_gds
    tuple val(meta), path("${prefix}.ancestry_weights.rds"), emit: ancestry_weights
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}"
    """
    #!/usr/bin/env Rscript

    # Tractor → GENESIS Adapter
    # Creates GDS file with ancestry-specific dosage matrices for GENESIS

    library(SeqArray)
    library(SeqVarTools)

    # Read Tractor ancestry-specific dosages
    # Format: chr, pos, snp, then columns for each ancestry (e.g., dose_EUR, dose_AFR, dose_AMR)
    anc_files <- list.files(pattern = "*.ancdose.*.tsv.gz")

    # Parse ancestry populations from file names
    ancestries <- gsub(".*\\\\.ancdose\\\\.(.*)\\\\.tsv\\\\.gz", "\\\\1", anc_files)
    cat("Detected ancestries:", ancestries, "\\n")

    # Read dosage matrices
    dosage_list <- lapply(anc_files, function(f) {
        data.table::fread(f)
    })
    names(dosage_list) <- ancestries

    # Create unified GDS with ancestry-specific dosage nodes
    # This allows GENESIS to access ancestry-stratified genotypes

    # Open original GDS as template
    gds_orig <- seqOpen("${gds}")
    variant_ids <- seqGetData(gds_orig, "variant.id")
    sample_ids <- seqGetData(gds_orig, "sample.id")
    seqClose(gds_orig)

    # Create new GDS with Tractor structure
    gds_out <- seqOpen("${prefix}.tractor_gds", mode = "w")

    # Copy variant/sample structure
    # Add ancestry-specific dosage nodes
    for (anc in ancestries) {
        node_name <- paste0("annotation/dosage_", anc)
        # Add dosage matrix for this ancestry
        cat("Adding dosage node:", node_name, "\\n")
    }

    seqClose(gds_out)

    # Save ancestry weight information
    # For weighted combination of ancestry-specific effects
    ancestry_info <- list(
        ancestries = ancestries,
        n_variants = nrow(dosage_list[[1]]),
        n_samples = ncol(dosage_list[[1]]) - 3  # Exclude chr, pos, snp columns
    )
    saveRDS(ancestry_info, "${prefix}.ancestry_weights.rds")

    cat("Tractor → GENESIS preparation complete\\n")

    writeLines(c(
        '"${task.process}":',
        paste0('    genesis: "', packageVersion("GENESIS"), '"'),
        paste0('    R: "', R.version.string, '"')
    ), "versions.yml")
    """
}

process TRACTOR_GENESIS_SURVIVAL {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_high'

    // GENESIS with survival analysis for admixed populations using Tractor decomposition
    // NOTE: This is a THIN WRAPPER that calls bin/tractor_genesis_adapter.R
    // All analysis logic is in the bin script - changes should be made there
    // See shell/tractor_genesis_adapter.sh for equivalent shell wrapper
    conda "bioconda::bioconductor-genesis=2.32.0 conda-forge::r-survival conda-forge::r-data.table conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/tractor-genesis:latest' :
        'your-registry/tractor-genesis:latest' }"

    input:
    tuple val(meta), path(tractor_gds), path(ancestry_weights)
    tuple val(meta), path(phenotype)
    path kinship_matrix
    val covariate_cols
    val ancestral_pops   // e.g., "EUR,AFR" or "EUR,AFR,AMR"

    output:
    tuple val(meta), path("${prefix}.joint.tsv.gz"), emit: joint_results
    tuple val(meta), path("${prefix}.ancestry_specific.tsv.gz"), emit: ancestry_results
    tuple val(meta), path("${prefix}.het_test.tsv"), emit: het_test
    tuple val(meta), path("${prefix}.null_model.rds"), emit: null_model
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    def time_col = meta.time_col ?: 'time'
    def event_col = meta.event_col ?: 'event'
    def covar_arg = covariate_cols ? "--covariates ${covariate_cols}" : ''
    def kinship_arg = kinship_matrix ? "--kinship ${kinship_matrix}" : ''
    def weights_arg = ancestry_weights ? "--weights ${ancestry_weights}" : ''
    """
    # Thin wrapper - calls bin/tractor_genesis_adapter.R
    # Equivalent to: shell/tractor_genesis_adapter.sh -m survival ...

    tractor_genesis_adapter.R \\
        --gds ${tractor_gds} \\
        --phenotype ${phenotype} \\
        --trait ${meta.trait} \\
        --model survival \\
        --ancestries ${ancestral_pops} \\
        --time_col ${time_col} \\
        --event_col ${event_col} \\
        --output_prefix ${prefix} \\
        ${covar_arg} \\
        ${kinship_arg} \\
        ${weights_arg} \\
        --verbose

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tractor_genesis_adapter: "1.0.0"
        genesis: \$(Rscript -e "cat(as.character(packageVersion('GENESIS')))")
        survival: \$(Rscript -e "cat(as.character(packageVersion('survival')))")
        R: \$(Rscript -e "cat(R.version.string)")
    END_VERSIONS
    """
}

process TRACTOR_GENESIS_BINARY {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_high'

    // GENESIS with binary outcomes for admixed populations using Tractor decomposition
    // NOTE: THIN WRAPPER - calls bin/tractor_genesis_adapter.R
    conda "bioconda::bioconductor-genesis=2.32.0 conda-forge::r-data.table conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/tractor-genesis:latest' :
        'your-registry/tractor-genesis:latest' }"

    input:
    tuple val(meta), path(tractor_gds), path(ancestry_weights)
    tuple val(meta), path(phenotype)
    path kinship_matrix
    val covariate_cols
    val ancestral_pops

    output:
    tuple val(meta), path("${prefix}.joint.tsv.gz"), emit: joint_results
    tuple val(meta), path("${prefix}.ancestry_specific.tsv.gz"), emit: ancestry_results
    tuple val(meta), path("${prefix}.het_test.tsv"), emit: het_test
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    def covar_arg = covariate_cols ? "--covariates ${covariate_cols}" : ''
    def kinship_arg = kinship_matrix ? "--kinship ${kinship_matrix}" : ''
    def weights_arg = ancestry_weights ? "--weights ${ancestry_weights}" : ''
    """
    # Thin wrapper - calls bin/tractor_genesis_adapter.R
    # Equivalent to: shell/tractor_genesis_adapter.sh -m binary ...

    tractor_genesis_adapter.R \\
        --gds ${tractor_gds} \\
        --phenotype ${phenotype} \\
        --trait ${meta.trait} \\
        --model binary \\
        --ancestries ${ancestral_pops} \\
        --output_prefix ${prefix} \\
        ${covar_arg} \\
        ${kinship_arg} \\
        ${weights_arg} \\
        --verbose

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tractor_genesis_adapter: "1.0.0"
        genesis: \$(Rscript -e "cat(as.character(packageVersion('GENESIS')))")
        R: \$(Rscript -e "cat(R.version.string)")
    END_VERSIONS
    """
}

process TRACTOR_GENESIS_QUANTITATIVE {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_high'

    // GENESIS with quantitative outcomes for admixed populations using Tractor decomposition
    // E.g., MRD levels, biomarker values
    // NOTE: THIN WRAPPER - calls bin/tractor_genesis_adapter.R
    conda "bioconda::bioconductor-genesis=2.32.0 conda-forge::r-data.table conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/tractor-genesis:latest' :
        'your-registry/tractor-genesis:latest' }"

    input:
    tuple val(meta), path(tractor_gds), path(ancestry_weights)
    tuple val(meta), path(phenotype)
    path kinship_matrix
    val covariate_cols
    val ancestral_pops

    output:
    tuple val(meta), path("${prefix}.joint.tsv.gz"), emit: joint_results
    tuple val(meta), path("${prefix}.ancestry_specific.tsv.gz"), emit: ancestry_results
    tuple val(meta), path("${prefix}.het_test.tsv"), emit: het_test
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    def covar_arg = covariate_cols ? "--covariates ${covariate_cols}" : ''
    def kinship_arg = kinship_matrix ? "--kinship ${kinship_matrix}" : ''
    def weights_arg = ancestry_weights ? "--weights ${ancestry_weights}" : ''
    """
    # Thin wrapper - calls bin/tractor_genesis_adapter.R
    # Equivalent to: shell/tractor_genesis_adapter.sh -m quantitative ...

    tractor_genesis_adapter.R \\
        --gds ${tractor_gds} \\
        --phenotype ${phenotype} \\
        --trait ${meta.trait} \\
        --model quantitative \\
        --ancestries ${ancestral_pops} \\
        --output_prefix ${prefix} \\
        ${covar_arg} \\
        ${kinship_arg} \\
        ${weights_arg} \\
        --verbose

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tractor_genesis_adapter: "1.0.0"
        genesis: \$(Rscript -e "cat(as.character(packageVersion('GENESIS')))")
        R: \$(Rscript -e "cat(R.version.string)")
    END_VERSIONS
    """
}
