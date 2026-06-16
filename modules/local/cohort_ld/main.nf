// Cohort-Specific LD Calculation Module
// Calculates LD matrices from cohort data, stratified by user-defined ancestry groups
//
// FLEXIBLE DESIGN:
// - Not hardcoded to specific ancestry inference tools (GRAF-ANC, ADMIXTURE, etc.)
// - Users define strata via ancestry_file with sample→group mapping
// - Can combine groups (e.g., LAT1+LAT2 → LATINO)
// - Supports any grouping variable (ancestry, site, cohort, etc.)

process COHORT_LD_CALCULATE {
    tag "$meta.id - ${ancestry_groups}"
    label 'process_high'

    // NOTE: THIN WRAPPER - calls bin/calculate_cohort_ld.R
    conda "bioconda::bioconductor-seqarray conda-forge::r-data.table conda-forge::r-matrix conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/cohort-ld:latest' :
        'your-registry/cohort-ld:latest' }"

    input:
    tuple val(meta), path(gds)
    path ancestry_file           // TSV with sample_id, ancestry columns
    val ancestry_col             // Column name in ancestry_file
    val ancestry_groups          // Comma-separated groups to process (or "ALL")
    val combine_groups           // Optional: "NEW=OLD1+OLD2" format

    output:
    tuple val(meta), path("${prefix}.*.ld.*"), emit: ld_files
    tuple val(meta), path("${prefix}.*.variants.tsv"), emit: variant_info
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def anc_file_arg = ancestry_file ? "--ancestry_file ${ancestry_file}" : ''
    def anc_col_arg = ancestry_col ? "--ancestry_col ${ancestry_col}" : ''
    def groups_arg = ancestry_groups ? "--groups ${ancestry_groups}" : ''
    def combine_arg = combine_groups ? "--combine_groups '${combine_groups}'" : ''
    def window_kb = task.ext.window_kb ?: 1000
    def format = task.ext.ld_format ?: 'ldstore'
    """
    calculate_cohort_ld.R \\
        --gds ${gds} \\
        ${anc_file_arg} \\
        ${anc_col_arg} \\
        ${groups_arg} \\
        ${combine_arg} \\
        --window_kb ${window_kb} \\
        --format ${format} \\
        --threads ${task.cpus} \\
        --output_prefix ${prefix} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cohort_ld: "1.0.0"
        R: \$(Rscript -e "cat(R.version.string)")
    END_VERSIONS
    """
}

process COHORT_LD_REGION {
    tag "$meta.id - $region"
    label 'process_medium'

    // Calculate LD for a specific region (used for fine-mapping, coloc)
    conda "bioconda::bioconductor-seqarray conda-forge::r-data.table conda-forge::r-matrix conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/cohort-ld:latest' :
        'your-registry/cohort-ld:latest' }"

    input:
    tuple val(meta), path(gds)
    path ancestry_file
    val ancestry_col
    val ancestry_groups
    val region                   // "chr:start-end" format

    output:
    tuple val(meta), path("${prefix}.*.ld.rds"), emit: ld_matrix
    tuple val(meta), path("${prefix}.*.variants.tsv"), emit: variant_info
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    def anc_file_arg = ancestry_file ? "--ancestry_file ${ancestry_file}" : ''
    def anc_col_arg = ancestry_col ? "--ancestry_col ${ancestry_col}" : ''
    def groups_arg = ancestry_groups ? "--groups ${ancestry_groups}" : ''
    """
    calculate_cohort_ld.R \\
        --gds ${gds} \\
        ${anc_file_arg} \\
        ${anc_col_arg} \\
        ${groups_arg} \\
        --region ${region} \\
        --format coloc \\
        --output_prefix ${prefix} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cohort_ld: "1.0.0"
        R: \$(Rscript -e "cat(R.version.string)")
    END_VERSIONS
    """
}

process COHORT_LD_FROM_BED {
    tag "$meta.id"
    label 'process_high'

    // Calculate LD for multiple regions from BED file
    conda "bioconda::bioconductor-seqarray conda-forge::r-data.table conda-forge::r-matrix conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/cohort-ld:latest' :
        'your-registry/cohort-ld:latest' }"

    input:
    tuple val(meta), path(gds)
    path ancestry_file
    val ancestry_col
    val ancestry_groups
    path regions_bed             // BED file with regions

    output:
    tuple val(meta), path("${prefix}.*.ld.rds"), emit: ld_matrices
    tuple val(meta), path("${prefix}.*.variants.tsv"), emit: variant_info
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def anc_file_arg = ancestry_file ? "--ancestry_file ${ancestry_file}" : ''
    def anc_col_arg = ancestry_col ? "--ancestry_col ${ancestry_col}" : ''
    def groups_arg = ancestry_groups ? "--groups ${ancestry_groups}" : ''
    """
    calculate_cohort_ld.R \\
        --gds ${gds} \\
        ${anc_file_arg} \\
        ${anc_col_arg} \\
        ${groups_arg} \\
        --regions_file ${regions_bed} \\
        --format coloc \\
        --threads ${task.cpus} \\
        --output_prefix ${prefix} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cohort_ld: "1.0.0"
        R: \$(Rscript -e "cat(R.version.string)")
    END_VERSIONS
    """
}
