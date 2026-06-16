// Multi-Method Colocalization Workflow
// TIER 1: coloc.susie (per QTL type) → TIER 2: HyPrColoc → TIER 3: OPERA

process COLOC_TIER1 {
    tag "${meta.id} - ${qtl_type}"
    label 'process_medium'

    // TIER 1: coloc.susie per QTL type individually
    conda "bioconda::r-coloc conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/coloc:latest' :
        'your-registry/coloc:latest' }"

    input:
    tuple val(meta), path(gwas)
    tuple val(meta), path(qtl)
    tuple val(meta), path(ld_matrix)
    val qtl_type
    val gwas_n
    val qtl_n

    output:
    tuple val(meta), path("${prefix}.tier1.tsv"), emit: results
    tuple val(meta), path("${prefix}.tier1.candidates.tsv"), emit: candidates, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}.${qtl_type}"
    def ld_arg = ld_matrix ? "--ld ${ld_matrix}" : ''
    """
    run_colocalization.R \\
        --gwas ${gwas} \\
        --qtl ${qtl} \\
        ${ld_arg} \\
        --gwas_n ${gwas_n} \\
        --qtl_n ${qtl_n} \\
        --qtl_types ${qtl_type} \\
        --tier2_threshold 0.5 \\
        --output_prefix ${prefix} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coloc: \$(Rscript -e "cat(as.character(packageVersion('coloc')))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}

process COLOC_TIER2_HYPRCOLOC {
    tag "${meta.id}"
    label 'process_medium'

    // TIER 2: HyPrColoc for multi-trait colocalization
    conda "r-hyprcoloc conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/hyprcoloc:latest' :
        'your-registry/hyprcoloc:latest' }"

    input:
    tuple val(meta), path(tier1_candidates)
    tuple val(meta), path(gwas)
    tuple val(meta), path(qtl_files)
    tuple val(meta), path(ld_matrix)

    output:
    tuple val(meta), path("${prefix}.tier2.hyprcoloc.tsv"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(hyprcoloc)

    # Load candidates from Tier 1
    candidates <- fread("${tier1_candidates}")

    cat("TIER 2: HyPrColoc\\n")
    cat("Candidates:", nrow(candidates), "\\n")

    # Run HyPrColoc on candidates
    # ... implementation

    cat("HyPrColoc implementation - placeholder\\n")

    # Write results
    fwrite(candidates, "${prefix}.tier2.hyprcoloc.tsv", sep = "\\t")

    writeLines(c(
        '"${task.process}":',
        paste0('    hyprcoloc: "', packageVersion("hyprcoloc"), '"')
    ), "versions.yml")
    """
}

process COLOC_TIER3_OPERA {
    tag "${meta.id}"
    label 'process_medium'

    // TIER 3: OPERA (SMR-based, NOT causal mediation)
    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/opera:latest' :
        'your-registry/opera:latest' }"

    input:
    tuple val(meta), path(tier2_results)
    tuple val(meta), path(gwas)
    tuple val(meta), path(qtl_files)

    output:
    tuple val(meta), path("${prefix}.tier3.opera.tsv"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    #!/usr/bin/env Rscript

    library(data.table)

    cat("TIER 3: OPERA\\n")
    cat("Note: OPERA uses SMR coding but is NOT causal mediation\\n")
    cat("Tests for shared causal variants between GWAS and QTL\\n")

    # OPERA implementation
    # ... placeholder

    tier2 <- fread("${tier2_results}")
    fwrite(tier2, "${prefix}.tier3.opera.tsv", sep = "\\t")

    writeLines(c(
        '"${task.process}":',
        '    opera: "1.0.0"'
    ), "versions.yml")
    """
}

process COLOC_COMBINE_TIERS {
    tag "${meta.id}"
    label 'process_low'

    // Combine results across all tiers
    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/coloc:latest' :
        'your-registry/coloc:latest' }"

    input:
    tuple val(meta), path(tier1_results)
    tuple val(meta), path(tier2_results), optional: true
    tuple val(meta), path(tier3_results), optional: true

    output:
    tuple val(meta), path("${prefix}.coloc_summary.tsv"), emit: summary
    tuple val(meta), path("${prefix}.coloc_final.tsv"), emit: final
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    #!/usr/bin/env Rscript

    library(data.table)

    # Combine all tier results
    tier1 <- rbindlist(lapply(list.files(".", pattern = "tier1.tsv"), fread), fill = TRUE)

    # Add tier annotations
    tier1[, tier := "tier1"]
    tier1[, evidence := ifelse(PP4 >= 0.8, "strong", ifelse(PP4 >= 0.5, "moderate", "weak"))]

    # Write summary
    fwrite(tier1, "${prefix}.coloc_summary.tsv", sep = "\\t")

    # Final results (high-confidence only)
    final <- tier1[PP4 >= 0.8]
    fwrite(final, "${prefix}.coloc_final.tsv", sep = "\\t")

    cat("Colocalization complete\\n")
    cat("  Total tested:", nrow(tier1), "\\n")
    cat("  High-confidence (PP4 >= 0.8):", nrow(final), "\\n")

    writeLines(c(
        '"${task.process}":',
        '    coloc_combine: "1.0.0"'
    ), "versions.yml")
    """
}
