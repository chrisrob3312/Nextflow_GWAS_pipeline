// Multi-Method Colocalization Workflow
// PARALLEL ANALYSIS (each method handles its own MTC):
// 1. coloc.susie (per QTL type)
// 2. HyPrColoc (multi-trait, ALL QTL types)
// 3. OPERA (multi-QTL SMR, ALL QTLs with built-in MTC)

process COLOC_SUSIE {
    tag "${meta.id} - ${qtl_type}"
    label 'process_medium'

    // coloc.susie: Fine-mapping based colocalization per QTL type
    conda "bioconda::r-coloc bioconda::r-susier conda-forge::r-data.table"
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
    tuple val(meta), path("${prefix}.coloc.tsv"), emit: results
    tuple val(meta), path("${prefix}.coloc.significant.tsv"), emit: significant, optional: true
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
        --methods coloc_susie \\
        --output_prefix ${prefix} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coloc: \$(Rscript -e "cat(as.character(packageVersion('coloc')))" 2>/dev/null || echo "unknown")
        susieR: \$(Rscript -e "cat(as.character(packageVersion('susieR')))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}

process HYPRCOLOC {
    tag "${meta.id}"
    label 'process_medium'

    // HyPrColoc: Multi-trait colocalization across ALL QTL types
    conda "r-hyprcoloc conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/hyprcoloc:latest' :
        'your-registry/hyprcoloc:latest' }"

    input:
    tuple val(meta), path(gwas)
    tuple val(meta), path(qtl_dir)
    tuple val(meta), path(ld_matrix)
    val qtl_types  // ALL types: eqtl,sqtl,pqtl,mqtl,caqtl,hqtl

    output:
    tuple val(meta), path("${prefix}.hyprcoloc.tsv"), emit: results
    tuple val(meta), path("${prefix}.hyprcoloc.significant.tsv"), emit: significant, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    run_colocalization.R \\
        --gwas ${gwas} \\
        --qtl_dir ${qtl_dir} \\
        --ld ${ld_matrix} \\
        --qtl_types ${qtl_types} \\
        --methods hyprcoloc \\
        --include_all_qtl_types true \\
        --output_prefix ${prefix} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hyprcoloc: \$(Rscript -e "cat(as.character(packageVersion('hyprcoloc')))" 2>/dev/null || echo "unknown")
    END_VERSIONS
    """
}

process OPERA {
    tag "${meta.id}"
    label 'process_high'

    // OPERA: Multi-QTL SMR analysis - runs on ALL QTLs with built-in MTC
    // NOT causal mediation - tests for shared causal variants
    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/opera:latest' :
        'your-registry/opera:latest' }"

    input:
    tuple val(meta), path(gwas)
    tuple val(meta), path(qtl_dir)
    tuple val(meta), path(geno_files)
    val qtl_types  // ALL types for simultaneous testing

    output:
    tuple val(meta), path("${prefix}.opera.tsv"), emit: results
    tuple val(meta), path("${prefix}.opera.significant.tsv"), emit: significant, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    run_colocalization.R \\
        --gwas ${gwas} \\
        --qtl_dir ${qtl_dir} \\
        --geno ${geno_files[0].baseName} \\
        --qtl_types ${qtl_types} \\
        --methods opera \\
        --run_opera true \\
        --include_all_qtl_types true \\
        --opera_fdr 0.05 \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        opera: "1.0.0"
    END_VERSIONS
    """
}

process COLOC_COMBINE {
    tag "${meta.id}"
    label 'process_low'

    // Combine results from parallel methods
    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/coloc:latest' :
        'your-registry/coloc:latest' }"

    input:
    tuple val(meta), path(coloc_results)
    tuple val(meta), path(hyprcoloc_results), stageAs: 'hyprcoloc/*'
    tuple val(meta), path(opera_results), stageAs: 'opera/*'

    output:
    tuple val(meta), path("${prefix}.coloc_combined.tsv"), emit: combined
    tuple val(meta), path("${prefix}.coloc_final.tsv"), emit: final
    tuple val(meta), path("${prefix}.method_comparison.tsv"), emit: comparison
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    #!/usr/bin/env Rscript

    library(data.table)

    # Load coloc.susie results
    coloc_files <- list.files(".", pattern = "coloc\\\\.tsv\$", full.names = TRUE)
    coloc_all <- rbindlist(lapply(coloc_files, fread), fill = TRUE)
    coloc_all[, method := "coloc_susie"]

    # Load HyPrColoc results
    hypr_files <- list.files("hyprcoloc", pattern = "\\\\.tsv\$", full.names = TRUE)
    if (length(hypr_files) > 0) {
        hypr_all <- rbindlist(lapply(hypr_files, fread), fill = TRUE)
        hypr_all[, method := "hyprcoloc"]
    } else {
        hypr_all <- data.table()
    }

    # Load OPERA results
    opera_files <- list.files("opera", pattern = "\\\\.tsv\$", full.names = TRUE)
    if (length(opera_files) > 0) {
        opera_all <- rbindlist(lapply(opera_files, fread), fill = TRUE)
        opera_all[, method := "opera"]
    } else {
        opera_all <- data.table()
    }

    # Combine all methods
    all_results <- rbindlist(list(coloc_all, hypr_all, opera_all), fill = TRUE)

    # Consensus scoring
    if ("gene" %in% names(all_results)) {
        consensus <- all_results[, .(
            n_methods = .N,
            methods = paste(unique(method), collapse = ";"),
            max_pp4 = max(PP4, na.rm = TRUE),
            min_smr_p = min(smr_p, na.rm = TRUE)
        ), by = .(gene, qtl_type)]

        # High confidence: multiple methods agree
        consensus[, confidence := fifelse(
            n_methods >= 2 & (max_pp4 >= 0.8 | min_smr_p < 0.05),
            "high",
            fifelse(max_pp4 >= 0.5, "medium", "low")
        )]
    }

    # Write combined results
    fwrite(all_results, "${prefix}.coloc_combined.tsv", sep = "\\t")

    # Final high-confidence results
    if (exists("consensus")) {
        final <- consensus[confidence %in% c("high", "medium")]
        fwrite(final, "${prefix}.coloc_final.tsv", sep = "\\t")
    } else {
        fwrite(all_results, "${prefix}.coloc_final.tsv", sep = "\\t")
    }

    # Method comparison
    comparison <- all_results[, .(
        n_genes = length(unique(gene)),
        n_significant = sum(PP4 >= 0.8 | smr_p < 0.05, na.rm = TRUE)
    ), by = .(method, qtl_type)]
    fwrite(comparison, "${prefix}.method_comparison.tsv", sep = "\\t")

    cat("Colocalization summary:\\n")
    cat("  Total gene-QTL pairs:", nrow(all_results), "\\n")
    cat("  Methods used:", paste(unique(all_results\$method), collapse = ", "), "\\n")

    writeLines(c(
        '"${task.process}":',
        '    coloc_combine: "1.0.0"'
    ), "versions.yml")
    """
}

// Workflow to run all colocalization methods in parallel
workflow COLOCALIZATION {
    take:
    gwas           // tuple: meta, gwas_sumstats
    qtl_dir        // tuple: meta, qtl_directory (all QTL types)
    ld_matrix      // tuple: meta, ld_matrix
    geno_files     // tuple: meta, [bed, bim, fam]
    qtl_types_str  // string: comma-separated QTL types

    main:
    // Parse QTL types
    qtl_types = qtl_types_str.tokenize(',')

    // Run coloc.susie per QTL type (parallel across types)
    coloc_results = COLOC_SUSIE(
        gwas,
        qtl_dir.map { meta, dir -> [meta, file("${dir}/*.harmonized.tsv.gz")] },
        ld_matrix,
        qtl_types,
        params.gwas_n,
        params.qtl_n
    )

    // Run HyPrColoc on all QTL types simultaneously
    hyprcoloc_results = HYPRCOLOC(
        gwas,
        qtl_dir,
        ld_matrix,
        qtl_types_str
    )

    // Run OPERA on all QTL types (built-in MTC)
    opera_results = OPERA(
        gwas,
        qtl_dir,
        geno_files,
        qtl_types_str
    )

    // Combine all results
    combined = COLOC_COMBINE(
        coloc_results.results.collect(),
        hyprcoloc_results.results,
        opera_results.results
    )

    emit:
    coloc      = coloc_results.results
    hyprcoloc  = hyprcoloc_results.results
    opera      = opera_results.results
    combined   = combined.combined
    final      = combined.final
    comparison = combined.comparison
    versions   = combined.versions
}
