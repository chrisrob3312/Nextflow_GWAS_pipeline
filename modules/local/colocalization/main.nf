// Multi-Method Colocalization (thin wrappers around bin/run_colocalization.R)
// PARALLEL methods, each with its own multiple-testing handling:
//   coloc.susie per QTL type | HyPrColoc across ALL QTL types | OPERA on ALL QTLs
// One COLOC_RUN task per GWAS (stratum x trait, or meta) x QTL megaset;
// parallelism comes from the many GWAS inputs. COLOC_COMBINE builds the
// per-trait consensus across strata and methods.

process COLOC_RUN {
    tag "${meta.id ?: meta.trait} - ${meta.ancestry ?: meta.analysis_type}"
    label 'process_medium'

    conda "bioconda::r-coloc bioconda::r-susier conda-forge::r-data.table conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/coloc:latest' :
        'your-registry/coloc:latest' }"

    input:
    tuple val(meta), path(gwas), path(qtl_files, stageAs: 'qtl_in/*'), val(qtl_types), path(ld_matrix)
    val gwas_n
    val qtl_n
    val p1
    val p2
    val p12

    output:
    tuple val(meta), path("${prefix}.coloc_susie.tsv"),          emit: coloc_susie
    tuple val(meta), path("${prefix}.hyprcoloc.tsv"),            emit: hyprcoloc, optional: true
    tuple val(meta), path("${prefix}.opera.tsv"),                emit: opera, optional: true
    tuple val(meta), path("${prefix}.opera.significant.tsv"),    emit: opera_significant, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id ?: meta.trait}.${meta.ancestry ?: meta.analysis_type ?: 'meta'}.${meta.trait}"
    def ld_arg = ld_matrix ? "--ld ${ld_matrix}" : ''
    def s_arg  = meta.case_prop ? "--gwas_s ${meta.case_prop}" : ''
    def gtype  = meta.binary ? 'cc' : 'quant'
    // stage pooled QTL files as <type>.harmonized.tsv.gz (what the script expects)
    def stage = [qtl_files instanceof List ? qtl_files : [qtl_files], qtl_types instanceof List ? qtl_types : [qtl_types]]
        .transpose().collect { f, t -> "ln -sf \$(readlink -f ${f}) qtl/${t}.harmonized.tsv.gz" }.join('\n    ')
    """
    mkdir -p qtl
    ${stage}

    run_colocalization.R \\
        --gwas ${gwas} \\
        --gwas_n ${gwas_n} \\
        --gwas_type ${gtype} \\
        ${s_arg} \\
        --qtl_dir qtl \\
        --qtl_n ${qtl_n} \\
        --qtl_types ${(qtl_types instanceof List ? qtl_types : [qtl_types]).join(',')} \\
        --methods coloc_susie,hyprcoloc,opera \\
        --run_opera \\
        ${ld_arg} \\
        --p1 ${p1} --p2 ${p2} --p12 ${p12} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coloc: \$(Rscript -e "cat(as.character(packageVersion('coloc')))" 2>/dev/null || echo unknown)
        susieR: \$(Rscript -e "cat(as.character(packageVersion('susieR')))" 2>/dev/null || echo unknown)
        hyprcoloc: \$(Rscript -e "cat(tryCatch(as.character(packageVersion('hyprcoloc')), error=function(e) 'not installed'))")
    END_VERSIONS
    """
}

process COLOC_COMBINE {
    tag "${meta.trait}"
    label 'process_low'

    // Per-trait consensus across strata (EUR, AAC, LAT1, LAT2, ..., meta) and methods
    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/coloc:latest' :
        'your-registry/coloc:latest' }"

    input:
    tuple val(meta), path(result_files, stageAs: 'results/*')

    output:
    tuple val(meta), path("${prefix}.coloc_combined.tsv"),   emit: combined
    tuple val(meta), path("${prefix}.coloc_consensus.tsv"),  emit: consensus
    tuple val(meta), path("${prefix}.method_comparison.tsv"), emit: comparison
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.trait}"
    """
    #!/usr/bin/env Rscript
    library(data.table)

    load_set <- function(pattern, method) {
        fs <- list.files("results", pattern = pattern, full.names = TRUE)
        if (length(fs) == 0) return(data.table())
        rbindlist(lapply(fs, function(f) {
            d <- fread(f); d[, method := method]
            d[, source_file := basename(f)]
            d[, stratum := sub("^[^.]+\\\\.([^.]+)\\\\..*\$", "\\\\1", basename(f))]
            d
        }), fill = TRUE)
    }
    cs <- load_set("coloc_susie\\\\.tsv\$", "coloc_susie")
    hy <- load_set("hyprcoloc\\\\.tsv\$",  "hyprcoloc")
    op <- load_set("opera\\\\.tsv\$",      "opera")
    all <- rbindlist(list(cs, hy, op), fill = TRUE)
    if (nrow(all) == 0) all <- data.table(gene = character(), method = character())

    if (!"PP4" %in% names(all)) all[, PP4 := NA_real_]
    if (!"smr_p" %in% names(all)) all[, smr_p := NA_real_]
    if (!"posterior_prob" %in% names(all)) all[, posterior_prob := NA_real_]
    if (!"qtl_type" %in% names(all)) all[, qtl_type := NA_character_]

    fwrite(all, "${prefix}.coloc_combined.tsv", sep = "\\t")

    # Consensus per gene x QTL type: how many methods / strata support it
    cons <- all[!is.na(gene), .(
        n_methods = length(unique(method)),
        methods = paste(sort(unique(method)), collapse = ";"),
        n_strata = length(unique(stratum)),
        strata = paste(sort(unique(stratum)), collapse = ";"),
        max_pp4 = suppressWarnings(max(PP4, na.rm = TRUE)),
        max_hyprcoloc_pp = suppressWarnings(max(posterior_prob, na.rm = TRUE)),
        min_smr_p = suppressWarnings(min(smr_p, na.rm = TRUE))
    ), by = .(gene, qtl_type)]
    cons[is.infinite(max_pp4), max_pp4 := NA]; cons[is.infinite(max_hyprcoloc_pp), max_hyprcoloc_pp := NA]; cons[is.infinite(min_smr_p), min_smr_p := NA]
    cons[, confidence := fifelse(n_methods >= 2 & (max_pp4 >= 0.8 | max_hyprcoloc_pp >= 0.8 | min_smr_p < 0.05), "high",
                         fifelse(max_pp4 >= 0.5 | max_hyprcoloc_pp >= 0.5 | min_smr_p < 0.05, "medium", "low"))]
    cons[is.na(confidence), confidence := "low"]
    setorder(cons, -n_methods, -max_pp4, na.last = TRUE)
    fwrite(cons, "${prefix}.coloc_consensus.tsv", sep = "\\t")

    comp <- all[, .(n_genes = length(unique(gene)),
                    n_significant = sum(PP4 >= 0.8 | posterior_prob >= 0.8 | smr_p < 0.05, na.rm = TRUE)),
                by = .(method, stratum, qtl_type)]
    fwrite(comp, "${prefix}.method_comparison.tsv", sep = "\\t")
    cat("Combined", nrow(all), "rows;", nrow(cons), "gene x QTL-type consensus entries;",
        sum(cons\$confidence == "high"), "high confidence\\n")

    writeLines(c('"${task.process}":', '    coloc_combine: "2.0.0"'), "versions.yml")
    """
}
