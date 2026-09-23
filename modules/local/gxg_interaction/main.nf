// GxG (Epistasis) Interaction Testing on GWAS Hits
// END OF PIPELINE: pairwise SNP x SNP tests overall and within ancestry strata,
// plus ancestry-modification of the interaction (Cochran's Q + 3-way LRT)
// and optional local-ancestry-aware epistasis from Tractor dosages.

process GXG_SELECT_HITS {
    tag "${meta.trait}"
    label 'process_low'

    conda "conda-forge::r-data.table conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-r:latest' :
        'your-registry/gwas-r:latest' }"

    input:
    tuple val(meta), path(sumstats_files)
    path known_loci        // default known-risk-loci TSV, or [] to disable
    path custom_variants   // user GRCh38 list (rsID / chr:pos / chr:pos:ref:alt), or []
    val p_threshold
    val max_hits

    output:
    tuple val(meta), path("${prefix}.gxg_hits.txt"), emit: snp_list
    tuple val(meta), path("${prefix}.gxg_hits.tsv"), emit: hits_table
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.trait}"
    def known_arg = known_loci ? "--known_loci ${known_loci}" : ''
    """
    #!/usr/bin/env Rscript
    library(data.table)

    files <- strsplit("${sumstats_files}", " ")[[1]]
    known <- if (nzchar("${known_loci}")) fread("${known_loci}") else data.table()
    p_thr <- ${p_threshold}
    max_hits <- ${max_hits}

    find_col <- function(dt, c) { h <- intersect(c, names(dt)); if (length(h)) h[1] else NA }
    hits <- rbindlist(lapply(files, function(f) {
        ss <- fread(f)
        idc <- find_col(ss, c("ID","SNP","variant_id","rsid"))
        pc  <- find_col(ss, c("P_JOINT","P","PVAL","pvalue","P_META"))
        lp  <- find_col(ss, c("LOG10P"))
        if (is.na(idc)) return(NULL)
        p <- if (!is.na(pc)) ss[[pc]] else if (!is.na(lp)) 10^(-ss[[lp]]) else NA
        data.table(ID = as.character(ss[[idc]]), P = p, SOURCE = basename(f))[!is.na(P) & P < p_thr]
    }), fill = TRUE)

    if (nrow(known) > 0) {
        kc <- find_col(known, c("snp","SNP","ID","rsid"))
        hits <- rbind(hits, data.table(ID = as.character(known[[kc]]), P = 0,
                      SOURCE = paste0("known:", known\$gene)), fill = TRUE)
    }
    # Custom GRCh38 list (rsID / chr:pos / chr:pos:ref:alt) - resolved to genotype IDs in GXG_TEST
    if (nzchar("${custom_variants}") && file.exists("${custom_variants}")) {
        cv <- trimws(readLines("${custom_variants}")); cv <- cv[nzchar(cv) & !grepl("^#", cv)]
        cv <- sapply(strsplit(cv, "[ \\t]+"), `[`, 1)
        hits <- rbind(hits, data.table(ID = cv, P = -1, SOURCE = "custom"), fill = TRUE)
        cat("Custom variants:", length(cv), "\\n")
    }
    hits <- hits[order(P)][, .SD[1], by = ID]
    if (nrow(hits) > max_hits) hits <- hits[1:max_hits]    # custom (P=-1) and known (P=0) rank first

    writeLines(hits\$ID, "${prefix}.gxg_hits.txt")
    fwrite(hits, "${prefix}.gxg_hits.tsv", sep = "\\t")
    cat("Selected", nrow(hits), "hit SNPs for GxG testing\\n")

    writeLines(c('"${task.process}":', paste0('    r-data.table: "', packageVersion("data.table"), '"')), "versions.yml")
    """
}

process GXG_TEST {
    tag "${meta.trait} - ${stratum}"
    label 'process_high'

    // One task per stratum (POOLED, EUR, AAC, LAT1, LAT2, EAS, SAS, OTHER)
    conda "conda-forge::r-data.table conda-forge::r-optparse conda-forge::r-survival bioconda::plink2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-r:latest' :
        'your-registry/gwas-r:latest' }"

    input:
    tuple val(meta), path(snp_list), path(bed), path(bim), path(fam), path(phenotype)
    each stratum
    path known_loci
    path tractor_files, stageAs: 'tractor/*'
    val covariates
    val ancestry_col
    val min_stratum_n

    output:
    tuple val(meta), val(stratum), path("${prefix}.gxg.*.tsv"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.trait}.${stratum}"
    def model = meta.survival ? 'survival' : (meta.binary ? 'binary' : 'quantitative')
    def surv_args = meta.survival ? "--time_col ${meta.time_col} --event_col ${meta.event_col}" : ''
    def cov_arg = covariates ? "--covariates ${covariates}" : ''
    def known_arg = known_loci ? "--known_loci ${known_loci}" : ''
    def tractor_arg = tractor_files ? "--tractor_prefix tractor/${meta.tractor_prefix}" : ''
    """
    run_gxg_interaction.R \\
        --snp_list ${snp_list} \\
        ${known_arg} \\
        --geno ${bed.baseName} \\
        --phenotype ${phenotype} \\
        --trait ${meta.trait} \\
        --model ${model} \\
        ${surv_args} \\
        ${cov_arg} \\
        --ancestry_col ${ancestry_col} \\
        --stratum ${stratum} \\
        --min_stratum_n ${min_stratum_n} \\
        --prune_mode ${params.gxg_prune_mode ?: 'conditional'} \\
        --cond_p_threshold ${params.gxg_cond_p_threshold ?: 1e-4} \\
        --min_distance_kb ${params.gxg_min_distance_kb ?: 1000} \\
        --max_pair_r2 ${params.gxg_max_pair_r2 ?: 0.2} \\
        ${tractor_arg} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-survival: \$(Rscript -e "cat(as.character(packageVersion('survival')))")
        plink2: \$(plink2 --version 2>&1 | head -1 | awk '{print \$2}' || echo unknown)
    END_VERSIONS
    """
}

process GXG_ANCESTRY_HETEROGENEITY {
    tag "${meta.trait}"
    label 'process_medium'

    // Runs the full script once (all strata) to get Cochran's Q across strata
    // and the pooled 3-way g1:g2:ancestry LRT. Cheap relative to GXG_TEST
    // because only the pair-level heterogeneity models are added.
    conda "conda-forge::r-data.table conda-forge::r-optparse conda-forge::r-survival bioconda::plink2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-r:latest' :
        'your-registry/gwas-r:latest' }"

    input:
    tuple val(meta), path(snp_list), path(bed), path(bim), path(fam), path(phenotype)
    path known_loci
    val covariates
    val ancestry_col
    val min_stratum_n

    output:
    tuple val(meta), path("${prefix}.gxg.ancestry_heterogeneity.tsv"), emit: heterogeneity, optional: true
    tuple val(meta), path("${prefix}.gxg.all.tsv"), emit: all_strata
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.trait}"
    def model = meta.survival ? 'survival' : (meta.binary ? 'binary' : 'quantitative')
    def surv_args = meta.survival ? "--time_col ${meta.time_col} --event_col ${meta.event_col}" : ''
    def cov_arg = covariates ? "--covariates ${covariates}" : ''
    def known_arg = known_loci ? "--known_loci ${known_loci}" : ''
    """
    run_gxg_interaction.R \\
        --snp_list ${snp_list} \\
        ${known_arg} \\
        --geno ${bed.baseName} \\
        --phenotype ${phenotype} \\
        --trait ${meta.trait} \\
        --model ${model} \\
        ${surv_args} \\
        ${cov_arg} \\
        --ancestry_col ${ancestry_col} \\
        --min_stratum_n ${min_stratum_n} \\
        --prune_mode ${params.gxg_prune_mode ?: 'conditional'} \\
        --cond_p_threshold ${params.gxg_cond_p_threshold ?: 1e-4} \\
        --min_distance_kb ${params.gxg_min_distance_kb ?: 1000} \\
        --max_pair_r2 ${params.gxg_max_pair_r2 ?: 0.2} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        run_gxg_interaction: "1.1.0"
    END_VERSIONS
    """
}

process GXG_COMBINE {
    tag "${meta.trait}"
    label 'process_low'

    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gwas-r:latest' :
        'your-registry/gwas-r:latest' }"

    input:
    tuple val(meta), path(stratum_results, stageAs: 'strata/*')
    tuple val(meta2), path(heterogeneity)

    output:
    tuple val(meta), path("${prefix}.gxg_summary.tsv"), emit: summary
    tuple val(meta), path("${prefix}.gxg_top_interactions.tsv"), emit: top
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.trait}"
    """
    #!/usr/bin/env Rscript
    library(data.table)

    files <- list.files("strata", pattern = "gxg\\\\.(POOLED|EUR|AAC|LAT1|LAT2|EAS|SAS|OTHER)\\\\.tsv\$",
                        full.names = TRUE, recursive = TRUE)
    res <- rbindlist(lapply(files, fread), fill = TRUE)
    het <- if (file.exists("${heterogeneity}")) fread("${heterogeneity}") else data.table()

    # Wide table: one row per pair, interaction P per stratum + ancestry heterogeneity
    wide <- dcast(res, SNP1 + SNP2 + CHR1 + POS1 + CHR2 + POS2 ~ analysis,
                  value.var = c("BETA_INT", "P_LRT", "N"))
    if (nrow(het) > 0) wide <- merge(wide, het[, .(SNP1, SNP2, strongest_stratum, strongest_p,
                                                    P_HET_INT, I2_INT, P_3WAY_LRT, P_3WAY_FDR)],
                                     by = c("SNP1", "SNP2"), all.x = TRUE)
    if ("P_LRT_POOLED" %in% names(wide)) wide <- wide[order(P_LRT_POOLED)]
    fwrite(wide, "${prefix}.gxg_summary.tsv", sep = "\\t")

    # Top interactions: FDR < 0.05 in any stratum OR ancestry-modified (3-way FDR < 0.10)
    top <- res[P_FDR < 0.05]
    if (nrow(het) > 0) {
        mod <- het[P_3WAY_FDR < 0.10, .(SNP1, SNP2)]
        top <- rbind(top, res[mod, on = c("SNP1", "SNP2")], fill = TRUE)
    }
    top <- unique(top)[order(P_LRT)]
    fwrite(top, "${prefix}.gxg_top_interactions.tsv", sep = "\\t")

    cat("GxG summary:", nrow(wide), "pairs;", nrow(top), "top interaction rows\\n")
    writeLines(c('"${task.process}":', '    gxg_combine: "1.0.0"'), "versions.yml")
    """
}
