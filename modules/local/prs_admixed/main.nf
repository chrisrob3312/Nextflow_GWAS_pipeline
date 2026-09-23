// Multi-Method PRS for Admixed Populations (thin wrappers around bin/calculate_prs_admixed.R)
// Primary: PRS-CSx | Adjuncts: GAUDI, DiscoDivas, SDPR_admix, MUSSEL, PROSPER
// Plus: local-ancestry PARTIAL scores (DiscoDivas/GAUDI idea: the part of each
// person's score carried on EUR / AFR / AMR haplotypes) and validation that is
// reported overall AND within each ancestry stratum.
// Aligned with the SLURM path (slurm/submit_prs_array.sh), same script.

process PRS_METHOD {
    tag "${meta.trait} - ${method}"
    label 'process_high'

    conda "conda-forge::r-data.table conda-forge::r-optparse conda-forge::python conda-forge::numpy conda-forge::scipy bioconda::plink2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/prs-admixed:latest' :
        'your-registry/prs-admixed:latest' }"

    input:
    tuple val(meta), val(ancestries), path(sumstats_files, stageAs: 'ss_in/*'), path(bed), path(bim), path(fam), path(local_ancestry), path(phenotype)
    each method
    path ld_ref

    output:
    tuple val(meta), val(method), path("${prefix}.${method}*"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.trait}"
    def anc_list = ancestries instanceof List ? ancestries : [ancestries]
    def ss_list  = sumstats_files instanceof List ? sumstats_files : [sumstats_files]
    // stage as <ancestry>.sumstats.gz (what PRS-CSx / SDPR expect)
    def stage = [ss_list, anc_list].transpose().collect { f, a -> "ln -sf \$(readlink -f ${f}) sumstats/${a}.sumstats.gz" }.join('\n    ')
    def first_ss = "sumstats/${anc_list[0]}.sumstats.gz"
    def la_arg   = local_ancestry ? "--local_ancestry ${local_ancestry}" : ''
    def ld_arg   = ld_ref ? "--ld_ref ${ld_ref}" : ''
    def pheno_arg = phenotype ? "--phenotype ${phenotype} --trait ${meta.trait}" : ''
    """
    mkdir -p sumstats
    ${stage}

    calculate_prs_admixed.R \\
        --method ${method} \\
        --sumstats ${first_ss} \\
        --sumstats_dir sumstats \\
        --sumstats_pattern '{ancestry}.sumstats.gz' \\
        ${ld_arg} \\
        --geno ${bed.baseName} \\
        ${la_arg} \\
        --ancestries ${anc_list.join(',')} \\
        ${pheno_arg} \\
        --gwas_type ${meta.gwas_type ?: 'joint'} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        calculate_prs_admixed: "2.0.0"
        method: "${method}"
        plink2: \$(plink2 --version 2>&1 | head -1 | awk '{print \$2}' || echo unknown)
    END_VERSIONS
    """
}

process PRS_LA_PARTIAL {
    tag "${meta.trait}"
    label 'process_medium'

    // Local-ancestry partial PRS: PRS_anc,i = sum_j w_anc,j * Dose_anc,ij
    // Weights = ancestry-specific effects (Tractor-GENESIS BETA_<anc>, or
    // PRS-CSx per-ancestry posterior weights); dosages = Tractor ancdose files.
    conda "conda-forge::r-data.table conda-forge::r-optparse"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/prs-admixed:latest' :
        'your-registry/prs-admixed:latest' }"

    input:
    tuple val(meta), val(ancestries), path(tractor_sumstats), path(prscsx_weights, stageAs: 'prscsx/*'), path(tractor_dosages, stageAs: 'tractor/*'), path(phenotype)

    output:
    tuple val(meta), path("${prefix}.la_partial.scores.tsv"),  emit: scores
    tuple val(meta), path("${prefix}.la_partial.summary.tsv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.trait}"
    def anc_list = ancestries instanceof List ? ancestries : [ancestries]
    def dos = tractor_dosages instanceof List ? tractor_dosages : [tractor_dosages]
    def tractor_prefix = "tractor/" + dos[0].name.replaceAll(/\.(ancdose|dosage)\..*$/, '')
    // Prefer PRS-CSx shrunk per-ancestry posterior weights; fall back to Tractor-GENESIS betas
    def weights_arg = prscsx_weights ? "prscsx" : "${tractor_sumstats}"
    """
    calculate_prs_admixed.R \\
        --method la_partial \\
        --weights ${weights_arg} \\
        --tractor_prefix ${tractor_prefix} \\
        --ancestries ${anc_list.join(',')} \\
        --phenotype ${phenotype} \\
        --trait ${meta.trait} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        --verbose

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        calculate_prs_admixed: "2.0.0"
    END_VERSIONS
    """
}

process PRS_VALIDATE_STRATIFIED {
    tag "${meta.trait}"
    label 'process_low'

    // Truth basis = the phenotype. AUC / R2 / C-index OVERALL and WITHIN each
    // ancestry stratum, with disparity flags; picks the best method.
    conda "conda-forge::r-data.table conda-forge::r-optparse conda-forge::r-proc conda-forge::r-survival"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/prs-admixed:latest' :
        'your-registry/prs-admixed:latest' }"

    input:
    tuple val(meta), path(score_files, stageAs: 'scores/*'), path(phenotype)
    val ancestry_col

    output:
    tuple val(meta), path("${prefix}.validation.tsv"), emit: validation
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.trait}"
    def model = meta.survival ? 'survival' : (meta.binary ? 'binary' : 'quantitative')
    def surv = meta.survival ? "--time_col ${meta.time_col} --event_col ${meta.event_col}" : ''
    """
    # the script looks for *.sscore in the working directory
    for f in scores/*.sscore scores/*.scores.tsv; do [ -e "\$f" ] && ln -sf "\$f" ./; done

    calculate_prs_admixed.R \\
        --method validate \\
        --phenotype ${phenotype} \\
        --trait ${meta.trait} \\
        --trait_type ${model} \\
        ${surv} \\
        --ancestry_col ${ancestry_col} \\
        --validate --stratify_validation \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        calculate_prs_admixed: "2.0.0"
    END_VERSIONS
    """
}

process PRS_COMBINE_METHODS {
    tag "${meta.trait}"
    label 'process_low'

    conda "conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/prs-admixed:latest' :
        'your-registry/prs-admixed:latest' }"

    input:
    tuple val(meta), path(validation)

    output:
    tuple val(meta), path("${prefix}.prs_comparison.tsv"), emit: comparison
    tuple val(meta), path("${prefix}.best_prs.tsv"),       emit: best_prs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.trait}"
    """
    #!/usr/bin/env Rscript
    library(data.table)
    v <- fread("${validation}")
    overall <- v[stratum == "OVERALL"][order(-value)]
    strat   <- v[stratum != "OVERALL"]

    # Best overall, and per stratum; flag ancestry disparity (> 0.10)
    per_stratum <- strat[, .SD[which.max(value)], by = stratum][, .(stratum, best_method = method, value)]
    disparity <- strat[, .(range = max(value, na.rm = TRUE) - min(value, na.rm = TRUE),
                           worst_stratum = stratum[which.min(value)]), by = method]
    comp <- merge(overall[, .(method, metric, overall = value, n)], disparity, by = "method", all.x = TRUE)
    comp[, disparity_flag := !is.na(range) & range > 0.10]
    setorder(comp, -overall)
    fwrite(comp, "${prefix}.prs_comparison.tsv", sep = "\\t")

    best <- rbind(data.table(stratum = "OVERALL", best_method = overall\$method[1], value = overall\$value[1]), per_stratum)
    fwrite(best, "${prefix}.best_prs.tsv", sep = "\\t")
    cat("Best overall:", overall\$method[1], "\\n"); print(best)
    writeLines(c('"${task.process}":', '    prs_combine: "2.0.0"'), "versions.yml")
    """
}
