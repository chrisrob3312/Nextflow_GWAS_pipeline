// Heritability Estimation Module
// Multiple methods appropriate for admixed populations
// cov-LDSC (preferred for admixed), GREML, BOLT-REML

process HERITABILITY_LDSC {
    tag "${meta.id} - ${meta.ancestry}"
    label 'process_low'

    conda "bioconda::ldsc conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/ldsc:latest' :
        'your-registry/ldsc:latest' }"

    input:
    tuple val(meta), path(sumstats)
    path ld_scores                   // LD score files prefix

    output:
    tuple val(meta), path("${prefix}.hsq"), emit: h2
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    """
    estimate_heritability.R \\
        --method ldsc \\
        --sumstats ${sumstats} \\
        --ld_scores ${ld_scores} \\
        --output_prefix ${prefix} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ldsc: \$(ldsc.py --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}

process HERITABILITY_COV_LDSC {
    tag "${meta.id}"
    label 'process_medium'

    // cov-LDSC: Covariate-stratified LDSC for admixed populations
    conda "bioconda::ldsc conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/ldsc:latest' :
        'your-registry/ldsc:latest' }"

    input:
    tuple val(meta), path(sumstats)
    path ld_scores
    path ancestry_file
    val ancestries                   // Comma-separated list

    output:
    tuple val(meta), path("${prefix}.heritability.rds"), emit: h2
    tuple val(meta), path("${prefix}*.log"), emit: logs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    """
    estimate_heritability.R \\
        --method cov_ldsc \\
        --sumstats ${sumstats} \\
        --ld_scores ${ld_scores} \\
        --ancestry_file ${ancestry_file} \\
        --ancestries ${ancestries} \\
        --output_prefix ${prefix} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cov_ldsc: "1.0.0"
    END_VERSIONS
    """
}

process HERITABILITY_GREML {
    tag "${meta.id} - ${meta.ancestry}"
    label 'process_high'

    // GCTA-GREML: Gold standard, individual-level data
    conda "bioconda::gcta"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gcta:latest' :
        'your-registry/gcta:latest' }"

    input:
    tuple val(meta), path(grm_files)  // .grm.bin, .grm.N.bin, .grm.id
    tuple val(meta), path(phenotype)
    val trait
    val covariates

    output:
    tuple val(meta), path("${prefix}.hsq"), emit: h2
    tuple val(meta), path("${prefix}.heritability.rds"), emit: results
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    def grm_prefix = grm_files[0].baseName.replaceAll(/\\.grm.*/, '')
    def covar_arg = covariates ? "--covariates ${covariates}" : ''
    """
    estimate_heritability.R \\
        --method greml \\
        --grm ${grm_prefix} \\
        --phenotype ${phenotype} \\
        --trait ${trait} \\
        ${covar_arg} \\
        --output_prefix ${prefix} \\
        --threads ${task.cpus} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gcta: \$(gcta64 --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}

process HERITABILITY_ANCESTRY_STRATIFIED {
    tag "${meta.id}"
    label 'process_high'

    // Run GREML separately per ancestry, then combine
    conda "bioconda::gcta conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/gcta:latest' :
        'your-registry/gcta:latest' }"

    input:
    tuple val(meta), path(plink_files)
    tuple val(meta), path(phenotype)
    path ancestry_file
    val trait
    val ancestries
    val covariates

    output:
    tuple val(meta), path("${prefix}.combined_h2.tsv"), emit: combined
    tuple val(meta), path("${prefix}.*.hsq"), emit: per_ancestry
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.trait}"
    def covar_arg = covariates ? "--covariates ${covariates}" : ''
    """
    #!/bin/bash
    set -euo pipefail

    # Run GREML for each ancestry
    for anc in \$(echo "${ancestries}" | tr ',' ' '); do
        echo "Processing ancestry: \$anc"

        # Extract samples for this ancestry
        awk -v anc="\$anc" '\$2 == anc {print \$1, \$1}' ${ancestry_file} > \${anc}.samples.txt

        if [[ \$(wc -l < \${anc}.samples.txt) -lt 100 ]]; then
            echo "  Skipping \$anc: insufficient samples"
            continue
        fi

        # Subset PLINK files
        plink2 --bfile ${plink_files[0].baseName} \\
            --keep \${anc}.samples.txt \\
            --make-bed \\
            --out ${prefix}.\${anc}

        # Compute GRM
        gcta64 --bfile ${prefix}.\${anc} \\
            --make-grm \\
            --out ${prefix}.\${anc} \\
            --thread-num ${task.cpus}

        # Run GREML
        estimate_heritability.R \\
            --method greml \\
            --grm ${prefix}.\${anc} \\
            --phenotype ${phenotype} \\
            --trait ${trait} \\
            ${covar_arg} \\
            --output_prefix ${prefix}.\${anc} \\
            --threads ${task.cpus}
    done

    # Combine results
    echo -e "ancestry\\th2\\tse\\tn" > ${prefix}.combined_h2.tsv
    for f in ${prefix}.*.hsq; do
        anc=\$(basename \$f | sed 's/${prefix}.//; s/.hsq//')
        h2=\$(awk '/V\\(G\\)\\/Vp/ {print \$2}' \$f)
        se=\$(awk '/V\\(G\\)\\/Vp/ {print \$3}' \$f)
        n=\$(awk '/^n/ {print \$2}' \$f)
        echo -e "\${anc}\\t\${h2}\\t\${se}\\t\${n}" >> ${prefix}.combined_h2.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gcta: \$(gcta64 --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """
}

process GENETIC_CORRELATION {
    tag "${meta.id1} vs ${meta.id2}"
    label 'process_medium'

    // Cross-ancestry or cross-trait genetic correlation
    conda "bioconda::ldsc conda-forge::r-data.table"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://your-registry/ldsc:latest' :
        'your-registry/ldsc:latest' }"

    input:
    tuple val(meta), path(sumstats1), path(sumstats2)
    path ld_scores
    val method                       // ldsc, popcorn, s_ldxr

    output:
    tuple val(meta), path("${prefix}.rg.tsv"), emit: rg
    tuple val(meta), path("${prefix}.rg.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id1}_vs_${meta.id2}"
    """
    estimate_heritability.R \\
        --method ldsc \\
        --sumstats ${sumstats1} \\
        --sumstats2 ${sumstats2} \\
        --ld_scores ${ld_scores} \\
        --estimate_rg \\
        --rg_method ${method} \\
        --output_prefix ${prefix} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rg_method: "${method}"
    END_VERSIONS
    """
}
