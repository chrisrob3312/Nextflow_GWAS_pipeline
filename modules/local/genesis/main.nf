process GENESIS_NULL_MODEL {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_high'

    // GENESIS is an R/Bioconductor package optimized for admixed populations
    conda "bioconda::bioconductor-genesis=2.32.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-genesis:2.32.0' :
        'quay.io/biocontainers/bioconductor-genesis:2.32.0' }"

    input:
    tuple val(meta), path(gds), path(phenotype)
    val covariate_cols
    path kinship_matrix  // GRM for relatedness adjustment

    output:
    tuple val(meta), path("${prefix}.null_model.rds"), emit: null_model
    tuple val(meta), path("${prefix}.log"), emit: log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    def covar_arg = covariate_cols ? "covars = c('${covariate_cols.split(',').join(\"', '\")}')" : "covars = NULL"
    def family = meta.binary ? "binomial" : "gaussian"
    """
    #!/usr/bin/env Rscript

    library(GENESIS)
    library(GWASTools)
    library(SeqArray)
    library(SeqVarTools)

    # Load phenotype data (TSV; sample ID column sample_id / IID -> scanID for GENESIS)
    pheno <- read.table("${phenotype}", header=TRUE, sep="\\t", stringsAsFactors=FALSE, check.names=FALSE)
    idc <- intersect(c("sample_id", "IID", "ID", "sample"), names(pheno))[1]
    if (is.na(idc)) stop("Phenotype file needs a sample_id or IID column")
    pheno\$scanID <- as.character(pheno[[idc]])
    pheno <- pheno[, c("scanID", setdiff(names(pheno), "scanID"))]

    # Load kinship / GRM if provided (PC-Relate GRM from GENESIS_PCAIR_PCRELATE, or user file)
    kinship <- NULL
    if (file.exists("${kinship_matrix}")) {
        kinship <- readRDS("${kinship_matrix}")
        common <- intersect(rownames(kinship), pheno\$scanID)
        pheno <- pheno[pheno\$scanID %in% common, ]
        kinship <- kinship[pheno\$scanID, pheno\$scanID]
    }

    # Create AnnotatedDataFrame
    scanAnnot <- ScanAnnotationDataFrame(pheno)

    # Fit null model
    # GENESIS handles:
    # - Population structure via PCs
    # - Relatedness via kinship/GRM
    # - Admixed populations with heterogeneous ancestry
    nullmod <- fitNullModel(
        scanAnnot,
        outcome = "${meta.trait}",
        ${covar_arg},
        cov.mat = kinship,
        family = "${family}",
        verbose = TRUE
    )

    # Save null model
    saveRDS(nullmod, "${prefix}.null_model.rds")

    # Log model summary
    sink("${prefix}.log")
    print(summary(nullmod))
    sink()

    # Version info
    writeLines(c(
        '"${task.process}":',
        paste0('    genesis: "', packageVersion("GENESIS"), '"'),
        paste0('    R: "', R.version.string, '"')
    ), "versions.yml")
    """
}

process GENESIS_ASSOC {
    tag "$meta.id - $meta.ancestry - $meta.trait"
    label 'process_high'

    conda "bioconda::bioconductor-genesis=2.32.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-genesis:2.32.0' :
        'quay.io/biocontainers/bioconductor-genesis:2.32.0' }"

    input:
    tuple val(meta), path(gds), path(null_model)
    val test_type  // 'Score', 'Wald', 'BinomiRare', 'CMP'

    output:
    tuple val(meta), path("${prefix}.genesis.assoc.tsv.gz"), emit: summary_stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.${meta.ancestry}.${meta.trait}"
    """
    #!/usr/bin/env Rscript

    library(GENESIS)
    library(SeqArray)
    library(SeqVarTools)

    # Open GDS file
    gds <- seqOpen("${gds}")

    # Load null model
    nullmod <- readRDS("${null_model}")

    # Create SeqVarData object
    seqData <- SeqVarData(gds)

    # Create iterator for genome-wide scan
    iterator <- SeqVarBlockIterator(seqData, verbose=TRUE)

    # Run association test
    # GENESIS supports multiple test types optimized for different scenarios:
    # - Score: fast, good for common variants
    # - Wald: provides effect sizes
    # - BinomiRare: for rare variants in case-control
    # - CMP: robust to population stratification in admixed
    assoc <- assocTestSingle(
        iterator,
        nullmod,
        test = "${test_type}",
        verbose = TRUE
    )

    # Format results
    results <- as.data.frame(assoc)

    # Write results
    write.table(
        results,
        gzfile("${prefix}.genesis.assoc.tsv.gz"),
        sep = "\\t",
        row.names = FALSE,
        quote = FALSE
    )

    seqClose(gds)

    # Version info
    writeLines(c(
        '"${task.process}":',
        paste0('    genesis: "', packageVersion("GENESIS"), '"')
    ), "versions.yml")
    """
}
