#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ancestry-Aware GWAS Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A comprehensive Nextflow pipeline for multi-ancestry genome-wide association studies
    with fine-mapping, meta-analysis, PRS calculation, and functional annotation

    Key Features:
    - GRAF-ANC ancestry inference (9 groups: EUR, AFR, EAS, SAS, AMR, MID, AAC, AHI, HET)
    - Parallelized GWAS by ancestry and trait
    - MR-MEGA/MR-MEGA-env meta-analysis across ancestries
    - PolyFun+SuSIE within-ancestry fine-mapping
    - MG-FLASH-FM multi-ancestry fine-mapping (related traits)
    - SuSIE-ME multi-ancestry fine-mapping (single/unrelated traits)
    - Multiple PRS methods (PRS-CSx, PRS-CS, GAUDI, LDpred2)
    - Ancestry-specific heritability estimation
    - Functional annotation (MAGMA, FUMA, LAVA, FLAMES)
    - Colocalization analysis with diverse QTL datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT HELP / VERSION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.help) {
    helpMessage()
    exit 0
}

if (params.version) {
    println "Ancestry-Aware GWAS Pipeline v${workflow.manifest.version}"
    exit 0
}

def helpMessage() {
    log.info """
    =========================================
    Ancestry-Aware GWAS Pipeline v${workflow.manifest.version}
    =========================================

    Usage:
        nextflow run main.nf [options]

    Mandatory Arguments:
        --input                 Path to samplesheet with genotype data
        --phenotype_file        Path to phenotype file
        --phenotype_cols        Phenotype column names (comma-separated)
        --outdir                Output directory

    Ancestry Options:
        --run_ancestry_inference  Run ancestry inference (default: false, can use pre-computed)
        --ancestry_calls_file     Pre-computed ancestry calls file
        --ancestry_method         Ancestry inference method [graf-anc, admixture, pca] (default: graf-anc)

    GWAS Options:
        --gwas_tool             GWAS tool [regenie, saige, bolt-lmm, plink2, genesis] (default: regenie)
        --gwas_model            Genetic model [additive, dominant, recessive] (default: additive)
        --kinship_matrix        Pre-computed kinship matrix (for GENESIS)

    Tractor Options (Local Ancestry-Aware GWAS):
        --run_tractor           Run Tractor for admixed populations (default: false)
        --tractor_aac_pops      Ancestral populations for AAC (default: EUR,AFR)
        --tractor_lat_pops      Ancestral populations for Latino (default: EUR,AFR,NAT)
        --local_ancestry_files  Local ancestry MSP files (RFMix format)

    Survival Analysis:
        --survival_analysis     Run survival GWAS with SPA-Cox (default: false)
        --time_col              Time-to-event column name
        --event_col             Event indicator column name

    Analysis Options:
        --meta_analysis         Run MR-MEGA meta-analysis [true/false] (default: true)
        --fine_mapping          Run fine-mapping [true/false] (default: true)
        --prs_analysis          Run PRS analysis [true/false] (default: true)
        --heritability          Run heritability analysis [true/false] (default: true)
        --colocalization        Run colocalization [true/false] (default: true)
        --functional_annotation Run functional annotation [true/false] (default: true)

    For full parameter list, see nextflow.config or documentation.
    """.stripIndent()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check mandatory parameters
def checkMandatoryParams() {
    def errors = []
    if (!params.input) errors << "Missing mandatory parameter: --input"
    if (!params.phenotype_file) errors << "Missing mandatory parameter: --phenotype_file"
    if (!params.phenotype_cols) errors << "Missing mandatory parameter: --phenotype_cols"

    if (errors) {
        log.error "Parameter validation failed:\n" + errors.join("\n")
        exit 1
    }
}

checkMandatoryParams()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT WORKFLOWS AND MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Subworkflows
include { INPUT_CHECK         } from './subworkflows/input_check'
include { QC_WORKFLOW         } from './subworkflows/qc'
include { ANCESTRY_INFERENCE  } from './subworkflows/ancestry'
include { GWAS_WORKFLOW       } from './subworkflows/gwas'
include { META_ANALYSIS       } from './subworkflows/meta_analysis'
include { FINE_MAPPING        } from './subworkflows/fine_mapping'
include { COLOCALIZATION      } from './subworkflows/colocalization'
include { PRS_WORKFLOW        } from './subworkflows/prs'
include { HERITABILITY        } from './subworkflows/heritability'
include { FUNCTIONAL_ANNOT    } from './subworkflows/functional'
include { VISUALIZATION       } from './subworkflows/visualization'
include { REPORTING           } from './subworkflows/reporting'
include { GXG_INTERACTION     } from './subworkflows/gxg'

// Utility modules
include { SOFTWARE_VERSIONS   } from './modules/local/software_versions'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {
    // Create channels
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Parse input samplesheet
    INPUT_CHECK(
        file(params.input)
    )
    ch_genotypes = INPUT_CHECK.out.genotypes
    ch_phenotypes = Channel.fromPath(params.phenotype_file)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Parse traits and assign each one a Tractor-GENESIS model so every trait
    // goes through the same adapter (MRD -> binary; OS, relapse -> survival)
    def survival_traits = params.survival_traits ? params.survival_traits.split(',') as List : []
    def surv_time_cols  = params.survival_time_cols  ? params.survival_time_cols.split(',')  as List : []
    def surv_event_cols = params.survival_event_cols ? params.survival_event_cols.split(',') as List : []

    ch_traits = Channel.of(params.phenotype_cols.split(','))
        .map { trait ->
            def is_binary   = params.binary_traits?.split(',')?.contains(trait) ?: false
            def is_survival = survival_traits.contains(trait)
            def idx         = survival_traits.indexOf(trait)
            def time_col    = is_survival ? (idx < surv_time_cols.size()  ? surv_time_cols[idx]  : "${trait}_time")   : null
            def event_col   = is_survival ? (idx < surv_event_cols.size() ? surv_event_cols[idx] : "${trait}_status") : null
            def model       = is_survival ? 'survival' : (is_binary ? 'binary' : 'quantitative')
            [trait: trait, binary: is_binary, survival: is_survival, model: model,
             time_col: time_col, event_col: event_col]
        }

    // =========================================================================
    // QUALITY CONTROL
    // =========================================================================
    QC_WORKFLOW(
        ch_genotypes,
        params.maf_threshold,
        params.hwe_threshold,
        params.geno_missing,
        params.mind_missing
    )
    ch_qc_genotypes = QC_WORKFLOW.out.genotypes
    ch_versions = ch_versions.mix(QC_WORKFLOW.out.versions)

    // =========================================================================
    // ANCESTRY INFERENCE
    // =========================================================================
    ANCESTRY_INFERENCE(
        ch_qc_genotypes,
        params.ancestry_method,
        params.graf_reference ? file(params.graf_reference) : [],
        params.local_ancestry,
        params.local_ancestry_method
    )
    ch_ancestry_calls = ANCESTRY_INFERENCE.out.ancestry_calls
    ch_stratified_genotypes = ANCESTRY_INFERENCE.out.stratified_genotypes
    ch_local_ancestry = ANCESTRY_INFERENCE.out.local_ancestry
    ch_versions = ch_versions.mix(ANCESTRY_INFERENCE.out.versions)

    // Filter ancestry groups with sufficient sample size
    ch_ancestry_filtered = ch_stratified_genotypes
        .filter { meta, bed, bim, fam ->
            def n_samples = fam.countLines()
            if (n_samples < params.min_ancestry_n) {
                log.warn "Skipping ancestry ${meta.ancestry}: only ${n_samples} samples (minimum: ${params.min_ancestry_n})"
                return false
            }
            return true
        }

    // Create ancestry x trait combinations for parallel GWAS
    ch_gwas_inputs = ch_ancestry_filtered
        .combine(ch_traits)
        .map { meta, bed, bim, fam, trait_info ->
            def new_meta = meta + trait_info
            [new_meta, bed, bim, fam]
        }

    // =========================================================================
    // GWAS BY ANCESTRY GROUP (PARALLELIZED)
    // =========================================================================
    // Prepare phenotype channel with metadata
    ch_phenotypes_with_meta = ch_phenotypes
        .map { pheno -> [[id: 'phenotypes'], pheno] }

    GWAS_WORKFLOW(
        ch_gwas_inputs,
        ch_phenotypes_with_meta,
        params.covariate_cols,
        params.gwas_tool,
        params.gwas_model,
        params.kinship_matrix ? file(params.kinship_matrix) : [],
        params.run_tractor,
        ch_local_ancestry,
        params.tractor_aac_pops,
        params.tractor_lat_pops,
        params.survival_analysis,
        params.time_col,
        params.event_col
    )
    ch_gwas_results = GWAS_WORKFLOW.out.summary_stats
    ch_gwas_filtered = GWAS_WORKFLOW.out.filtered_results
    ch_tractor_results = GWAS_WORKFLOW.out.tractor_results
    ch_survival_results = GWAS_WORKFLOW.out.survival_results
    ch_versions = ch_versions.mix(GWAS_WORKFLOW.out.versions)

    // =========================================================================
    // VISUALIZATION (Manhattan, QQ plots - run before fine-mapping per spec)
    // =========================================================================
    if (params.visualization) {
        VISUALIZATION(
            ch_gwas_results,
            ch_gwas_filtered,
            params.gene_annotation ? file(params.gene_annotation) : []
        )
        ch_multiqc_files = ch_multiqc_files.mix(VISUALIZATION.out.plots)
        ch_versions = ch_versions.mix(VISUALIZATION.out.versions)
    }

    // =========================================================================
    // META-ANALYSIS ACROSS ANCESTRIES (MR-MEGA)
    // =========================================================================
    if (params.meta_analysis) {
        // Group GWAS results by trait for meta-analysis
        ch_meta_inputs = ch_gwas_results
            .map { meta, sumstats ->
                [[trait: meta.trait, binary: meta.binary], meta.ancestry, sumstats]
            }
            .groupTuple(by: 0)
            .map { trait_meta, ancestries, sumstats_list ->
                [trait_meta, ancestries, sumstats_list]
            }

        META_ANALYSIS(
            ch_meta_inputs,
            params.meta_method,
            params.mr_mega_env,
            params.env_variable
        )
        ch_meta_results = META_ANALYSIS.out.meta_results
        ch_heterogeneity = META_ANALYSIS.out.heterogeneity
        ch_versions = ch_versions.mix(META_ANALYSIS.out.versions)
    }

    // =========================================================================
    // FINE-MAPPING
    // =========================================================================
    if (params.fine_mapping) {
        // Determine which fine-mapping method to use based on trait relationships
        ch_within_ancestry_fm = ch_gwas_filtered
            .map { meta, sumstats ->
                [meta, sumstats, params.within_ancestry_fm]
            }

        // For multi-ancestry fine-mapping, choose method based on trait relatedness
        ch_multi_ancestry_fm_method = params.related_traits ?
            params.multi_ancestry_fm :
            params.single_trait_fm

        // Group results by trait for multi-ancestry fine-mapping
        ch_multi_anc_inputs = ch_gwas_filtered
            .map { meta, sumstats ->
                [[trait: meta.trait], meta.ancestry, sumstats]
            }
            .groupTuple(by: 0)

        FINE_MAPPING(
            ch_within_ancestry_fm,
            ch_multi_anc_inputs,
            ch_multi_ancestry_fm_method,
            params.ld_reference_dir ? file(params.ld_reference_dir) : [],
            params.polyfun_annotations ? file(params.polyfun_annotations) : [],
            params.fm_window_kb,
            params.fm_max_causal
        )
        ch_credible_sets = FINE_MAPPING.out.credible_sets
        ch_fm_results = FINE_MAPPING.out.fine_mapping_results
        ch_versions = ch_versions.mix(FINE_MAPPING.out.versions)
    }

    // =========================================================================
    // COLOCALIZATION ANALYSIS
    // =========================================================================
    if (params.colocalization) {
        // Prepare QTL datasets channel
        ch_qtl_data = Channel.empty()
        if (params.eqtl_datasets) {
            ch_qtl_data = ch_qtl_data.mix(
                Channel.fromPath(params.eqtl_datasets)
                    .splitCsv(header: true)
                    .map { row -> [row.name, 'eqtl', file(row.path)] }
            )
        }
        if (params.custom_qtl_datasets) {
            ch_qtl_data = ch_qtl_data.mix(
                Channel.fromPath(params.custom_qtl_datasets)
                    .splitCsv(header: true)
                    .map { row -> [row.name, row.type, file(row.path)] }
            )
        }

        COLOCALIZATION(
            params.fine_mapping ? ch_fm_results : ch_gwas_filtered,
            ch_qtl_data,
            params.coloc_method,
            params.coloc_p1,
            params.coloc_p2,
            params.coloc_p12
        )
        ch_coloc_results = COLOCALIZATION.out.coloc_results
        ch_versions = ch_versions.mix(COLOCALIZATION.out.versions)
    }

    // =========================================================================
    // HERITABILITY ANALYSIS
    // =========================================================================
    if (params.heritability) {
        HERITABILITY(
            ch_gwas_results,
            params.meta_analysis ? ch_meta_results : Channel.empty(),
            params.h2_method,
            params.ancestry_specific_h2,
            params.gxe_heritability,
            params.genetic_correlation,
            params.ld_reference_dir ? file(params.ld_reference_dir) : []
        )
        ch_h2_results = HERITABILITY.out.h2_estimates
        ch_rg_results = HERITABILITY.out.genetic_correlations
        ch_versions = ch_versions.mix(HERITABILITY.out.versions)
    }

    // =========================================================================
    // PRS ANALYSIS
    // =========================================================================
    if (params.prs_analysis) {
        // Prepare inputs for different PRS methods
        ch_prs_inputs = params.meta_analysis ?
            ch_meta_results.map { meta, results -> [meta, results, 'meta'] } :
            ch_gwas_results.map { meta, results -> [meta, results, meta.ancestry] }

        PRS_WORKFLOW(
            ch_prs_inputs,
            ch_stratified_genotypes,
            ch_local_ancestry,
            params.prs_methods.split(',') as List,
            params.prscsx_reference ? file(params.prscsx_reference) : [],
            params.ld_reference_dir ? file(params.ld_reference_dir) : [],
            params.validation_cohort ? file(params.validation_cohort) : [],
            params.prs_validation,
            params.prs_best_method
        )
        ch_prs_scores = PRS_WORKFLOW.out.prs_scores
        ch_prs_validation = PRS_WORKFLOW.out.validation_results
        ch_prs_comparison = PRS_WORKFLOW.out.method_comparison
        ch_versions = ch_versions.mix(PRS_WORKFLOW.out.versions)
    }

    // =========================================================================
    // FUNCTIONAL ANNOTATION
    // =========================================================================
    if (params.functional_annotation) {
        // Use meta-analysis results if available, otherwise use per-ancestry results
        ch_annot_inputs = params.meta_analysis ? ch_meta_results : ch_gwas_filtered

        FUNCTIONAL_ANNOT(
            ch_annot_inputs,
            params.annotation_tools.split(',') as List,
            params.magma_geneset ? file(params.magma_geneset) : [],
            params.magma_annotation ? file(params.magma_annotation) : [],
            params.lava_partition ? file(params.lava_partition) : [],
            params.cell_type_analysis,
            params.cell_type_datasets ? file(params.cell_type_datasets) : []
        )
        ch_functional_results = FUNCTIONAL_ANNOT.out.annotation_results
        ch_versions = ch_versions.mix(FUNCTIONAL_ANNOT.out.versions)
    }

    // =========================================================================
    // GxG (EPISTASIS) INTERACTION TESTING ON GWAS HITS - END OF PIPELINE
    // =========================================================================
    // Takes hits from ALL GWAS variants (per-ancestry, Tractor, meta-analysis)
    // plus known leukemia risk loci; tests pairwise SNP x SNP interactions in
    // the pooled cohort and within each ancestry stratum, then tests whether
    // the interaction itself differs by ancestry background.
    if (params.gxg_interaction) {
        ch_gxg_sumstats = ch_gwas_results
            .mix(params.run_tractor ? ch_tractor_results : Channel.empty())
            .mix(params.meta_analysis ? ch_meta_results : Channel.empty())
            .mix(params.survival_analysis ? ch_survival_results : Channel.empty())
            .map { meta, ss -> [meta + [survival: meta.survival ?: false,
                                        time_col: meta.time_col ?: params.time_col,
                                        event_col: meta.event_col ?: params.event_col], ss] }

        GXG_INTERACTION(
            ch_gxg_sumstats,
            ch_qc_genotypes,
            ch_phenotypes_with_meta,
            params.run_tractor ? ch_local_ancestry : Channel.empty(),
            (params.gxg_use_default_loci && params.gxg_known_loci) ? file(params.gxg_known_loci) : [],
            params.gxg_custom_variants ? file(params.gxg_custom_variants) : [],
            params.gxg_strata.split(',') as List,
            params.gxg_p_threshold,
            params.gxg_max_hits,
            params.gxg_covariates ?: params.covariate_cols,
            params.gxg_ancestry_col,
            params.min_stratum_n
        )
        ch_gxg_summary = GXG_INTERACTION.out.summary
        ch_gxg_top = GXG_INTERACTION.out.top
        ch_versions = ch_versions.mix(GXG_INTERACTION.out.versions)
    }

    // =========================================================================
    // COLLECT SOFTWARE VERSIONS
    // =========================================================================
    SOFTWARE_VERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
    ch_multiqc_files = ch_multiqc_files.mix(SOFTWARE_VERSIONS.out.yml)

    // =========================================================================
    // GENERATE FINAL REPORT
    // =========================================================================
    REPORTING(
        ch_multiqc_files.collect(),
        ch_gwas_results.collect{ it[1] },
        params.meta_analysis ? ch_meta_results.collect{ it[1] } : Channel.empty().collect(),
        params.fine_mapping ? ch_fm_results.collect{ it[1] } : Channel.empty().collect(),
        params.heritability ? ch_h2_results.collect{ it[1] } : Channel.empty().collect(),
        params.prs_analysis ? ch_prs_comparison.collect{ it[1] } : Channel.empty().collect()
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION HANDLERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    log.info ""
    log.info "========================================="
    log.info "Pipeline completed successfully!"
    log.info "========================================="
    log.info "Started     : ${workflow.start}"
    log.info "Completed   : ${workflow.complete}"
    log.info "Duration    : ${workflow.duration}"
    log.info "Results     : ${params.outdir}"
    log.info ""
}

workflow.onError {
    log.error "Pipeline failed. Please check logs for details."
}
