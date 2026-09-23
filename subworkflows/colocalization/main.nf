/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLOCALIZATION SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Runs the parallel multi-method colocalization (bin/run_colocalization.R:
    coloc.susie per QTL type | HyPrColoc across QTL types | OPERA on all QTLs)
    on BOTH:
      1. ancestry-stratified GWAS (population-specific regulation)
      2. meta-analysis / POOLED GWAS (power for shared causal variants)
    against the pooled QTL megaset (eQTL, sQTL, pQTL, mQTL, caQTL, hQTL),
    then builds a per-trait consensus and the shared-vs-divergent analysis.

    Aligned with the SLURM path (slurm/submit_full_pipeline.sh step 4), which
    calls the same script.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PREPARE_QTL_DATA           } from '../../modules/local/prepare_qtl_data'
include { POOL_QTL_DATASETS          } from '../../modules/local/prepare_qtl_data'
include { COLOC_RUN                  } from '../../modules/local/colocalization'
include { COLOC_COMBINE              } from '../../modules/local/colocalization'
include { COLOC_DIVERGENCE_ANALYSIS  } from '../../modules/local/coloc_divergence'
include { COLOC_ANCESTRY_HEATMAP     } from '../../modules/local/coloc_divergence'

workflow COLOCALIZATION {
    take:
    ch_ancestry_gwas   // channel: [ meta(trait, ancestry, binary), sumstats ]  ancestry-stratified (incl. Tractor-GENESIS)
    ch_meta_gwas       // channel: [ meta(trait), sumstats ]                    meta-analysis / POOLED
    ch_qtl_data        // channel: [ name, type, qtl_file ]
    ch_cohort_ld       // channel: [ meta, ld_rds ] cohort-specific LD, or empty
    gwas_n             // value: GWAS sample size (fallback when meta.n is absent)
    qtl_n              // value: QTL sample size
    p1                 // value: prior for trait association
    p2                 // value: prior for QTL association
    p12                // value: prior for colocalization

    main:
    ch_versions = Channel.empty()

    // ------------------------------------------------------------------
    // QTL megaset: standardize each source, pool by QTL type
    // ------------------------------------------------------------------
    PREPARE_QTL_DATA(ch_qtl_data)
    ch_versions = ch_versions.mix(PREPARE_QTL_DATA.out.versions)

    ch_qtl_by_type = PREPARE_QTL_DATA.out.standardized
        .map { name, type, qtl_file -> [type, name, qtl_file] }
        .groupTuple(by: 0)

    POOL_QTL_DATASETS(ch_qtl_by_type)
    ch_versions = ch_versions.mix(POOL_QTL_DATASETS.out.versions)

    // All pooled QTL types as ONE list -> every GWAS sees the whole megaset
    ch_qtl_megaset = POOL_QTL_DATASETS.out.pooled
        .toList()
        .map { rows -> [rows.collect { it[1] }, rows.collect { it[0] }] }   // [files, types]

    // ------------------------------------------------------------------
    // GWAS inputs: ancestry-stratified + meta, tagged
    // ------------------------------------------------------------------
    ch_gwas_all = ch_ancestry_gwas
        .map { meta, ss -> [meta + [analysis_type: 'ancestry_stratified'], ss] }
        .mix(ch_meta_gwas.map { meta, ss -> [meta + [analysis_type: 'meta_analysis', ancestry: meta.ancestry ?: 'META'], ss] })

    // Optional cohort LD keyed by ancestry (empty -> ABF instead of SuSiE)
    ch_ld_by_anc = ch_cohort_ld.map { meta, ld -> [meta.ancestry ?: 'ALL', ld] }
    ch_coloc_input = ch_gwas_all
        .map { meta, ss -> [meta.ancestry ?: 'ALL', meta, ss] }
        .combine(ch_ld_by_anc.ifEmpty { ['__none__', []] }, by: 0)
        .map { anc, meta, ss, ld -> [meta, ss, ld] }
        .mix(  // GWAS with no matching LD entry still run (ABF)
            ch_gwas_all.map { meta, ss -> [meta, ss, []] }
        )
        .unique { it[0] }
        .combine(ch_qtl_megaset)
        .map { meta, ss, ld, qtl_files, qtl_types -> [meta, ss, qtl_files, qtl_types, ld] }

    COLOC_RUN(ch_coloc_input, gwas_n, qtl_n, p1, p2, p12)
    ch_versions = ch_versions.mix(COLOC_RUN.out.versions)

    ch_all_results = COLOC_RUN.out.coloc_susie
        .mix(COLOC_RUN.out.hyprcoloc)
        .mix(COLOC_RUN.out.opera)

    // ------------------------------------------------------------------
    // Per-trait consensus across strata and methods
    // ------------------------------------------------------------------
    ch_by_trait = ch_all_results
        .map { meta, f -> [[trait: meta.trait], f] }
        .groupTuple(by: 0)

    COLOC_COMBINE(ch_by_trait)
    ch_versions = ch_versions.mix(COLOC_COMBINE.out.versions)

    // ------------------------------------------------------------------
    // Shared vs divergent (ancestry-stratified vs meta) on coloc.susie PP4
    // ------------------------------------------------------------------
    ch_cs_ancestry = COLOC_RUN.out.coloc_susie
        .filter { meta, f -> meta.analysis_type == 'ancestry_stratified' }
        .map { meta, f -> [meta.trait, f] }
        .groupTuple()
    ch_cs_meta = COLOC_RUN.out.coloc_susie
        .filter { meta, f -> meta.analysis_type == 'meta_analysis' }
        .map { meta, f -> [meta.trait, f] }

    ch_divergence_input = ch_cs_ancestry.join(ch_cs_meta)

    COLOC_DIVERGENCE_ANALYSIS(
        ch_divergence_input.map { it[1] },
        ch_divergence_input.map { it[2] },
        ch_divergence_input.map { it[0] },
        0.8
    )
    ch_versions = ch_versions.mix(COLOC_DIVERGENCE_ANALYSIS.out.versions)

    ch_heatmap_input = COLOC_RUN.out.coloc_susie.map { meta, f -> [meta.trait, f] }.groupTuple()
    COLOC_ANCESTRY_HEATMAP(
        ch_heatmap_input.map { trait, files -> files },
        ch_heatmap_input.map { trait, files -> trait }
    )
    ch_versions = ch_versions.mix(COLOC_ANCESTRY_HEATMAP.out.versions)

    emit:
    coloc_results      = ch_all_results                           // channel: [ meta, results ] (all methods)
    coloc_susie        = COLOC_RUN.out.coloc_susie                // channel: [ meta, coloc_susie.tsv ]
    hyprcoloc          = COLOC_RUN.out.hyprcoloc                  // channel: [ meta, hyprcoloc.tsv ]
    opera              = COLOC_RUN.out.opera                      // channel: [ meta, opera.tsv ]
    combined           = COLOC_COMBINE.out.combined               // channel: [ meta(trait), combined ]
    consensus          = COLOC_COMBINE.out.consensus              // channel: [ meta(trait), consensus ]
    method_comparison  = COLOC_COMBINE.out.comparison
    shared_coloc       = COLOC_DIVERGENCE_ANALYSIS.out.shared
    divergent_coloc    = COLOC_DIVERGENCE_ANALYSIS.out.divergent
    meta_only_coloc    = COLOC_DIVERGENCE_ANALYSIS.out.meta_only
    divergence_summary = COLOC_DIVERGENCE_ANALYSIS.out.summary
    ancestry_heatmap   = COLOC_ANCESTRY_HEATMAP.out.heatmap
    pooled_qtl         = POOL_QTL_DATASETS.out.pooled             // channel: [ type, pooled_qtl ]
    versions           = ch_versions
}
