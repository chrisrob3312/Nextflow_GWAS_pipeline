/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FINE_MAPPING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Within-ancestry : regions around GWAS signals -> PolyFun+SuSiE (when
                      annotations are given) or SuSiE-RSS with LD from a PLINK
                      reference (cohort-specific per ancestry, or a panel)
    Multi-ancestry  : MG-FLASH-FM (related traits) or SuSiE-ME
    Credible sets from every method are merged per GWAS.
    FINEMAP and CAVIAR are not bundled: selecting them raises an error.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { EXTRACT_REGIONS } from '../../modules/local/fine_mapping_utils'
include { SUSIE_RSS       } from '../../modules/local/fine_mapping_utils'
include { CREDIBLE_SETS   } from '../../modules/local/fine_mapping_utils'
include { POLYFUN_SUSIE   } from '../../modules/local/polyfun_susie'
include { MG_FLASH_FM     } from '../../modules/local/mg_flash_fm'
include { SUSIE_ME        } from '../../modules/local/susie_me'

workflow FINE_MAPPING {
    take:
    ch_within_ancestry_input   // channel: [ meta, sumstats, method ]
    ch_multi_ancestry_input    // channel: [ [trait], [ancestries], [sumstats_files] ]
    multi_ancestry_method      // value: 'mg-flash-fm' or 'susie-me'
    ld_reference_dir           // file: directory with <ancestry>.bed/bim/fam (or [])
    polyfun_annotations        // file: PolyFun annotation files (or [])
    window_kb                  // value
    max_causal                 // value

    main:
    ch_versions = Channel.empty()
    ch_within_results = Channel.empty()
    ch_multi_results = Channel.empty()
    ch_cs = Channel.empty()

    // ------------------------------------------------------------------
    // Within-ancestry
    // ------------------------------------------------------------------
    ch_within_ancestry_input
        .map { meta, ss, method -> [meta + [fm_method: method], ss] }
        .set { ch_within }

    EXTRACT_REGIONS(ch_within, window_kb, params.fm_region_p ?: 5e-8)
    ch_versions = ch_versions.mix(EXTRACT_REGIONS.out.versions)

    EXTRACT_REGIONS.out.region_sumstats
        .branch {
            polyfun: it[0].fm_method == 'polyfun-susie' && polyfun_annotations
            susie:   it[0].fm_method in ['susie', 'polyfun-susie']   // polyfun without annotations -> SuSiE
            other:   true
        }
        .set { ch_fm }

    ch_fm.other.subscribe { meta, files ->
        error "within_ancestry_fm = '${meta.fm_method}' is not bundled. Use 'polyfun-susie' (needs --polyfun_annotations) or 'susie'."
    }

    if (polyfun_annotations) {
        POLYFUN_SUSIE(ch_fm.polyfun, ld_reference_dir, max_causal)
        ch_within_results = ch_within_results.mix(POLYFUN_SUSIE.out.results)
        ch_cs = ch_cs.mix(POLYFUN_SUSIE.out.credible_sets)
        ch_versions = ch_versions.mix(POLYFUN_SUSIE.out.versions)
    }

    SUSIE_RSS(ch_fm.susie, ld_reference_dir, max_causal)
    ch_within_results = ch_within_results.mix(SUSIE_RSS.out.results)
    ch_cs = ch_cs.mix(SUSIE_RSS.out.credible_sets)
    ch_versions = ch_versions.mix(SUSIE_RSS.out.versions)

    // ------------------------------------------------------------------
    // Multi-ancestry
    // ------------------------------------------------------------------
    ch_multi = ch_multi_ancestry_input.filter { trait_meta, ancestries, ss -> ancestries.size() >= 2 }

    if (multi_ancestry_method == 'mg-flash-fm') {
        MG_FLASH_FM(ch_multi, ld_reference_dir, window_kb, max_causal)
        ch_multi_results = MG_FLASH_FM.out.results
        ch_cs = ch_cs.mix(MG_FLASH_FM.out.credible_sets)
        ch_versions = ch_versions.mix(MG_FLASH_FM.out.versions)
    } else if (multi_ancestry_method == 'susie-me') {
        SUSIE_ME(ch_multi, ld_reference_dir, max_causal)
        ch_multi_results = SUSIE_ME.out.results
        ch_cs = ch_cs.mix(SUSIE_ME.out.credible_sets)
        ch_versions = ch_versions.mix(SUSIE_ME.out.versions)
    } else {
        error "multi_ancestry_fm = '${multi_ancestry_method}' is not bundled. Use 'mg-flash-fm' or 'susie-me'."
    }

    // ------------------------------------------------------------------
    // Credible sets merged per trait
    // ------------------------------------------------------------------
    ch_cs_by_trait = ch_cs
        .map { meta, f -> [[trait: meta.trait], f] }
        .groupTuple(by: 0)

    CREDIBLE_SETS(ch_cs_by_trait)
    ch_versions = ch_versions.mix(CREDIBLE_SETS.out.versions)

    emit:
    fine_mapping_results = ch_within_results.mix(ch_multi_results)   // channel: [ meta, results ]
    credible_sets        = CREDIBLE_SETS.out.credible_sets            // channel: [ meta(trait), credible_sets.tsv ]
    regions              = EXTRACT_REGIONS.out.regions                // channel: [ meta, regions.tsv ]
    within_ancestry      = ch_within_results
    multi_ancestry       = ch_multi_results
    versions             = ch_versions
}
