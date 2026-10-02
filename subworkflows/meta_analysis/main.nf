/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    META_ANALYSIS SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Across-strata meta-analysis per trait:
      'mr-mega' : MR-MEGA (axes of genetic variation; + MR-MEGA-env optional)
      'ivw'     : bin/meta_analyze_strata.R (random-effects IVW per ancestry
                  background; keeps the Tractor-GENESIS structure; same as
                  SLURM step 2)
    METAL / METASOFT are not bundled: selecting them raises an error.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MR_MEGA     } from '../../modules/local/mr_mega'
include { MR_MEGA_ENV } from '../../modules/local/mr_mega_env'
include { META_STRATA } from '../../modules/local/meta_strata'

workflow META_ANALYSIS {
    take:
    ch_gwas_grouped    // channel: [ [trait, binary], [ancestries], [sumstats_files] ]
    meta_method        // value: 'mr-mega' | 'ivw'
    run_mr_mega_env    // boolean
    env_variable       // value: environmental variable column (MR-MEGA-env)

    main:
    ch_versions = Channel.empty()
    ch_meta_results = Channel.empty()
    ch_heterogeneity = Channel.empty()
    ch_env_results = Channel.empty()

    ch_multi = ch_gwas_grouped.filter { trait_meta, ancestries, ss -> ancestries.size() >= 2 }

    if (meta_method == 'mr-mega') {
        MR_MEGA(ch_multi)
        ch_meta_results = MR_MEGA.out.meta_results
        ch_heterogeneity = MR_MEGA.out.heterogeneity
        ch_versions = ch_versions.mix(MR_MEGA.out.versions)

        if (run_mr_mega_env && env_variable) {
            MR_MEGA_ENV(ch_multi, env_variable)
            ch_env_results = MR_MEGA_ENV.out.meta_results
            ch_versions = ch_versions.mix(MR_MEGA_ENV.out.versions)
        }
    } else if (meta_method == 'ivw') {
        META_STRATA(ch_multi, params.meta_ivw_model ?: 'random')
        ch_meta_results = META_STRATA.out.meta_results
        ch_versions = ch_versions.mix(META_STRATA.out.versions)
    } else {
        error "meta_method = '${meta_method}' is not bundled. Use 'mr-mega' or 'ivw'."
    }

    emit:
    meta_results     = ch_meta_results      // channel: [ meta(trait), meta_sumstats ]
    heterogeneity    = ch_heterogeneity     // channel: [ meta(trait), heterogeneity ] (MR-MEGA only)
    ancestry_effects = ch_heterogeneity
    env_results      = ch_env_results       // channel: [ meta(trait), mr_mega_env ]
    versions         = ch_versions
}
