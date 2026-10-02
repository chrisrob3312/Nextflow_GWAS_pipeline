/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HERITABILITY SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LDSC h2 per stratum and for the meta/POOLED GWAS (munged ONCE), cross-
    ancestry genetic correlation (LDSC rg on every stratum pair), and a
    combined comparison table.
    GCTA-GREML / BOLT-REML are not bundled: selecting them raises an error.
    GxE heritability: use MR-MEGA-env (meta_analysis) rather than a separate
    GxE REML; a warning is logged if gxe_heritability is requested here.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LDSC_MUNGE    } from '../../modules/local/ldsc_munge'
include { LDSC_H2       } from '../../modules/local/ldsc_h2'
include { LDSC_RG       } from '../../modules/local/ldsc_rg'
include { H2_COMPARISON } from '../../modules/local/h2_comparison'

workflow HERITABILITY {
    take:
    ch_gwas_results        // channel: [ meta(trait, ancestry), sumstats ]
    ch_meta_results        // channel: [ meta(trait), meta_sumstats ]  (may be empty)
    h2_method              // value: 'ldsc'
    ancestry_specific_h2   // boolean
    gxe_heritability       // boolean (logged only)
    genetic_correlation    // boolean
    ld_reference_dir       // file: LD score directory

    main:
    ch_versions = Channel.empty()

    if (h2_method != 'ldsc') {
        error "h2_method = '${h2_method}' is not bundled. Use 'ldsc' (cov-LDSC weights for admixed strata via bin/estimate_heritability.R on the SLURM path)."
    }
    if (gxe_heritability) {
        log.warn "gxe_heritability = true: no GxE REML is bundled; MR-MEGA-env (meta_analysis) provides the environment-interaction test."
    }

    // Munge once: per-stratum GWAS (if requested) + meta/POOLED tagged META
    ch_to_munge = (ancestry_specific_h2 ? ch_gwas_results : Channel.empty())
        .mix(ch_meta_results.map { meta, ss -> [meta + [ancestry: meta.ancestry ?: 'META'], ss] })

    LDSC_MUNGE(ch_to_munge)
    ch_versions = ch_versions.mix(LDSC_MUNGE.out.versions)

    LDSC_H2(LDSC_MUNGE.out.munged, ld_reference_dir)
    ch_h2 = LDSC_H2.out.results
    ch_versions = ch_versions.mix(LDSC_H2.out.versions)

    // Cross-ancestry rg on every pair of strata within a trait
    ch_rg = Channel.empty()
    if (genetic_correlation) {
        ch_pairs = LDSC_MUNGE.out.munged
            .filter { meta, m -> meta.ancestry != 'META' }
            .map { meta, m -> [[trait: meta.trait], meta.ancestry, m] }
            .groupTuple(by: 0)
            .filter { t, ancs, files -> ancs.size() >= 2 }
            .flatMap { t, ancs, files ->
                def pairs = []
                for (int i = 0; i < ancs.size() - 1; i++)
                    for (int j = i + 1; j < ancs.size(); j++)
                        pairs << [t + [anc1: ancs[i], anc2: ancs[j], id: "${t.trait}.${ancs[i]}_vs_${ancs[j]}"], files[i], files[j]]
                pairs
            }
        LDSC_RG(ch_pairs, ld_reference_dir)
        ch_rg = LDSC_RG.out.results
        ch_versions = ch_versions.mix(LDSC_RG.out.versions)
    }

    H2_COMPARISON(
        ch_h2.map { m, f -> f }.collect(),
        ch_rg.map { m, f -> f }.collect().ifEmpty([])
    )
    ch_versions = ch_versions.mix(H2_COMPARISON.out.versions)

    emit:
    h2_estimates         = ch_h2                          // channel: [ meta, h2_results.tsv ]
    genetic_correlations = ch_rg                          // channel: [ meta(anc1, anc2), rg_results.tsv ]
    gxe_estimates        = Channel.empty()
    comparison           = H2_COMPARISON.out.comparison   // channel: [ h2_comparison.tsv ]
    versions             = ch_versions
}
