/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GXG_INTERACTION SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    End-of-pipeline epistasis testing on GWAS hits.

    1. Select hits from all GWAS outputs for the trait (Tractor, standard,
       meta-analysis) plus known leukemia risk loci
    2. Test every pairwise SNP x SNP interaction:
         - POOLED cohort
         - within each ancestry stratum (EUR, AAC, LAT1, LAT2, EAS/SAS if N>=30, OTHER)
       (one task per stratum -> SLURM array)
    3. Ancestry modification: Cochran's Q across strata + pooled 3-way LRT
    4. Combine into one summary table per trait
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GXG_SELECT_HITS            } from '../../modules/local/gxg_interaction'
include { GXG_TEST                   } from '../../modules/local/gxg_interaction'
include { GXG_ANCESTRY_HETEROGENEITY } from '../../modules/local/gxg_interaction'
include { GXG_COMBINE                } from '../../modules/local/gxg_interaction'

workflow GXG_INTERACTION {
    take:
    ch_sumstats        // channel: [ meta(trait, ancestry, binary, ...), sumstats ]  (all GWAS types)
    ch_genotypes       // channel: [ meta, bed, bim, fam ]  (full QC'd cohort, NOT stratified)
    ch_phenotypes      // channel: [ meta, phenotype_file ]
    ch_tractor_files   // channel: [ meta, tractor_dosage_files ] or empty
    known_loci         // file or []
    strata             // list: ['POOLED','EUR','AAC','LAT1','LAT2','EAS','SAS','OTHER']
    p_threshold        // numeric: hit selection threshold
    max_hits           // integer
    covariates         // string
    ancestry_col       // string
    min_stratum_n      // integer

    main:
    ch_versions = Channel.empty()

    // Gather every sumstats file for a trait into one list
    ch_by_trait = ch_sumstats
        .map { meta, ss -> [[trait: meta.trait, binary: meta.binary ?: false,
                             survival: meta.survival ?: false,
                             time_col: meta.time_col ?: params.time_col,
                             event_col: meta.event_col ?: params.event_col], ss] }
        .groupTuple(by: 0)

    GXG_SELECT_HITS(ch_by_trait, known_loci, p_threshold, max_hits)
    ch_versions = ch_versions.mix(GXG_SELECT_HITS.out.versions)

    // Full-cohort genotypes + phenotype attached to each trait
    ch_geno_single = ch_genotypes.first()
    ch_pheno_single = ch_phenotypes.map { m, p -> p }.first()

    ch_test_input = GXG_SELECT_HITS.out.snp_list
        .combine(ch_geno_single)
        .combine(ch_pheno_single)
        .map { meta, snps, gmeta, bed, bim, fam, pheno -> [meta, snps, bed, bim, fam, pheno] }

    ch_tractor = ch_tractor_files.map { m, f -> f }.collect().ifEmpty([])

    // One task per stratum (POOLED + each ancestry) -> array-friendly
    GXG_TEST(
        ch_test_input,
        strata,
        known_loci,
        ch_tractor,
        covariates,
        ancestry_col,
        min_stratum_n
    )
    ch_versions = ch_versions.mix(GXG_TEST.out.versions)

    // Ancestry modification of interaction effects (needs all strata in one process)
    GXG_ANCESTRY_HETEROGENEITY(
        ch_test_input,
        known_loci,
        covariates,
        ancestry_col,
        min_stratum_n
    )
    ch_versions = ch_versions.mix(GXG_ANCESTRY_HETEROGENEITY.out.versions)

    // Combine per trait
    ch_strata_results = GXG_TEST.out.results
        .map { meta, stratum, files -> [meta, files] }
        .groupTuple(by: 0)
        .map { meta, files -> [meta, files.flatten()] }

    ch_het = GXG_ANCESTRY_HETEROGENEITY.out.heterogeneity
        .ifEmpty { [[trait: 'none'], file('NO_FILE')] }

    GXG_COMBINE(ch_strata_results, ch_het)
    ch_versions = ch_versions.mix(GXG_COMBINE.out.versions)

    emit:
    hits            = GXG_SELECT_HITS.out.hits_table
    stratum_results = GXG_TEST.out.results
    heterogeneity   = GXG_ANCESTRY_HETEROGENEITY.out.heterogeneity
    summary         = GXG_COMBINE.out.summary
    top             = GXG_COMBINE.out.top
    versions        = ch_versions
}
