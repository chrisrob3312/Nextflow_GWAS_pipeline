/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GWAS_WORKFLOW SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ancestry-stratified GWAS using multiple tools
    Supports: REGENIE, SAIGE, BOLT-LMM, PLINK2, GENESIS, SPA-Cox
    Admixed-optimized: GENESIS, SAIGE, Tractor (local ancestry-aware)

    Tractor combines LAT1+LAT2 for Latino analysis (3-way: EUR,AFR,NAT)
    Tractor runs AAC separately (2-way: EUR,AFR)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REGENIE_STEP1       } from '../../modules/local/regenie_step1'
include { REGENIE_STEP2       } from '../../modules/local/regenie_step2'
include { SAIGE_STEP1         } from '../../modules/local/saige_step1'
include { SAIGE_STEP2         } from '../../modules/local/saige_step2'
include { BOLT_LMM            } from '../../modules/local/bolt_lmm'
include { PLINK2_GWAS         } from '../../modules/local/plink2_gwas'
include { GENESIS_NULL_MODEL  } from '../../modules/local/genesis'
include { GENESIS_ASSOC       } from '../../modules/local/genesis'
include { SPA_COX_STEP1       } from '../../modules/local/spa_cox'
include { SPA_COX_STEP2       } from '../../modules/local/spa_cox'
include { TRACTOR_EXTRACT_TRACTS } from '../../modules/local/tractor'
include { TRACTOR_GENESIS     } from '../../modules/local/tractor_genesis_adapter'
include { GWAS_FILTER         } from '../../modules/local/gwas_filter'
include { CLUMP_REGIONS       } from '../../modules/local/clump_regions'
include { PLINK2_TO_GDS       } from '../../modules/local/plink2_to_gds'
include { PLINK_TO_VCF        } from '../../modules/local/format_conversion'
include { MERGE_PLINK_FILES   } from '../../modules/local/format_conversion'
include { VALIDATE_GWAS_INPUT } from '../../modules/local/validation'

workflow GWAS_WORKFLOW {
    take:
    ch_genotypes       // channel: [ meta, bed, bim, fam ] (meta includes ancestry, trait)
    ch_phenotypes      // channel: [ meta, phenotype_file ]
    covariate_cols     // value: covariate column names
    gwas_tool          // value: 'regenie', 'saige', 'bolt-lmm', 'plink2', 'genesis'
    gwas_model         // value: 'additive', 'dominant', 'recessive'
    kinship_matrix     // path: kinship/GRM matrix (optional, for GENESIS)
    run_tractor        // boolean: run Tractor for admixed populations
    ch_local_ancestry  // channel: [ meta, msp_file ] local ancestry calls (for Tractor)
    tractor_aac_pops   // value: ancestral populations for AAC (e.g., "EUR,AFR")
    tractor_lat_pops   // value: ancestral populations for Latino (e.g., "EUR,AFR,NAT")
    survival_analysis  // boolean: run survival/time-to-event analysis
    time_col           // value: time column for survival analysis
    event_col          // value: event column for survival analysis

    main:
    ch_versions = Channel.empty()
    ch_gwas_results = Channel.empty()
    ch_tractor_results = Channel.empty()
    ch_tractor_dosages = Channel.empty()
    ch_survival_results = Channel.empty()

    // =========================================================================
    // INPUT VALIDATION
    // =========================================================================

    VALIDATE_GWAS_INPUT(
        ch_genotypes,
        ch_phenotypes.map { meta, pheno -> pheno }.first()
    )
    ch_versions = ch_versions.mix(VALIDATE_GWAS_INPUT.out.versions)

    // Combine genotypes with phenotypes
    ch_gwas_input = ch_genotypes.combine(ch_phenotypes.map { meta, pheno -> pheno })

    // =========================================================================
    // STANDARD GWAS TOOLS
    // =========================================================================

    if (gwas_tool == 'regenie') {
        // REGENIE Step 1: Fit null model
        // Group by ancestry for Step 1 (fit once per ancestry)
        ch_step1_input = ch_gwas_input
            .map { meta, bed, bim, fam, pheno ->
                [[id: meta.id, ancestry: meta.ancestry], bed, bim, fam, pheno]
            }
            .groupTuple(by: 0)
            .map { meta, beds, bims, fams, phenos ->
                [meta, beds[0], bims[0], fams[0], phenos[0]]
            }

        REGENIE_STEP1(
            ch_step1_input,
            covariate_cols
        )
        ch_versions = ch_versions.mix(REGENIE_STEP1.out.versions)

        // REGENIE Step 2: Association testing (per trait)
        ch_step2_input = ch_gwas_input
            .map { meta, bed, bim, fam, pheno ->
                [[id: meta.id, ancestry: meta.ancestry], meta, bed, bim, fam, pheno]
            }
            .combine(REGENIE_STEP1.out.predictions.map { meta, pred -> [meta, pred] }, by: 0)
            .map { key, meta, bed, bim, fam, pheno, predictions ->
                [meta, bed, bim, fam, pheno, predictions]
            }

        REGENIE_STEP2(
            ch_step2_input,
            covariate_cols,
            gwas_model
        )
        ch_gwas_results = REGENIE_STEP2.out.summary_stats
        ch_versions = ch_versions.mix(REGENIE_STEP2.out.versions)

    } else if (gwas_tool == 'saige') {
        // SAIGE - optimized for admixed populations and case-control imbalance
        ch_step1_input = ch_gwas_input
            .map { meta, bed, bim, fam, pheno ->
                [[id: meta.id, ancestry: meta.ancestry, trait: meta.trait, binary: meta.binary ?: false],
                 bed, bim, fam, pheno]
            }

        SAIGE_STEP1(
            ch_step1_input,
            covariate_cols
        )
        ch_versions = ch_versions.mix(SAIGE_STEP1.out.versions)

        ch_step2_input = ch_step1_input
            .map { meta, bed, bim, fam, pheno -> [meta, bed, bim, fam] }
            .join(SAIGE_STEP1.out.model)

        SAIGE_STEP2(
            ch_step2_input,
            gwas_model
        )
        ch_gwas_results = SAIGE_STEP2.out.summary_stats
        ch_versions = ch_versions.mix(SAIGE_STEP2.out.versions)

    } else if (gwas_tool == 'genesis') {
        // GENESIS - optimized for admixed populations with complex relatedness
        // Requires GDS format and optionally kinship matrix

        // Convert PLINK to GDS
        PLINK2_TO_GDS(
            ch_genotypes
        )
        ch_versions = ch_versions.mix(PLINK2_TO_GDS.out.versions)

        // Prepare input for GENESIS null model
        ch_genesis_input = PLINK2_TO_GDS.out.gds
            .combine(ch_phenotypes.map { meta, pheno -> pheno })

        // Fit null model
        GENESIS_NULL_MODEL(
            ch_genesis_input,
            covariate_cols,
            kinship_matrix ?: []
        )
        ch_versions = ch_versions.mix(GENESIS_NULL_MODEL.out.versions)

        // Run association
        ch_assoc_input = PLINK2_TO_GDS.out.gds
            .join(GENESIS_NULL_MODEL.out.null_model)

        GENESIS_ASSOC(
            ch_assoc_input,
            'Score'  // Test type: Score, Wald, BinomiRare
        )
        ch_gwas_results = GENESIS_ASSOC.out.summary_stats
        ch_versions = ch_versions.mix(GENESIS_ASSOC.out.versions)

    } else if (gwas_tool == 'bolt-lmm') {
        BOLT_LMM(
            ch_gwas_input,
            covariate_cols,
            gwas_model
        )
        ch_gwas_results = BOLT_LMM.out.summary_stats
        ch_versions = ch_versions.mix(BOLT_LMM.out.versions)

    } else if (gwas_tool == 'plink2') {
        PLINK2_GWAS(
            ch_gwas_input,
            covariate_cols,
            gwas_model
        )
        ch_gwas_results = PLINK2_GWAS.out.summary_stats
        ch_versions = ch_versions.mix(PLINK2_GWAS.out.versions)
    }

    // =========================================================================
    // TRACTOR - Local Ancestry-Aware GWAS for Admixed Populations
    // =========================================================================
    // - AAC (African American): 2-way admixture (EUR, AFR)
    // - Latino (LAT1 + LAT2 COMBINED): 3-way admixture (EUR, AFR, NAT)
    //   LAT1 = Mexican/Central American, LAT2 = Caribbean/South American
    //   These are combined for better statistical power in Tractor analysis
    // Decomposes genetic effects by ancestral origin

    // ALL traits go through Tractor-GENESIS (bin/tractor_genesis_adapter.R) so
    // MRD (binary), relapse and OS (time-to-event) are directly comparable:
    // one conditional score test, kinship in the null model for every trait.
    // meta.model / meta.time_col / meta.event_col are set per trait in main.nf.
    //
    // Each Tractor process is invoked ONCE on a mixed AAC + Latino channel
    // (DSL2 forbids calling the same process twice without an alias).

    if (run_tractor && ch_local_ancestry) {
        // =====================================================================
        // AFRICAN AMERICAN (AAC) - 2-way admixture (EUR, AFR)
        // =====================================================================
        ch_aac = ch_genotypes
            .filter { meta, bed, bim, fam -> meta.ancestry == 'AAC' }
            .join(ch_local_ancestry)
            .map { meta, bed, bim, fam, la ->
                [meta + [tractor_pops: tractor_aac_pops, tractor_group: 'AAC',
                         ref_ancestry: params.tractor_ref_ancestry ?: 'EUR'], bed, bim, fam, la]
            }

        // =====================================================================
        // LATINO (LAT1 + LAT2 COMBINED) - 3-way admixture (EUR, AFR, NAT/AMR)
        // =====================================================================
        // LAT1 and LAT2 are merged for the Tractor decomposition (power for the
        // 3-way model); stratified LAT1/LAT2 results come from the standard
        // per-ancestry GWAS above and from the SLURM per-stratum scripts.
        ch_latino_separate = ch_genotypes
            .filter { meta, bed, bim, fam -> meta.ancestry in ['LAT1', 'LAT2', 'AHI'] }
            .join(ch_local_ancestry)

        ch_latino_for_merge = ch_latino_separate
            .map { meta, bed, bim, fam, la ->
                def merge_key = meta.id.replaceAll(/\.(LAT1|LAT2|AHI)$/, '')
                [[merge_id: merge_key, trait: meta.trait], meta, bed, bim, fam, la]
            }
            .groupTuple(by: 0)
            .map { merge_meta, metas, beds, bims, fams, las ->
                def base = metas[0].findAll { k, v -> !(k in ['id', 'ancestry']) }
                def new_meta = base + [id: "${merge_meta.merge_id}.LATINO", ancestry: 'LATINO',
                                       original_ancestries: metas.collect { it.ancestry }.join(','),
                                       tractor_pops: tractor_lat_pops, tractor_group: 'LATINO',
                                       ref_ancestry: params.tractor_ref_ancestry ?: 'EUR']
                [new_meta, beds, bims, fams, las]
            }

        ch_latino_for_merge
            .branch {
                single:   it[1].size() == 1
                multiple: it[1].size() > 1
            }
            .set { ch_latino_branched }

        ch_latino_single = ch_latino_branched.single
            .map { meta, beds, bims, fams, las -> [meta, beds[0], bims[0], fams[0], las[0]] }

        MERGE_PLINK_FILES(
            ch_latino_branched.multiple.map { meta, beds, bims, fams, las -> [meta, beds, bims, fams] }
        )
        ch_latino_merged_la = ch_latino_branched.multiple.map { meta, beds, bims, fams, las -> [meta, las] }

        ch_latino_combined = ch_latino_single.mix(
            MERGE_PLINK_FILES.out.merged
                .join(ch_latino_merged_la)
                .map { meta, bed, bim, fam, las -> [meta, bed, bim, fam, las[0]] }  // LA files combined upstream
        )

        // =====================================================================
        // AAC + LATINO -> one channel -> each process called once
        // =====================================================================
        ch_tractor_groups = ch_aac.mix(ch_latino_combined)

        PLINK_TO_VCF(
            ch_tractor_groups.map { meta, bed, bim, fam, la -> [meta, bed, bim, fam] },
            'tractor'
        )

        ch_tractor_extract_input = PLINK_TO_VCF.out.vcf
            .join(ch_tractor_groups.map { meta, bed, bim, fam, la -> [meta, la] })
            .map { meta, vcf, vcf_idx, la_files -> [meta, vcf, la_files] }

        // Ancestry-specific dosages + haplotype counts (populations read from meta.tractor_pops)
        TRACTOR_EXTRACT_TRACTS(
            ch_tractor_extract_input,
            tractor_lat_pops
        )
        ch_versions = ch_versions.mix(TRACTOR_EXTRACT_TRACTS.out.versions.first())
        ch_tractor_dosages = TRACTOR_EXTRACT_TRACTS.out.ancestry_dosages   // reused by PRS (LA partial scores) and GxG

        // Tractor-GENESIS for every trait: binary, quantitative and time-to-event
        ch_tractor_genesis_input = TRACTOR_EXTRACT_TRACTS.out.ancestry_dosages
            .join(TRACTOR_EXTRACT_TRACTS.out.haplotype_counts)
            .combine(ch_phenotypes.map { meta, pheno -> pheno })

        TRACTOR_GENESIS(
            ch_tractor_genesis_input,
            kinship_matrix,
            covariate_cols
        )
        ch_tractor_results = ch_tractor_results.mix(TRACTOR_GENESIS.out.sumstats)
        ch_versions = ch_versions.mix(TRACTOR_GENESIS.out.versions.first())
    }

    // =========================================================================
    // SURVIVAL / TIME-TO-EVENT ANALYSIS (SPA-Cox)
    // =========================================================================

    if (survival_analysis && time_col && event_col) {
        // Filter to traits with time-to-event data
        ch_survival_input = ch_gwas_input

        SPA_COX_STEP1(
            ch_survival_input,
            covariate_cols,
            time_col,
            event_col
        )
        ch_versions = ch_versions.mix(SPA_COX_STEP1.out.versions)

        ch_cox_step2_input = ch_survival_input
            .map { meta, bed, bim, fam, pheno -> [meta, bed, bim, fam] }
            .join(SPA_COX_STEP1.out.model)

        SPA_COX_STEP2(
            ch_cox_step2_input
        )
        ch_survival_results = SPA_COX_STEP2.out.summary_stats
        ch_versions = ch_versions.mix(SPA_COX_STEP2.out.versions)
    }

    // =========================================================================
    // FILTER AND CLUMP RESULTS
    // =========================================================================

    // Combine all GWAS results
    ch_all_gwas = ch_gwas_results
        .mix(ch_tractor_results)
        .mix(ch_survival_results)

    // Only run filtering if we have results
    ch_all_gwas
        .ifEmpty { log.warn "No GWAS results to filter" }
        .set { ch_gwas_to_filter }

    GWAS_FILTER(
        ch_gwas_to_filter
    )
    ch_versions = ch_versions.mix(GWAS_FILTER.out.versions)

    // Identify independent signals by clumping
    CLUMP_REGIONS(
        GWAS_FILTER.out.significant,
        ch_genotypes.map { meta, bed, bim, fam -> [meta.ancestry, bed, bim, fam] }.unique()
    )
    ch_versions = ch_versions.mix(CLUMP_REGIONS.out.versions)

    emit:
    summary_stats     = ch_gwas_results                    // channel: [ meta, sumstats ]
    tractor_results   = ch_tractor_results                 // channel: [ meta, tractor_genesis.tsv.gz ]
    tractor_dosages   = ch_tractor_dosages                 // channel: [ meta, ancdose files ] (PRS LA-partial, GxG)
    survival_results  = ch_survival_results                // channel: [ meta, survival_sumstats ]
    filtered_results  = GWAS_FILTER.out.filtered           // channel: [ meta, filtered_sumstats ]
    significant       = GWAS_FILTER.out.significant        // channel: [ meta, sig_variants ]
    clumped           = CLUMP_REGIONS.out.clumped          // channel: [ meta, clumped_results ]
    versions          = ch_versions                        // channel: [ versions.yml ]
}
