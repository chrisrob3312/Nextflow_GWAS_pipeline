/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRS_WORKFLOW SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Multi-method PRS for admixed populations (bin/calculate_prs_admixed.R):
      primary PRS-CSx; adjuncts GAUDI, DiscoDivas, SDPR_admix, MUSSEL, PROSPER
      (no clumping + thresholding).
    Local-ancestry PARTIAL scores from Tractor dosages x ancestry-specific
    effects (Tractor-GENESIS BETA_<anc>): the part of each person's score
    carried on EUR / AFR / AMR haplotypes.
    Validation on the phenotype, OVERALL and WITHIN each ancestry stratum,
    with disparity flags; best method overall and per stratum.
    Aligned with the SLURM path (slurm/submit_prs_array.sh), same script.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PRS_METHOD              } from '../../modules/local/prs_admixed'
include { PRS_LA_PARTIAL          } from '../../modules/local/prs_admixed'
include { PRS_VALIDATE_STRATIFIED } from '../../modules/local/prs_admixed'
include { PRS_COMBINE_METHODS     } from '../../modules/local/prs_admixed'

workflow PRS_WORKFLOW {
    take:
    ch_sumstats          // channel: [ meta(trait, ancestry, binary, survival, ...), sumstats ]  ancestry-specific GWAS
    ch_genotypes         // channel: [ meta, bed, bim, fam ]   full QC'd cohort (scoring + validation)
    ch_local_ancestry    // channel: [ meta, msp ] or empty
    ch_phenotype         // channel: [ meta, phenotype_with_pcs ]
    ch_tractor_sumstats  // channel: [ meta(trait), tractor_genesis.tsv.gz ] ancestry-specific weights, or empty
    ch_tractor_dosages   // channel: [ meta, ancdose files ] or empty
    prs_methods          // list: e.g. ['prs_csx','gaudi','disco_divas','sdpr_admix','mussel','prosper']
    ld_ref               // file: LD reference dir (PRS-CSx / SDPR), or []
    run_validation       // boolean
    ancestry_col         // string

    main:
    ch_versions = Channel.empty()

    ch_geno_single  = ch_genotypes.first()
    ch_pheno_single = ch_phenotype.map { m, p -> p }.first()
    ch_la_single    = ch_local_ancestry.map { m, f -> f }.first().ifEmpty([])

    // One record per trait: ancestry list + matching sumstats list (for PRS-CSx / SDPR)
    ch_by_trait = ch_sumstats
        .filter { meta, ss -> meta.ancestry && meta.ancestry != 'META' }
        .map { meta, ss -> [[trait: meta.trait, binary: meta.binary ?: false, survival: meta.survival ?: false,
                             time_col: meta.time_col, event_col: meta.event_col], meta.ancestry, ss] }
        .groupTuple(by: 0)
        .map { meta, ancs, files -> [meta, ancs, files] }

    ch_method_input = ch_by_trait
        .combine(ch_geno_single)
        .combine(ch_la_single)
        .combine(ch_pheno_single)
        .map { meta, ancs, files, gmeta, bed, bim, fam, la, pheno -> [meta, ancs, files, bed, bim, fam, la, pheno] }

    // ------------------------------------------------------------------
    // All methods in parallel (each method = one task per trait)
    // ------------------------------------------------------------------
    PRS_METHOD(ch_method_input, prs_methods, ld_ref)
    ch_versions = ch_versions.mix(PRS_METHOD.out.versions)

    ch_scores = PRS_METHOD.out.results
        .map { meta, method, files -> [meta + [method: method], files] }

    // ------------------------------------------------------------------
    // Local-ancestry partial scores (Tractor dosages x Tractor-GENESIS betas)
    // ------------------------------------------------------------------
    // Weights: PRS-CSx shrunk per-ancestry posterior weights when that method
    // ran (preferred), else Tractor-GENESIS betas thresholded on P_JOINT
    ch_prscsx_weights = ch_scores
        .filter { meta, files -> meta.method == 'prs_csx' }
        .map { meta, files -> [[trait: meta.trait], (files instanceof List ? files : [files]).findAll { it.name =~ /weights\.tsv/ }] }

    ch_la_input = ch_tractor_sumstats
        .map { meta, ss -> [[trait: meta.trait], meta.tractor_pops ?: params.tractor_lat_pops, ss] }
        .join(ch_prscsx_weights, remainder: true)
        .map { meta, ancs, ss, w -> [meta, ancs, ss, w ?: []] }
        .combine(ch_tractor_dosages.map { m, f -> f }.collect().map { [it] })
        .combine(ch_pheno_single)
        .map { meta, ancs, ss, w, dos, pheno -> [meta, ancs.tokenize(','), ss, w, dos, pheno] }

    PRS_LA_PARTIAL(ch_la_input)
    ch_versions = ch_versions.mix(PRS_LA_PARTIAL.out.versions)

    // ------------------------------------------------------------------
    // Validation: overall + per stratum, then method comparison
    // ------------------------------------------------------------------
    ch_validation = Channel.empty()
    ch_comparison = Channel.empty()
    ch_best       = Channel.empty()
    if (run_validation) {
        ch_val_input = ch_scores
            .map { meta, files -> [[trait: meta.trait, binary: meta.binary, survival: meta.survival,
                                    time_col: meta.time_col, event_col: meta.event_col], files] }
            .groupTuple(by: 0)
            .map { meta, files -> [meta, files.flatten()] }
            .combine(ch_pheno_single)

        PRS_VALIDATE_STRATIFIED(ch_val_input, ancestry_col)
        ch_validation = PRS_VALIDATE_STRATIFIED.out.validation
        ch_versions = ch_versions.mix(PRS_VALIDATE_STRATIFIED.out.versions)

        PRS_COMBINE_METHODS(ch_validation)
        ch_comparison = PRS_COMBINE_METHODS.out.comparison
        ch_best = PRS_COMBINE_METHODS.out.best_prs
        ch_versions = ch_versions.mix(PRS_COMBINE_METHODS.out.versions)
    }

    emit:
    prs_scores         = ch_scores                        // channel: [ meta(trait, method), files ]
    la_partial         = PRS_LA_PARTIAL.out.scores        // channel: [ meta(trait), la_partial.scores.tsv ]
    validation_results = ch_validation                    // channel: [ meta(trait), validation.tsv ]
    method_comparison  = ch_comparison                    // channel: [ meta(trait), prs_comparison.tsv ]
    best_method        = ch_best                          // channel: [ meta(trait), best_prs.tsv ]
    versions           = ch_versions
}
