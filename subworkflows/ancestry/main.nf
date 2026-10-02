/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANCESTRY_INFERENCE SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Two ways in:
      run_ancestry_inference = true  -> GRAF-ANC is run on the QC'd cohort
      run_ancestry_inference = false -> --ancestry_calls_file (sample_id, ancestry)
    Either way the calls are normalised (GRAF codes -> labels) and the cohort is
    split into one bed/bim/fam per ancestry group with N >= min_ancestry_n.
    Local ancestry is NOT inferred here: supply RFMix MSP files with
    --local_ancestry_files (glob). They are attached to every stratum so the
    GWAS subworkflow can hand them to Tractor.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GRAF_ANCESTRY          } from '../../modules/local/graf_ancestry'
include { ANCESTRY_CALLS_PREPARE } from '../../modules/local/ancestry_stratify'
include { ANCESTRY_STRATIFY      } from '../../modules/local/ancestry_stratify'

workflow ANCESTRY_INFERENCE {
    take:
    ch_genotypes            // channel: [ meta, bed, bim, fam ]  QC'd cohort
    run_inference           // boolean: run GRAF-ANC here
    ancestry_method         // value: 'graf-anc' | 'graf-pop' | 'graf' (others not bundled)
    graf_reference          // file or []
    ancestry_calls_file     // file or []: pre-computed calls (sample_id, ancestry)
    min_ancestry_n          // integer: minimum samples per stratum
    run_local_ancestry      // boolean: Tractor needs local ancestry
    local_ancestry_files    // string glob or null: RFMix MSP files

    main:
    ch_versions = Channel.empty()

    // ------------------------------------------------------------------
    // Ancestry calls
    // ------------------------------------------------------------------
    if (run_inference) {
        def m = (ancestry_method ?: 'graf-anc').toLowerCase()
        if (!(m in ['graf', 'graf-anc', 'graf-pop'])) {
            error "ancestry_method '${ancestry_method}' is not bundled with this pipeline. Use 'graf-anc' or supply --ancestry_calls_file."
        }
        GRAF_ANCESTRY(ch_genotypes, graf_reference)
        ch_raw_calls = GRAF_ANCESTRY.out.ancestry_calls
        ch_versions = ch_versions.mix(GRAF_ANCESTRY.out.versions)
    } else {
        if (!ancestry_calls_file) {
            error "run_ancestry_inference = false requires --ancestry_calls_file (columns: sample_id, ancestry)."
        }
        ch_raw_calls = ch_genotypes.map { meta, bed, bim, fam -> [meta, file(ancestry_calls_file)] }
    }

    ANCESTRY_CALLS_PREPARE(ch_raw_calls)
    ch_ancestry_calls = ANCESTRY_CALLS_PREPARE.out.ancestry_calls
    ch_versions = ch_versions.mix(ANCESTRY_CALLS_PREPARE.out.versions)

    // ------------------------------------------------------------------
    // Stratify: one [meta + ancestry, bed, bim, fam] per group
    // ------------------------------------------------------------------
    ANCESTRY_STRATIFY(ch_genotypes.join(ch_ancestry_calls), min_ancestry_n)
    ch_versions = ch_versions.mix(ANCESTRY_STRATIFY.out.versions)

    ch_stratified = ANCESTRY_STRATIFY.out.stratified_genotypes
        .flatMap { meta, beds, bims, fams ->
            def bed_list = beds instanceof List ? beds : [beds]
            def bim_list = bims instanceof List ? bims : [bims]
            def fam_list = fams instanceof List ? fams : [fams]
            bed_list.collect { bed ->
                def anc = bed.name.replaceFirst(/\.bed$/, '').tokenize('.').last()
                def bim = bim_list.find { it.name == bed.name.replaceFirst(/\.bed$/, '.bim') }
                def fam = fam_list.find { it.name == bed.name.replaceFirst(/\.bed$/, '.fam') }
                [meta + [id: "${meta.id}.${anc}", ancestry: anc], bed, bim, fam]
            }
        }

    // ------------------------------------------------------------------
    // Local ancestry (pre-computed MSP files attached to every stratum)
    // ------------------------------------------------------------------
    ch_local_ancestry = Channel.empty()
    if (run_local_ancestry) {
        if (!local_ancestry_files) {
            error "local_ancestry = true (needed for Tractor) but no --local_ancestry_files given. RFMix / G-Nomix are not bundled; run them first and pass the MSP files."
        }
        ch_msp = Channel.fromPath(local_ancestry_files, checkIfExists: true).collect()
        ch_local_ancestry = ch_stratified
            .combine(ch_msp.map { [it] })
            .map { meta, bed, bim, fam, msp -> [meta, msp] }
    }

    emit:
    ancestry_calls       = ch_ancestry_calls                     // channel: [ meta, sample_id/ancestry file ]
    ancestry_counts      = ANCESTRY_CALLS_PREPARE.out.counts     // channel: [ meta, counts ]
    stratified_genotypes = ch_stratified                         // channel: [ meta(ancestry), bed, bim, fam ]
    local_ancestry       = ch_local_ancestry                     // channel: [ meta(ancestry), [msp files] ]
    versions             = ch_versions
}
