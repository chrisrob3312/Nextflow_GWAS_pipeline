/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VISUALIZATION SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Lambda GC, Manhattan and QQ per GWAS, regional plots at significant loci,
    Miami plots of each ancestry stratum against the reference stratum.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LAMBDA_GC      } from '../../modules/local/lambda_gc'
include { MANHATTAN_PLOT } from '../../modules/local/manhattan_plot'
include { QQ_PLOT        } from '../../modules/local/qq_plot'
include { REGIONAL_PLOT  } from '../../modules/local/regional_plot'
include { MIAMI_PLOT     } from '../../modules/local/regional_plots_enhanced'

workflow VISUALIZATION {
    take:
    ch_gwas_results       // channel: [ meta, sumstats ]
    ch_significant        // channel: [ meta, sig_variants ]
    gene_annotation       // file or []

    main:
    ch_versions = Channel.empty()
    ch_plots = Channel.empty()

    LAMBDA_GC(ch_gwas_results)
    ch_versions = ch_versions.mix(LAMBDA_GC.out.versions)

    ch_with_lambda = ch_gwas_results.join(LAMBDA_GC.out.lambda)

    MANHATTAN_PLOT(ch_with_lambda)
    ch_plots = ch_plots.mix(MANHATTAN_PLOT.out.plots.map { m, f -> f })
    ch_versions = ch_versions.mix(MANHATTAN_PLOT.out.versions)

    QQ_PLOT(ch_with_lambda)
    ch_plots = ch_plots.mix(QQ_PLOT.out.plots.map { m, f -> f })
    ch_versions = ch_versions.mix(QQ_PLOT.out.versions)

    if (gene_annotation) {
        REGIONAL_PLOT(ch_significant, gene_annotation)
        ch_plots = ch_plots.mix(REGIONAL_PLOT.out.plots.map { m, f -> f })
        ch_versions = ch_versions.mix(REGIONAL_PLOT.out.versions)
    }

    // Miami: each stratum vs the reference stratum (EUR if present, else first)
    ch_for_miami = ch_gwas_results
        .filter { meta, ss -> meta.ancestry }
        .map { meta, ss -> [[trait: meta.trait], meta.ancestry, ss] }
        .groupTuple(by: 0)
        .filter { t, ancs, files -> ancs.size() >= 2 }
        .flatMap { t, ancs, files ->
            def ref = ancs.indexOf('EUR') >= 0 ? ancs.indexOf('EUR') : 0
            def pairs = []
            for (int i = 0; i < ancs.size(); i++) {
                if (i != ref) pairs << [t + [anc1: ancs[ref], anc2: ancs[i], id: "${t.trait}.${ancs[ref]}_vs_${ancs[i]}"], files[ref], files[i]]
            }
            pairs
        }

    MIAMI_PLOT(ch_for_miami, 'ancestry', 'top', 'bottom')   // labels replaced by meta.anc1 / meta.anc2 in the process
    ch_plots = ch_plots.mix(MIAMI_PLOT.out.png.map { m, f -> f })
    ch_versions = ch_versions.mix(MIAMI_PLOT.out.versions)

    emit:
    plots     = ch_plots                 // channel: plot files
    lambda_gc = LAMBDA_GC.out.lambda     // channel: [ meta, lambda ]
    versions  = ch_versions
}
