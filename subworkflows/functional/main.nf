/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONAL_ANNOT SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    'magma'   : SNP -> gene annotation, gene-level test, gene-set test (optional)
    'annovar' : variant annotation (ANNOVAR), RegulomeDB / HaploReg optional
    FUMA, LAVA, FLAMES and cell-type enrichment are not bundled: listing them
    in annotation_tools raises an error.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAGMA_ANNOTATE      } from '../../modules/local/magma_annotate'
include { FUNCTIONAL_SUMMARY  } from '../../modules/local/magma_annotate'
include { MAGMA_GENE          } from '../../modules/local/magma_gene'
include { ANNOVAR_ANNOTATE    } from '../../modules/local/enhanced_annotation'

workflow FUNCTIONAL_ANNOT {
    take:
    ch_gwas_results       // channel: [ meta, sumstats ]
    annotation_tools      // list
    magma_geneset         // file or [] (gene-set analysis not bundled as a process yet; logged)
    magma_annotation      // file: MAGMA gene.loc
    lava_partition        // file or [] (LAVA not bundled)
    cell_type_analysis    // boolean (not bundled)
    cell_type_datasets    // file or []

    main:
    ch_versions = Channel.empty()
    ch_results = Channel.empty()

    def supported = ['magma', 'annovar']
    def unknown = annotation_tools.findAll { !(it in supported) }
    if (unknown) error "annotation_tools contains ${unknown}: only ${supported} are bundled."
    if (cell_type_analysis) log.warn "cell_type_analysis = true: cell-type enrichment is not bundled; skipped."
    if (lava_partition) log.warn "lava_partition given but LAVA is not bundled; skipped."

    // ------------------------------------------------------------------
    // MAGMA
    // ------------------------------------------------------------------
    if ('magma' in annotation_tools) {
        if (!magma_annotation) error "annotation_tools includes 'magma' but magma_annotation (gene.loc) is not set."
        if (!params.magma_reference) error "annotation_tools includes 'magma' but magma_reference (PLINK prefix dir) is not set."

        def magma_window = params.magma_window ?: '35,10'
        MAGMA_ANNOTATE(ch_gwas_results, magma_annotation, magma_window)
        ch_versions = ch_versions.mix(MAGMA_ANNOTATE.out.versions)

        MAGMA_GENE(MAGMA_ANNOTATE.out.annotated, file(params.magma_reference), params.gwas_n)
        ch_results = ch_results.mix(MAGMA_GENE.out.genes.map { meta, f -> [meta + [tool: 'magma_gene'], f] })
        ch_versions = ch_versions.mix(MAGMA_GENE.out.versions)

        if (magma_geneset) log.warn "magma_geneset given: gene-set analysis process is not bundled yet; run 'magma --gene-results <prefix>.genes.raw --set-annot <geneset>' on the .genes.raw outputs."
    }

    // ------------------------------------------------------------------
    // ANNOVAR
    // ------------------------------------------------------------------
    if ('annovar' in annotation_tools) {
        if (!params.annovar_db) error "annotation_tools includes 'annovar' but annovar_db is not set."
        def build = (params.genome_build ?: 'GRCh38') == 'GRCh37' ? 'hg19' : 'hg38'
        ANNOVAR_ANNOTATE(ch_gwas_results, file(params.annovar_db), build, ['refGene', 'gnomad', 'clinvar', 'cadd'])
        ch_results = ch_results.mix(ANNOVAR_ANNOTATE.out.annotated.map { meta, f -> [meta + [tool: 'annovar'], f] })
        ch_versions = ch_versions.mix(ANNOVAR_ANNOTATE.out.versions)
    }

    FUNCTIONAL_SUMMARY(ch_results.map { m, f -> f }.collect().ifEmpty([]))
    ch_versions = ch_versions.mix(FUNCTIONAL_SUMMARY.out.versions)

    emit:
    annotation_results = ch_results                        // channel: [ meta(tool), results ]
    summary            = FUNCTIONAL_SUMMARY.out.summary    // channel: [ functional_summary.tsv ]
    versions           = ch_versions
}
