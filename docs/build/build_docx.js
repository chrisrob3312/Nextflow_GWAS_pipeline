// Build docs/GWAS_pipeline_workflow.docx from the rendered diagrams.
// Usage (from docs/build): node render_diagrams.js && node build_docx.js ../GWAS_pipeline_workflow.docx
const fs = require('fs');
const {
    Document, Packer, Paragraph, TextRun, ImageRun, HeadingLevel, AlignmentType,
    Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle, PageOrientation, LevelFormat
} = require('docx');

const sizes = JSON.parse(fs.readFileSync('diagrams.json', 'utf8'));
const OUT = process.argv[2] || '../GWAS_pipeline_workflow.docx';

const LETTER = { width: 12240, height: 15840 };
const MARGIN = 1440;
const PORTRAIT_W_PX = 624, PORTRAIT_H_PX = 780;
const LANDSCAPE_W_PX = 864, LANDSCAPE_H_PX = 540;

function img(name, maxW, maxH) {
    const { w, h } = sizes[name];
    const s = Math.min(maxW / w, maxH / h, 1);
    return new Paragraph({
        alignment: AlignmentType.CENTER, spacing: { before: 120, after: 120 },
        children: [new ImageRun({ type: 'png', data: fs.readFileSync(`${name}.png`),
            transformation: { width: Math.round(w * s), height: Math.round(h * s) } })]
    });
}
const H1 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_1, spacing: { before: 240, after: 120 }, children: [new TextRun(t)] });
const H2 = (t) => new Paragraph({ heading: HeadingLevel.HEADING_2, spacing: { before: 200, after: 80 }, children: [new TextRun(t)] });
const P  = (t, opts = {}) => new Paragraph({ spacing: { after: 120 }, children: [new TextRun({ text: t, ...opts })] });
const bullet = (t) => new Paragraph({ numbering: { reference: 'bullets', level: 0 }, spacing: { after: 60 }, children: [new TextRun(t)] });
const caption = (t) => new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 200 }, children: [new TextRun({ text: t, italics: true, size: 18, color: '555555' })] });

const thin = { style: BorderStyle.SINGLE, size: 4, color: 'BBBBBB' };
const borders = { top: thin, bottom: thin, left: thin, right: thin };
function table(headers, rows, widths) {
    const total = widths.reduce((a, b) => a + b, 0);
    const cell = (text, w, head) => new TableCell({
        width: { size: w, type: WidthType.DXA }, borders,
        shading: head ? { type: ShadingType.CLEAR, fill: 'E7EEF5', color: 'auto' } : undefined,
        margins: { top: 60, bottom: 60, left: 90, right: 90 },
        children: [new Paragraph({ children: [new TextRun({ text, bold: head, size: 17, font: 'Calibri' })] })]
    });
    return new Table({
        width: { size: total, type: WidthType.DXA }, columnWidths: widths,
        rows: [new TableRow({ tableHeader: true, children: headers.map((h, i) => cell(h, widths[i], true)) }),
               ...rows.map(r => new TableRow({ children: r.map((c, i) => cell(c, widths[i], false)) }))]
    });
}

const portrait  = { page: { size: LETTER, margin: { top: MARGIN, bottom: MARGIN, left: MARGIN, right: MARGIN } } };
const landscape = { page: { size: { ...LETTER, orientation: PageOrientation.LANDSCAPE }, margin: { top: MARGIN, bottom: MARGIN, left: MARGIN, right: MARGIN } } };

const inputRows = [
    ['Samplesheet', 'input', 'CSV with header. Required column sample_id (one row per genotype set), plus ONE of bed,bim,fam / pgen,pvar,psam / bgen,sample[,bgi] / vcf[,vcf_index]; optional cohort', 'Validated by bin/check_samplesheet.py; non-PLINK1 sets are converted in QC'],
    ['Phenotype file', 'phenotype_file', 'TSV with header. FID and IID must be the first two columns (REGENIE, SAIGE, SPA-Cox require it; the R scripts accept IID or sample_id). IID must match the genotype sample IDs', 'Nothing is inferred from the file'],
    ['Traits', 'phenotype_cols', 'comma list of column names', 'Traits are exactly these columns. binary_traits (0/1) and survival_traits (needs <trait>_time, <trait>_status) pick the model; anything else is quantitative'],
    ['Covariates', 'covariate_cols', 'comma list of column names', 'Covariates are exactly these columns. List PC1..PCn here; with run_pcair = true they are replaced by PC-AiR PCs before any model runs'],
    ['Ancestry labels', 'gxg_ancestry_col (default GRAF_ANC)', 'a column of the phenotype file', 'GRAF-ANC codes or labels; used by GxG and PRS validation to stratify'],
    ['Ancestry calls', 'ancestry_calls_file (when run_ancestry_inference = false)', 'TSV or CSV: sample_id, ancestry (codes or labels)', 'Normalised by ANCESTRY_CALLS_PREPARE; defines the per-ancestry strata (N >= min_ancestry_n)'],
    ['Local ancestry', 'local_ancestry_files', 'glob of RFMix .msp.tsv files', 'Required for Tractor; RFMix / G-Nomix are not bundled'],
    ['Kinship / GRM', 'none needed', 'produced by GENESIS_PCAIR_PCRELATE', 'kinship_matrix is read only when run_pcair = false'],
    ['LD reference', 'ld_reference_dir', 'directory with <ancestry>.bed/bim/fam (cohort-specific from calculate_cohort_ld.R, or a panel)', 'SuSiE-RSS, PRS, colocalization'],
    ['QTL datasets', 'custom_qtl_datasets / eqtl_datasets', 'CSV: name,type,path[,ancestry,tissue]', 'plus the curated megaset via download_qtl_datasets.R']
];

const moduleRows = [
    ['INPUT_CHECK', 'samplesheet CSV', 'validated genotype channel', 'QC'],
    ['QC_WORKFLOW', 'bed / bim / fam', 'QC\'d bed / bim / fam, QC report', 'ANCESTRY, PC-AiR, GxG (full cohort)'],
    ['GENESIS_PCAIR_PCRELATE (bin/genesis_pcair_pcrelate.R)', 'QC\'d full cohort, phenotype', 'PC-AiR PCs; PC-Relate GRM (.pcrelate.grm.rds); kinship; unrelated set; phenotype with PC1..PCn replaced by PC-AiR PCs (FID, IID order kept)', 'GWAS null models, Tractor-GENESIS, GxG, PRS validation: every model uses the same GRM and PCs'],
    ['ANCESTRY_INFERENCE (GRAF-ANC or ancestry_calls_file -> ANCESTRY_CALLS_PREPARE -> ANCESTRY_STRATIFY)', 'QC\'d genotypes; calls file; RFMix MSP files', 'normalised calls (sample_id, ancestry); one bed/bim/fam per ancestry with N >= min_ancestry_n; MSP files attached to every stratum', 'GWAS, Tractor, PRS'],
    ['GWAS_WORKFLOW (standard)', 'stratified genotypes x trait, phenotype, covariates, GRM', '<id>.<ancestry>.<trait>.sumstats.gz; filtered sumstats; significant variants', 'META, FM, H2, VIS, GxG'],
    ['PLINK_TO_VCF -> TRACTOR_EXTRACT_TRACTS', 'AAC and merged LATINO genotypes + MSP; meta.tractor_pops', '<id>.ancdose.<k>.tsv.gz and <id>.hapcount.<k>.tsv.gz per ancestry k', 'TRACTOR_GENESIS, PRS (LA partial), GxG'],
    ['TRACTOR_GENESIS (bin/tractor_genesis_adapter.R)', 'dosages + hapcounts, phenotype, GRM, covariates; meta.model', '<prefix>.tractor_genesis.tsv.gz (sorted by P_JOINT), .genomic_order.tsv.gz, .top_hits.tsv, .significant.tsv, .heterogeneous.tsv, .<anc>_specific.tsv, .null_model.rds', 'META, GxG, PRS'],
    ['META_ANALYSIS (MR_MEGA or META_STRATA = bin/meta_analyze_strata.R)', 'sumstats grouped by trait across strata', 'MR-MEGA results + heterogeneity, or POOLED.sumstats.gz (IVW per ancestry background)', 'FM, H2, PRS, FUNC, GxG'],
    ['FINE_MAPPING (EXTRACT_REGIONS -> SUSIE_RSS / POLYFUN_SUSIE; MG_FLASH_FM / SUSIE_ME; CREDIBLE_SETS)', 'filtered sumstats, ld_reference_dir (<ancestry>.bed)', 'per-region results, PIPs, merged credible sets per trait', 'COLOC'],
    ['COLOCALIZATION (COLOC_RUN = bin/run_colocalization.R)', 'every ancestry-stratified GWAS (standard + Tractor-GENESIS) and the meta / POOLED GWAS; QTL megaset pooled by type; optional cohort LD', 'per GWAS: .coloc_susie.tsv, .hyprcoloc.tsv, .opera.tsv; per trait: coloc_combined.tsv, coloc_consensus.tsv, shared / divergent / meta-only tables, PP4 heatmap', 'REPORTING'],
    ['HERITABILITY (LDSC_MUNGE once -> LDSC_H2, LDSC_RG -> H2_COMPARISON)', 'sumstats per stratum + POOLED, LD scores', 'h2 per stratum and POOLED, cross-ancestry rg per pair, comparison table', 'REPORTING'],
    ['PRS_WORKFLOW (PRS_METHOD = bin/calculate_prs_admixed.R)', 'ancestry-specific sumstats per trait, LD reference, full-cohort genotypes, local ancestry, phenotype with PC-AiR PCs; Tractor dosages + PRS-CSx weights (or Tractor betas) for the partial scores', 'per method: weights and .sscore; <trait>.la_partial.scores.tsv; <trait>.validation.tsv (OVERALL + per stratum); prs_comparison.tsv, best_prs.tsv', 'REPORTING'],
    ['FUNCTIONAL_ANNOT (MAGMA_ANNOTATE -> MAGMA_GENE; ANNOVAR_ANNOTATE; FUNCTIONAL_SUMMARY)', 'meta or filtered sumstats, gene.loc, MAGMA LD reference, ANNOVAR db', 'gene-level results, variant annotations, summary', 'REPORTING'],
    ['GXG_INTERACTION (bin/run_gxg_interaction.R)', 'all sumstats for the trait, known loci, custom list, full-cohort genotypes, phenotype, optional Tractor dosages', 'gxg.<stratum>.tsv, gxg.all.tsv, gxg.significant.tsv, gxg.hits_pruned.tsv, gxg.ancestry_heterogeneity.tsv, gxg.local_ancestry.tsv, gxg_summary.tsv, gxg_top_interactions.tsv', 'REPORTING'],
    ['REPORTING', 'everything above', 'HTML report', 'user']
];

const customRows = [
    ['Traits', 'phenotype_cols', 'TRAITS_LIST', 'any phenotype columns'],
    ['Binary traits', 'binary_traits', 'everything not in SURVIVAL_TRAITS', '0 / 1'],
    ['Time-to-event traits', 'survival_traits (+ survival_time_cols, survival_event_cols)', 'SURVIVAL_TRAITS', 'default columns <trait>_time, <trait>_status'],
    ['Covariates', 'covariate_cols, gxg_covariates', 'COVARIATES', 'any phenotype columns'],
    ['Ancestry column', 'gxg_ancestry_col', 'ANCESTRY_COL', 'GRAF-ANC codes or labels'],
    ['Strata rules', 'min_ancestry_n, min_stratum_n, bin/ancestry_config.R', 'STRATA_LIST, bin/ancestry_config.R', 'N < 30 pools to OTHER'],
    ['Known loci', 'gxg_known_loci, gxg_use_default_loci', 'KNOWN_LOCI, GXG_USE_DEFAULT_LOCI', 'default table is a starting point; verify positions'],
    ['Custom variants', 'gxg_custom_variants', 'GXG_CUSTOM_VARIANTS', 'GRCh38, one per line'],
    ['LD-correlated hits', 'gxg_prune_mode, gxg_cond_p_threshold, gxg_max_pair_r2, gxg_min_distance_kb', 'GXG_PRUNE_MODE, GXG_COND_P', 'conditional (default): tested conditional on kept hits, kept if cond P < 1e-4; GWAS > custom > known'],
    ['Tractor', 'tractor_aac_pops, tractor_lat_pops, tractor_ref_ancestry, tractor_mac_min', 'ANCESTRIES case in submit_gwas_array.sh', ''],
    ['Optional steps', 'meta_method (mr-mega | ivw), within_ancestry_fm (susie | polyfun-susie), h2_method (ldsc), annotation_tools (magma, annovar)', '', 'other method names raise an explicit error']
];

const doc = new Document({
    creator: 'Nextflow GWAS pipeline',
    title: 'Multi-ancestry leukemia GWAS pipeline: workflow',
    styles: {
        default: { document: { run: { font: 'Calibri', size: 22 } } },
        paragraphStyles: [
            { id: 'Heading1', name: 'Heading 1', basedOn: 'Normal', next: 'Normal', quickFormat: true,
              run: { size: 32, bold: true, color: '1F3864', font: 'Calibri' }, paragraph: { spacing: { before: 240, after: 120 }, outlineLevel: 0 } },
            { id: 'Heading2', name: 'Heading 2', basedOn: 'Normal', next: 'Normal', quickFormat: true,
              run: { size: 26, bold: true, color: '2E5C8A', font: 'Calibri' }, paragraph: { spacing: { before: 200, after: 80 }, outlineLevel: 1 } }
        ]
    },
    numbering: { config: [{ reference: 'bullets', levels: [{ level: 0, format: LevelFormat.BULLET, text: '•', alignment: AlignmentType.LEFT, style: { paragraph: { indent: { left: 720, hanging: 360 } } } }] }] },
    sections: [
        { properties: portrait, children: [
            new Paragraph({ heading: HeadingLevel.TITLE, spacing: { after: 80 }, children: [new TextRun({ text: 'Multi-ancestry leukemia GWAS pipeline', bold: true, size: 40, color: '1F3864' })] }),
            P('Workflow: what happens to the data, module by module', { size: 26, color: '555555' }),
            P('Every box is a Nextflow subworkflow or process and every arrow is a file handed from one module to the next. The same steps exist as SLURM array scripts under slurm/, which call the same bin/ scripts (section 6). Ancestry strata are EUR, AAC, LAT1 and LAT2, with EAS and SAS analysed separately when they reach 30 samples and pooled into OTHER otherwise.'),
            H1('0. Inputs: what the pipeline expects at the start'),
            table(['Input', 'Param', 'Format', 'Notes'], inputRows, [1500, 2000, 3260, 2600]),
            P('Example phenotype file header: FID  IID  OS_time  OS_status  DFS_time  DFS_status  MRD  age  sex  PC1 ... PC5  GRAF_ANC, with phenotype_cols = OS,DFS,MRD; survival_traits = OS,DFS; binary_traits = MRD; covariate_cols = age,sex,PC1,PC2,PC3,PC4,PC5.')
        ] },
        { properties: portrait, children: [
            H1('1. End-to-end overview'),
            H2('Part A: genotypes to summary statistics'),
            img('overview_a', PORTRAIT_W_PX, PORTRAIT_H_PX - 140),
            caption('Every trait goes through Tractor-GENESIS (green) as well as the standard per-stratum GWAS, so MRD, DFS and OS are comparable. PC-AiR PCs and the PC-Relate GRM feed every model.')
        ] },
        { properties: landscape, children: [
            H2('Part B: summary statistics to results'),
            img('overview_b', LANDSCAPE_W_PX, LANDSCAPE_H_PX),
            caption('GxG interaction testing is the last analysis step and reads all summary statistics produced upstream.')
        ] },
        { properties: portrait, children: [
            H1('2. Module inputs and outputs'),
            table(['Module', 'Reads', 'Writes', 'Feeds'], moduleRows, [2100, 2500, 3160, 1600])
        ] },
        { properties: landscape, children: [
            H1('3. Tractor-GENESIS: one engine for every trait'),
            P('Every trait goes through the same conditional score test. The trait\'s model, set in main.nf from binary_traits and survival_traits, only changes where the projection comes from: the GENESIS mixed-model null for binary and quantitative traits, or a Cox null with a kinship frailty and martingale residuals for OS and DFS. DFS (relapse or death) and OS are single-event outcomes, so no competing-risk model is needed.'),
            img('tractor', LANDSCAPE_W_PX, LANDSCAPE_H_PX - 120),
            caption('D = [LA terms | Dose terms]. The dose block is conditioned on the local-ancestry block, giving the joint k-df test, per-ancestry effects and a covariance-aware heterogeneity test.')
        ] },
        { properties: portrait, children: [
            H1('4. GxG interaction (end of pipeline)'),
            img('gxg', PORTRAIT_W_PX, PORTRAIT_H_PX - 40),
            caption('Models: logistic for binary traits, linear for quantitative, Cox for survival. Each pair reports Wald and LRT p-values with Bonferroni over pairs and BH FDR.'),
            H2('How hits in LD are resolved: conditional testing'),
            bullet('Hits are ordered GWAS hits by p-value first, then the custom list, then the default known loci.'),
            bullet('Walking down that order, a hit whose r2 with an already-kept hit on the same chromosome exceeds gxg_max_pair_r2 (0.2) is not dropped. It is tested conditional on all its kept LD partners in the pooled cohort (y ~ hit + partners + covariates + stratum, likelihood-ratio test).'),
            bullet('If its conditional p-value is below gxg_cond_p_threshold (1e-4) it is kept as an independent signal; otherwise it is recorded as explained by its proxy. hits_pruned.tsv lists status, proxy, r2 and conditional p for every hit.'),
            bullet('This is an individual-level stepwise conditional and joint analysis, the analysis GCTA-COJO approximates from summary statistics and reference LD. It runs on the hit list, not genome-wide.'),
            bullet('Proximity alone never drops a variant: nearby independent signals both stay, they are just not paired with each other. gxg_prune_mode = variant keeps only the strongest per cluster; pair keeps everything and skips LD pairs.')
        ] },
        { properties: portrait, children: [
            H1('5. Traits, covariates, strata: what is customisable'),
            table(['Setting', 'Nextflow', 'SLURM', 'Notes'], customRows, [1700, 3000, 2560, 2100]),
            H1('6. SLURM orchestration (slurm/submit_full_pipeline.sh)'),
            img('slurm', PORTRAIT_W_PX, 400),
            caption('Step 6 needs only the combined summary statistics, so it runs in parallel with PRS and colocalization.'),
            H1('7. Dual structure: both paths call the same scripts'),
            P('Every analysis step has one implementation in bin/ that both the Nextflow modules and the SLURM scripts call, so the two paths cannot drift apart.'),
            table(['Step', 'Script', 'Nextflow process', 'SLURM'], [
                ['PC-AiR + PC-Relate', 'genesis_pcair_pcrelate.R', 'GENESIS_PCAIR_PCRELATE', 'submit_pcair.sh (step 0)'],
                ['Tractor-GENESIS GWAS', 'tractor_genesis_adapter.R', 'TRACTOR_GENESIS', 'submit_gwas_array.sh'],
                ['Strata meta-analysis', 'meta_analyze_strata.R', 'META_STRATA (meta_method = ivw)', 'step 2'],
                ['PRS: 6 methods, LA partial, stratified validation', 'calculate_prs_admixed.R', 'PRS_METHOD, PRS_LA_PARTIAL, PRS_VALIDATE_STRATIFIED', 'submit_prs_array.sh'],
                ['Colocalization: coloc.susie, HyPrColoc, OPERA', 'run_colocalization.R', 'COLOC_RUN, COLOC_COMBINE', 'step 4'],
                ['GxG interaction', 'run_gxg_interaction.R', 'GXG_SELECT_HITS, GXG_TEST, GXG_ANCESTRY_HETEROGENEITY, GXG_COMBINE', 'submit_gxg_array.sh']
            ], [2400, 2200, 2760, 2000]),
            H1('8. What comes out per ancestry'),
            P('Tractor-GENESIS reports, for every variant and trait: P_JOINT (any effect) and, for each ancestral background, BETA, SE, Z, P and MAC, plus P_HET and I2 for whether the effects differ across backgrounds. The per-stratum SLURM runs add stratum-level results (EUR, AAC, LAT1, LAT2, and so on) on top, and meta_analyze_strata.R pools each background across strata. PRS adds the local-ancestry partial scores. GxG reports every interaction pooled, per stratum, and for ancestry modification.')
        ] }
    ]
});

Packer.toBuffer(doc).then(buf => { fs.writeFileSync(OUT, buf); console.log('wrote', OUT, buf.length, 'bytes'); });
