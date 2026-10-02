// Render the workflow diagrams to PNG with Mermaid in headless Chromium.
// Usage (from docs/build):  npm install && node render_diagrams.js
// Chromium: Playwright's bundle (PLAYWRIGHT_BROWSERS_PATH) or CHROMIUM_PATH env.
const { chromium } = require('playwright-core');
const fs = require('fs');
const mermaidPath = require.resolve('mermaid/dist/mermaid.min.js');
const exe = process.env.CHROMIUM_PATH || '/opt/pw-browsers/chromium';

const CLASSES = `
    classDef in fill:#fff4d6,stroke:#b58a00,color:#222;
    classDef core fill:#e8f0fb,stroke:#3b6ea5,color:#222;
    classDef tr fill:#e6f5ea,stroke:#2e7d4f,color:#222;
    classDef end_ fill:#f3e8f7,stroke:#7a3e9d,color:#222;
`;

const diagrams = {
overview_a: { wrap: 320, code: `
flowchart TD
    IN["INPUTS<br/>samplesheet: PLINK bed / bim / fam<br/>phenotype TSV: FID, IID, traits, covariates, GRAF_ANC<br/>ancestry calls &nbsp;|&nbsp; local ancestry MSP (RFMix)"]
    QC["QC_WORKFLOW<br/>MAF, HWE, missingness filters"]
    PCA["GENESIS_PCAIR_PCRELATE<br/>PC-AiR ancestry PCs + PC-Relate GRM<br/>used by EVERY model (GWAS, Tractor, GxG, PRS)"]
    ANC["ANCESTRY_INFERENCE<br/>GRAF-ANC or pre-computed calls -> strata: EUR, AAC, LAT1, LAT2,<br/>EAS / SAS if N >= 30, else OTHER"]
    GWAS["GWAS_WORKFLOW (standard)<br/>per stratum x trait<br/>REGENIE / SAIGE / GENESIS / SPA-Cox"]
    TR["TRACTOR-GENESIS<br/>AAC (EUR, AFR) + LATINO (EUR, AFR, AMR)<br/>ALL traits: MRD binary, OS and DFS time-to-event"]
    META["META_ANALYSIS<br/>MR-MEGA or IVW across strata -> POOLED sumstats"]
    VIS["VISUALIZATION<br/>lambda GC, Manhattan, QQ, regional, Miami"]
    NEXT(["-> Part B: summary statistics to results"])
    IN --> QC --> ANC
    QC --> PCA
    ANC -->|stratified bed / bim / fam| GWAS
    ANC -->|stratified genotypes + MSP| TR
    PCA -->|GRM + PC covariates| GWAS
    PCA -->|GRM + PC covariates| TR
    IN -->|phenotype| GWAS
    IN -->|phenotype| TR
    GWAS -->|sumstats per stratum x trait| META
    TR -->|tractor_genesis.tsv.gz per group x trait| META
    GWAS --> VIS
    META --> NEXT
    GWAS -.->|filtered sumstats, all sumstats| NEXT
    ${CLASSES}
    class IN in; class QC,ANC,GWAS,META,VIS core; class TR,PCA tr; class NEXT end_;
` },

overview_b: { wrap: 300, font: 20, code: `
flowchart TD
    PREV(["From Part A: GWAS + Tractor-GENESIS + POOLED sumstats"])
    FM["FINE_MAPPING<br/>regions -> SuSiE-RSS / PolyFun+SuSiE within ancestry<br/>SuSiE-ME / MG-FLASH-FM multi-ancestry"]
    COLOC["COLOCALIZATION (parallel methods)<br/>coloc.susie | HyPrColoc | OPERA<br/>on the QTL megaset (eQTL, sQTL, pQTL, mQTL, caQTL, hQTL)"]
    H2["HERITABILITY<br/>LDSC h2 per stratum + POOLED, cross-ancestry rg"]
    PRS["PRS_WORKFLOW<br/>PRS-CSx + GAUDI, DiscoDivas, SDPR_admix, MUSSEL, PROSPER<br/>local-ancestry partial scores; validation overall AND per stratum"]
    FUNC["FUNCTIONAL_ANNOT<br/>MAGMA, ANNOVAR"]
    GXG["GXG_INTERACTION (end of pipeline)<br/>hits + known loci + custom list -> conditional selection<br/>pairwise SNP x SNP: pooled + per stratum<br/>ancestry modification of the interaction"]
    REP["REPORTING<br/>MultiQC + summary report"]
    PREV -->|filtered + POOLED sumstats| FM
    FM -->|credible sets| COLOC
    PREV -->|sumstats + meta| H2
    PREV -->|ancestry-specific sumstats<br/>+ genotypes, local ancestry| PRS
    PREV -->|meta sumstats| FUNC
    PREV -->|ALL sumstats: standard, Tractor, meta, survival<br/>+ full-cohort genotypes + phenotype| GXG
    FM --> REP
    COLOC --> REP
    H2 --> REP
    PRS --> REP
    FUNC --> REP
    GXG --> REP
    ${CLASSES}
    class PREV in; class FM,COLOC,H2,PRS,FUNC core; class GXG,REP end_;
` },

tractor: { wrap: 220, font: 20, code: `
flowchart LR
    subgraph IN["Inputs"]
        D["Tractor dosages<br/>Dose_EUR, Dose_AFR, Dose_AMR"]
        H["Tractor haplotype counts<br/>LA_AFR, LA_AMR (EUR = reference)"]
        P["phenotype + covariates (PC-AiR PCs)"]
        K["PC-Relate GRM"]
    end
    subgraph NULL["Null model (fitted once per trait)"]
        NB["binary / quantitative<br/>GENESIS fitNullModel<br/>kinship in Sigma"]
        NS["survival (OS, DFS)<br/>coxme: Cox + frailty on 2 x kinship<br/>Breslow baseline -> martingale residuals"]
    end
    E["Per-variant CONDITIONAL SCORE TEST<br/>D = [ LA | Dose ]<br/>U = D'Py &nbsp;&nbsp; M = D'PD<br/>Dose block conditioned on LA block<br/>joint k-df test &nbsp;|&nbsp; per-ancestry beta, SE, P<br/>covariance-aware heterogeneity &nbsp;|&nbsp; MAC gate"]
    O["tractor_genesis.tsv.gz<br/>P_JOINT, BETA/SE/P per ancestry,<br/>MAC per ancestry, P_HET, I2"]
    P --> NB
    K --> NB
    P --> NS
    K --> NS
    NB -->|Py, cholSigmaInv, CX, CXCXI| E
    NS -->|martingale residuals,<br/>risk-set information| E
    D --> E
    H --> E
    E --> O
    ${CLASSES}
    class D,H,P,K in; class NB core; class NS tr; class E,O end_;
` },

gxg: { wrap: 320, code: `
flowchart TD
    S["All sumstats for the trait<br/>standard + Tractor-GENESIS + meta + survival"]
    KL["Default known leukemia loci<br/>assets/known_leukemia_risk_loci.tsv<br/>(gxg_use_default_loci: on / off)"]
    CV["Custom GRCh38 list<br/>rsID | chr:pos | chr:pos:ref:alt<br/>(gxg_custom_variants)"]
    SEL["GXG_SELECT_HITS<br/>P < gxg_p_threshold, cap gxg_max_hits"]
    G["Full-cohort genotypes<br/>plink2 --export A on the hit list; IDs resolved allele-aware"]
    PR["Independent-signal selection (gxg_prune_mode = conditional)<br/>walk GWAS hits by P > custom > known; a hit in LD (r2 > 0.2) with kept hits<br/>is tested CONDITIONAL on them (y ~ hit + partners + covariates, LRT)<br/>kept if cond P < 1e-4; hits_pruned.tsv records status, proxy, cond P"]
    PAIRS["Pair list<br/>skip same-chromosome pairs < 1 Mb apart"]
    T0["POOLED<br/>y ~ g1 + g2 + g1:g2 + covariates + stratum"]
    T1["Per stratum (N >= 30)<br/>EUR, AAC, LAT1, LAT2, EAS, SAS, OTHER"]
    HET["Ancestry modification<br/>Cochran's Q / I2 across strata<br/>pooled 3-way LRT g1:g2:stratum"]
    LA["Optional Tractor dosages<br/>Dose_anc(SNP1) x Dose_anc(SNP2)"]
    OUT["gxg.stratum.tsv, gxg.all.tsv, gxg.significant.tsv<br/>gxg.ancestry_heterogeneity.tsv, gxg.local_ancestry.tsv<br/>gxg_summary.tsv (wide), gxg_top_interactions.tsv"]
    S --> SEL
    KL --> SEL
    CV --> SEL
    SEL --> G --> PR --> PAIRS
    PAIRS --> T0
    PAIRS --> T1
    T1 --> HET
    T0 --> HET
    PAIRS --> LA
    T0 --> OUT
    T1 --> OUT
    HET --> OUT
    LA --> OUT
    ${CLASSES}
    class S,KL,CV in; class SEL,G,PR,PAIRS core; class T0,T1,HET,LA tr; class OUT end_;
` },

slurm: { wrap: 320, code: `
flowchart TD
    S0["STEP 0 &nbsp; submit_pcair.sh<br/>PC-AiR PCs + PC-Relate GRM (once per cohort)"]
    S1["STEP 1 &nbsp; submit_gwas_array.sh<br/>traits x strata, 22-chromosome arrays<br/>Tractor-GENESIS for every trait"]
    S2["STEP 2 &nbsp; combine per trait<br/>chr* -> stratum.sumstats.gz<br/>meta_analyze_strata.R -> POOLED.sumstats.gz"]
    S3["STEP 3 &nbsp; submit_prs_array.sh<br/>6 PRS methods per trait + LA partial scores"]
    S4["STEP 4 &nbsp; colocalization array<br/>6 QTL types x 3 methods"]
    S5["STEP 5 &nbsp; PRS validation<br/>overall + per stratum"]
    S6["STEP 6 &nbsp; submit_gxg_array.sh<br/>per trait: POOLED + 7 strata,<br/>then ancestry-heterogeneity job"]
    R["results/gwas/trait/stratum/ &nbsp;|&nbsp; results/prs/trait/<br/>results/coloc/ &nbsp;|&nbsp; results/gxg/trait/"]
    S0 --> S1
    S1 --> S2
    S2 --> S3 --> S5
    S2 --> S4
    S2 --> S6
    S5 --> R
    S4 --> R
    S6 --> R
    ${CLASSES}
    class S0 tr; class S1,S2,S3,S4,S5,S6 core; class R end_;
` }
};

(async () => {
    const browser = await chromium.launch({ executablePath: exe, args: ['--no-sandbox'] });
    const page = await browser.newPage({ viewport: { width: 2400, height: 2400 }, deviceScaleFactor: 2 });
    await page.setContent('<html><body style="margin:0;background:#fff"><div id="d" style="display:inline-block;padding:8px"></div></body></html>');
    await page.addScriptTag({ path: mermaidPath });
    const meta = {};
    for (const [name, { code, wrap, font }] of Object.entries(diagrams)) {
        const size = await page.evaluate(async ({ code, wrap, font }) => {
            mermaid.initialize({
                startOnLoad: false, theme: 'neutral', securityLevel: 'loose',
                flowchart: { htmlLabels: true, useMaxWidth: false, nodeSpacing: 30, rankSpacing: 44, curve: 'basis', wrappingWidth: wrap },
                themeVariables: { fontSize: (font || 16) + 'px', fontFamily: 'Helvetica, Arial, sans-serif' }
            });
            const { svg } = await mermaid.render('m' + Math.random().toString(36).slice(2), code);
            document.getElementById('d').innerHTML = svg;
            const r = document.querySelector('#d svg').getBoundingClientRect();
            return { w: Math.ceil(r.width), h: Math.ceil(r.height) };
        }, { code, wrap, font });
        await (await page.$('#d')).screenshot({ path: `${name}.png` });
        meta[name] = size;
        console.log(name, size);
    }
    fs.writeFileSync('diagrams.json', JSON.stringify(meta));
    await browser.close();
})().catch(e => { console.error(e); process.exit(1); });
