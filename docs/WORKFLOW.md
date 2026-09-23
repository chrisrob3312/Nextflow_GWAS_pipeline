# Pipeline workflow: what happens to the data, module by module

This is the current data flow for the multi-ancestry leukemia GWAS pipeline.
Boxes are Nextflow subworkflows or processes; every arrow is a file handed from
one module to the next. The same steps exist as SLURM array scripts under
`slurm/` (section 6), which call the same `bin/` scripts.

## 1. End-to-end overview

```mermaid
flowchart TD
    IN["INPUTS<br/>samplesheet (PLINK bed/bim/fam)<br/>phenotype TSV (sample_id, traits, covariates, GRAF_ANC)<br/>kinship RDS<br/>local ancestry MSP (RFMix)"]
    QC["QC_WORKFLOW<br/>MAF / HWE / missingness filters"]
    ANC["ANCESTRY_INFERENCE<br/>GRAF-ANC -> strata<br/>EUR, AAC, LAT1, LAT2, EAS, SAS (N>=30) else OTHER"]
    GWAS["GWAS_WORKFLOW<br/>per stratum x trait<br/>REGENIE / SAIGE / GENESIS / SPA-Cox"]
    TR["TRACTOR-GENESIS<br/>AAC (EUR,AFR) + LATINO (EUR,AFR,AMR)<br/>ALL traits: MRD binary, OS/DFS time-to-event"]
    META["META_ANALYSIS<br/>MR-MEGA across strata -> POOLED"]
    VIS["VISUALIZATION<br/>Manhattan, QQ, Miami"]
    FM["FINE_MAPPING<br/>PolyFun+SuSiE (within) / SuSiE-ME, MG-FLASH-FM (multi)"]
    COLOC["COLOCALIZATION (parallel)<br/>coloc.susie | HyPrColoc | OPERA<br/>on the QTL megaset"]
    H2["HERITABILITY<br/>cov-LDSC, GCTA, local-ancestry h2"]
    PRS["PRS_WORKFLOW<br/>PRS-CSx + GAUDI, DiscoDivas, SDPR_admix, MUSSEL, PROSPER<br/>validation: overall AND per stratum"]
    FUNC["FUNCTIONAL_ANNOT<br/>MAGMA, FUMA, LAVA, FLAMES"]
    GXG["GXG_INTERACTION (end of pipeline)<br/>hits + known loci + custom list -> LD prune<br/>pairwise SNP x SNP: pooled + per stratum<br/>ancestry modification of the interaction"]
    REP["REPORTING<br/>MultiQC + summary report"]

    IN --> QC --> ANC
    ANC -->|"stratified bed/bim/fam"| GWAS
    ANC -->|"stratified genotypes + MSP"| TR
    IN -->|"phenotype, kinship"| GWAS
    IN -->|"phenotype, kinship"| TR
    GWAS -->|"sumstats per stratum x trait"| META
    TR -->|"tractor_genesis.tsv.gz per group x trait"| META
    GWAS --> VIS
    GWAS -->|"filtered sumstats"| FM
    META -->|"POOLED sumstats"| FM
    FM -->|"credible sets"| COLOC
    GWAS --> H2
    META --> H2
    META -->|"ancestry-specific sumstats"| PRS
    ANC -->|"stratified genotypes, local ancestry"| PRS
    META --> FUNC
    GWAS -->|"all sumstats (standard + Tractor + meta + survival)"| GXG
    QC -->|"full-cohort bed/bim/fam"| GXG
    IN -->|"phenotype (traits, covariates, GRAF_ANC)"| GXG
    TR -.->|"optional: ancestry dosages"| GXG
    META --> REP
    FM --> REP
    H2 --> REP
    PRS --> REP
    GXG --> REP
```

## 2. Module inputs and outputs

| Module (subworkflow / process) | Reads | Writes | Feeds |
|---|---|---|---|
| `INPUT_CHECK` | samplesheet CSV | validated genotype channel | QC |
| `QC_WORKFLOW` | bed/bim/fam | QC'd bed/bim/fam, QC report | ANCESTRY, GxG (full cohort) |
| `ANCESTRY_INFERENCE` | QC'd genotypes, GRAF reference, MSP | `ancestry_calls.tsv`, stratified bed/bim/fam per group, local ancestry channel | GWAS, Tractor, PRS |
| `GWAS_WORKFLOW` (standard) | stratified genotypes x trait, phenotype, covariates, kinship | `<id>.<ancestry>.<trait>.sumstats.gz`, filtered sumstats, significant variants | META, FM, H2, VIS, GxG |
| `PLINK_TO_VCF` -> `TRACTOR_EXTRACT_TRACTS` | AAC and merged LATINO genotypes + MSP; `meta.tractor_pops` | `<id>.ancdose.<k>.tsv.gz`, `<id>.hapcount.<k>.tsv.gz` per ancestry k | TRACTOR_GENESIS |
| `TRACTOR_GENESIS` (`bin/tractor_genesis_adapter.R`) | dosages + hapcounts, phenotype, kinship, covariates; `meta.model` | `<prefix>.tractor_genesis.tsv.gz` (sorted by P_JOINT), `.genomic_order.tsv.gz`, `.top_hits.tsv`, `.significant.tsv`, `.heterogeneous.tsv`, `.<anc>_specific.tsv`, `.null_model.rds` | META, GxG |
| `META_ANALYSIS` | sumstats grouped by trait across strata | MR-MEGA results, heterogeneity table; SLURM path: `bin/meta_analyze_strata.R` -> `POOLED.sumstats.gz` | FM, H2, PRS, FUNC, GxG |
| `FINE_MAPPING` | filtered sumstats, LD reference / cohort LD | credible sets, PIPs | COLOC |
| `COLOCALIZATION` | credible sets or filtered sumstats, QTL megaset (`assets/qtl_datasets/curated_qtl_sources.yml` via `bin/download_qtl_datasets.R`, harmonized by `bin/harmonize_qtl.R`), cohort LD | `*.coloc.tsv`, `*.hyprcoloc.tsv`, `*.opera.tsv`, combined + consensus tables | REP |
| `HERITABILITY` | sumstats, meta results, LD reference | h2 per stratum, local-ancestry h2, genetic correlations | REP |
| `PRS_WORKFLOW` (`bin/calculate_prs_admixed.R` on the SLURM path) | ancestry-specific sumstats, LD reference, stratified genotypes, local ancestry, phenotype | weights and `.sscore` per method, `<prefix>.validation.tsv` (OVERALL + one row per stratum, AUC / R2 / C-index, disparity flags) | REP |
| `FUNCTIONAL_ANNOT` | meta or filtered sumstats | MAGMA / FUMA / LAVA / FLAMES outputs | REP |
| `GXG_INTERACTION` (`bin/run_gxg_interaction.R`) | see section 4 | see section 4 | REP |
| `REPORTING` | everything above | HTML report | user |

## 3. Tractor-GENESIS: one engine for every trait

Every trait goes through the same conditional score test so MRD, DFS and OS
are comparable. The trait's `model` (set in `main.nf` from `binary_traits` and
`survival_traits`) only changes where the projection comes from.

```mermaid
flowchart LR
    subgraph inputs
        D["Tractor dosages<br/>Dose_EUR, Dose_AFR, Dose_AMR"]
        H["Tractor hapcounts<br/>LA_AFR, LA_AMR (EUR = reference)"]
        P["phenotype + covariates"]
        K["kinship"]
    end
    subgraph null["null model (fitted once)"]
        NB["binary / quantitative<br/>GENESIS fitNullModel<br/>kinship in Sigma"]
        NS["survival (OS, DFS)<br/>coxme: Cox + frailty on 2*kinship<br/>Breslow baseline -> martingale residuals"]
    end
    E["per-variant conditional score test<br/>D = [LA | Dose]<br/>U = D'Py, M = D'PD<br/>Dose block conditioned on LA block<br/>joint k-df test, per-ancestry beta/SE,<br/>covariance-aware heterogeneity, MAC gate"]
    O["tractor_genesis.tsv.gz<br/>P_JOINT, BETA/SE/P per ancestry, MAC per ancestry, P_HET, I2"]
    P --> NB
    K --> NB
    P --> NS
    K --> NS
    NB -->|"Py, cholSigmaInv, CX, CXCXI"| E
    NS -->|"martingale residuals, risk-set information"| E
    D --> E
    H --> E
    E --> O
```

## 4. GxG interaction (end of pipeline)

```mermaid
flowchart TD
    S["all sumstats for the trait<br/>standard + Tractor-GENESIS + meta + survival"]
    KL["default known leukemia loci<br/>assets/known_leukemia_risk_loci.tsv<br/>(gxg_use_default_loci)"]
    CV["custom GRCh38 list<br/>rsID | chr:pos | chr:pos:ref:alt<br/>(gxg_custom_variants)"]
    SEL["GXG_SELECT_HITS<br/>P < gxg_p_threshold, cap gxg_max_hits<br/>priority: custom > known > GWAS"]
    G["full-cohort genotypes<br/>plink2 --export A on the hit list"]
    PR["LD prune (gxg_prune_mode = variant)<br/>walk custom > known > best P;<br/>drop a hit if r2 > 0.2 with a kept hit on the same chr<br/>-> hits_pruned.tsv records the proxy"]
    PAIRS["pair list<br/>skip same-chr pairs < 1 Mb apart"]
    T0["POOLED<br/>y ~ g1 + g2 + g1:g2 + covariates + stratum"]
    T1["per stratum (N >= 30)<br/>EUR, AAC, LAT1, LAT2, EAS, SAS, OTHER"]
    HET["ancestry modification<br/>Cochran's Q / I2 across strata<br/>pooled 3-way LRT g1:g2:stratum"]
    LA["optional Tractor dosages<br/>Dose_anc(SNP1) x Dose_anc(SNP2)"]
    OUT["gxg.<stratum>.tsv, gxg.all.tsv, gxg.significant.tsv<br/>gxg.ancestry_heterogeneity.tsv, gxg.local_ancestry.tsv<br/>gxg_summary.tsv (wide), gxg_top_interactions.tsv"]
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
```

Models: logistic for binary traits, linear for quantitative, Cox for survival.
Each pair reports Wald and LRT p-values with Bonferroni over pairs and BH FDR.

## 5. Traits, covariates, strata: what is customisable

| Setting | Nextflow | SLURM | Notes |
|---|---|---|---|
| Traits | `phenotype_cols` | `TRAITS_LIST` | any phenotype columns |
| Binary traits | `binary_traits` | everything not in `SURVIVAL_TRAITS` | 0/1 |
| Time-to-event traits | `survival_traits` (+ `survival_time_cols`, `survival_event_cols`) | `SURVIVAL_TRAITS` | default columns `<trait>_time`, `<trait>_status` |
| Covariates | `covariate_cols`, `gxg_covariates` | `COVARIATES` | any phenotype columns |
| Ancestry column | `gxg_ancestry_col` | `ANCESTRY_COL` | GRAF-ANC codes or labels |
| Strata rules | `min_stratum_n`, `bin/ancestry_config.R` | `STRATA_LIST`, `bin/ancestry_config.R` | N < 30 pools to OTHER |
| Known loci | `gxg_known_loci`, `gxg_use_default_loci` | `KNOWN_LOCI`, `GXG_USE_DEFAULT_LOCI` | default table is a starting point; verify positions |
| Custom variants | `gxg_custom_variants` | `GXG_CUSTOM_VARIANTS` | GRCh38, one per line |
| Pair filtering | `gxg_prune_mode`, `gxg_max_pair_r2`, `gxg_min_distance_kb` | `GXG_PRUNE_MODE` | variant = keep strongest per LD cluster |
| Tractor | `tractor_aac_pops`, `tractor_lat_pops`, `tractor_ref_ancestry`, `tractor_mac_min` | `ANCESTRIES` case in `submit_gwas_array.sh` | |

## 6. SLURM orchestration (`slurm/submit_full_pipeline.sh`)

```mermaid
flowchart TD
    S1["STEP 1  submit_gwas_array.sh<br/>traits x strata, 22-chromosome arrays<br/>Tractor-GENESIS for every trait"]
    S2["STEP 2  combine per trait<br/>chr* -> <stratum>.sumstats.gz<br/>bin/meta_analyze_strata.R -> POOLED.sumstats.gz"]
    S3["STEP 3  submit_prs_array.sh<br/>6 PRS methods per trait"]
    S4["STEP 4  colocalization array<br/>6 QTL types x 3 methods"]
    S5["STEP 5  PRS validation<br/>overall + per stratum"]
    S6["STEP 6  submit_gxg_array.sh<br/>per trait: POOLED + 7 strata,<br/>then ancestry-heterogeneity job"]
    S1 --> S2
    S2 --> S3 --> S5
    S2 --> S4
    S2 --> S6
```

Results land under `results/gwas/<trait>/<stratum>/`, `results/prs/<trait>/`,
`results/coloc/`, `results/gxg/<trait>/`.

## 7. Known drift between the two paths

The Nextflow PRS and colocalization subworkflows (`subworkflows/prs`,
`subworkflows/colocalization`) still call the older per-method modules
(`prs_csx`, `prs_cs`, `gaudi`, `ldpred2`, single-method coloc). The SLURM path
calls the newer scripts (`bin/calculate_prs_admixed.R` with six methods and
stratified validation; `bin/run_colocalization.R` with parallel coloc.susie,
HyPrColoc and OPERA). The GWAS, Tractor-GENESIS and GxG steps are aligned on
both paths. Aligning PRS and colocalization is the next piece of dual-structure
work.
