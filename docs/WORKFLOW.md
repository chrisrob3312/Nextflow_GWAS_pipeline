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

    PCA["GENESIS_PCAIR_PCRELATE<br/>PC-AiR ancestry PCs + PC-Relate GRM<br/>used by every model"]
    IN --> QC --> ANC
    QC --> PCA
    ANC -->|"stratified bed/bim/fam"| GWAS
    ANC -->|"stratified genotypes + MSP"| TR
    PCA -->|"GRM + PC covariates"| GWAS
    PCA -->|"GRM + PC covariates"| TR
    PCA -->|"PC covariates"| GXG
    IN -->|"phenotype"| GWAS
    IN -->|"phenotype"| TR
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
| `QC_WORKFLOW` | bed/bim/fam | QC'd bed/bim/fam, QC report | ANCESTRY, PC-AiR, GxG (full cohort) |
| `GENESIS_PCAIR_PCRELATE` (`bin/genesis_pcair_pcrelate.R`) | QC'd full cohort, phenotype | PC-AiR PCs (`.pcair.pcs.tsv`), PC-Relate GRM (`.pcrelate.grm.rds`), kinship, unrelated set, phenotype with PC1..PCn replaced | GWAS null models, Tractor-GENESIS, GxG, PRS validation (every model uses the same GRM and PCs) |
| `ANCESTRY_INFERENCE` | QC'd genotypes, GRAF reference, MSP | `ancestry_calls.tsv`, stratified bed/bim/fam per group, local ancestry channel | GWAS, Tractor, PRS |
| `GWAS_WORKFLOW` (standard) | stratified genotypes x trait, phenotype, covariates, kinship | `<id>.<ancestry>.<trait>.sumstats.gz`, filtered sumstats, significant variants | META, FM, H2, VIS, GxG |
| `PLINK_TO_VCF` -> `TRACTOR_EXTRACT_TRACTS` | AAC and merged LATINO genotypes + MSP; `meta.tractor_pops` | `<id>.ancdose.<k>.tsv.gz`, `<id>.hapcount.<k>.tsv.gz` per ancestry k | TRACTOR_GENESIS |
| `TRACTOR_GENESIS` (`bin/tractor_genesis_adapter.R`) | dosages + hapcounts, phenotype, kinship, covariates; `meta.model` | `<prefix>.tractor_genesis.tsv.gz` (sorted by P_JOINT), `.genomic_order.tsv.gz`, `.top_hits.tsv`, `.significant.tsv`, `.heterogeneous.tsv`, `.<anc>_specific.tsv`, `.null_model.rds` | META, GxG |
| `META_ANALYSIS` | sumstats grouped by trait across strata | MR-MEGA results, heterogeneity table; SLURM path: `bin/meta_analyze_strata.R` -> `POOLED.sumstats.gz` | FM, H2, PRS, FUNC, GxG |
| `FINE_MAPPING` | filtered sumstats, LD reference / cohort LD | credible sets, PIPs | COLOC |
| `COLOCALIZATION` (`COLOC_RUN` = `bin/run_colocalization.R`) | every ancestry-stratified GWAS (standard + Tractor-GENESIS) and the meta/POOLED GWAS; QTL megaset pooled by type (`curated_qtl_sources.yml` via `download_qtl_datasets.R`, `harmonize_qtl.R`); optional cohort LD | per GWAS: `.coloc_susie.tsv`, `.hyprcoloc.tsv`, `.opera.tsv`; per trait: `coloc_combined.tsv`, `coloc_consensus.tsv`, shared / divergent / meta-only tables, PP4 heatmap | REP |
| `HERITABILITY` | sumstats, meta results, LD reference | h2 per stratum, local-ancestry h2, genetic correlations | REP |
| `PRS_WORKFLOW` (`PRS_METHOD` = `bin/calculate_prs_admixed.R`) | ancestry-specific sumstats per trait, LD reference, full-cohort genotypes, local ancestry, phenotype with PC-AiR PCs | per method: weights and `.sscore`; `<trait>.la_partial.scores.tsv` (Tractor dosages x Tractor-GENESIS betas: the part of each score on EUR / AFR / AMR haplotypes); `<trait>.validation.tsv` (OVERALL + one row per stratum, AUC / R2 / C-index); `prs_comparison.tsv`, `best_prs.tsv` (best overall and per stratum, disparity flags) | REP |
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
    S0["STEP 0  submit_pcair.sh<br/>PC-AiR PCs + PC-Relate GRM (once)"]
    S1["STEP 1  submit_gwas_array.sh<br/>traits x strata, 22-chromosome arrays<br/>Tractor-GENESIS for every trait"]
    S0 --> S1
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

## 7. Dual structure: both paths call the same scripts

Every analysis step now has one implementation in `bin/` that both the
Nextflow modules and the SLURM scripts call:

| Step | Script | Nextflow process | SLURM |
|---|---|---|---|
| PC-AiR + PC-Relate | `genesis_pcair_pcrelate.R` | `GENESIS_PCAIR_PCRELATE` | `submit_pcair.sh` (step 0) |
| Tractor-GENESIS GWAS | `tractor_genesis_adapter.R` | `TRACTOR_GENESIS` | `submit_gwas_array.sh` |
| Strata meta-analysis | `meta_analyze_strata.R` | (MR-MEGA in Nextflow) | step 2 |
| PRS (6 methods, LA partial, stratified validation) | `calculate_prs_admixed.R` | `PRS_METHOD`, `PRS_LA_PARTIAL`, `PRS_VALIDATE_STRATIFIED` | `submit_prs_array.sh` |
| Colocalization (coloc.susie, HyPrColoc, OPERA) | `run_colocalization.R` | `COLOC_RUN`, `COLOC_COMBINE` | step 4 |
| GxG interaction | `run_gxg_interaction.R` | `GXG_*` | `submit_gxg_array.sh` |

The older single-method modules (`prs_csx`, `prs_cs`, `gaudi`, `ldpred2`,
`coloc`, `ecaviar`, `fastenloc`) remain in `modules/local/` but are no longer
wired into the workflow.

## 8. What comes out per ancestry

Tractor-GENESIS reports, for every variant and trait: `P_JOINT` (any effect),
and for each ancestral background `BETA_<anc>`, `SE_<anc>`, `Z_<anc>`,
`P_<anc>`, `MAC_<anc>`, plus `P_HET` and `I2` for whether the effects differ
across backgrounds. The per-stratum SLURM runs add stratum-level results
(EUR, AAC, LAT1, LAT2, ...) on top, and `meta_analyze_strata.R` pools each
background across strata. PRS adds the local-ancestry partial scores, and
GxG reports every interaction pooled, per stratum, and for ancestry
modification.
