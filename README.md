# Ancestry-Aware GWAS Pipeline

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A comprehensive Nextflow pipeline for multi-ancestry genome-wide association studies (GWAS) with fine-mapping, meta-analysis, polygenic risk score calculation, and functional annotation. **Optimized for admixed populations** with Tractor, GENESIS, and SAIGE support.

## Overview

This pipeline addresses the critical need for ancestry-aware genetic analysis by:

1. **Ancestry-stratified GWAS**: Runs separate GWAS for each major ancestry group (EUR, AFR, EAS, SAS, AMR, MID, AAC, AHI, LAT1, LAT2, HET)
2. **Admixed-optimized analysis**: GENESIS, SAIGE, and **Tractor** for local ancestry-aware GWAS
3. **Multi-ancestry meta-analysis**: Uses MR-MEGA to model ancestry-correlated heterogeneity
4. **Within-ancestry fine-mapping**: PolyFun+SuSIE for annotation-prioritized fine-mapping
5. **Multi-ancestry fine-mapping**: MG-FLASH-FM for related traits, SuSIE-ME for single traits
6. **Ancestry-aware PRS**: Multiple methods including PRS-CSx, GAUDI, and local ancestry-informed approaches
7. **Heritability analysis**: Ancestry-specific h² and cross-ancestry genetic correlations
8. **Functional annotation**: MAGMA, FUMA, LAVA, FLAMES with cell-type specificity
9. **Diverse QTL colocalization**: Curated datasets beyond GTEx for blood/immune, mQTL, caQTL, and more
10. **Survival analysis**: SPA-Cox for time-to-event GWAS

## Key Features for Admixed Populations

### Tractor - Local Ancestry-Aware GWAS
Tractor decomposes genetic effects by ancestral origin in admixed individuals:
- **African American (AAC)**: 2-way admixture (EUR, AFR)
- **Latino (LAT1+LAT2 combined)**: 3-way admixture (EUR, AFR, NAT)
  - LAT1 = Mexican/Central American
  - LAT2 = Caribbean/South American
  - Combined for better statistical power

### GENESIS - Admixed Population GWAS
R/Bioconductor package optimized for:
- Complex population structure
- Related individuals (with kinship matrix)
- Heterogeneous ancestry

### SAIGE
Optimized for:
- Case-control imbalance
- Rare variants
- Large biobank data

## Pipeline Architecture

```
┌─────────────────┐
│   Input Data    │
│  (Genotypes +   │
│   Phenotypes)   │
└────────┬────────┘
         │
         ▼
┌─────────────────┐     ┌──────────────────────┐
│   QC + Ancestry │────▶│  Optional Ancestry   │
│    (GRAF-ANC)   │     │  Inference (or use   │
│                 │     │  pre-computed calls) │
└────────┬────────┘     └──────────────────────┘
         │
    ┌────┴────┬────────┬────────┬────────┬────────┬────────┬────────┐
    ▼         ▼        ▼        ▼        ▼        ▼        ▼        ▼
┌───────┐ ┌───────┐ ┌───────┐ ┌───────┐ ┌───────┐ ┌───────┐ ┌───────┐ ┌───────┐
│  EUR  │ │  AFR  │ │  EAS  │ │  SAS  │ │  AMR  │ │  MID  │ │  AAC  │ │ LATINO│
│ GWAS  │ │ GWAS  │ │ GWAS  │ │ GWAS  │ │ GWAS  │ │ GWAS  │ │Tractor│ │Tractor│
└───┬───┘ └───┬───┘ └───┬───┘ └───┬───┘ └───┬───┘ └───┬───┘ └───┬───┘ └───┬───┘
    │         │        │        │        │        │        │        │
    └────┬────┴────────┴────────┴────────┴────────┴────────┴────────┘
         │
         ▼
┌─────────────────┐
│   MR-MEGA       │
│  Meta-Analysis  │
└────────┬────────┘
         │
    ┌────┴─────────────┐
    ▼                  ▼
┌─────────────┐  ┌─────────────┐
│ Within-Anc  │  │ Multi-Anc   │
│ Fine-Map    │  │ Fine-Map    │
│ PolyFun+    │  │ MG-FLASH-FM │
│ SuSIE       │  │ or SuSIE-ME │
└──────┬──────┘  └──────┬──────┘
       │                │
       └────────┬───────┘
                ▼
┌───────────────────────────────────────────────────────────┐
│                    Downstream Analysis                     │
├─────────────┬─────────────┬──────────────┬────────────────┤
│     PRS     │ Heritability│ Colocalization│   Functional   │
│  PRS-CSx    │    LDSC     │ Diverse QTL  │     MAGMA      │
│  PRS-CS     │    GCTA     │   Datasets   │     FUMA       │
│   GAUDI     │  Cross-Anc  │              │     LAVA       │
│  LDpred2    │     rg      │              │    FLAMES      │
└─────────────┴─────────────┴──────────────┴────────────────┘
```

## Quick Start

```bash
# Basic run with REGENIE
nextflow run main.nf \
  --input samplesheet.csv \
  --phenotype_file phenotypes.txt \
  --phenotype_cols "trait1,trait2" \
  --outdir results \
  -profile docker

# With pre-computed ancestry calls (recommended)
nextflow run main.nf \
  --input samplesheet.csv \
  --phenotype_file phenotypes.txt \
  --phenotype_cols "WBC,Lymphocytes,Monocytes" \
  --covariate_cols "sex,age_at_diagnosis,cd19_pos_pct,storage_days,PC1,PC2,PC3,PC4,PC5" \
  --run_ancestry_inference false \
  --ancestry_calls_file ancestry_calls.tsv \
  --gwas_tool genesis \
  --kinship_matrix grm.rds \
  --outdir results \
  -profile singularity

# For admixed populations with Tractor
nextflow run main.nf \
  --input samplesheet.csv \
  --phenotype_file phenotypes.txt \
  --phenotype_cols "trait1" \
  --run_tractor true \
  --local_ancestry_files local_ancestry.msp \
  --tractor_aac_pops "EUR,AFR" \
  --tractor_lat_pops "EUR,AFR,NAT" \
  --outdir results \
  -profile singularity

# For survival/time-to-event analysis
nextflow run main.nf \
  --input samplesheet.csv \
  --phenotype_file phenotypes.txt \
  --phenotype_cols "leukemia" \
  --survival_analysis true \
  --time_col "time_to_event" \
  --event_col "event_status" \
  --outdir results \
  -profile singularity
```

## Running on SLURM without Nextflow

The same `bin/` scripts also run as plain SLURM jobs from `slurm/`. The site settings to fill in
are in [docs/HPC_SITE.md](docs/HPC_SITE.md); the setup and check steps are in the combined runbook
(the harmonization repo's `genotyping_pipeline/benchmarking/SMOKE_RUNBOOK.md`, part "GWAS
pipeline"); the run record is described in [docs/WORKFLOW.md](docs/WORKFLOW.md), section 9. In short:

1. **Environments.** Build the four conda envs in `envs/` (`gwas-r`, `gwas-py`, `lai`, `mrmega`)
   into `$GWAS_ENVS_DIR`, run `bin/envs_tsv.sh` to record the sha256 of each yml in `envs.tsv`, and
   check them with `sbatch ... bin/check_envs.sh`. Every job starts with
   `source bin/activate_env.sh <env>`, which also puts `TMPDIR` on disk (`$SLURM_TMPDIR`, else
   `$GWAS_SCRATCH/tmp`), never on a tmpfs `/tmp`.
2. **Site settings.** `SLURM_PARTITION` and `SLURM_ACCOUNT` are required; the launcher refuses to
   start without them. No site path is written in any script.
   ```bash
   export SLURM_PARTITION=<partition> SLURM_ACCOUNT=<account> GWAS_ENVS_DIR=<dir> GWAS_SCRATCH=<dir>
   ```
3. **Launch**, from any directory; the dry run prints every sbatch line and submits nothing:
   ```bash
   slurm/submit_full_pipeline.sh --dry-run
   slurm/submit_full_pipeline.sh
   ```
   Step 0 ingests the population structure the handoff bundle ships (`HANDOFF_DIR`: cohort and
   within-stratum PCs, the GRM; `STRUCTURE_SOURCE=compute` recomputes it with PC-AiR / PC-Relate,
   as a comparator, only when asked). Step 0, the Tractor-GENESIS GWAS (traits × strata × 22
   chromosomes) and the combine + meta-analysis step always run. PRS, colocalization, PRS validation and GxG run only
   with `RUN_PRS=true`, `RUN_COLOC=true`, `RUN_VALIDATION=true` (needs `RUN_PRS`) and
   `RUN_GXG=true`; every launch prints the ones it deferred, e.g.
   `deferred: PRS COLOC VALIDATION GXG (RUN_*=false)`.
   Traits are `mrd` (binary), `efs` and `os` (time-to-event), the handoff contract's columns: each
   survival trait reads its own `<trait>_time` and `<trait>_event` (override one trait with
   `TIME_COL_<trait>` / `EVENT_COL_<trait>`); the launcher checks every column before queuing.
   The number-changing switches (all recorded in the run record) include `MIN_STRATUM_N`, the
   smallest stratum analysed (default 50, decided 10 Nov 2026), and `MAC_MIN`, the minor allele
   count an ancestry's dose term needs (default 10; the handoff contract says 20, still open).
4. **Ancestries per stratum.** Each stratum needs `data/tractor/<stratum>/ancestries.txt`, the
   ancestry labels its Tractor files carry; GxG's pooled task and PRS read the cohort-wide
   `data/tractor/ancestries.txt` ([docs/HANDOFF_CONTRACT.md](docs/HANDOFF_CONTRACT.md)). Tractor
   files are read by ancestry label, or through `data/tractor/flare_index_map.tsv` when numbered;
   never by their position in the list. The strata are the bundle's (`ancestry/strata.tsv`; set
   `STRATA_LIST` to run fewer). With `data/ancestry/strata.tsv` present, the adapter
   checks every sample against its stratum; a stratum under `MIN_STRATUM_N` analysable samples is
   skipped, and the run record says so.
5. **Run record.** Each launch writes `results/run_record/<run_id>/`: `params.tsv` (every
   number-changing switch) at launch, a row per task as it runs, and `trace.tsv` (with MaxRSS and
   Elapsed from `sacct`, and the adapter's sample and variant counts) and `summary.md` once the final
   job, which waits on all the others, has run.
6. **Switches.** Every number-changing switch has one default and a domain (`RR_SWITCHES`,
   `bin/run_record.sh`); a value outside it is refused before anything is queued, and `params.tsv`
   records every value. `RUN_LABEL=<name>` writes a sensitivity arm under `results/<name>/`;
   `PC_SCALE` and `MIX_RETAIN` are two such arms ([docs/WORKFLOW.md](docs/WORKFLOW.md) section 9).
7. **QC gates.** Steps file QC verdicts and per-join sample counts into the run record's `qc/`;
   each GWAS task requires the structure step's `pcair_pc_merge`, the combine job the adapter's
   join checks. A failed **or missing** verdict stops the job (exit 3);
   `QC_GATE_OVERRIDE=<id>` proceeds on purpose and is recorded
   ([docs/WORKFLOW.md](docs/WORKFLOW.md) section 10).

Logs go to `<repo>/logs/` (absolute paths; `LOG_DIR` to change). The wrappers under `shell/` are
debug-only, not the execution path. `tests/test_chassis.sh` checks all of this locally with fake
`sbatch`, `sacct` and `conda`.

## Input Requirements

### Samplesheet

Create a CSV file with genotype file paths:

```csv
sample_id,cohort,bed,bim,fam
study1,cohort_A,/path/to/study1.bed,/path/to/study1.bim,/path/to/study1.fam
study2,cohort_B,/path/to/study2.bed,/path/to/study2.bim,/path/to/study2.fam
```

Supported formats: PLINK1 (bed/bim/fam), PLINK2 (pgen/pvar/psam), BGEN, VCF

### Phenotype File

Tab-separated file with sample IDs and phenotype values:

```
IID	trait1	trait2	age	sex	PC1	PC2	PC3	PC4	PC5
sample1	25.3	170.2	45	1	0.01	-0.02	0.005	0.001	-0.003
sample2	28.1	165.8	38	2	-0.01	0.03	-0.002	0.004	0.001
```

### Pre-computed Ancestry Calls (Optional but Recommended)

If you have pre-computed ancestry calls, provide a TSV file:

```
IID	ancestry
sample1	EUR
sample2	AFR
sample3	AAC
sample4	LAT1
```

## Key Parameters

### Ancestry Options

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_ancestry_inference` | `false` | Run ancestry inference in pipeline |
| `--ancestry_calls_file` | `null` | Pre-computed ancestry calls file |
| `--ancestry_method` | `graf-anc` | Method: graf-anc, admixture, pca |
| `--ancestry_groups` | `EUR,AFR,EAS,SAS,AMR,MID,AAC,AHI,LAT1,LAT2,HET` | Ancestry groups |
| `--min_ancestry_n` | `100` | Minimum samples per ancestry |

### GWAS Tools

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--gwas_tool` | `regenie` | Tool: regenie, saige, bolt-lmm, plink2, **genesis** |
| `--gwas_model` | `additive` | Model: additive, dominant, recessive |
| `--kinship_matrix` | `null` | Pre-computed kinship matrix (for GENESIS) |
| `--p_threshold` | `5e-8` | Genome-wide significance threshold |

### Tractor (Local Ancestry-Aware GWAS)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--run_tractor` | `false` | Run Tractor for admixed populations |
| `--tractor_groups` | `AAC,LATINO` | Groups for Tractor (LATINO = LAT1+LAT2) |
| `--tractor_aac_pops` | `EUR,AFR` | Ancestral populations for African American |
| `--tractor_lat_pops` | `EUR,AFR,NAT` | Ancestral populations for Latino |
| `--local_ancestry_files` | `null` | Local ancestry MSP files (RFMix format) |

### Survival Analysis (SPA-Cox)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--survival_analysis` | `false` | Run time-to-event GWAS |
| `--time_col` | `null` | Time-to-event column name |
| `--event_col` | `null` | Event indicator column (1=event, 0=censored) |

### Meta-Analysis

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--meta_method` | `mr-mega` | Method: mr-mega, metal, metasoft |
| `--mr_mega_env` | `false` | Include environmental variable |

### Fine-Mapping

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--within_ancestry_fm` | `polyfun-susie` | Within-ancestry method |
| `--multi_ancestry_fm` | `mg-flash-fm` | Multi-ancestry method for related traits |
| `--single_trait_fm` | `susie-me` | Multi-ancestry method for single traits |
| `--related_traits` | `true` | Whether traits are related |

### Colocalization & QTL Datasets

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--use_curated_qtl` | `true` | Use curated diverse QTL datasets |
| `--qtl_ancestry_matching` | `true` | Match QTL datasets to GWAS ancestry |
| `--custom_qtl_datasets` | `null` | Path to custom QTL samplesheet |

**Curated QTL datasets include:**
- **Blood/Immune**: OneK1K (14 immune cell types), DICE, BLUEPRINT, eQTLGen
- **African ancestry**: GENOA (mQTL/eQTL), AFGR (caQTL)
- **Multi-ethnic**: MESA (EUR, AFR, HISP, EAS), TOPMed
- **Methylation**: GENOA mQTL, GoDMC, BLUEPRINT mQTL
- **Chromatin accessibility**: DICE caQTL, BLUEPRINT hQTL
- **Single-cell**: OneK1K, HCA Immune Atlas

### PRS Methods

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--prs_methods` | `prs-csx,prs-cs,gaudi,ldpred2` | PRS methods to run |
| `--prs_validation` | `true` | Validate PRS in target cohort |
| `--prs_best_method` | `true` | Determine best method per ancestry |

## Ancestry Groups

| Code | Description | Tractor Support |
|------|-------------|-----------------|
| EUR | European | - |
| AFR | African | - |
| EAS | East Asian | - |
| SAS | South Asian | - |
| AMR | Native American | - |
| MID | Middle Eastern | - |
| AAC | African American (admixed) | Yes (2-way: EUR,AFR) |
| AHI | American Hispanic (admixed) | Yes (3-way: EUR,AFR,NAT) |
| LAT1 | Latino Type 1 (Mexican/Central American) | Combined as LATINO |
| LAT2 | Latino Type 2 (Caribbean/South American) | Combined as LATINO |
| HET | Heterogeneous/Other | - |

## Output Structure

```
results/
├── ancestry/
│   ├── graf/                    # GRAF ancestry calls
│   └── stratified/              # Ancestry-stratified genotypes
├── gwas/
│   ├── {ancestry}/{trait}/      # Per-ancestry GWAS results
│   └── tractor/                 # Tractor results (AAC, LATINO)
│       ├── AAC/                 # 2-way admixture results
│       └── LATINO/              # 3-way admixture results (LAT1+LAT2)
├── survival/
│   └── {ancestry}/{trait}/      # SPA-Cox survival GWAS
├── meta_analysis/
│   └── {trait}/                 # MR-MEGA meta-analysis
├── fine_mapping/
│   ├── {ancestry}/{trait}/      # Within-ancestry fine-mapping
│   └── multi_ancestry/{trait}/  # Multi-ancestry fine-mapping
├── prs/
│   └── {trait}/                 # PRS weights and scores
├── heritability/
│   └── {ancestry}/{trait}/      # h² estimates
├── functional/
│   └── {trait}/                 # Functional annotations
├── colocalization/
│   └── {trait}/                 # Colocalization results
├── plots/
│   └── {trait}/                 # Manhattan, QQ, regional plots
└── reports/
    ├── multiqc_report.html      # MultiQC report
    └── ancestry_gwas_report.html # Custom pipeline report
```

## Methods Description

### Tractor
Local ancestry-aware GWAS that decomposes genetic effects by ancestral origin:
- Extracts ancestry-specific haplotype dosages from local ancestry calls
- Tests for ancestry-specific effects and heterogeneity
- Particularly powerful for identifying ancestry-specific associations

### GENESIS
R/Bioconductor package for:
- Mixed model association testing accounting for population structure and relatedness
- Score, Wald, and BinomiRare tests
- PC-AiR for ancestry inference in related samples
- PC-Relate for kinship estimation in structured populations

### SPA-Cox
Saddlepoint approximation for Cox proportional hazards GWAS:
- Efficient for large-scale survival analysis
- Handles time-to-event outcomes (e.g., disease onset, death)
- Accounts for censoring

### MR-MEGA Meta-Analysis
Meta-regression of genetic association data:
- Models allelic effects as a function of axes of genetic variation
- Accounts for ancestry-correlated heterogeneity
- Detects residual heterogeneity beyond ancestry

### PolyFun+SuSIE Fine-Mapping
Annotation-informed fine-mapping:
- PolyFun: Functional annotation-informed prior probabilities
- SuSIE: Sum of Single Effects regression for multiple causal variants
- 95% credible set construction

### MG-FLASH-FM & SuSIE-ME
Multi-ancestry fine-mapping:
- **MG-FLASH-FM**: For related traits (shares information across traits)
- **SuSIE-ME**: For single/unrelated traits

### PRS-CSx & GAUDI
Cross-population PRS:
- **PRS-CSx**: Coupled continuous shrinkage priors across ancestries
- **GAUDI**: Local ancestry-informed PRS for admixed individuals

## Citation

If you use this pipeline, please cite:

```
[Pipeline citation - to be added]

Key methods:
- REGENIE: Mbatchou et al. (2021) Nature Genetics
- SAIGE: Zhou et al. (2018) Nature Genetics
- GENESIS: Gogarten et al. (2019) Bioinformatics
- Tractor: Atkinson et al. (2021) Nature Genetics
- MR-MEGA: Magi et al. (2017) Nature Genetics
- PolyFun: Weissbrod et al. (2020) Nature Genetics
- SuSIE: Wang et al. (2020) J R Stat Soc
- PRS-CSx: Ruan et al. (2022) Nature Genetics
- GAUDI: Marnetto et al. (2020) Nature Communications
```

## Requirements

- Nextflow >= 23.04.0
- Container runtime: Docker, Singularity, or Apptainer
- OR Conda/Mamba for local execution
- For the SLURM path without Nextflow: Slurm and conda (or mamba), with the envs in `envs/`
  (see "Running on SLURM without Nextflow" above)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please read our [contributing guidelines](CONTRIBUTING.md) first.

## Support

- [Open an issue](https://github.com/your-org/ancestry-aware-gwas/issues)
- [Documentation](https://your-org.github.io/ancestry-aware-gwas/)
