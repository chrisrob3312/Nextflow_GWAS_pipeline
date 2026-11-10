# Summary statistics: what the GWAS writes

This file describes the tables the Tractor-GENESIS step writes
(`results/gwas/<trait>/<stratum>/chr<N>.tractor_genesis.*`, combined into
`results/gwas/<trait>/<stratum>.sumstats.gz` and pooled into `POOLED.ivw_descriptive.sumstats.gz`).

## Columns of a stratum's table

Every effect is per **ALT** allele (`ALT` is the effect allele, `REF` the other). Order as written:

| column | what |
|---|---|
| `STRATUM` | the stratum (`--stratum`); the first column |
| `CHR`, `POS`, `ID`, `REF`, `ALT` | the variant, from the Tractor files (`ID` = `chr:pos:ref:alt`) |
| `N` | samples in the test (complete cases; sporadic missing doses are mean-imputed) |
| `MAC_STRATUM` | minor allele count in the stratum, `min(ALT, 2N - ALT)`; below `MAC_MIN` nothing is tested |
| `EAF` | ALT allele frequency in the stratum, `ALT / 2N` |
| `BETA_TOTAL`, `SE_TOTAL`, `P_TOTAL` | the **total-dose 1-df score test**: the summed ALT dose of all ancestries, no local-ancestry term, the same null model (CP-19). One effect per variant and stratum, the standard single-effect GWAS test; what MR-MEGA pools. With one ancestry it is the conditional test itself (`BETA_TOTAL` = `BETA_<anc>`) |
| `STAT_JOINT`, `DF_JOINT`, `P_JOINT` | the joint test of all testable dose terms, given local ancestry (`DF_JOINT` = number of testable backgrounds; with one, `P_JOINT` = `P_<anc>`) |
| `MAC_<anc>` | minor allele count on that background, `min(ALT_anc, haplotypes_anc - ALT_anc)`; its dose term is tested only at `MAC_MIN` or more |
| `EAF_<anc>` | ALT frequency on that background, `ALT_anc / haplotypes_anc` (`NA` with no haplotypes) |
| `BETA_<anc>`, `SE_<anc>`, `Z_<anc>`, `P_<anc>` | that background's dose effect, conditional on local ancestry and the other doses (one-step score estimate) |
| `Q_HET`, `DF_HET`, `P_HET`, `I2` | covariance-aware heterogeneity of the background effects (`DF_HET` = testable backgrounds - 1) |

`BETA` columns are log odds ratios (binary `mrd`), log hazard ratios (`efs`, `os`) or trait units
(quantitative), from a one-step score estimate: accurate near zero, conservative for large effects.

## The MR-MEGA input (`bin/sumstats_to_mrmega.py`)

One file per stratum, `results/mrmega/<trait>/<S>.mrmega_in.tsv`, MR-MEGA's column names:

| column | from |
|---|---|
| `MARKERNAME` | `ID` |
| `CHROMOSOME`, `POSITION` | `CHR` (without `chr`), `POS` |
| `EA`, `NEA` | `ALT`, `REF` |
| `EAF`, `N` | `EAF`, `N` |
| `BETA`, `SE` | `BETA_TOTAL`, `SE_TOTAL` (rows without them are left out); MR-MEGA runs with `--qt` |

Its `##` lines name the sources (`beta_source=BETA_TOTAL`, `se_source=SE_TOTAL`, `ea=ALT`,
`nea=REF`, the input file, the stratum) and how many rows passed the orientation check: with the
bundle's `.bim` (`--variants`), every row found there (by chromosome, position and allele pair)
must have `NEA` = its ref and `EA` = its alt, or the conversion stops. MR-MEGA is given copies
without the `##` lines. `slurm/submit_mrmega.sh` writes `<trait>.mrmega.result.tsv` (K >= 3 strata;
`##` lines: method, pc_axes, K, strata, beta_source) or, for K = 2, `<trait>.fixed_effect.tsv`
("pc_axes=0: fixed-effect meta-analysis, not MR-MEGA").

## Genotype dose: hard calls in every stratum, including k = 1

One rule for every stratum (decision M-06):

- The allele dose an ancestry term carries is the **hard-call** count of ALT alleles on that
  ancestry's haplotypes, from the phased genotypes (`GT`) that local ancestry was called on: the
  Tractor `dosage` files, 0, 1 or 2 per sample, alongside the `hapcount` files (the number of
  haplotypes of that ancestry, 0, 1 or 2).
- This holds for a stratum with one ancestry (k = 1, e.g. EUR) too: its single dose term is the
  hard-call genotype from the same phased calls. It does **not** switch to the imputed dosage
  (`FORMAT/DS` of the dosage VCFs), so a k = 1 stratum's effect sizes are on the same scale as the
  ancestry-specific effects of a k = 2 or 3 stratum, and the pooled meta-analysis combines like
  with like.
- The imputed dosages (`FORMAT/DS`, and the per-haplotype `FORMAT/HDS`) stay in the handoff bundle
  for reporting and for any later dosage-based analysis; no step of this pipeline reads them into a
  model.
- Consequence to state in methods: imputation uncertainty is not carried into the effect
  estimates (a hard call is used where the imputed dosage was, say, 1.4). Quality filtering
  happened upstream with MagicalRsq-X (`docs/HANDOFF_CONTRACT.md`, rule 5).

## One sample key

The `.fam` IID is the only sample key (`docs/HANDOFF_CONTRACT.md`, rule 1): VCF header sample IDs,
`.fam` IID and `phenotypes.tsv` IID / `sample_id` must agree. `bin/check_handoff_bundle.sh` refuses a
bundle where they do not, with `sample ID mismatch: <file> vs fam: ...`, and names PLINK's FID_IID
VCF export (sample names such as `F_S001_S001`) when that is the cause. With `--qc-dir` (or
`GWAS_QC_DIR` set) it also files `handoff_sample_ids` (PASS or FAIL) and a sample-accounting row.

## `## key=value` lines

A table may start with lines `## key=value` (no tab in them) before its column header:

| key | written by | value |
|---|---|---|
| `stratum`, `trait`, `model`, `n` | the adapter, every table | the stratum, the trait, the model (binary, quantitative, survival) and the samples in the null model |
| `method` | the adapter, every table | the tests (conditional Tractor score test + total-dose 1-df test) and the null model actually fitted (GENESIS logistic / linear, with or without the GRM; Cox with coxme frailty on the GRM, or Cox without it) |
| `covariates` | the adapter, every table; kept once by the combine step; `meta_analyze_strata.R` | the covariate list the null model used, global-ancestry columns included (`none` when there are none); `covariates_<stratum>` per stratum in the pooled file when the strata differ |
| `strata` | `meta_analyze_strata.R` (pooled file) | the strata pooled |
| `strata_skipped` | `meta_analyze_strata.R` (pooled file) | the strata left out because they were skipped below `MIN_STRATUM_N` (`none` when none) |
| `backgrounds` | `meta_analyze_strata.R` (pooled file) | the ancestry backgrounds pooled (from the `BETA_<anc>` columns) |

`data.table::fread` skips these lines by itself; `bin/sumstats_io.R` reads them
(`sumstats_meta()`) and writes them (`fwrite_sumstats()`); `bin/filter_gwas_results.py` skips them.
