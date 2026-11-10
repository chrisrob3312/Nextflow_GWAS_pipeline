# Summary statistics: what the GWAS writes

This file describes the tables the Tractor-GENESIS step writes
(`results/gwas/<trait>/<stratum>/chr<N>.tractor_genesis.*`, combined into
`results/gwas/<trait>/<stratum>.sumstats.gz` and pooled into `POOLED.sumstats.gz`).

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
| `covariates` | the adapter, every table; kept once by the combine step; `meta_analyze_strata.R` | the covariate list the null model used, global-ancestry columns included (`none` when there are none); `covariates_<stratum>` per stratum in the pooled file when the strata differ |

`data.table::fread` skips these lines by itself; `bin/sumstats_io.R` reads them
(`sumstats_meta()`) and writes them (`fwrite_sumstats()`); `bin/filter_gwas_results.py` skips them.
