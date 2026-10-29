# Handoff contract — harmonization pipeline → GWAS pipeline

**One directory, one manifest, no guessing.** The harmonization pipeline writes everything the GWAS
pipeline needs into `handoff/<run_id>/`, lists it in `manifest.tsv`, and `check_handoff_bundle.sh`
refuses the bundle if anything is missing, empty, mis-shaped, or mis-keyed. Every file below is named
exactly as written; the chromosome placeholder `{N}` is 1–22 with no `chr` prefix in the file name
(contigs *inside* VCFs are `chr`-prefixed, hg38).

| File | Role downstream | Format and required content | Produced by |
|---|---|---|---|
| `manifest.tsv` | the contract itself | tab-separated, header `path sha256 bytes role`; one row per file below; header lines `# run_id`, `# upstream_git_sha`, `# imputation_server` (topmed \| michigan \| anvil), `# duplicates_removed` (yes \| no), `# build` (hg38) | `make_handoff_manifest.sh` |
| `genotypes/cohort.hardcall.nomaf.bed/.bim/.fam` | the GWAS set: stratum keep-files, the per-stratum MAC gate, GrafAnc audit | PLINK 1 trio; hg38; QC-complete, **no MAF filter** — the minor-allele-count gate is applied per stratum downstream; `.fam` IID is the sample key used everywhere else; `.bim` IDs `chr:pos:ref:alt` | Module 6 `11_final_datasets` (`*_final_nomaf.*`) |
| `genotypes/cohort.hardcall.maf01.bed/.bim/.fam` | the same QC-complete set at cohort MAF ≥ 1 %, for anyone wanting a common-variant dataset without re-running QC | PLINK 1 trio; identical `.fam` to the nomaf set; a subset of its sites | Module 6 `11_final_datasets` (`*_final_maf01.*`) |
| `genotypes/chr{N}.phased.dose.vcf.gz` + `.csi` | FLARE input (phased GT); SeqVarGDS source; Tractor extraction | bgzipped VCF, hg38, `chr{N}` contigs, phased `GT` (`|`), `FORMAT/DS` present, INFO carries the per-service quality field and the MagicalRsq-X score as `INFO/MRSQX`; sample order identical to `.fam` | Module 6 per-chromosome filtered VCFs |
| `ancestry/ancestry_calls.tsv` | stratification, per-stratum MAC gate, the label on every result, the 21 continuous scores as covariates or figures | tab-separated, 24 columns exactly: `sample_id AncGroupID MajorGroup GD1 GD2 GD3 EA1 EA2 EA3 EA4 AF1 AF2 AF3 EU1 EU2 EU3 SA1 SA2 IC1 IC2 IC3 Pe Pf Pa`. `AncGroupID` is one of GrafAnc 1.0.0's 46 three-digit codes; `MajorGroup` is derived from its hundreds digit (1 AFR, 2 MEN, 3 EUR, 4 SAS, 5 EAS, 6 AMR, 7 OCN, 8 MIX) — GrafAnc emits no text group itself. One row per `.fam` IID, no duplicates. Made by `grafanc_to_calls.sh`, which also survives GrafAnc's multi-round output (later rounds are *appended* with a repeated header) and strips an `FID_` prefix from `Sample` if that is what the `.fam` matches | Module 7 `01_grafanc`, from `*_grafanc_results.txt` |
| `ancestry/grafanc_raw.txt` | audit only | GrafAnc's own 27-column output, untouched: `Sample #SNPs GD1-3 EA1-4 AF1-3 EU1-3 SA1-2 IC1-3 Pe Pf Pa RawPe RawPf RawPa AncGroupID` (tab-delimited with a header, inferred from the author's own reader) | Module 7 |
| `ancestry/strata.tsv` | the analysis unit for every downstream result | tab-separated, header `sample_id stratum`; one row per `.fam` IID; strata already bundled by the rule in `strata_rule.txt`, MIX retained as its own stratum | Module 7, `ASSIGN_STRATA` |
| `ancestry/strata_rule.txt` | so a stratum set is reproducible | the floor, the bundling rule applied, and the list of groups bundled or dropped | Module 7 |
| `structure/cohort.pcair.pcs.tsv` | the whole-cohort PCs for the pooled analysis and as covariates | header `sample_id PC1 … PCn`; one row per `.fam` IID | Module 6/7 PC-AiR, cohort scale |
| `structure/strata/<stratum>.pcair.pcs.tsv` | the within-stratum PCs for each stratified GENESIS run | one file per stratum in `strata.tsv`, same header, rows = that stratum's members | Module 6/7 PC-AiR, stratum scale |
| `structure/cohort.pcrelate.grm.rds`, `structure/cohort.pcrelate.kinship.rds` | the relatedness term in every null model | R matrices with **dimnames = `.fam` IID**; KING-robust divergence used throughout the PC-AiR / PC-Relate iteration | Module 6/7 PC-Relate |
| `local_ancestry/flare/chr{N}.anc.vcf.gz` + `.csi` | Tractor extraction (`extract_tracts_flare.py`) | FLARE output per chromosome, **indexed by upstream** (FLARE does not index); samples identical to `.fam` | Module 7 `FLARE_LAI` |
| `local_ancestry/flare_panels.txt` | provenance of the LAI reference | the exact FLARE `ref-panel=` sample map used: `sample<TAB>label` | Module 7 (copied from `params.flare_panel_file`) |
| `local_ancestry/ancestry_labels.txt` | the ancestry list the adapter is given (`--ancestries`) | one label per line, **alphabetical** — the order the Tractor-GENESIS adapter reads and reports in, so a column is always identifiable by name | derived from the panel file |
| `local_ancestry/flare_index_map.tsv` | the bridge from FLARE's numbering to those labels | header `anc_index label`; rows `0,1,2…` in the panel file's order of first appearance, which is how FLARE numbers `anc0…ancK-1`; the extraction step renames by this map so every downstream file is keyed by label | derived from the panel file |
| `phenotypes/phenotypes.tsv` | null model, Cox branch, covariates | tab-separated, header `FID IID sample_id time event <traits…> <covariates…>`; `sample_id == IID`; `time` numeric, `event` 0/1; **no free-text columns** (no `comments`, `cyto_text`, `death_cause`…, no `mrn`, no `dob` — derive age upstream) | the study, placed here by the harmonization run |
| `phenotypes/covariates.txt` | what the model adjusts for | one column name per line, every name present in `phenotypes.tsv` | the study |
| `qc/removed_duplicates.tsv` | answers "were duplicates removed upstream" | header `kept_iid removed_iid kinship reason`; may be empty but must exist; consistent with `# duplicates_removed` | Module 6 kinship step |
| `qc/sample_accounting.tsv` | n at every join | header `step n_in n_out n_removed reason` from first intake to this bundle | every Module |
| `qc/qc_assertions.tsv` | the upstream verdicts | `id verdict n message`; a bundle with no PASS rows is refused | `qc_assert()` |
| `run_record/params.tsv`, `run_record/trace.tsv` | provenance | the upstream run's resolved parameters and execution trace | Nextflow run |

## Two bundles, one contract

The bundle is delivered twice. **`--stage cohort`** (January, at M-COHORT) carries everything except
`local_ancestry/` — including `structure/` and `strata.tsv` — and is enough for the per-stratum and
pooled GENESIS runs. January is **first-pass** imputation; the second pass lands in the March re-run. **`--stage tractor`** (March, at M-COHORT-FINAL) is the same bundle regenerated on second-pass
imputation, plus the FLARE calls, the panel file, the alphabetical label list and the index map — what
the Tractor extraction consumes. The checker takes `--stage`; the default
is `tractor`, so an unqualified call demands everything.

**On GrafAnc's `Sample` column.** The documentation does not say whether it carries the PLINK IID or
`FID_IID`. Rule 1 holds either way because the normaliser maps it onto the `.fam` IID and refuses a row it
cannot place; confirm the behaviour on the first real run and record it in `strata_rule.txt`.

## Five rules

1. **One key.** The sample key is the `.fam` IID. `phenotypes.tsv` `sample_id`, every VCF header,
   `ancestry_calls.tsv`, and the FLARE VCFs must carry exactly that set, in the same order where order
   matters (VCFs). The checker fails on the first mismatch.
2. **Two hard-call sets, one sample set.** `nomaf` is the GWAS input; `maf01` is the same QC-complete
   cohort at MAF ≥ 1 % for anyone who wants a common-variant dataset without re-running QC. Same `.fam`,
   byte-exact. The checker refuses a bundle where they differ.
2b. **The bundle is complete or it is not a bundle.** `manifest.tsv` lists every file; the checker fails on
   the first listed path that is missing or zero bytes, and on any `.fam` IID absent from the phenotypes.
3. **Structure is computed once, upstream, on the whole cohort — and shipped.** PC-AiR and PC-Relate run
   in the harmonization pipeline at both scales (cohort and within stratum), because the stratified runs
   depend on quantities that need everyone together. The bundle carries the PCs at both scales and the
   GRM; downstream ingests them and re-derives one stratum only to *verify*, never to replace. Strata are
   assigned upstream too. Downstream computes only what is per-run: the MAC gate and the Tractor extraction.
4. **Quality filtering happened upstream, once per imputation, with MagicalRsq-X.** The bundle carries
   the score in `INFO/MRSQX` so the downstream can report it; it does not re-filter on raw `R2`/`DR2`.
5. **Nothing identifiable beyond the sample ID.** Study sample IDs are fine. Free-text clinical fields,
   MRNs and dates of birth are refused by the checker by column name.

## How to use it

```
# upstream, at the end of a harmonization run
make_handoff_manifest.sh handoff/<run_id> --run-id <run_id> --git-sha <sha> --server topmed --dups yes
# downstream, before anything runs
check_handoff_bundle.sh handoff/<run_id>          # exit 0 = accepted; non-zero names the first failure
```
