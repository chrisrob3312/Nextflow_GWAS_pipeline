# Handoff contract — harmonization pipeline → GWAS pipeline

**One directory, one manifest, no guessing.** The harmonization pipeline writes everything the GWAS
pipeline needs into `handoff/<run_id>/`, lists it in `manifest.tsv`, and `check_handoff_bundle.sh`
refuses the bundle if anything is missing, empty, mis-shaped, or mis-keyed. Every file below is named
exactly as written; the chromosome placeholder `{N}` is 1–22 with no `chr` prefix in the file name
(contigs *inside* VCFs are `chr`-prefixed, hg38).

| File | Role downstream | Format and required content | Produced by |
|---|---|---|---|
| `manifest.tsv` | the contract itself | tab-separated, header `path sha256 bytes role`; one row per file below; header lines `# run_id`, `# upstream_git_sha`, `# imputation_server` (topmed \| michigan \| anvil), `# duplicates_removed` (yes \| no), `# build` (hg38) | `make_handoff_manifest.sh` |
| `genotypes/cohort.hardcall.bed/.bim/.fam` | PC-AiR / PC-Relate, GrafAnc audit, stratum keep-files | PLINK 1 trio; hg38; the **nomaf** set (the 1 % MAF filter is descriptive and is applied downstream only as a report); `.fam` IID is the sample key used everywhere else; `.bim` IDs `chr:pos:ref:alt` | Module 6 `11_final_datasets` (`*_final_nomaf.*`) |
| `genotypes/chr{N}.phased.dose.vcf.gz` + `.csi` | FLARE input (phased GT); SeqVarGDS source; Tractor extraction | bgzipped VCF, hg38, `chr{N}` contigs, phased `GT` (`|`), `FORMAT/DS` present, INFO carries the per-service quality field and the MagicalRsq-X score as `INFO/MRSQX`; sample order identical to `.fam` | Module 6 per-chromosome filtered VCFs |
| `ancestry/ancestry_calls.tsv` | stratification (`ANCESTRY_STRATIFY`), per-stratum MAC gate, labels on every result | tab-separated, header `sample_id AncGroupID MajorGroup GD1 GD2 GD3`; `MajorGroup` = floor(AncGroupID/100) as the label AFR/MEN/EUR/SAS/EAS/AMR/OCN/MIX; one row per `.fam` IID, no extras | Module 7 `01_grafanc`, normalised from `*_grafanc_results.txt` |
| `ancestry/grafanc_raw.txt` | audit only | GrafAnc's own output, untouched | Module 7 |
| `local_ancestry/flare/chr{N}.anc.vcf.gz` + `.csi` | Tractor extraction (`extract_tracts_flare.py`) | FLARE output per chromosome, **indexed by upstream** (FLARE does not index); samples identical to `.fam` | Module 7 `FLARE_LAI` |
| `local_ancestry/flare_panels.txt` | ancestry index → label order | the exact FLARE `ref-panel=` sample map used: `sample<TAB>label`; the **order of first appearance** of labels defines `anc0, anc1, …` | Module 7 (copied from `params.flare_panel_file`) |
| `local_ancestry/ancestry_labels.txt` | the same order, explicitly | one label per line, in `anc0…ancK-1` order | derived from the panel file |
| `phenotypes/phenotypes.tsv` | null model, Cox branch, covariates | tab-separated, header `FID IID sample_id time event <traits…> <covariates…>`; `sample_id == IID`; `time` numeric, `event` 0/1; **no free-text columns** (no `comments`, `cyto_text`, `death_cause`…, no `mrn`, no `dob` — derive age upstream) | the study, placed here by the harmonization run |
| `phenotypes/covariates.txt` | what the model adjusts for | one column name per line, every name present in `phenotypes.tsv` | the study |
| `qc/removed_duplicates.tsv` | answers "were duplicates removed upstream" | header `kept_iid removed_iid kinship reason`; may be empty but must exist; consistent with `# duplicates_removed` | Module 6 kinship step |
| `qc/sample_accounting.tsv` | n at every join | header `step n_in n_out n_removed reason` from first intake to this bundle | every Module |
| `qc/qc_assertions.tsv` | the upstream verdicts | `id verdict n message`; a bundle with no PASS rows is refused | `qc_assert()` |
| `run_record/params.tsv`, `run_record/trace.tsv` | provenance | the upstream run's resolved parameters and execution trace | Nextflow run |

## Two bundles, one contract

The bundle is delivered twice. **`--stage cohort`** (January, at M-COHORT) carries everything except
`local_ancestry/`, and is enough for PC-AiR, PC-Relate, stratification and the per-stratum GENESIS
runs. **`--stage tractor`** (March, at M-COHORT-FINAL) adds the FLARE calls, the panel file and the
label order, and is what the Tractor extraction consumes. The checker takes `--stage`; the default
is `tractor`, so an unqualified call demands everything.

## Five rules

1. **One key.** The sample key is the `.fam` IID. `phenotypes.tsv` `sample_id`, every VCF header,
   `ancestry_calls.tsv`, and the FLARE VCFs must carry exactly that set, in the same order where order
   matters (VCFs). The checker fails on the first mismatch.
2. **The bundle is complete or it is not a bundle.** `manifest.tsv` lists every file; the checker fails on
   the first listed path that is missing or zero bytes, and on any `.fam` IID absent from the phenotypes.
3. **Nothing is computed twice.** Upstream delivers hard calls, dosage VCFs, ancestry calls, local
   ancestry, and the duplicate decision. Downstream computes PC-AiR, PC-Relate, strata, the MAC gate, and
   the Tractor extraction. If upstream also writes PCs or a GRM, they are *not* in this bundle.
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
