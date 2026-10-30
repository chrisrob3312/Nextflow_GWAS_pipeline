#!/usr/bin/env bash
# Build a small, internally consistent stub handoff bundle for exercising bin/check_handoff_bundle.sh.
# usage: make_stub_bundle.sh <out_dir> [--chromosomes 1-22|21,22]
# Nothing in it is real data: 8 samples, 6 unfiltered sites of which 4 pass MAF >= 1 %. The dosage and
# FLARE VCFs are placeholder bytes (bcftools-level checks only run where bcftools is installed), and the
# .rds files are placeholders for the same reason. What the stub DOES exercise is every structural rule in
# the contract: the one key, the two hard-call sets, the subset rule, the 24-column ancestry table, strata
# and per-stratum PCs, the alphabetical label list and the panel-order index map, the phenotype header
# for MRD/EFS/OS, the covariate list, the QC verdicts and the manifest.
set -euo pipefail
OUT="${1:?out_dir}"; shift || true
CHROM_SPEC="1-22"
while [[ $# -gt 0 ]]; do case "$1" in --chromosomes) CHROM_SPEC="$2"; shift 2;; *) echo "unknown arg $1" >&2; exit 1;; esac; done
CHROMS=(); for part in ${CHROM_SPEC//,/ }; do if [[ "$part" =~ ^([0-9]+)-([0-9]+)$ ]]; then CHROMS+=($(seq "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}")); elif [[ "$part" =~ ^[0-9]+$ ]]; then CHROMS+=("$part"); else echo "bad --chromosomes value: $part" >&2; exit 1; fi; done
rm -rf "$OUT"; mkdir -p "$OUT"/{genotypes,ancestry,structure/strata,local_ancestry/flare,phenotypes,qc,run_record}
cd "$OUT"

SAMPLES=(S001 S002 S003 S004 S005 S006 S007 S008)
# .fam: FID IID father mother sex pheno — IID is the key
for s in "${SAMPLES[@]}"; do printf 'F_%s\t%s\t0\t0\t1\t-9\n' "$s" "$s"; done > genotypes/cohort.hardcall.nomaf.fam
cp genotypes/cohort.hardcall.nomaf.fam genotypes/cohort.hardcall.maf01.fam
# .bim: chr id cm pos a1 a2 ; nomaf has 6 sites, maf01 the 4 common ones (a subset)
printf '1\tchr1:10000:A:G\t0\t10000\tG\tA\n1\tchr1:20000:C:T\t0\t20000\tT\tC\n1\tchr1:30000:G:A\t0\t30000\tA\tG\n2\tchr2:10000:T:C\t0\t10000\tC\tT\n2\tchr2:20000:A:C\t0\t20000\tC\tA\n2\tchr2:30000:G:T\t0\t30000\tT\tG\n' > genotypes/cohort.hardcall.nomaf.bim
grep -v -e 'chr1:30000' -e 'chr2:30000' genotypes/cohort.hardcall.nomaf.bim > genotypes/cohort.hardcall.maf01.bim
# .bed: 3 magic bytes + ceil(8/4)=2 bytes per variant
mkbed() { local nv=$1 out=$2; { printf '\x6c\x1b\x01'; head -c $((2 * nv)) /dev/zero; } > "$out"; }
mkbed 6 genotypes/cohort.hardcall.nomaf.bed; mkbed 4 genotypes/cohort.hardcall.maf01.bed
for n in "${CHROMS[@]}"; do printf 'placeholder dosage vcf chr%s\n' "$n" | gzip -c > "genotypes/chr${n}.phased.dose.vcf.gz"; printf 'csi' > "genotypes/chr${n}.phased.dose.vcf.gz.csi"; done

# ancestry: 24 columns; codes from GrafAnc 1.0.0; MajorGroup = hundreds digit
{ printf 'sample_id\tAncGroupID\tMajorGroup'; for c in GD1 GD2 GD3 EA1 EA2 EA3 EA4 AF1 AF2 AF3 EU1 EU2 EU3 SA1 SA2 IC1 IC2 IC3 Pe Pf Pa; do printf '\t%s' "$c"; done; printf '\n'
  codes=(300 300 301 100 101 600 600 800); groups=(EUR EUR EUR AFR AFR AMR AMR MIX)
  for i in "${!SAMPLES[@]}"; do printf '%s\t%s\t%s' "${SAMPLES[$i]}" "${codes[$i]}" "${groups[$i]}"; for k in $(seq 1 21); do printf '\t0.%02d' "$((i + k))"; done; printf '\n'; done; } > ancestry/ancestry_calls.tsv
{ printf 'Sample\t#SNPs\tGD1\tGD2\tGD3\tEA1\tEA2\tEA3\tEA4\tAF1\tAF2\tAF3\tEU1\tEU2\tEU3\tSA1\tSA2\tIC1\tIC2\tIC3\tPe\tPf\tPa\tRawPe\tRawPf\tRawPa\tAncGroupID\n'
  for s in "${SAMPLES[@]}"; do printf '%s\t100000' "$s"; for k in $(seq 1 24); do printf '\t0.1'; done; printf '\t300\n'; done; } > ancestry/grafanc_raw.txt
{ printf 'sample_id\tstratum\n'; strata=(EUR EUR EUR AFR AFR AMR AMR MIX); for i in "${!SAMPLES[@]}"; do printf '%s\t%s\n' "${SAMPLES[$i]}" "${strata[$i]}"; done; } > ancestry/strata.tsv
printf 'MajorGroup from GrafAnc AncGroupID hundreds digit; floor 50 per stratum with bundling; MIX retained (stub)\n' > ancestry/strata_rule.txt

# structure: cohort PCs, per-stratum PCs, placeholder rds
{ printf 'sample_id\tPC1\tPC2\tPC3\n'; for s in "${SAMPLES[@]}"; do printf '%s\t0.01\t0.02\t0.03\n' "$s"; done; } > structure/cohort.pcair.pcs.tsv
for S in $(awk -F'\t' 'NR>1{print $2}' ancestry/strata.tsv | sort -u); do
  { printf 'sample_id\tPC1\tPC2\n'; awk -F'\t' -v s="$S" 'NR>1 && $2==s {printf "%s\t0.01\t0.02\n", $1}' ancestry/strata.tsv; } > "structure/strata/${S}.pcair.pcs.tsv"; done
printf 'placeholder rds' > structure/cohort.pcrelate.grm.rds; printf 'placeholder rds' > structure/cohort.pcrelate.kinship.rds

# local ancestry: panel file in first-appearance order EUR, AFR, AMR -> index map 0,1,2; labels sorted
printf 'REF001\tEUR\nREF002\tAFR\nREF003\tEUR\nREF004\tAMR\nREF005\tAFR\n' > local_ancestry/flare_panels.txt
printf 'AFR\nAMR\nEUR\n' > local_ancestry/ancestry_labels.txt
printf 'anc_index\tlabel\n0\tEUR\n1\tAFR\n2\tAMR\n' > local_ancestry/flare_index_map.tsv
for n in "${CHROMS[@]}"; do printf 'placeholder flare vcf chr%s\n' "$n" | gzip -c > "local_ancestry/flare/chr${n}.anc.vcf.gz"; printf 'csi' > "local_ancestry/flare/chr${n}.anc.vcf.gz.csi"; done

# phenotypes: MRD (binary), EFS and OS (time, event), four covariates
{ printf 'FID\tIID\tsample_id\tmrd\tefs_time\tefs_event\tos_time\tos_event\tprotocol\tsex\tage_at_diagnosis\tcd19_pos_pct\tstorage_days\n'
  i=0; protos=(P1 P2 P3); for s in "${SAMPLES[@]}"; do printf 'F_%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\n' "$s" "$s" "$s" $((i % 2)) "$((300 + i * 10)).5" $((i % 2)) "$((400 + i * 10))" $(( (i + 1) % 2 )) "${protos[$((i % 3))]}" $((1 + i % 2)) $((5 + i)) "9$i.5" $((100 + i)); i=$((i + 1)); done; } > phenotypes/phenotypes.tsv
printf 'sex\nage_at_diagnosis\ncd19_pos_pct\nstorage_days\n' > phenotypes/covariates.txt

# qc and run record
printf 'kept_iid\tremoved_iid\tkinship\treason\n' > qc/removed_duplicates.tsv   # header-only: nothing removed
printf 'step\tn_in\tn_out\tn_removed\treason\n01_intake\t8\t8\t0\t-\n11_final_datasets\t8\t8\t0\t-\n' > qc/sample_accounting.tsv
printf 'id\tverdict\tn\tmessage\nM6_callrate_95\tPASS\t8\tall samples >= 0.95\nM6_mrsqx_applied\tPASS\t6\tMagicalRsq-X applied after each pass\n' > qc/qc_assertions.tsv
printf 'param\tvalue\nrun_lai\ttrue\nlai_methods\tflare\nimputation_quality_method\tmagicalrsqx\n' > run_record/params.tsv
printf 'task_id\tname\tstatus\n1\tMODULE1:harmonize\tCOMPLETED\n' > run_record/trace.tsv

# manifest in the documented format
role() { case "$1" in
  genotypes/cohort.hardcall.nomaf.*) echo hardcall_nomaf;; genotypes/cohort.hardcall.maf01.*) echo hardcall_maf01;; genotypes/chr*.phased.dose.vcf.gz) echo dosage_vcf;;
  genotypes/chr*.phased.dose.vcf.gz.csi) echo dosage_vcf_index;; ancestry/ancestry_calls.tsv) echo ancestry_calls;; ancestry/strata.tsv) echo strata;; ancestry/strata_rule.txt) echo strata_rule;;
  structure/strata/*.pcair.pcs.tsv) echo stratum_pcs;; structure/cohort.pcair.pcs.tsv) echo cohort_pcs;; structure/cohort.pcrelate.*.rds) echo structure_rds;;
  ancestry/grafanc_raw.txt) echo grafanc_raw;; local_ancestry/flare/chr*.anc.vcf.gz) echo flare_anc;;
  local_ancestry/flare/chr*.anc.vcf.gz.csi) echo flare_anc_index;; local_ancestry/flare_panels.txt) echo flare_panels;;
  local_ancestry/ancestry_labels.txt) echo ancestry_labels;; local_ancestry/flare_index_map.tsv) echo flare_index_map;; phenotypes/phenotypes.tsv) echo phenotypes;;
  phenotypes/covariates.txt) echo covariates;; qc/removed_duplicates.tsv) echo removed_duplicates;;
  qc/sample_accounting.tsv) echo sample_accounting;; qc/qc_assertions.tsv) echo qc_assertions;;
  run_record/*) echo run_record;; *) echo other;; esac; }
{ printf '# run_id\tstub-run\n# upstream_git_sha\t0000000\n# imputation_server\ttopmed\n# duplicates_removed\tno\n# build\thg38\n'
  printf 'path\tsha256\tbytes\trole\n'
  find genotypes ancestry structure local_ancestry phenotypes qc run_record -type f | LC_ALL=C sort | while read -r p; do
    printf '%s\t%s\t%s\t%s\n' "$p" "$(sha256sum "$p" | cut -d' ' -f1)" "$(stat -c %s "$p")" "$(role "$p")"; done; } > manifest.tsv
echo "stub bundle written to $OUT (chromosomes $CHROM_SPEC)"
