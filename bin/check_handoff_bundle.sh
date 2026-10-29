#!/usr/bin/env bash
# Refuse a handoff bundle that is incomplete, mis-shaped, or mis-keyed. Exits non-zero naming the FIRST
# failure, so the message is always one thing to fix. Run downstream before anything else.
# usage: check_handoff_bundle.sh <bundle_dir> [--stage cohort|tractor] [--no-sha]
#   --stage cohort   the January bundle: genotypes, ancestry calls, phenotypes, QC, run record (no local ancestry)
#   --stage tractor  the March bundle: everything above PLUS local_ancestry/ (default)
#   --no-sha         skip the slow checksum pass
set -euo pipefail
B="${1:?bundle dir}"; shift; STAGE=tractor; NOSHA=""
while [[ $# -gt 0 ]]; do case "$1" in --stage) STAGE="$2"; shift 2;; --no-sha) NOSHA="--no-sha"; shift;; *) echo "unknown arg $1" >&2; exit 1;; esac; done
[[ "$STAGE" =~ ^(cohort|tractor)$ ]] || { echo "--stage must be cohort|tractor" >&2; exit 1; }
fail() { echo "HANDOFF REFUSED: $*" >&2; exit 2; }
need() { [[ -s "$B/$1" ]] || fail "$1 is missing or empty"; }
cd "$B" 2>/dev/null || fail "$B is not a directory"

# 1. the manifest and its header fields
need manifest.tsv
for k in run_id upstream_git_sha imputation_server duplicates_removed build; do
  grep -qP "^# ${k}\t\S" manifest.tsv || fail "manifest.tsv lacks header '# ${k}'"
done
SERVER=$(grep -P '^# imputation_server\t' manifest.tsv | cut -f2); DUPS=$(grep -P '^# duplicates_removed\t' manifest.tsv | cut -f2)
[[ "$SERVER" =~ ^(topmed|michigan|anvil)$ ]] || fail "imputation_server '$SERVER' is not topmed|michigan|anvil"
[[ "$DUPS" =~ ^(yes|no)$ ]] || fail "duplicates_removed '$DUPS' is not yes|no"

# 2. every listed file exists and is non-empty; checksum unless --no-sha
grep -v '^#' manifest.tsv | tail -n +2 | while IFS=$'\t' read -r p sha bytes role; do
  [[ -s "$p" ]] || fail "listed in manifest but missing or empty: $p"
  [[ "$NOSHA" == "--no-sha" ]] || [[ "$(sha256sum "$p" | cut -d' ' -f1)" == "$sha" ]] || fail "checksum mismatch: $p"
done

# 3. the fixed file set
for f in genotypes/cohort.hardcall.bed genotypes/cohort.hardcall.bim genotypes/cohort.hardcall.fam \
         ancestry/ancestry_calls.tsv ancestry/grafanc_raw.txt phenotypes/phenotypes.tsv phenotypes/covariates.txt \
         qc/sample_accounting.tsv qc/qc_assertions.tsv run_record/params.tsv run_record/trace.tsv; do need "$f"; done
[[ -e qc/removed_duplicates.tsv ]] || fail "qc/removed_duplicates.tsv must exist (may be header-only)"
for n in $(seq 1 22); do need "genotypes/chr${n}.phased.dose.vcf.gz"; need "genotypes/chr${n}.phased.dose.vcf.gz.csi"; done
if [[ "$STAGE" == tractor ]]; then
  need local_ancestry/flare_panels.txt; need local_ancestry/ancestry_labels.txt
  for n in $(seq 1 22); do need "local_ancestry/flare/chr${n}.anc.vcf.gz"; need "local_ancestry/flare/chr${n}.anc.vcf.gz.csi"; done
fi

# 4. the PLINK trio is internally consistent: bed bytes == 3 + ceil(n_samples/4) * n_variants
NS=$(wc -l < genotypes/cohort.hardcall.fam); NV=$(wc -l < genotypes/cohort.hardcall.bim)
EXP=$(( 3 + ((NS + 3) / 4) * NV )); GOT=$(stat -c %s genotypes/cohort.hardcall.bed)
[[ "$GOT" -eq "$EXP" ]] || fail "cohort.hardcall.bed is $GOT bytes; $NS samples x $NV variants implies $EXP"

# 5. one key: .fam IID == phenotypes sample_id == ancestry_calls sample_id == every VCF header
awk '{print $2}' genotypes/cohort.hardcall.fam | LC_ALL=C sort > .fam_iid
H=$(head -1 phenotypes/phenotypes.tsv)
for c in FID IID sample_id time event; do grep -qP "(^|\t)${c}(\t|$)" <<<"$H" || fail "phenotypes.tsv lacks column '$c'"; done
for bad in mrn dob comments cyto_text oth_features death_cause relapse_site_oth; do
  grep -qiP "(^|\t)${bad}(\t|$)" <<<"$H" && fail "phenotypes.tsv carries forbidden column '$bad'"; done
SC=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="sample_id") print i}' phenotypes/phenotypes.tsv)
awk -F'\t' -v c="$SC" 'NR>1{print $c}' phenotypes/phenotypes.tsv | LC_ALL=C sort > .ph_id
cmp -s .fam_iid .ph_id || fail "phenotypes.tsv sample_id set differs from .fam IID set (first diff: $(comm -3 .fam_iid .ph_id | head -1))"
head -1 ancestry/ancestry_calls.tsv | grep -qP '^sample_id\tAncGroupID\tMajorGroup\tGD1\tGD2\tGD3' || fail "ancestry_calls.tsv header must start sample_id AncGroupID MajorGroup GD1 GD2 GD3"
awk -F'\t' 'NR>1{print $1}' ancestry/ancestry_calls.tsv | LC_ALL=C sort > .anc_id
cmp -s .fam_iid .anc_id || fail "ancestry_calls.tsv sample set differs from .fam IID set"
awk -F'\t' 'NR>1 && $3 !~ /^(AFR|MEN|EUR|SAS|EAS|AMR|OCN|MIX)$/ {print; exit 1}' ancestry/ancestry_calls.tsv >/dev/null \
  || fail "ancestry_calls.tsv MajorGroup has a value outside AFR|MEN|EUR|SAS|EAS|AMR|OCN|MIX"
if command -v bcftools >/dev/null; then
  for n in 1 22; do   # first and last are enough to catch a wrong sample set; full check is a run-time job
    bcftools query -l "genotypes/chr${n}.phased.dose.vcf.gz" | LC_ALL=C sort | cmp -s - .fam_iid || fail "chr${n}.phased.dose.vcf.gz samples differ from .fam IID set"
    bcftools view -h "genotypes/chr${n}.phased.dose.vcf.gz" | grep -q '^##FORMAT=<ID=DS' || fail "chr${n}.phased.dose.vcf.gz has no FORMAT/DS"
    bcftools view -h "genotypes/chr${n}.phased.dose.vcf.gz" | grep -q '^##INFO=<ID=MRSQX' || fail "chr${n}.phased.dose.vcf.gz has no INFO/MRSQX (MagicalRsq-X score)"
  done
fi

# 6. covariates named exist; ancestry labels match the panel order; duplicates answer is consistent
while read -r cv; do [[ -z "$cv" ]] && continue; grep -qP "(^|\t)${cv}(\t|$)" <<<"$H" || fail "covariate '$cv' not in phenotypes.tsv"; done < phenotypes/covariates.txt
if [[ "$STAGE" == tractor ]]; then
  awk '{print $2}' local_ancestry/flare_panels.txt | awk '!seen[$0]++' > .panel_order
  cmp -s .panel_order local_ancestry/ancestry_labels.txt || fail "ancestry_labels.txt is not the order of first appearance in flare_panels.txt"
fi
NR_DUP=$(($(wc -l < qc/removed_duplicates.tsv) - 1))
[[ "$DUPS" == "yes" || "$NR_DUP" -eq 0 ]] || fail "duplicates_removed=no but removed_duplicates.tsv lists $NR_DUP removals"
grep -qP '\tPASS\t' qc/qc_assertions.tsv || fail "qc_assertions.tsv has no PASS rows — 'did not run' is not 'passed'"
rm -f .fam_iid .ph_id .anc_id .panel_order
echo "HANDOFF ACCEPTED ($STAGE stage): $B  ($NS samples, $NV hardcall variants, server=$SERVER, duplicates_removed=$DUPS)"
