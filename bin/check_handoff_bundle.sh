#!/usr/bin/env bash
# Refuse a handoff bundle that is incomplete, mis-shaped, or mis-keyed. Exits non-zero naming the FIRST
# failure, so the message is always one thing to fix. Run downstream before anything else.
# usage: check_handoff_bundle.sh <bundle_dir> [--no-sha] [--chromosomes 1-22|21,22]
#   One bundle, complete: second-pass imputation, QC, ancestry, strata, structure, local ancestry, phenotypes.
#   --no-sha              skip the slow checksum pass
#   --chromosomes LIST    which chromosomes the bundle must carry (default 1-22); a module-ladder
#                         bundle built on two chromosomes is checked with --chromosomes 21,22
set -euo pipefail
B="${1:?bundle dir}"; shift; NOSHA=""
CHROM_SPEC="1-22"
while [[ $# -gt 0 ]]; do case "$1" in --no-sha) NOSHA="--no-sha"; shift;; --chromosomes) CHROM_SPEC="$2"; shift 2;; *) echo "unknown arg $1" >&2; exit 1;; esac; done
CHROMS=(); for part in ${CHROM_SPEC//,/ }; do if [[ "$part" =~ ^([0-9]+)-([0-9]+)$ ]]; then CHROMS+=($(seq "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}")); elif [[ "$part" =~ ^[0-9]+$ ]]; then CHROMS+=("$part"); else echo "bad --chromosomes value: $CHROM_SPEC" >&2; exit 1; fi; done
[[ ${#CHROMS[@]} -gt 0 ]] || { echo "empty --chromosomes" >&2; exit 1; }
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
for set in nomaf maf01; do for ext in bed bim fam; do need "genotypes/cohort.hardcall.${set}.${ext}"; done; done
for f in ancestry/ancestry_calls.tsv ancestry/grafanc_raw.txt ancestry/strata.tsv ancestry/strata_rule.txt \
         structure/cohort.pcair.pcs.tsv structure/cohort.pcrelate.grm.rds structure/cohort.pcrelate.kinship.rds \
         phenotypes/phenotypes.tsv phenotypes/covariates.txt \
         qc/sample_accounting.tsv qc/qc_assertions.tsv run_record/params.tsv run_record/trace.tsv; do need "$f"; done
[[ -e qc/removed_duplicates.tsv ]] || fail "qc/removed_duplicates.tsv must exist (may be header-only)"
for n in "${CHROMS[@]}"; do need "genotypes/chr${n}.phased.dose.vcf.gz"; need "genotypes/chr${n}.phased.dose.vcf.gz.csi"; done
need local_ancestry/flare_panels.txt; need local_ancestry/ancestry_labels.txt
for n in "${CHROMS[@]}"; do need "local_ancestry/flare/chr${n}.anc.vcf.gz"; need "local_ancestry/flare/chr${n}.anc.vcf.gz.csi"; done

# 4. both PLINK trios are internally consistent: bed bytes == 3 + ceil(n_samples/4) * n_variants,
#    and they carry the same samples: maf01 (the GWAS input, cohort MAF >= 1 %) is the same people as nomaf
#    and a subset of its sites
for set in nomaf maf01; do
  ns=$(wc -l < genotypes/cohort.hardcall.$set.fam); nv=$(wc -l < genotypes/cohort.hardcall.$set.bim)
  exp=$(( 3 + ((ns + 3) / 4) * nv )); got=$(stat -c %s genotypes/cohort.hardcall.$set.bed)
  [[ "$got" -eq "$exp" ]] || fail "cohort.hardcall.$set.bed is $got bytes; $ns samples x $nv variants implies $exp"
done
cmp -s <(awk '{print $2}' genotypes/cohort.hardcall.nomaf.fam) <(awk '{print $2}' genotypes/cohort.hardcall.maf01.fam) \
  || fail "nomaf and maf01 .fam files list different samples"
awk '{print $2}' genotypes/cohort.hardcall.maf01.bim | LC_ALL=C sort > .bim_maf01
awk '{print $2}' genotypes/cohort.hardcall.nomaf.bim | LC_ALL=C sort > .bim_nomaf
EXTRA=$(LC_ALL=C comm -23 .bim_maf01 .bim_nomaf | head -1)
[[ -z "$EXTRA" ]] || fail "maf01 carries a site absent from nomaf ($EXTRA): maf01 must be a subset of the unfiltered set"
NS=$(wc -l < genotypes/cohort.hardcall.maf01.fam); NV=$(wc -l < genotypes/cohort.hardcall.maf01.bim); NV_ALL=$(wc -l < genotypes/cohort.hardcall.nomaf.bim)

# 5. one key: .fam IID (the maf01 set is the GWAS input; nomaf carries the same .fam) == phenotypes sample_id
#    == ancestry_calls sample_id == every VCF header
awk '{print $2}' genotypes/cohort.hardcall.maf01.fam | LC_ALL=C sort > .fam_iid
H=$(head -1 phenotypes/phenotypes.tsv)
for c in FID IID sample_id mrd efs_time efs_event os_time os_event; do grep -qP "(^|\t)${c}(\t|$)" <<<"$H" || fail "phenotypes.tsv lacks column '$c' (the GWAS runs MRD, EFS and OS)"; done
# event columns are 0/1 or NA; time columns numeric or NA
for pair in efs os; do
  tc=$(awk -F'\t' -v n="${pair}_time" 'NR==1{for(i=1;i<=NF;i++) if($i==n) print i}' phenotypes/phenotypes.tsv)
  ec=$(awk -F'\t' -v n="${pair}_event" 'NR==1{for(i=1;i<=NF;i++) if($i==n) print i}' phenotypes/phenotypes.tsv)
  awk -F'\t' -v t="$tc" -v e="$ec" 'NR>1 && $e !~ /^(0|1|NA)$/ {exit 1}' phenotypes/phenotypes.tsv || fail "phenotypes.tsv ${pair}_event has a value other than 0, 1 or NA"
  awk -F'\t' -v t="$tc" 'NR>1 && $t !~ /^(NA|[0-9]+(\.[0-9]+)?)$/ {exit 1}' phenotypes/phenotypes.tsv || fail "phenotypes.tsv ${pair}_time has a non-numeric value"
done
for bad in mrn dob comments cyto_text oth_features death_cause relapse_site_oth; do
  grep -qiP "(^|\t)${bad}(\t|$)" <<<"$H" && fail "phenotypes.tsv carries forbidden column '$bad'"; done
SC=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="sample_id") print i}' phenotypes/phenotypes.tsv)
awk -F'\t' -v c="$SC" 'NR>1{print $c}' phenotypes/phenotypes.tsv | LC_ALL=C sort > .ph_id
cmp -s .fam_iid .ph_id || fail "phenotypes.tsv sample_id set differs from .fam IID set (first diff: $(comm -3 .fam_iid .ph_id | head -1))"
ANC_HDR=$'sample_id\tAncGroupID\tMajorGroup\tGD1\tGD2\tGD3\tEA1\tEA2\tEA3\tEA4\tAF1\tAF2\tAF3\tEU1\tEU2\tEU3\tSA1\tSA2\tIC1\tIC2\tIC3\tPe\tPf\tPa'
[[ "$(head -1 ancestry/ancestry_calls.tsv)" == "$ANC_HDR" ]] || fail "ancestry_calls.tsv header is not the 24-column GrafAnc-normalised header (sample_id AncGroupID MajorGroup GD1-3 EA1-4 AF1-3 EU1-3 SA1-2 IC1-3 Pe Pf Pa)"
awk -F'\t' 'NR>1{print $1}' ancestry/ancestry_calls.tsv | LC_ALL=C sort > .anc_id
[[ "$(sort -u .anc_id | wc -l)" -eq "$(wc -l < .anc_id)" ]] || fail "ancestry_calls.tsv has a duplicated sample_id (GrafAnc appends later rounds; the normaliser must not)"
cmp -s .fam_iid .anc_id || fail "ancestry_calls.tsv sample set differs from .fam IID set"
# the 46 AncGroupID codes GrafAnc 1.0.0 can emit, and the hundreds digit -> MajorGroup rule
GRAF_CODES=" 100 101 102 103 104 105 106 107 108 200 201 202 203 300 301 302 303 304 305 306 307 308 400 401 402 403 404 405 500 501 502 503 504 505 506 507 508 509 510 511 600 601 602 603 700 800 "
awk -F'\t' -v codes="$GRAF_CODES" 'BEGIN{split("AFR MEN EUR SAS EAS AMR OCN MIX",G," ")}
  NR>1 { if (index(codes, " " $2 " ")==0) {print "code " $2 " on " $1; exit 1}
         if ($3 != G[int($2/100)]) {print "MajorGroup " $3 " does not match AncGroupID " $2 " on " $1; exit 1} }' ancestry/ancestry_calls.tsv \
  || fail "ancestry_calls.tsv: $(awk -F'\t' -v codes="$GRAF_CODES" 'BEGIN{split("AFR MEN EUR SAS EAS AMR OCN MIX",G," ")} NR>1 { if (index(codes, " " $2 " ")==0) {print "AncGroupID " $2 " is not a GrafAnc code (sample " $1 ")"; exit} if ($3 != G[int($2/100)]) {print "MajorGroup " $3 " does not match AncGroupID " $2 " (sample " $1 ")"; exit} }' ancestry/ancestry_calls.tsv)"
head -1 ancestry/grafanc_raw.txt | grep -qP '^Sample\t' || fail "grafanc_raw.txt does not start with GrafAnc's header (first column Sample)"
if command -v bcftools >/dev/null; then
  for n in "${CHROMS[0]}" "${CHROMS[-1]}"; do   # first and last are enough to catch a wrong sample set; full check is a run-time job
    bcftools query -l "genotypes/chr${n}.phased.dose.vcf.gz" | LC_ALL=C sort | cmp -s - .fam_iid || fail "chr${n}.phased.dose.vcf.gz samples differ from .fam IID set"
    bcftools view -h "genotypes/chr${n}.phased.dose.vcf.gz" | grep -q '^##FORMAT=<ID=DS' || fail "chr${n}.phased.dose.vcf.gz has no FORMAT/DS"
    bcftools view -h "genotypes/chr${n}.phased.dose.vcf.gz" | grep -q '^##INFO=<ID=MRSQX' || fail "chr${n}.phased.dose.vcf.gz has no INFO/MRSQX (MagicalRsq-X score)"
  done
fi

# 5b. strata and structure: every sample has a stratum; cohort PCs and one PC file per stratum; GRM dimnames
head -1 ancestry/strata.tsv | grep -qP '^sample_id\tstratum' || fail "strata.tsv header must start sample_id<TAB>stratum"
awk -F'\t' 'NR>1{print $1}' ancestry/strata.tsv | LC_ALL=C sort > .st_id
cmp -s .fam_iid .st_id || fail "strata.tsv sample set differs from .fam IID set"
head -1 structure/cohort.pcair.pcs.tsv | grep -qP '^sample_id\tPC1(\t|$)' || fail "cohort.pcair.pcs.tsv header must start sample_id<TAB>PC1"
awk -F'\t' 'NR>1{print $1}' structure/cohort.pcair.pcs.tsv | LC_ALL=C sort > .pc_id
cmp -s .fam_iid .pc_id || fail "cohort.pcair.pcs.tsv sample set differs from .fam IID set"
awk -F'\t' 'NR>1{print $2}' ancestry/strata.tsv | LC_ALL=C sort -u | while read -r S; do
  [[ -s "structure/strata/${S}.pcair.pcs.tsv" ]] || fail "no per-stratum PCs for stratum '$S' (structure/strata/${S}.pcair.pcs.tsv)"
  awk -F'\t' -v s="$S" 'NR>1 && $2==s {print $1}' ancestry/strata.tsv | LC_ALL=C sort > .st_members
  awk -F'\t' 'NR>1{print $1}' "structure/strata/${S}.pcair.pcs.tsv" | LC_ALL=C sort | cmp -s - .st_members \
    || fail "structure/strata/${S}.pcair.pcs.tsv samples differ from stratum '$S' members in strata.tsv"
done
if command -v Rscript >/dev/null; then
  Rscript --vanilla -e 'ids<-readLines(".fam_iid"); for(f in c("structure/cohort.pcrelate.grm.rds","structure/cohort.pcrelate.kinship.rds")){m<-readRDS(f); if(is.null(rownames(m))||!identical(sort(rownames(m)),sort(ids))) {cat(f,"\n"); quit(status=3)}}' \
    || fail "a structure .rds has no dimnames or its dimnames differ from the .fam IID set"
fi

# 6. covariates named exist; ancestry labels are ALPHABETICAL and the index map is explicit; duplicates answer is consistent
while read -r cv; do [[ -z "$cv" ]] && continue; grep -qP "(^|\t)${cv}(\t|$)" <<<"$H" || fail "covariate '$cv' not in phenotypes.tsv"; done < phenotypes/covariates.txt
if true; then
  need local_ancestry/flare_index_map.tsv
  awk '{print $2}' local_ancestry/flare_panels.txt | LC_ALL=C sort -u > .panel_labels
  LC_ALL=C sort -c local_ancestry/ancestry_labels.txt 2>/dev/null || fail "ancestry_labels.txt is not in alphabetical order"
  cmp -s .panel_labels local_ancestry/ancestry_labels.txt || fail "ancestry_labels.txt is not the sorted set of labels in flare_panels.txt"
  head -1 local_ancestry/flare_index_map.tsv | grep -qP '^anc_index\tlabel' || fail "flare_index_map.tsv header must be anc_index<TAB>label"
  awk '{print $2}' local_ancestry/flare_panels.txt | awk '!seen[$0]++' > .first_seen
  awk -F'\t' 'NR>1{print $2}' local_ancestry/flare_index_map.tsv | cmp -s - .first_seen \
    || fail "flare_index_map.tsv labels are not in the panel file's order of first appearance (that is how FLARE numbers anc0..ancK-1)"
  awk -F'\t' 'NR>1 && $1 != NR-2 {exit 1}' local_ancestry/flare_index_map.tsv || fail "flare_index_map.tsv anc_index must be 0,1,2,... in row order"
fi
NR_DUP=$(($(wc -l < qc/removed_duplicates.tsv) - 1))
[[ "$DUPS" == "yes" || "$NR_DUP" -eq 0 ]] || fail "duplicates_removed=no but removed_duplicates.tsv lists $NR_DUP removals"
grep -qP '\tPASS\t' qc/qc_assertions.tsv || fail "qc_assertions.tsv has no PASS rows — 'did not run' is not 'passed'"
rm -f .fam_iid .ph_id .anc_id .st_id .pc_id .st_members .panel_labels .first_seen .bim_maf01 .bim_nomaf
echo "HANDOFF ACCEPTED: $B  ($NS samples, $NV GWAS-input variants at MAF >= 1 %, $NV_ALL unfiltered, chromosomes $CHROM_SPEC, server=$SERVER, duplicates_removed=$DUPS)"
