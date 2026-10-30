#!/usr/bin/env bash
# Refuse a handoff bundle that is incomplete, mis-shaped, or mis-keyed. Exits non-zero naming the FIRST
# failure, so the message is always one thing to fix. Run downstream before anything else.
# usage: check_handoff_bundle.sh <bundle_dir> [--no-sha] [--chromosomes 1-22|21,22] [--allow-missing-tools]
#   One bundle, complete: second-pass imputation, QC, ancestry, strata, structure, local ancestry, phenotypes.
#   --no-sha              skip the slow checksum pass
#   --chromosomes LIST    which chromosomes the bundle must carry (default 1-22); a module-ladder
#                         bundle built on two chromosomes is checked with --chromosomes 21,22
#   --allow-missing-tools accept a bundle even when bcftools / Rscript are absent: their checks are then
#                         SKIPPED and named in the ACCEPTED line. Without it a missing tool is a refusal,
#                         because "did not run" must never read as "passed".
set -euo pipefail
B="${1:?bundle dir}"; shift; NOSHA=""
CHROM_SPEC="1-22"; ALLOW_MISSING=0
while [[ $# -gt 0 ]]; do case "$1" in
  --no-sha) NOSHA="--no-sha"; shift;;
  --chromosomes) [[ $# -ge 2 ]] || { echo "--chromosomes needs a value (e.g. 1-22 or 21,22)" >&2; exit 1; }; CHROM_SPEC="$2"; shift 2;;
  --allow-missing-tools) ALLOW_MISSING=1; shift;;
  *) echo "unknown arg $1" >&2; exit 1;; esac; done
CHROMS=(); for part in ${CHROM_SPEC//,/ }; do
  if [[ "$part" =~ ^([0-9]+)-([0-9]+)$ ]]; then CHROMS+=($(seq "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"))
  elif [[ "$part" =~ ^[0-9]+$ ]]; then CHROMS+=("$part")
  else echo "bad --chromosomes value: $CHROM_SPEC" >&2; exit 1; fi; done
[[ ${#CHROMS[@]} -gt 0 ]] || { echo "empty --chromosomes" >&2; exit 1; }
for n in "${CHROMS[@]}"; do [[ "$n" -ge 1 && "$n" -le 22 ]] || { echo "--chromosomes: $n is outside 1-22 (the contract carries autosomes only)" >&2; exit 1; }; done
fail() { echo "HANDOFF REFUSED: $*" >&2; exit 2; }
[[ -d "$B" ]] || fail "$B is not a directory"
B="$(cd "$B" && pwd)"                       # absolute, so a relative bundle path works after the cd below
T="$(mktemp -d)"; trap 'rm -rf "$T"' EXIT  # scratch lives outside the bundle; a refused or read-only bundle stays clean
cd "$B"
: > "$T/required"
need() { echo "$1" >> "$T/required"; [[ -s "$1" ]] || fail "$1 is missing or empty"; }
SKIPPED=()
hdr() { awk -F'\t' -v k="# $1" '$1==k {print $2; exit}' manifest.tsv; }

# 1. the manifest and its header fields (the manifest is not listed in itself)
[[ -s manifest.tsv ]] || fail "manifest.tsv is missing or empty"
for f in run_id upstream_git_sha imputation_server duplicates_removed build; do
  [[ -n "$(hdr "$f")" ]] || fail "manifest.tsv lacks header field '# $f'"; done
SERVER=$(hdr imputation_server); DUPS=$(hdr duplicates_removed)
[[ "$SERVER" =~ ^(topmed|michigan|anvil)$ ]] || fail "imputation_server '$SERVER' is not topmed|michigan|anvil"
[[ "$DUPS" =~ ^(yes|no)$ ]] || fail "duplicates_removed '$DUPS' is not yes|no"
grep -q $'^path\tsha256\tbytes\trole$' manifest.tsv || fail "manifest.tsv column header must be: path sha256 bytes role"

# 2. every listed file exists, is non-empty, has the listed size; checksum unless --no-sha
grep -v '^#' manifest.tsv | tail -n +2 | while IFS=$'\t' read -r p sha bytes role; do
  [[ -s "$p" ]] || fail "listed in manifest but missing or empty: $p"
  [[ "$(stat -c %s "$p")" == "$bytes" ]] || fail "size mismatch: $p (manifest says $bytes bytes, found $(stat -c %s "$p"))"
  [[ "$NOSHA" == "--no-sha" ]] || [[ "$(sha256sum "$p" | cut -d' ' -f1)" == "$sha" ]] || fail "checksum mismatch: $p"
done

# 3. the fixed file set, and every required path must be LISTED in the manifest
for set in maf01 nomaf; do for ext in bed bim fam; do need "genotypes/cohort.hardcall.${set}.${ext}"; done; done
for f in ancestry/ancestry_calls.tsv ancestry/grafanc_raw.txt ancestry/strata.tsv ancestry/strata_rule.txt \
         structure/cohort.pcair.pcs.tsv structure/cohort.pcrelate.grm.rds structure/cohort.pcrelate.kinship.rds \
         phenotypes/phenotypes.tsv phenotypes/covariates.txt \
         qc/removed_duplicates.tsv qc/sample_accounting.tsv qc/qc_assertions.tsv run_record/params.tsv run_record/trace.tsv; do need "$f"; done
for n in "${CHROMS[@]}"; do need "genotypes/chr${n}.phased.dose.vcf.gz"; need "genotypes/chr${n}.phased.dose.vcf.gz.csi"; done
need local_ancestry/flare_panels.txt; need local_ancestry/ancestry_labels.txt; need local_ancestry/flare_index_map.tsv
for n in "${CHROMS[@]}"; do need "local_ancestry/flare/chr${n}.anc.vcf.gz"; need "local_ancestry/flare/chr${n}.anc.vcf.gz.csi"; done
while read -r p; do
  awk -F'\t' -v p="$p" '!/^#/ && $1==p {f=1} END{exit !f}' manifest.tsv || fail "required file is not listed in manifest.tsv: $p"
done < "$T/required"
awk -F'\t' '!/^#/ && $1 ~ /^genotypes\/cohort\.hardcall\.maf01\./ && $4 != "hardcall_maf01" {bad=1}
            !/^#/ && $1 ~ /^genotypes\/cohort\.hardcall\.nomaf\./ && $4 != "hardcall_nomaf" {bad=1} END{exit bad+0}' manifest.tsv \
  || fail "manifest role of a hard-call file is not hardcall_maf01 / hardcall_nomaf (see the contract's manifest row)"

# 4. both PLINK trios are internally consistent: bed bytes == 3 + ceil(n_samples/4) * n_variants; the two .fam
#    files are byte-identical (one sample set, one order); maf01 (the GWAS input, cohort MAF >= 1 %) is a subset
#    of nomaf's sites
for set in maf01 nomaf; do
  ns=$(wc -l < genotypes/cohort.hardcall.$set.fam); nv=$(wc -l < genotypes/cohort.hardcall.$set.bim)
  exp=$(( 3 + ((ns + 3) / 4) * nv )); got=$(stat -c %s genotypes/cohort.hardcall.$set.bed)
  [[ "$got" -eq "$exp" ]] || fail "cohort.hardcall.$set.bed is $got bytes; $ns samples x $nv variants implies $exp"
done
cmp -s genotypes/cohort.hardcall.nomaf.fam genotypes/cohort.hardcall.maf01.fam \
  || fail "nomaf and maf01 .fam files differ (the contract requires one byte-identical .fam)"
awk '{print $2}' genotypes/cohort.hardcall.maf01.bim | LC_ALL=C sort > "$T/bim_maf01"
awk '{print $2}' genotypes/cohort.hardcall.nomaf.bim | LC_ALL=C sort > "$T/bim_nomaf"
LC_ALL=C comm -23 "$T/bim_maf01" "$T/bim_nomaf" > "$T/extra"
NEXTRA=$(wc -l < "$T/extra")
[[ "$NEXTRA" -eq 0 ]] || fail "maf01 carries $NEXTRA site(s) absent from nomaf (first: $(head -1 "$T/extra")): maf01 must be a subset of the unfiltered set"
NS=$(wc -l < genotypes/cohort.hardcall.maf01.fam); NV=$(wc -l < genotypes/cohort.hardcall.maf01.bim); NV_ALL=$(wc -l < genotypes/cohort.hardcall.nomaf.bim)

# 5. one key: .fam IID (maf01 and nomaf share it) == phenotypes sample_id == IID == ancestry_calls sample_id == VCF samples
awk '{print $2}' genotypes/cohort.hardcall.maf01.fam > "$T/fam_order"
LC_ALL=C sort "$T/fam_order" > "$T/fam_iid"
H=$(head -1 phenotypes/phenotypes.tsv)
for c in FID IID sample_id mrd efs_time efs_event os_time os_event protocol; do
  grep -qP "(^|\t)${c}(\t|$)" <<<"$H" || fail "phenotypes.tsv lacks column '$c' (the GWAS runs MRD, EFS and OS; protocol is the env-MR-MEGA covariate)"; done
col() { awk -F'\t' -v n="$1" 'NR==1{for(i=1;i<=NF;i++) if($i==n) {print i; exit}}' phenotypes/phenotypes.tsv; }
IC=$(col IID); SC=$(col sample_id); PC=$(col protocol)
awk -F'\t' -v i="$IC" -v s="$SC" 'NR>1 && $i != $s {print $i; exit 1}' phenotypes/phenotypes.tsv > "$T/badrow" \
  || fail "phenotypes.tsv row where IID != sample_id (IID $(cat "$T/badrow")); the contract keys everything on one id"
awk -F'\t' -v p="$PC" 'NR>1 && ($p == "" || $p ~ /^[[:space:]]+$/) {exit 1}' phenotypes/phenotypes.tsv \
  || fail "phenotypes.tsv has a blank protocol (use NA where unknown; the value is a short protocol code, never free text)"
for pair in efs os; do
  tc=$(col "${pair}_time"); ec=$(col "${pair}_event")
  awk -F'\t' -v e="$ec" 'NR>1 && $e !~ /^(0|1|NA)$/ {exit 1}' phenotypes/phenotypes.tsv || fail "phenotypes.tsv ${pair}_event has a value other than 0, 1 or NA"
  awk -F'\t' -v t="$tc" 'NR>1 && $t !~ /^(NA|[0-9]+(\.[0-9]+)?)$/ {exit 1}' phenotypes/phenotypes.tsv || fail "phenotypes.tsv ${pair}_time has a non-numeric value"
done
for bad in mrn dob comments cyto_text oth_features death_cause relapse_site_oth; do
  grep -qiP "(^|\t)${bad}(\t|$)" <<<"$H" && fail "phenotypes.tsv carries forbidden column '$bad'"; done
awk -F'\t' -v c="$SC" 'NR>1{print $c}' phenotypes/phenotypes.tsv | LC_ALL=C sort > "$T/ph_id"
cmp -s "$T/fam_iid" "$T/ph_id" || fail "phenotypes.tsv sample_id set differs from .fam IID set (first diff: $(LC_ALL=C comm -3 "$T/fam_iid" "$T/ph_id" | head -1))"
ANC_HDR=$'sample_id\tAncGroupID\tMajorGroup\tGD1\tGD2\tGD3\tEA1\tEA2\tEA3\tEA4\tAF1\tAF2\tAF3\tEU1\tEU2\tEU3\tSA1\tSA2\tIC1\tIC2\tIC3\tPe\tPf\tPa'
[[ "$(head -1 ancestry/ancestry_calls.tsv)" == "$ANC_HDR" ]] || fail "ancestry_calls.tsv header is not the 24-column GrafAnc-normalised header (sample_id AncGroupID MajorGroup GD1-3 EA1-4 AF1-3 EU1-3 SA1-2 IC1-3 Pe Pf Pa)"
awk -F'\t' 'NR>1{print $1}' ancestry/ancestry_calls.tsv | LC_ALL=C sort > "$T/anc_id"
[[ "$(sort -u "$T/anc_id" | wc -l)" -eq "$(wc -l < "$T/anc_id")" ]] || fail "ancestry_calls.tsv has a duplicated sample_id (GrafAnc appends later rounds; the normaliser must not)"
cmp -s "$T/fam_iid" "$T/anc_id" || fail "ancestry_calls.tsv sample set differs from .fam IID set"
# the 46 AncGroupID codes GrafAnc 1.0.0 can emit, and the hundreds digit -> MajorGroup rule
GRAF_CODES=" 100 101 102 103 104 105 106 107 108 200 201 202 203 300 301 302 303 304 305 306 307 308 400 401 402 403 404 405 500 501 502 503 504 505 506 507 508 509 510 511 600 601 602 603 700 800 "
MSG=$(awk -F'\t' -v codes="$GRAF_CODES" 'BEGIN{split("AFR MEN EUR SAS EAS AMR OCN MIX",G," ")}
  NR>1 { if (index(codes, " " $2 " ")==0) {print "AncGroupID " $2 " is not a GrafAnc code (sample " $1 ")"; exit}
         if ($3 != G[int($2/100)]) {print "MajorGroup " $3 " does not match AncGroupID " $2 " (sample " $1 ")"; exit} }' ancestry/ancestry_calls.tsv)
[[ -z "$MSG" ]] || fail "ancestry_calls.tsv: $MSG"
head -1 ancestry/grafanc_raw.txt | grep -qP '^Sample\t' || fail "grafanc_raw.txt does not start with GrafAnc's header (first column Sample)"
if command -v bcftools >/dev/null; then
  for n in "${CHROMS[0]}" "${CHROMS[-1]}"; do   # first and last are enough to catch a wrong sample set or order; the full check is a run-time job
    bcftools query -l "genotypes/chr${n}.phased.dose.vcf.gz" | cmp -s - "$T/fam_order" || fail "chr${n}.phased.dose.vcf.gz samples differ from the .fam IID list or are in a different order"
    bcftools view -h "genotypes/chr${n}.phased.dose.vcf.gz" | grep -q '^##FORMAT=<ID=DS' || fail "chr${n}.phased.dose.vcf.gz has no FORMAT/DS"
    bcftools view -h "genotypes/chr${n}.phased.dose.vcf.gz" | grep -q '^##INFO=<ID=MRSQX' || fail "chr${n}.phased.dose.vcf.gz has no INFO/MRSQX (MagicalRsq-X score)"
  done
else
  [[ "$ALLOW_MISSING" -eq 1 ]] || fail "bcftools is not on PATH, so the VCF sample-order, FORMAT/DS and INFO/MRSQX checks cannot run (pass --allow-missing-tools to accept the bundle with those checks SKIPPED)"
  SKIPPED+=("VCF sample order / FORMAT/DS / INFO/MRSQX (no bcftools)")
fi

# 5b. strata and structure: every sample has a stratum; cohort PCs and one PC file per stratum; GRM dimnames
head -1 ancestry/strata.tsv | grep -qP '^sample_id\tstratum(\t|$)' || fail "strata.tsv header must start sample_id<TAB>stratum"
awk -F'\t' 'NR>1{print $1}' ancestry/strata.tsv | LC_ALL=C sort > "$T/st_id"
cmp -s "$T/fam_iid" "$T/st_id" || fail "strata.tsv sample set differs from .fam IID set"
head -1 structure/cohort.pcair.pcs.tsv | grep -qP '^sample_id\tPC1(\t|$)' || fail "cohort.pcair.pcs.tsv header must start sample_id<TAB>PC1"
awk -F'\t' 'NR>1{print $1}' structure/cohort.pcair.pcs.tsv | LC_ALL=C sort > "$T/pc_id"
cmp -s "$T/fam_iid" "$T/pc_id" || fail "cohort.pcair.pcs.tsv sample set differs from .fam IID set"
awk -F'\t' 'NR>1{print $2}' ancestry/strata.tsv | LC_ALL=C sort -u > "$T/strata"
while read -r S; do
  [[ -s "structure/strata/${S}.pcair.pcs.tsv" ]] || fail "no per-stratum PCs for stratum '$S' (structure/strata/${S}.pcair.pcs.tsv)"
  awk -F'\t' -v s="$S" 'NR>1 && $2==s {print $1}' ancestry/strata.tsv | LC_ALL=C sort > "$T/st_members"
  awk -F'\t' 'NR>1{print $1}' "structure/strata/${S}.pcair.pcs.tsv" | LC_ALL=C sort | cmp -s - "$T/st_members" \
    || fail "structure/strata/${S}.pcair.pcs.tsv samples differ from stratum '$S' members in strata.tsv"
done < "$T/strata"
if command -v Rscript >/dev/null; then
  FAMIDS="$T/fam_iid" Rscript --vanilla -e 'ids<-readLines(Sys.getenv("FAMIDS")); for(f in c("structure/cohort.pcrelate.grm.rds","structure/cohort.pcrelate.kinship.rds")){m<-readRDS(f); if(is.null(rownames(m))||!identical(sort(rownames(m)),sort(ids))) {cat(f,"\n"); quit(status=3)}}' \
    || fail "a structure .rds has no dimnames or its dimnames differ from the .fam IID set"
else
  [[ "$ALLOW_MISSING" -eq 1 ]] || fail "Rscript is not on PATH, so the GRM/kinship dimnames check cannot run (pass --allow-missing-tools to accept the bundle with it SKIPPED)"
  SKIPPED+=("GRM/kinship dimnames (no Rscript)")
fi

# 6. covariates named exist (whole-field match); labels alphabetical; index map explicit; QC files shaped as the contract says
while read -r cv || [[ -n "$cv" ]]; do
  [[ -z "$cv" ]] && continue
  awk -F'\t' -v c="$cv" 'NR==1{for(i=1;i<=NF;i++) if($i==c) f=1} END{exit !f}' phenotypes/phenotypes.tsv || fail "covariate '$cv' not in phenotypes.tsv"
done < phenotypes/covariates.txt
awk '{print $2}' local_ancestry/flare_panels.txt | LC_ALL=C sort -u > "$T/panel_labels"
LC_ALL=C sort -c local_ancestry/ancestry_labels.txt 2>/dev/null || fail "ancestry_labels.txt is not in alphabetical order"
cmp -s "$T/panel_labels" local_ancestry/ancestry_labels.txt || fail "ancestry_labels.txt is not the sorted set of labels in flare_panels.txt"
head -1 local_ancestry/flare_index_map.tsv | grep -qP '^anc_index\tlabel(\t|$)' || fail "flare_index_map.tsv header must be anc_index<TAB>label"
awk '{print $2}' local_ancestry/flare_panels.txt | awk '!seen[$0]++' > "$T/first_seen"
awk -F'\t' 'NR>1{print $2}' local_ancestry/flare_index_map.tsv | cmp -s - "$T/first_seen" \
  || fail "flare_index_map.tsv labels are not in the panel file's order of first appearance (that is how FLARE numbers anc0..ancK-1)"
awk -F'\t' 'NR>1 && $1 != NR-2 {exit 1}' local_ancestry/flare_index_map.tsv || fail "flare_index_map.tsv anc_index must be 0,1,2,... in row order"
head -1 qc/removed_duplicates.tsv | grep -qP '^kept_iid\tremoved_iid\tkinship\treason(\t|$)' || fail "removed_duplicates.tsv header must be kept_iid removed_iid kinship reason (header-only when nothing was removed)"
NR_DUP=$(($(wc -l < qc/removed_duplicates.tsv) - 1))
[[ "$DUPS" == "yes" || "$NR_DUP" -eq 0 ]] || fail "duplicates_removed=no but removed_duplicates.tsv lists $NR_DUP removals"
head -1 qc/sample_accounting.tsv | grep -qP '^step\tn_in\tn_out\tn_removed\treason(\t|$)' || fail "sample_accounting.tsv header must be step n_in n_out n_removed reason"
head -1 qc/qc_assertions.tsv | grep -qP '^id\tverdict\tn\tmessage(\t|$)' || fail "qc_assertions.tsv header must be id verdict n message"
FAILS=$(awk -F'\t' 'NR>1 && $2=="FAIL" {print $1}' qc/qc_assertions.tsv | paste -sd, -)
[[ -z "$FAILS" ]] || fail "qc_assertions.tsv carries FAIL verdict(s): $FAILS — an upstream check failed; the bundle is not analysis-ready"
awk -F'\t' 'NR>1 && $2 !~ /^(PASS|FAIL|INFO)$/ {exit 1}' qc/qc_assertions.tsv || fail "qc_assertions.tsv has a verdict other than PASS, FAIL or INFO"
[[ "$(awk -F'\t' 'NR>1 && $2=="PASS"' qc/qc_assertions.tsv | wc -l)" -ge 1 ]] || fail "qc_assertions.tsv has no PASS verdict — 'did not run' is not 'passed'"
NOTE=""; [[ ${#SKIPPED[@]} -eq 0 ]] || NOTE="; SKIPPED: $(IFS='; '; echo "${SKIPPED[*]}")"
echo "HANDOFF ACCEPTED: $B  ($NS samples, $NV GWAS-input variants at MAF >= 1 %, $NV_ALL unfiltered, chromosomes $CHROM_SPEC, server=$SERVER, duplicates_removed=$DUPS$NOTE)"
