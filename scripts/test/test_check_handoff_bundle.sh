#!/usr/bin/env bash
# Exercise bin/check_handoff_bundle.sh: an intact stub is accepted; each deliberate break is refused
# with its own message. usage: scripts/test/test_check_handoff_bundle.sh [work_dir]
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"
CHECK="$REPO/bin/check_handoff_bundle.sh"; MK="$HERE/make_stub_bundle.sh"
W="${1:-$(mktemp -d)}"; mkdir -p "$W"; pass=0; failn=0; skipn=0
ok()   { echo "  ok   $1"; pass=$((pass + 1)); }
bad()  { echo "  FAIL $1"; failn=$((failn + 1)); }
TOOLS=(); command -v bcftools >/dev/null && command -v Rscript >/dev/null || TOOLS=(--allow-missing-tools)
expect_accept() { # name, bundle, args...
  local name=$1 b=$2; shift 2
  if out=$("$CHECK" "$b" "${TOOLS[@]}" "$@" 2>&1); then ok "$name: accepted"; else bad "$name: refused — $out"; fi; }
expect_refuse() { # name, bundle, pattern, args...
  local name=$1 b=$2 pat=$3; shift 3
  out=$("$CHECK" "$b" "${TOOLS[@]}" "$@" 2>&1); rc=$?
  if [[ $rc -ne 0 && "$out" == *"$pat"* ]]; then ok "$name: refused with '$pat'"; else bad "$name: rc=$rc, output: $out"; fi; }
remanifest() { # rewrite checksums after a deliberate edit so the break under test is the first failure
  ( cd "$1" && { grep '^#' manifest.tsv; printf 'path\tsha256\tbytes\trole\n'
      dirs=(); for d in genotypes ancestry structure local_ancestry phenotypes qc run_record; do [[ -d "$d" ]] && dirs+=("$d"); done
      find "${dirs[@]}" -type f | LC_ALL=C sort | while read -r p; do
        r=other; case "$p" in genotypes/cohort.hardcall.maf01.*) r=hardcall_maf01;; genotypes/cohort.hardcall.nomaf.*) r=hardcall_nomaf;; esac
        printf '%s\t%s\t%s\t%s\n' "$p" "$(sha256sum "$p" | cut -d' ' -f1)" "$(stat -c %s "$p")" "$r"; done; } > manifest.new && mv manifest.new manifest.tsv ); }
fresh() { "$MK" "$W/$1" "${@:2}" >/dev/null; echo "$W/$1"; }

b=$(fresh intact);                                         expect_accept "intact bundle" "$b"
b=$(fresh intact_nosha);                                   expect_accept "intact bundle, --no-sha" "$b" --no-sha
b=$(fresh ladder --chromosomes 21,22);                     expect_accept "ladder bundle on chr21-22 with --chromosomes 21,22" "$b" --chromosomes 21,22
b=$(fresh ladder_default --chromosomes 21,22);             expect_refuse "ladder bundle checked as 1-22" "$b" "chr1.phased.dose.vcf.gz is missing or empty"
b=$(fresh extra_site); printf '3\tchr3:5:A:T\t0\t5\tT\tA\n' >> "$b/genotypes/cohort.hardcall.maf01.bim"; printf '\0\0' >> "$b/genotypes/cohort.hardcall.maf01.bed"; remanifest "$b"
                                                           expect_refuse "maf01 site absent from nomaf" "$b" "site(s) absent from nomaf"
b=$(fresh bed_short); truncate -s 10 "$b/genotypes/cohort.hardcall.maf01.bed"; remanifest "$b"
                                                           expect_refuse "truncated maf01 .bed" "$b" "cohort.hardcall.maf01.bed is 10 bytes"
b=$(fresh fam_differ); sed -i 's/^F_S008\tS008/F_S008\tS999/' "$b/genotypes/cohort.hardcall.nomaf.fam"; remanifest "$b"
                                                           expect_refuse "nomaf and maf01 .fam differ" "$b" ".fam files differ"
b=$(fresh old_pheno); sed -i '1s/efs_time\tefs_event\tos_time\tos_event/time\tevent\tos_time\tos_event/' "$b/phenotypes/phenotypes.tsv"; remanifest "$b"
                                                           expect_refuse "old two-column time/event header" "$b" "lacks column 'efs_time'"
b=$(fresh event2); awk -F'\t' 'BEGIN{OFS="\t"} NR==2{$8=2} {print}' "$b/phenotypes/phenotypes.tsv" > "$b/p.tmp" && mv "$b/p.tmp" "$b/phenotypes/phenotypes.tsv"; remanifest "$b"
                                                           expect_refuse "os_event = 2" "$b" "os_event has a value other than 0, 1 or NA"
b=$(fresh mrn); sed -i '1s/$/\tmrn/' "$b/phenotypes/phenotypes.tsv"; sed -i '2,$s/$/\t1/' "$b/phenotypes/phenotypes.tsv"; remanifest "$b"
                                                           expect_refuse "forbidden column mrn" "$b" "forbidden column 'mrn'"
b=$(fresh no_lai); rm -r "$b/local_ancestry"; remanifest "$b"
                                                           expect_refuse "bundle without local_ancestry/" "$b" "local_ancestry/flare_panels.txt is missing or empty"
b=$(fresh labels_unsorted); printf 'EUR\nAFR\nAMR\n' > "$b/local_ancestry/ancestry_labels.txt"; remanifest "$b"
                                                           expect_refuse "ancestry_labels.txt not alphabetical" "$b" "not in alphabetical order"
b=$(fresh index_order); printf 'anc_index\tlabel\n0\tAFR\n1\tEUR\n2\tAMR\n' > "$b/local_ancestry/flare_index_map.tsv"; remanifest "$b"
                                                           expect_refuse "flare_index_map.tsv not panel order" "$b" "order of first appearance"
b=$(fresh bad_code); sed -i '2s/\t300\tEUR/\t999\tEUR/' "$b/ancestry/ancestry_calls.tsv"; remanifest "$b"
                                                           expect_refuse "AncGroupID 999" "$b" "is not a GrafAnc code"
b=$(fresh wrong_group); sed -i '2s/\t300\tEUR/\t300\tAFR/' "$b/ancestry/ancestry_calls.tsv"; remanifest "$b"
                                                           expect_refuse "MajorGroup disagrees with code" "$b" "does not match AncGroupID"
b=$(fresh no_pass); printf 'id\tverdict\tn\tmessage\nM6_note\tINFO\t8\tinformational only\n' > "$b/qc/qc_assertions.tsv"; remanifest "$b"
                                                           expect_refuse "no PASS in qc_assertions (INFO only)" "$b" "no PASS verdict"
b=$(fresh dups_contradict); printf 'S002\tS001\tduplicate\n' >> "$b/qc/removed_duplicates.tsv"; remanifest "$b"
                                                           expect_refuse "duplicates_removed=no with removals listed" "$b" "duplicates_removed=no but"
b=$(fresh covariate_missing); printf 'not_a_column\n' >> "$b/phenotypes/covariates.txt"; remanifest "$b"
                                                           expect_refuse "covariate not a phenotype column" "$b" "covariate 'not_a_column' not in"
b=$(fresh checksum_drift); sed -i 's/flare/flarf/' "$b/run_record/params.tsv"
                                                           expect_refuse "checksum drift (same size)" "$b" "checksum mismatch"
b=$(fresh size_drift); printf 'x' >> "$b/run_record/params.tsv"
                                                           expect_refuse "size drift" "$b" "size mismatch"

if ! command -v bcftools >/dev/null || ! command -v Rscript >/dev/null; then
  b=$(fresh tools_default); out=$("$CHECK" "$b" 2>&1); rc=$?
  if [[ $rc -ne 0 && "$out" == *"is not on PATH"* ]]; then ok "missing tools refused by default (no --allow-missing-tools)"; else bad "missing tools: rc=$rc $out"; fi
fi
# M-06: one sample key; PLINK's FID_IID VCF naming is named as the cause
if command -v bcftools >/dev/null; then
  rehead() { local b=$1 names=$2; for v in "$b"/genotypes/chr*.phased.dose.vcf.gz; do
      bcftools reheader -s "$names" -o "$v.new" "$v" && mv "$v.new" "$v" && bcftools index -f -c "$v"; done; remanifest "$b"; }
  b=$(fresh fid_iid); awk '{print $1 "_" $2}' "$b/genotypes/cohort.hardcall.maf01.fam" > "$W/fid_iid.txt"; rehead "$b" "$W/fid_iid.txt"
  expect_refuse "VCF header uses FID_IID" "$b" "sample ID mismatch: chr1.phased.dose.vcf.gz vs fam: the VCF names its samples FID_IID (e.g. F_S001_S001, .fam IID S001)"
  out=$("$CHECK" "$b" "${TOOLS[@]}" --qc-dir "$W/qc_fid" 2>&1)
  grep -q $'^handoff_sample_ids\tFAIL\t8\tsample ID mismatch: chr1.phased.dose.vcf.gz vs fam' "$W"/qc_fid/assertions.d/*.tsv 2>/dev/null \
    && grep -q $'^handoff_sample_ids\tcohort\tNA\tNA\t8\t0\t8\t' "$W"/qc_fid/accounting.d/*.tsv 2>/dev/null \
    && ok "--qc-dir: the FID_IID refusal is filed as handoff_sample_ids FAIL with an accounting row 8 -> 0" || bad "qc-dir fail rows: $(cat "$W"/qc_fid/*/*.tsv 2>&1)"
  b=$(fresh vcf_order); awk '{print $2}' "$b/genotypes/cohort.hardcall.maf01.fam" | sed '1{h;d};2{G}' > "$W/swapped.txt"; rehead "$b" "$W/swapped.txt"
  expect_refuse "VCF samples in another order" "$b" "sample ID mismatch: chr1.phased.dose.vcf.gz vs fam: same samples, different order (first difference at line 1)"
  b=$(fresh ids_agree); out=$("$CHECK" "$b" "${TOOLS[@]}" --qc-dir "$W/qc_ok" 2>&1); rc=$?
  [[ $rc -eq 0 ]] && grep -q $'^handoff_sample_ids\tPASS\t8\tVCF samples (chr1, chr22) == .fam IID in order == phenotypes.tsv IID == sample_id$' "$W"/qc_ok/assertions.d/*.tsv \
    && grep -q $'^handoff_sample_ids\tcohort\tNA\tNA\t8\t8\t0\t' "$W"/qc_ok/accounting.d/*.tsv \
    && ok "IDs agree: accepted; --qc-dir gets handoff_sample_ids PASS and an accounting row 8 -> 8" || bad "ids agree: rc=$rc $out"
fi
b=$(fresh pheno_id); sed -i 's/^F_S008\tS008\tS008\t/F_S008\tS999\tS999\t/' "$b/phenotypes/phenotypes.tsv"; remanifest "$b"
                                                           expect_refuse "phenotype sample_id not a .fam IID" "$b" "sample ID mismatch: phenotypes.tsv vs fam: the sample_id set differs from the .fam IID set"
b=$(fresh rel_path); ( cd "$W" && "$CHECK" rel_path "${TOOLS[@]}" >/dev/null 2>&1 ) && ok "relative bundle path accepted" || bad "relative bundle path refused"
b=$(fresh many_extra); awk 'BEGIN{for(i=1;i<=200000;i++) printf "3\tchr3:%d:A:T\t0\t%d\tT\tA\n", i, i}' >> "$b/genotypes/cohort.hardcall.maf01.bim"; head -c $((2 * 200000)) /dev/zero >> "$b/genotypes/cohort.hardcall.maf01.bed"; remanifest "$b"
                                                           expect_refuse "200,000 maf01 sites absent from nomaf (no SIGPIPE death)" "$b" "maf01 carries 200000 site(s) absent from nomaf"
b=$(fresh no_protocol); awk -F'\t' 'BEGIN{OFS="\t"} {$9=""; print}' "$b/phenotypes/phenotypes.tsv" | sed '1s/\t\t/\t/' > "$b/p.tmp"; sed -i '1s/os_event\tsex/os_event\tsex/' "$b/p.tmp"; awk -F'\t' 'BEGIN{OFS="\t"} NR==1{$9="dropme"} {print}' "$b/p.tmp" > "$b/phenotypes/phenotypes.tsv"; remanifest "$b"
                                                           expect_refuse "phenotypes without a protocol column" "$b" "lacks column 'protocol'"
b=$(fresh blank_protocol); awk -F'\t' 'BEGIN{OFS="\t"} NR==3{$9=""} {print}' "$b/phenotypes/phenotypes.tsv" > "$b/p.tmp" && mv "$b/p.tmp" "$b/phenotypes/phenotypes.tsv"; remanifest "$b"
                                                           expect_refuse "blank protocol value" "$b" "blank protocol"
b=$(fresh iid_ne_sid); awk -F'\t' 'BEGIN{OFS="\t"} NR==2{$3="S002"} NR==3{$3="S001"} {print}' "$b/phenotypes/phenotypes.tsv" > "$b/p.tmp" && mv "$b/p.tmp" "$b/phenotypes/phenotypes.tsv"; remanifest "$b"
                                                           expect_refuse "IID differs from sample_id on a row" "$b" "IID != sample_id"
b=$(fresh fail_verdict); printf 'M6_het\tFAIL\t8\tthree samples beyond 3 SD\n' >> "$b/qc/qc_assertions.tsv"; remanifest "$b"
                                                           expect_refuse "a FAIL verdict in qc_assertions" "$b" "FAIL verdict(s): M6_het"
b=$(fresh unlisted); grep -v $'^ancestry/strata_rule.txt\t' "$b/manifest.tsv" > "$b/m.tmp" && mv "$b/m.tmp" "$b/manifest.tsv"
                                                           expect_refuse "required file missing from the manifest" "$b" "not listed in manifest.tsv: ancestry/strata_rule.txt"
b=$(fresh dup_header); printf 'sample_id\tkept_of\treason\n' > "$b/qc/removed_duplicates.tsv"; remanifest "$b"
                                                           expect_refuse "removed_duplicates.tsv with the wrong header" "$b" "header must be kept_iid removed_iid kinship reason"
b=$(fresh strata_prefix); sed -i '1s/stratum/stratumX/' "$b/ancestry/strata.tsv"; remanifest "$b"
                                                           expect_refuse "strata header column named stratumX" "$b" "strata.tsv header must start"
b=$(fresh cov_regex); printf 'age.at.diagnosis\n' >> "$b/phenotypes/covariates.txt"; remanifest "$b"
                                                           expect_refuse "covariate name with regex dots is not a column" "$b" "covariate 'age.at.diagnosis' not in"
b=$(fresh cov_nonl); printf 'not_a_column' >> "$b/phenotypes/covariates.txt"; remanifest "$b"
                                                           expect_refuse "last covariate without a trailing newline" "$b" "covariate 'not_a_column' not in"
b=$(fresh bad_role); sed -i $'s/\thardcall_maf01$/\thardcall/' "$b/manifest.tsv"
                                                           expect_refuse "wrong manifest role on a hard-call row" "$b" "manifest role of a hard-call file"
b=$(fresh clean_on_refuse); printf 'x' >> "$b/run_record/params.tsv"; "$CHECK" "$b" "${TOOLS[@]}" >/dev/null 2>&1 || true
                                                           if ls -A "$b" | grep -q '^\.'; then bad "refused bundle left scratch dotfiles"; else ok "refused bundle left no scratch files"; fi

# Per-haplotype dosages (FORMAT/HDS). These edit the stub's real dosage VCFs, so they need bcftools; without it
# the stub's VCFs are placeholders, the checker SKIPS its VCF checks, and these cases are reported as skipped.
vcf_filter() { # bundle, chr, filter command...: VCF text through the filter, re-bgzipped and re-indexed in place
  local v="$1/genotypes/chr$2.phased.dose.vcf.gz"; shift 2
  bcftools view "$v" | "$@" | bcftools view -Oz -o "$v.tmp" - && mv "$v.tmp" "$v" && bcftools index -f -c "$v"; }
s003_first() { # replace DS:HDS of sample S003 (column 12) on the first record; GT kept
  awk -v cell="$1" 'BEGIN { FS = OFS = "\t" } !/^#/ && !done { split($12, a, ":"); $12 = a[1] ":" cell; done = 1 } { print }'; }
HDS_CASES=("dosage VCF without FORMAT/HDS (chr1)" "dosage VCF without FORMAT/HDS on a middle chromosome (chr12)" "FORMAT/HDS declared Number=1"
           "FORMAT/HDS inconsistent with FORMAT/DS (chr22)" "FORMAT/HDS missing where FORMAT/DS is present" "FORMAT/HDS off FORMAT/DS by exactly 0.01 (rounding)")
if command -v bcftools >/dev/null; then
  b=$(fresh no_hds); vcf_filter "$b" 1 bcftools annotate -x FORMAT/HDS -; remanifest "$b"
                                                           expect_refuse "${HDS_CASES[0]}" "$b" "chr1.phased.dose.vcf.gz has no FORMAT/HDS header line"
  b=$(fresh no_hds_mid); vcf_filter "$b" 12 bcftools annotate -x FORMAT/HDS -; remanifest "$b"
                                                           expect_refuse "${HDS_CASES[1]}" "$b" "chr12.phased.dose.vcf.gz has no FORMAT/HDS header line"
  b=$(fresh hds_number1); v="$b/genotypes/chr1.phased.dose.vcf.gz"
  bcftools view -h "$v" | sed 's/^##FORMAT=<ID=HDS,Number=2,/##FORMAT=<ID=HDS,Number=1,/' > "$W/hds_number1.hdr"
  bcftools reheader -h "$W/hds_number1.hdr" -o "$v.tmp" "$v" && mv "$v.tmp" "$v" && bcftools index -f -c "$v"; remanifest "$b"
                                                           expect_refuse "${HDS_CASES[2]}" "$b" "chr1.phased.dose.vcf.gz declares FORMAT/HDS with Number=1, not Number=2"
  b=$(fresh hds_vs_ds); vcf_filter "$b" 22 s003_first '0.500:0.250,0.750'; remanifest "$b"
                                                           expect_refuse "${HDS_CASES[3]}" "$b" "chr22.phased.dose.vcf.gz: FORMAT/HDS disagrees with FORMAT/DS (|HDS1+HDS2-DS| must be <= 0.01) at chr22:10000 sample S003 (HDS 0.25,0.75 sums to 1, DS 0.5)"
  b=$(fresh hds_missing); vcf_filter "$b" 1 s003_first '0.500:.'; remanifest "$b"
                                                           expect_refuse "${HDS_CASES[4]}" "$b" "at chr1:10000 sample S003 (DS 0.5, HDS .: not two haplotype values against one DS)"
  b=$(fresh hds_rounding); vcf_filter "$b" 1 s003_first '0.500:0.250,0.240'; remanifest "$b"
                                                           expect_accept "${HDS_CASES[5]}" "$b"
else
  for c in "${HDS_CASES[@]}"; do echo "  skip $c: needs bcftools"; skipn=$((skipn + 1)); done
fi
echo "passed $pass, failed $failn$( [[ $skipn -gt 0 ]] && echo ", skipped $skipn (need bcftools)")"; [[ $failn -eq 0 ]]
