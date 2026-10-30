#!/usr/bin/env bash
# Exercise bin/check_handoff_bundle.sh: an intact stub is accepted; each deliberate break is refused
# with its own message. usage: scripts/test/test_check_handoff_bundle.sh [work_dir]
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"
CHECK="$REPO/bin/check_handoff_bundle.sh"; MK="$HERE/make_stub_bundle.sh"
W="${1:-$(mktemp -d)}"; mkdir -p "$W"; pass=0; failn=0
ok()   { echo "  ok   $1"; pass=$((pass + 1)); }
bad()  { echo "  FAIL $1"; failn=$((failn + 1)); }
expect_accept() { # name, bundle, args...
  local name=$1 b=$2; shift 2
  if out=$("$CHECK" "$b" "$@" 2>&1); then ok "$name: accepted"; else bad "$name: refused — $out"; fi; }
expect_refuse() { # name, bundle, pattern, args...
  local name=$1 b=$2 pat=$3; shift 3
  out=$("$CHECK" "$b" "$@" 2>&1); rc=$?
  if [[ $rc -ne 0 && "$out" == *"$pat"* ]]; then ok "$name: refused with '$pat'"; else bad "$name: rc=$rc, output: $out"; fi; }
remanifest() { # rewrite checksums after a deliberate edit so the break under test is the first failure
  ( cd "$1" && { grep '^#' manifest.tsv; printf 'path\tsha256\tbytes\trole\n'
      dirs=(); for d in genotypes ancestry structure local_ancestry phenotypes qc run_record; do [[ -d "$d" ]] && dirs+=("$d"); done
      find "${dirs[@]}" -type f | LC_ALL=C sort | while read -r p; do
        printf '%s\t%s\t%s\t%s\n' "$p" "$(sha256sum "$p" | cut -d' ' -f1)" "$(stat -c %s "$p")" "x"; done; } > manifest.new && mv manifest.new manifest.tsv ); }
fresh() { "$MK" "$W/$1" "${@:2}" >/dev/null; echo "$W/$1"; }

b=$(fresh intact);                                         expect_accept "intact bundle" "$b"
b=$(fresh intact_nosha);                                   expect_accept "intact bundle, --no-sha" "$b" --no-sha
b=$(fresh ladder --chromosomes 21,22);                     expect_accept "ladder bundle on chr21-22 with --chromosomes 21,22" "$b" --chromosomes 21,22
b=$(fresh ladder_default --chromosomes 21,22);             expect_refuse "ladder bundle checked as 1-22" "$b" "chr1.phased.dose.vcf.gz is missing or empty"
b=$(fresh extra_site); printf '3\tchr3:5:A:T\t0\t5\tT\tA\n' >> "$b/genotypes/cohort.hardcall.maf01.bim"; printf '\0\0' >> "$b/genotypes/cohort.hardcall.maf01.bed"; remanifest "$b"
                                                           expect_refuse "maf01 site absent from nomaf" "$b" "maf01 carries a site absent from nomaf"
b=$(fresh bed_short); truncate -s 10 "$b/genotypes/cohort.hardcall.maf01.bed"; remanifest "$b"
                                                           expect_refuse "truncated maf01 .bed" "$b" "cohort.hardcall.maf01.bed is 10 bytes"
b=$(fresh fam_differ); sed -i 's/^F_S008\tS008/F_S008\tS999/' "$b/genotypes/cohort.hardcall.nomaf.fam"; remanifest "$b"
                                                           expect_refuse "nomaf and maf01 .fam differ" "$b" "list different samples"
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
b=$(fresh no_pass); printf 'id\tverdict\tn\tmessage\nM6_callrate_95\tFAIL\t8\tx\n' > "$b/qc/qc_assertions.tsv"; remanifest "$b"
                                                           expect_refuse "no PASS in qc_assertions" "$b" "no PASS rows"
b=$(fresh dups_contradict); printf 'S002\tS001\tduplicate\n' >> "$b/qc/removed_duplicates.tsv"; remanifest "$b"
                                                           expect_refuse "duplicates_removed=no with removals listed" "$b" "duplicates_removed=no but"
b=$(fresh covariate_missing); printf 'not_a_column\n' >> "$b/phenotypes/covariates.txt"; remanifest "$b"
                                                           expect_refuse "covariate not a phenotype column" "$b" "covariate 'not_a_column' not in"
b=$(fresh checksum_drift); printf 'x' >> "$b/run_record/params.tsv"
                                                           expect_refuse "checksum drift" "$b" "checksum mismatch"
echo "passed $pass, failed $failn"; [[ $failn -eq 0 ]]
