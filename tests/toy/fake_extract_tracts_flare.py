#!/usr/bin/env python3
"""
A stand-in for Tractor's scripts/extract_tracts_flare.py, for tests/test_chassis.sh only
(the fixture copies it into the fake lai env as extract_tracts_flare.py).

Same interface and output as the real script as this repo uses it (--vcf, --num-ancs,
--output-dir, --compress-output; per ancestry index k, <vcf name>.anc<k>.hapcount.txt and
<vcf name>.anc<k>.dosage.txt, header CHROM POS ID REF ALT <samples>). The real script runs on the
cluster (the verification notes). FAKE_EXTRACT_DOSAGE_ONLY=1 makes it write dosage files only,
the failure the wrapper must catch.
"""
import argparse
import gzip
import os

ap = argparse.ArgumentParser()
ap.add_argument("--vcf", required=True)
ap.add_argument("--num-ancs", type=int, required=True)
ap.add_argument("--output-dir", default=".")
ap.add_argument("--compress-output", action="store_true")
a = ap.parse_args()
base = os.path.basename(a.vcf)
for suf in (".gz", ".vcf"):
    if base.endswith(suf):
        base = base[: -len(suf)]
kinds = ["dosage"] if os.environ.get("FAKE_EXTRACT_DOSAGE_ONLY") == "1" else ["hapcount", "dosage"]
opn = (lambda p: gzip.open(p + ".gz", "wt")) if a.compress_output else (lambda p: open(p, "w"))
outs = {(k, kind): opn(os.path.join(a.output_dir, f"{base}.anc{k}.{kind}.txt")) for k in range(a.num_ancs) for kind in kinds}
with gzip.open(a.vcf, "rt") as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#CHROM"):
            for o in outs.values():
                o.write("CHROM\tPOS\tID\tREF\tALT\t" + "\t".join(f[9:]) + "\n")
            continue
        fmt = f[8].split(":")
        igt, i1, i2 = fmt.index("GT"), fmt.index("AN1"), fmt.index("AN2")
        hap = [[0] * len(f[9:]) for _ in range(a.num_ancs)]
        dos = [[0] * len(f[9:]) for _ in range(a.num_ancs)]
        for j, cell in enumerate(f[9:]):
            v = cell.split(":")
            g = v[igt].split("|")
            for allele, anc in ((g[0], int(v[i1])), (g[1], int(v[i2]))):
                hap[anc][j] += 1
                dos[anc][j] += allele == "1"
        for k in range(a.num_ancs):
            for kind, vals in (("hapcount", hap[k]), ("dosage", dos[k])):
                if (k, kind) in outs:
                    outs[(k, kind)].write("\t".join(f[:5]) + "\t" + "\t".join(str(int(x)) for x in vals) + "\n")
for o in outs.values():
    o.close()
