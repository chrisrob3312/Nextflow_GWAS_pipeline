#!/usr/bin/env python3
"""
la_tractor_check.py: check a chromosome's Tractor files against what the adapter needs (CP-16)

    la_tractor_check.py --prefix data/tractor/cohort/chr22 --labels AFR,AMR,EUR
                        [--vcf chr22.anc.vcf.gz --index_map flare_index_map.tsv --rows 200]

For every label, <prefix>.hapcount.<label>.txt.gz and <prefix>.dosage.<label>.txt.gz must exist with
the header CHROM POS ID REF ALT <samples>, the same samples and the same variants in every file; on
every row and sample the hapcounts sum to 2 across the labels and each dosage is between 0 and its
hapcount. With --vcf (the FLARE output) and --index_map, the first --rows records are recomputed from
GT and AN1 / AN2 (FLARE's per-haplotype ancestry indices, mapped to labels) and must agree, so a
wrong index -> label mapping is caught, not just a malformed file. Standard library only.
"""

import argparse
import gzip
import sys

INFO = ["CHROM", "POS", "ID", "REF", "ALT"]


def die(msg):
    sys.exit(f"ERROR: {msg}")


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--vcf")
    ap.add_argument("--index_map")
    ap.add_argument("--rows", type=int, default=200)
    a = ap.parse_args()
    labels = a.labels.split(",")
    files = {(lab, kind): f"{a.prefix}.{kind}.{lab}.txt.gz" for lab in labels for kind in ("hapcount", "dosage")}
    handles = {}
    header = None
    for key, path in files.items():
        try:
            fh = opener(path)
        except FileNotFoundError:
            die(f"{path} does not exist (every ancestry needs a hapcount AND a dosage file)")
        h = fh.readline().rstrip("\n").split("\t")
        if h[:5] != INFO:
            die(f"{path}: header starts {' '.join(h[:5])}, not CHROM POS ID REF ALT")
        if header is None:
            header = h
        elif h != header:
            die(f"{path}: samples differ from {files[(labels[0], 'hapcount')]}")
        handles[key] = fh
    samples = header[5:]

    expect = {}
    flare_map = None
    if a.vcf:
        if not a.index_map:
            die("--vcf needs --index_map")
        with open(a.index_map) as fh:
            rows = [l.rstrip("\n").split("\t") for l in fh][1:]
        idx2lab = {r[0]: r[1] for r in rows}
        flare_map = None
        with opener(a.vcf) as fh:
            n = 0
            for line in fh:
                if line.startswith("##ANCESTRY=<"):
                    # FLARE names its indices itself: <EUR=0,AFR=1,...>. The index map must say the same.
                    flare_map = {kv.split("=")[1]: kv.split("=")[0] for kv in line.strip()[len("##ANCESTRY=<"):-1].split(",")}
                    if flare_map != idx2lab:
                        die(f"{a.index_map} maps {sorted(idx2lab.items())}, but the FLARE VCF's own ##ANCESTRY line says "
                            f"{sorted(flare_map.items())}: the files would be labelled with the wrong ancestries")
                    continue
                if line.startswith("##"):
                    continue
                f = line.rstrip("\n").split("\t")
                if line.startswith("#CHROM"):
                    if f[9:] != samples:
                        die(f"{a.vcf}: samples differ from the Tractor files'")
                    continue
                fmt = f[8].split(":")
                try:
                    igt, ia1, ia2 = fmt.index("GT"), fmt.index("AN1"), fmt.index("AN2")
                except ValueError:
                    die(f"{a.vcf}: FORMAT {f[8]} lacks GT, AN1 or AN2 (FLARE's per-haplotype ancestry)")
                hap = {lab: [0] * len(samples) for lab in labels}
                dos = {lab: [0] * len(samples) for lab in labels}
                for j, cell in enumerate(f[9:]):
                    v = cell.split(":")
                    gt = v[igt].replace("/", "|").split("|")
                    for allele, anc in ((gt[0], v[ia1]), (gt[1], v[ia2])):
                        lab = idx2lab.get(anc)
                        if lab is None:
                            die(f"{a.vcf}: ancestry index {anc} is not in {a.index_map}")
                        hap[lab][j] += 1
                        dos[lab][j] += 1 if allele == "1" else 0
                expect[n] = (f[0], f[1], hap, dos)
                n += 1
                if n >= a.rows:
                    break

    nrow = 0
    while True:
        lines = {k: h.readline() for k, h in handles.items()}
        empty = [k for k, l in lines.items() if not l]
        if empty:
            if len(empty) != len(lines):
                die(f"{files[empty[0]]} ends before the other files (row {nrow + 1})")
            break
        parts = {k: l.rstrip("\n").split("\t") for k, l in lines.items()}
        ref = parts[(labels[0], "hapcount")]
        for k, p in parts.items():
            if p[:5] != ref[:5]:
                die(f"{files[k]} row {nrow + 1}: {':'.join(p[:5])} vs {':'.join(ref[:5])} (files out of step)")
        hsum = [0] * len(samples)
        for lab in labels:
            hc = parts[(lab, "hapcount")][5:]
            ds = parts[(lab, "dosage")][5:]
            for j in range(len(samples)):
                h, d = int(hc[j]), int(ds[j])
                if not 0 <= d <= h:
                    die(f"{ref[0]}:{ref[1]} sample {samples[j]}: {lab} dosage {d} outside 0..hapcount {h}")
                hsum[j] += h
            if nrow in expect:
                c, p_, hap, dos = expect[nrow]
                if (c, p_) != (ref[0], ref[1]):
                    die(f"row {nrow + 1}: Tractor files at {ref[0]}:{ref[1]}, the VCF at {c}:{p_}")
                if [int(x) for x in hc] != hap[lab] or [int(x) for x in ds] != dos[lab]:
                    die(f"{ref[0]}:{ref[1]}: the {lab} hapcount / dosage differ from GT + AN1 / AN2 in the VCF "
                        "(an index -> label mapping that does not follow the index map?)")
        bad = [samples[j] for j in range(len(samples)) if hsum[j] != 2]
        if bad:
            die(f"{ref[0]}:{ref[1]}: hapcounts sum to {hsum[samples.index(bad[0])]} for {bad[0]} (must be 2)")
        nrow += 1
    for h in handles.values():
        h.close()
    note = ""
    if a.vcf:
        note = f"; first {len(expect)} records agree with {a.vcf}"
        note += "; the index map matches its ##ANCESTRY line" if flare_map else "; NOTE: no ##ANCESTRY line in the VCF, so the index map was not checked against FLARE's own"
    print(f"ok: {a.prefix}: {nrow} variants x {len(samples)} samples, labels {','.join(labels)}; hapcounts sum to 2, "
          f"0 <= dosage <= hapcount everywhere" + note)


if __name__ == "__main__":
    main()
