#!/usr/bin/env python3
"""
la_tractor_layout.py: the cohort's Tractor files -> one stratum's (CP-16). Standard library only.

shares: how much of each ancestry a stratum's haplotypes carry (deciding its ancestries.txt)
    la_tractor_layout.py shares --cohort_prefixes data/tractor/cohort/chr1,...,chr22 --labels AFR,AMR,EUR
                         --samples results/strata/S.samples.txt [--every 50]
  prints "label<TAB>share" for each label: its haplotypes / all haplotypes of the stratum's samples,
  over every --every-th variant of the given chromosomes.

split: one chromosome of one stratum
    la_tractor_layout.py split --cohort_prefix data/tractor/cohort/chr22 --labels AFR,AMR,EUR
                         --keep AFR,EUR --samples results/strata/S.samples.txt --out_prefix data/tractor/S/chr22
                         [--fold_into EUR]
  Columns: CHROM POS ID REF ALT, then the stratum's samples (in samples.txt order).
  Two or more kept ancestries: <out>.hapcount.<L>.txt.gz and <out>.dosage.<L>.txt.gz for each kept
  label. A label not kept is folded into --fold_into (its haplotypes and alleles added there), so
  every sample's hapcounts still sum to 2 and no allele is lost.
  One kept ancestry (k = 1): one file, <out>.ancdose.<L>.tsv.gz, the hard-call ALT count summed over
  every ancestry (the genotype), and no hapcount file (the adapter takes 2 haplotypes per sample;
  the same hard calls as the admixed strata, docs/SUMSTATS_SCHEMA.md).
"""
import argparse
import gzip
import sys

INFO = 5


def die(msg):
    sys.exit(f"ERROR: {msg}")


def read_samples(path):
    with open(path) as fh:
        s = [l.strip() for l in fh if l.strip()]
    if not s:
        die(f"{path} lists no samples")
    return s


def open_set(prefix, labels, kinds=("hapcount", "dosage")):
    hs = {}
    for lab in labels:
        for kind in kinds:
            path = f"{prefix}.{kind}.{lab}.txt.gz"
            try:
                fh = gzip.open(path, "rt")
            except FileNotFoundError:
                die(f"{path} does not exist (run slurm/submit_la_extract.sh first)")
            hs[(lab, kind)] = (fh, fh.readline().rstrip("\n").split("\t"), path)
    return hs


def shares(a):
    labels = a.labels.split(",")
    keep = set(read_samples(a.samples))
    tot = {lab: 0 for lab in labels}
    for prefix in a.cohort_prefixes.split(","):
        hs = open_set(prefix, labels, ("hapcount",))
        for (lab, _), (fh, hdr, path) in hs.items():
            idx = [i for i, s in enumerate(hdr) if i >= INFO and s in keep]
            if not idx:
                die(f"{path}: none of the stratum's samples")
            for n, line in enumerate(fh):
                if n % a.every:
                    continue
                f = line.rstrip("\n").split("\t")
                tot[lab] += sum(int(f[i]) for i in idx)
            fh.close()
    all_h = sum(tot.values())
    if all_h == 0:
        die("no haplotypes counted")
    for lab in labels:
        print(f"{lab}\t{tot[lab] / all_h:.6f}")


def split(a):
    labels = a.labels.split(",")
    keep_labels = a.keep.split(",")
    for k in keep_labels:
        if k not in labels:
            die(f"--keep {k} is not one of --labels {a.labels}")
    samples = read_samples(a.samples)
    hs = open_set(a.cohort_prefix, labels)
    hdr0 = hs[(labels[0], "hapcount")][1]
    pos = {s: i for i, s in enumerate(hdr0)}
    absent = [s for s in samples if s not in pos]
    if absent:
        die(f"{hs[(labels[0], 'hapcount')][2]} lacks {len(absent)} of the stratum's samples, e.g. {absent[0]}")
    for key, (fh, hdr, path) in hs.items():
        if hdr != hdr0:
            die(f"{path}: other samples than {hs[(labels[0], 'hapcount')][2]}")
    idx = [pos[s] for s in samples]
    single = len(keep_labels) == 1
    fold = a.fold_into or (keep_labels[0] if single else None)
    if not single and fold not in keep_labels:
        die(f"--fold_into {fold} must be one of the kept labels {a.keep}")
    head = "\t".join(hdr0[:INFO] + samples) + "\n"
    if single:
        outs = {"ancdose": gzip.open(f"{a.out_prefix}.ancdose.{keep_labels[0]}.tsv.gz", "wt")}
    else:
        outs = {(lab, kind): gzip.open(f"{a.out_prefix}.{kind}.{lab}.txt.gz", "wt") for lab in keep_labels for kind in ("hapcount", "dosage")}
    for o in outs.values():
        o.write(head)
    n = 0
    while True:
        lines = {k: v[0].readline() for k, v in hs.items()}
        if not any(lines.values()):
            break
        if not all(lines.values()):
            die(f"the cohort files end at different rows (row {n + 1})")
        parts = {k: l.rstrip("\n").split("\t") for k, l in lines.items()}
        info = parts[(labels[0], "hapcount")][:INFO]
        for k, p in parts.items():
            if p[:INFO] != info:
                die(f"{hs[k][2]} row {n + 1}: {':'.join(p[:INFO])} vs {':'.join(info)} (files out of step)")
        vals = {k: [int(parts[k][i]) for i in idx] for k in parts}
        if single:
            dose = [sum(vals[(lab, "dosage")][j] for lab in labels) for j in range(len(idx))]
            outs["ancdose"].write("\t".join(info + [str(x) for x in dose]) + "\n")
        else:
            for kind in ("hapcount", "dosage"):
                acc = {lab: list(vals[(lab, kind)]) for lab in keep_labels}
                for lab in labels:
                    if lab not in keep_labels:
                        acc[fold] = [x + y for x, y in zip(acc[fold], vals[(lab, kind)])]
                for lab in keep_labels:
                    outs[(lab, kind)].write("\t".join(info + [str(x) for x in acc[lab]]) + "\n")
        n += 1
    for o in outs.values():
        o.close()
    kinds = "one ancdose file (k = 1: summed hard-call dose, no hapcount)" if single else \
        f"hapcount + dosage for {','.join(keep_labels)}" + (f" (other labels folded into {fold})" if set(labels) - set(keep_labels) else "")
    print(f"{a.out_prefix}: {n} variants x {len(samples)} samples; {kinds}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="mode", required=True)
    s1 = sub.add_parser("shares")
    s1.add_argument("--cohort_prefixes", required=True)
    s1.add_argument("--labels", required=True)
    s1.add_argument("--samples", required=True)
    s1.add_argument("--every", type=int, default=50)
    s2 = sub.add_parser("split")
    s2.add_argument("--cohort_prefix", required=True)
    s2.add_argument("--labels", required=True)
    s2.add_argument("--keep", required=True)
    s2.add_argument("--samples", required=True)
    s2.add_argument("--out_prefix", required=True)
    s2.add_argument("--fold_into", default=None)
    a = ap.parse_args()
    shares(a) if a.mode == "shares" else split(a)


if __name__ == "__main__":
    main()
