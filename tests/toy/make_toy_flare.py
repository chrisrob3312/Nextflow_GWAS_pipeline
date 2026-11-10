#!/usr/bin/env python3
"""
make_toy_flare.py: a toy FLARE local-ancestry VCF (what FLARE writes per chromosome)

    make_toy_flare.py <out.anc.vcf.gz> <samples.txt|N> <n_variants> <panel_order> [chrom] [seed] [pure=S1,S2,...]

panel_order: the labels in FLARE's index order (the panel file's order of first appearance), e.g.
EUR,AFR,AMR -> AN1 / AN2 = 0 for EUR, 1 for AFR, 2 for AMR. Each haplotype carries tracts of one
ancestry; GT is phased. pure=...: samples whose both haplotypes are the first label everywhere
(a one-ancestry stratum). Standard library only.
"""
import gzip
import random
import sys

a = sys.argv[1:]
if len(a) < 4:
    sys.exit(__doc__)
out, who, m, order = a[0], a[1], int(a[2]), a[3].split(",")
chrom = a[4] if len(a) > 4 and "=" not in a[4] else "chr22"
seed = int(a[5]) if len(a) > 5 and "=" not in a[5] else 1
pure = set()
for x in a[4:]:
    if x.startswith("pure="):
        pure = set(x[5:].split(","))
random.seed(seed)
try:
    samples = [l.strip() for l in open(who) if l.strip()]
except FileNotFoundError:
    samples = [f"S{i:04d}" for i in range(1, int(who) + 1)]
k = len(order)
freq = [[min(0.9, max(0.05, random.uniform(0.1, 0.6) + random.gauss(0, 0.1))) for _ in range(k)] for _ in range(m)]


def tracts(sample):
    if sample in pure:
        return [0] * m
    la, pos = [], 0
    while pos < m:
        anc = random.randrange(k)
        n = 1 + int(random.expovariate(1 / 25))
        la.extend([anc] * n)
        pos += n
    return la[:m]


hap = {s: (tracts(s), tracts(s)) for s in samples}
with gzip.open(out, "wt") as fh:
    fh.write("##fileformat=VCFv4.2\n")
    fh.write("##ANCESTRY=<" + ",".join(f"{lab}={i}" for i, lab in enumerate(order)) + ">\n")
    fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fh.write('##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="Ancestry of first haplotype">\n')
    fh.write('##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Ancestry of second haplotype">\n')
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
    for v in range(m):
        pos = 10000 + 1000 * v
        cells = []
        for s in samples:
            a1, a2 = hap[s][0][v], hap[s][1][v]
            g1 = int(random.random() < freq[v][a1])
            g2 = int(random.random() < freq[v][a2])
            cells.append(f"{g1}|{g2}:{a1}:{a2}")
        fh.write(f"{chrom}\t{pos}\t{chrom}:{pos}:A:G\tA\tG\t.\tPASS\t.\tGT:AN1:AN2\t" + "\t".join(cells) + "\n")
print(f"toy FLARE VCF {out}: {len(samples)} samples, {m} variants, indices {', '.join(f'{i}={l}' for i, l in enumerate(order))}")
