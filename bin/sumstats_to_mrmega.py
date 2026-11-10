#!/usr/bin/env python3
"""
sumstats_to_mrmega.py: one stratum's Tractor-GENESIS sumstats -> MR-MEGA input (CP-19)

    sumstats_to_mrmega.py --input results/gwas/os/EUR.sumstats.gz --output EUR.mrmega_in.tsv
                          [--variants HANDOFF/genotypes/cohort.hardcall.maf01.bim]

MR-MEGA pools ONE effect per variant and cohort. Each stratum's is the total-dose 1-df test
(BETA_TOTAL, SE_TOTAL: the adapter's standard single-effect test; docs/SUMSTATS_SCHEMA.md).
Output columns, MR-MEGA's names (BETA / SE: run MR-MEGA with --qt):
    MARKERNAME CHROMOSOME POSITION EA NEA EAF N BETA SE
EA is the ALT allele and NEA the REF allele: every effect in this pipeline is per ALT allele.
Rows without BETA_TOTAL / SE_TOTAL (below the MAC gate, or a failed test) are left out.
The file starts with "## key=value" lines naming the source columns (beta_source=BETA_TOTAL);
MR-MEGA itself must be given a copy without them (slurm/submit_mrmega.sh strips them).

--variants: the handoff bundle's variant table (.bim, IDs chr:pos:ref:alt). Every output row whose
MARKERNAME is in it must have NEA == its ref and EA == its alt; a row the other way round stops the
run (the round-trip orientation check).

    sumstats_to_mrmega.py --fixed_effect --inputs A.mrmega_in.tsv,B.mrmega_in.tsv --output pooled.tsv

K = 2 strata: an inverse-variance fixed-effect meta-analysis of the two (MR-MEGA needs K >= 3 for
one axis); the header says "pc_axes=0: fixed-effect meta-analysis, not MR-MEGA".
"""

import argparse
import gzip
import math
import sys

import pandas as pd


def meta_lines(path):
    opener = gzip.open if path.endswith(".gz") else open
    out = []
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            out.append(line.rstrip("\n"))
    return out


def read_table(path):
    return pd.read_csv(path, sep="\t", skiprows=len(meta_lines(path)),
                       compression="gzip" if path.endswith(".gz") else None,
                       dtype={"CHR": str, "ID": str, "REF": str, "ALT": str, "MARKERNAME": str, "EA": str, "NEA": str})


def write_with_meta(df, path, meta):
    with open(path, "w") as fh:
        for m in meta:
            fh.write(m + "\n")
        df.to_csv(fh, sep="\t", index=False, na_rep="NA")


def convert(args):
    ss = read_table(args.input)
    need = ["CHR", "POS", "ID", "REF", "ALT", "N", "EAF", "BETA_TOTAL", "SE_TOTAL"]
    missing = [c for c in need if c not in ss.columns]
    if missing:
        sys.exit(f"ERROR: {args.input} lacks {', '.join(missing)} (written by the adapter since CP-19)")
    ok = ss["BETA_TOTAL"].notna() & ss["SE_TOTAL"].notna() & (ss["SE_TOTAL"] > 0)
    out = pd.DataFrame({
        "MARKERNAME": ss.loc[ok, "ID"],
        "CHROMOSOME": ss.loc[ok, "CHR"].str.replace("^chr", "", regex=True),
        "POSITION": ss.loc[ok, "POS"],
        "EA": ss.loc[ok, "ALT"],
        "NEA": ss.loc[ok, "REF"],
        "EAF": ss.loc[ok, "EAF"],
        "N": ss.loc[ok, "N"],
        "BETA": ss.loc[ok, "BETA_TOTAL"],
        "SE": ss.loc[ok, "SE_TOTAL"],
    })
    checked = 0
    if args.variants:
        bim = pd.read_csv(args.variants, sep=r"\s+", header=None, usecols=[1], names=["ID"], dtype=str)
        parts = bim["ID"].str.split(":", expand=True)
        if parts.shape[1] < 4:
            sys.exit(f"ERROR: {args.variants}: variant IDs are not chr:pos:ref:alt (e.g. {bim['ID'].iloc[0]})")
        # keyed by chromosome, position and the unordered allele pair, so a variant written the other
        # way round (ref and alt swapped) is found, and caught
        def key(c, pos, a1, a2):
            return f"{str(c).replace('chr', '')}:{pos}:{':'.join(sorted([str(a1), str(a2)]))}"
        table = {key(c, p_, r, a): (r, a) for c, p_, r, a in zip(parts[0], parts[1], parts[2], parts[3])}
        keys = [key(c, p_, ea, nea) for c, p_, ea, nea in zip(out["CHROMOSOME"], out["POSITION"], out["EA"], out["NEA"])]
        known = pd.Series([k in table for k in keys], index=out.index)
        bad = pd.Series([k in table and (nea != table[k][0] or ea != table[k][1])
                         for k, ea, nea in zip(keys, out["EA"], out["NEA"])], index=out.index)
        if bad.any():
            i = bad.idxmax(); k = keys[out.index.get_loc(i)]; r = out.loc[i]
            sys.exit(f"ERROR: orientation: {int(bad.sum())} row(s) do not match {args.variants}, e.g. {r['MARKERNAME']}: "
                     f"EA {r['EA']} / NEA {r['NEA']}, the variant table has ref {table[k][0]} / alt {table[k][1]}")
        checked = int(known.sum())
        print(f"orientation: {checked} of {len(out)} rows checked against {args.variants}: EA = alt, NEA = ref on all; "
              f"{len(out) - checked} not in the variant table")
    src = {m.split("=", 1)[0][3:].strip(): m.split("=", 1)[1] for m in meta_lines(args.input) if "=" in m}
    meta = ["## beta_source=BETA_TOTAL", "## se_source=SE_TOTAL", "## ea=ALT", "## nea=REF",
            f"## input={args.input}", f"## stratum={src.get('stratum', 'unknown')}",
            f"## orientation_checked={checked}" + (f" (against {args.variants})" if args.variants else " (no --variants)")]
    write_with_meta(out, args.output, meta)
    print(f"{args.output}: {len(out)} variants (of {len(ss)}; {int((~ok).sum())} without BETA_TOTAL / SE_TOTAL left out)")


def fixed_effect(args):
    files = args.inputs.split(",")
    if len(files) != 2:
        sys.exit("ERROR: --fixed_effect is for K = 2 strata: give exactly two --inputs")
    tabs = [read_table(f) for f in files]
    a, b = tabs
    m = a.merge(b, on="MARKERNAME", how="outer", suffixes=("_1", "_2"))
    both = m["BETA_1"].notna() & m["BETA_2"].notna()
    flip = both & ((m["EA_1"] != m["EA_2"]) | (m["NEA_1"] != m["NEA_2"]))
    if flip.any():
        sys.exit(f"ERROR: {int(flip.sum())} variant(s) with other EA / NEA in the two strata, e.g. {m.loc[flip, 'MARKERNAME'].iloc[0]}")
    w1 = 1 / m["SE_1"] ** 2
    w2 = 1 / m["SE_2"] ** 2
    w1 = w1.fillna(0)
    w2 = w2.fillna(0)
    sw = w1 + w2
    beta = (m["BETA_1"].fillna(0) * w1 + m["BETA_2"].fillna(0) * w2) / sw
    se = 1 / sw ** 0.5
    z = beta / se
    p = z.abs().map(lambda x: math.erfc(x / math.sqrt(2)) if pd.notna(x) else float("nan"))
    take = lambda c: m[c + "_1"].combine_first(m[c + "_2"])
    out = pd.DataFrame({"MARKERNAME": m["MARKERNAME"], "CHROMOSOME": take("CHROMOSOME"), "POSITION": take("POSITION"),
                        "EA": take("EA"), "NEA": take("NEA"), "N": m["N_1"].fillna(0) + m["N_2"].fillna(0),
                        "BETA": beta, "SE": se, "P": p, "N_STRATA": both.astype(int) + 1})
    out = out.sort_values("P", na_position="last")
    meta = ["## pc_axes=0: fixed-effect meta-analysis, not MR-MEGA",
            "## method=inverse-variance fixed effect of BETA_TOTAL across K=2 strata",
            f"## strata={args.strata or ','.join(files)}", "## beta_source=BETA_TOTAL"]
    write_with_meta(out, args.output, meta)
    print(f"{args.output}: {len(out)} variants, {int(both.sum())} in both strata (fixed effect, K=2)")


def main():
    ap = argparse.ArgumentParser(description="Tractor-GENESIS sumstats -> MR-MEGA input; or a K=2 fixed-effect meta-analysis")
    ap.add_argument("--input", help="one stratum's sumstats (results/gwas/<trait>/<S>.sumstats.gz)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--variants", default=None, help="the bundle's .bim (IDs chr:pos:ref:alt) for the orientation check")
    ap.add_argument("--fixed_effect", action="store_true", help="pool two MR-MEGA inputs (K = 2)")
    ap.add_argument("--inputs", default=None, help="with --fixed_effect: two MR-MEGA input files, comma-separated")
    ap.add_argument("--strata", default=None, help="with --fixed_effect: the two strata, for the header")
    args = ap.parse_args()
    if args.fixed_effect:
        fixed_effect(args)
    elif args.input:
        convert(args)
    else:
        ap.error("--input (or --fixed_effect --inputs A,B) is required")


if __name__ == "__main__":
    main()
