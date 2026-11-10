#!/usr/bin/env python3
"""
Filter GWAS results and extract significant variants.

Outputs: the suggestive variants (--output, with the LOCUS each belongs to), the genome-wide
significant ones (--significant), and a locus table (--loci; M-05). Loci are distance-based, not
LD-based: per stratum (STRATUM column, if any), the variant with the lowest P is an index variant
and every not-yet-assigned suggestive variant within --window-kb of it on its chromosome joins its
locus; repeat with the lowest remaining P. The locus table says so in its "##" header lines.
"""

import argparse
import gzip
import pandas as pd
from pathlib import Path


def count_meta_lines(path: str) -> int:
    """Leading "## key=value" lines of a sumstats table (e.g. "## covariates=..."), skipped on read."""
    opener = gzip.open if path.endswith(".gz") else open
    n = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            n += 1
    return n


def detect_columns(df: pd.DataFrame, p_col: str = None) -> dict:
    """Detect column names for different GWAS tools.

    The P column: --p-col when given (exactly, else ignoring case; an error when absent); without
    it P_JOINT (the Tractor joint test) comes first, then the usual single-P names (CP-14).
    """
    column_mapping = {
        "chr": None,
        "pos": None,
        "snp": None,
        "p": None,
        "beta": None,
        "se": None,
        "af": None,
    }

    # Common column name patterns
    chr_patterns = ["CHR", "CHROM", "#CHROM", "chromosome", "chr"]
    pos_patterns = ["POS", "BP", "GENPOS", "position", "pos"]
    snp_patterns = ["SNP", "ID", "RSID", "variant_id", "MarkerName"]
    p_patterns = ["P_JOINT", "P", "PVALUE", "P_VALUE", "LOG10P", "p.value", "P-value"]
    beta_patterns = ["BETA", "EFFECT", "B", "beta", "Effect"]
    se_patterns = ["SE", "STDERR", "se", "StdErr"]
    af_patterns = ["A1FREQ", "AF", "MAF", "EAF", "Freq", "freq"]

    for col in df.columns:
        col_upper = col.upper()
        if column_mapping["chr"] is None and any(
            p.upper() == col_upper for p in chr_patterns
        ):
            column_mapping["chr"] = col
        if column_mapping["pos"] is None and any(
            p.upper() == col_upper for p in pos_patterns
        ):
            column_mapping["pos"] = col
        if column_mapping["snp"] is None and any(
            p.upper() == col_upper for p in snp_patterns
        ):
            column_mapping["snp"] = col
        if column_mapping["p"] is None and any(
            p.upper() == col_upper for p in p_patterns
        ):
            column_mapping["p"] = col
        if column_mapping["beta"] is None and any(
            p.upper() == col_upper for p in beta_patterns
        ):
            column_mapping["beta"] = col
        if column_mapping["se"] is None and any(
            p.upper() == col_upper for p in se_patterns
        ):
            column_mapping["se"] = col
        if column_mapping["af"] is None and any(
            p.upper() == col_upper for p in af_patterns
        ):
            column_mapping["af"] = col

    # P column by priority, not by position in the file
    upper = {}
    for col in df.columns:
        upper.setdefault(col.upper(), col)
    if p_col:
        if p_col in df.columns:
            column_mapping["p"] = p_col
        elif p_col.upper() in upper:
            column_mapping["p"] = upper[p_col.upper()]
        else:
            p_like = [c for c in df.columns if c.upper().startswith("P_") or c.upper() in ("P", "PVALUE")]
            raise ValueError(f"--p-col {p_col} is not a column (P columns: {', '.join(p_like)})")
    else:
        column_mapping["p"] = next((upper[x.upper()] for x in p_patterns if x.upper() in upper), None)

    return column_mapping


def collapse_loci(df, chr_col, pos_col, id_col, p_col, window_bp, stratum_col=None):
    """Distance-based loci (not LD-based). Returns (df with LOCUS / LOCUS_INDEX / IS_INDEX, loci)."""
    df = df.copy()
    df["LOCUS"] = ""
    df["LOCUS_INDEX"] = ""
    df["IS_INDEX"] = False
    rows = []
    groups = df.groupby(stratum_col, sort=True) if stratum_col else [(None, df)]
    for stratum, g in groups:
        g = g.sort_values(p_col, kind="mergesort")
        left = list(g.index)
        k = 0
        while left:
            i = left[0]
            k += 1
            chrom, pos = str(df.at[i, chr_col]), int(df.at[i, pos_col])
            members = [j for j in left if str(df.at[j, chr_col]) == chrom and abs(int(df.at[j, pos_col]) - pos) <= window_bp]
            locus = f"{stratum}:{k}" if stratum_col else str(k)
            index_id = str(df.at[i, id_col]) if id_col else f"{chrom}:{pos}"
            df.loc[members, "LOCUS"] = locus
            df.loc[members, "LOCUS_INDEX"] = index_id
            df.at[i, "IS_INDEX"] = True
            m_pos = [int(df.at[j, pos_col]) for j in members]
            row = {} if not stratum_col else {stratum_col: stratum}
            row.update({"LOCUS": locus, chr_col: df.at[i, chr_col], "INDEX_ID": index_id, "INDEX_" + pos_col: pos,
                        "START": min(m_pos), "END": max(m_pos), "N_VARIANTS": len(members), p_col: df.at[i, p_col]})
            rows.append(row)
            taken = set(members)
            left = [j for j in left if j not in taken]
    loci = pd.DataFrame(rows)
    return df, loci


def filter_gwas(
    input_path: str,
    output_path: str,
    significant_path: str,
    p_threshold: float = 5e-8,
    suggestive_threshold: float = 1e-5,
    p_col_requested: str = None,
    loci_path: str = None,
    window_kb: float = 500,
) -> None:
    """Filter GWAS results."""

    # Read input (handle gzipped files)
    skip = count_meta_lines(input_path)
    if input_path.endswith(".gz"):
        df = pd.read_csv(input_path, sep="\t", compression="gzip", skiprows=skip)
    else:
        df = pd.read_csv(input_path, sep="\t", skiprows=skip)

    # Detect columns
    cols = detect_columns(df, p_col_requested)

    if cols["p"] is None:
        raise ValueError("Could not detect P-value column (looked for P_JOINT, P, PVALUE, ...; or pass --p-col)")
    print(f"P column: {cols['p']}")

    # Convert P-values if needed (handle LOG10P)
    p_col = cols["p"]
    if "LOG10" in p_col.upper():
        df["P_converted"] = 10 ** (-df[p_col])
        p_col = "P_converted"

    # Filter by suggestive threshold for output
    df_filtered = df[df[p_col] <= suggestive_threshold].copy()

    # Extract genome-wide significant
    df_significant = df[df[p_col] <= p_threshold].copy()

    # Sort by P-value
    df_filtered = df_filtered.sort_values(p_col)
    df_significant = df_significant.sort_values(p_col)

    # Distance-based loci over the suggestive variants (M-05)
    stratum_col = "STRATUM" if "STRATUM" in df.columns else None
    window_bp = int(round(window_kb * 1000))
    if cols["chr"] is None or cols["pos"] is None:
        raise ValueError("Could not detect the chromosome / position columns needed for the locus table")
    df_filtered, loci = collapse_loci(df_filtered, cols["chr"], cols["pos"], cols["snp"], p_col, window_bp, stratum_col)

    # Write outputs
    df_filtered.to_csv(output_path, sep="\t", index=False, compression="gzip")
    df_significant.to_csv(significant_path, sep="\t", index=False)
    if loci_path is None:
        loci_path = (significant_path[: -len(".significant.tsv")] if significant_path.endswith(".significant.tsv")
                     else significant_path) + ".loci.tsv"
    with open(loci_path, "w") as fh:
        fh.write(f"## method=distance-based, not LD-based: a locus is every suggestive variant within {window_kb:g} kb "
                 f"of its index variant on the same chromosome; index variants taken lowest P first"
                 f"{', per stratum' if stratum_col else ''}\n")
        fh.write(f"## window_kb={window_kb:g}\n## p_col={p_col}\n## p_threshold={suggestive_threshold:g}\n")
        if len(loci):
            loci.to_csv(fh, sep="\t", index=False)
        else:
            fh.write("\t".join(([stratum_col] if stratum_col else []) + ["LOCUS", cols["chr"], "INDEX_ID", "INDEX_" + cols["pos"],
                                "START", "END", "N_VARIANTS", p_col]) + "\n")

    print(f"Total variants: {len(df)}")
    print(f"Suggestive (P < {suggestive_threshold}): {len(df_filtered)}")
    print(f"Genome-wide significant (P < {p_threshold}): {len(df_significant)}")
    print(f"Loci (distance-based, {window_kb:g} kb window, not LD-based): {len(loci)} -> {loci_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Filter GWAS results"
    )
    parser.add_argument("--input", required=True, help="Input summary statistics")
    parser.add_argument("--output", required=True, help="Filtered output")
    parser.add_argument("--significant", required=True, help="Significant variants")
    parser.add_argument(
        "--p-threshold", type=float, default=5e-8, help="P-value threshold"
    )
    parser.add_argument(
        "--suggestive-threshold",
        type=float,
        default=1e-5,
        help="Suggestive P-value threshold",
    )
    parser.add_argument(
        "--loci",
        default=None,
        help="Locus table [default: --significant with .significant.tsv replaced by .loci.tsv]",
    )
    parser.add_argument(
        "--window-kb",
        type=float,
        default=500,
        help="Distance-based locus window either side of each index variant, in kb [default: 500]",
    )
    parser.add_argument(
        "--p-col",
        default=None,
        help="P column to filter on [default: P_JOINT, else P / PVALUE / P_VALUE / LOG10P]",
    )
    args = parser.parse_args()

    filter_gwas(
        args.input,
        args.output,
        args.significant,
        args.p_threshold,
        args.suggestive_threshold,
        args.p_col,
        args.loci,
        args.window_kb,
    )


if __name__ == "__main__":
    main()
