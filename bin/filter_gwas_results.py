#!/usr/bin/env python3
"""
Filter GWAS results and extract significant variants.
"""

import argparse
import gzip
import pandas as pd
from pathlib import Path


def detect_columns(df: pd.DataFrame) -> dict:
    """Detect column names for different GWAS tools."""
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
    # P_JOINT first so Tractor-GENESIS files use the joint test, not a per-ancestry P_<anc>
    p_patterns = ["P_JOINT", "P", "PVALUE", "P_VALUE", "PVAL", "LOG10P", "p.value", "P-value", "Score.pval", "P_META", "P_BOLT_LMM"]
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

    return column_mapping


def filter_gwas(
    input_path: str,
    output_path: str,
    significant_path: str,
    p_threshold: float = 5e-8,
    suggestive_threshold: float = 1e-5,
) -> None:
    """Filter GWAS results."""

    # Read input (handle gzipped files)
    if input_path.endswith(".gz"):
        df = pd.read_csv(input_path, sep="\t", compression="gzip")
    else:
        df = pd.read_csv(input_path, sep="\t")

    # Detect columns
    cols = detect_columns(df)

    if cols["p"] is None:
        raise ValueError("Could not detect P-value column")

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

    # Write outputs
    df_filtered.to_csv(output_path, sep="\t", index=False, compression="gzip")
    df_significant.to_csv(significant_path, sep="\t", index=False)

    print(f"Total variants: {len(df)}")
    print(f"Suggestive (P < {suggestive_threshold}): {len(df_filtered)}")
    print(f"Genome-wide significant (P < {p_threshold}): {len(df_significant)}")


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
    args = parser.parse_args()

    filter_gwas(
        args.input,
        args.output,
        args.significant,
        args.p_threshold,
        args.suggestive_threshold,
    )


if __name__ == "__main__":
    main()
