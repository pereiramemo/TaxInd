#!/usr/bin/env python3
"""
Compute alpha diversity metrics from a rarefied abundance matrix.

Metrics:
  - Observed richness (S_obs): number of species with count > 0
  - Chao1 richness (S_chao1): estimated total richness based on singletons and doubletons
  - Shannon diversity (H'): -sum(p_i * ln(p_i))
  - Simpson diversity (D): 1 - sum(p_i^2)
  - Pielou's evenness (J'): H' / ln(S)

Input matrix format (CSV or TSV):
  - Rows: samples
  - Columns: species / OTUs / ASVs
  - First column: sample identifiers
  - Values: integer counts
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import sys

import numpy as np
import pandas as pd

###############################################################################
# 2. Define auxiliary functions
###############################################################################

def shannon_diversity(counts: np.ndarray) -> float:
    """Shannon diversity index (H') using natural log."""
    counts = counts[counts > 0]
    if counts.size == 0:
        return 0.0
    proportions = counts / counts.sum()
    return float(-np.sum(proportions * np.log(proportions)))


def simpson_diversity(counts: np.ndarray) -> float:
    """Simpson diversity index (D) = 1 - sum(p_i^2). Ranges from 0 (no diversity) to 1 (max diversity)."""
    counts = counts[counts > 0]
    if counts.size == 0:
        return 0.0
    proportions = counts / counts.sum()
    return float(1 - np.sum(proportions ** 2))


def observed_richness(counts: np.ndarray) -> int:
    """Observed species richness (S_obs): number of taxa with at least one count."""
    return int(np.sum(counts > 0))


def chao1(counts: np.ndarray) -> float:
    """
    Chao1 richness estimator.

    S_chao1 = S_obs + (f1^2 / (2 * f2))

    where f1 = number of singletons (count == 1)
          f2 = number of doubletons (count == 2)

    When f2 == 0: S_chao1 = S_obs + f1 * (f1 - 1) / 2 (bias-corrected form).
    """
    counts = counts[counts > 0].astype(int)
    s_obs = counts.size
    f1 = int(np.sum(counts == 1))
    f2 = int(np.sum(counts == 2))
    if f2 > 0:
        return s_obs + (f1 ** 2) / (2 * f2)
    else:
        return s_obs + f1 * (f1 - 1) / 2


def pielou(h: float, s: int) -> float:
    """Pielou's evenness (J') = H' / ln(S). Returns NaN when S < 2."""
    if s < 2:
        return float("nan")
    return h / np.log(s)


def compute_diversity(matrix: pd.DataFrame) -> pd.DataFrame:
    results = []
    for sample, row in matrix.iterrows():
        counts = row.values.astype(float)
        s_obs = observed_richness(counts)
        s_chao1 = chao1(counts)
        h = shannon_diversity(counts)
        d = simpson_diversity(counts)
        j = pielou(h, s_obs)
        results.append({
            "sample": sample,
            "observed_richness": s_obs,
            "chao1": s_chao1,
            "shannon_diversity": h,
            "simpson_diversity": d,
            "pielou_evenness": j,
        })
    return pd.DataFrame(results).set_index("sample")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute alpha diversity metrics from a rarefied abundance matrix."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Rarefied abundance matrix (CSV or TSV). Rows = samples, columns = species.",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Output file path. Defaults to stdout.",
    )
    parser.add_argument(
        "--sep",
        default=None,
        help="Column separator. Auto-detected from file extension if not provided (csv -> ',', tsv/txt -> '\\t').",
    )
    parser.add_argument(
        "--index-col",
        default=0,
        type=int,
        help="Column index to use as row (sample) labels (default: 0).",
    )
    return parser.parse_args()


def detect_sep(path: str) -> str:
    if path.endswith(".tsv") or path.endswith(".txt"):
        return "\t"
    return ","

###############################################################################
# 3. Define the Main function
###############################################################################

def main():
    args = parse_args()
    sep = args.sep if args.sep is not None else detect_sep(args.input)

    try:
        matrix = pd.read_csv(args.input, sep=sep, index_col=args.index_col)
    except FileNotFoundError:
        sys.exit(f"Error: file not found: {args.input}")
    except Exception as exc:
        sys.exit(f"Error reading input file: {exc}")

    if matrix.empty:
        sys.exit("Error: input matrix is empty.")

    results = compute_diversity(matrix)

    out_str = results.to_csv(sep="\t", float_format="%.6f", na_rep="NaN")

    if args.output:
        with open(args.output, "w") as fh:
            fh.write(out_str)
            print(f"Diversity metrics written to: {args.output}")
    else:
        sys.stdout.write(out_str)

###############################################################################
# 4. Run Main
###############################################################################

if __name__ == "__main__":
    main()
