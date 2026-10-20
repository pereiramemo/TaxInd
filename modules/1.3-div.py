#!/usr/bin/env python3
"""
Compute Bray-Curtis dissimilarity between samples from a rarefied abundance matrix.

Workflow:
  1. Read abundance matrix (rows = samples, columns = ASVs/species)
  2. Compute pairwise Bray-Curtis dissimilarity between all samples
  3. Output: pairwise dissimilarity matrix as TSV

Formula:
  BC(i, j) = sum(|x_ik - x_jk|) / (sum(x_ik) + sum(x_jk))

  where x_ik is the count of taxon k in sample i.
  Values range from 0 (identical) to 1 (no taxa in common).

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


def detect_sep(path: str) -> str:
    if path.endswith(".tsv") or path.endswith(".txt"):
        return "\t"
    return ","


def bray_curtis(u: np.ndarray, v: np.ndarray) -> float:
    """
    Bray-Curtis dissimilarity between two samples.

    BC(i, j) = sum(|u_k - v_k|) / (sum(u_k) + sum(v_k))

    Returns 0.0 if both samples are empty.
    """
    denom = u.sum() + v.sum()
    if denom == 0:
        return 0.0
    return float(np.sum(np.abs(u - v)) / denom)


def compute_bray_curtis(matrix: pd.DataFrame) -> pd.DataFrame:
    """
    Compute pairwise Bray-Curtis dissimilarity between all samples.

    Returns a symmetric DataFrame (samples x samples) with zeros on the diagonal.
    """
    samples = list(matrix.index)
    n = len(samples)
    data = matrix.values.astype(float)

    result = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            bc = bray_curtis(data[i], data[j])
            result[i, j] = bc
            result[j, i] = bc

    return pd.DataFrame(result, index=samples, columns=samples)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute pairwise Bray-Curtis dissimilarity from a rarefied abundance matrix."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Rarefied abundance matrix (CSV or TSV). Rows = samples, columns = species/ASVs.",
    )
    parser.add_argument(
        "-o", "--output",
        default="bray_curtis.tsv",
        help="Output TSV file for the pairwise dissimilarity matrix (default: bray_curtis.tsv).",
    )
    parser.add_argument(
        "--sep",
        default=None,
        help="Column separator. Auto-detected from file extension if not provided.",
    )
    parser.add_argument(
        "--index-col",
        default=0,
        type=int,
        help="Column index to use as row (sample) labels (default: 0).",
    )
    return parser.parse_args()

###############################################################################
# 3. Define the Main function
###############################################################################


def main():
    args = parse_args()
    sep = args.sep if args.sep is not None else detect_sep(args.input)

    # Load abundance matrix
    try:
        matrix = pd.read_csv(args.input, sep=sep, index_col=args.index_col)
    except FileNotFoundError:
        sys.exit(f"Error: file not found: {args.input}")
    except Exception as exc:
        sys.exit(f"Error reading input file: {exc}")

    if matrix.empty:
        sys.exit("Error: input matrix is empty.")

    if len(matrix) < 2:
        sys.exit("Error: at least two samples are required to compute pairwise dissimilarity.")

    # Compute and write Bray-Curtis dissimilarity matrix
    result = compute_bray_curtis(matrix)
    result.to_csv(args.output, sep="\t", float_format="%.6f")
    print(f"Bray-Curtis dissimilarity matrix written to: {args.output}")

###############################################################################
# 4. Run Main
###############################################################################

if __name__ == "__main__":
    main()
