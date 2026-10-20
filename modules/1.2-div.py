#!/usr/bin/env python3
"""
Compute Faith's phylogenetic diversity from a rarefied abundance matrix.

Workflow:
  1. Read abundance matrix (rows = samples, columns = ASVs/sequences)
  2. Optionally read sequences from a FASTA file (headers must match column names).
     If no FASTA is provided, column names are treated as the sequences themselves.
  3. Align sequences with MAFFT (--auto)
  4. Build a phylogenetic tree with FastTree (-nt, JC69 model)
  5. Compute Faith's PD, Phylogenetic Entropy (Hp), and Rao's Quadratic Diversity (Q) per sample
  6. Output: Newick tree file + phylogenetic diversity metrics per sample as TSV

Input matrix format (CSV or TSV):
  - Rows: samples
  - Columns: ASVs / sequences
  - First column: sample identifiers
  - Values: integer counts
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from skbio import DNA, TreeNode
from skbio import io as skbio_io
from skbio.diversity.alpha import faith_pd as skbio_faith_pd

###############################################################################
# 2. Define auxiliary functions
###############################################################################


def detect_sep(path: str) -> str:
    if path.endswith(".tsv") or path.endswith(".txt"):
        return "\t"
    return ","


def read_fasta(path: str) -> dict:
    """Read a FASTA file and return {header: sequence}."""
    return {
        str(seq.metadata["id"]): str(seq)
        for seq in skbio_io.read(path, format="fasta", constructor=DNA)
    }


def run_mafft(sequences: dict, tmp_dir: str) -> tuple:
    """
    Align sequences with MAFFT (--auto).

    Uses short internal IDs (ASV-1, ASV-2, ...) in the FASTA to avoid issues
    with long or special-character column names as headers.

    Returns (aligned_fasta_path, id_map) where id_map maps safe_id -> original_id.
    """
    id_map = {f"ASV-{i+1}": orig_id for i, orig_id in enumerate(sequences)}
    reverse_map = {v: k for k, v in id_map.items()}

    input_fasta = Path(tmp_dir) / "input.fasta"
    aligned_fasta = Path(tmp_dir) / "aligned.fasta"

    with open(input_fasta, "w") as fh:
        for orig_id, seq in sequences.items():
            fh.write(f">{reverse_map[orig_id]}\n{seq}\n")

    result = subprocess.run(
        ["mafft", "--auto", str(input_fasta)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT error:\n{result.stderr}")

    with open(aligned_fasta, "w") as fh:
        fh.write(result.stdout)

    return str(aligned_fasta), id_map, result.stdout


def run_fasttree(aligned_fasta: str) -> str:
    """
    Build a nucleotide phylogenetic tree with FastTree (JC69 model).

    Tries 'FastTree' first, then 'fasttree' as a fallback.
    Returns the Newick string.
    """
    result = subprocess.run(
        ["FastTree", "-nt", aligned_fasta],
        capture_output=True,
        text=True,
    )
    if result.returncode == 0:
        return result.stdout.strip()
    raise RuntimeError(f"FastTree error:\n{result.stderr}")


def compute_faith_pd(matrix: pd.DataFrame, tree: TreeNode) -> pd.Series:
    """Compute Faith's PD for each sample using scikit-bio."""
    otu_ids = list(matrix.columns)
    results = {}
    for sample, row in matrix.iterrows():
        counts = row.values.astype(int)
        results[sample] = skbio_faith_pd(counts, otu_ids, tree)
    return pd.Series(results, name="faith_pd")


def phylogenetic_entropy(counts: np.ndarray, otu_ids: list, tree: TreeNode) -> float:
    """
    Phylogenetic entropy (Hp) from Allen et al. 2009.

    Hp = -sum over all branches b of T: l(b) * p(b) * ln(p(b))

    where:
      l(b) = length of branch b
      p(b) = proportion of individuals descending from branch b

    When all species are equally distinct (uniform branch lengths),
    Hp reduces to twice the Shannon index.
    """
    total = counts.sum()
    if total == 0:
        return 0.0
    abund = {otu_id: count / total for otu_id, count in zip(otu_ids, counts)}
    hp = 0.0
    for node in tree.traverse():
        if node.is_root():
            continue
        length = node.length
        if length is None or length == 0:
            continue
        p_b = sum(abund.get(tip.name, 0.0) for tip in node.tips(include_self=True))
        if p_b > 0:
            hp += length * p_b * np.log(p_b)
    return -hp


def compute_hp(matrix: pd.DataFrame, tree: TreeNode) -> pd.Series:
    """Compute phylogenetic entropy (Hp) for each sample."""
    otu_ids = list(matrix.columns)
    results = {}
    for sample, row in matrix.iterrows():
        counts = row.values.astype(float)
        results[sample] = phylogenetic_entropy(counts, otu_ids, tree)
    return pd.Series(results, name="hp")


def rao_quadratic_diversity(counts: np.ndarray, D: np.ndarray) -> float:
    """
    Rao's quadratic diversity (Q) from Rao (1982).

    Q = sum_i sum_j d_ij * p_i * p_j = p^T * D * p

    where:
      d_ij = phylogenetic distance between species i and j
      p_i, p_j = relative abundances of species i and j
    """
    total = counts.sum()
    if total == 0:
        return 0.0
    proportions = counts / total
    return float(proportions @ D @ proportions)


def compute_rao(matrix: pd.DataFrame, tree: TreeNode) -> pd.Series:
    """Compute Rao's quadratic diversity (Q) for each sample."""
    otu_ids = list(matrix.columns)

    # Compute pairwise tip distances once for all samples
    dist_matrix = tree.tip_tip_distances()
    dm_ids = list(dist_matrix.ids)
    idx = [dm_ids.index(otu_id) for otu_id in otu_ids]
    D = dist_matrix.data[np.ix_(idx, idx)]

    results = {}
    for sample, row in matrix.iterrows():
        counts = row.values.astype(float)
        results[sample] = rao_quadratic_diversity(counts, D)
    return pd.Series(results, name="rao_q")


def parse_args(): 
    parser = argparse.ArgumentParser(
        description="Compute Faith's phylogenetic diversity from a rarefied abundance matrix."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Rarefied abundance matrix (CSV or TSV). Rows = samples, columns = ASVs/sequences.",
    )
    parser.add_argument(
        "-f", "--fasta",
        default=None,
        help=(
            "FASTA file with sequences. Headers must match column names of the abundance matrix. "
            "If not provided, column names are treated as the sequences themselves."
        ),
    )
    parser.add_argument(
        "--tree-output",
        default="tree.nwk",
        help="Output path for the Newick tree (default: tree.nwk).",
    )
    parser.add_argument(
        "-o", "--output",
        default="faith_pd.tsv",
        help="Output TSV file for Faith PD values (default: faith_pd.tsv).",
    )
    parser.add_argument(
        "--alignment-output",
        default=None,
        help="If provided, save the MAFFT alignment as a FASTA file at this path.",
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

    # Resolve sequences
    if args.fasta:
        try:
            all_seqs = read_fasta(args.fasta)
        except FileNotFoundError:
            sys.exit(f"Error: FASTA file not found: {args.fasta}")

        missing = [col for col in matrix.columns if col not in all_seqs]
        if missing:
            sys.exit(
                f"Error: {len(missing)} column(s) not found in FASTA. "
                f"First missing: {missing[:5]}"
            )
        sequences = {col: all_seqs[col] for col in matrix.columns}
    else:
        # Column names are the sequences themselves (typical ASV table format)
        sequences = {col: col for col in matrix.columns}

    # Align with MAFFT and build tree with FastTree
    with tempfile.TemporaryDirectory() as tmp_dir:
        aligned_fasta, id_map, aligned_content = run_mafft(sequences, tmp_dir)
        newick_str = run_fasttree(aligned_fasta)

    # Optionally save the alignment with original headers
    if args.alignment_output:
        with open(args.alignment_output, "w") as fh:
            for line in aligned_content.splitlines():
                if line.startswith(">"):
                    safe_id = line[1:]
                    fh.write(f">{id_map[safe_id]}\n")
                else:
                    fh.write(f"{line}\n")
        print(f"Alignment written to: {args.alignment_output}")

    # Rename tree tips from safe IDs back to original column names
    tree = TreeNode.read([newick_str])
    for tip in tree.tips():
        tip.name = id_map[tip.name]

    # Root the tree at the midpoint (required by faith_pd)
    tree = tree.root_at_midpoint()

    # Write Newick tree
    tree.write(args.tree_output)
    print(f"Tree written to:     {args.tree_output}")

    # Compute and write phylogenetic diversity metrics
    results = pd.concat([
        compute_faith_pd(matrix, tree),
        compute_hp(matrix, tree),
        compute_rao(matrix, tree),
    ], axis=1)
    results.to_csv(args.output, sep="\t", float_format="%.6f")
    print(f"Phylogenetic diversity written to: {args.output}")

###############################################################################
# 4. Run Main
###############################################################################

if __name__ == "__main__":
    main()
