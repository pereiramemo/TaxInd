# taxind

Compute alpha and phylogenetic diversity metrics from rarefied abundance matrices.

## Installation

### 1. Clone the repository

```bash
git clone <repo-url>
cd taxind
```

### 2. Create and activate the conda environment

```bash
conda env create -f environment.yml
conda activate taxind_env
```

## Modules

### 1.1-div.py — Alpha diversity (richness, Shannon, Simpson, Pielou)

Computes five alpha diversity metrics from a rarefied abundance matrix.

**Metrics:**

| Column | Description |
|--------|-------------|
| `observed_richness` | Number of taxa with abundance > 0 (S_obs) |
| `chao1` | Chao1 estimated total richness based on singletons and doubletons |
| `shannon_diversity` | Shannon diversity index (H'), natural log |
| `simpson_diversity` | Simpson diversity index (D = 1 - sum(p_i²)); ranges from 0 to 1 |
| `pielou_evenness` | Pielou's evenness (J' = H' / ln(S_obs)); NaN when S_obs < 2 |

**Usage:**

```bash
python modules/1.1-div.py -i <abundance_matrix> [-o output.tsv]
```

**Arguments:**

| Argument | Description |
|----------|-------------|
| `-i / --input` | Rarefied abundance matrix (CSV or TSV). Required. |
| `-o / --output` | Output file path. Defaults to stdout. |
| `--sep` | Column separator. Auto-detected from file extension if not provided. |
| `--index-col` | Column index to use as sample labels (default: 0). |

**Example:**

```bash
python modules/1.1-div.py -i examples/abundance_matrix.csv
python modules/1.1-div.py -i examples/abundance_matrix.csv -o results/alpha_div.tsv
```

---

### 1.2-div.py — Phylogenetic diversity (Faith's PD, Phylogenetic Entropy, Rao's Q)

Computes phylogenetic diversity metrics from a rarefied abundance matrix by aligning
sequences with MAFFT, building a phylogenetic tree with FastTree, and computing metrics
using scikit-bio.

**Workflow:**
1. Read abundance matrix (rows = samples, columns = ASVs/sequences)
2. Resolve sequences from a FASTA file or directly from column names
3. Align sequences with MAFFT (`--auto`)
4. Build a phylogenetic tree with FastTree (`-nt`, JC69 model)
5. Root the tree at the midpoint
6. Compute Faith's PD, Phylogenetic Entropy (Hp), and Rao's Quadratic Diversity (Q) per sample

**Metrics:**

| Column | Description |
|--------|-------------|
| `faith_pd` | Faith's phylogenetic diversity: total branch length of the minimum spanning tree connecting all present taxa |
| `hp` | Phylogenetic entropy (Allen et al. 2009): Hp = -∑ ℓ(b) · p(b) · ln(p(b)) over all branches; weights Shannon entropy by evolutionary distinctness |
| `rao_q` | Rao's quadratic diversity (Rao 1982): Q = ∑∑ d_ij · p_i · p_j; mean pairwise phylogenetic distance weighted by relative abundances |

**Usage:**

```bash
python modules/1.2-div.py -i <abundance_matrix> [-f sequences.fasta] [options]
```

**Arguments:**

| Argument | Description |
|----------|-------------|
| `-i / --input` | Rarefied abundance matrix (CSV or TSV). Required. |
| `-f / --fasta` | FASTA file with sequences. Headers must match column names. If not provided, column names are treated as sequences. |
| `--tree-output` | Output path for the Newick tree (default: `tree.nwk`). |
| `-o / --output` | Output TSV for phylogenetic diversity metrics (default: `faith_pd.tsv`). |
| `--alignment-output` | If provided, save the MAFFT alignment as a FASTA file at this path. |
| `--sep` | Column separator. Auto-detected from file extension if not provided. |
| `--index-col` | Column index to use as sample labels (default: 0). |

**Example:**

```bash
python modules/1.2-div.py -i examples/abundance_matrix.csv -f examples/sequences.fasta
python modules/1.2-div.py -i examples/abundance_matrix.csv -f examples/sequences.fasta \
    --tree-output results/tree.nwk \
    -o results/phylo_div.tsv \
    --alignment-output results/aligned.fasta
```

## Input format

A CSV or TSV file where:
- First column: sample identifiers
- Remaining columns: species / OTU / ASV counts (integers)

An example matrix is provided in `examples/abundance_matrix.csv`.

## References

- **Chao1 richness estimator**: Chao, A. (1984). Nonparametric estimation of the number of classes in a population. *Scandinavian Journal of Statistics*, 11(4), 265–270.

- **Phylogenetic entropy (Hp)**: Allen, B., Kon, M., & Bar-Yam, Y. (2009). A new phylogenetic diversity measure generalizing the Shannon index and its application to phyllostomid bats. *The American Naturalist*, 174(2), 236–243. https://doi.org/10.1086/600101

- **Faith's phylogenetic diversity**: Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. *Biological Conservation*, 61(1), 1–10. https://doi.org/10.1016/0006-3207(92)91201-3

- **Rao's quadratic diversity**: Rao, C. R. (1982). Diversity and dissimilarity coefficients: a unified approach. *Theoretical Population Biology*, 21(1), 24–43. https://doi.org/10.1016/0040-5809(82)90004-1
# TaxInd
