# PhyloDist

**PhyloDist** is a Python library for comparing phylogenetic trees using classical and information-theoretic tree distance and similarity measures.

The library is designed to:

- Reproduce results from **TreeDist (R)**
- Support both **rooted** and **unrooted** trees
- Provide both **distance** and **similarity** measures
- Support **normalised scores in [0,1]**
- Provide a clean and extensible Python interface
- Enable validation against published phylogenetic metrics

---

# Features

Supported comparison methods:

| Method | Description | Return Type | TreeDist (R) Equivalent |
|--------|-------------|-------------|--------------------------|
| `rf` | Robinson–Foulds distance | Distance | `RobinsonFoulds(tree1, tree2, normalize = FALSE)` |
| `jrf` | Jaccard Robinson–Foulds distance | Distance | `JaccardRobinsonFoulds(tree1, tree2, normalize = FALSE)` |
| `icrf` | Information-corrected RF | Similarity (TreeDist-equivalent distance also available) | `InfoRobinsonFoulds(tree1, tree2, normalize = FALSE)` |
| `nye` | Nye similarity | Similarity | `NyeSimilarity(tree1, tree2, normalize = FALSE)` |
| `msi` | Matching Split Information | Similarity | `MatchingSplitInfo(tree1, tree2, normalize = FALSE)` |
| `msid` | Matching Split Information Distance | Distance | `MatchingSplitInfoDistance(tree1, tree2, normalize = FALSE)` |
| `spi` | Shared Phylogenetic Information | Similarity | `SharedPhylogeneticInfo(tree1, tree2, normalize = FALSE)` |
| `pid` | Phylogenetic Information Distance | Distance | `PhylogeneticInfoDistance(tree1, tree2, normalize = FALSE)` |
| `mci` | Mutual Clustering Information | Similarity | `MutualClusteringInfo(tree1, tree2, normalize = FALSE)` |
| `cid` | Clustering Information Distance | Distance | `ClusteringInfoDist(tree1, tree2, normalize = FALSE)` |

---

# Installation

## Method 1 — Install directly from GitHub

```bash
pip install git+https://github.com/OrthoFinder/PhyloDist.git
```

## Method 2 - Clone the repository

```bash
git clone https://github.com/OrthoFinder/PhyloDist.git
cd PhyloDist
pip install .
```

## Dependencies

- Python ≥ 3.11
- numpy
- scipy
- ete4

# Basic Usage
```python
from phylodist import PhyloDist, Tree

# Trees must be valid Newick strings or file paths.
tree1 = "tree1.txt"
tree2 = "tree2.txt"
# or 
# tree1 = "((A,B),(C,D));"
# tree2 = "((A,C),(B,D));"

# Create tree objects
t1 = Tree(
    tree1,
    is_rooted=False,   # Set rootedness here
    parser=1           # Parser value follows ete4 conventions
)

t2 = Tree(
    tree2,
    is_rooted=False,
    parser=1
)

# Create PhyloDist object
grf = PhyloDist(t1, t2)

# Compare trees
result = grf.compare(
    method="mci",
    rooted=False,     # Force comparison mode
    normalise=True    # Scale result to [0,1]
)

print(result.score)
```

# Scalability Notes

Unlike **TreeDist (R)**, which is often used on trees of approximately **≤1000 leaves** in practice, **PhyloDist** does not impose a fixed upper limit on tree size.

However, many tree comparison methods implemented here rely on bipartite matching algorithms (Hungarian algorithm), whose computational complexity grows rapidly with tree size.

Typical complexity:

| Component | Complexity |
|----------|-------------|
| Split comparison | O(n²) |
| Matching (Hungarian algorithm) | O(n³) |
| Memory usage | O(n²) |

Where:

- `n` = number of splits (roughly proportional to number of leaves)

---

## Practical Recommendations

For best performance:

- ≤ **500 leaves** — fast
- **500–2000 leaves** — moderate runtime
- ≥ **2000 leaves** — potentially slow
- ≥ **5000 leaves** — large memory usage expected

Performance depends strongly on:

- number of taxa
- tree balance
- chosen metric
- available system memory

---

## Important Note

Large trees are supported, but runtime may increase substantially.  
Users working with large phylogenies are encouraged to:

- test on smaller subsets first
- monitor memory usage
- consider parallel workflows for batch comparisons

# Important Interpretation Notes

Many tree comparison metrics implemented here are information-based. For these metrics:

- Similarity values increase as trees become more similar.
- Distance values increase as trees become more different.
- Raw information-based scores (similarity and distance) depend on tree size.
- These values should not be compared across datasets without normalization.
- Use `normalise=True` when comparing trees of different sizes.

Metrics especially affected include:

- icrf
- msi / msid
- spi / pid
- mci / cid


# Reference
- Martin R Smith, Information theoretic generalized Robinson–Foulds metrics for comparing phylogenetic trees, Bioinformatics, Volume 36, Issue 20, October 2020, Pages 5007–5013, [![DOI](https://img.shields.io/badge/DOI-10.1093%2Fbioinformatics%2Fbtaa614-blue)](https://doi.org/10.1093/bioinformatics/btaa614)