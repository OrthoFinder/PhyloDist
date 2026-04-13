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

## Features

Supported comparison methods:

| Method | Description | Return Type | TreeDist (R) Equivalent |
|--------|-------------|-------------|--------------------------|
| `rf` | Robinson–Foulds distance | Distance | `RobinsonFoulds(tree1, tree2, normalise = FALSE)` |
| `jrf` | Jaccard Robinson–Foulds distance | Distance | `JaccardRobinsonFoulds(tree1, tree2, normalise = FALSE)` |
| `icrf` | Information-corrected RF | Similarity (TreeDist-equivalent distance also available) | `InfoRobinsonFoulds(tree1, tree2, normalise = FALSE)` |
| `nye` | Nye similarity | Similarity | `NyeSimilarity(tree1, tree2, normalise = FALSE)` |
| `msi` | Matching Split Information | Similarity | `MatchingSplitInfo(tree1, tree2, normalise = FALSE)` |
| `msid` | Matching Split Information Distance | Distance | `MatchingSplitInfoDistance(tree1, tree2, normalise = FALSE)` |
| `spi` | Shared Phylogenetic Information | Similarity | `SharedPhylogeneticInfo(tree1, tree2, normalise = FALSE)` |
| `pid` | Phylogenetic Information Distance | Distance | `PhylogeneticInfoDistance(tree1, tree2, normalise = FALSE)` |
| `mci` | Mutual Clustering Information | Similarity | `MutualClusteringInfo(tree1, tree2, normalise = FALSE)` |
| `cid` | Clustering Information Distance | Distance | `ClusteringInfoDist(tree1, tree2, normalise = FALSE)` |

---

## Installation

### Method 1 - Install directly from GitHub

```bash
pip install git+https://github.com/OrthoFinder/PhyloDist.git
```

### Method 2 - Clone the repository

```bash
git clone https://github.com/OrthoFinder/PhyloDist.git
cd PhyloDist
pip install .
```

### Dependencies

- Python ≥ 3.11
- numpy
- scipy
- ete4

## Basic Usage
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

## Command Line Usage

PhyloDist can be run in several ways depending on how it is installed.

---

#### Method 1 - Installed Command Line Tool
After installation, you can run:
```bash
phylodist tree1.txt tree2.txt \
    --method mci \
    --tree1-unrooted \
    --tree2-unrooted \
    --compare-unrooted \
    --normalise
```
#### Method 2 - Run Using the Standalone Script
After downloading the repo locally, you can run:
```bash
python phylodist.py tree1.txt tree2.txt\
    --method mci \
    --tree1-unrooted \
    --tree2-unrooted \
    --compare-unrooted \
    --normalise
```

### Command Line Options
- Required Arguments

| Argument | Description                       |
| -------- | --------------------------------- |
| `tree1`  | First tree file or Newick string  |
| `tree2`  | Second tree file or Newick string |

- General Options

| Option          | Type   | Default    | Description                        |
| --------------- | ------ | ---------- | ---------------------------------- |
| `--method`      | string | `rf`       | Comparison method                  |
| `--normalise`   | flag   | `False`    | Return normalised score in `[0,1]` |
| `--return-type` | string | `distance` | Return `distance` or `similarity`  |
| `--k`           | float  | `1.0`      | Exponent parameter for JRF         |
| `--parser`      | int    | `1`        | ETE parser format                  |

- Rootedness Options

| Option               | Description               |
| -------------------- | ------------------------- |
| `--tree1-rooted`     | Mark tree1 as rooted      |
| `--tree1-unrooted`   | Mark tree1 as unrooted    |
| `--tree2-rooted`     | Mark tree2 as rooted      |
| `--tree2-unrooted`   | Mark tree2 as unrooted    |
| `--compare-rooted`   | Force rooted comparison   |
| `--compare-unrooted` | Force unrooted comparison |

- Output Options

| Option             | Description                       |
| ------------------ | --------------------------------- |
| `--verbose`        | Print diagnostic information      |
| `--print-matching` | Display matched split/clade pairs |
| `--json`           | Output full result as JSON        |
- Available Methods

| Method | Description                         | Output                 |
| ------ | ----------------------------------- | ---------------------- |
| `rf`   | Robinson–Foulds                     | Distance               |
| `jrf`  | Jaccard Robinson–Foulds             | Distance               |
| `icrf` | Information-corrected RF            | Similarity or Distance |
| `nye`  | Nye similarity                      | Similarity             |
| `msi`  | Matching Split Information          | Similarity             |
| `msid` | Matching Split Information Distance | Distance               |
| `spi`  | Shared Phylogenetic Information     | Similarity             |
| `pid`  | Phylogenetic Information Distance   | Distance               |
| `mci`  | Mutual Clustering Information       | Similarity             |
| `cid`  | Clustering Information Distance     | Distance               |

- Get Help
```bash
phylodist --help
```


## Scalability Notes

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

### Practical Recommendations

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

### Important Note

Large trees are supported, but runtime may increase substantially.  
Users working with large phylogenies are encouraged to:

- test on smaller subsets first
- monitor memory usage
- consider parallel workflows for batch comparisons

## Important Interpretation Notes

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


## Reference
- Martin R Smith, Information theoretic generalized Robinson–Foulds metrics for comparing phylogenetic trees, Bioinformatics, Volume 36, Issue 20, October 2020, Pages 5007–5013, [![DOI](https://img.shields.io/badge/DOI-10.1093%2Fbioinformatics%2Fbtaa614-blue)](https://doi.org/10.1093/bioinformatics/btaa614)