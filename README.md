# geneSCOPE: gene Spatial Co-Occurrence of Pairwise Expression

[![R](https://img.shields.io/badge/R-%3E%3D4.4.1-blue.svg)](https://cran.r-project.org/)
[![License: GPL-3](https://img.shields.io/badge/License-GPL%203-yellow.svg)](https://opensource.org/licenses/GPL-3.0)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

**geneSCOPE** is an R package for comprehensive spatial gene expression
analysis using Lee's L statistics. It provides a complete workflow for
analyzing spatially variable and spatially co-occurring genes in spatial
transcriptomics data, particularly Xenium, CosMx, and Visium.

## Features

- **Spatial Statistics**: Compute Lee's L statistics for spatial gene co-occurrence analysis
- **Gene Network Analysis**: Build and visualize gene co-occurrence networks with spatial context
- **Gene Clustering**: Identify spatially co-occurring gene modules using clustering algorithms
- **Quality Control & Visualization**: Plot-ready QC and exploratory views for grids, single cells, and modules
- **High Performance**: C++ implementations with R fallback for large spatial datasets
- **Flexible Workflow**: Modular design allowing customization for different analysis needs

## Installation

### Option 1 (Recommended): Conda-heavy install

Use conda to install and lock dependencies, including R >= 4.4.1:

```bash
conda create -n genescope -c conda-forge -c bioconda \
  python=3.11 pyarrow \
  "r-base>=4.4.1" r-devtools r-xml2 r-sf r-data.table r-dplyr r-foreach \
  r-ggplot2 r-ggraph r-ggrepel r-igraph r-scales r-tidyr r-arrow \
  r-rcpparmadillo r-rcppeigen r-future.apply r-ggforce r-ggnewscale \
  r-spdep r-rhpcblasctl r-hdf5r bioconductor-rhdf5
conda activate genescope
```

Install geneSCOPE inside the environment:

```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("CoooRossa/geneSCOPE")
```

### Option 2: Pure R install

Ensure your system R is >= 4.4.1, then run:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

install.packages(c(
  "data.table", "dplyr", "foreach", "ggplot2", "ggraph", "ggrepel",
  "igraph", "scales", "sf", "tidyr", "arrow", "RcppArmadillo", "RcppEigen",
  "future.apply", "ggforce", "ggnewscale", "spdep", "rhpcBLASctl",
  "devtools", "hdf5r"
))

# Bioconductor dependencies (optional but recommended for full functionality)
# Note: These dependencies are OPTIONAL and only needed for specific advanced features.
# You can install them as needed based on your use case.
BiocManager::install(c(
  "rhdf5",           # HDF5 file support
  "glmGamPoi",       # Gamma-Poisson GLM for count data
  "STRINGdb",        # STRING database integration
  "dendsort"         # Dendrogram sorting
))

# CRAN optional dependencies (for specific advanced features)
install.packages(c(
  "dendextend",      # Dendrogram manipulation
  "fastcluster",     # Fast hierarchical clustering
  "gclus",           # Clustering utilities
  "ggfx",            # ggplot2 effects
  "seriation",       # Seriation algorithms
  "bigmemory",       # Large matrix support
  "jpeg"             # Image support
))

devtools::install_github("CoooRossa/geneSCOPE")
```

### Optional Dependencies

Some advanced features in geneSCOPE depend on the following optional packages.
Install them only as needed.

**Statistical Computing & Clustering**:
- `bigmemory` (CRAN): Large matrix support for handling ultra-large datasets
- `fastcluster` (CRAN): Fast hierarchical clustering algorithms
- `glmGamPoi` (Bioconductor): Gamma-Poisson GLM models for count data

**Visualization & Graphics**:
- `ggfx` (CRAN): ggplot2 effect layers for advanced visualization
- `jpeg` (CRAN): JPEG image read/write support

**Network & Database**:
- `STRINGdb` (Bioconductor): STRING database interface for gene functional annotation

**Dendrogram & Seriation**:
- `dendextend` (CRAN): Dendrogram manipulation and visualization
- `dendsort` (Bioconductor): Dendrogram sorting algorithms
- `seriation` (CRAN): Seriation algorithms for matrix rearrangement
- `gclus` (CRAN): Clustering utility collection

**File Format Support**:
- `rhdf5` (Bioconductor): HDF5 file format support for large-scale data storage

`pyarrow` is optional but recommended for robust Xenium fallback ingest:

```bash
python3 -m pip install pyarrow
```

## Data-Type Workflows

The examples below follow the same exposed, script-like style as the paper
trial P5 workflow. Replace placeholder paths with your own directories.

### Native C++ Backends on macOS

On Darwin/macOS, geneSCOPE starts with native C++ backends turned off and uses
R fallback paths to avoid platform-specific crashes. To run the C++ paths for
spatial weights, Lee's L, permutations, correlation, and related native
helpers, set the global option below:

```r
options(geneSCOPE.disable_native_all = FALSE)
```

### Xenium

```r
library(geneSCOPE)

xenium_path <- "<XENIUM_DATA_DIR>"
xenium_roi <- "<ROI_COORDINATE_CSV>"  # optional
grid_um <- 30
grid_name <- paste0("grid", grid_um)

# 1) Create SCOPE object from raw data
xenium <- createSCOPE(
    data_dir = xenium_path,
    grid_length = c(grid_um),
    seg_type = "cell",
    coord_file = xenium_roi,
    ncores = 32
)

# 2) Attach single-cell matrix
xenium <- addSingleCells(
    scope_obj = xenium,
    xenium_dir = xenium_path,
    prefer = "xenium"
)

# 3) Normalize single-cell expression
xenium <- normalizeSingleCells(
    scope_obj = xenium,
    input_layer = "counts",
    output_layer = "logCPM",
    scale_factor = 1e4
)

# 4) Normalize grid-level molecule counts
xenium <- normalizeMoleculesInGrid(
    scope_obj = xenium,
    grid_name = grid_name
)

# 5) Build spatial weight matrix
xenium <- computeWeights(
    scope_obj = xenium,
    grid_name = grid_name
)

# 6) Compute Lee's L statistics
xenium <- computeL(
    scope_obj = xenium,
    grid_name = grid_name,
    use_bigmemory = FALSE,
    ncores = 32
)

# 7) Correlate expression at single-cell level
xenium <- computeCorrelation(
    scope_obj = xenium,
    level = "cell",
    layer = "logCPM",
    method = "pearson",
    blocksize = 2000,
    ncores = 32
)

# 8) Fit empirical L vs R curve
curve_name <- paste0("LR_curve_", grid_um)
xenium <- computeLvsRCurve(
    scope_obj = xenium,
    level = "cell",
    grid_name = grid_name,
    ncores = 32,
    downsample = 0.05,
    k_max = 2000,
    n_strata = 1000,
    min_rel_width = 0.15,
    widen_span = 0.1,
    curve_name = curve_name
)

# 9) Cluster genes into spatial modules
cluster_col <- paste0("q95_res0.1_grid", grid_um, "_log1p_freq0.95")
xenium <- clusterGenes(
    scope_obj = xenium,
    grid_name = grid_name,
    L_min = 0,
    algo = "leiden",
    resolution = 0.10,
    pct_min = "q95",
    cluster_name = cluster_col,
    graph_slot_name = cluster_col,
    use_log1p_weight = TRUE,
    use_consensus = TRUE,
    consensus_thr = 0.95,
    n_restart = 1000
)
```

Optional plotting and inspection steps from the P5 workflow:

```r
# L vs R plot
p_lvsr <- plotLvsR(
    scope_obj = xenium,
    grid_name = grid_name,
    pear_level = "cell",
    delta_top_n = 0,
    flip = TRUE
)

# Network and dendrogram
p_network <- plotNetwork(
    scope_obj = xenium,
    lee_stats_layer = "LeeStats_Xz",
    grid_name = grid_name,
    use_consensus_graph = TRUE,
    graph_slot_name = cluster_col,
    cluster_vec = cluster_col,
    show_sign = TRUE,
    drop_isolated = TRUE,
    L_min = 0.11,
    L_min_neg = 0.11
)

p_dendro <- plotDendroNetwork(
    scope_obj = xenium,
    lee_stats_layer = "LeeStats_Xz",
    grid_name = grid_name,
    use_consensus_graph = TRUE,
    graph_slot_name = cluster_col,
    cluster_vec = cluster_col,
    tree_mode = "radial"
)

# Density maps
xenium <- computeDensity(
    scope_obj = xenium,
    grid_name = grid_name,
    layer_name = "counts",
    normalize_method = "none",
    density_name = "CEACAM5",
    genes = "CEACAM5"
)

plotDensity(
    scope_obj = xenium,
    density1_name = "CEACAM5",
    seg_type = "cell",
    grid_name = grid_name
)

# Top L-vs-R pairs
top_pairs <- getTopLvsR(
    scope_obj = xenium,
    grid_name = grid_name,
    pear_level = "cell",
    L_range = c(0.0, 1),
    top_n = 100,
    ncores = 32,
    direction = "largest",
    pval_mode = "uniform",
    curve_layer = curve_name,
    CI_rule = "remove_within"
)
```

### CosMx

Optional pre-processing for faster reads and ROI generation:

```bash
conda create -n cosmx-py311 -c conda-forge python=3.11 pandas numpy pyarrow matplotlib
conda activate cosmx-py311
```

```bash
python path/to/geneSCOPE/cosmx_flatfiles_to_parquet.py \
  --root "<COSMX_ROOT_DIR>" \
  --pixel-size-um 0.120280945 \
  --build-transcripts --build-cells --build-segmentation \
  --overwrite

python path/to/geneSCOPE/auto_rois_from_fovs.py \
  --root "<COSMX_ROOT_DIR>" \
  --method graph --k 0 --min-fovs 3 \
  --eps-um 800 \
  --margin-um 5 \
  --out-dir "<OUTPUT_ROI_DIR>"
```

```r
library(geneSCOPE)

cosmx_path <- "<COSMX_DATA_DIR>"

cosmx <- createSCOPE(
    data_dir = cosmx_path,
    prefer = "cosmx",
    grid_length = c(100),
    seg_type = "cell",
    ncores = 32,
    flip_y = FALSE
)

cosmx <- addSingleCells(
    scope_obj = cosmx,
    data_dir = cosmx_path
)

cosmx <- normalizeSingleCells(
    scope_obj = cosmx,
    input_layer = "counts",
    output_layer = "logCPM",
    scale_factor = 1e4
)

cosmx <- normalizeMoleculesInGrid(
    scope_obj = cosmx,
    grid_name = "grid100"
)

cosmx <- computeWeights(
    scope_obj = cosmx,
    grid_name = "grid100"
)

cosmx <- computeL(
    scope_obj = cosmx,
    grid_name = "grid100",
    use_bigmemory = FALSE,
    ncores = 32
)

cosmx <- clusterGenes(
    scope_obj = cosmx,
    grid_name = "grid100",
    resolution = 0.1,
    pct_min = "q95",
    cluster_name = "q95_r0.1_grid100",
    graph_slot_name = "q95_r0.1_grid100",
    consensus_thr = 0.95,
    n_restart = 100,
    mode = "safe_sequential"
)
```

CosMx relaxed coarse-grid settings such as `pct_min = "q0"` or
`use_significance = FALSE` should be treated as exploratory unless validated
against a smaller ROI or block-wise full-grid run.

### Visium

Prerequisites:

```r
install.packages(c("Seurat", "hdf5r", "sctransform"))
```

```r
library(geneSCOPE)

visium_path <- "<VISIUM_DATA_DIR>"

visium <- createSCOPE(
    data_dir = visium_path,
    prefer = "visium",
    sctransform = TRUE,
    flip_y = FALSE
)

visium <- computeWeights(
    scope_obj = visium,
    grid_name = "grid55"
)

visium <- computeL_visium(
    scope_obj = visium,
    grid_name = "grid55",
    norm_layer = "SCT",
    use_idelta = TRUE,
    S_target = 300,
    mem_limit_GB = 8,
    ncores = 32
)

visium <- clusterGenes(
    scope_obj = visium,
    grid_name = "grid55",
    lee_stats_layer = "LeeStats_SCT",
    resolution = 0.1,
    pct_min = "q95",
    cluster_name = "q95_r0.1_grid55_SCT",
    graph_slot_name = "q95_r0.1_grid55_SCT",
    consensus_thr = 0.95,
    n_restart = 100,
    mode = "safe_sequential"
)

visium <- computeCorrelation(
    scope_obj = visium,
    level = "grid",
    layer = "SCT",
    method = "pearson",
    blocksize = 2000,
    ncores = 32
)
```

For the P5-style Xenium workflow, `grid_length = 30` is the standard starting
point.

## Workflow Overview

The geneSCOPE analysis pipeline consists of three main parts:

### Part 1: Basic Data Processing

- **Object Construction**: Create SCOPE object from Xenium/CosMx/Visium data
- **Data Integration**: Add single-cell data and normalize expression
- **Quality Control**: Optional preliminary visualizations

### Part 2: Spatial Analysis

- **Spatial Statistics**: Compute Lee's L statistics for spatial gene relationships
- **Correlation Analysis**: Calculate gene expression correlations
- **Curve Fitting**: Optional L vs R curve analysis for validation

### Part 3: Gene Clustering and Networks

- **Gene Clustering**: Identify spatially co-expressed gene modules
- **Network Construction**: Build gene co-expression networks
- **Visualization**: Create module footprints and density maps

## Viewing Function Parameters

To view the parameters of any geneSCOPE function, use R's `args()` function:

```r
# View parameters for createSCOPE
args(createSCOPE)

# View parameters for clusterGenes
args(clusterGenes)

# View parameters for computeWeights
args(computeWeights)
```

For detailed parameter descriptions, use `?function_name`:

```r
# View detailed documentation for createSCOPE
?createSCOPE

# View detailed documentation for clusterGenes
?clusterGenes
```

### Important Notes on Parameters

- **Ellipsis (`...`) arguments**: Some functions accept `...` for backward compatibility, but this should not be used as the primary entry point. Always use explicit named parameters for clarity and reliability.
- **P5 workflow recommended parameters**: For the P5 Xenium workflow, use:
  - `clusterGenes(resolution = 0.05-0.1)` for robust clustering results
  - `consensus_thr = 0.95` for consensus graph building
  - `n_restart = 1000` for stable clustering

## Key Functions

| Function | Description |
|----------|-------------|
| `createSCOPE()` | Create SCOPE object from spatial data |
| `addSingleCells()` | Attach single-cell matrices |
| `normalizeMoleculesInGrid()` | Normalize grid-level molecule counts |
| `normalizeSingleCells()` | Normalize single-cell expression |
| `computeWeights()` | Calculate spatial weight matrix |
| `computeL()` | Calculate Lee's L spatial statistics |
| `computeCorrelation()` | Calculate gene expression correlations |
| `computeLvsRCurve()` | Fit empirical Lee's L vs Pearson-r curve |
| `getTopLvsR()` | Identify top gene pairs with spatial-correlation differences |
| `clusterGenes()` | Perform spatial gene clustering |
| `computeDensity()` | Compute per-grid gene/module density |
| `plotNetwork()` | Plot gene co-occurrence networks |
| `plotDendroNetwork()` | Plot dendrogram-style module networks |

## Notes

- Lee's L and Pearson correlation currently use C++ and R implementations; no
  Python compute backend is shipped.
- `computeCorrelation()` currently supports Pearson correlation.
- `computeL(..., approximate_q = TRUE)` is available for low-permutation,
  large-panel screening. It should not replace high-permutation confirmatory
  inference. When active, the Lee's L stats metadata records
  `approximate_q_status`, `approximate_q_pi0_hat`, `approximate_q_warning`,
  and any `approximate_q_fallback` reason.
- `addSingleCells(platform=)` is deprecated. Use
  `addSingleCells(prefer = "xenium")`, `prefer = "cosmx"`, or
  `prefer = "auto"`; conflicting `platform=` and `prefer=` values fail fast.
- Developers who regenerate Rcpp attributes must run
  `Rscript tools/update-rcpp-exports.R` before committing. The script restores
  Darwin spatial guards in both generated R and C++ export wrappers and
  classifies every native export.
- q-guard or `k_perms` fallback clustering should be interpreted cautiously.
- The formal permutation backend is the package-internal C++ path in
  `src/2.LeeL.cpp`; external `sourceCpp()` functions are not used by the
  production pipeline.

## Citation

If you use geneSCOPE in your research, please cite:

```text
geneSCOPE: gene Spatial Co-Occurrence of Pairwise Expression
Shicheng Zhang, Koichi Saeki, Hiroshi Haeno
bioRxiv 2025.12.03.691993; doi: https://doi.org/10.64898/2025.12.03.691993
```

## License

This project is licensed under the GPL-3 License - see the [LICENSE](LICENSE) file for details.

## Support

- **Issues**: Report bugs or request features via [GitHub Issues](../../issues)
- **Email**: Contact the maintainers at [0322704@ed.tus.ac.jp]

---

**Note**: This package is under active development. Some features may change in
future versions. Please check the [releases page](../../releases) for the
latest stable version.
