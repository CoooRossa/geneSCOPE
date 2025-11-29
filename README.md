# geneSCOPE: gene spatial correlation of pairwise expression

[![R](https://img.shields.io/badge/R-%3E%3D4.4.1-blue.svg)](https://cran.r-project.org/)
[![License: GPL-3](https://img.shields.io/badge/License-GPL%203-yellow.svg)](https://opensource.org/licenses/GPL-3.0)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

**geneSCOPE** is an R package for comprehensive spatial gene expression analysis using Lee's L statistics. It provides a complete workflow for analyzing spatially resolved transcriptomics data, particularly designed for spatial transcriptomics datasets such as Xenium.

## Features

- **Spatial Statistics**: Compute Lee's L statistics for spatial gene co-expression analysis
- **Gene Network Analysis**: Build and visualize gene co-expression networks with spatial context
- **Gene Clustering**: Identify spatially co-expressed gene modules using advanced clustering algorithms
- **Quality Control**: Comprehensive visualization tools for data validation and exploration
- **High Performance**: Optimized C++ implementations for large-scale spatial datasets
- **Flexible Workflow**: Modular design allowing customization for different analysis needs

## Installation

### Prerequisites

You can install the required dependencies using either conda or R's built-in package manager.

#### Option 1: Using Conda (Recommended)

Create a new conda environment with all required packages in one step:

```bash
conda create -n genescope -c conda-forge -c bioconda \
  "r-base>=4.4.1" \
  r-data.table r-dplyr r-foreach r-ggplot2 r-ggraph r-ggrepel \
  r-igraph r-scales r-sf r-tidyr r-devtools r-arrow \
  r-rcpparmadillo r-rcppeigen r-future.apply r-ggforce \
  r-ggnewscale r-spdep r-rhpcblasctl

# Activate the environment
conda activate genescope
```

#### Option 2: Using R Package Manager

Launch R and run the following commands:

```r
# Install BiocManager (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install CRAN packages
install.packages(c(
  "rhpcBLASctl",
  "arrow",
  "data.table",
  "dplyr",
  "foreach",
  "future.apply",
  "ggforce",
  "ggnewscale",
  "ggplot2",
  "ggraph",
  "ggrepel",
  "igraph",
  "scales",
  "sf",
  "spdep",
  "tidyr",
  "RcppArmadillo",
  "RcppEigen",
  "devtools"
))

# Install Bioconductor packages
BiocManager::install(c(
  "rhdf5"
))
```

### Install geneSCOPE

Once dependencies are installed, install geneSCOPE from GitHub:

```r
# Install geneSCOPE
devtools::install_github("CoooRossa/geneSCOPE")
```

## Quick Start

Here's a basic example of how to use geneSCOPE for spatial gene expression analysis:

```r
library(geneSCOPE)

# Load example data path
data_path <- "/path/to/your/xenium/data"

# Step 1: Create SCOPE object
scope_obj <- createSCOPE(
    xenium_dir = data_path,
    grid_length = c(30),
    seg_type = "cell",
    ncores = 8
)

# Step 2: Add single cell data (Optional but recommended)
scope_obj <- addSingleCells(
    scope_obj = scope_obj,
    xenium_dir = data_path
)

# Step 3: Data normalization
scope_obj <- normalizeMoleculesInGrid(
    scope_obj = scope_obj,
    grid_name = "grid30"
)

scope_obj <- normalizeSingleCells(scope_obj)

# Step 4: Compute spatial weights and Lee's L
scope_obj <- computeWeights(
    scope_obj = scope_obj,
    grid_name = "grid30"
)

scope_obj <- computeL(
    scope_obj = scope_obj,
    grid_name = "grid30",
    ncores = 8
)

# Step 5: Gene clustering
scope_obj <- clusterGenes(
    scope_obj = scope_obj,
    grid_name = "grid30",
    pct_min = "q95.0",
    cluster_name = "spatial_clusters"
)

# Step 6: Network visualization
plotNetwork(
    scope_obj = scope_obj,
    grid_name = "grid30",
    cluster_vec = "spatial_clusters",
    pct_min = "q95.0"
)

# Step 7: Compute L vs R and identify top gene pairs
topLvsR <- getTopLvsR(
    scope_obj = scope_obj,
    pear_level = "cell",
    grid_name = "grid30",
    cor_method = "pearson",
    ncores = 8,
    top_n = 500,
    direction = "largest"
)
```

## Workflow Overview

The geneSCOPE analysis pipeline consists of three main parts:

### Part 1: Basic Data Processing
- **Object Construction**: Create SCOPE object from Xenium data
- **Data Integration**: Add single cell data and normalize expression
- **Quality Control**: Optional preliminary visualizations

### Part 2: Spatial Analysis
- **Spatial Statistics**: Compute Lee's L statistics for spatial gene relationships
- **Correlation Analysis**: Calculate gene expression correlations
- **Curve Fitting**: Optional L vs R curve analysis for validation

### Part 3: Gene Clustering and Networks
- **Gene Clustering**: Identify spatially co-expressed gene modules
- **Network Construction**: Build gene co-expression networks
- **Visualization**: Create publication-ready network plots and dendrograms

## Key Functions

| Function | Description |
|----------|-------------|
| `createSCOPE()` | Create SCOPE object from spatial data |
| `computeWeights()` | Calculate spatial weight matrix |
| `computeL()` | Calculate Lee's L spatial statistics |
| `computeCorrelation()` | Calculate gene expression correlations |
| `getTopLvsR()` | Identify top gene pairs with spatial correlation differences |
| `clusterGenes()` | Perform spatial gene clustering |

## Data Format

geneSCOPE works with Xenium spatial transcriptomics data format, including:

- Gene expression count matrices
- Spatial coordinate information
- Cell segmentation data

## Documentation

For detailed usage examples and function references, see:

- [Complete Workflow Tutorial](vignettes/geneSCOPE_basic_workflow.html)
- [Function Reference](man/)
- [GitHub Issues](../../issues) for questions and bug reports

## Citation

If you use geneSCOPE in your research, please cite:

```
```

## License

This project is licensed under the GPL-3 License - see the [LICENSE](LICENSE) file for details.

## Support

- **Issues**: Report bugs or request features via [GitHub Issues](../../issues)
- **Discussions**: Join the community discussion in [GitHub Discussions](../../discussions)
- **Email**: Contact the maintainers at [0322704@ed.tus.ac.jp]

---

**Note**: This package is under active development. Some features may change in future versions. Please check the [releases page](../../releases) for the latest stable version.
