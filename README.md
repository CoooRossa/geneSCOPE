# geneSCOPE

> **Gene Spatial Correlation Of Pairwise Expression**

[![R build status](https://github.com/CoooRossa/geneSCOPE/workflows/R-CMD-check/badge.svg)](https://github.com/CoooRossa/geneSCOPE/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Tools to compute and visualize spatial correlation of pairwise gene expression, including Lee's L statistics, permutation-based inference, spatial weights, and network-based visualizations over grid-aggregated spatial transcriptomics data.

## Installation

### Quick Install (Recommended)

```r
# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("CoooRossa/geneSCOPE")
```

### Prerequisites

- **R** ≥ 4.1.0
- **C++17** compiler (GCC ≥ 9 or Clang)
- **OpenMP** support (recommended for performance)

#### Platform-specific setup:

**macOS:**
```bash
# Install Xcode Command Line Tools
xcode-select --install

# Install OpenMP (via Homebrew)
brew install libomp
```

**Linux:**
```bash
# Ubuntu/Debian
sudo apt install build-essential libomp-dev

# CentOS/RHEL
sudo yum install gcc-c++ libgomp-devel
```

### HPC/Cluster Installation

For HPC environments where dependencies may fail, use our specialized installers:

```r
# Step 1: Install dependencies (HPC-optimized)
source("https://raw.githubusercontent.com/CoooRossa/geneSCOPE/main/install_hpc_dependencies.R")

# Step 2: Install geneSCOPE
source("https://raw.githubusercontent.com/CoooRossa/geneSCOPE/main/install_genescope_hpc.R")
```

See [HPC_INSTALL.md](HPC_INSTALL.md) for detailed troubleshooting.

### Manual Installation with Dependencies

If you encounter dependency issues on regular systems:

```r
# Download and run dependency installer
source("https://raw.githubusercontent.com/CoooRossa/geneSCOPE/main/install_dependencies.R")

# Then install the package
devtools::install_github("CoooRossa/geneSCOPE")
```

## Quick Start

```r
library(geneSCOPE)

# Load your spatial transcriptomics data
# Create coordinate object
coord_obj <- createCoordObj(coordinates, gene_expression)

# Compute spatial weights
coord_obj <- computeSpatialWeights(coord_obj, method = "queen")

# Calculate Lee's L statistics
coord_obj <- addLeeStats(coord_obj, genes = c("GENE1", "GENE2"))

# Visualize results
plotLeeLDistribution(coord_obj)
plotNetworkGene(coord_obj, gene = "GENE1")
```

## Key Features

- **Spatial Statistics**: Lee's L, Moran's I, and custom spatial correlation metrics
- **Network Analysis**: Gene co-expression networks with spatial constraints
- **Visualization**: Interactive plots for spatial patterns and gene networks
- **High Performance**: C++17 backend with OpenMP parallelization
- **Cross-platform**: Optimized for macOS, Linux, and HPC environments

## Documentation

- **Function Reference**: Run `help(package = "geneSCOPE")` in R
- **Vignettes**: `browseVignettes("geneSCOPE")`
- **Examples**: See `/examples` directory

## Citation

If you use geneSCOPE in your research, please cite:

```
Your Paper Title (2024)
Authors et al.
Journal Name
DOI: xxx
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## Support

- **Issues**: [GitHub Issues](https://github.com/CoooRossa/geneSCOPE/issues)
- **Discussions**: [GitHub Discussions](https://github.com/CoooRossa/geneSCOPE/discussions)
- **Email**: your.email@domain.com
