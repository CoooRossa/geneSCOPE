# geneSCOPE Quick Start Guide

## 🚀 Installation

### Standard Installation
```r
# Install from GitHub (when repository is public)
devtools::install_github("CoooRossa/geneSCOPE")
```

### HPC Installation (Recommended for Clusters)
```r
# For HPC environments with dependency issues
source("https://raw.githubusercontent.com/CoooRossa/geneSCOPE/main/emergency_hpc_install.R")
```

## 📚 Basic Usage

### Load the Package
```r
library(geneSCOPE)
# [geneSCOPE] Gene Spatial Correlation Of Pairwise Expression
# [geneSCOPE] Version: 0.0.0.9000
# [geneSCOPE] Detected 32 CPU cores
```

### Available Functions
```r
# View all available functions
ls("package:geneSCOPE")

# Get help for the package
help(package = "geneSCOPE")

# Example core functions:
# - createCoordObj()        # Create coordinate objects
# - computeSpatialWeights() # Calculate spatial weights
# - computeDensity()        # Gene density analysis
# - plotDensity()          # Visualization functions
# - clusterGenes()         # Gene clustering
# - plotNetworkGenes()     # Network visualization
```

## 🔧 Core Workflow Example

```r
# 1. Create coordinate object
coord_obj <- createCoordObj(x_coords, y_coords, gene_data)

# 2. Compute spatial weights
weights <- computeSpatialWeights(coord_obj)

# 3. Analyze gene density
density_result <- computeDensity(coord_obj, gene_name)

# 4. Visualize results
plotDensity(density_result)

# 5. Cluster genes based on spatial patterns
clusters <- clusterGenes(coord_obj, gene_list)

# 6. Network analysis and visualization
plotNetworkGenes(clusters)
```

## 🎯 Key Features

- **Spatial Gene Expression Analysis**: Analyze spatial patterns in gene expression
- **Network-based Clustering**: Cluster genes based on spatial relationships
- **High-Performance Computing**: Optimized for large datasets with parallel processing
- **Flexible Visualization**: Multiple plotting functions for different analysis types
- **Cross-platform Support**: Works on macOS, Linux, and Windows
- **HPC Compatible**: Special installation tools for cluster environments

## 🆘 Troubleshooting

### Installation Issues
```r
# If standard installation fails, try:
source("https://raw.githubusercontent.com/CoooRossa/geneSCOPE/main/diagnose_hpc.R")
# This will identify specific issues in your environment
```

### Getting Help
```r
# View function documentation
?createCoordObj
?computeDensity

# Get package information
packageVersion("geneSCOPE")
citation("geneSCOPE")
```

## 📈 Performance Tips

### For Large Datasets
```r
# Control thread usage
library(geneSCOPE)
# Package automatically optimizes thread usage

# Monitor memory usage
gc()  # Garbage collection
```

### HPC Optimization
```r
# The package automatically:
# - Detects available CPU cores
# - Manages BLAS threading conflicts
# - Optimizes memory usage for large datasets
```

## 📖 Documentation

- **Full Documentation**: `help(package = "geneSCOPE")`
- **HPC Installation**: See `HPC_INSTALL.md`
- **Troubleshooting**: See `HPC_TROUBLESHOOTING.md`
- **Function Reference**: All functions documented with examples

## ✨ What's New in This Version

- ✅ **51 Functions Available**: Complete spatial analysis toolkit
- ✅ **HPC Optimized**: Special installation tools for cluster environments
- ✅ **Cross-platform**: Works on all major operating systems
- ✅ **Thread Management**: Automatic optimization for parallel processing
- ✅ **Clean Documentation**: Warning-free, professional documentation

---

**geneSCOPE** - Gene Spatial Correlation Of Pairwise Expression  
Version: 0.0.0.9000 | Ready for Production Use 🚀
