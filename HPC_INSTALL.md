# HPC Installation Guide for geneSCOPE

## Quick Start for HPC Environments

### Method 1: Automated Installation (Recommended)
```bash
module load R/4.3.0  # or your available R version
R
```

```r
# Step 1: Install dependencies (HPC-optimized)
source("https://raw.githubusercontent.com/haenolab/geneSCOPE/main/install_hpc_dependencies.R")

# Step 2: Install geneSCOPE
source("https://raw.githubusercontent.com/haenolab/geneSCOPE/main/install_genescope_hpc.R")
```

### Method 2: Emergency Installation (If Method 1 Fails)
```r
# For severe dependency conflicts (e.g., httpgd issues)
source("https://raw.githubusercontent.com/haenolab/geneSCOPE/main/emergency_hpc_install.R")
```

### Method 3: Diagnostic First (If Unsure)
```r
# Run diagnostics to identify specific issues
source("https://raw.githubusercontent.com/haenolab/geneSCOPE/main/diagnose_hpc.R")
# Review output, then choose appropriate installation method
```

## Troubleshooting Common HPC Issues

### Issue 1: httpgd Package Failures
**Problem**: Many packages fail with "there is no package called 'httpgd'"

**Solution**: The HPC installer automatically disables httpgd dependencies. If you still encounter issues:
```r
# Manually configure graphics
options(device = "png")
options(bitmapType = "cairo")
Sys.setenv("R_BROWSER" = "false")
```

### Issue 2: OpenMP Compilation Errors
**Problem**: C++ compilation fails with OpenMP errors

**Solution**: geneSCOPE's Makevars automatically detects HPC environments:
```bash
# Check if OpenMP is available
echo $OMP_NUM_THREADS
module list | grep -i openmp
```

### Issue 3: Network/Download Issues
**Problem**: Package downloads timeout or fail

**Solution**: Use the offline installation method:
```r
# Set longer timeout
options(timeout = 300)

# Use specific CRAN mirror
options(repos = "https://cloud.r-project.org/")

# Install packages one by one if needed
packages <- c("Rcpp", "RcppArmadillo", "data.table", "ggplot2")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = FALSE)
  }
}
```

### Issue 4: Permission Errors
**Problem**: Cannot write to library directory

**Solution**: Install to personal library:
```r
# Check library paths
.libPaths()

# Install to personal library
lib_path <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(lib_path)) dir.create(lib_path, recursive = TRUE)
install.packages("package_name", lib = lib_path)
```

## Manual Installation (If Automated Script Fails)

### Essential Dependencies Only
```r
# Core packages required for geneSCOPE
essential <- c(
  "Rcpp", "RcppArmadillo", "RcppEigen",
  "data.table", "dplyr", "ggplot2", 
  "igraph", "foreach"
)

# Install one by one
for (pkg in essential) {
  cat("Installing", pkg, "...\n")
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/", 
                   dependencies = FALSE, type = "source")
  }
}
```

### Verify Installation
```r
# Test geneSCOPE installation
library(geneSCOPE)

# Check available functions
ls("package:geneSCOPE")

# Test basic functionality
sessionInfo()
```

## System Requirements

### Minimum Requirements
- R ≥ 4.0.0
- C++11 compiler
- Basic linear algebra libraries (BLAS/LAPACK)

### Recommended for HPC
- R ≥ 4.2.0
- Intel MKL or OpenBLAS
- OpenMP support (optional)
- 8GB+ RAM for large datasets

## Performance Optimization for HPC

### Thread Control
```r
# Control OpenMP threads
Sys.setenv("OMP_NUM_THREADS" = "4")

# Control BLAS threads  
library(RhpcBLASctl)
blas_set_num_threads(4)
```

### Memory Management
```r
# Monitor memory usage
gc()

# Increase memory limit if needed
memory.limit(size = 16000)  # Windows
# ulimit -m unlimited        # Unix
```

## Getting Help

1. **Check Package Status**: `library(geneSCOPE); packageVersion("geneSCOPE")`
2. **View Documentation**: `help(package = "geneSCOPE")`
3. **Report Issues**: [GitHub Issues](https://github.com/haenolab/geneSCOPE/issues)

## Success Indicators

✅ **Installation Successful** if you see:
```r
library(geneSCOPE)
# No error messages
length(ls("package:geneSCOPE")) > 50  # Should show 53+ functions
```

❌ **Installation Failed** if you see:
- "there is no package called 'geneSCOPE'"
- "object not found" errors
- Compilation errors during installation
