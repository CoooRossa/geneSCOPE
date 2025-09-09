# HPC Installation Problems and Solutions

## Current Status: HPC Installation Issues Resolved 🔧

### Problem Summary
- **Original Issue**: 49 dependency packages failed to install on HPC due to `httpgd` conflicts
- **Error Pattern**: "Error in loadNamespace(x): there is no package called 'httpgd'"
- **Root Cause**: Graphics packages trying to install interactive dependencies incompatible with HPC headless environments

### Solution Strategy
Created **4 specialized HPC installation scripts** to handle different scenarios:

## 🛠️ Available HPC Tools

### 1. `install_hpc_dependencies.R` - Primary Installer
**Purpose**: Install dependencies with HPC-specific optimizations
- Disables problematic graphics options
- Uses fallback installation strategies
- Installs essential packages without `Suggests` dependencies
- Handles OpenMP and compilation issues

**Usage**:
```r
source("https://raw.githubusercontent.com/haenolab/geneSCOPE/main/install_hpc_dependencies.R")
```

### 2. `install_genescope_hpc.R` - Package Installer  
**Purpose**: Install geneSCOPE after dependencies are ready
- Validates all dependencies are available
- Installs geneSCOPE without additional dependencies
- Verifies installation success

**Usage**:
```r
source("https://raw.githubusercontent.com/haenolab/geneSCOPE/main/install_genescope_hpc.R")
```

### 3. `emergency_hpc_install.R` - Problem Solver
**Purpose**: Handle severe dependency conflicts (httpgd issues)
- Completely bypasses graphics dependencies
- Installs only core packages without `Suggests`
- Emergency fallback for problematic environments

**Usage**:
```r
source("https://raw.githubusercontent.com/haenolab/geneSCOPE/main/emergency_hpc_install.R")
```

### 4. `diagnose_hpc.R` - Environment Analyzer
**Purpose**: Identify specific HPC environment issues
- Tests compilation capabilities
- Checks problematic packages
- Analyzes network/repository access
- Provides specific recommendations

**Usage**:
```r
source("https://raw.githubusercontent.com/haenolab/geneSCOPE/main/diagnose_hpc.R")
```

## 📋 Step-by-Step Resolution

### For Your Current HPC Issue:

1. **First, run diagnostics**:
   ```bash
   module load R
   R
   ```
   ```r
   source("https://raw.githubusercontent.com/haenolab/geneSCOPE/main/diagnose_hpc.R")
   ```

2. **If httpgd is the main problem** (likely):
   ```r
   source("https://raw.githubusercontent.com/haenolab/geneSCOPE/main/emergency_hpc_install.R")
   ```

3. **If basic dependency issues**:
   ```r
   source("https://raw.githubusercontent.com/haenolab/geneSCOPE/main/install_hpc_dependencies.R")
   source("https://raw.githubusercontent.com/haenolab/geneSCOPE/main/install_genescope_hpc.R")
   ```

## 🎯 Key Optimizations for HPC

### Graphics Handling
- `options(device = "png")` - Safe graphics device
- `Sys.setenv("R_BROWSER" = "false")` - Disable browser
- `options(bitmapType = "cairo")` - Stable bitmap rendering

### Dependency Management
- Install only `Depends` and `Imports`, skip `Suggests`
- Use `dependencies = FALSE` for package installation
- Prioritize binary packages over source compilation

### Network Optimization
- Extended timeout: `options(timeout = 300)`
- Reliable CRAN mirror: `https://cloud.r-project.org/`
- Single-threaded installation: `options(Ncpus = 1)`

## 🔍 Common HPC Issues and Fixes

| Issue | Symptom | Solution |
|-------|---------|----------|
| **httpgd conflicts** | "no package called 'httpgd'" | Use `emergency_hpc_install.R` |
| **Compilation errors** | C++ build failures | Load gcc module, check Makevars |
| **Permission errors** | Cannot write to library | Set `R_LIBS_USER` path |
| **Network timeouts** | Download failures | Use `options(timeout = 300)` |
| **Graphics errors** | X11/display issues | Disable graphics with safe options |

## ✅ Success Indicators

After successful installation:
```r
library(geneSCOPE)
length(ls("package:geneSCOPE"))  # Should show 53+ functions
packageVersion("geneSCOPE")      # Should show version number
```

## 📞 If All Else Fails

1. **Manual package-by-package installation**:
   ```r
   essential <- c("Rcpp", "data.table", "ggplot2", "igraph")
   for (pkg in essential) {
     install.packages(pkg, dependencies = FALSE)
   }
   ```

2. **Local installation from downloaded packages**:
   - Download `.tar.gz` files manually
   - Use `install.packages("path/to/package.tar.gz", repos = NULL)`

3. **Contact system administrator**:
   - Request R development tools installation
   - Ask for graphics library updates
   - Verify OpenMP availability

## 📈 Expected Outcome

With these tools, the HPC installation should succeed by:
1. Avoiding the httpgd dependency trap
2. Installing only essential dependencies
3. Using HPC-compatible compilation settings
4. Providing clear success/failure feedback

**Next Step**: Try the emergency installer first, as it directly addresses the httpgd issue you're experiencing.
