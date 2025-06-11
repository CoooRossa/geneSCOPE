# FG²CLI 
## A Fast Geometric Gene Clustering using Lee’s L and Idelta

**Spatial gene‑co‑distribution analysis toolkit (R/C++17, OpenMP)**  
Supported platforms: **macOS** (Intel & Apple Silicon) **or Linux** (x86‑64, ARM64)

***
## 1 Prerequisites

| Component | macOS | Linux |
|-----------|-------|-------|
| R | ≥ 4.4 (download from <https://cran.r‑project.org>) | ≥ 4.4 (system repo or CRAN binary) |
| OpenMP runtime | `brew install libomp` | `apt install libomp-dev` \| `dnf install libomp` |
| Build tools | Xcode CLT (`xcode‑select --install`) | GCC ≥ 9 (or Clang + libstdc++17) |

### Optional but recommended
* **Homebrew** on macOS – <https://brew.sh>
* **pkg‑config** for easier configure detection


### Configure Makevars (two options)

<details>
<summary><strong>Option A — quick one‑liner</strong></summary>

Append the three needed lines in one shot:

```bash
cat >> ~/.R/Makevars <<'EOF'
CXX_STD      = CXX17
PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS) -I$(brew --prefix libomp)/include
PKG_LIBS     = $(SHLIB_OPENMP_CXXFLAGS) -L$(brew --prefix libomp)/lib -lomp
EOF
```
*(Open a new terminal or run `R CMD config CXX` to confirm the change.)*

</details>

<details>
<summary><strong>Option B — manual edit (fine‑grained)</strong></summary>

1. Open the file:
   ```bash
   nano -w ~/.R/Makevars   # or vi / code / gedit
   ```
2. Add / adjust these lines (avoid duplicates):
   ```makefile
   CXX_STD = CXX17
   PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS) -I$(brew --prefix libomp)/include
   PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS) -L$(brew --prefix libomp)/lib -lomp
   ```
3. **Linux users:** replace `$(brew --prefix libomp)` with `/usr` or your distro‑specific prefix if headers & libs live elsewhere.

</details>

***
## 2 Installation

### Install from GitHub
```r
# a) remotes
install.packages("remotes")
remotes::install_git("https://github.com/CoooRossa/FG2CLI.git", ref = "main")

# b) devtools
install.packages("devtools")     # one‑time
library(devtools)
install_github("CoooRossa/FG2CLI", ref = "main", dependencies = TRUE)
```

> **Tip (macOS & Linux):** Before starting R set
> ```bash
> export OPENBLAS_NUM_THREADS=1   # Keep BLAS single‑threaded
> export OMP_NUM_THREADS=8        # Adjust to <= physical cores
> ```

***
## 3 Troubleshooting

| Error message | Cause | Quick fix |
|---------------|-------|-----------|
| `symbol not found '___kmpc_*'` | OpenMP runtime not linked | Install `libomp` and confirm `Makevars` paths |
| `clang: error: unsupported option '-fopenmp'` | Missing Xcode CLT or GCC OpenMP | macOS: `xcode‑select --install` · Linux: `apt install build-essential libomp-dev` |
| R session aborts when loading FG²CLI | BLAS & OpenMP oversubscribe threads | `OPENBLAS_NUM_THREADS=1` and set reasonable `OMP_NUM_THREADS` |
| `omp.h` not found (Apple Silicon) | Wrong include path | Use `/opt/homebrew/opt/libomp/include` in `Makevars` |
| Compilation fails on Linux | Old GCC/Clang without C++17 | Upgrade compiler (GCC ≥ 9) or install `g++‑11` etc. |

If problems persist, open an issue and attach the full compiler log plus `sessionInfo()`.  
Happy clustering!