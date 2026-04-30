# geneSCOPE NEWS

## v1.1.0 (2026-04-30)

### Major Changes

- **C++-first compute pipeline**: Core statistics (Lee's L, Pearson correlation, spatial weights) now use C++ by default with automatic R fallback. No manual backend selection required.
- **Segmentation-hybrid grid construction**: Xenium grid construction uses `grid_keep_mode = "segmentation_hybrid"` by default, retaining segmentation-supported bins as the primary analysis unit with strongly-supported molecule-only fallback bins.
- **Darwin-safe native dispatch**: macOS/Darwin spatial grid/listw helpers use R fallback by default (safe default); Lee's L and Pearson C++ paths remain enabled.
- **Approximate q-values**: `computeL()` supports low-permutation q-value approximation via `approximate_q = TRUE` for large-panel screening (not a replacement for adequate permutations).

### New Features

- `clusterGenes()` q-distribution guard with automatic `k_perms` fallback when similarity spread is too narrow.
- `clusterGenes(mode = "safe_sequential")` for large Visium/CosMx graphs.
- `computeLvsRCurve()` for empirical Lee's L vs Pearson-r curve fitting.
- `getTopLvsR()` with Delta (Lee's L - Pearson r) exploratory output.
- `computeDensity()` for per-grid gene/module density maps.
- Visium workflow with SCT normalization support (`computeL_visium()`).
- CosMx workflow with flexible grid and clustering options.
- Verbose logging system with START/DONE markers and configurable log levels.
- Style-B variance correction for Lee's L statistics with `weight_style` metadata persistence.
- Backend dispatch contract across `computeL()`, `computeCorrelation()`, and `compareLeesL()`.

### Bug Fixes

- Fixed grid alignment in `normalizeMoleculesInGrid()` ensuring Xz rows match grid_info.
- Fixed FDR calculation to use global Benjamini-Hochberg correction.
- Fixed duplicate function definitions across multiple source files.
- Fixed hex kernel coordinate conversion boundary conditions.
- Fixed Visium `norm_layer` forwarding in `computeL_visium()`.
- Fixed permutation retry loop to terminate cleanly on structural failures.
- Fixed consensus sparse-matrix construction in `clusterGenes()`.
- Fixed `compareLeesL()` support percentage gating with missing values.
- Fixed block permutation to use within-block row shuffles.


