#' @title Visium/SCTransform Optional Module
#' @description Adapter layer for Visium spatial transcriptomics data processing.
#'   This module requires Seurat, SeuratObject, and glmGamPoi packages.
#'   If not installed, functions will throw informative errors.
#' @name optional-visium
#' @keywords internal
NULL

#' Check if Visium dependencies are available
#' @return Logical indicating if all Visium dependencies are installed
#' @keywords internal
.visium_deps_available <- function() {
  requireNamespace("Seurat", quietly = TRUE) &&
  requireNamespace("SeuratObject", quietly = TRUE)
}

#' Check if glmGamPoi is available (for SCTransform)
#' @return Logical indicating if glmGamPoi is installed
#' @keywords internal
.glmgampoi_available <- function() {
  requireNamespace("glmGamPoi", quietly = TRUE)
}

#' Require Visium dependencies or stop with informative error
#' @param require_glmgam_poi Whether to also require glmGamPoi (default FALSE)
#' @return Invisible TRUE if available
#' @keywords internal
.require_visium_or_stop <- function(require_glmgam_poi = FALSE) {
  if (!.visium_deps_available()) {
    stop(
      "Visium processing requires Seurat and SeuratObject packages.\n",
      "Install with:\n",
      "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
      "  BiocManager::install(c('Seurat', 'SeuratObject'))"
    )
  }
  if (isTRUE(require_glmgam_poi) && !.glmgampoi_available()) {
    stop(
      "SCTransform with glmGamPoi backend requires the glmGamPoi package.\n",
      "Install with:\n",
      "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
      "  BiocManager::install('glmGamPoi')\n",
      "Alternatively, use sctransform = FALSE in createSCOPE()."
    )
  }
  invisible(TRUE)
}

#' Load Visium data using Seurat (adapter)
#' @description Wrapper around Seurat::Load10X_Spatial with dependency checking
#' @param data_dir Path to Visium data directory
#' @param ... Additional arguments passed to Seurat::Load10X_Spatial
#' @return Seurat object
#' @keywords internal
.load_visium_data <- function(data_dir, ...) {
  .require_visium_or_stop()
  Seurat::Load10X_Spatial(data.dir = data_dir, ...)
}

#' Run SCTransform on Seurat object (adapter)
#' @description Wrapper around Seurat::SCTransform with optional glmGamPoi backend
#' @param seurat_obj Seurat object
#' @param use_glmgam_poi Whether to use glmGamPoi backend (default FALSE)
#' @param ... Additional arguments passed to Seurat::SCTransform
#' @return Seurat object with SCT assay
#' @keywords internal
.run_sctransform <- function(seurat_obj, use_glmgam_poi = FALSE, ...) {
  .require_visium_or_stop(require_glmgam_poi = use_glmgam_poi)
  
  if (use_glmgam_poi) {
    Seurat::SCTransform(seurat_obj, method = "glmGamPoi", ...)
  } else {
    Seurat::SCTransform(seurat_obj, ...)
  }
}

#' Get assay data from Seurat object (adapter)
#' @description Wrapper around SeuratObject::GetAssayData
#' @param seurat_obj Seurat object
#' @param assay Assay name (default "SCT")
#' @param slot Slot name (default "data")
#' @return Matrix of assay data
#' @keywords internal
.get_seurat_assay_data <- function(seurat_obj, assay = "SCT", slot = "data") {
  .require_visium_or_stop()
  SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = slot)
}

#' Set default assay in Seurat object (adapter)
#' @description Wrapper around Seurat::DefaultAssay
#' @param seurat_obj Seurat object
#' @param assay Assay name
#' @return Seurat object with updated default assay
#' @keywords internal
.set_seurat_default_assay <- function(seurat_obj, assay) {
  .require_visium_or_stop()
  Seurat::DefaultAssay(seurat_obj) <- assay
  seurat_obj
}

#' Compute Lee's L for Visium data
#' @description
#' Wrapper around computeL() with Visium-specific defaults. The selected source
#' layer is centered and scaled into an internal analysis layer so Lee
#' L is computed on deviations rather than translation-dependent logCPM values.
#' @param scope_obj A scope_object containing Visium data
#' @param grid_name Name of the grid layer to use
#' @param norm_layer Name of the normalized source layer (default "logCPM").
#'   The wrapper centers and scales every gene before Lee L is computed.
#' @param ncores Number of cores to use
#' @param lee_stats_layer_name Custom name for the Lee stats layer
#' @param perms Number of permutations (default 0)
#' @param backend Backend to use ("auto", "r", "cpp")
#' @param verbose Logical; whether to print progress messages
#' @param cache_inputs Cache the z-scored layer inputs (default `FALSE`). The
#'   internal layer itself is retained so later provenance and permutation
#'   checks can reproduce the observed statistic.
#' @param ... Additional arguments passed to computeL()
#' @return Updated scope_object with Lee's L statistics
#' @export
computeL_visium <- function(
    scope_obj,
    grid_name,
    norm_layer = "logCPM",
    ncores = 1L,
    lee_stats_layer_name = NULL,
    perms = 0L,
    backend = "auto",
    verbose = FALSE,
    ...,
    cache_inputs = FALSE) {
  
  # Check that the grid exists
  if (!grid_name %in% names(scope_obj@grid)) {
    stop(sprintf("Grid '%s' not found in scope_obj", grid_name))
  }
  
  # Check that the norm_layer exists
  grid_data <- scope_obj@grid[[grid_name]]
  if (!norm_layer %in% names(grid_data)) {
    stop(sprintf("Layer '%s' not found in grid '%s'", norm_layer, grid_name))
  }
  
  # Lee's L is defined on deviations.  A raw logCPM matrix is translation
  # dependent in the native kernel, so create an internal z-scored analysis
  # layer while retaining the public source-layer/stats-layer naming contract.
  X <- as.matrix(grid_data[[norm_layer]])
  if (!is.numeric(X) || any(!is.finite(X))) {
    stop("Visium norm_layer must be a finite numeric matrix.", call. = FALSE)
  }
  Xz <- scale(X, center = TRUE, scale = TRUE)
  Xz[!is.finite(Xz)] <- 0
  dimnames(Xz) <- dimnames(X)
  analysis_layer <- paste0(".__geneSCOPE_visium_Xz_", norm_layer)
  scope_obj@grid[[grid_name]][[analysis_layer]] <- Xz

  # Determine the stats layer name
  if (is.null(lee_stats_layer_name)) {
    lee_stats_layer_name <- paste0("LeeStats_", norm_layer)
  }
  
  # Retain the exact centred/scaled layer. Downstream Delta permutations and
  # provenance validation must use the same observations that produced L.
  computeL(
    scope_obj = scope_obj,
    grid_name = grid_name,
    norm_layer = analysis_layer,
    ncores = ncores,
    lee_stats_layer_name = lee_stats_layer_name,
    perms = perms,
    backend = backend,
    cache_inputs = cache_inputs,
    verbose = verbose,
    ...
  )
}
