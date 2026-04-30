#' @title Advanced Visualization Optional Module
#' @description Adapter layer for advanced visualization features.
#'   Requires ComplexHeatmap, circlize, ggfx, seriation, dendextend, dendsort, gclus.
#' @name optional-advanced-viz
#' @keywords internal
NULL

#' Check if ComplexHeatmap is available
#' @return Logical
#' @keywords internal
.complexheatmap_available <- function() {
  requireNamespace("ComplexHeatmap", quietly = TRUE) &&
  requireNamespace("circlize", quietly = TRUE)
}

#' Check if ggfx is available
#' @return Logical
#' @keywords internal
.ggfx_available <- function() {
  requireNamespace("ggfx", quietly = TRUE)
}

#' Check if seriation is available
#' @return Logical
#' @keywords internal
.seriation_available <- function() {
  requireNamespace("seriation", quietly = TRUE)
}

#' Check if dendextend is available
#' @return Logical
#' @keywords internal
.dendextend_available <- function() {
  requireNamespace("dendextend", quietly = TRUE)
}

#' Check if dendsort is available
#' @return Logical
#' @keywords internal
.dendsort_available <- function() {
  requireNamespace("dendsort", quietly = TRUE)
}

#' Check if gclus is available
#' @return Logical
#' @keywords internal
.gclus_available <- function() {
  requireNamespace("gclus", quietly = TRUE)
}

#' Require ComplexHeatmap or stop
#' @return Invisible TRUE
#' @keywords internal
.require_complexheatmap_or_stop <- function() {
  if (!.complexheatmap_available()) {
    stop(
      "Advanced heatmap visualization requires ComplexHeatmap and circlize.\n",
      "Install with:\n",
      "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
      "  BiocManager::install('ComplexHeatmap')\n",
      "  install.packages('circlize')"
    )
  }
  invisible(TRUE)
}

#' Create advanced heatmap (adapter)
#' @description Wrapper around ComplexHeatmap::Heatmap
#' @param matrix Matrix to plot
#' @param ... Additional arguments passed to Heatmap
#' @return Heatmap object
#' @keywords internal
.create_advanced_heatmap <- function(matrix, ...) {
  .require_complexheatmap_or_stop()
  ComplexHeatmap::Heatmap(matrix, ...)
}

#' Create color ramp for heatmaps (adapter)
#' @description Wrapper around circlize::colorRamp2
#' @param breaks Break points
#' @param colors Colors for each break
#' @return Color mapping function
#' @keywords internal
.create_color_ramp <- function(breaks, colors) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("circlize package required for colorRamp2")
  }
  circlize::colorRamp2(breaks, colors)
}

#' Apply ggfx effects (adapter)
#' @description Apply visual effects using ggfx
#' @param plot ggplot object
#' @param effect_type Type of effect ("blend", "reference", etc.)
#' @param ... Additional arguments for the effect
#' @return Modified ggplot object
#' @keywords internal
.apply_ggfx_effect <- function(plot, effect_type = "blend", ...) {
  if (!.ggfx_available()) {
    warning("ggfx not available, returning original plot")
    return(plot)
  }
  
  if (effect_type == "blend") {
    ggfx::with_blend(plot, ...)
  } else if (effect_type == "reference") {
    ggfx::as_reference(plot, ...)
  } else {
    plot
  }
}

#' Seriate matrix (adapter)
#' @description Reorder matrix using seriation
#' @param mat Matrix to seriate
#' @param method Seriation method (default "PCA")
#' @return Seriated matrix
#' @keywords internal
.seriate_matrix <- function(mat, method = "PCA") {
  if (!.seriation_available()) {
    warning("seriation not available, returning original matrix")
    return(mat)
  }
  
  order <- seriation::seriate(mat, method = method)
  mat[seriation::get_order(order), seriation::get_order(order)]
}

#' Rotate dendrogram (adapter)
#' @description Rotate dendrogram using dendextend
#' @param dend Dendrogram object
#' @param ... Additional arguments
#' @return Rotated dendrogram
#' @keywords internal
.rotate_dendrogram <- function(dend, ...) {
  if (!.dendextend_available()) {
    warning("dendextend not available, returning original dendrogram")
    return(dend)
  }
  dendextend::rotate(dend, ...)
}

#' Sort dendrogram (adapter)
#' @description Sort dendrogram using dendsort
#' @param dend Dendrogram object
#' @return Sorted dendrogram
#' @keywords internal
.sort_dendrogram <- function(dend) {
  if (!.dendsort_available()) {
    warning("dendsort not available, returning original dendrogram")
    return(dend)
  }
  dendsort::dendsort(dend)
}

#' Reorder hclust (adapter)
#' @description Reorder hierarchical clustering using gclus
#' @param hclust_obj hclust object
#' @param dist Distance matrix
#' @return Reordered hclust object
#' @keywords internal
.reorder_hclust <- function(hclust_obj, dist) {
  if (!.gclus_available()) {
    warning("gclus not available, returning original hclust")
    return(hclust_obj)
  }
  gclus::reorder.hclust(hclust_obj, dist)
}
