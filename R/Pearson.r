# ─────────────────────────────────────────────────────────────
#  geneCorrelation()  v2  –  works with cells list slot
# ─────────────────────────────────────────────────────────────
#' @title Gene-gene correlation across grid cells or single cells
#'
#' @description
#'   Computes a gene × gene correlation matrix on either
#'   **grid‐level** data (\code{level = "grid"}) or
#'   **cell‐level** data (\code{level = "cell"}).  
#'   * Grid mode → result stored in
#'     \code{coordObj@grid[[grid_name]][[cor_layer_name]]}.  
#'   * Cell mode → result stored in
#'     \code{coordObj@cells[[cor_layer_name]]}.
#'
#' @param coordObj      A \code{CoordObj}.
#' @param level         "grid" (default) or "cell".
#' @param grid_name     Grid layer to use when \code{level = "grid"}.
#' @param layer         *Grid mode*: "Xz" or "counts".  
#'                      *Cell mode*: name of existing layer inside
#'                      \code{@cells} (default "counts").
#' @param genes         Vector of genes to keep; NULL = all.
#' @param method        "pearson", "spearman", or "kendall".
#' @param blocksize     Number of genes per block (memory control).
#' @param as_sparse     Convert small-coef entries (<\code{thr}) to 0?
#' @param thr           Threshold used when \code{as_sparse = TRUE}.
#' @param use_parallel  Use \code{parallel::mclapply}?
#' @param ncores        Cores for parallel.
#' @param cor_layer_name Name to save; default paste0(method, "_cor").
#'
#' @return Modified \code{CoordObj} invisibly.
#' @export
geneCorrelation <- function(coordObj,
                            level           = c("grid", "cell"),
                            grid_name       = "grid_lenGrid16",
                            layer           = "counts",   # ◀︎ for cell mode
                            genes           = NULL,
                            method          = "pearson",
                            blocksize       = 2000,
                            as_sparse       = TRUE,
                            thr             = 0.3,
                            use_parallel    = FALSE,
                            ncores          = 15,
                            cor_layer_name  = NULL)
{
  level  <- match.arg(level)
  stopifnot(method %in% c("pearson", "spearman", "kendall"))

  ## ---------- fetch matrix ---------------------------------------------
  if (level == "grid") {

    g <- coordObj@grid[[grid_name]]
    if (is.null(g)) stop("Grid layer '", grid_name, "' not found.")

    mat <- switch(layer,
      "Xz"     = g$Xz,
      "counts" = {
        requireNamespace("Matrix"); requireNamespace("data.table")
        dt <- data.table::as.data.table(g$counts)[count > 0]
        dt[, rn := match(grid_id, g$grid_info$grid_id)]
        i <- dt$rn
        gene_levels <- unique(dt$gene)
        j <- match(dt$gene, gene_levels)
        Matrix::sparseMatrix(
          i = i, j = j, x = dt$count,
          dims = c(nrow(g$grid_info), length(gene_levels)),
          dimnames = list(g$grid_info$grid_id, gene_levels)
        )
      },
      stop("Unsupported layer for grid level: ", layer)
    )

  } else {                                                # ---- cell mode
    stopifnot(is.list(coordObj@cells),
              layer %in% names(coordObj@cells))
    mat <- Matrix::t(coordObj@cells[[layer]])             # rows = cells
  }

  ## ---------- gene filter ----------------------------------------------
  if (!is.null(genes)) {
    keep <- intersect(genes, colnames(mat))
    if (!length(keep)) stop("None of the requested genes were found.")
    if (length(keep) < length(genes))
      warning("Dropped ", length(setdiff(genes, keep)),
              " genes absent from matrix.")
    mat <- mat[, keep, drop = FALSE]
  }

  nr <- nrow(mat); nc <- ncol(mat)
  message("[", level, "]  Observations: ", nr, "  |  Genes: ", nc)

  ## ---------- helper (block correlation) -------------------------------
  corBlock <- function(j1, j2 = NULL) {
    m1 <- as.matrix(mat[, j1, drop = FALSE])      # dense, numeric
    m2 <- if (is.null(j2)) m1 else as.matrix(mat[, j2, drop = FALSE])
    stats::cor(m1, m2, method = method)
  }

  ## ---------- compute ---------------------------------------------------
  if (is.na(blocksize) || nc <= blocksize) {

    C <- corBlock(seq_len(nc))

  } else {

    idx <- split(seq_len(nc), ceiling(seq_len(nc) / blocksize))
    n   <- length(idx)

    C <- Matrix::bandSparse(nc, nc, k = 0,
                            diagonals = list(rep(1, nc)),
                            dimnames  = list(colnames(mat), colnames(mat)))

    run_pair <- function(pair) {
      a <- pair[1]; b <- pair[2]
      cab <- corBlock(idx[[a]], idx[[b]])
      list(a = a, b = b, cab = cab)
    }

    pairs <- do.call(rbind,
                     lapply(seq_along(idx),
                            function(a) cbind(a, b = a:length(idx))))
    pairlist <- split(t(pairs), seq_len(nrow(pairs)))

    res <- if (use_parallel && ncores > 1) {
      parallel::mclapply(pairlist, run_pair,
                         mc.cores = ncores, mc.preschedule = FALSE)
    } else {
      lapply(pairlist, run_pair)
    }

    for (r in res) {
      ja <- idx[[r$a]]; jb <- idx[[r$b]]
      if (r$a == r$b) {
        C[ja, jb] <- r$cab
      } else {
        C[ja, jb] <- r$cab
        C[jb, ja] <- t(r$cab)
      }
    }
  }

  ## ---------- sparsify --------------------------------------------------
  if (as_sparse) {
    C[abs(C) < thr] <- 0
    C <- as(C, "dgCMatrix")
  }

  ## ---------- write back -----------------------------------------------
  if (is.null(cor_layer_name))
    cor_layer_name <- paste0(method, "_cor")

  if (level == "grid") {
    coordObj@grid[[grid_name]][[cor_layer_name]] <- C
  } else {
    coordObj@cells[[cor_layer_name]] <- C
  }

  invisible(coordObj)
}


#' @title Compute Δ = Lee’s L − Pearson r with permutation p-values (+QC)
#'
#' @description
#'   Calculates the difference Δ = L − r either on a grid layer
#'   (**level = "grid"**) or on the full single-cell matrix
#'   (**level = "cell"**), performs permutation testing via the
#'   C++ backend \code{leeL_minusR_perm_ge_cpp()}, adds empirical
#'   p-values, FDR, classification, and QC metrics, and writes the
#'   result back into the \code{CoordObj}.
#'
#' @inheritParams computeDeltaLee
#' @param level           "grid" (default) or "cell".
#' @param grid_name       Grid layer name (grid mode only).
#' @param lee_stats_layer Name of the Lee’s L sub-layer / attribute.
#' @param cor_layer_name  Name of the Pearson correlation sub-layer / attribute.
#' @param perms           Number of permutations.
#' @param fdr_method      FDR adjustment method.
#' @param seed            RNG seed.
#'
#' @return If successful, the modified \code{CoordObj}.
#' @export
computeDeltaLee <- function(coordObj,
                            level            = c("grid", "cell"),
                            grid_name        = "grid_lenGrid50",
                            lee_stats_layer  = "LeeStats_Xz",
                            cor_layer_name   = "pearson_Xz",
                            use_delta = FALSE,
                            delta_threshold = 0.1,
                            r_threshold     = 0.1
                            )
{
  level <- match.arg(level)
  stopifnot(grid_name %in% names(coordObj@grid))
  g <- coordObj@grid[[grid_name]]
  L_mat <- g[[lee_stats_layer]]$L
  ## ---------- A. fetch matrices ----------------------------------------
  if (level == "grid") {
    R_mat <- g[[cor_layer_name]]
    stopifnot(!is.null(R_mat))
  } else {                           # ---- level == "cell" --------------
    R_mat  <- coordObj@cells[[cor_layer_name]]
    stopifnot(!is.null(R_mat))
  }

  ## Type safety
  if (inherits(L_mat, "Matrix")) L_mat <- as.matrix(L_mat)
  if (inherits(R_mat, "Matrix")) R_mat <- as.matrix(R_mat)

  ## ---------- B. intersect gene sets -----------------------------------

  common <- intersect(colnames(L_mat), colnames(R_mat))

  if (!length(common))
    stop("No overlapping genes between L and r matrices.")

  L_mat  <- L_mat[common, common, drop = FALSE]
  R_mat  <- R_mat[common, common, drop = FALSE]
  Delta  <- L_mat - R_mat                                # observed Δ

  ## ---------- E. classification ----------------------------------------
  class_mat <- matrix("Other",
    nrow = nrow(Delta), ncol = ncol(Delta),
    dimnames = dimnames(Delta)
  )
  if(!use_delta){
  class_mat[(L_mat > 0) & (R_mat > 0)] <- "Intra-cellular"
  class_mat[(L_mat > 0) & (R_mat < r_threshold)] <- "Inter-cellular"
  } else {
  class_mat[(L_mat > 0) & (R_mat > 0)] <- "Intra-cellular"
  class_mat[(L_mat > 0) & (Delta > delta_threshold)] <- "Inter-cellular"
  }

  ## ---------- G. assemble & write back ---------------------------------
  out <- list(
    Delta = Delta,
    class = class_mat
  )

  if (level == "grid") {
    coordObj@grid[[grid_name]][["LeeMinusR_bygrid_stats"]] <- out
  } else {  # cell level
    coordObj@grid[[grid_name]][["LeeMinusR_bycell_stats"]] <- out
  }
  coordObj
}