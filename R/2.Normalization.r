#' 2.Normalization.r (2025-06-27)
#' @title norMPG - Row normalization + Column Z-score (legacy numeric compatibility, memory-efficient compromise)
#' @description
#'   Generates a "grid × gene" matrix → row normalization → **dense conversion**; reproduces legacy results exactly using \code{base::scale()} with sample SD. If the matrix element count exceeds \code{mem_dense_limit} it is still coerced to dense (with a warning).
#' @param coordObj A \code{CoordObj} that already contains a grid layer.
#' @param grid_name Target grid sub-layer; if \code{NULL} the sole sub-layer is used automatically.
#' @param keep_zero_grids Logical. Whether to keep grids whose counts are all zero.
#' @param row_norm Logical. If \code{TRUE} (default) perform row normalization first.
#' @param zero_var_to_zero Logical. If \code{TRUE} (default) set columns with variance ≤ \code{var_eps} to zero.
#' @param var_eps Numeric. Threshold below which variance is considered "near zero."
#' @param mem_dense_limit Numeric. Threshold for `nrow * ncol`; if the product exceeds this value the sparse matrix is coerced to dense anyway (default \code{2e11}, roughly 16 GB of RAM for a \code{double} matrix).
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#' @return The modified \code{CoordObj} containing a new \code{$Xz} slot.
#' @details
#' The function materializes a sparse \code{dgCMatrix} of grid-level transcript
#' counts, optionally normalizes each row to library size, and then performs a
#' column-wise standardization identical to \code{base::scale()} with sample
#' standard deviation (\eqn{n-1}). Conversion to a dense matrix (\code{Xd}) is
#' required to replicate legacy results; \code{mem_dense_limit} protects against
#' inadvertently allocating very large dense matrices.
#'
#' Zero-variance columns are optionally replaced by zeros because
#' \code{scale()} returns \code{NA} for such columns in R < 4.4 and drops
#' the \code{scaled:scale} attribute from R 4.4 onwards. Setting
#' \code{var_eps} to a positive tolerance can aggressively zero
#' near-constant genes, reducing downstream file sizes.
#' @examples
#' \dontrun{
#' ## Compute Z-score expression for 25 µm grid, retaining empty grids
#' P5.coord <- norMPG(
#'   P5.coord,
#'   grid_name       = "25um",
#'   keep_zero_grids = TRUE,
#'   mem_dense_limit = 5e8 # 500 M elements ~ 4 GB RAM
#' )
#' }
#' @importFrom Matrix sparseMatrix rowSums
#' @importFrom data.table as.data.table
#' @export
norMPG <- function(coordObj,
                   grid_name = NULL,
                   keep_zero_grids = FALSE,
                   row_norm = TRUE,
                   zero_var_to_zero = TRUE,
                   var_eps = 0,
                   mem_dense_limit = 2e11,
                   verbose = TRUE) {
  ## ---- 0. Layer check -------------------------------------------------
  if (verbose) message("[geneSCOPE::norMPG] Starting normalization process...")

  stopifnot(!is.null(coordObj@grid))
  grid_layer_name <- if (is.null(grid_name)) {
    nms <- names(coordObj@grid)
    if (length(nms) == 1L) {
      nms
    } else {
      stop("Multiple sub-layers; specify `grid_name`.")
    }
  } else {
    as.character(grid_name)
  }

  if (verbose) message("[geneSCOPE::norMPG] Processing grid layer: ", grid_layer_name)

  g <- coordObj@grid[[grid_layer_name]]
  stopifnot(all(c("counts", "grid_info") %in% names(g)))

  ## ---- 1. Sparse counts → dgCMatrix ----------------------------------
  if (verbose) message("[geneSCOPE::norMPG] Building sparse count matrix...")
  dt <- data.table::as.data.table(g$counts)

  grids <- if (keep_zero_grids) {
    sort(unique(g$grid_info$grid_id))
  } else {
    sort(unique(dt$grid_id))
  }
  genes <- sort(unique(dt$gene))

  if (verbose) message("[geneSCOPE::norMPG] Matrix dimensions: ", length(grids), " grids × ", length(genes), " genes")

  X <- Matrix::sparseMatrix(
    i = match(dt$grid_id, grids),
    j = match(dt$gene, genes),
    x = dt$count,
    dims = c(length(grids), length(genes)),
    dimnames = list(grids, genes)
  )

  ## ---- 2. Row normalization ------------------------------------------
  if (row_norm) {
    if (verbose) message("[geneSCOPE::norMPG] Applying row normalization...")
    lib <- Matrix::rowSums(X)
    lib[lib == 0] <- 1
    X@x <- X@x / lib[X@i + 1L]
  }

  ## ---- 3. Dense center + sample SD -----------------------------------
  elems <- prod(dim(X))
  if (verbose) message("[geneSCOPE::norMPG] Converting to dense matrix for standardization...")
  if (elems > mem_dense_limit) {
    warning(
      "Matrix size ", format(elems, scientific = FALSE),
      " > mem_dense_limit; coercing to dense may use much RAM."
    )
  }

  Xd <- as.matrix(X) # Dense copy
  if (verbose) message("[geneSCOPE::norMPG] Applying column-wise standardization...")
  Xd <- scale(Xd, center = TRUE, scale = TRUE) # base::scale → sample SD

  Xd[is.na(Xd)] <- 0 # Constant columns NA → 0

  ## ---- 4. Zero near-variance columns as requested --------------------
  if (zero_var_to_zero) {
    if (verbose) message("[geneSCOPE::norMPG] Processing zero-variance columns...")
    csd <- attr(Xd, "scaled:scale")
    if (is.null(csd)) { # R 4.4 scale() no longer attaches attributes, can re-compute
      csd <- apply(Xd, 2, sd)
    }
    Xd[, csd <= var_eps] <- 0
  }

  ## ---- 5. Save & return ----------------------------------------------
  g$Xz <- Xd
  coordObj@grid[[grid_layer_name]] <- g

  if (verbose) message("[geneSCOPE::norMPG] Normalization completed, Xz layer stored")

  invisible(coordObj)
}

#' @title Add log-CPM normalized matrix to a CoordObj
#' @description
#'   Reads the specified layer (default \code{"counts"}) from
#'   \code{coordObj@cells}, converts to counts-per-\code{scale_factor}
#'   (default 10,000), applies \code{log1p}, and stores the result as
#'   a new layer (default \code{"logCPM"}). The modified object is
#'   returned invisibly so you can overwrite the original variable.
#' @param coordObj A \code{CoordObj}.
#' @param input_layer Character. Name of the existing layer to read
#'                    (must be a \code{dgCMatrix} inside
#'                    \code{coordObj@cells}); default \code{"counts"}.
#' @param output_layer Character. Name of the layer to write to
#'                     (will overwrite if exists); default \code{"logCPM"}.
#' @param scale_factor Numeric. Library size each cell is scaled to
#'                     before log1p; default 1e4.
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#' @return The modified \code{CoordObj} (invisibly).
#' @details
#' Counts per million (CPM) are computed by dividing each cell's counts by its
#' total library size, multiplying by \code{scale_factor}, and finally applying
#' \code{log1p}. Because the matrix is stored in compressed sparse column (CSC)
#' form, these operations touch only the non-zero entries, preserving sparsity.
#'
#' The resulting layer is functionally equivalent to Seurat's
#' \code{NormalizeData(method = "LogNormalize")}, but this implementation keeps
#' all steps sparse and avoids introducing Seurat as a dependency. The new
#' matrix is stored (or overwritten) in \code{coordObj@cells[[output_layer]]};
#' no other slots are modified.
#' @examples
#' \dontrun{
#' ## Add log-CPM layer and inspect a housekeeping gene
#' P5.coord <- normalizeCellsCPMlog(P5.coord)
#' }
#' @import Matrix
#' @export
normalizeCellsCPMlog <- function(coordObj,
                                 input_layer = "counts",
                                 output_layer = "logCPM",
                                 scale_factor = 1e4,
                                 verbose = TRUE) {
  if (verbose) message("[geneSCOPE::normalizeCellsCPMlog] Starting log-CPM normalization...")

  stopifnot(
    inherits(coordObj, "CoordObj"),
    is.list(coordObj@cells),
    input_layer %in% names(coordObj@cells)
  )

  mat <- coordObj@cells[[input_layer]]
  if (!inherits(mat, "dgCMatrix")) {
    stop("Layer '", input_layer, "' is not a dgCMatrix.")
  }

  if (verbose) message("[geneSCOPE::normalizeCellsCPMlog] Computing library sizes...")
  lib <- Matrix::colSums(mat) # Library size per column
  zerocell <- lib == 0
  lib[zerocell] <- 1 # Avoid division by zero

  if (verbose) message("[geneSCOPE::normalizeCellsCPMlog] Applying CPM transformation and log1p...")
  ## Vectorized: match corresponding column factors for each non-zero element
  col_rep <- rep(seq_len(ncol(mat)), diff(mat@p))
  sf_vec <- (scale_factor / lib)[col_rep]
  out <- mat # Copy sparse structure
  out@x <- log1p(out@x * sf_vec)

  if (any(zerocell)) {
    if (verbose) warning(sum(zerocell), " cells have total count = 0; they remain 0 in log-CPM.")
  }

  coordObj@cells[[output_layer]] <- out

  if (verbose) message("[geneSCOPE::normalizeCellsCPMlog] Log-CPM normalization completed, stored in @cells$", output_layer)

  invisible(coordObj)
}
