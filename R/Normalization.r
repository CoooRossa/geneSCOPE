#' @title Construct and Save a Z-score Matrix for Grid × Gene (with optional zero-filling)
#' @description
#'   Based on `coordObj@grid[[grid_name]]$counts`, construct a “grid × gene” matrix,
#'   perform row normalization (so that each grid’s row sums to 1), then apply column-wise
#'   Z-score standardization. Optionally retain grid cells with all-zero counts. The result
#'   is saved under `coordObj@grid[[grid_name]]$Xz`.
#'
#' @param coordObj         A list containing `coordObj@grid[[grid_name]]$counts` and `grid_info`.
#' @param grid_name        Character. Name of the grid layer to read. If NULL, attempts to auto-detect
#'                         a unique sublayer; otherwise must be specified explicitly.
#' @param keep_zero_grids  Logical. Whether to retain grid cells that have zero counts for all genes
#'                         (default TRUE).
#'                         - TRUE: Keep all grid cells, even if their `counts` are all zero.
#'                         - FALSE: Keep only grid cells with at least one gene count > 0.
#'
#' @return The modified `coordObj`, with a new entry:
#'   - `coordObj@grid[[grid_name]]$Xz`: the standardized matrix (rows = grid cells, columns = genes).
#'
#' @import data.table
#' @import Matrix
#' @export
norMPG <- function(coordObj, grid_name = NULL, keep_zero_grids = FALSE) {
  # 0. Validate that coordObj@grid exists and determine grid_layer_name
  stopifnot(!is.null(coordObj@grid))


  if (is.null(grid_name)) {
    sub_layers <- names(coordObj@grid)
    if (length(sub_layers) == 1L) {
      grid_layer_name <- sub_layers[[1]]
    } else {
      stop("Multiple sub-layers found under coordObj@grid; please specify grid_name explicitly.")
    }
  } else {
    grid_layer_name <- as.character(grid_name)
  }
  if (!(grid_layer_name %in% names(coordObj@grid))) {
    stop("Specified grid_name does not exist: ", grid_layer_name)
  }
  grid_layer <- coordObj@grid[[grid_layer_name]]
  stopifnot("counts"    %in% names(grid_layer),
            "grid_info" %in% names(grid_layer))

  # 1. Extract counts_dt and ensure it is a data.table
  counts_dt <- grid_layer$counts
  if (!inherits(counts_dt, "data.table")) {
    counts_dt <- data.table::as.data.table(counts_dt)
  }
  # 2. Determine which grid IDs to retain based on keep_zero_grids
  if (keep_zero_grids) {
    # Retain all grid cells defined in grid_info
    grid_ids_all <- sort(unique(grid_layer$grid_info$grid_id))
  } else {
    # Retain only grid cells that appear in counts_dt (at least one gene count > 0)
    grid_ids_all <- sort(unique(counts_dt$grid_id))
    if (length(grid_ids_all) == 0L) {
      stop("No non-zero grid cells to retain; check keep_zero_grids or counts data.")
    }
  }

  # 3. Extract all genes
  genes_all <- sort(unique(counts_dt$gene))

  # 4. Build a sparse matrix P_mat: rows = genes, columns = retained grid cells
  P_mat <- Matrix::sparseMatrix(
    i        = match(counts_dt$gene,    genes_all),
    j        = match(counts_dt$grid_id, grid_ids_all),
    x        = counts_dt$count,
    dims     = c(length(genes_all), length(grid_ids_all)),
    dimnames = list(genes_all, grid_ids_all)
  )

  # 5. Transpose to obtain “grid × gene” sparse matrix X_full (zero-fill automatically)
  X_full <- Matrix::t(P_mat)
  rownames(X_full) <- grid_ids_all
  colnames(X_full) <- genes_all

  # 6. Row-normalize so each row (grid cell) sums to 1
  #    Convert to dense matrix for normalization
  X_mat <- as.matrix(X_full)
  #    Compute row sums
  row_sums <- rowSums(X_mat)
  #    For rows with sum == 0, set denominator to 1 to avoid division by zero
  zero_rows <- which(row_sums == 0)
  if (length(zero_rows) > 0L) {
    row_sums[zero_rows] <- 1
  }
  P_norm <- X_mat / row_sums
  # 7. Column-wise Z-score standardization
  Xz <- scale(P_norm, center = TRUE, scale = TRUE)
  #    Replace any NA (from zero-variance columns) with 0
  Xz[is.na(Xz)] <- 0
  # 8. Save the standardized matrix under coordObj@grid[[grid_layer_name]]$Xz
  coordObj@grid[[grid_layer_name]]$Xz <- Xz

  return(coordObj)
}
