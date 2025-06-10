#' @title Compute and Store Raw Iδ Values (Supports Specifying Grid Sublayer and Normalization)
#' @description
#'   Read the grid×gene counts (or normalized Z‐score matrix) from a specified sublayer in `coordObj@grid`,  
#'   construct a “gene × grid” sparse or dense matrix, and call the C++ function `idelta_cpp` to compute raw Iδ (no further normalization or log transforms).  
#'   The results are saved under `coordObj@grid[[grid_name]]$iDeltaStats`, and the raw Iδ values are also inserted into `coordObj@meta.data`.
#'
#' @param coordObj   A `coordObj` that already contains outputs of `countMPG`, `normalizeMPG` (if using normalized data), and `computeSpatialWeights`.  
#' @param grid_name  Character (optional). Name of the grid sublayer to use (e.g., `"grid_lenGrid50"`).  
#'                   If NULL and `coordObj@grid` has exactly one sublayer, that sublayer is used. Otherwise, this parameter is required.  
#' @param normalized Logical (default FALSE). Whether to use the normalized matrix:  
#'                   - FALSE: use `coordObj@grid[[grid_name]]$counts` (raw counts).  
#'                   - TRUE: use `coordObj@grid[[grid_name]]$Xz` (Z‐score normalized matrix).  
#'
#' @return The modified `coordObj`, where:  
#'   • `coordObj@grid[[grid_name]]$iDeltaStats` is a list containing:  
#'       - `genes`: Character vector of gene IDs.  
#'       - `delta_raw`: Named numeric vector of each gene’s raw Iδ value.  
#'   • A new column named `"<grid_name>_iDelta_raw"` (or `"<grid_name>_iDelta_nor"` if `normalized = TRUE`) is added to `coordObj@meta.data`, storing each gene’s raw Iδ.  
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom FG2CLI idelta_cpp
#' @export
computeIDeltaMetrics <- function(coordObj,
                                 grid_name  = NULL) {
  # —— 0. Determine grid sublayer name —— 
  if (is.null(coordObj@grid) || length(coordObj@grid) == 0) {
    stop("coordObj@grid is empty; cannot compute Iδ.")
  }
  if (is.null(grid_name)) {
    sub_layers <- names(coordObj@grid)
    if (length(sub_layers) == 1L) {
      grid_layer_name <- sub_layers[[1]]
    } else {
      stop("coordObj@grid contains multiple sublayers; please specify `grid_name` explicitly.")
    }
  } else {
    grid_layer_name <- as.character(grid_name)
  }
  if (!(grid_layer_name %in% names(coordObj@grid))) {
    stop("Specified `grid_name` '", grid_layer_name, "' does not exist in coordObj@grid.")
  }

  # —— 1. Extract either raw counts or the Z‐score matrix —— 
  grid_layer <- coordObj@grid[[grid_layer_name]]
    if (is.null(grid_layer$counts)) {
      stop("coordObj@grid[['", grid_layer_name, "']]$counts is NULL; please run countMPG first.")
    }
    counts_dt <- grid_layer$counts
    # List of all genes and cells
    genes <- sort(unique(counts_dt$gene))
    cells <- sort(unique(counts_dt$grid_id))
    # Build a sparseMatrix: rows = genes, columns = cells
    mat <- Matrix::sparseMatrix(
      i        = match(counts_dt$gene,    genes),
      j        = match(counts_dt$grid_id, cells),
      x        = counts_dt$count,
      dims     = c(length(genes), length(cells)),
      dimnames = list(genes, cells)
    )
  
  # —— 2. Call C++ `idelta_cpp` on the “gene × grid” matrix —— 
  # idelta_cpp expects a matrix with rows = genes, columns = cells
  delta_raw <- FG2CLI::idelta_cpp(mat)
  # Ensure the returned vector is named
  names(delta_raw) <- rownames(mat)

  # —— 3. Save results into coordObj@grid[[grid_layer_name]]$iDeltaStats —— 
  coordObj@grid[[grid_layer_name]]$iDeltaStats <- list(
    genes     = rownames(mat),
    delta_raw = delta_raw
  )

  # —— 4. Write `delta_raw` into coordObj@meta.data as a new column —— 
    # If using raw counts, suffix = "_iDelta_raw"
    meta_colname <- paste0(grid_layer_name, "_iDelta_raw")

  if (is.null(coordObj@meta.data)) {
    coordObj@meta.data <- data.frame(row.names = rownames(mat))
  }
  # Ensure all genes are present as row names in meta.data
  missing_genes <- setdiff(rownames(mat), rownames(coordObj@meta.data))
  if (length(missing_genes) > 0L) {
    na_mat <- matrix(
      NA,
      nrow = length(missing_genes),
      ncol = ncol(coordObj@meta.data),
      dimnames = list(missing_genes, colnames(coordObj@meta.data))
    )
    coordObj@meta.data <- rbind(coordObj@meta.data, as.data.frame(na_mat))
  }

  coordObj@meta.data[rownames(mat), meta_colname] <- delta_raw

  return(coordObj)
}
