#' @title Compute and Save a Spatial Weights Matrix (Queen Adjacency) for a Specified Grid Layer
#' @description
#'   Construct an undirected Queen adjacency sparse matrix based on
#'   `coordObj@grid[[grid_name]]$grid_info`, 并保存到
#'   `coordObj@grid[[grid_name]]$W`，供后续 Lee’s L 计算使用。
#'
#' @param coordObj  A `coordObj` object containing grid information.
#' @param grid_name Character. Name of the grid layer to read. If NULL, attempts to auto-detect
#'                 a unique sublayer; otherwise must be specified explicitly.
#' #' @return The modified `coordObj`, with a new entry:
#'  - `coordObj@grid[[grid_name]]$W`: the spatial weights matrix (Queen adjacency).
#'
#' @importFrom Matrix Matrix
#' @importFrom spdep cell2nb nb2listw listw2mat
#' @export
computeSpatialWeights <- function(coordObj, grid_name = NULL) {
  # 0. Validate that coordObj@grid exists and determine grid_layer_name
  stopifnot(!is.null(coordObj@grid))
  
  # Determine the grid layer name
  # If grid_name is NULL, try to auto-detect a unique sublayer
  # If multiple sub-layers exist, require explicit specification
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
  stopifnot(grid_layer_name %in% names(coordObj@grid))
  # Ensure the specified grid layer has the required components
    stopifnot("grid_info"  %in% names(coordObj@grid[[grid_layer_name]]))
  stopifnot("xbins_eff"  %in% names(coordObj@grid[[grid_layer_name]]))
  stopifnot("ybins_eff"  %in% names(coordObj@grid[[grid_layer_name]]))
  
  # 1. Extract grid_info, xbins_eff, and ybins_eff
  #    from the specified grid layer
  grid_info  <- coordObj@grid[[grid_layer_name]]$grid_info
  xbins_eff  <- coordObj@grid[[grid_layer_name]]$xbins_eff
  ybins_eff  <- coordObj@grid[[grid_layer_name]]$ybins_eff
  
  # 2. Validate grid_info
  nbins <- max(xbins_eff, ybins_eff)
  if (xbins_eff != ybins_eff) {
    warning(sprintf(
      "Grid xbins (%d) and ybins (%d) are not equal; using max(%d, %d) = %d as nbins.",
      xbins_eff, ybins_eff, nbins
    ))
  }
  
  # Create a matrix to hold grid IDs
  # Initialize a matrix with NA_character_ to hold grid IDs
  # This matrix will have dimensions nbins x nbins
  mat_id <- matrix(NA_character_, nrow = nbins, ncol = nbins)
  
  # Fill the matrix with grid IDs based on grid_info
  # Loop through each row of grid_info to populate mat_id
  # Ensure that gx_i and gy_i are within the bounds of nbins
  for (i in seq_len(nrow(grid_info))) {
    gx_i   <- grid_info$gx[i]
    gy_i   <- grid_info$gy[i]
    # Ensure gx_i and gy_i are within the bounds of nbins
    if (gx_i <= nbins && gy_i <= nbins) {
      mat_id[gx_i, gy_i] <- grid_info$grid_id[i]
    }
  }
  
  # Convert the matrix to a vector of grid IDs
  # This will flatten the matrix row-wise into a single vector
  # This vector will contain all grid IDs in the order they appear in the matrix
  # This is useful for creating the spatial weights matrix later
  cells_all <- as.vector(t(mat_id))
  
  # Create the spatial weights matrix using Queen adjacency
  # Use the spdep package to create a Queen adjacency matrix
  # The cell2nb function creates a neighbor list based on the grid structure
  # The nb2listw function converts the neighbor list to a listw object
  # The listw2mat function converts the listw object to a sparse matrix
  # The resulting matrix W0 will have dimensions nbins x nbins
  nb0 <- spdep::cell2nb(nbins, nbins, type = "queen")
  lw0 <- spdep::nb2listw(nb0, style = "B")
  W0  <- Matrix::Matrix(spdep::listw2mat(lw0), sparse = TRUE)
  
  # Set the row and column names of the spatial weights matrix
  dimnames(W0) <- list(cells_all, cells_all)
  
  # Filter the spatial weights matrix to retain only valid grid cells
  # Extract the grid IDs from grid_info
  valid_grids <- grid_info$grid_id
  # Ensure that valid_grids is a character vector
  W <- W0[valid_grids, valid_grids, drop = FALSE]
  
  # Assign the spatial weights matrix to the coordObj
  # This will add the W matrix to the specified grid layer in coordObj
  coordObj@grid[[grid_layer_name]]$W <- W
  return(coordObj)
}