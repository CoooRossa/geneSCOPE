## ---- S4 class & methods ----
if (!methods::isClass("CoordObj")) {
  methods::setClass(
    "CoordObj",
    slots = c(
      coord     = "list", # centroids, segmentation, …
      grid      = "list", # grid layers
      meta.data = "data.frame", # miscellaneous gene / grid meta
      cells     = "list",
      stats     = "list", # statistics layers
      density   = "list" # density layers
    )
  )
}

if (is.null(methods::selectMethod("initialize", "CoordObj", optional = TRUE))) {
  methods::setMethod(
    "initialize", "CoordObj",
    function(.Object,
             coord = list(),
             grid = list(),
             meta.data = data.frame(),
             cells = list(),
             stats = list(),
             density = list()) {
      ## slot assignment ----
      .Object@coord <- coord
      .Object@grid <- grid
      .Object@meta.data <- meta.data
      .Object@cells <- cells
      .Object@stats <- stats
      .Object@density <- density

      ## basic validation – avoid downstream function errors ----
      if (nrow(.Object@meta.data) > 0 && is.null(rownames(.Object@meta.data))) {
        warning(
          "meta.data has no row-names; downstream functions assume gene names. ",
          "Consider setting row-names now."
        )
      }
      .Object
    }
  )

  ## concise object printout (PDF suggestion) ----
  methods::setMethod(
    "show", "CoordObj",
    function(object) {
      cat("CoordObj with:\n")
      cat("  coord:", length(object@coord), "elements\n")
      cat("  grid:", length(object@grid), "layers\n")
      cat("  cells:", length(object@cells), "matrices\n")
      cat("  stats:", length(object@stats), "result sets\n")
      cat("  density:", length(object@density), "density tables\n")
      cat("  meta.data:", nrow(object@meta.data), "genes x", ncol(object@meta.data), "features\n")
    }
  )

  ## summary() – quick overview ----
  methods::setMethod(
    "summary", signature(object = "CoordObj"),
    function(object, ...) {
      cat("=== CoordObj Summary ===\n\n")

      # Coordinate data
      if (length(object@coord) > 0) {
        cat("Coordinate data:\n")
        for (nm in names(object@coord)) {
          if (is.data.frame(object@coord[[nm]])) {
            cat(" ", nm, ":", nrow(object@coord[[nm]]), "rows\n")
          } else {
            cat(" ", nm, ":", length(object@coord[[nm]]), "elements\n")
          }
        }
        cat("\n")
      }

      # Grid layers
      if (length(object@grid) > 0) {
        cat("Grid layers:\n")
        for (nm in names(object@grid)) {
          g <- object@grid[[nm]]
          if ("grid_info" %in% names(g)) {
            cat(" ", nm, ":", nrow(g$grid_info), "grid cells\n")
          }
        }
        cat("\n")
      }

      # Cell data
      if (length(object@cells) > 0) {
        cat("Cell matrices:\n")
        for (nm in names(object@cells)) {
          mat <- object@cells[[nm]]
          if (is.matrix(mat) || inherits(mat, "Matrix")) {
            cat(" ", nm, ":", nrow(mat), "x", ncol(mat), "\n")
          }
        }
        cat("\n")
      }

      # Stats
      if (length(object@stats) > 0) {
        cat("Analysis results:\n")
        for (nm in names(object@stats)) {
          cat(" ", nm, ":", length(object@stats[[nm]]), "result sets\n")
        }
        cat("\n")
      }

      # Meta data
      if (nrow(object@meta.data) > 0) {
        cat("Gene metadata:", nrow(object@meta.data), "genes\n")
        if (ncol(object@meta.data) > 0) {
          cat("  Columns:", paste(colnames(object@meta.data), collapse = ", "), "\n")
        }
      }
    }
  )
}

## flipCoordinates ----
#' @title Flip Coordinates
#' @description
#' Flips the y-coordinates of a set of points or shapes, typically to adjust for different coordinate system origins.
#'
#' @param data A data.frame or data.table containing at least a column \code{y} with y-coordinates to flip.
#' @param y_max Numeric value of the maximum Y coordinate (e.g., image height) used as the reference for flipping.
#'
#' @return The same data object with the \code{y} coordinates flipped (i.e., each \code{y} replaced by \code{y_max - y}).
#'
#' @details
#' Author: FG Team. Created: 2025-07-05.
flipCoordinates <- function(data, y_max) {
  stopifnot(is.numeric(y_max) && length(y_max) == 1)
  if (!"y" %in% names(data)) {
    stop("Input data must have a 'y' column for coordinates.")
  }
  # If data is a data.table, modify in place for efficiency
  if (data.table::is.data.table(data)) {
    data[, y := y_max - y]
  } else {
    data$y <- y_max - data$y
  }
  return(data)
}

#' @noRd
.flipCoordinates <- function(data, y_max) {
  stopifnot(is.numeric(y_max) && length(y_max) == 1)
  if (!"y" %in% names(data)) {
    stop("Input data must have a 'y' column for coordinates.")
  }
  # If data is a data.table, modify in place for efficiency
  if (data.table::is.data.table(data)) {
    data[, y := y_max - y]
  } else {
    data$y <- y_max - data$y
  }
  return(data)
}

## clipPointsToPolygon ----
#' @title Clip Points to Polygon
#' @description
#' Filters a set of points, returning only those that lie inside a given polygon region. Operates in chunks to handle large datasets efficiently.
#'
#' @param points A data.frame or data.table of points with columns \code{x} and \code{y}.
#' @param polygon An \code{sfc} polygon or \code{sf} polygon object defining the region of interest.
#' @param chunk_size Integer. Number of points to process per chunk (default 5e5).
#' @param ncores Integer. Number of cores for parallel processing (default 1).
#'
#' @return A data.table (or data.frame) containing only the points from the input that fall within the polygon.
#'
#' @details
#' This function splits the points into chunks to avoid excessive memory usage during spatial queries. It leverages \pkg{sf} for point-in-polygon checks and can use parallel processing on non-Windows platforms for speed.
#'
#' Author: FG Team. Created: 2025-07-05.
## ---------------------------------------------------------------------
##  Clip points to polygon  – bug-fix: ensure 'polygon' is sfc / sf
## ---------------------------------------------------------------------
clipPointsToPolygon <- function(points, polygon,
                                chunk_size = 5e5, ncores = 1) {
  if (is.null(polygon) || nrow(points) == 0) {
    return(points)
  }

  ## -------- 1. Force sfc (as before) ---------------------------------------
  if (inherits(polygon, "sf")) {
    polygon <- sf::st_geometry(polygon)
  }

  if (!inherits(polygon, "sfc")) {
    if (all(vapply(polygon, inherits, logical(1), "sfg"))) {
      polygon <- sf::st_sfc(polygon)
    } else {
      stop("`polygon` must be an sf/sfc object (or list of sfg).")
    }
  }

  ## ❶—— If length 0, return empty result to avoid crash ------------------------------
  if (length(polygon) == 0) {
    return(points[0])
  } # empty data.table / data.frame

  ## (Optional) Only union if more than one feature; prevents EMPTY
  if (length(polygon) > 1) {
    polygon <- sf::st_union(polygon)
  }

  ## ❷—— Check again after union
  if (length(polygon) == 0) {
    return(points[0])
  }

  ## -------- 2. Preprocessing & parallel clipping (as before) -----------------------------
  stopifnot(all(c("x", "y") %in% names(points)))

  dt_pts <- data.table::as.data.table(points)
  crs_use <- sf::st_crs(polygon) %||% NA

  idx_split <- split(
    seq_len(nrow(dt_pts)),
    ceiling(seq_len(nrow(dt_pts)) / chunk_size)
  )

  worker <- function(idx) {
    sub <- dt_pts[idx]
    inside <- lengths(sf::st_within(
      sf::st_as_sf(sub,
        coords = c("x", "y"),
        crs = crs_use, remove = FALSE
      ),
      polygon
    )) > 0
    sub[inside]
  }

  res_list <- if (ncores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(idx_split, worker, mc.cores = ncores)
  } else if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(sf))
    parallel::clusterExport(cl, c("dt_pts", "polygon", "crs_use", "worker"),
      envir = environment()
    )
    parallel::parLapply(cl, idx_split, worker)
  } else {
    lapply(idx_split, worker)
  }

  data.table::rbindlist(res_list)
}

#' @noRd
.clipPointsToPolygon <- function(points, polygon,
                                 chunk_size = 5e5, ncores = 1) {
  if (is.null(polygon) || nrow(points) == 0) {
    return(points)
  }

  ## -------- 1. Force sfc (as before) ---------------------------------------
  if (inherits(polygon, "sf")) {
    polygon <- sf::st_geometry(polygon)
  }

  if (!inherits(polygon, "sfc")) {
    if (all(vapply(polygon, inherits, logical(1), "sfg"))) {
      polygon <- sf::st_sfc(polygon)
    } else {
      stop("`polygon` must be an sf/sfc object (or list of sfg).")
    }
  }

  ## ❶—— If length 0, return empty result to avoid crash ------------------------------
  if (length(polygon) == 0) {
    return(points[0])
  } # empty data.table / data.frame

  ## (Optional) Only union if more than one feature; prevents EMPTY
  if (length(polygon) > 1) {
    polygon <- sf::st_union(polygon)
  }

  ## ❷—— Check again after union
  if (length(polygon) == 0) {
    return(points[0])
  }

  ## -------- 2. Preprocessing & parallel clipping (as before) -----------------------------
  stopifnot(all(c("x", "y") %in% names(points)))

  dt_pts <- data.table::as.data.table(points)
  crs_use <- sf::st_crs(polygon) %||% NA

  idx_split <- split(
    seq_len(nrow(dt_pts)),
    ceiling(seq_len(nrow(dt_pts)) / chunk_size)
  )

  worker <- function(idx) {
    sub <- dt_pts[idx]
    inside <- lengths(sf::st_within(
      sf::st_as_sf(sub,
        coords = c("x", "y"),
        crs = crs_use, remove = FALSE
      ),
      polygon
    )) > 0
    sub[inside]
  }

  res_list <- if (ncores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(idx_split, worker, mc.cores = ncores)
  } else if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(sf))
    parallel::clusterExport(cl, c("dt_pts", "polygon", "crs_use", "worker"),
      envir = environment()
    )
    parallel::parLapply(cl, idx_split, worker)
  } else {
    lapply(idx_split, worker)
  }

  data.table::rbindlist(res_list)
}

## processSegmentation ----
#' @title Process Segmentation Boundaries
#' @description
#' Reads a segmentation boundaries file (e.g., cell or nucleus boundaries) and processes it into coordinate data and polygon geometries.
#'
#' @param seg_file Path to the segmentation boundaries file (Parquet format).
#' @param tag Character tag to label the segmentation type (e.g., "cell" or "nucleus").
#' @param flip_y Logical. Whether to flip the y-coordinates (default \code{FALSE}). If \code{TRUE}, argument \code{y_max} must be provided.
#' @param y_max Numeric. The maximum Y value for flipping coordinates (required if \code{flip_y = TRUE}).
#' @param keep_cells Optional vector of cell IDs. If provided, only segmentation data for these cells will be retained.
#' @param ncores Integer. Number of cores for parallel polygon processing (default 1).
#'
#' @return A list with two components: \code{points} (a data.table of segmentation points with columns \code{cell}, \code{x}, \code{y}, \code{label_id}) and \code{polygons} (an \code{sfc} object of polygon geometries).
#'
#' @details
#' This function requires the \pkg{arrow} package to read Parquet files. After reading the segmentation data, it flips the coordinates if needed, filters by cell IDs if provided, and then constructs polygon geometries for each segmented object.
#'
#' Author: FG Team. Created: 2025-07-05.
processSegmentation <- function(seg_file, tag = "cell",
                                flip_y = FALSE, y_max = NULL,
                                keep_cells = NULL, ncores = 1) {
  if (!file.exists(seg_file)) {
    stop("Segmentation file not found: ", seg_file)
  }

  seg_dt <- data.table::as.data.table(
    arrow::read_parquet(seg_file)
  )[, .(cell = cell_id, x = vertex_x, y = vertex_y, label_id)]

  if (flip_y) {
    if (is.null(y_max)) {
      stop("y_max must be provided when flip_y = TRUE.")
    }
    seg_dt[, y := y_max - y]
  }

  if (!is.null(keep_cells)) {
    seg_dt <- seg_dt[cell %in% keep_cells]
  }

  ## Build polygons (ensure closure) – PDF suggestion
  split_idx <- split(seq_len(nrow(seg_dt)), seg_dt$label_id)
  build_poly <- function(idx) {
    sub <- seg_dt[idx, .(x, y)]
    sub <- sub[complete.cases(sub)] # ① Remove NA
    if (nrow(sub) < 3) {
      return(NULL)
    }

    # ② Ensure double storage
    coords <- as.matrix(sub)
    storage.mode(coords) <- "double"

    # Close polygon if not already closed
    if (!all(coords[1, ] == coords[nrow(coords), ])) {
      coords <- rbind(coords, coords[1, ])
    }

    # ③ try-catch, return NULL on failure
    tryCatch(sf::st_polygon(list(coords)), error = function(e) {
      message("label ", sub$label_id[1], " invalid: ", conditionMessage(e))
      NULL
    })
  }

  poly_list <- if (ncores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(split_idx, build_poly, mc.cores = ncores)
  } else if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(sf))
    parallel::clusterExport(cl, c("seg_dt", "build_poly"),
      envir = environment()
    )
    parallel::parLapply(cl, split_idx, build_poly)
  } else {
    lapply(split_idx, build_poly)
  }

  polygons <- sf::st_sfc(poly_list[!sapply(poly_list, is.null)])
  list(points = seg_dt, polygons = polygons)
}

#' @noRd
.processSegmentation <- function(seg_file, tag = "cell",
                                 flip_y = FALSE, y_max = NULL,
                                 keep_cells = NULL, ncores = 1) {
  if (!file.exists(seg_file)) {
    stop("Segmentation file not found: ", seg_file)
  }

  seg_dt <- data.table::as.data.table(
    arrow::read_parquet(seg_file)
  )[, .(cell = cell_id, x = vertex_x, y = vertex_y, label_id)]

  if (flip_y) {
    if (is.null(y_max)) {
      stop("y_max must be provided when flip_y = TRUE.")
    }
    seg_dt[, y := y_max - y]
  }

  if (!is.null(keep_cells)) {
    seg_dt <- seg_dt[cell %in% keep_cells]
  }

  ## Build polygons (ensure closure) – PDF suggestion
  split_idx <- split(seq_len(nrow(seg_dt)), seg_dt$label_id)
  build_poly <- function(idx) {
    sub <- seg_dt[idx, .(x, y)]
    sub <- sub[complete.cases(sub)] # ① Remove NA
    if (nrow(sub) < 3) {
      return(NULL)
    }

    # ② Ensure double storage
    coords <- as.matrix(sub)
    storage.mode(coords) <- "double"

    # Close polygon if not already closed
    if (!all(coords[1, ] == coords[nrow(coords), ])) {
      coords <- rbind(coords, coords[1, ])
    }

    # ③ try-catch, return NULL on failure
    tryCatch(sf::st_polygon(list(coords)), error = function(e) {
      message("label ", sub$label_id[1], " invalid: ", conditionMessage(e))
      NULL
    })
  }

  poly_list <- if (ncores > 1 && .Platform$OS.type != "windows") {
    parallel::mclapply(split_idx, build_poly, mc.cores = ncores)
  } else if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(sf))
    parallel::clusterExport(cl, c("seg_dt", "build_poly"),
      envir = environment()
    )
    parallel::parLapply(cl, split_idx, build_poly)
  } else {
    lapply(split_idx, build_poly)
  }

  polygons <- sf::st_sfc(poly_list[!sapply(poly_list, is.null)])
  list(points = seg_dt, polygons = polygons)
}


## checkGridContent ----
#' @noRd
.checkGridContent <- function(coordObj, grid_name) {
  # Ensure the specified grid layer exists
  grid_layer <- .selectGridLayer(coordObj, grid_name)
  # Check for grid_info and counts presence
  if (is.null(grid_layer$grid_info)) {
    stop("Grid layer '", grid_name, "' does not contain grid_info.")
  }
  if (is.null(grid_layer$counts)) {
    stop("Grid layer '", grid_name, "' does not contain counts data.")
  }
  # Required columns in grid_info
  required_info_cols <- c("grid_id", "xmin", "xmax", "ymin", "ymax")
  if (!all(required_info_cols %in% colnames(grid_layer$grid_info))) {
    stop("grid_info is missing required columns: ", paste(setdiff(required_info_cols, colnames(grid_layer$grid_info)), collapse = ", "))
  }
  # Required columns in counts data
  required_count_cols <- c("grid_id", "gene", "count")
  if (!all(required_count_cols %in% colnames(grid_layer$counts))) {
    stop("counts data is missing required columns: ", paste(setdiff(required_count_cols, colnames(grid_layer$counts)), collapse = ", "))
  }
  return(invisible(TRUE))
}

## getGeneSubset ----
#' @title Get Gene Subset
#' @description
#' Determines a subset of gene names based on either an explicit list or a cluster membership from the \code{CoordObj}'s metadata.
#'
#' @param coordObj A \code{CoordObj} that contains a \code{meta.data} slot with gene metadata (e.g., clustering results).
#' @param genes Character vector of gene names to select (optional).
#' @param cluster_col Character. Column name in \code{coordObj@meta.data} that contains cluster identifiers (optional).
#' @param cluster_num Value indicating which cluster to select (must correspond to an entry in \code{cluster_col}, optional).
#'
#' @return A character vector of gene names constituting the subset.
#'
#' @details
#' You must specify either \code{genes} directly, or provide both \code{cluster_col} and \code{cluster_num} to define the subset via clustering. If a cluster is specified, the gene names are taken as the row names of \code{coordObj@meta.data} where the given cluster column equals \code{cluster_num}.
#'
# filepath: /Users/haenolab/Documents/FG2CIL Paper/Code/Package/FG2CLI/R/zzz_classes.r

#' @noRd
.getGeneSubset <- function(coordObj, genes = NULL, cluster_col = NULL, cluster_num = NULL) {
  if (!is.null(genes)) {
    # Use the provided gene list directly
    sel_genes <- unique(as.character(genes))
  } else if (!is.null(cluster_col) && !is.null(cluster_num)) {
    # Ensure meta.data exists and contains the cluster column
    if (is.null(coordObj@meta.data) || !(cluster_col %in% colnames(coordObj@meta.data))) {
      stop("coordObj@meta.data not found or does not contain column: ", cluster_col)
    }
    # Get genes belonging to the specified cluster
    sel_genes <- rownames(coordObj@meta.data)[coordObj@meta.data[[cluster_col]] == cluster_num]
    if (length(sel_genes) == 0) {
      stop("No genes found in meta.data where ", cluster_col, " == ", cluster_num, ".")
    }
  } else {
    stop("Please specify either `genes` or both `cluster_col` and `cluster_num`.")
  }
  return(sel_genes)
}

## getLeeMatrix ----
#' @title Get Lee's L Matrix
#' @description
#' Retrieves the Lee’s L spatial correlation matrix from a specified grid layer of a \code{CoordObj}. If a Pearson correlation matrix is present for the same layer, the result will be aligned to have the same gene set.
#'
#' @param coordObj A \code{CoordObj} containing grid layers with computed Lee’s L statistics.
#' @param grid_name Character. Name of the grid layer containing Lee’s L results. If \code{NULL}, a unique grid layer will be used by default (if only one exists).
#' @param lee_layer Character. Name of the sub-layer that stores Lee’s L statistics (default \code{"LeeStats_Xz"}).
#'
#' @return A numeric matrix of Lee’s L values, with row names and column names as gene identifiers.
#'
#' @details
#' The returned matrix is typically symmetric and has genes in both rows and columns. If a Pearson correlation matrix is available in the same grid layer (usually stored as \code{pearson_cor}), the function will restrict the Lee’s L matrix to the intersecting gene set for direct comparability.
#'
#' Author: FG Team. Created: 2025-07-05.
## ======================================================================
##  Helper – getLeeMatrix  (auto-match LeeStats layer; compatible with @stats and @grid)
## ======================================================================
getLeeMatrix <- function(coordObj,
                         grid_name = NULL,
                         lee_layer = NULL) {
  ## ---- 0. Select grid sub-layer ------------------------------------------------
  g_layer <- .selectGridLayer(coordObj, grid_name)
  if (is.null(grid_name)) { # Write back the actual name
    grid_name <- names(coordObj@grid)[
      vapply(coordObj@grid, identical, logical(1), g_layer)
    ]
  }

  ## ---- 1. Auto-detect LeeStats layer name ----------------------------------------
  if (is.null(lee_layer)) {
    cand <- character(0)
    if (!is.null(coordObj@stats[[grid_name]])) {
      cand <- names(coordObj@stats[[grid_name]])
    }
    cand <- c(cand, names(g_layer))
    cand <- unique(cand[grepl("^LeeStats_", cand)])

    if (length(cand) == 0L) {
      stop(
        "No layer starting with 'LeeStats_' found for grid '",
        grid_name, "'."
      )
    }
    if (length(cand) == 1L) {
      lee_layer <- cand
    } else if ("LeeStats_Xz" %in% cand) { # Prefer default naming
      lee_layer <- "LeeStats_Xz"
    } else {
      stop(
        "Multiple LeeStats layers detected (",
        paste(cand, collapse = ", "),
        "); please specify `lee_layer` explicitly."
      )
    }
  }

  ## ---- 2. Search @stats → @grid in order ------------------------------------
  leeStat <- NULL
  if (!is.null(coordObj@stats[[grid_name]]) &&
    !is.null(coordObj@stats[[grid_name]][[lee_layer]])) {
    leeStat <- coordObj@stats[[grid_name]][[lee_layer]]
  }

  if (is.null(leeStat) && !is.null(g_layer[[lee_layer]])) {
    leeStat <- g_layer[[lee_layer]]
  }

  if (is.null(leeStat) || is.null(leeStat$L)) {
    stop(
      "Layer '", lee_layer, "' in grid '", grid_name,
      "' does not contain a valid Lee’s L matrix."
    )
  }

  Lmat <- leeStat$L

  ## ---- 3. If Pearson correlation matrix exists, take intersection ---------------------------------
  if (!is.null(g_layer$pearson_cor)) {
    common <- intersect(rownames(Lmat), rownames(g_layer$pearson_cor))
    Lmat <- Lmat[common, common, drop = FALSE]
  }
  Lmat
}

#' @noRd
.getLeeMatrix <- function(coordObj,
                          grid_name = NULL,
                          lee_layer = NULL) {
  ## ---- 0. Select grid sub-layer ------------------------------------------------
  g_layer <- .selectGridLayer(coordObj, grid_name)
  if (is.null(grid_name)) { # Write back the actual name
    grid_name <- names(coordObj@grid)[
      vapply(coordObj@grid, identical, logical(1), g_layer)
    ]
  }

  ## ---- 1. Auto-detect LeeStats layer name ----------------------------------------
  if (is.null(lee_layer)) {
    cand <- character(0)
    if (!is.null(coordObj@stats[[grid_name]])) {
      cand <- names(coordObj@stats[[grid_name]])
    }
    cand <- c(cand, names(g_layer))
    cand <- unique(cand[grepl("^LeeStats_", cand)])

    if (length(cand) == 0L) {
      stop(
        "No layer starting with 'LeeStats_' found for grid '",
        grid_name, "'."
      )
    }
    if (length(cand) == 1L) {
      lee_layer <- cand
    } else if ("LeeStats_Xz" %in% cand) { # Prefer default naming
      lee_layer <- "LeeStats_Xz"
    } else {
      stop(
        "Multiple LeeStats layers detected (",
        paste(cand, collapse = ", "),
        "); please specify `lee_layer` explicitly."
      )
    }
  }

  ## ---- 2. Search @stats → @grid in order ------------------------------------
  leeStat <- NULL
  if (!is.null(coordObj@stats[[grid_name]]) &&
    !is.null(coordObj@stats[[grid_name]][[lee_layer]])) {
    leeStat <- coordObj@stats[[grid_name]][[lee_layer]]
  }

  if (is.null(leeStat) && !is.null(g_layer[[lee_layer]])) {
    leeStat <- g_layer[[lee_layer]]
  }

  if (is.null(leeStat) || is.null(leeStat$L)) {
    stop(
      "Layer '", lee_layer, "' in grid '", grid_name,
      "' does not contain a valid Lee’s L matrix."
    )
  }

  Lmat <- leeStat$L

  ## ---- 3. If Pearson correlation matrix exists, take intersection ---------------------------------
  if (!is.null(g_layer$pearson_cor)) {
    common <- intersect(rownames(Lmat), rownames(g_layer$pearson_cor))
    Lmat <- Lmat[common, common, drop = FALSE]
  }
  Lmat
}

#' @title Get Pearson Correlation Matrix
#' @description
#'   Retrieves the gene-gene Pearson correlation matrix either at the grid
#'   level (aggregated grids) or the single-cell level from a \code{CoordObj}.
#'   If a Lee’s L matrix is present, the result is subset to the intersecting
#'   gene set for comparability.
#'
#' @param coordObj   A \code{CoordObj} with correlation results.
#' @param grid_name  Grid sub-layer name (ignored when \code{level = "cell"}).
#' @param level      \code{"grid"} (Default) or \code{"cell"}.
#'
#' @return A numeric correlation matrix with genes in both rows / columns.
#' @export
getPearsonMatrix <- function(coordObj,
                             grid_name = NULL,
                             level = c("grid", "cell")) {
  level <- match.arg(level)

  ## ---------- 1. Set layer name & final target ------------------------------------
  corr_name <- "pearson_cor"
  f_cell_suf <- "_cell" # single-cell layer suffix

  ## ---------- 2. Get matrix --------------------------------------------------
  if (level == "grid") {
    g_layer <- .selectGridLayer(coordObj, grid_name)
    if (is.null(grid_name)) {
      grid_name <- names(coordObj@grid)[
        vapply(coordObj@grid, identical, logical(1), g_layer)
      ]
    }

    ##   2a. New version: @stats[[grid_name]]
    rmat <- if (!is.null(coordObj@stats[[grid_name]]) &&
      !is.null(coordObj@stats[[grid_name]][[corr_name]])) {
      coordObj@stats[[grid_name]][[corr_name]]
    } else {
      NULL
    }

    ##   2b. Fallback: @grid[[grid_name]]
    if (is.null(rmat) && !is.null(g_layer[[corr_name]])) {
      rmat <- g_layer[[corr_name]]
    }

    if (is.null(rmat)) {
      stop("Pearson matrix not found for grid layer '", grid_name, "'.")
    }

    ## Lee’s L alignment
    Lmat <- tryCatch(
      .getLeeMatrix(coordObj, grid_name = grid_name),
      error = function(e) NULL
    )
    if (!is.null(Lmat)) {
      common <- intersect(rownames(rmat), rownames(Lmat))
      rmat <- rmat[common, common, drop = FALSE]
    }
  } else { # ---------- single-cell ----------

    ##   2a. New version: @stats[["cell"]]
    rmat <- if (!is.null(coordObj@stats[["cell"]]) &&
      !is.null(coordObj@stats[["cell"]][[paste0(corr_name, f_cell_suf)]])) {
      coordObj@stats[["cell"]][[paste0(corr_name, f_cell_suf)]]
    } else {
      NULL
    }

    ##   2b. Fallback: @cells$pearson_cor
    if (is.null(rmat) && !is.null(coordObj@cells[[corr_name]])) {
      rmat <- coordObj@cells[[corr_name]]
    }

    if (is.null(rmat)) {
      stop("Pearson matrix at cell level not found in coordObj.")
    }

    ## Lee’s L alignment (if single-cell layer exists)
    Lmat <- NULL
    if (!is.null(coordObj@stats[["cell"]])) {
      layer_L <- intersect(
        grep("^LeeStats_", names(coordObj@stats[["cell"]]),
          value = TRUE
        ),
        paste0("LeeStats_Xz", f_cell_suf)
      )
      if (length(layer_L)) {
        Lmat <- coordObj@stats[["cell"]][[layer_L[1]]]$L
      }
    }
    if (is.null(Lmat) && !is.null(coordObj@cells$LeeStats_Xz)) {
      Lmat <- coordObj@cells$LeeStats_Xz$L
    }

    if (!is.null(Lmat)) {
      common <- intersect(rownames(rmat), rownames(Lmat))
      rmat <- rmat[common, common, drop = FALSE]
    }
  }

  return(rmat)
}

#' @noRd
.getPearsonMatrix <- function(coordObj,
                              grid_name = NULL,
                              level = c("grid", "cell")) {
  level <- match.arg(level)

  ## ---------- 1. Set layer name & final target ------------------------------------
  corr_name <- "pearson_cor"
  f_cell_suf <- "_cell" # single-cell layer suffix

  ## ---------- 2. Get matrix --------------------------------------------------
  if (level == "grid") {
    g_layer <- .selectGridLayer(coordObj, grid_name)
    if (is.null(grid_name)) {
      grid_name <- names(coordObj@grid)[
        vapply(coordObj@grid, identical, logical(1), g_layer)
      ]
    }

    ##   2a. New version: @stats[[grid_name]]
    rmat <- if (!is.null(coordObj@stats[[grid_name]]) &&
      !is.null(coordObj@stats[[grid_name]][[corr_name]])) {
      coordObj@stats[[grid_name]][[corr_name]]
    } else {
      NULL
    }

    ##   2b. Fallback: @grid[[grid_name]]
    if (is.null(rmat) && !is.null(g_layer[[corr_name]])) {
      rmat <- g_layer[[corr_name]]
    }

    if (is.null(rmat)) {
      stop("Pearson matrix not found for grid layer '", grid_name, "'.")
    }

    ## Lee’s L alignment
    Lmat <- tryCatch(
      .getLeeMatrix(coordObj, grid_name = grid_name),
      error = function(e) NULL
    )
    if (!is.null(Lmat)) {
      common <- intersect(rownames(rmat), rownames(Lmat))
      rmat <- rmat[common, common, drop = FALSE]
    }
  } else { # ---------- single-cell ----------

    ##   2a. New version: @stats[["cell"]]
    rmat <- if (!is.null(coordObj@stats[["cell"]]) &&
      !is.null(coordObj@stats[["cell"]][[paste0(corr_name, f_cell_suf)]])) {
      coordObj@stats[["cell"]][[paste0(corr_name, f_cell_suf)]]
    } else {
      NULL
    }

    ##   2b. Fallback: @cells$pearson_cor
    if (is.null(rmat) && !is.null(coordObj@cells[[corr_name]])) {
      rmat <- coordObj@cells[[corr_name]]
    }

    if (is.null(rmat)) {
      stop("Pearson matrix at cell level not found in coordObj.")
    }

    ## Lee’s L alignment (if single-cell layer exists)
    Lmat <- NULL
    if (!is.null(coordObj@stats[["cell"]])) {
      layer_L <- intersect(
        grep("^LeeStats_", names(coordObj@stats[["cell"]]),
          value = TRUE
        ),
        paste0("LeeStats_Xz", f_cell_suf)
      )
      if (length(layer_L)) {
        Lmat <- coordObj@stats[["cell"]][[layer_L[1]]]$L
      }
    }
    if (is.null(Lmat) && !is.null(coordObj@cells$LeeStats_Xz)) {
      Lmat <- coordObj@cells$LeeStats_Xz$L
    }

    if (!is.null(Lmat)) {
      common <- intersect(rownames(rmat), rownames(Lmat))
      rmat <- rmat[common, common, drop = FALSE]
    }
  }

  return(rmat)
}

#' @noRd
.leeL_perm_block <- function(Xz, W, L_ref,
                             block_id,
                             perms = 999,
                             block_size = 64,
                             n_threads = 1) {
  stopifnot(
    is.matrix(Xz), inherits(W, "dgCMatrix"),
    length(block_id) == nrow(Xz)
  )

  RhpcBLASctl::blas_set_num_threads(1)
  ngen <- ncol(Xz)
  geCnt <- matrix(0, ngen, ngen)
  done <- 0L

  ## ---- Pre-group row indices by block ----
  split_rows <- split(seq_along(block_id), block_id)
  blk_keys <- names(split_rows)
  n_blk <- length(split_rows)

  while (done < perms) {
    bsz <- min(block_size, perms - done)

    idx_mat <- replicate(bsz,
      {
        new_order <- sample(n_blk) # shuffle block order
        unlist(split_rows[new_order], use.names = FALSE)
      },
      simplify = "matrix"
    ) - 1L # 0-based for C++

    storage.mode(idx_mat) <- "integer"
    geCnt <- geCnt + lee_perm(Xz, W, idx_mat, L_ref, n_threads)
    done <- done + bsz
  }
  (geCnt + 1) / (perms + 1)
}

#' @noRd
.assign_block_id <- function(grid_info, block_side = 8) {
  # gx / gy start at 1
  bx <- (grid_info$gx - 1L) %/% block_side
  by <- (grid_info$gy - 1L) %/% block_side
  # merge into a single integer id
  max_by <- max(by)
  block_id <- bx * (max_by + 1L) + by + 1L
  block_id
}

#' @noRd
.pick_grid_layer <- function(coordObj, grid_name = NULL) {
  stopifnot(!is.null(coordObj@grid))
  if (is.null(grid_name)) {
    sub_layers <- names(coordObj@grid)
    if (length(sub_layers) != 1L) {
      stop("coordObj@grid contains multiple sub-layers; please specify grid_name explicitly")
    }
    sub_layers[[1]]
  } else {
    as.character(grid_name)
  }
}

#' @noRd
# .computeLeeL <- function(coordObj,
#                          grid_name = NULL,
#                          norm_layer = "Xz",
#                          genes = NULL,
#                          within = TRUE,
#                          ncore = 1,
#                          mem_limit_GB = 16,
#                          chunk_size = 256L,
#                          use_bigmemory = TRUE,
#                          backing_path = tempdir()) {
#   if (use_bigmemory && !requireNamespace("bigmemory", quietly = TRUE)) {
#     stop("Package 'bigmemory' is required.")
#   }

#   ## ---- 0. Get grid layer (new helper) ---------------------------------
#   g_layer <- selectGridLayer(coordObj, grid_name)
#   # Fill back layer name string, for later writing back to stats
#   if (is.null(grid_name)) {
#     grid_name <- names(coordObj@grid)[
#       vapply(coordObj@grid, identical, logical(1), g_layer)
#     ]
#   }

#   ## ---- 1. Extract expression matrix and weight matrix ----------------------------------
#   Xz_full <- as.matrix(g_layer[[norm_layer]])
#   ord <- match(g_layer$grid_info$grid_id, rownames(Xz_full))
#   Xz_full <- Xz_full[ord, , drop = FALSE]
#   W <- g_layer$W[ord, ord]


#   all_genes <- colnames(Xz_full)
#   idx_keep <- if (is.null(genes)) {
#     seq_along(all_genes)
#   } else {
#     m <- match(genes, all_genes, nomatch = 0L)
#     if (any(m == 0L)) stop("Some genes were not found in column names")
#     m
#   }

#   n_g <- length(idx_keep)
#   bytes_L <- n_g^2 * 8 / 1024^3 # in GB
#   need_stream <- use_bigmemory && bytes_L > mem_limit_GB

#   ## ======================================================
#   ## =============  A. One‑shot computation (fits in RAM) ===============
#   ## ======================================================
#   if (!need_stream) {
#     RhpcBLASctl::blas_set_num_threads(1)
#     Sys.setenv(OMP_NUM_THREADS = ncore)
#     L_full <- leeL_cpp_cache(Xz_full, W, n_threads = ncore)

#     if (within) {
#       Lmat <- L_full[idx_keep, idx_keep, drop = FALSE]
#       Xuse <- Xz_full[, idx_keep, drop = FALSE]
#     } else {
#       Lmat <- L_full[idx_keep, , drop = FALSE]
#       Xuse <- Xz_full[, idx_keep, drop = FALSE]
#     }

#     ## ======================================================
#     ## =============  B. Chunked computation with file mapping ================
#     ## ======================================================
#   } else {
#     message(sprintf(
#       "Estimated matrix %.1f GB > limit %.1f GB; switching to chunked, streamed write …",
#       bytes_L, mem_limit_GB
#     ))

#     bm_file <- file.path(
#       backing_path,
#       sprintf("LeeL_%s.bin", grid_name)
#     )
#     bm_desc <- file.path(
#       backing_path,
#       sprintf("LeeL_%s.desc", grid_name)
#     )

#     ## *** Key change: init=NULL, dimnames=NULL ***
#     L_bm <- bigmemory::filebacked.big.matrix(
#       nrow = n_g, ncol = n_g, type = "double",
#       backingfile = basename(bm_file),
#       descriptorfile = basename(bm_desc),
#       backingpath = backing_path,
#       init = NULL, # ← do not initialize full matrix
#       dimnames = NULL,
#       shared = TRUE
#     )

#     RhpcBLASctl::blas_set_num_threads(1)
#     Sys.setenv(OMP_NUM_THREADS = ncore)

#     ## ---- Process column blocks and write ----
#     for (start in seq(1L, n_g, by = chunk_size)) {
#       idx_chunk <- idx_keep[start:min(n_g, start + chunk_size - 1L)]
#       L_block <- leeL_cpp_cols(
#         Xz_full, W,
#         cols0 = as.integer(idx_chunk - 1L),
#         n_threads = ncore
#       )

#       ## Write column block + symmetric columns
#       L_bm[, idx_chunk] <- L_block
#       L_bm[idx_chunk, ] <- t(L_block)

#       rm(L_block)
#       gc(verbose = FALSE)
#     }

#     ## ---- Add gene names (only two vectors) ----
#     dimnames(L_bm) <- list(
#       all_genes[idx_keep],
#       all_genes[idx_keep]
#     )

#     Lmat <- L_bm
#     Xuse <- Xz_full[, idx_keep, drop = FALSE]
#   }

#   dimnames(Lmat) <- list(
#     row = colnames(Xuse),
#     col = colnames(Xuse)
#   )

#   list(
#     Lmat      = Lmat,
#     X_used    = Xuse,
#     X_full    = Xz_full,
#     cells     = rownames(Xz_full),
#     W         = W,
#     grid_info = g_layer$grid_info,
#     grid_name = grid_name # ← Additional return, convenient for addLeeStats
#   )
# }
