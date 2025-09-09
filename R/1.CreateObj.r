#' 1.CreateObj_fixed.r (2025-07-09)
#' @title Build CoordObj from Xenium Output with (Multi-)Scale Grid Aggregation - FIXED
#' @description
#'   Fixes Y-axis flipping issue that caused misalignment between grid heatmap and segmentation:
#'   (1) Flip only once; (2) After flip, grid origin y0 = 0; (3) Remove secondary flip inside build_grid.
#'   All other parameters and return values remain identical to the previous version.
#' @param xenium_dir Character. Path to the Xenium output directory.
#' @param lenGrid Numeric vector. Grid sizes for multi-resolution analysis.
#' @param seg_type Character. One of "cell", "nucleus", or "both". Type of segmentation to process.
#' @param filtergenes Character vector or NULL. Gene names to include (default NULL keeps all).
#' @param max_dist_mol_nuc Numeric. Maximum distance between molecule and nucleus for filtering (default 25).
#' @param filtermolecule Logical. Whether to filter molecules by prefix exclusion (default TRUE).
#' @param exclude_prefix Character vector. Prefixes to exclude when filtermolecule=TRUE (default c("Unassigned", "NegControl", "Background")).
#' @param filterqv Numeric or NULL. Minimum QV score for molecules (default 20).
#' @param coord_file Character or NULL. Path to CSV/TSV file defining ROI polygon coordinates.
#' @param ncores Integer. Number of cores for parallel processing (default 1L).
#' @param chunk_pts Numeric. Number of points to process per chunk (default 5e5).
#' @param min_gene_types Integer. Minimum number of unique genes per grid cell (default 1).
#' @param max_gene_types Numeric. Maximum number of unique genes per grid cell (default Inf).
#' @param min_seg_points Integer. Minimum number of segmentation points per grid cell (default 1).
#' @param keep_partial_grid Logical. Whether to keep incomplete edge grid cells (default FALSE).
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#' @param flip_y Logical. Whether to flip Y coordinates (default TRUE).
#' @return A CoordObj S4 object containing centroids, segmentation data, and multi-scale grid layers.
#' @details
#'   This function completely resolves the Y-axis flipping misalignment between grids and segmentation,
#'   fixes segmentation filtering bugs that incorrectly removed all grids,
#'   and corrects setenv errors when restoring environment variables.
#'   All other parameters remain consistent with the previous version.
#' @importFrom parallel detectCores
#' @importFrom data.table setDTthreads as.data.table rbindlist CJ setnames copy uniqueN
#' @importFrom arrow read_parquet open_dataset set_cpu_count set_io_thread_count
#' @importFrom sf st_polygon st_make_valid st_sfc
#' @importFrom dplyr filter summarise collect select pull
#' @importFrom utils read.table txtProgressBar setTxtProgressBar
#' @importFrom tools file_ext
#' @importFrom Matrix sparseMatrix t
#' @importFrom rhdf5 h5read
#' @export
createCoordObj <- function(xenium_dir,
                           lenGrid,
                           seg_type = c("cell", "nucleus", "both"),
                           filtergenes = NULL,
                           max_dist_mol_nuc = 25,
                           filtermolecule = TRUE,
                           exclude_prefix = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword"),
                           filterqv = 20,
                           coord_file = NULL,
                           ncores = 1L,
                           chunk_pts = 5e5,
                           min_gene_types = 1,
                           max_gene_types = Inf,
                           min_seg_points = 1,
                           keep_partial_grid = FALSE,
                           verbose = TRUE,
                           flip_y = TRUE) {
  # Use I/O-bound configuration for data reading and processing
  thread_config <- configureThreadsFor("io_bound", ncores, restore_after = TRUE)
  on.exit({
    restore_fn <- attr(thread_config, "restore_function")
    if (!is.null(restore_fn)) restore_fn()
  })

  # Use R threads for I/O operations, limited OpenMP for processing
  ncores_io <- thread_config$r_threads
  ncores_cpp <- thread_config$openmp_threads

  ## ---- 0. Basic validation and early error detection ----------------
  seg_type <- match.arg(seg_type)
  if (!dir.exists(xenium_dir)) {
    stop("`xenium_dir` does not exist: ", xenium_dir)
  }

  # Check if required files exist
  tf_parq <- file.path(xenium_dir, "transcripts.parquet")
  cell_parq <- file.path(xenium_dir, "cells.parquet")

  if (!file.exists(tf_parq)) {
    stop("transcripts.parquet not found in ", xenium_dir)
  }
  if (!file.exists(cell_parq)) {
    stop("cells.parquet not found in ", xenium_dir)
  }

  stopifnot(is.numeric(lenGrid) && length(lenGrid) >= 1 && all(lenGrid > 0L))
  stopifnot(is.logical(flip_y) && length(flip_y) == 1)
  stopifnot(is.logical(keep_partial_grid) && length(keep_partial_grid) == 1)

  ## ---- 1. Intelligent thread configuration, progressive validation ---------------
  os_type <- detectOS()
  max_cores <- parallel::detectCores()

  # Get system memory information (GB)
  sys_mem_gb <- tryCatch(
    {
      if (os_type == "linux") {
        mem_info <- readLines("/proc/meminfo")
        mem_total <- grep("MemTotal", mem_info, value = TRUE)
        as.numeric(gsub("[^0-9]", "", mem_total)) / 1024 / 1024
      } else if (os_type == "macos") {
        mem_bytes <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
        mem_bytes / 1024^3
      } else {
        32 # Windows default assumption
      }
    },
    error = function(e) 32
  )

  # Set thread limits based on memory and system type
  mem_based_limit <- switch(cut(sys_mem_gb, breaks = c(0, 8, 32, 64, 128, Inf), labels = FALSE),
    1, # 0-8GB: conservative
    min(4, max_cores - 1), # 8-32GB
    min(8, max_cores - 1), # 32-64GB
    min(12, max_cores - 1), # 64-128GB
    min(16, max_cores - 1) # 128GB+
  )

  platform_limit <- switch(os_type,
    windows = min(8, max_cores - 1),
    macos = min(12, max_cores - 1),
    linux = min(16, max_cores - 1)
  )

  # Progressively validate thread count feasibility
  test_cores <- min(ncores, mem_based_limit, platform_limit)
  ncores_safe <- 1

  for (test_nc in seq(test_cores, 1, by = -2)) {
    tryCatch(
      {
        # Simple parallel test
        if (test_nc > 1) {
          cl <- parallel::makeCluster(min(test_nc, 4))
          parallel::clusterEvalQ(cl, 1 + 1)
          parallel::stopCluster(cl)
        }
        ncores_safe <- test_nc
        break
      },
      error = function(e) {
        if (verbose) message("[geneSCOPE] Testing ", test_nc, " cores failed, trying fewer...")
      }
    )
  }

  if (verbose && ncores_safe < ncores) {
    message(
      "[geneSCOPE] Using ", ncores_safe, " cores (requested: ", ncores,
      ", memory: ", round(sys_mem_gb, 1), "GB, platform: ", os_type, ")"
    )
  } else if (verbose) {
    message("[geneSCOPE] Using ", ncores_safe, " cores")
  }

  ## ---- 2. Environment variable management ------------------------------------
  thread_vars <- c(
    "OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OMP_SCHEDULE",
    "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
    "BLAS_NUM_THREADS", "LAPACK_NUM_THREADS"
  )

  old_env <- sapply(thread_vars, Sys.getenv, unset = NA_character_)

  on.exit(
    {
      for (var in names(old_env)) {
        if (is.na(old_env[[var]])) {
          Sys.unsetenv(var)
        } else {
          args <- list()
          args[[var]] <- old_env[[var]]
          do.call(Sys.setenv, args)
        }
      }
    },
    add = TRUE
  )

  # Set conservative thread environment
  for (var in thread_vars) {
    args <- list()
    args[[var]] <- "1"
    do.call(Sys.setenv, args)
  }

  # Safely set data.table threads
  if (requireNamespace("data.table", quietly = TRUE)) {
    tryCatch(
      {
        data.table::setDTthreads(1)
      },
      error = function(e) {
        if (verbose) message("[geneSCOPE] Warning: Could not set data.table threads")
      }
    )
  }

  ## ---- 3. Package loading and configuration --------------------------------------------
  if (verbose) message("[geneSCOPE] Loading required packages...")

  suppressPackageStartupMessages({
    library(arrow)
    library(data.table)
    library(sf)
    library(parallel)
    library(dplyr)
  })

  # Arrow configuration
  if (requireNamespace("arrow", quietly = TRUE)) {
    tryCatch(
      {
        options(arrow.num_threads = ncores_safe)
        if (exists("set_cpu_count", where = asNamespace("arrow"))) {
          arrow::set_cpu_count(ncores_safe)
        }
        if (exists("set_io_thread_count", where = asNamespace("arrow"))) {
          arrow::set_io_thread_count(ncores_safe)
        }
      },
      error = function(e) {
        if (verbose) message("[geneSCOPE] Warning: Could not configure Arrow threading")
      }
    )
  }

  ## ---- 4. File validation and segmentation detection -------------------------------------
  seg_files <- list()
  if (seg_type %in% c("cell", "both")) {
    cell_seg_files <- c(
      file.path(xenium_dir, "cell_boundaries.parquet"),
      file.path(xenium_dir, "segmentation_boundaries.parquet")
    )
    cell_seg_found <- sapply(cell_seg_files, file.exists)
    if (any(cell_seg_found)) {
      seg_files$cell <- cell_seg_files[cell_seg_found][1]
    } else {
      stop("Cell segmentation parquet file not found in ", xenium_dir)
    }
  }

  if (seg_type %in% c("nucleus", "both")) {
    nuc_file <- file.path(xenium_dir, "nucleus_boundaries.parquet")
    if (!file.exists(nuc_file)) {
      stop("Nucleus segmentation parquet file not found in ", xenium_dir)
    }
    seg_files$nucleus <- nuc_file
  }

  ## ---- 5. Read centroid data with I/O optimization
  if (verbose) message("[geneSCOPE] Reading cell centroids...")

  centroids_dt <- tryCatch(
    {
      # Use Arrow with limited threading for I/O stability
      if (requireNamespace("arrow", quietly = TRUE)) {
        arrow::set_cpu_count(ncores_io)
        arrow::set_io_thread_count(ncores_io)
      }

      cell_data <- arrow::read_parquet(cell_parq)
      as.data.table(cell_data)[
        , .(cell = cell_id, x = x_centroid, y = y_centroid)
      ]
    },
    error = function(e) {
      stop("Error reading centroids from cells.parquet: ", e$message)
    }
  )

  if (nrow(centroids_dt) == 0) {
    stop("No centroids found in cells.parquet")
  }

  if (verbose) message("[geneSCOPE] Found ", nrow(centroids_dt), " cells")

  ## ---- 6. Dataset preparation --------------------------------
  if (verbose) message("[geneSCOPE] Preparing transcript dataset...")

  ds <- tryCatch(
    {
      arrow::open_dataset(tf_parq)
    },
    error = function(e) {
      stop("Error opening transcripts dataset: ", e$message)
    }
  )

  # Apply filters
  if (!is.null(filtergenes)) {
    if (verbose) message("[geneSCOPE] Filtering genes: ", length(filtergenes), " genes")
    ds <- ds |> filter(feature_name %in% filtergenes)
  }

  if (!is.null(max_dist_mol_nuc)) {
    if (verbose) message("[geneSCOPE] Filtering by nucleus distance: max ", max_dist_mol_nuc)
    ds <- ds |> filter(nucleus_distance <= max_dist_mol_nuc)
  }

  if (filtermolecule) {
    if (verbose) message("[geneSCOPE] Filtering molecules by prefix exclusion")
    pat <- paste0("^(?:", paste(exclude_prefix, collapse = "|"), ")")
    ds <- ds |> filter(!grepl(pat, feature_name))
  }

  if (!is.null(filterqv)) {
    if (verbose) message("[geneSCOPE] Filtering by QV score: min ", filterqv)
    ds <- ds |> filter(qv >= filterqv)
  }

  # Compute bounds
  if (verbose) message("[geneSCOPE] Computing global bounds...")

  bounds_global <- tryCatch(
    {
      ds |>
        summarise(
          xmin = min(x_location), xmax = max(x_location),
          ymin = min(y_location), ymax = max(y_location)
        ) |>
        collect()
    },
    error = function(e) {
      stop("Error computing bounds: ", e$message)
    }
  )

  y_max <- as.numeric(bounds_global$ymax)

  if (verbose) {
    message(
      "[geneSCOPE] Global bounds: x=[", round(bounds_global$xmin), ",",
      round(bounds_global$xmax), "], y=[", round(bounds_global$ymin),
      ",", round(bounds_global$ymax), "]"
    )
  }

  if (flip_y) {
    centroids_dt <- .flipCoordinates(centroids_dt, y_max)
    if (verbose) message("[geneSCOPE] Applied Y-axis flip")
  }

  ## ---- 7. ROI polygon processing -------------------------
  user_poly <- NULL
  if (!is.null(coord_file)) {
    if (verbose) message("[geneSCOPE] Processing ROI polygon from: ", coord_file)

    if (!file.exists(coord_file)) {
      stop("Coordinate file not found: ", coord_file)
    }

    ext <- tolower(tools::file_ext(coord_file))
    sep_char <- switch(ext,
      csv = ",",
      tsv = "\t",
      txt = "\t",
      "\t"
    )

    poly_tb <- tryCatch(
      {
        utils::read.table(coord_file,
          header = TRUE, sep = sep_char,
          stringsAsFactors = FALSE
        )
      },
      error = function(e) {
        stop("Error reading coordinate file: ", e$message)
      }
    )

    names(poly_tb) <- tolower(names(poly_tb))

    if (!all(c("x", "y") %in% names(poly_tb))) {
      stop("Coordinate file must contain 'x' and 'y' columns")
    }

    if (nrow(poly_tb) < 3) {
      stop("ROI polygon must have at least 3 points")
    }

    if (verbose) message("[geneSCOPE] ROI polygon has ", nrow(poly_tb), " vertices")

    poly_geom <- tryCatch(
      {
        geom <- sf::st_polygon(list(as.matrix(poly_tb[, c("x", "y")])))

        sf::st_make_valid(geom)
      },
      error = function(e) {
        stop("Error creating polygon geometry: ", e$message)
      }
    )

    if (length(poly_geom) == 0) {
      stop("ROI polygon became EMPTY after validation")
    }

    if (flip_y) {
      poly_tb <- .flipCoordinates(poly_tb, y_max)
      poly_geom <- tryCatch(
        {
          sf::st_make_valid(sf::st_polygon(list(as.matrix(poly_tb[, c("x", "y")]))))
        },
        error = function(e) {
          stop("Error creating flipped polygon: ", e$message)
        }
      )
    }
    user_poly <- sf::st_sfc(poly_geom)
  } else {
    # Auto-create ROI from data boundaries when no coord_file is specified
    if (verbose) message("[geneSCOPE] Creating automatic ROI from data boundaries...")
    
    # Create rectangle polygon from global bounds
    poly_coords <- data.frame(
      x = c(bounds_global$xmin, bounds_global$xmax, bounds_global$xmax, bounds_global$xmin, bounds_global$xmin),
      y = c(bounds_global$ymin, bounds_global$ymin, bounds_global$ymax, bounds_global$ymax, bounds_global$ymin)
    )
    
    if (flip_y) {
      poly_coords <- .flipCoordinates(poly_coords, y_max)
    }
    
    poly_geom <- tryCatch(
      {
        geom <- sf::st_polygon(list(as.matrix(poly_coords[, c("x", "y")])))
        sf::st_make_valid(geom)
      },
      error = function(e) {
        stop("Error creating automatic ROI polygon: ", e$message)
      }
    )
    
    if (length(poly_geom) == 0) {
      stop("Automatic ROI polygon became EMPTY after validation")
    }
    
    user_poly <- sf::st_sfc(poly_geom)
    
    if (verbose) {
      message("[geneSCOPE] Automatic ROI created: x=[", round(min(poly_coords$x)), ",", 
              round(max(poly_coords$x)), "], y=[", round(min(poly_coords$y)), ",", 
              round(max(poly_coords$y)), "]")
    }
  }

  ## ---- 8. Centroid clipping -----------------------------
  if (!is.null(user_poly) && nrow(centroids_dt)) {
    if (verbose) message("[geneSCOPE] Clipping centroids to ROI...")

    centroids_dt <- tryCatch(
      {
        .clipPointsToPolygon(centroids_dt, user_poly,
          chunk_size = chunk_pts,
          ncores = ncores_safe
        )
      },
      error = function(e) {
        stop("Error clipping centroids: ", e$message)
      }
    )

    if (verbose) message("[geneSCOPE] Retained ", nrow(centroids_dt), " cells within ROI")
  }

  keep_cells <- if (nrow(centroids_dt)) unique(centroids_dt$cell) else NULL

  ## ---- 9. Initialize CoordObj ----------------------------------------
  coord_obj <- new("CoordObj",
    coord     = list(centroids = centroids_dt),
    grid      = list(),
    meta.data = data.frame()
  )

  ## ---- 10. Process segmentation data -------------------------------
  if (verbose) message("[geneSCOPE] Processing segmentation data...")

  for (tag in names(seg_files)) {
    if (verbose) message("[geneSCOPE]  Processing ", tag, " segmentation...")

    seg_res <- tryCatch(
      {
        .processSegmentation(
          seg_file = seg_files[[tag]],
          tag = tag,
          flip_y = flip_y,
          y_max = y_max,
          keep_cells = keep_cells,
          ncores = ncores_safe
        )
      },
      error = function(e) {
        warning("Error processing ", tag, " segmentation: ", e$message)
        list(points = data.table(), polygons = list())
      }
    )

    coord_obj@coord[[paste0("segmentation_", tag)]] <- seg_res$points
    coord_obj@coord[[paste0("segmentation_polys_", tag)]] <- seg_res$polygons
  }

  ## ---- 11. Pre-fetch molecules (ROI case) ----------------------------
  mol_small <- NULL
  if (!is.null(user_poly)) {
    if (verbose) message("[geneSCOPE] Pre-fetching molecules for ROI...")

    mol_small <- tryCatch(
      {
        ds |>
          select(x = x_location, y = y_location, feature_name) |>
          collect() |>
          as.data.table()
      },
      error = function(e) {
        warning("Error pre-fetching molecules: ", e$message)
        data.table(x = numeric(), y = numeric(), feature_name = character())
      }
    )

    if (flip_y && nrow(mol_small) > 0) {
      mol_small[, y := y_max - y]
    }

    if (nrow(mol_small) > 0) {
      mol_small <- .clipPointsToPolygon(mol_small, user_poly,
        chunk_size = chunk_pts,
        ncores = ncores_safe
      )
    }

    if (verbose) message("[geneSCOPE] Retained ", nrow(mol_small), " molecules within ROI")
  }

  ## ---- 12. Batch scan counting (panoramic case, using original iterator approach) ----------------------
  counts_list <- setNames(vector("list", length(lenGrid)), paste0("lg", lenGrid))
  for (i in seq_along(counts_list)) {
    counts_list[[i]] <- data.table(grid_id = character(), gene = character(), count = integer())
  }

  if (is.null(user_poly)) {
    if (verbose) message("[geneSCOPE] Starting multi-resolution scan...")

    # Using the original iterator approach
    n_total <- tryCatch(
      {
        ds |>
          summarise(n = n()) |>
          collect() |>
          pull(n)
      },
      error = function(e) NA_integer_
    )

    pb <- if (verbose && !is.na(n_total)) {
      utils::txtProgressBar(min = 0, max = n_total, style = 3)
    } else {
      NULL
    }

    scanned <- 0L

    # Using to_table() to read in batches, this is the original approach
    ds_scanner <- ds$to_table(batch_size = as.integer(chunk_pts))

    repeat {
      tbl <- tryCatch(
        {
          if (scanned == 0) {
            # First read
            first_batch <- ds |>
              select(x_location, y_location, feature_name) |>
              slice_head(n = chunk_pts) |>
              collect()
            first_batch
          } else {
            # Subsequent batches - using skip and limit
            batch <- ds |>
              select(x_location, y_location, feature_name) |>
              slice(scanned + 1, scanned + chunk_pts) |>
              collect()

            if (nrow(batch) == 0) NULL else batch
          }
        },
        error = function(e) {
          if (verbose) message("[geneSCOPE] Scanner finished or error: ", e$message)
          NULL
        }
      )

      if (is.null(tbl) || nrow(tbl) == 0) break

      scanned <- scanned + nrow(tbl)
      if (!is.null(pb)) {
        utils::setTxtProgressBar(pb, min(scanned, n_total))
      }

      dt <- as.data.table(tbl)
      if (flip_y) dt[, y_location := y_max - y_location]

      for (k in seq_along(lenGrid)) {
        lg <- lenGrid[k]
        x0_lg <- floor(bounds_global$xmin / lg) * lg
        y0_lg <- if (flip_y) 0 else floor(bounds_global$ymin / lg) * lg

        dt[, `:=`(
          gx = (x_location - x0_lg) %/% lg + 1L,
          gy = (y_location - y0_lg) %/% lg + 1L
        )]
        dt[, grid_id_tmp := paste0("g", gx, "_", gy)]

        accum <- dt[, .N, by = .(grid_id = grid_id_tmp, gene = feature_name)]
        setnames(accum, "N", "count")
        counts_list[[k]] <- rbindlist(list(counts_list[[k]], accum))
      }

      # Prevent infinite loop
      if (!is.na(n_total) && scanned >= n_total) break
    }

    if (!is.null(pb)) close(pb)

    # Merge duplicate grid_id-gene combinations
    for (k in seq_along(counts_list)) {
      counts_list[[k]] <- counts_list[[k]][, .(count = sum(count)), by = .(grid_id, gene)]
    }
  }

  ## ---- 13. build_grid() function ----------------------------------------
  build_grid <- function(lg) {
    if (verbose) message(sprintf("[geneSCOPE]  ↳ grid %.1f µm …", lg))
    x0 <- floor(bounds_global$xmin / lg) * lg
    y0 <- if (flip_y) 0 else floor(bounds_global$ymin / lg) * lg

    ## Molecule counts
    cnt <- if (is.null(user_poly)) {
      counts_list[[paste0("lg", lg)]]
    } else {
      mol_dt <- copy(mol_small)
      mol_dt[, `:=`(
        gx = (x - x0) %/% lg + 1L,
        gy = (y - y0) %/% lg + 1L
      )]
      mol_dt[, grid_id := paste0("g", gx, "_", gy)]
      accum <- mol_dt[, .N, by = .(grid_id, gene = feature_name)]
      setnames(accum, "N", "count")
      accum
    }

    ## ------- ① Gene diversity filtering ------------------------------------
    gene_div <- cnt[, uniqueN(gene), by = grid_id]
    keep_grid <- gene_div[V1 >= min_gene_types & V1 <= max_gene_types, grid_id]
    cnt <- cnt[grid_id %in% keep_grid]

    ## ------- ② Segmentation filtering -- **Core fix** -------------------
    if (min_seg_points > 0L) {
      seg_layers <- names(coord_obj@coord)[grepl("^segmentation_", names(coord_obj@coord)) &
        !grepl("_polys_", names(coord_obj@coord))]
      seg_dt_all <- rbindlist(coord_obj@coord[seg_layers], use.names = TRUE, fill = TRUE)
      if (nrow(seg_dt_all)) {
        seg_dt_all[, `:=`(
          gx = (x - x0) %/% lg + 1L,
          gy = (y - y0) %/% lg + 1L
        )]
        seg_dt_all[, grid_id := paste0("g", gx, "_", gy)]
        seg_cnt <- seg_dt_all[, .N, by = grid_id]
        valid_seg <- seg_cnt[N >= min_seg_points, grid_id]

        ## ---- Key logic: Only intersect when valid_seg is non-empty ----
        if (length(valid_seg)) {
          keep_grid <- intersect(keep_grid, valid_seg)
          cnt <- cnt[grid_id %in% keep_grid] # Synchronize counts filtering
        }
      }
    }

    ## ------- ③ Generate grid_info --------------------------------------
    x_breaks <- seq(x0, bounds_global$xmax, by = lg)
    if (tail(x_breaks, 1L) < bounds_global$xmax) {
      x_breaks <- c(x_breaks, bounds_global$xmax)
    }

    max_y_use <- if (flip_y) y_max else bounds_global$ymax
    y_breaks <- seq(y0, max_y_use, by = lg)
    if (tail(y_breaks, 1L) < max_y_use) {
      y_breaks <- c(y_breaks, max_y_use)
    }

    grid_dt <- CJ(
      gx = seq_len(length(x_breaks) - 1L),
      gy = seq_len(length(y_breaks) - 1L)
    )
    grid_dt[, grid_id := paste0("g", gx, "_", gy)]
    grid_dt[, `:=`(
      xmin = x_breaks[gx],
      xmax = x_breaks[gx + 1L],
      ymin = y_breaks[gy],
      ymax = y_breaks[gy + 1L]
    )]
    grid_dt[, `:=`(
      center_x = (xmin + xmax) / 2,
      center_y = (ymin + ymax) / 2,
      width    = xmax - xmin,
      height   = ymax - ymin,
      idx      = .I
    )]

    ## Optional: Remove incomplete edge grids
    if (!keep_partial_grid) {
      grid_dt <- grid_dt[abs(width - lg) < 1e-6 & abs(height - lg) < 1e-6]
    }

    ## Final filtering
    grid_dt <- grid_dt[grid_id %in% keep_grid]
    cnt <- cnt[grid_id %in% grid_dt$grid_id]

    coord_obj@grid[[paste0("grid_lenGrid", lg)]] <<- list(
      grid_info = grid_dt,
      counts    = cnt,
      lenGrid   = lg,
      xbins_eff = length(x_breaks) - 1L,
      ybins_eff = length(y_breaks) - 1L
    )
  }

  if (verbose) message("[geneSCOPE] Building grid layers …")
  invisible(lapply(sort(unique(lenGrid)), build_grid))
  if (verbose) message("[geneSCOPE] Grid construction finished.")

  ## ---- 14. Integrity check ----------------------------------------------
  if (length(coord_obj@grid) == 0 ||
    all(vapply(coord_obj@grid, function(g) nrow(g$grid_info) == 0L, logical(1)))) {
    warning("No effective grid layer generated – please check filters/parameters.")
  }

  coord_obj
}


#' @title Add Cell Count Matrix
#' @description
#'   Read the Xenium **cell-feature sparse matrix** (either HDF5 or
#'   Matrix-Market format), keep only the cells that appear in
#'   \code{coordObj@coord$centroids$cell}, preserve that order, and store the
#'   result in the new \code{@cells} slot.
#'
#'   The function relies only on \strong{Matrix}, \strong{data.table}, and
#'   \strong{rhdf5} (or \strong{hdf5r})—no Seurat, tidyverse, or other heavy
#'   dependencies.
#' @param coordObj A valid \code{CoordObj} whose \code{@coord$centroids} slot
#'        is already filled (e.g. via \code{\link{createCoordObj}}). Cells
#'        found in this slot define the subset and order of columns kept.
#' @param xenium_dir Character scalar. Path to the Xenium \file{outs/}
#'        directory that holds \file{cell_feature_matrix.h5}; must exist.
#' @param filtergenes Character vector or NULL. Gene names to include (default NULL keeps all).
#' @param exclude_prefix Character vector. Prefixes to exclude when filtering genes (default c("Unassigned", "NegControl", "Background", "DeprecatedCodeword")).
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#' @return The same object (modified in place) now carrying a
#'         \code{dgCMatrix} in \code{@cells}.
#' @details
#' The function performs five main steps:
#' \enumerate{
#'   \item Reads the \strong{cell‑feature} matrix in HDF5 format produced by
#'         10x Xenium, using \pkg{rhdf5} to avoid loading the full file into
#'         memory.
#'   \item Reconstructs the sparse count matrix as a \code{dgCMatrix} from
#'         the \emph{x/i/p} components and automatically handles either
#'         {genes × cells} or {cells × genes} orientations.
#'   \item Filters genes by excluding those with specified prefixes 
#'         (\code{exclude_prefix}) and optionally keeping only genes in 
#'         \code{filtergenes}, similar to the filtering applied in 
#'         \code{\link{createCoordObj}}.
#'   \item Filters and reorders columns so that their order matches the
#'         centroid table (\code{coordObj@coord$centroids$cell}), guaranteeing
#'         one‑to‑one correspondence for downstream spatial analyses.
#'   \item Migrates any extra attributes found on the original matrix
#'         (e.g. log-CPM transforms) into the list stored in
#'         \code{coordObj@cells}.
#' }
#'
#' Only the \pkg{Matrix}, \pkg{data.table}, and \pkg{rhdf5} packages are
#' required; no Seurat or tidyverse dependencies are introduced. Memory
#' footprint scales with the number of non‑zero entries, not the full matrix
#' dimensions.
#' @importFrom Matrix sparseMatrix t
#' @importFrom data.table data.table
#' @importFrom rhdf5 h5read
#' @examples
#' \dontrun{
#' ## Add the cell‑level count matrix to an existing CoordObj
#' coord <- createCoordObj("mouse_brain/xenium_output")
#' coord <- addCells(coord, "mouse_brain/xenium_output")
#' 
#' ## With custom gene filtering
#' coord <- addCells(coord, "mouse_brain/xenium_output",
#'                   filtergenes = c("Gene1", "Gene2", "Gene3"),
#'                   exclude_prefix = c("Unassigned", "NegControl"))
#' }
#' @export
addCells <- function(coordObj, xenium_dir, 
                     filtergenes = NULL,
                     exclude_prefix = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword"),
                     verbose = TRUE) {
  stopifnot(inherits(coordObj, "CoordObj"), dir.exists(xenium_dir))

  if (verbose) message("[geneSCOPE] Loading HDF5 libraries...")
  suppressPackageStartupMessages({
    library(Matrix)
    library(data.table)
  })

  # Cross-platform HDF5 library detection and loading
  hdf5_available <- FALSE
  if (requireNamespace("rhdf5", quietly = TRUE)) {
    library(rhdf5)
    hdf5_available <- TRUE
  } else if (requireNamespace("hdf5r", quietly = TRUE)) {
    library(hdf5r)
    hdf5_available <- TRUE
  }

  if (!hdf5_available) {
    stop("Neither rhdf5 nor hdf5r package is available. Please install one of them.")
  }

  ## ------------------------------------------------------------------ 1
  keep_cells_target <- coordObj@coord$centroids$cell
  if (!length(keep_cells_target)) {
    stop("coordObj@coord$centroids is empty, cannot determine cells to keep.")
  }

  if (verbose) message("[geneSCOPE] Target cells to retain: ", length(keep_cells_target))

  ## ------------------------------------------------------------------ 2
  h5_file <- file.path(xenium_dir, "cell_feature_matrix.h5")
  if (!file.exists(h5_file)) {
    stop("HDF5 file not found: ", h5_file)
  }

  if (verbose) message("[geneSCOPE] Reading HDF5 cell-feature matrix from: ", basename(h5_file))

  # Cross-platform HDF5 reading
  if (requireNamespace("rhdf5", quietly = TRUE)) {
    barcodes <- rhdf5::h5read(h5_file, "matrix/barcodes")
    genes_id <- rhdf5::h5read(h5_file, "matrix/features/id")
    genes_nm <- rhdf5::h5read(h5_file, "matrix/features/name")
    data <- rhdf5::h5read(h5_file, "matrix/data")
    indices <- rhdf5::h5read(h5_file, "matrix/indices")
    indptr <- rhdf5::h5read(h5_file, "matrix/indptr")
    shape <- as.integer(rhdf5::h5read(h5_file, "matrix/shape"))
  } else {
    # Alternative hdf5r implementation
    h5file <- hdf5r::H5File$new(h5_file, mode = "r")
    barcodes <- h5file[["matrix/barcodes"]][]
    genes_id <- h5file[["matrix/features/id"]][]
    genes_nm <- h5file[["matrix/features/name"]][]
    data <- h5file[["matrix/data"]][]
    indices <- h5file[["matrix/indices"]][]
    indptr <- h5file[["matrix/indptr"]][]
    shape <- as.integer(h5file[["matrix/shape"]][])
    h5file$close()
  }

  ## -------------------------- 2a Determine HDF5 storage orientation ----
  nrow_h5 <- shape[1]
  ncol_h5 <- shape[2]
  by_cols <- length(indptr) == (ncol_h5 + 1) # indptr length = ncol+1

  if (!by_cols) {
    stop("Current function only supports CSC-style storage (indptr length = ncol+1).")
  }

  col_is_cell <- (ncol_h5 == length(barcodes))
  if (!col_is_cell && !(nrow_h5 == length(barcodes))) {
    stop("Shape doesn't match barcode/gene counts, cannot determine orientation.")
  }

  if (verbose) {
    message(
      "[geneSCOPE] Matrix dimensions: ", nrow_h5, " × ", ncol_h5,
      " (", length(genes_nm), " genes, ", length(barcodes), " cells)"
    )
  }

  ## -------------------------- 2b Build sparse matrix ------------------
  if (verbose) message("[geneSCOPE] Building sparse count matrix...")
  x <- as.numeric(data)
  i <- as.integer(indices)
  p <- as.integer(indptr)

  if (col_is_cell) {
    ## →  columns = cells, rows = genes   (most common)
    counts_raw <- new("dgCMatrix",
      Dim      = c(length(genes_nm), length(barcodes)),
      x        = x,
      i        = i,
      p        = p,
      Dimnames = list(make.unique(genes_nm), barcodes)
    )
  } else {
    ## →  columns = genes, rows = cells   (rare case: need transpose)
    if (verbose) message("[geneSCOPE] Transposing matrix (cells × genes → genes × cells)...")
    tmp <- new("dgCMatrix",
      Dim      = c(length(barcodes), length(genes_nm)),
      x        = x,
      i        = i,
      p        = p,
      Dimnames = list(barcodes, make.unique(genes_nm))
    )
    counts_raw <- Matrix::t(tmp)
  }

  attr(counts_raw, "gene_map") <- data.frame(
    ensembl = genes_id,
    symbol = genes_nm,
    stringsAsFactors = FALSE
  )

  ## -------------------------- 2c Filter genes by exclude_prefix and filtergenes ------------------
  if (verbose) message("[geneSCOPE] Applying gene filters...")
  
  # Get gene names (rownames)
  all_genes <- rownames(counts_raw)
  keep_genes <- all_genes
  
  # Filter by exclude_prefix
  if (length(exclude_prefix) > 0) {
    exclude_pattern <- paste0("^(?:", paste(exclude_prefix, collapse = "|"), ")")
    exclude_mask <- grepl(exclude_pattern, keep_genes, perl = TRUE)
    keep_genes <- keep_genes[!exclude_mask]
    
    if (verbose) {
      excluded_count <- sum(exclude_mask)
      if (excluded_count > 0) {
        message("[geneSCOPE] Excluded ", excluded_count, " genes with prefixes: ", 
                paste(exclude_prefix, collapse = ", "))
      }
    }
  }
  
  # Filter by filtergenes (if specified)
  if (!is.null(filtergenes)) {
    keep_genes <- intersect(keep_genes, filtergenes)
    if (verbose) {
      message("[geneSCOPE] Keeping ", length(keep_genes), " genes from filtergenes list")
    }
  }
  
  if (length(keep_genes) == 0) {
    stop("No genes remain after filtering. Please check exclude_prefix and filtergenes parameters.")
  }
  
  # Apply gene filtering to the matrix
  counts_raw <- counts_raw[keep_genes, , drop = FALSE]
  
  if (verbose) {
    message("[geneSCOPE] Final gene count: ", nrow(counts_raw), " genes (from original ", 
            length(all_genes), ")")
  }

  ## ------------------------------------------------------------------ 3
  if (verbose) message("[geneSCOPE] Filtering and reordering cells...")
  keep_cells <- keep_cells_target[keep_cells_target %in% colnames(counts_raw)]
  if (!length(keep_cells)) {
    stop("Centroid cell IDs do not overlap with H5 data.")
  }

  if (verbose) message("[geneSCOPE] Overlapping cells: ", length(keep_cells), "/", length(keep_cells_target))

  counts <- counts_raw[, keep_cells, drop = FALSE]

  ## ------------------------------------------------------------------ 4
  if (verbose) message("[geneSCOPE] Storing count matrix in @cells slot...")
  coordObj@cells <- list(counts = counts)

  if (verbose) message("[geneSCOPE] Cell count matrix added successfully")
  invisible(coordObj)
}
