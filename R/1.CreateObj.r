#' 1.CreateObj_fixed.r (2025-07-09)
#' @title Build scope_object from Xenium Output with (Multi-)Scale Grid Aggregation - FIXED
#' @description
#'   Fixes Y-axis flipping issue that caused misalignment between grid heatmap and segmentation:
#'   (1) Flip only once; (2) After flip, grid origin y0 = 0; (3) Remove secondary flip inside build_grid.
#'   All other parameters and return values remain identical to the previous version.
#' @param xenium_dir Character. Path to the Xenium output directory.
#' @param grid_length Numeric vector. Grid sizes for multi-resolution analysis.
#' @param seg_type Character. One of "cell", "nucleus", or "both". Type of segmentation to process.
#' @param filter_genes Character vector or NULL. Gene names to include (default NULL keeps all).
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
#' @return A scope_object S4 object containing centroids, segmentation data, and multi-scale grid layers.
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
createSCOPE_xenium <- function(xenium_dir,
                        grid_length,
                        seg_type = c("cell", "nucleus", "both"),
                        filter_genes = NULL,
                        max_dist_mol_nuc = 25,
                        filtermolecule = TRUE,
                        exclude_prefix = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword",
                                           "SystemControl", "Negative"),
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

  stopifnot(is.numeric(grid_length) && length(grid_length) >= 1 && all(grid_length > 0L))
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
        mem_kb <- as.numeric(gsub("[^0-9]", "", mem_total))
        mem_gb <- mem_kb / 1024 / 1024
        # Add bounds checking to catch unrealistic values - raised threshold for HPC systems
        if (mem_gb > 50000) {
          warning("Detected unusually high memory: ", round(mem_gb, 1), "GB. This may indicate a parsing error.")
        }
        mem_gb
      } else if (os_type == "macos") {
        mem_bytes <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
        mem_bytes / 1024^3
      } else {
        32 # Windows default assumption
      }
    },
    error = function(e) 32
  )

  # Set thread limits based on memory and system type (HPC-optimized)
  mem_based_limit <- if (sys_mem_gb > 500) {
    # HPC environment with high memory - use aggressive scaling
    min(ncores, max_cores - 2)
  } else {
    switch(cut(sys_mem_gb, breaks = c(0, 8, 32, 64, 128, Inf), labels = FALSE),
      1, # 0-8GB: conservative
      min(4, max_cores - 1), # 8-32GB
      min(8, max_cores - 1), # 32-64GB
      min(16, max_cores - 1), # 64-128GB
      min(ncores, max_cores - 2) # 128GB+ (typically HPC)
    )
  }

  platform_limit <- switch(os_type,
    windows = min(8, max_cores - 1),
    macos = min(12, max_cores - 1),
    linux = if (sys_mem_gb > 100) {
      # HPC Linux environment - much more aggressive scaling
      min(ncores, max_cores - 2)
    } else {
      min(32, max_cores - 1)
    }
  )

  # Progressively validate thread count feasibility (HPC-optimized)
  test_cores <- min(ncores, mem_based_limit, platform_limit)
  ncores_safe <- 1

  # For HPC environments, use more aggressive testing
  test_step <- if (sys_mem_gb > 100 || max_cores > 32) 8 else 2

  for (test_nc in seq(test_cores, 1, by = -test_step)) {
    tryCatch(
      {
        # Simplified parallel test for HPC
        if (test_nc > 1) {
          if (sys_mem_gb > 100) {
            # In HPC environments with >100GB RAM, trust the system more
            ncores_safe <- test_nc
            break
          } else {
            # Standard test for smaller systems
            cl <- parallel::makeCluster(min(test_nc, 4))
            parallel::clusterEvalQ(cl, 1 + 1)
            parallel::stopCluster(cl)
          }
        }
        ncores_safe <- test_nc
        break
      },
      error = function(e) {
        if (verbose) message("[geneSCOPE::createSCOPE xenium] Testing ", test_nc, " cores failed, trying fewer")
      }
    )
  }

  if (verbose) {
    # Format memory display appropriately
    mem_display <- if (sys_mem_gb >= 1024) {
      paste0(round(sys_mem_gb / 1024, 1), "TB")
    } else {
      paste0(round(sys_mem_gb, 1), "GB")
    }

    if (ncores_safe < ncores) {
      message("[geneSCOPE::createSCOPE xenium] Core configuration: using ", ncores_safe, "/", ncores, " cores (", mem_display, " RAM, ", os_type, ")")
    } else {
      message("[geneSCOPE::createSCOPE xenium] Core configuration: using ", ncores_safe, " cores")
    }
  }

  ## ---- 2. Environment variable management ------------------------------------
  thread_vars <- c(
    "OMP_NUM_THREADS", "OMP_THREAD_LIMIT",
    "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
    "BLAS_NUM_THREADS", "LAPACK_NUM_THREADS"
  )

  old_env <- sapply(thread_vars, Sys.getenv, unset = NA_character_)
  old_schedule <- Sys.getenv("OMP_SCHEDULE", unset = NA_character_)

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
      # Restore OMP_SCHEDULE separately
      if (is.na(old_schedule)) {
        Sys.unsetenv("OMP_SCHEDULE")
      } else {
        Sys.setenv(OMP_SCHEDULE = old_schedule)
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

  # Set valid OpenMP schedule only for Linux (where it's more stable)
  # macOS OpenMP can be problematic with schedule settings
  if (os_type == "linux") {
    Sys.setenv(OMP_SCHEDULE = "static")
  }

  # Safely set data.table threads
  if (requireNamespace("data.table", quietly = TRUE)) {
    tryCatch(
      {
        data.table::setDTthreads(1)
      },
      error = function(e) {
        if (verbose) message("[geneSCOPE::createSCOPE xenium] Warning: Could not set data.table threads")
      }
    )
  }

  ## ---- 3. Package loading and configuration --------------------------------------------
  if (verbose) message("[geneSCOPE::createSCOPE xenium] Initializing data processing environment")

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
        if (verbose) message("[geneSCOPE::createSCOPE xenium]  Warning: Could not configure Arrow threading")
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
  if (verbose) message("[geneSCOPE::createSCOPE xenium] Loading cell centroids and transcripts")

  centroids_dt <- tryCatch(
    {
      # Use Arrow with limited threading for I/O stability
      if (requireNamespace("arrow", quietly = TRUE)) {
        arrow::set_cpu_count(ncores_io)
        arrow::set_io_thread_count(ncores_io)
      }

      cell_data <- arrow::read_parquet(cell_parq)
      cd <- as.data.table(cell_data)
      # CosMx: cell_id may repeat across FOVs; if a 'fov' column exists,
      # build a unique key 'cell_id_fov'.
      if ("fov" %in% names(cd)) {
        cd[, cell := paste0(as.character(cell_id), "_", as.character(fov))]
      } else {
        cd[, cell := as.character(cell_id)]
      }
      cd[, .(cell = cell, x = x_centroid, y = y_centroid)]
    },
    error = function(e) {
      stop("Error reading centroids from cells.parquet: ", e$message)
    }
  )

  if (nrow(centroids_dt) == 0) {
    stop("No centroids found in cells.parquet")
  }

  if (verbose) message("[geneSCOPE::createSCOPE xenium]  Found ", nrow(centroids_dt), " cells")

  ## ---- 6. Dataset preparation --------------------------------
  ds <- tryCatch(
    {
      arrow::open_dataset(tf_parq)
    },
    error = function(e) {
      stop("Error opening transcripts dataset: ", e$message)
    }
  )

  # Apply filters
  if (!is.null(filter_genes)) {
    if (verbose) message("[geneSCOPE::createSCOPE xenium]  Filtering genes: ", length(filter_genes), " genes")
    ds <- ds |> filter(feature_name %in% filter_genes)
  }

  if (!is.null(max_dist_mol_nuc)) {
    if (verbose) message("[geneSCOPE::createSCOPE xenium]  Filtering by nucleus distance: max ", max_dist_mol_nuc)
    ds <- ds |> filter(nucleus_distance <= max_dist_mol_nuc)
  }

  if (filtermolecule) {
    if (verbose) message("[geneSCOPE::createSCOPE xenium]  Filtering molecules by prefix exclusion")
    pat <- paste0("^(?:", paste(exclude_prefix, collapse = "|"), ")")
    ds <- ds |> filter(!grepl(pat, feature_name))
  }

  if (!is.null(filterqv)) {
    if (verbose) message("[geneSCOPE::createSCOPE xenium]  Filtering by QV score: min ", filterqv)
    ds <- ds |> filter(qv >= filterqv)
  }

  # Compute bounds
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
      "[geneSCOPE::createSCOPE xenium] Data bounds: x=[", round(bounds_global$xmin), ",",
      round(bounds_global$xmax), "], y=[", round(bounds_global$ymin),
      ",", round(bounds_global$ymax), "]"
    )
  }

  if (flip_y) {
    centroids_dt <- .flipCoordinates(centroids_dt, y_max)
    if (verbose) message("[geneSCOPE::createSCOPE xenium]  Applied Y-axis flip")
  }

  ## ---- 7. ROI polygon processing -------------------------
  user_poly <- NULL
  if (!is.null(coord_file)) {
    if (verbose) message("[geneSCOPE::createSCOPE xenium] Processing ROI polygon from file")

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

    if (verbose) message("[geneSCOPE::createSCOPE xenium]  Polygon has ", nrow(poly_tb), " vertices")

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
    if (verbose) message("[geneSCOPE::createSCOPE xenium] Creating automatic ROI from data boundaries")

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
      message(
        "[geneSCOPE::createSCOPE xenium] Automatic ROI created: x=[", round(min(poly_coords$x)), ",",
        round(max(poly_coords$x)), "], y=[", round(min(poly_coords$y)), ",",
        round(max(poly_coords$y)), "]"
      )
    }
  }

  ## ---- 8. Centroid clipping -----------------------------
  if (!is.null(user_poly) && nrow(centroids_dt)) {
    if (verbose) message("[geneSCOPE::createSCOPE xenium] Clipping centroids to ROI...")

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

    if (verbose) message("[geneSCOPE::createSCOPE xenium] Retained ", nrow(centroids_dt), " cells within ROI")
  }

  keep_cells <- if (nrow(centroids_dt)) unique(centroids_dt$cell) else NULL

  ## ---- 9. Initialize scope_object ----------------------------------------
  scope_obj <- new("scope_object",
    coord     = list(centroids = centroids_dt),
    grid      = list(),
    meta.data = data.frame()
  )

  ## ---- 10. Process segmentation data -------------------------------
  if (verbose) message("[geneSCOPE::createSCOPE xenium] Processing segmentation data...")

  for (tag in names(seg_files)) {
    if (verbose) message("[geneSCOPE::createSCOPE xenium] Processing ", tag, " segmentation...")

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

    scope_obj@coord[[paste0("segmentation_", tag)]] <- seg_res$points
    scope_obj@coord[[paste0("segmentation_polys_", tag)]] <- seg_res$polygons
  }

  ## ---- 11. Pre-fetch molecules (ROI case) ----------------------------
  mol_small <- NULL
  if (!is.null(user_poly)) {
    if (verbose) message("[geneSCOPE::createSCOPE xenium] Pre-fetching molecules for ROI...")

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

    if (verbose) message("[geneSCOPE::createSCOPE xenium] Retained ", nrow(mol_small), " molecules within ROI")
  }

  ## ---- 12. Batch scan counting (panoramic case, using original iterator approach) ----------------------
  counts_list <- setNames(vector("list", length(grid_length)), paste0("lg", grid_length))
  for (i in seq_along(counts_list)) {
    counts_list[[i]] <- data.table(grid_id = character(), gene = character(), count = integer())
  }

  if (is.null(user_poly)) {
    if (verbose) message("[geneSCOPE::createSCOPE xenium] Starting multi-resolution scan...")

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
          if (verbose) message("[geneSCOPE::createSCOPE xenium] Scanner finished or error: ", e$message)
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

      for (k in seq_along(grid_length)) {
        lg <- grid_length[k]
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
    if (verbose) message(sprintf("[geneSCOPE::createSCOPE xenium] Grid %.1f µm …", lg))
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
      seg_layers <- names(scope_obj@coord)[grepl("^segmentation_", names(scope_obj@coord)) &
        !grepl("_polys_", names(scope_obj@coord))]
      seg_dt_all <- rbindlist(scope_obj@coord[seg_layers], use.names = TRUE, fill = TRUE)
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

    scope_obj@grid[[paste0("grid", lg)]] <<- list(
      grid_info = grid_dt,
      counts = cnt,
      grid_length = lg,
      xbins_eff = length(x_breaks) - 1L,
      ybins_eff = length(y_breaks) - 1L
    )
  }

  if (verbose) message("[geneSCOPE::createSCOPE xenium] Building grid layers …")
  invisible(lapply(sort(unique(grid_length)), build_grid))
  if (verbose) message("[geneSCOPE::createSCOPE xenium] Grid construction finished.")

  ## ---- 14. Integrity check ----------------------------------------------
  if (length(scope_obj@grid) == 0 ||
    all(vapply(scope_obj@grid, function(g) nrow(g$grid_info) == 0L, logical(1)))) {
    warning("No effective grid layer generated – please check filters/parameters.")
  }

  scope_obj
}


#' @title Build scope_object from CosMx (AtomX flat files) with Sopa-like conventions
#' @description
#'   Normalizes CosMx AtomX flat files (flatFiles) into the parquet layout
#'   expected by \code{createSCOPE_xenium} (transcripts.parquet, cells.parquet,
#'   segmentation_boundaries.parquet, etc.), then delegates to
#'   \code{createSCOPE_xenium} to build a scope_object identical in structure
#'   to the Xenium pipeline.
#'
#'   The function first checks whether the target parquet files already exist;
#'   if not, it attempts to synthesize them from
#'   \code{flatFiles/*_{tx|polygons|fov_positions}_file.csv(.gz)}. Only
#'   \code{transcripts.parquet} is guaranteed to be produced; missing
#'   \code{cells.parquet}/segmentation parquet files are handled by
#'   \code{seg_type} fallback.
#'
#' @param cosmx_root Character. Root directory of a CosMx project; must contain
#'        \file{flatFiles/} and optionally \file{DecodedFiles/}.
#' @param grid_length Numeric vector. Grid edge length(s) (in the same
#'        coordinate unit/pixels) forwarded to \code{createSCOPE_xenium}.
#' @param seg_type Character. One of \code{"cell"}, \code{"nucleus"}, or
#'        \code{"both"}. Automatically downgraded if a segmentation parquet is
#'        missing.
#' @param dataset_id Character or NULL. If NULL, inferred from
#'        \file{flatFiles/*/*_tx_file.csv(.gz)} filename.
#' @param pixel_size_um Numeric. Pixel size (µm/px) used when unit conversion is
#'        required; default \code{0.120280945}.
#' @param verbose Logical. Print progress messages.
#' @param ... Additional arguments forwarded to \code{createSCOPE_xenium}
#'        (e.g. \code{filter_genes}, \code{filterqv}, \code{flip_y}).
#' @return A \code{scope_object} whose structure matches the output of
#'         \code{createSCOPE_xenium}.
#' @importFrom data.table fread as.data.table
#' @importFrom arrow write_parquet
#' @export
createSCOPE_cosmx <- function(cosmx_root,
                              grid_length,
                              seg_type = c("cell", "nucleus", "both"),
                              dataset_id = NULL,
                              pixel_size_um = 0.120280945,
                              verbose = TRUE,
                              build_from_flatfiles = TRUE,
                              build_cells_from_polygons = TRUE,
                              ...) {
  stopifnot(dir.exists(cosmx_root))
  seg_type <- match.arg(seg_type)

  if (verbose) message("[geneSCOPE::createSCOPE cosmx] Root = ", cosmx_root)

  # 1) Target parquet paths
  tf_parq   <- file.path(cosmx_root, "transcripts.parquet")
  cell_parq <- file.path(cosmx_root, "cells.parquet")
  seg_cell  <- file.path(cosmx_root, "cell_boundaries.parquet")
  seg_any   <- file.path(cosmx_root, "segmentation_boundaries.parquet")
  nuc_parq  <- file.path(cosmx_root, "nucleus_boundaries.parquet")

  # Working directory: if we need to synthesize any inputs, use a temp dir
  work_dir <- cosmx_root
  need_temp <- FALSE

  # 2) If transcripts.parquet is missing, synthesize from flatFiles
  if (!file.exists(tf_parq) && build_from_flatfiles) {
    if (verbose) message("[geneSCOPE::createSCOPE cosmx] transcripts.parquet not found; generating from flatFiles …")

    # Locate flatFiles and infer dataset ID if needed
    ff_dir <- file.path(cosmx_root, "flatFiles")
    if (!dir.exists(ff_dir)) {
      stop("flatFiles/ not found under ", cosmx_root)
    }
    # Search for *_tx_file.csv(.gz)
    tx_candidates <- c(
      Sys.glob(file.path(ff_dir, "*/*_tx_file.csv.gz")),
      Sys.glob(file.path(ff_dir, "*/*_tx_file.csv"))
    )
    if (length(tx_candidates) == 0L) {
      stop("No *_tx_file.csv(.gz) found under ", ff_dir)
    }
    tx_path <- tx_candidates[[1]]

    if (is.null(dataset_id)) {
      # Infer dataset ID from filename: External_xxx_tx_file.csv.gz → External_xxx
      base <- basename(tx_path)
      dataset_id <- sub("_tx_file\\.csv(\\.gz)?$", "", base)
    }

    # Read tx CSV. Expected minimum columns: fov, x_global_px, y_global_px, target
    suppressPackageStartupMessages({
      library(data.table)
      library(arrow)
    })
    if (verbose) message("[geneSCOPE::createSCOPE cosmx] Reading ", tx_path)
    tx <- tryCatch({ data.table::fread(tx_path) }, error = function(e) stop("Failed to read ", tx_path, ": ", e$message))
    cn <- tolower(names(tx)); names(tx) <- cn
    headerless <- all(grepl("^v[0-9]+$", cn))

    if (!headerless && all(c("x_global_px","y_global_px","target") %in% cn)) {
      xg <- tx[["x_global_px"]]; yg <- tx[["y_global_px"]]; tgt <- tx[["target"]]
      qv <- if ("qv" %in% cn) tx[["qv"]] else Inf
      nd <- if ("nucleus_distance" %in% cn) tx[["nucleus_distance"]] else 0
    } else {
      # headerless (positional mapping) or missing expected names
      if (ncol(tx) < 9L) stop("headerless tx_file has insufficient columns (<9)")
      if (verbose) message("[geneSCOPE::createSCOPE cosmx] tx_file headerless – using positional mapping col6/col7/col9 (and col8 as nucleus_distance)")
      xg <- tx[[6]]; yg <- tx[[7]]; tgt <- tx[[9]]
      qv <- if (ncol(tx) >= 10 && is.numeric(tx[[10]])) tx[[10]] else Inf
      nd <- if (ncol(tx) >= 8 && is.numeric(tx[[8]])) tx[[8]] else 0
    }
    tx_out <- data.table::data.table(
      x_location = as.numeric(xg),
      y_location = as.numeric(yg),
      feature_name = as.character(tgt),
      qv = as.numeric(qv),
      nucleus_distance = as.numeric(nd)
    )

    # choose working directory
    work_dir <- file.path(tempdir(), paste0("cosmx_scope_", as.integer(Sys.time())))
    if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
    need_temp <- TRUE
    tf_parq <- file.path(work_dir, "transcripts.parquet")
    if (verbose) message("[geneSCOPE::createSCOPE cosmx] Writing ", tf_parq)
    arrow::write_parquet(tx_out, tf_parq)
  } else if (verbose) {
    message("[geneSCOPE::createSCOPE cosmx] Using existing transcripts.parquet")
  }

  # 3) Segmentation and centroids: downgrade seg_type if segmentation is missing
  have_cell_seg <- file.exists(seg_cell) || file.exists(seg_any)
  have_nuc_seg  <- file.exists(nuc_parq)

  if (seg_type %in% c("cell","both") && !have_cell_seg) {
    if (verbose) message("[geneSCOPE::createSCOPE cosmx] Cell segmentation parquet not found; downgrade seg_type")
    seg_type <- if (have_nuc_seg) "nucleus" else "cell" # Prefer nucleus; otherwise pass 'cell' and let downstream handle absence
  }
  if (seg_type %in% c("nucleus","both") && !have_nuc_seg) {
    if (verbose) message("[geneSCOPE::createSCOPE cosmx] Nucleus segmentation parquet not found; downgrade seg_type")
    seg_type <- if (have_cell_seg) "cell" else "nucleus"
  }

  # 3.1) If cells.parquet is missing and polygon derivation is allowed, generate a temporary cells.parquet
  if (!file.exists(cell_parq) && build_from_flatfiles && build_cells_from_polygons) {
    # find polygons in flatFiles
    poly_path <- NULL
    ff_dir <- file.path(cosmx_root, "flatFiles")
    cand <- c(Sys.glob(file.path(ff_dir, "*/*-polygons.csv.gz")), Sys.glob(file.path(ff_dir, "*/*_polygons.csv.gz")))
    if (length(cand)) poly_path <- cand[[1]]
    if (!is.null(poly_path)) {
      if (verbose) message("[geneSCOPE::createSCOPE cosmx] cells.parquet not found; deriving centroids from polygons: ", basename(poly_path))
      pol <- data.table::fread(poly_path)
      cn <- tolower(names(pol)); names(pol) <- cn
      headerless <- all(grepl("^v[0-9]+$", cn))
      if (headerless) {
        # Map: col1=fov, col2=cell_id, col6=x_global_px, col7=y_global_px
        if (ncol(pol) < 7L) stop("headerless polygons has insufficient columns (<7)")
        fov <- pol[[1]]; cid <- pol[[2]]; xg <- pol[[6]]; yg <- pol[[7]]
      } else {
        # try common lowercase names
        fx <- if ("fov" %in% cn) pol[["fov"]] else pol[[1]]
        cx <- if ("cellid" %in% cn) pol[["cellid"]] else if ("cell_id" %in% cn) pol[["cell_id"]] else pol[[2]]
        xg <- if ("x_global_px" %in% cn) pol[["x_global_px"]] else pol[[which.max(cn=="x")]]
        yg <- if ("y_global_px" %in% cn) pol[["y_global_px"]] else pol[[which.max(cn=="y")]]
        fov <- fx; cid <- cx
      }
      dt <- data.table::data.table(fov=as.integer(fov), cell_id=as.character(cid), x=as.numeric(xg), y=as.numeric(yg))
      cent <- dt[, .(x_centroid = mean(x, na.rm=TRUE), y_centroid = mean(y, na.rm=TRUE)), by=.(cell_id, fov)]
      # write to work_dir
      if (!need_temp) { work_dir <- file.path(tempdir(), paste0("cosmx_scope_", as.integer(Sys.time()))); dir.create(work_dir, recursive = TRUE); need_temp <- TRUE }
      cell_parq <- file.path(work_dir, "cells.parquet")
      if (verbose) message("[geneSCOPE::createSCOPE cosmx] Writing ", cell_parq)
      arrow::write_parquet(cent, cell_parq)
    } else if (verbose) {
      message("[geneSCOPE::createSCOPE cosmx] polygons not found; cannot derive cells.parquet")
    }
  }

  # 3.2) If using a temporary working directory, symlink/copy existing segmentation parquets into it so downstream code sees a single root
  if (need_temp) {
    link_or_copy <- function(src, dst) {
      if (file.exists(src) && !file.exists(dst)) {
        ok <- FALSE
        suppressWarnings({ ok <- file.symlink(src, dst) })
        if (!ok) invisible(file.copy(src, dst, overwrite = FALSE))
      }
    }
    link_or_copy(file.path(cosmx_root, "segmentation_boundaries.parquet"), file.path(work_dir, "segmentation_boundaries.parquet"))
    link_or_copy(file.path(cosmx_root, "cell_boundaries.parquet"),         file.path(work_dir, "cell_boundaries.parquet"))
    link_or_copy(file.path(cosmx_root, "nucleus_boundaries.parquet"),      file.path(work_dir, "nucleus_boundaries.parquet"))
  }

  # 4) Delegate to createSCOPE_xenium (flip_y defaults to TRUE to match Xenium coordinates)
  if (verbose) message("[geneSCOPE::createSCOPE cosmx] Delegating to createSCOPE_xenium() …")
  scope_obj <- createSCOPE_xenium(
    xenium_dir = work_dir,
    grid_length = grid_length,
    seg_type = seg_type,
    verbose = verbose,
    ...
  )

  # Attach platform/dataset metadata
  if (nrow(scope_obj@meta.data) > 0L) {
    scope_obj@meta.data$platform <- rep("CosMx", nrow(scope_obj@meta.data))
    scope_obj@meta.data$dataset  <- rep(if (is.null(dataset_id)) NA_character_ else dataset_id,
                                        nrow(scope_obj@meta.data))
  } else {
    # When meta.data is empty, record identifiers in @stats to avoid data.frame length mismatch
    scope_obj@stats$platform <- "CosMx"
    scope_obj@stats$dataset  <- if (is.null(dataset_id)) NA_character_ else dataset_id
  }
  scope_obj
}

#' @title Build scope_object from 10x Visium spot-level outputs
#' @description
#'   Reads standard 10x Visium outputs (matrix + spatial metadata) and constructs
#'   a \code{scope_object} where the sole grid layer matches the Visium spot size.
#'   Coordinates are expressed in microns using \code{scalefactors_json.json},
#'   and only in-tissue spots are retained by default.
#' @param visium_dir Character. Path to the Visium run (typically the \code{outs/} directory)
#'        containing \code{spatial/} and the feature-barcode matrix folder.
#' @param use_filtered Logical. If \code{TRUE} (default), prefer
#'        \code{filtered_feature_bc_matrix}; otherwise fall back to raw counts.
#' @param include_in_tissue Logical. Keep only spots marked with
#'        \code{in_tissue == 1} (default \code{TRUE}).
#' @param grid_multiplier Numeric scalar. Multiplier applied to the inferred spot
#'        diameter when defining the grid width/height (default 1, i.e. spot-sized grid).
#' @param flip_y Logical. Whether to flip Y so that the origin is bottom-left
#'        (default \code{FALSE}; Visium images are usually top-left origin).
#' @param verbose Logical. Print progress messages (default \code{TRUE}).
#' @return A \code{scope_object} containing spot centroids, a single grid layer,
#'         and a spot-level count matrix stored in \code{@cells$counts}.
#' @details
#'   This helper bypasses \code{\link{createSCOPE}} because Visium provides
#'   spot-level counts instead of molecule coordinates. Grid tiles are centered
#'   on Visium spot centroids and their edge length equals the spot diameter
#'   reported in \code{scalefactors_json.json} multiplied by \code{grid_multiplier}.
#' @export
createSCOPE_visium <- function(visium_dir,
                               use_filtered = TRUE,
                               include_in_tissue = TRUE,
                               grid_multiplier = 1,
                               flip_y = FALSE,
                               verbose = TRUE) {
  stopifnot(dir.exists(visium_dir))

  if (!is.numeric(grid_multiplier) || length(grid_multiplier) != 1L || !is.finite(grid_multiplier)) {
    stop("grid_multiplier must be a single finite numeric value.")
  }
  if (grid_multiplier <= 0) {
    stop("grid_multiplier must be positive.")
  }

  suppressPackageStartupMessages({
    library(Matrix)
    library(data.table)
  })

  candidate_dirs <- unique(file.path(visium_dir, c(
    if (isTRUE(use_filtered)) "filtered_feature_bc_matrix",
    if (!isTRUE(use_filtered)) "raw_feature_bc_matrix",
    "filtered_feature_bc_matrix",
    "raw_feature_bc_matrix"
  )))
  matrix_dir <- candidate_dirs[dir.exists(candidate_dirs)][1]

  if (is.na(matrix_dir)) {
    stop("Could not find filtered_feature_bc_matrix/ or raw_feature_bc_matrix/ under ", visium_dir)
  }
  if (verbose) message("[geneSCOPE::createSCOPE visium] Using matrix directory: ", basename(matrix_dir))

  matrix_file <- file.path(matrix_dir, "matrix.mtx.gz")
  if (!file.exists(matrix_file)) matrix_file <- file.path(matrix_dir, "matrix.mtx")
  feature_file <- file.path(matrix_dir, "features.tsv.gz")
  if (!file.exists(feature_file)) feature_file <- file.path(matrix_dir, "features.tsv")
  barcode_file <- file.path(matrix_dir, "barcodes.tsv.gz")
  if (!file.exists(barcode_file)) barcode_file <- file.path(matrix_dir, "barcodes.tsv")

  if (!file.exists(matrix_file)) {
    stop("matrix.mtx(.gz) not found in ", matrix_dir)
  }
  if (!file.exists(feature_file) || !file.exists(barcode_file)) {
    stop("Required features.tsv(.gz) or barcodes.tsv(.gz) missing in ", matrix_dir)
  }

  counts_mt <- Matrix::readMM(if (grepl("\\.gz$", matrix_file, ignore.case = TRUE)) gzfile(matrix_file) else matrix_file)
  # Use CsparseMatrix to avoid Matrix::DeprecatedCoerce warnings
  counts_mat <- as(counts_mt, "CsparseMatrix")
  rm(counts_mt)

  features_dt <- data.table::fread(feature_file, header = FALSE)
  barcodes_dt <- data.table::fread(barcode_file, header = FALSE)
  if (ncol(features_dt) < 2) {
    stop("features.tsv must contain at least two columns (gene id, gene name).")
  }
  gene_ids <- as.character(features_dt[[1]])
  gene_names <- make.unique(as.character(features_dt[[2]]))
  feature_type <- if (ncol(features_dt) >= 3) as.character(features_dt[[3]]) else rep(NA_character_, length(gene_names))
  barcodes <- as.character(barcodes_dt[[1]])
  dimnames(counts_mat) <- list(gene_names, barcodes)

  gene_meta <- data.frame(
    feature_id = gene_ids,
    feature_type = feature_type,
    stringsAsFactors = FALSE,
    row.names = gene_names
  )

  spatial_dir <- file.path(visium_dir, "spatial")
  if (!dir.exists(spatial_dir)) {
    stop("spatial/ directory not found under ", visium_dir)
  }

  pos_candidates <- c(
    file.path(spatial_dir, "tissue_positions_list.csv"),
    file.path(spatial_dir, "tissue_positions.csv")
  )
  pos_path <- pos_candidates[file.exists(pos_candidates)][1]
  if (is.na(pos_path)) {
    stop("Could not find tissue_positions_list.csv or tissue_positions.csv under ", spatial_dir)
  }
  pos_raw <- data.table::fread(pos_path)
  if (verbose) message("[geneSCOPE::createSCOPE visium] Loaded ", nrow(pos_raw), " spot positions")

  standardize_positions <- function(dt) {
    if (!data.table::is.data.table(dt)) dt <- data.table::as.data.table(dt)
    if (ncol(dt) < 6) {
      stop("Visium position file must contain at least six columns.")
    }
    if (all(grepl("^V[0-9]+$", names(dt)))) {
      cols <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres")
      data.table::setnames(dt, cols[seq_len(ncol(dt))])
    } else {
      lower <- tolower(names(dt))
      rename_first <- function(target, aliases) {
        idx <- which(lower %in% aliases)
        if (length(idx)) {
          data.table::setnames(dt, idx[1], target)
          lower[idx[1]] <<- target
        }
      }
      rename_first("barcode", c("barcode", "barcodes"))
      rename_first("in_tissue", c("in_tissue", "intissue"))
      rename_first("array_row", c("array_row", "row", "spot_row", "arrayrow"))
      rename_first("array_col", c("array_col", "col", "spot_col", "arraycol"))
      rename_first("pxl_col_in_fullres", c("pxl_col_in_fullres", "pixel_x", "pxl_col", "imagecol", "col_coord"))
      rename_first("pxl_row_in_fullres", c("pxl_row_in_fullres", "pixel_y", "pxl_row", "imagerow", "row_coord"))
    }
    required <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres")
    missing_cols <- setdiff(required, names(dt))
    if (length(missing_cols)) {
      stop("Position file missing required columns: ", paste(missing_cols, collapse = ", "))
    }
    dt[, `:=`(
      barcode = as.character(barcode),
      in_tissue = as.integer(in_tissue),
      array_row = as.integer(array_row),
      array_col = as.integer(array_col),
      pxl_col_in_fullres = as.numeric(pxl_col_in_fullres),
      pxl_row_in_fullres = as.numeric(pxl_row_in_fullres)
    )]
    dt
  }

  pos_dt <- standardize_positions(pos_raw)
  if (isTRUE(include_in_tissue)) {
    pos_dt <- pos_dt[in_tissue == 1L]
    if (verbose) message("[geneSCOPE::createSCOPE visium] Retained ", nrow(pos_dt), " in-tissue spots")
  }
  if (nrow(pos_dt) == 0L) {
    stop("No spots remain after applying the in_tissue filter.")
  }

  col_idx <- match(pos_dt$barcode, colnames(counts_mat))
  if (anyNA(col_idx)) {
    missing <- pos_dt$barcode[is.na(col_idx)]
    warning("Removing ", length(missing), " spots absent from the count matrix.")
    pos_dt <- pos_dt[!is.na(col_idx)]
    col_idx <- col_idx[!is.na(col_idx)]
  }
  if (!length(col_idx)) {
    stop("No overlap between Visium positions and the count matrix barcodes.")
  }

  counts_mat <- counts_mat[, col_idx, drop = FALSE]
  colnames(counts_mat) <- pos_dt$barcode

  nonzero_genes <- Matrix::rowSums(counts_mat) > 0
  if (!any(nonzero_genes)) {
    stop("All genes have zero counts across the selected spots.")
  }
  if (any(!nonzero_genes)) {
    counts_mat <- counts_mat[nonzero_genes, , drop = FALSE]
  }
  gene_meta <- gene_meta[rownames(counts_mat), , drop = FALSE]

  scalefactor_candidates <- c(
    file.path(spatial_dir, "scalefactors_json.json"),
    file.path(spatial_dir, "scalefactors_json_fullres.json")
  )
  scalefactor_path <- scalefactor_candidates[file.exists(scalefactor_candidates)][1]
  if (is.na(scalefactor_path)) {
    stop("scalefactors_json.json not found under ", spatial_dir)
  }
  scalefactors <- jsonlite::fromJSON(scalefactor_path)
  microns_per_pixel <- scalefactors$microns_per_pixel
  spot_diameter_um <- scalefactors$spot_diameter_microns
  spot_diameter_px <- scalefactors$spot_diameter_fullres

  # Derive missing pieces when possible
  if (is.null(spot_diameter_um) && !is.null(spot_diameter_px) && !is.null(microns_per_pixel)) {
    spot_diameter_um <- as.numeric(spot_diameter_px) * as.numeric(microns_per_pixel)
  }
  if (is.null(microns_per_pixel) && !is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
    microns_per_pixel <- as.numeric(spot_diameter_um) / as.numeric(spot_diameter_px)
  }

  # Fallback for outputs lacking microns_per_pixel and spot_diameter_microns
  if (is.null(microns_per_pixel) && is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
    default_spot_um <- getOption("geneSCOPE.default_visium_spot_um", 55)
    microns_per_pixel <- as.numeric(default_spot_um) / as.numeric(spot_diameter_px)
    spot_diameter_um  <- as.numeric(default_spot_um)
    if (isTRUE(verbose)) {
      message(
        "[geneSCOPE::createSCOPE visium] microns_per_pixel missing; assumed spot_diameter_microns = ",
        default_spot_um, " um and computed microns_per_pixel = ",
        format(round(microns_per_pixel, 9), scientific = FALSE)
      )
    }
  }

  if (is.null(microns_per_pixel)) {
    stop("scalefactors file missing microns_per_pixel; cannot convert coordinates.")
  }
  if (is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
    spot_diameter_um <- as.numeric(spot_diameter_px) * as.numeric(microns_per_pixel)
  }
  if (is.null(spot_diameter_um)) {
    stop("scalefactors file missing spot diameter information.")
  }

  microns_per_pixel <- as.numeric(microns_per_pixel)
  base_spot_um <- as.numeric(spot_diameter_um)
  grid_size_um <- base_spot_um * as.numeric(grid_multiplier)

  pos_dt[, `:=`(
    x_um = pxl_col_in_fullres * microns_per_pixel,
    y_um = pxl_row_in_fullres * microns_per_pixel
  )]

  if (isTRUE(flip_y)) {
    y_max <- max(pos_dt$y_um, na.rm = TRUE)
    pos_dt[, y_um := y_max - y_um]
  }

  centroids_dt <- data.table::data.table(
    cell = pos_dt$barcode,
    x = pos_dt$x_um,
    y = pos_dt$y_um,
    array_row = pos_dt$array_row,
    array_col = pos_dt$array_col,
    in_tissue = pos_dt$in_tissue
  )

  unique_cols <- sort(unique(pos_dt$array_col))
  unique_rows <- sort(unique(pos_dt$array_row))
  col_map <- setNames(seq_along(unique_cols), unique_cols)
  row_map <- setNames(seq_along(unique_rows), unique_rows)

  grid_dt <- data.table::data.table(
    grid_id = pos_dt$barcode,
    gx = as.integer(col_map[as.character(pos_dt$array_col)]),
    gy = as.integer(row_map[as.character(pos_dt$array_row)]),
    xmin = pos_dt$x_um - grid_size_um / 2,
    xmax = pos_dt$x_um + grid_size_um / 2,
    ymin = pos_dt$y_um - grid_size_um / 2,
    ymax = pos_dt$y_um + grid_size_um / 2,
    center_x = pos_dt$x_um,
    center_y = pos_dt$y_um,
    width = rep(grid_size_um, nrow(pos_dt)),
    height = rep(grid_size_um, nrow(pos_dt)),
    idx = seq_len(nrow(pos_dt)),
    array_row = pos_dt$array_row,
    array_col = pos_dt$array_col
  )

  xbins_eff <- length(unique(grid_dt$gx))
  ybins_eff <- length(unique(grid_dt$gy))

  counts_tbl <- Matrix::summary(counts_mat)
  counts_dt <- data.table::data.table(
    grid_id = colnames(counts_mat)[counts_tbl$j],
    gene = rownames(counts_mat)[counts_tbl$i],
    count = as.numeric(counts_tbl$x)
  )
  counts_dt <- counts_dt[count > 0]

  gene_meta$platform <- "Visium"

  fmt_len <- formatC(grid_size_um, format = "fg", digits = 6)
  fmt_len <- trimws(fmt_len)
  fmt_len <- sub("\\.?0+$", "", fmt_len)
  grid_name <- paste0("grid", fmt_len)

  scope_obj <- new("scope_object",
    coord = list(centroids = centroids_dt),
    grid = list(),
    meta.data = gene_meta,
    cells = list(counts = counts_mat),
    stats = list(),
    density = list()
  )

  ## --- infer hex orientation (flat-top vs pointy-top) -----------------
  infer_hex_orientation <- function(cent, tol_deg = 15) {
    cn <- c("x","y","array_row","array_col")
    if (!all(cn %in% names(cent))) return(list(orientation = NA_character_, horiz_frac = NA_real_, vert_frac = NA_real_, tol_deg = tol_deg))
    dt <- data.table::as.data.table(cent)[, .(x, y, array_row = as.integer(array_row), array_col = as.integer(array_col))]
    # helper to join neighbor at offset (dr, dc)
    get_pairs <- function(dr, dc) {
      a <- dt
      b <- data.table::copy(dt)
      b[, `:=`(array_row_nb = array_row - dr, array_col_nb = array_col - dc, x_nb = x, y_nb = y)]
      m <- merge(a, b, by.x = c("array_row", "array_col"), by.y = c("array_row_nb", "array_col_nb"), all = FALSE)
      if (!nrow(m)) return(NULL)
      m[, .(dx = x_nb - x, dy = y_nb - y)]
    }
    edges <- data.table::rbindlist(list(
      get_pairs(0, 1),   # E
      get_pairs(1, 0),   # S (vertical-ish in pointy-top)
      get_pairs(1, 1),   # SE/NW
      get_pairs(1, -1)   # SW/NE
    ), use.names = TRUE, fill = TRUE)
    if (is.null(edges) || !nrow(edges)) return(list(orientation = NA_character_, horiz_frac = NA_real_, vert_frac = NA_real_, tol_deg = tol_deg))
    theta <- abs(atan2(edges$dy, edges$dx) * 180 / pi)
    horiz <- pmin(theta, 180 - theta) <= tol_deg
    vert  <- abs(theta - 90) <= tol_deg
    hfrac <- mean(horiz, na.rm = TRUE)
    vfrac <- mean(vert,  na.rm = TRUE)
    list(orientation = ifelse(hfrac >= vfrac, "flat-top", "pointy-top"), horiz_frac = hfrac, vert_frac = vfrac, tol_deg = tol_deg)
  }
  ori <- infer_hex_orientation(centroids_dt, tol_deg = 15)

  scope_obj@grid[[grid_name]] <- list(
    grid_info = grid_dt,
    counts = counts_dt,
    grid_length = grid_size_um,
    xbins_eff = xbins_eff,
    ybins_eff = ybins_eff,
    spot_diameter_um = base_spot_um,
    microns_per_pixel = microns_per_pixel,
    visium_orientation = ori$orientation,
    visium_orientation_horiz_frac = ori$horiz_frac,
    visium_orientation_vert_frac  = ori$vert_frac,
    visium_orientation_tol_deg    = ori$tol_deg
  )

  if (verbose) {
    message(
      "[geneSCOPE::createSCOPE visium] scope_object created with ",
      nrow(centroids_dt), " spots and grid_length = ", round(grid_size_um, 3), " µm"
    )
  }

  scope_obj
}

#' @title Auto-detect data type and create scope_object
#' @description
#'   Wrapper that inspects a data directory and dispatches to one of:
#'   \code{createSCOPE_xenium()}, \code{createSCOPE_cosmx()}, or
#'   \code{.seurat_visium_to_scope()} based on characteristic files.
#'   For all platforms (Xenium/CosMx/Visium), a \code{scope_object} is returned.
#' @param data_dir Character. Directory containing a Xenium, CosMx, or Visium run.
#'        As a convenience, you may also pass \code{xenium_dir}, \code{cosmx_root},
#'        or \code{visium_dir} via \code{...}; the first non-NULL will be used.
#' @param prefer Character. One of \code{"auto"}, \code{"xenium"}, \code{"cosmx"},
#'        \code{"visium"}. When \code{"auto"}, detection order is Xenium → Visium → CosMx.
#' @param return_both Logical. For Visium only, return list(scope_obj, seurat) when TRUE.
#' @param ... Additional parameters forwarded to the chosen builder. Note that
#'        Xenium/CosMx require \code{grid_length}; Visium does not.
#' @return A scope_object (or a list for Visium when \code{return_both=TRUE}).
#' @export
createSCOPE <- function(data_dir = NULL,
                        prefer = c("auto", "xenium", "cosmx", "visium"),
                        verbose = TRUE,
                        sctransform = TRUE,
                        ...) {
  dots <- list(...)
  if (is.null(dots$verbose)) dots$verbose <- verbose
  # Backward-compatible aliasing
  if (is.null(data_dir)) {
    if (!is.null(dots$xenium_dir)) data_dir <- dots$xenium_dir
    else if (!is.null(dots$visium_dir)) data_dir <- dots$visium_dir
    else if (!is.null(dots$cosmx_root)) data_dir <- dots$cosmx_root
  }
  if (is.null(data_dir)) stop("Please provide data_dir (or xenium_dir/visium_dir/cosmx_root).")
  if (!dir.exists(data_dir)) stop("Directory not found: ", data_dir)

  prefer <- match.arg(prefer)

  is_xenium <- function(d) {
    file.exists(file.path(d, "transcripts.parquet")) ||
      file.exists(file.path(d, "cells.parquet"))
  }
  is_visium <- function(d) {
    sp <- file.path(d, "spatial")
    dir.exists(sp) && (
      file.exists(file.path(d, "filtered_feature_bc_matrix.h5")) ||
      dir.exists(file.path(d, "filtered_feature_bc_matrix")) ||
      dir.exists(file.path(d, "raw_feature_bc_matrix")) ||
      file.exists(file.path(sp, "scalefactors_json.json"))
    )
  }
  is_cosmx <- function(d) {
    ff <- file.path(d, "flatFiles")
    dec <- file.path(d, "DecodedFiles")
    dir.exists(ff) || dir.exists(dec)
  }

  pick <- switch(prefer,
    xenium = "xenium",
    visium = "visium",
    cosmx  = "cosmx",
    auto   = {
      if (is_xenium(data_dir)) "xenium"
      else if (is_visium(data_dir)) "visium"
      else if (is_cosmx(data_dir)) "cosmx"
      else stop("Unable to detect data type under ", data_dir,
                 "; expected Xenium (parquet), Visium (spatial/ + matrix), or CosMx (flatFiles/).")
    }
  )

  # Clean conflicting aliases in dots
  dots$xenium_dir <- NULL; dots$visium_dir <- NULL; dots$cosmx_root <- NULL

  if (pick == "xenium") {
    if (isTRUE(verbose)) message("[geneSCOPE::createSCOPE xenium] Dispatching to createSCOPE_xenium() …")
    if (is.null(dots$grid_length)) stop("grid_length is required for Xenium. Pass grid_length= ...")
    dots$xenium_dir <- data_dir
    return(do.call(createSCOPE_xenium, dots))
  }
  if (pick == "cosmx") {
    if (isTRUE(verbose)) message("[geneSCOPE::createSCOPE cosmx] Dispatching to createSCOPE_cosmx() …")
    if (is.null(dots$grid_length)) stop("grid_length is required for CosMx. Pass grid_length= ...")
    dots$cosmx_root <- data_dir
    return(do.call(createSCOPE_cosmx, dots))
  }

  # Visium → choose Seurat-first (SCT) or direct builder based on sctransform
  dots$visium_dir <- data_dir
  if (isTRUE(sctransform)) {
    if (isTRUE(verbose)) message("[geneSCOPE::createSCOPE visium] Dispatching to Seurat pipeline (SCTransform=TRUE) …")
    # Ensure the downstream function gets the flag
    if (is.null(dots$run_sctransform)) dots$run_sctransform <- TRUE
    res <- do.call(.seurat_visium_to_scope, dots)
    res
  } else {
    if (isTRUE(verbose)) message("[geneSCOPE::createSCOPE visium] Dispatching to createSCOPE_visium (SCTransform=FALSE) …")
    obj <- do.call(createSCOPE_visium, dots)
    obj
  }
}
#' @title Seurat-first Visium pipeline to scope_object
#' @description
#' Read a 10x Visium run directory using Seurat, optionally run SCTransform,
#' then construct and return a \code{scope_object}. Raw counts are kept in
#' \code{@cells$counts}; SCTransform residuals are written to
#' \code{@cells$SCT} if computed.
#'
#' This path is useful if you prefer to rely on Seurat's IO and preprocessing
#' steps first, then switch to geneSCOPE for grid-level workflows.
#'
#' @param visium_dir Path to the Visium run (the Space Ranger outs directory
#'   that contains \code{spatial/} and a feature-barcode matrix).
#' @param run_sctransform Logical; run SCTransform on the Seurat object
#'   (default TRUE).
#' @param sct_return_only_var_genes Logical; pass to SCTransform (default FALSE
#'   to keep all genes in scale.data).
#' @param vars_to_regress Optional character vector passed to SCTransform.
#' @param glmGamPoi Logical; use glmGamPoi backend if available (default TRUE).
#' @param include_in_tissue Keep only spots with \code{in_tissue == 1}
#'   (default TRUE).
#' @param grid_multiplier Numeric; multiplier on spot diameter for grid length
#'   (default 1).
#' @param flip_y Logical; flip Y coordinates to bottom-left origin (default FALSE).
#' @param verbose Logical; print progress (default TRUE).
#'
#' @return A \code{scope_object} with:
#'   - \code{@cells$counts}: raw counts (dgCMatrix)
#'   - \code{@cells$SCT}: SCTransform scale.data (if run)
#'   - one grid layer sized to the Visium spot diameter * \code{grid_multiplier}
#' @noRd
.seurat_visium_to_scope <- function(visium_dir,
                                   run_sctransform = TRUE,
                                   sct_return_only_var_genes = FALSE,
                                   vars_to_regress = NULL,
                                   glmGamPoi = TRUE,
                                   include_in_tissue = TRUE,
                                   grid_multiplier = 1,
                                   flip_y = FALSE,
                                   verbose = TRUE) {
  stopifnot(dir.exists(visium_dir))
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for .seurat_visium_to_scope(). Please install it.")
  }

  if (!is.numeric(grid_multiplier) || length(grid_multiplier) != 1L || !is.finite(grid_multiplier) || grid_multiplier <= 0) {
    stop("grid_multiplier must be a single positive finite numeric value.")
  }

  suppressPackageStartupMessages({
    library(Matrix)
    library(data.table)
  })

  # 1) Read with Seurat -----------------------------------------------------
  if (isTRUE(verbose)) message("[geneSCOPE::createSCOPE visium] Calling Seurat::Load10X_Spatial() …")
  # Prefer H5 if present, otherwise let Seurat handle defaults
  h5 <- file.path(visium_dir, "filtered_feature_bc_matrix.h5")
  has_h5 <- file.exists(h5)
  seu <- tryCatch({
    if (has_h5) {
      Seurat::Load10X_Spatial(data.dir = visium_dir, filename = basename(h5))
    } else {
      Seurat::Load10X_Spatial(data.dir = visium_dir)
    }
  }, error = function(e) {
    stop("Load10X_Spatial failed: ", conditionMessage(e))
  })

  # 2) Optional SCTransform -------------------------------------------------
  if (isTRUE(run_sctransform)) {
    if (isTRUE(verbose)) message("[geneSCOPE::createSCOPE visium] Calling Seurat::SCTransform() …")
    use_ggp <- isTRUE(glmGamPoi) && requireNamespace("glmGamPoi", quietly = TRUE)
    seu <- Seurat::SCTransform(
      object = seu,
      assay = "Spatial",
      new.assay.name = "SCT",
      return.only.var.genes = sct_return_only_var_genes,
      vars.to.regress = vars_to_regress,
      method = if (use_ggp) "glmGamPoi" else "poisson",
      verbose = isTRUE(verbose)
    )
    Seurat::DefaultAssay(seu) <- "SCT"
  }

  # 3) Extract matrices -----------------------------------------------------
  # Use layer=.. if available (SeuratObject >= 5) to avoid deprecation warnings; fallback to slot=..
  getAssayLayer <- function(obj, assay, which) {
    ga <- Seurat::GetAssayData
    fm <- try(formals(ga), silent = TRUE)
    if (!inherits(fm, "try-error") && "layer" %in% names(fm)) {
      ga(object = obj, assay = assay, layer = which)
    } else {
      ga(object = obj, assay = assay, slot = which)
    }
  }
  raw_counts <- getAssayLayer(seu, assay = "Spatial", which = "counts")
  if (!inherits(raw_counts, "dgCMatrix")) raw_counts <- as(raw_counts, "dgCMatrix")

  sct_mat <- NULL
  if (isTRUE(run_sctransform)) {
    # Pearson residuals (scale.data) by convention
    sct_mat <- getAssayLayer(seu, assay = "SCT", which = "scale.data")
  }

  # 4) Positions & scalefactors --------------------------------------------
  spatial_dir <- file.path(visium_dir, "spatial")
  if (!dir.exists(spatial_dir)) stop("spatial/ directory not found under ", visium_dir)

  pos_candidates <- c(file.path(spatial_dir, "tissue_positions_list.csv"),
                      file.path(spatial_dir, "tissue_positions.csv"))
  pos_path <- pos_candidates[file.exists(pos_candidates)][1]
  if (is.na(pos_path)) stop("Could not find tissue_positions_list.csv or tissue_positions.csv under ", spatial_dir)

  pos_raw <- data.table::fread(pos_path)
  # Standardize columns similar to createSCOPE_visium
  standardize_positions <- function(dt) {
    if (!data.table::is.data.table(dt)) dt <- data.table::as.data.table(dt)
    if (ncol(dt) < 6) stop("Visium position file must contain at least six columns.")
    if (all(grepl("^V[0-9]+$", names(dt)))) {
      cols <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres")
      data.table::setnames(dt, cols[seq_len(ncol(dt))])
    } else {
      lower <- tolower(names(dt))
      rename_first <- function(target, aliases) {
        idx <- which(lower %in% aliases)
        if (length(idx)) {
          data.table::setnames(dt, idx[1], target)
          lower[idx[1]] <<- target
        }
      }
      rename_first("barcode", c("barcode", "barcodes"))
      rename_first("in_tissue", c("in_tissue", "intissue"))
      rename_first("array_row", c("array_row", "row", "spot_row", "arrayrow"))
      rename_first("array_col", c("array_col", "col", "spot_col", "arraycol"))
      rename_first("pxl_col_in_fullres", c("pxl_col_in_fullres", "pixel_x", "pxl_col", "imagecol", "col_coord"))
      rename_first("pxl_row_in_fullres", c("pxl_row_in_fullres", "pixel_y", "pxl_row", "imagerow", "row_coord"))
    }
    required <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres")
    missing_cols <- setdiff(required, names(dt))
    if (length(missing_cols)) stop("Position file missing required columns: ", paste(missing_cols, collapse = ", "))
    dt[, `:=`(
      barcode = as.character(barcode),
      in_tissue = as.integer(in_tissue),
      array_row = as.integer(array_row),
      array_col = as.integer(array_col),
      pxl_col_in_fullres = as.numeric(pxl_col_in_fullres),
      pxl_row_in_fullres = as.numeric(pxl_row_in_fullres)
    )]
    dt
  }
  pos_dt <- standardize_positions(pos_raw)
  if (isTRUE(include_in_tissue)) pos_dt <- pos_dt[in_tissue == 1L]

  # Align with columns in raw_counts (barcodes in Seurat object)
  keep <- pos_dt$barcode %in% colnames(raw_counts)
  if (!all(keep)) pos_dt <- pos_dt[keep]
  if (nrow(pos_dt) == 0L) stop("No overlap between positions and Seurat barcodes.")

  scalefactor_candidates <- c(file.path(spatial_dir, "scalefactors_json.json"),
                              file.path(spatial_dir, "scalefactors_json_fullres.json"))
  sf_path <- scalefactor_candidates[file.exists(scalefactor_candidates)][1]
  if (is.na(sf_path)) stop("scalefactors_json.json not found under ", spatial_dir)
  sf <- jsonlite::fromJSON(sf_path)
  microns_per_pixel <- sf$microns_per_pixel
  spot_diameter_um <- sf$spot_diameter_microns
  spot_diameter_px <- sf$spot_diameter_fullres
  if (is.null(spot_diameter_um) && !is.null(spot_diameter_px) && !is.null(microns_per_pixel)) {
    spot_diameter_um <- as.numeric(spot_diameter_px) * as.numeric(microns_per_pixel)
  }
  if (is.null(microns_per_pixel) && !is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
    microns_per_pixel <- as.numeric(spot_diameter_um) / as.numeric(spot_diameter_px)
  }
  if (is.null(microns_per_pixel) && is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
    default_spot_um <- getOption("geneSCOPE.default_visium_spot_um", 55)
    microns_per_pixel <- as.numeric(default_spot_um) / as.numeric(spot_diameter_px)
    spot_diameter_um  <- as.numeric(default_spot_um)
    if (isTRUE(verbose)) message("[geneSCOPE::createSCOPE visium] microns_per_pixel missing; assumed spot_diameter_microns = ", default_spot_um,
                                 ", computed microns_per_pixel = ", format(round(microns_per_pixel, 9), scientific = FALSE))
  }
  if (is.null(microns_per_pixel)) stop("scalefactors missing microns_per_pixel; cannot convert coordinates.")
  if (is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
    spot_diameter_um <- as.numeric(spot_diameter_px) * as.numeric(microns_per_pixel)
  }
  if (is.null(spot_diameter_um)) stop("scalefactors missing spot diameter information.")

  microns_per_pixel <- as.numeric(microns_per_pixel)
  base_spot_um <- as.numeric(spot_diameter_um)
  grid_size_um <- base_spot_um * as.numeric(grid_multiplier)

  # Convert pixel coordinates to microns
  pos_dt[, `:=`(
    x_um = pxl_col_in_fullres * microns_per_pixel,
    y_um = pxl_row_in_fullres * microns_per_pixel
  )]
  if (isTRUE(flip_y)) {
    y_max <- max(pos_dt$y_um, na.rm = TRUE)
    pos_dt[, y_um := y_max - y_um]
  }

  # 5) Build scope_object ---------------------------------------------------
  # Gene metadata: try features.tsv(.gz), else fallback to rownames
  candidate_dirs <- unique(file.path(visium_dir, c("filtered_feature_bc_matrix", "raw_feature_bc_matrix")))
  matrix_dir <- candidate_dirs[dir.exists(candidate_dirs)][1]
  feature_file <- if (!is.na(matrix_dir)) file.path(matrix_dir, "features.tsv.gz") else NA
  if (!is.na(feature_file) && !file.exists(feature_file)) feature_file <- file.path(matrix_dir, "features.tsv")
  if (!is.na(feature_file) && file.exists(feature_file)) {
    features_dt <- data.table::fread(feature_file, header = FALSE)
    gene_ids <- as.character(features_dt[[1]])
    gene_names <- make.unique(as.character(features_dt[[2]]))
    feature_type <- if (ncol(features_dt) >= 3) as.character(features_dt[[3]]) else rep(NA_character_, length(gene_names))
    # Align to matrix rows
    idx <- match(gene_names, rownames(raw_counts))
    ok <- !is.na(idx)
    gene_meta <- data.frame(
      feature_id = gene_ids[ok],
      feature_type = feature_type[ok],
      platform = rep("Visium", sum(ok)),
      stringsAsFactors = FALSE,
      row.names = gene_names[ok]
    )
  } else {
    gene_meta <- data.frame(
      feature_id = rownames(raw_counts),
      feature_type = NA_character_,
      platform = rep("Visium", nrow(raw_counts)),
      stringsAsFactors = FALSE,
      row.names = rownames(raw_counts)
    )
  }

  # centroids data.table
  centroids_dt <- data.table::data.table(
    cell = pos_dt$barcode,
    x = pos_dt$x_um,
    y = pos_dt$y_um,
    array_row = pos_dt$array_row,
    array_col = pos_dt$array_col,
    in_tissue = pos_dt$in_tissue
  )

  # grid_info using array_row/array_col ordering
  unique_cols <- sort(unique(pos_dt$array_col))
  unique_rows <- sort(unique(pos_dt$array_row))
  col_map <- setNames(seq_along(unique_cols), unique_cols)
  row_map <- setNames(seq_along(unique_rows), unique_rows)

  grid_dt <- data.table::data.table(
    grid_id = pos_dt$barcode,
    gx = as.integer(col_map[as.character(pos_dt$array_col)]),
    gy = as.integer(row_map[as.character(pos_dt$array_row)]),
    xmin = pos_dt$x_um - grid_size_um / 2,
    xmax = pos_dt$x_um + grid_size_um / 2,
    ymin = pos_dt$y_um - grid_size_um / 2,
    ymax = pos_dt$y_um + grid_size_um / 2,
    center_x = pos_dt$x_um,
    center_y = pos_dt$y_um,
    width = rep(grid_size_um, nrow(pos_dt)),
    height = rep(grid_size_um, nrow(pos_dt)),
    idx = seq_len(nrow(pos_dt)),
    array_row = pos_dt$array_row,
    array_col = pos_dt$array_col
  )
  xbins_eff <- length(unique(grid_dt$gx))
  ybins_eff <- length(unique(grid_dt$gy))

  # Reorder matrices to pos_dt order (important downstream)
  col_idx <- match(pos_dt$barcode, colnames(raw_counts))
  raw_counts <- raw_counts[, col_idx, drop = FALSE]
  colnames(raw_counts) <- pos_dt$barcode

  if (!is.null(sct_mat)) {
    col_idx2 <- match(pos_dt$barcode, colnames(sct_mat))
    sct_mat <- sct_mat[, col_idx2, drop = FALSE]
    colnames(sct_mat) <- pos_dt$barcode
  }

  # grid counts long table from raw counts
  counts_tbl <- Matrix::summary(raw_counts)
  counts_dt <- data.table::data.table(
    grid_id = colnames(raw_counts)[counts_tbl$j],
    gene = rownames(raw_counts)[counts_tbl$i],
    count = as.numeric(counts_tbl$x)
  )
  counts_dt <- counts_dt[count > 0]

  # Assemble scope_object
  fmt_len <- formatC(grid_size_um, format = "fg", digits = 6)
  fmt_len <- trimws(fmt_len)
  fmt_len <- sub("\\.?0+$", "", fmt_len)
  grid_name <- paste0("grid", fmt_len)

  # Infer hex orientation (flat-top vs pointy-top) from adjacency angles
  infer_hex_orientation <- function(cent, tol_deg = 15) {
    cn <- c("x","y","array_row","array_col")
    if (!all(cn %in% names(cent))) return(list(orientation = NA_character_, horiz_frac = NA_real_, vert_frac = NA_real_, tol_deg = tol_deg))
    dt <- data.table::as.data.table(cent)[, .(x, y, array_row = as.integer(array_row), array_col = as.integer(array_col))]
    get_pairs <- function(dr, dc) {
      a <- dt
      # Build neighbour table with unique column names to avoid x/y suffixing after merge
      b <- dt[, .(array_row_nb = array_row - dr,
                  array_col_nb = array_col - dc,
                  x_nb = x,
                  y_nb = y)]
      m <- merge(a, b, by.x = c("array_row", "array_col"), by.y = c("array_row_nb", "array_col_nb"), all = FALSE)
      if (!nrow(m)) return(NULL)
      m[, .(dx = x_nb - x, dy = y_nb - y)]
    }
    edges <- data.table::rbindlist(list(
      get_pairs(0, 1), get_pairs(1, 0), get_pairs(1, 1), get_pairs(1, -1)
    ), use.names = TRUE, fill = TRUE)
    if (is.null(edges) || !nrow(edges)) return(list(orientation = NA_character_, horiz_frac = NA_real_, vert_frac = NA_real_, tol_deg = tol_deg))
    theta <- abs(atan2(edges$dy, edges$dx) * 180 / pi)
    horiz <- pmin(theta, 180 - theta) <= tol_deg
    vert  <- abs(theta - 90) <= tol_deg
    hfrac <- mean(horiz, na.rm = TRUE)
    vfrac <- mean(vert,  na.rm = TRUE)
    list(orientation = ifelse(hfrac >= vfrac, "flat-top", "pointy-top"), horiz_frac = hfrac, vert_frac = vfrac, tol_deg = tol_deg)
  }
  ori <- infer_hex_orientation(centroids_dt, tol_deg = 15)

  scope_obj <- new("scope_object",
    coord = list(centroids = centroids_dt),
    grid = list(),
    meta.data = gene_meta,
    cells = list(counts = raw_counts),
    stats = list(),
    density = list()
  )
  scope_obj@grid[[grid_name]] <- list(
    grid_info = grid_dt,
    counts = counts_dt,
    grid_length = grid_size_um,
    xbins_eff = xbins_eff,
    ybins_eff = ybins_eff,
    spot_diameter_um = base_spot_um,
    microns_per_pixel = microns_per_pixel,
    visium_orientation = ori$orientation,
    visium_orientation_horiz_frac = ori$horiz_frac,
    visium_orientation_vert_frac  = ori$vert_frac,
    visium_orientation_tol_deg    = ori$tol_deg
  )
  if (!is.null(sct_mat)) {
    scope_obj@grid[[grid_name]]$SCT <- t(sct_mat)
  }

  if (isTRUE(verbose)) {
    message("[geneSCOPE::createSCOPE visium] scope_object created with ", nrow(centroids_dt),
            " spots and grid_length = ", round(grid_size_um, 3), " um")
  }

  scope_obj
}
