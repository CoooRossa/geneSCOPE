#' Modules:
#'   * .initialize_dataset_builder_state() – normalize arguments and shared params
#'   * .configure_worker_threads()         – thread detection + safe concurrency caps
#'   * .limit_thread_environment()         – clamp BLAS/OpenMP env vars with restore
#'   * .load_dataset_dependencies()        – package initialization + Arrow tuning
#'   * .resolve_input_paths()              – discover required parquet inputs
#'   * .load_centroid_table()              – read/flip centroids with data.table
#'   * .open_transcript_dataset()          – Arrow dataset open + filter wiring
#'   * .apply_dataset_filters()            – gene/QV/distance filtering + bounds
#'   * .prepare_roi_geometry()             – ROI polygon ingestion / auto rectangle
#'   * .clip_points_within_roi()           – clip centroids to ROI with chunking
#'   * .initialize_output_object()         – initialise output container coordinates
#'   * .ingest_segmentation_geometries()   – segmentation parquet ingestion helper
#'   * .prefetch_roi_molecules()           – optional ROI molecule prefetch
#'   * .scan_molecule_batches()            – streaming panoramic counting fallback
#'   * .build_grid_layers()                – multi-resolution grid assembly
#'   * .finalize_output_object()           – post checks + metadata stamping
#'   * createSCOPE_xenium()             – public wrapper stitching the pipeline
#'
#' The CosMx/Visium helpers reuse the same staging pattern via lightweight
#' orchestrators to keep behaviour identical to 1.CreateObj.r (2025-07-09).
#' Initialize dataset builder state.
#'
#' Normalises user-facing arguments and returns the shared pipeline state list
#' consumed by downstream helpers.
#'
#' @param input_dir Directory containing Xenium/CosMx inputs.
#' @param grid_length Numeric vector of grid sizes (microns).
#' @param seg_type Segmentation strategy (`"cell"`, `"nucleus"`, `"both"`, `"none"`).
#' @param filter_genes Optional character vector of genes to retain.
#' @param max_dist_mol_nuc Maximum molecule-to-nucleus distance.
#' @param filtermolecule Logical flag indicating molecule filtering.
#' @param exclude_prefix Gene prefixes to exclude.
#' @param filterqv Minimum molecule quality threshold.
#' @param coord_file Optional ROI coordinate file.
#' @param ncores Requested core count.
#' @param chunk_pts Chunk size for streaming operations.
#' @param min_gene_types Minimum gene diversity per grid.
#' @param max_gene_types Maximum gene diversity per grid.
#' @param min_seg_points Minimum segmentation vertex count.
#' @param keep_partial_grid Logical flag to retain partial grids.
#' @param verbose Logical toggle for progress logs.
#' @param flip_y Logical flag for Y-axis flipping.
#' @param data_type Character label for the dataset type.
#' @keywords internal
.initialize_dataset_builder_state <- function(
    input_dir,
    grid_length,
    seg_type,
    filter_genes,
    max_dist_mol_nuc,
    filtermolecule,
    exclude_prefix,
    filterqv,
    coord_file,
    ncores,
    chunk_pts,
    min_gene_types,
    max_gene_types,
    min_seg_points,
    keep_partial_grid,
    verbose,
    flip_y,
    data_type = "xenium") {
    seg_type <- match.arg(seg_type, c("cell", "nucleus", "both", "none"))
    if (!dir.exists(input_dir)) {
        stop("`input_dir` does not exist: ", input_dir)
    }
    stopifnot(is.numeric(grid_length), length(grid_length) >= 1, all(grid_length > 0))
    stopifnot(is.logical(flip_y), length(flip_y) == 1)
    stopifnot(is.logical(keep_partial_grid), length(keep_partial_grid) == 1)

    list(
        params = list(
            input_dir = input_dir,
            grid_length = grid_length,
            seg_type = seg_type,
            filter_genes = filter_genes,
            max_dist_mol_nuc = max_dist_mol_nuc,
            filtermolecule = filtermolecule,
            exclude_prefix = exclude_prefix,
            filterqv = filterqv,
            coord_file = coord_file,
            ncores = ncores,
            chunk_pts = chunk_pts,
            min_gene_types = min_gene_types,
            max_gene_types = max_gene_types,
            min_seg_points = min_seg_points,
            keep_partial_grid = keep_partial_grid,
            verbose = verbose,
            flip_y = flip_y,
            data_type = data_type,
            log_prefix = sprintf("[geneSCOPE::build_scope/%s]", data_type)
        ),
        thread_context = NULL,
        env_restore = NULL,
        paths = NULL,
        datasets = list(),
        roi = list(),
        objects = list()
    )
}

#' Configure worker threads for dataset construction.
#'
#' Probes the system to select a safe threading configuration and stores thread
#' context within the pipeline state.
#'
#' @param state Pipeline state list returned by `.initialize_dataset_builder_state()`.
#' @return Updated state with thread context metadata.
#' @keywords internal
.configure_worker_threads <- function(state) {
    p <- state$params
    thread_config <- configureThreadsFor("io_bound", p$ncores, restore_after = TRUE)
    restore_fn <- attr(thread_config, "restore_function")
    ncores_io <- thread_config$r_threads
    ncores_cpp <- thread_config$openmp_threads
    log_prefix <- p$log_prefix

    os_type <- detectOS()
    max_cores <- parallel::detectCores()

    sys_mem_gb <- tryCatch({
        if (os_type == "linux") {
            mem_info <- readLines("/proc/meminfo")
            mem_total <- grep("MemTotal", mem_info, value = TRUE)
            mem_kb <- as.numeric(gsub("[^0-9]", "", mem_total))
            mem_gb <- mem_kb / 1024 / 1024
            if (mem_gb > 50000) warning("Detected unusually high memory: ", round(mem_gb, 1), "GB. This may indicate a parsing error.")
            mem_gb
        } else if (os_type == "macos") {
            mem_bytes <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
            mem_bytes / 1024^3
        } else {
            32
        }
    }, error = function(e) 32)

    mem_based_limit <- if (sys_mem_gb > 500) {
        min(p$ncores, max_cores - 2)
    } else {
        switch(cut(sys_mem_gb, breaks = c(0, 8, 32, 64, 128, Inf), labels = FALSE),
            1,
            min(4, max_cores - 1),
            min(8, max_cores - 1),
            min(16, max_cores - 1),
            min(p$ncores, max_cores - 2)
        )
    }

    platform_limit <- switch(os_type,
        windows = min(8, max_cores - 1),
        macos = min(12, max_cores - 1),
        linux = if (sys_mem_gb > 100) min(p$ncores, max_cores - 2) else min(32, max_cores - 1)
    )

    test_cores <- min(p$ncores, mem_based_limit, platform_limit)
    ncores_safe <- 1L
    test_step <- if (sys_mem_gb > 100 || max_cores > 32) 8 else 2

    for (test_nc in seq(test_cores, 1, by = -test_step)) {
        ok <- tryCatch({
            if (test_nc > 1) {
                if (sys_mem_gb > 100) {
                    TRUE
                } else {
                    cl <- parallel::makeCluster(min(test_nc, 4))
                    on.exit(parallel::stopCluster(cl), add = TRUE)
                    parallel::clusterEvalQ(cl, 1 + 1)
                    parallel::stopCluster(cl)
                }
            }
            TRUE
        }, error = function(e) {
            if (p$verbose) message(log_prefix, " Testing ", test_nc, " cores failed, trying fewer")
            FALSE
        })
        if (isTRUE(ok)) {
            ncores_safe <- test_nc
            break
        }
    }

    if (p$verbose) {
        mem_display <- if (sys_mem_gb >= 1024) paste0(round(sys_mem_gb / 1024, 1), "TB") else paste0(round(sys_mem_gb, 1), "GB")
        if (ncores_safe < p$ncores) {
            message(log_prefix, " Core configuration: using ", ncores_safe, "/", p$ncores, " cores (", mem_display, " RAM, ", os_type, ")")
        } else {
            message(log_prefix, " Core configuration: using ", ncores_safe, " cores")
        }
    }

    state$thread_context <- list(
        thread_config = thread_config,
        restore_fn = restore_fn,
        ncores_safe = ncores_safe,
        ncores_io = ncores_io,
        ncores_cpp = ncores_cpp,
        os_type = os_type,
        sys_mem_gb = sys_mem_gb
    )
    state
}

#' Clamp threaded environment variables.
#'
#' Temporarily adjusts BLAS/OpenMP environment settings to avoid oversubscription
#' and records a restore handler in the pipeline state.
#'
#' @param state Pipeline state with thread context metadata.
#' @return Updated state containing environment restoration callback.
#' @keywords internal
.limit_thread_environment <- function(state) {
    p <- state$params
    ctx <- state$thread_context
    log_prefix <- p$log_prefix

    thread_vars <- c(
        "OMP_NUM_THREADS", "OMP_THREAD_LIMIT",
        "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
        "BLAS_NUM_THREADS", "LAPACK_NUM_THREADS"
    )

    old_env <- sapply(thread_vars, Sys.getenv, unset = NA_character_)
    old_schedule <- Sys.getenv("OMP_SCHEDULE", unset = NA_character_)

    restore_env <- function() {
        for (var in names(old_env)) {
            if (is.na(old_env[[var]])) {
                Sys.unsetenv(var)
            } else {
                args <- list(); args[[var]] <- old_env[[var]]
                do.call(Sys.setenv, args)
            }
        }
        if (is.na(old_schedule)) {
            Sys.unsetenv("OMP_SCHEDULE")
        } else {
            Sys.setenv(OMP_SCHEDULE = old_schedule)
        }
    }

    for (var in thread_vars) {
        args <- list(); args[[var]] <- "1"
        do.call(Sys.setenv, args)
    }
    if (ctx$os_type == "linux") Sys.setenv(OMP_SCHEDULE = "static")

    if (requireNamespace("data.table", quietly = TRUE)) {
        tryCatch(data.table::setDTthreads(1), error = function(e) {
            if (p$verbose) message(log_prefix, " Warning: Could not set data.table threads")
        })
    }

    state$env_restore <- restore_env
    state
}

#' Load package dependencies and configure Arrow.
#'
#' Ensures required packages are available, tunes Arrow threading and appends
#' dependency metadata to the pipeline state.
#'
#' @param state Pipeline state list.
#' @return Updated state with dependency configuration.
#' @keywords internal
.load_dataset_dependencies <- function(state) {
    p <- state$params
    ctx <- state$thread_context
    log_prefix <- p$log_prefix

    if (p$verbose) message(log_prefix, " Initializing data processing environment")
    suppressPackageStartupMessages({
        library(arrow)
        library(data.table)
        library(sf)
        library(parallel)
        library(dplyr)
    })

    if (requireNamespace("arrow", quietly = TRUE)) {
        tryCatch({
            options(arrow.num_threads = ctx$ncores_safe)
            if (exists("set_cpu_count", where = asNamespace("arrow"))) arrow::set_cpu_count(ctx$ncores_safe)
            if (exists("set_io_thread_count", where = asNamespace("arrow"))) arrow::set_io_thread_count(ctx$ncores_safe)
        }, error = function(e) {
            if (p$verbose) message(log_prefix, " Warning: Could not configure Arrow threading")
        })
    }
    state
}

#' Resolve canonical input file paths.
#'
#' Validates the expected transcripts, cells and segmentation parquet locations
#' and stores them under `state$paths`.
#'
#' @param state Pipeline state list.
#' @return Updated state with resolved path entries.
#' @keywords internal
.resolve_input_paths <- function(state) {
    p <- state$params
    tf_parq <- file.path(p$input_dir, "transcripts.parquet")
    cell_parq <- file.path(p$input_dir, "cells.parquet")
    if (!file.exists(tf_parq)) stop("transcripts.parquet not found in ", p$input_dir)
    if (!file.exists(cell_parq)) stop("cells.parquet not found in ", p$input_dir)

    seg_files <- list()
    if (identical(p$seg_type, "none")) {
        # no segmentation required
    } else if (p$seg_type %in% c("cell", "both")) {
        cell_seg_files <- c(
            file.path(p$input_dir, "cell_boundaries.parquet"),
            file.path(p$input_dir, "segmentation_boundaries.parquet")
        )
        cell_seg_found <- sapply(cell_seg_files, file.exists)
        if (any(cell_seg_found)) {
            seg_files$cell <- cell_seg_files[cell_seg_found][1]
        } else {
            stop("Cell segmentation parquet file not found in ", p$input_dir)
        }
    }
    if (!identical(p$seg_type, "none") && p$seg_type %in% c("nucleus", "both")) {
        nuc_file <- file.path(p$input_dir, "nucleus_boundaries.parquet")
        if (!file.exists(nuc_file)) stop("Nucleus segmentation parquet file not found in ", p$input_dir)
        seg_files$nucleus <- nuc_file
    }

    state$paths <- list(
        tf_parq = tf_parq,
        cell_parq = cell_parq,
        seg_files = seg_files
    )
    state
}

#' Load centroid coordinates.
#'
#' Reads centroid tables, applies optional flips and stores them in
#' `state$datasets$centroids`.
#'
#' @param state Pipeline state list.
#' @return Updated state with centroid data.
#' @keywords internal
.load_centroid_table <- function(state) {
    p <- state$params
    ctx <- state$thread_context
    paths <- state$paths
    log_prefix <- p$log_prefix
    if (p$verbose) message(log_prefix, " Loading cell centroids and transcripts")

    centroids_dt <- tryCatch({
        if (requireNamespace("arrow", quietly = TRUE)) {
            arrow::set_cpu_count(ctx$ncores_io)
            arrow::set_io_thread_count(ctx$ncores_io)
        }
        cell_data <- arrow::read_parquet(paths$cell_parq)
        cd <- as.data.table(cell_data)
        if ("fov" %in% names(cd)) {
            cd[, cell := paste0(as.character(cell_id), "_", as.character(fov))]
        } else {
            cd[, cell := as.character(cell_id)]
        }
        cd[, .(cell = cell, x = x_centroid, y = y_centroid)]
    }, error = function(e) {
        stop("Error reading centroids from cells.parquet: ", e$message)
    })

    if (nrow(centroids_dt) == 0) stop("No centroids found in cells.parquet")
    if (p$verbose) message(log_prefix, "  Found ", nrow(centroids_dt), " cells")

    state$objects$centroids <- centroids_dt
    state
}

#' Open transcript Arrow dataset.
#'
#' Creates an Arrow dataset handle for transcripts and stores it for later
#' filtering and aggregation.
#'
#' @param state Pipeline state list.
#' @return Updated state with Arrow dataset handle.
#' @keywords internal
.open_transcript_dataset <- function(state) {
    paths <- state$paths
    ds <- tryCatch(arrow::open_dataset(paths$tf_parq), error = function(e) {
        stop("Error opening transcripts dataset: ", e$message)
    })
    state$datasets$ds <- ds
    state
}

#' Apply dataset-level filters.
#'
#' Executes gene, distance and quality filters against the transcript dataset
#' while computing global bounds.
#'
#' @param state Pipeline state list.
#' @return Updated state with filtered dataset references.
#' @keywords internal
.apply_dataset_filters <- function(state) {
    p <- state$params
    ds <- state$datasets$ds
    log_prefix <- p$log_prefix

    if (!is.null(p$filter_genes)) {
        if (p$verbose) message(log_prefix, "  Filtering genes: ", length(p$filter_genes), " genes")
        ds <- ds |> dplyr::filter(feature_name %in% p$filter_genes)
    }
    if (!is.null(p$max_dist_mol_nuc)) {
        if (p$verbose) message(log_prefix, "  Filtering by nucleus distance: max ", p$max_dist_mol_nuc)
        ds <- ds |> dplyr::filter(nucleus_distance <= p$max_dist_mol_nuc)
    }
    if (isTRUE(p$filtermolecule)) {
        if (p$verbose) message(log_prefix, "  Filtering molecules by prefix exclusion")
        pat <- paste0("^(?:", paste(p$exclude_prefix, collapse = "|"), ")")
        ds <- ds |> dplyr::filter(!grepl(pat, feature_name))
    }
    if (!is.null(p$filterqv)) {
        if (p$verbose) message(log_prefix, "  Filtering by QV score: min ", p$filterqv)
        ds <- ds |> dplyr::filter(qv >= p$filterqv)
    }

    bounds_global <- tryCatch({
        ds |>
            dplyr::summarise(
                xmin = min(x_location), xmax = max(x_location),
                ymin = min(y_location), ymax = max(y_location)
            ) |>
            dplyr::collect()
    }, error = function(e) {
        stop("Error computing bounds: ", e$message)
    })

    y_max <- as.numeric(bounds_global$ymax)
    if (p$verbose) {
        message(
            log_prefix, " Data bounds: x=[", round(bounds_global$xmin), ",",
            round(bounds_global$xmax), "], y=[", round(bounds_global$ymin),
            ",", round(bounds_global$ymax), "]"
        )
    }

    if (p$flip_y) {
        state$objects$centroids <- flip_y_coordinates(state$objects$centroids, y_max)
        if (p$verbose) message(log_prefix, "  Applied Y-axis flip")
    }

    state$datasets$ds <- ds
    state$bounds <- list(bounds_global = bounds_global, y_max = y_max)
    state
}

#' Prepare ROI geometry.
#'
#' Ingests explicit ROI polygons or synthesises defaults, storing geometry
#' structures inside the pipeline state.
#'
#' @param state Pipeline state list.
#' @return Updated state containing ROI metadata.
#' @keywords internal
.prepare_roi_geometry <- function(state) {
    p <- state$params
    ds <- state$datasets$ds
    bounds <- state$bounds$bounds_global
    y_max <- state$bounds$y_max
    verbose <- p$verbose
    flip_y <- p$flip_y
    chunk_pts <- p$chunk_pts
    log_prefix <- p$log_prefix

    user_poly <- NULL
    if (!is.null(p$coord_file)) {
        if (verbose) message(log_prefix, " Processing ROI polygon from file")
        if (!file.exists(p$coord_file)) stop("Coordinate file not found: ", p$coord_file)

        ext <- tolower(tools::file_ext(p$coord_file))
        sep_char <- switch(ext,
            csv = ",",
            tsv = "\t",
            txt = "\t",
            "\t"
        )
        poly_tb <- tryCatch(utils::read.table(p$coord_file, header = TRUE, sep = sep_char, stringsAsFactors = FALSE), error = function(e) {
            stop("Error reading coordinate file: ", e$message)
        })
        names(poly_tb) <- tolower(names(poly_tb))
        if (!all(c("x", "y") %in% names(poly_tb))) stop("Coordinate file must contain 'x' and 'y' columns")
        if (nrow(poly_tb) < 3) stop("ROI polygon must have at least 3 points")
        if (verbose) message(log_prefix, "  Polygon has ", nrow(poly_tb), " vertices")

        poly_geom <- tryCatch({
            geom <- sf::st_polygon(list(as.matrix(poly_tb[, c("x", "y")])))
            sf::st_make_valid(geom)
        }, error = function(e) stop("Error creating polygon geometry: ", e$message))
        if (length(poly_geom) == 0) stop("ROI polygon became EMPTY after validation")

        if (flip_y) {
            poly_tb <- flip_y_coordinates(poly_tb, y_max)
            poly_geom <- tryCatch(sf::st_make_valid(sf::st_polygon(list(as.matrix(poly_tb[, c("x", "y")])))), error = function(e) {
                stop("Error creating flipped polygon: ", e$message)
            })
        }
        user_poly <- sf::st_sfc(poly_geom)
    } else {
        if (verbose) message(log_prefix, " Creating automatic ROI from data boundaries")
        poly_coords <- data.frame(
            x = c(bounds$xmin, bounds$xmax, bounds$xmax, bounds$xmin, bounds$xmin),
            y = c(bounds$ymin, bounds$ymin, bounds$ymax, bounds$ymax, bounds$ymin)
        )
        if (flip_y) poly_coords <- flip_y_coordinates(poly_coords, y_max)
        poly_geom <- tryCatch({
            geom <- sf::st_polygon(list(as.matrix(poly_coords[, c("x", "y")])))
            sf::st_make_valid(geom)
        }, error = function(e) stop("Error creating automatic ROI polygon: ", e$message))
        if (length(poly_geom) == 0) stop("Automatic ROI polygon became EMPTY after validation")
        user_poly <- sf::st_sfc(poly_geom)
        if (verbose) {
            message(
                log_prefix, " Automatic ROI created: x=[", round(min(poly_coords$x)), ",",
                round(max(poly_coords$x)), "], y=[", round(min(poly_coords$y)), ",",
                round(max(poly_coords$y)), "]"
            )
        }
    }

    state$roi <- list(
        geometry = user_poly,
        chunk_pts = chunk_pts
    )
    state
}

#' Clip centroids to ROI boundaries.
#'
#' Filters centroid coordinates to the prepared ROI polygon, retaining only
#' points that fall inside the analysis area.
#'
#' @param state Pipeline state list.
#' @return Updated state with clipped centroids.
#' @keywords internal
.clip_points_within_roi <- function(state) {
    p <- state$params
    roi <- state$roi$geometry
    centroids_dt <- state$objects$centroids
    log_prefix <- p$log_prefix
    if (!is.null(roi) && nrow(centroids_dt)) {
        if (p$verbose) message(log_prefix, " Clipping centroids to ROI...")
        centroids_dt <- tryCatch(
            clip_points_to_region(centroids_dt, roi, chunk_size = p$chunk_pts, workers = state$thread_context$ncores_safe),
            error = function(e) stop("Error clipping centroids: ", e$message)
        )
        if (p$verbose) message(log_prefix, " Retained ", nrow(centroids_dt), " cells within ROI")
    }
    state$objects$centroids <- centroids_dt
    state$objects$keep_cells <- if (nrow(centroids_dt)) unique(centroids_dt$cell) else NULL
    state
}

#' Initialize the output scope object.
#'
#' Creates and seeds the `scope_object` container in the pipeline state using
#' preprocessed centroid data.
#'
#' @param state Pipeline state list.
#' @return Updated state with initialised scope object.
#' @keywords internal
.initialize_output_object <- function(state) {
    centroids_dt <- state$objects$centroids
    scope_obj <- new("scope_object",
        coord = list(centroids = centroids_dt),
        grid = list(),
        meta.data = data.frame()
    )
    state$objects$scope_obj <- scope_obj
    state
}

#' Ingest segmentation geometries.
#'
#' Loads segmentation parquet files, converts polygons and writes them into the
#' scope object stored in the pipeline state.
#'
#' @param state Pipeline state list.
#' @return Updated state with segmentation layers.
#' @keywords internal
.ingest_segmentation_geometries <- function(state) {
    p <- state$params
    seg_files <- state$paths$seg_files
    scope_obj <- state$objects$scope_obj
    keep_cells <- state$objects$keep_cells
    y_max <- state$bounds$y_max
    log_prefix <- p$log_prefix

    if (!length(seg_files)) {
        if (p$verbose) message(log_prefix, " No segmentation files detected; skipping segmentation ingestion")
        state$objects$scope_obj <- scope_obj
        return(state)
    }

    if (p$verbose) message(log_prefix, " Processing segmentation data...")
    for (tag in names(seg_files)) {
        if (p$verbose) message(log_prefix, " Processing ", tag, " segmentation...")
        seg_res <- tryCatch({
            build_segmentation_geometries(
                path = seg_files[[tag]],
                label = tag,
                flip = p$flip_y,
                y_max = y_max,
                keep_ids = keep_cells,
                workers = state$thread_context$ncores_safe
            )
        }, error = function(e) {
            warning("Error processing ", tag, " segmentation: ", e$message)
            list(points = data.table(), polygons = list())
        })
        scope_obj@coord[[paste0("segmentation_", tag)]] <- seg_res$points
        scope_obj@coord[[paste0("segmentation_polys_", tag)]] <- seg_res$polygons
    }
    state$objects$scope_obj <- scope_obj
    state
}

#' Prefetch molecules within ROI.
#'
#' Optionally collects molecules intersecting the ROI to accelerate downstream
#' aggregation steps.
#'
#' @param state Pipeline state list.
#' @return Updated state with prefetched molecule data when requested.
#' @keywords internal
.prefetch_roi_molecules <- function(state) {
    p <- state$params
    roi <- state$roi$geometry
    ds <- state$datasets$ds
    y_max <- state$bounds$y_max
    mol_small <- NULL
    log_prefix <- p$log_prefix
    if (!is.null(roi)) {
        if (p$verbose) message(log_prefix, " Pre-fetching molecules for ROI...")
        mol_small <- tryCatch({
            ds |>
                dplyr::select(x = x_location, y = y_location, feature_name) |>
                dplyr::collect() |>
                as.data.table()
        }, error = function(e) {
            warning("Error pre-fetching molecules: ", e$message)
            data.table(x = numeric(), y = numeric(), feature_name = character())
        })
        if (p$flip_y && nrow(mol_small) > 0) mol_small[, y := y_max - y]
        if (nrow(mol_small) > 0) {
            mol_small <- clip_points_to_region(mol_small, roi, chunk_size = p$chunk_pts, workers = state$thread_context$ncores_safe)
        }
        if (p$verbose) message(log_prefix, " Retained ", nrow(mol_small), " molecules within ROI")
    }
    state$objects$mol_small <- mol_small
    state
}

#' Perform streaming molecule scans.
#'
#' Iteratively scans transcripts in batches to compute summary statistics and
#' attaches results to the pipeline state.
#'
#' @param state Pipeline state list.
#' @return Updated state with molecule scan outputs.
#' @keywords internal
.scan_molecule_batches <- function(state) {
    p <- state$params
    ds <- state$datasets$ds
    bounds <- state$bounds$bounds_global
    y_max <- state$bounds$y_max
    counts_list <- setNames(vector("list", length(p$grid_length)), paste0("lg", p$grid_length))
    for (i in seq_along(counts_list)) {
        counts_list[[i]] <- data.table(grid_id = character(), gene = character(), count = integer())
    }
    log_prefix <- p$log_prefix
    if (is.null(state$roi$geometry)) {
        if (p$verbose) message(log_prefix, " Starting multi-resolution scan...")
        n_total <- tryCatch({
            ds |> dplyr::summarise(n = dplyr::n()) |> dplyr::collect() |> dplyr::pull(n)
        }, error = function(e) NA_integer_)
        pb <- if (p$verbose && !is.na(n_total)) utils::txtProgressBar(min = 0, max = n_total, style = 3) else NULL
        scanned <- 0L

        repeat {
            tbl <- tryCatch({
                if (scanned == 0) {
                    ds |>
                        dplyr::select(x_location, y_location, feature_name) |>
                        dplyr::slice_head(n = p$chunk_pts) |>
                        dplyr::collect()
                } else {
                    batch <- ds |>
                        dplyr::select(x_location, y_location, feature_name) |>
                        dplyr::slice(scanned + 1, scanned + p$chunk_pts) |>
                        dplyr::collect()
                    if (nrow(batch) == 0) NULL else batch
                }
            }, error = function(e) {
                if (p$verbose) message(log_prefix, " Scanner finished or error: ", e$message)
                NULL
            })
            if (is.null(tbl) || nrow(tbl) == 0) break

            scanned <- scanned + nrow(tbl)
            if (!is.null(pb)) utils::setTxtProgressBar(pb, min(scanned, n_total))

            dt <- as.data.table(tbl)
            if (p$flip_y) dt[, y_location := y_max - y_location]
            for (k in seq_along(p$grid_length)) {
                lg <- p$grid_length[k]
                x0_lg <- floor(bounds$xmin / lg) * lg
                y0_lg <- if (p$flip_y) 0 else floor(bounds$ymin / lg) * lg
                dt[, `:=`(
                    gx = (x_location - x0_lg) %/% lg + 1L,
                    gy = (y_location - y0_lg) %/% lg + 1L
                )]
                dt[, grid_id_tmp := paste0("g", gx, "_", gy)]
                accum <- dt[, .N, by = .(grid_id = grid_id_tmp, gene = feature_name)]
                setnames(accum, "N", "count")
                counts_list[[k]] <- rbindlist(list(counts_list[[k]], accum))
            }
            if (!is.na(n_total) && scanned >= n_total) break
        }
        if (!is.null(pb)) close(pb)
        for (k in seq_along(counts_list)) {
            counts_list[[k]] <- counts_list[[k]][, .(count = sum(count)), by = .(grid_id, gene)]
        }
    }
    state$objects$counts_list <- counts_list
    state
}

#' Build multi-resolution grid layers.
#'
#' Generates spatial grids across requested resolutions and writes them onto the
#' scope object.
#'
#' @param state Pipeline state list.
#' @return Updated state with populated grid layers.
#' @keywords internal
.build_grid_layers <- function(state) {
    p <- state$params
    scope_obj <- state$objects$scope_obj
    bounds <- state$bounds$bounds_global
    y_max <- state$bounds$y_max
    mol_small <- state$objects$mol_small
    counts_list <- state$objects$counts_list
    flip_y <- p$flip_y
    log_prefix <- p$log_prefix

    if (p$verbose) message(log_prefix, " Building grid layers …")
    lg_unique <- sort(unique(p$grid_length))

    for (lg in lg_unique) {
        if (p$verbose) message(sprintf("%s Grid %.1f µm …", log_prefix, lg))
        x0 <- floor(bounds$xmin / lg) * lg
        y0 <- if (flip_y) 0 else floor(bounds$ymin / lg) * lg

        cnt <- if (is.null(state$roi$geometry)) {
            counts_list[[paste0("lg", lg)]]
        } else {
            mol_dt <- copy(mol_small)
            if (nrow(mol_dt)) {
                mol_dt[, `:=`(
                    gx = (x - x0) %/% lg + 1L,
                    gy = (y - y0) %/% lg + 1L
                )]
                mol_dt[, grid_id := paste0("g", gx, "_", gy)]
                accum <- mol_dt[, .N, by = .(grid_id, gene = feature_name)]
                setnames(accum, "N", "count")
            } else {
                data.table(grid_id = character(), gene = character(), count = integer())
            }
        }

        gene_div <- cnt[, uniqueN(gene), by = grid_id]
        keep_grid <- gene_div[V1 >= p$min_gene_types & V1 <= p$max_gene_types, grid_id]
        cnt <- cnt[grid_id %in% keep_grid]

        if (p$min_seg_points > 0L) {
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
                valid_seg <- seg_cnt[N >= p$min_seg_points, grid_id]
                if (length(valid_seg)) {
                    keep_grid <- intersect(keep_grid, valid_seg)
                    cnt <- cnt[grid_id %in% keep_grid]
                }
            }
        }

        x_breaks <- seq(x0, bounds$xmax, by = lg)
        if (tail(x_breaks, 1L) < bounds$xmax) x_breaks <- c(x_breaks, bounds$xmax)
        max_y_use <- if (flip_y) y_max else bounds$ymax
        y_breaks <- seq(y0, max_y_use, by = lg)
        if (tail(y_breaks, 1L) < max_y_use) y_breaks <- c(y_breaks, max_y_use)

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
            width = xmax - xmin,
            height = ymax - ymin,
            idx = .I
        )]
        if (!p$keep_partial_grid) {
            grid_dt <- grid_dt[abs(width - lg) < 1e-6 & abs(height - lg) < 1e-6]
        }
        grid_dt <- grid_dt[grid_id %in% keep_grid]
        cnt <- cnt[grid_id %in% grid_dt$grid_id]

        scope_obj@grid[[paste0("grid", lg)]] <- list(
            grid_info = grid_dt,
            counts = cnt,
            grid_length = lg,
            xbins_eff = length(x_breaks) - 1L,
            ybins_eff = length(y_breaks) - 1L,
            histology = list(
                lowres = NULL,
                hires = NULL
            )
        )
    }
    if (p$verbose) message(log_prefix, " Grid construction finished.")
    state$objects$scope_obj <- scope_obj
    state
}

#' Finalise the scope object.
#'
#' Performs final integrity checks, stamps metadata and prepares the object for
#' return to the caller.
#'
#' @param state Pipeline state list.
#' @return Updated state with fully initialised `scope_object`.
#' @keywords internal
.finalize_output_object <- function(state) {
    scope_obj <- state$objects$scope_obj
    type_raw <- state$params$data_type
    platform_label <- switch(tolower(type_raw),
        xenium = "Xenium",
        cosmx = "CosMx",
        visium = "Visium",
        `default` = type_raw
    )
    if (length(scope_obj@grid) == 0 ||
        all(vapply(scope_obj@grid, function(g) nrow(g$grid_info) == 0L, logical(1)))) {
        warning("No effective grid layer generated – please check filters/parameters.")
    }
    if (nrow(scope_obj@meta.data) > 0L) {
        if ("platform" %in% names(scope_obj@meta.data)) {
            scope_obj@meta.data$platform[] <- platform_label
        } else {
            scope_obj@meta.data$platform <- rep(platform_label, nrow(scope_obj@meta.data))
        }
    } else {
        scope_obj@meta.data <- data.frame(platform = platform_label, stringsAsFactors = FALSE)
        rownames(scope_obj@meta.data) <- "__scope_platform__"
    }
    scope_obj@stats$platform <- platform_label
    state$objects$scope_obj <- scope_obj
    state
}

#' Build a scope_object from Xenium outputs.
#'
#' Executes the modular builder pipeline on a Xenium outs directory to produce
#' a fully initialised `scope_object`.
#'
#' @param input_dir Path to the Xenium outs directory.
#' @param grid_sizes Numeric vector of grid sizes (microns).
#' @param segmentation_strategy Segmentation strategy (`"cell"`, `"nucleus"`, `"both"`).
#' @param filter_genes Optional gene allowlist.
#' @param max_distance_to_nucleus Maximum molecule-to-nucleus distance.
#' @param filter_molecules Logical flag enabling molecule-level filtering.
#' @param exclude_prefixes Gene prefixes to exclude.
#' @param min_quality Minimum molecule quality value.
#' @param roi_path Optional ROI coordinate file.
#' @param num_workers Requested core count.
#' @param chunk_size Chunk size for streaming reads.
#' @param min_gene_types Minimum gene diversity per grid.
#' @param max_gene_types Maximum gene diversity per grid.
#' @param min_segmentation_points Minimum segmentation vertex count.
#' @param keep_partial_tiles Logical flag to retain partial tiles.
#' @param verbose Logical toggle for progress output.
#' @param flip_y Logical flag to flip the Y axis.
#' @param data_type Character label stored on the output object.
#' @return A populated `scope_object`.
#' @export
build_scope_from_xenium <- function(input_dir,
                                    grid_sizes,
                                    segmentation_strategy = c("cell", "nucleus", "both"),
                                    filter_genes = NULL,
                                    max_distance_to_nucleus = 25,
                                    filter_molecules = TRUE,
                                    exclude_prefixes = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword",
                                                         "SystemControl", "Negative"),
                                    min_quality = 20,
                                    roi_path = NULL,
                                    num_workers = 1L,
                                    chunk_size = 5e5,
                                    min_gene_types = 1,
                                    max_gene_types = Inf,
                                    min_segmentation_points = 1,
                                    keep_partial_tiles = FALSE,
                                    verbose = TRUE,
                                    flip_y = TRUE,
                                    data_type = "xenium") {

    state <- .initialize_dataset_builder_state(
        input_dir = input_dir,
        grid_length = grid_sizes,
        seg_type = segmentation_strategy,
        filter_genes = filter_genes,
        max_dist_mol_nuc = max_distance_to_nucleus,
        filtermolecule = filter_molecules,
        exclude_prefix = exclude_prefixes,
        filterqv = min_quality,
        coord_file = roi_path,
        ncores = num_workers,
        chunk_pts = chunk_size,
        min_gene_types = min_gene_types,
        max_gene_types = max_gene_types,
        min_seg_points = min_segmentation_points,
        keep_partial_grid = keep_partial_tiles,
        verbose = verbose,
        flip_y = flip_y,
        data_type = data_type
    )

    state <- .configure_worker_threads(state)
    if (!is.null(state$thread_context$restore_fn)) {
        on.exit(state$thread_context$restore_fn(), add = TRUE)
    }
    state <- .limit_thread_environment(state)
    if (!is.null(state$env_restore)) on.exit(state$env_restore(), add = TRUE)

    state <- .load_dataset_dependencies(state)
    state <- .resolve_input_paths(state)
    state <- .load_centroid_table(state)
    state <- .open_transcript_dataset(state)
    state <- .apply_dataset_filters(state)
    state <- .prepare_roi_geometry(state)
    state <- .clip_points_within_roi(state)
    state <- .initialize_output_object(state)
    state <- .ingest_segmentation_geometries(state)
    state <- .prefetch_roi_molecules(state)
    state <- .scan_molecule_batches(state)
    state <- .build_grid_layers(state)
    state <- .finalize_output_object(state)

    state$objects$scope_obj
}

#' @describeIn build_scope_from_xenium Backward-compatible wrapper retaining legacy argument names.
#' @export
createSCOPE_xenium <- function(xenium_dir,
                               grid_length,
                               seg_type = c("cell", "nucleus", "both"),
                               data_type = "xenium",
                               ...) {
    seg_type <- match.arg(seg_type)
    dots <- list(...)
    rename_map <- list(
        coord_file = "roi_path",
        ncores = "num_workers",
        chunk_pts = "chunk_size",
        min_seg_points = "min_segmentation_points",
        keep_partial_grid = "keep_partial_tiles",
        filtermolecule = "filter_molecules",
        exclude_prefix = "exclude_prefixes",
        filterqv = "min_quality",
        max_dist_mol_nuc = "max_distance_to_nucleus"
    )
    for (old in names(rename_map)) {
        if (!is.null(dots[[old]]) && is.null(dots[[rename_map[[old]]]])) {
            dots[[rename_map[[old]]]] <- dots[[old]]
        }
        dots[[old]] <- NULL
    }
    build_args <- c(
        list(
            input_dir = xenium_dir,
            grid_sizes = grid_length,
            segmentation_strategy = seg_type,
            data_type = data_type
        ),
        dots
    )
    do.call(build_scope_from_xenium, build_args)
}

#' Build a scope_object from CosMx outputs.
#'
#' Prepares CosMx flatFiles/parquet inputs (generating intermediates if needed)
#' and then delegates to the modular Xenium pipeline to construct a
#' `scope_object`.
#'
#' @param input_dir Path to the CosMx project root.
#' @param grid_sizes Numeric vector of grid sizes (microns).
#' @param segmentation_strategy Segmentation strategy (`"cell"`, `"nucleus"`, `"both"`, `"none"`).
#' @param dataset_id Optional dataset identifier override.
#' @param pixel_size_um Pixel size used for unit conversion.
#' @param verbose Logical toggle for progress logs.
#' @param allow_flatfile_generation Logical allowing parquet generation from flatFiles.
#' @param derive_cells_from_polygons Logical controlling cell derivation fallback.
#' @param data_type Character label stored on the output object.
#' @param ... Additional arguments forwarded to the builder stages.
#' @return A populated `scope_object`.
#' @importFrom arrow write_parquet
#' @export
build_scope_from_cosmx <- function(input_dir,
                                   grid_sizes,
                                   segmentation_strategy = c("cell", "nucleus", "both"),
                                   dataset_id = NULL,
                                   pixel_size_um = 0.120280945,
                                   verbose = TRUE,
                                   allow_flatfile_generation = TRUE,
                                   derive_cells_from_polygons = TRUE,
                                   data_type = "cosmx",
                                   ...) {
    stopifnot(dir.exists(input_dir))
    segmentation_strategy <- match.arg(segmentation_strategy, c("cell", "nucleus", "both", "none"))

    if (verbose) message("[geneSCOPE::build_scope_from_cosmx] Root = ", input_dir)

    transcripts_path <- file.path(input_dir, "transcripts.parquet")
    cells_path <- file.path(input_dir, "cells.parquet")
    primary_cell_boundary <- file.path(input_dir, "cell_boundaries.parquet")
    fallback_cell_boundary <- file.path(input_dir, "segmentation_boundaries.parquet")
    nucleus_boundary <- file.path(input_dir, "nucleus_boundaries.parquet")

    working_dir <- input_dir
    using_temp_dir <- FALSE

    grid_length <- grid_sizes
    seg_type <- segmentation_strategy
    build_from_flatfiles <- allow_flatfile_generation
    build_cells_from_polygons <- derive_cells_from_polygons
    tf_parq <- transcripts_path
    cell_parq <- cells_path
    seg_cell <- primary_cell_boundary
    seg_any <- fallback_cell_boundary
    nuc_parq <- nucleus_boundary
    work_dir <- working_dir
    need_temp <- using_temp_dir

    if (!file.exists(transcripts_path) && allow_flatfile_generation) {
        if (verbose) message("[geneSCOPE::build_scope_from_cosmx] transcripts.parquet not found; generating from flatFiles …")
        ff_dir <- file.path(input_dir, "flatFiles")
        if (!dir.exists(ff_dir)) stop("flatFiles/ not found under ", input_dir)
        tx_candidates <- c(
            Sys.glob(file.path(ff_dir, "*/*_tx_file.csv.gz")),
            Sys.glob(file.path(ff_dir, "*/*_tx_file.csv"))
        )
        if (length(tx_candidates) == 0L) {
            stop("No *_tx_file.csv(.gz) found under ", ff_dir)
        }
        tx_path <- tx_candidates[[1]]
        if (is.null(dataset_id)) {
            base <- basename(tx_path)
            dataset_id <- sub("_tx_file\\.csv(\\.gz)?$", "", base)
        }
        suppressPackageStartupMessages({
            library(data.table)
            library(arrow)
        })
        if (verbose) message("[geneSCOPE::build_scope_from_cosmx] Reading ", tx_path)
        tx <- tryCatch({ data.table::fread(tx_path) }, error = function(e) stop("Failed to read ", tx_path, ": ", e$message))
        cn <- tolower(names(tx)); names(tx) <- cn
        headerless <- all(grepl("^v[0-9]+$", cn))
        if (!headerless && all(c("x_global_px", "y_global_px", "target") %in% cn)) {
            xg <- tx[["x_global_px"]]; yg <- tx[["y_global_px"]]; tgt <- tx[["target"]]
            qv <- if ("qv" %in% cn) tx[["qv"]] else Inf
            nd <- if ("nucleus_distance" %in% cn) tx[["nucleus_distance"]] else 0
        } else {
            if (ncol(tx) < 9L) stop("headerless tx_file has insufficient columns (<9)")
            if (verbose) message("[geneSCOPE::build_scope_from_cosmx] tx_file headerless – using positional mapping col6/col7/col9 (and col8 as nucleus_distance)")
            xg <- tx[[6]]; yg <- tx[[7]]; tgt <- tx[[9]]
            qv <- if (ncol(tx) >= 10 && is.numeric(tx[[10]])) tx[[10]] else Inf
            nd <- if (ncol(tx) >= 8 && is.numeric(tx[[8]])) tx[[8]] else 0
        }
        tx_out <- data.table::data.table(
            x_location = as.numeric(xg) * pixel_size_um,
            y_location = as.numeric(yg) * pixel_size_um,
            feature_name = as.character(tgt),
            qv = as.numeric(qv),
            nucleus_distance = as.numeric(nd)
        )
        work_dir <- file.path(tempdir(), paste0("cosmx_scope_", as.integer(Sys.time())))
        if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
        need_temp <- TRUE
        tf_parq <- file.path(work_dir, "transcripts.parquet")
        if (verbose) message("[geneSCOPE::build_scope_from_cosmx] Writing ", tf_parq)
        arrow::write_parquet(tx_out, tf_parq)
    } else if (verbose) {
        message("[geneSCOPE::build_scope_from_cosmx] Using existing transcripts.parquet")
    }

    have_cell_seg <- file.exists(seg_cell) || file.exists(seg_any)
    have_nuc_seg <- file.exists(nuc_parq)

    if (seg_type %in% c("cell", "both") && !have_cell_seg) {
        if (verbose) message("[geneSCOPE::build_scope_from_cosmx] Cell segmentation parquet not found; downgrade segmentation strategy")
        seg_type <- if (have_nuc_seg) "nucleus" else "none"
    }
    if (seg_type %in% c("nucleus", "both") && !have_nuc_seg) {
        if (verbose) message("[geneSCOPE::build_scope_from_cosmx] Nucleus segmentation parquet not found; downgrade segmentation strategy")
        seg_type <- if (have_cell_seg) "cell" else "none"
    }

    if (!file.exists(cell_parq) && build_from_flatfiles && build_cells_from_polygons) {
        poly_path <- NULL
        ff_dir <- file.path(input_dir, "flatFiles")
        cand <- c(Sys.glob(file.path(ff_dir, "*/*-polygons.csv.gz")), Sys.glob(file.path(ff_dir, "*/*_polygons.csv.gz")))
        if (length(cand)) poly_path <- cand[[1]]
        if (!is.null(poly_path)) {
            if (verbose) message("[geneSCOPE::build_scope_from_cosmx] cells.parquet not found; deriving centroids from polygons: ", basename(poly_path))
            pol <- data.table::fread(poly_path)
            cn <- tolower(names(pol)); names(pol) <- cn
            headerless <- all(grepl("^v[0-9]+$", cn))
            if (headerless) {
                if (ncol(pol) < 7L) stop("headerless polygons has insufficient columns (<7)")
                fov <- pol[[1]]; cid <- pol[[2]]; xg <- pol[[6]]; yg <- pol[[7]]
            } else {
                fx <- if ("fov" %in% cn) pol[["fov"]] else pol[[1]]
                cx <- if ("cellid" %in% cn) pol[["cellid"]] else if ("cell_id" %in% cn) pol[["cell_id"]] else pol[[2]]
                xg <- if ("x_global_px" %in% cn) pol[["x_global_px"]] else pol[[which.max(cn == "x")]]
                yg <- if ("y_global_px" %in% cn) pol[["y_global_px"]] else pol[[which.max(cn == "y")]]
                fov <- fx; cid <- cx
            }
            dt <- data.table::data.table(fov = as.integer(fov), cell_id = as.character(cid),
                x = as.numeric(xg) * pixel_size_um,
                y = as.numeric(yg) * pixel_size_um)
            cent <- dt[, .(x_centroid = mean(x, na.rm = TRUE),
                y_centroid = mean(y, na.rm = TRUE)), by = .(cell_id, fov)]
            if (!need_temp) {
                work_dir <- file.path(tempdir(), paste0("cosmx_scope_", as.integer(Sys.time())))
                dir.create(work_dir, recursive = TRUE)
                need_temp <- TRUE
            }
            cell_parq <- file.path(work_dir, "cells.parquet")
            if (verbose) message("[geneSCOPE::build_scope_from_cosmx] Writing ", cell_parq)
            arrow::write_parquet(cent, cell_parq)
        } else if (verbose) {
            message("[geneSCOPE::build_scope_from_cosmx] polygons not found; cannot derive cells.parquet")
        }
    }

    if (need_temp) {
        link_or_copy <- function(src, dst) {
            if (file.exists(src) && !file.exists(dst)) {
                ok <- FALSE
                suppressWarnings({ ok <- file.symlink(src, dst) })
                if (!ok) invisible(file.copy(src, dst, overwrite = FALSE))
            }
        }
        link_or_copy(file.path(input_dir, "segmentation_boundaries.parquet"), file.path(work_dir, "segmentation_boundaries.parquet"))
        link_or_copy(file.path(input_dir, "cell_boundaries.parquet"), file.path(work_dir, "cell_boundaries.parquet"))
        link_or_copy(file.path(input_dir, "nucleus_boundaries.parquet"), file.path(work_dir, "nucleus_boundaries.parquet"))
    }

    if (verbose) message("[geneSCOPE::build_scope_from_cosmx] Delegating to build_scope_from_xenium() …")
    scope_obj <- build_scope_from_xenium(
        input_dir = work_dir,
        grid_sizes = grid_length,
        segmentation_strategy = seg_type,
        data_type = data_type,
        verbose = verbose,
        ...
    )

    if (nrow(scope_obj@meta.data) > 0L) {
        if ("platform" %in% names(scope_obj@meta.data)) {
            scope_obj@meta.data$platform[] <- "CosMx"
        } else {
            scope_obj@meta.data$platform <- rep("CosMx", nrow(scope_obj@meta.data))
        }
        if ("dataset" %in% names(scope_obj@meta.data)) {
            scope_obj@meta.data$dataset[] <- if (is.null(dataset_id)) NA_character_ else dataset_id
        } else {
            scope_obj@meta.data$dataset <- rep(if (is.null(dataset_id)) NA_character_ else dataset_id,
                                               nrow(scope_obj@meta.data))
        }
    } else {
        scope_obj@meta.data <- data.frame(
            platform = "CosMx",
            dataset = if (is.null(dataset_id)) NA_character_ else dataset_id,
            stringsAsFactors = FALSE
        )
        rownames(scope_obj@meta.data) <- "__scope_platform__"
    }
    scope_obj@stats$platform <- "CosMx"
    scope_obj@stats$dataset <- if (is.null(dataset_id)) NA_character_ else dataset_id
    scope_obj
}

#' @describeIn build_scope_from_cosmx Backward-compatible wrapper retaining legacy argument names.
#' @export
createSCOPE_cosmx <- function(cosmx_root,
                              grid_length,
                              seg_type = c("cell", "nucleus", "both"),
                              dataset_id = NULL,
                              pixel_size_um = 0.120280945,
                              verbose = TRUE,
                              build_from_flatfiles = TRUE,
                              build_cells_from_polygons = TRUE,
                              data_type = "cosmx",
                              ...) {
    seg_type <- match.arg(seg_type)
    dots <- list(...)
    rename_map <- list(
        coord_file = "roi_path",
        ncores = "num_workers",
        chunk_pts = "chunk_size",
        min_seg_points = "min_segmentation_points",
        keep_partial_grid = "keep_partial_tiles",
        filtermolecule = "filter_molecules",
        exclude_prefix = "exclude_prefixes",
        filterqv = "min_quality",
        max_dist_mol_nuc = "max_distance_to_nucleus"
    )
    for (old in names(rename_map)) {
        if (!is.null(dots[[old]]) && is.null(dots[[rename_map[[old]]]])) {
            dots[[rename_map[[old]]]] <- dots[[old]]
        }
        dots[[old]] <- NULL
    }
    build_args <- c(
        list(
            input_dir = cosmx_root,
            grid_sizes = grid_length,
            segmentation_strategy = seg_type,
            dataset_id = dataset_id,
            pixel_size_um = pixel_size_um,
            verbose = verbose,
            allow_flatfile_generation = build_from_flatfiles,
            derive_cells_from_polygons = build_cells_from_polygons,
            data_type = data_type
        ),
        dots
    )
    do.call(build_scope_from_cosmx, build_args)
}

#' @title Build scope_object from 10x Visium spot-level outputs
#' @description
#'   Reads standard 10x Visium outputs (matrix + spatial metadata) and constructs
#'   a \code{scope_object} where the sole grid layer matches the Visium spot size.
#'   Coordinates are expressed in microns using \code{scalefactors_json.json},
#'   and only in-tissue spots are retained by default.
#' @param input_dir Character. Path to the Visium run (typically the \code{outs/} directory)
#'        containing \code{spatial/} and the feature-barcode matrix folder.
#' @param use_filtered_matrix Logical. If \code{TRUE} (default), prefer
#'        \code{filtered_feature_bc_matrix}; otherwise fall back to raw counts.
#' @param include_in_tissue_spots Logical. Keep only spots marked with
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
build_scope_from_visium <- function(input_dir,
                                    use_filtered_matrix = TRUE,
                                    include_in_tissue_spots = TRUE,
                                    grid_multiplier = 1,
                                    flip_y = FALSE,
                                    verbose = TRUE) {
    stopifnot(dir.exists(input_dir))

    visium_dir <- input_dir
    use_filtered <- use_filtered_matrix
    include_in_tissue <- include_in_tissue_spots

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
    if (verbose) message("[geneSCOPE::build_scope_from_visium] Using matrix directory: ", basename(matrix_dir))

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
    if (verbose) message("[geneSCOPE::build_scope_from_visium] Loaded ", nrow(pos_raw), " spot positions")

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
        if (verbose) message("[geneSCOPE::build_scope_from_visium] Retained ", nrow(pos_dt), " in-tissue spots")
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

    if (is.null(spot_diameter_um) && !is.null(spot_diameter_px) && !is.null(microns_per_pixel)) {
        spot_diameter_um <- as.numeric(spot_diameter_px) * as.numeric(microns_per_pixel)
    }
    if (is.null(microns_per_pixel) && !is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
        microns_per_pixel <- as.numeric(spot_diameter_um) / as.numeric(spot_diameter_px)
    }
    if (is.null(microns_per_pixel) && is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
        default_spot_um <- getOption("geneSCOPE.default_visium_spot_um", 55)
        microns_per_pixel <- as.numeric(default_spot_um) / as.numeric(spot_diameter_px)
        spot_diameter_um <- as.numeric(default_spot_um)
    }
    if (is.null(spot_diameter_um)) stop("spot_diameter_microns is missing from scalefactors_json.json.")
    if (is.null(microns_per_pixel)) stop("microns_per_pixel is missing from scalefactors_json.json.")

    spot_diameter_um <- as.numeric(spot_diameter_um) * grid_multiplier
    microns_per_pixel <- as.numeric(microns_per_pixel)

    pos_dt[, `:=`(
        x = pxl_col_in_fullres * microns_per_pixel,
        y = pxl_row_in_fullres * microns_per_pixel
    )]
    if (isTRUE(flip_y)) {
        ymax <- max(pos_dt$y, na.rm = TRUE)
        pos_dt[, y := ymax - y]
    }

    scope_obj <- new("scope_object",
        coord = list(centroids = pos_dt[, .(cell = barcode, x, y)]),
        grid = list(),
        meta.data = data.frame(row.names = pos_dt$barcode)
    )
    scope_obj@meta.data$platform <- "Visium"
    scope_obj@stats$platform <- "Visium"
    scope_obj@stats$visium <- list(
        counts = counts_mat,
        genes = gene_meta,
        microns_per_pixel = microns_per_pixel,
        spot_diameter_um = spot_diameter_um,
        include_in_tissue = include_in_tissue,
        grid_multiplier = grid_multiplier
    )
    scope_obj
}

#' @describeIn build_scope_from_visium Backward-compatible wrapper retaining legacy argument names.
#' @export
createSCOPE_visium <- function(visium_dir,
                               use_filtered = TRUE,
                               include_in_tissue = TRUE,
                               grid_multiplier = 1,
                               flip_y = FALSE,
                               verbose = TRUE) {
    build_scope_from_visium(
        input_dir = visium_dir,
        use_filtered_matrix = use_filtered,
        include_in_tissue_spots = include_in_tissue,
        grid_multiplier = grid_multiplier,
        flip_y = flip_y,
        verbose = verbose
    )
}

# ---------------------------------------------------------------------------
# Visium via Seurat (SCTransform optional)
# ---------------------------------------------------------------------------
 
#' Build a scope_object via Seurat Visium workflow.
#'
#' Uses Seurat's Visium loader (optionally running SCTransform) and converts the
#' resulting object into a `scope_object`.
#'
#' @param input_dir Path to the Visium run directory.
#' @param run_sctransform Logical; whether to run SCTransform.
#' @param sct_return_only_var_genes Logical; passed to SCTransform.
#' @param vars_to_regress Variables to regress during SCTransform.
#' @param glmGamPoi Logical; use glmGamPoi for SCTransform.
#' @param include_in_tissue_spots Logical; keep only in-tissue spots.
#' @param grid_multiplier Grid multiplier applied to spot diameter.
#' @param flip_y Logical flag to flip Y axis.
#' @param verbose Logical toggle for progress output.
#' @param data_type Character label stored on the output object.
#' @return A list containing the `scope_object` and the Seurat object when SCTransform is run.
#' @keywords internal
build_scope_from_visium_seurat <- function(input_dir,
                                          run_sctransform = TRUE,
                                          sct_return_only_var_genes = FALSE,
                                          vars_to_regress = NULL,
                                          glmGamPoi = TRUE,
                                          include_in_tissue_spots = TRUE,
                                          grid_multiplier = 1,
                                          flip_y = FALSE,
                                          verbose = TRUE,
                                          data_type = "visium") {
    stopifnot(dir.exists(input_dir))
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' is required for build_scope_from_visium_seurat(). Please install it.")
    }
    if (!is.numeric(grid_multiplier) || length(grid_multiplier) != 1L || !is.finite(grid_multiplier) || grid_multiplier <= 0) {
        stop("grid_multiplier must be a single positive finite numeric value.")
    }

    log_prefix <- sprintf("[geneSCOPE::build_scope/%s]", data_type)

    suppressPackageStartupMessages({
        library(Matrix)
        library(data.table)
    })

    if (isTRUE(verbose)) message(log_prefix, " Calling Seurat::Load10X_Spatial() …")
    h5 <- file.path(input_dir, "filtered_feature_bc_matrix.h5")
    has_h5 <- file.exists(h5)
    seu <- tryCatch({
        if (has_h5) {
            Seurat::Load10X_Spatial(data.dir = input_dir, filename = basename(h5))
        } else {
            Seurat::Load10X_Spatial(data.dir = input_dir)
        }
    }, error = function(e) {
        stop("Load10X_Spatial failed: ", conditionMessage(e))
    })

    if (isTRUE(run_sctransform)) {
        if (isTRUE(verbose)) message(log_prefix, " Calling Seurat::SCTransform() …")
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
        sct_mat <- getAssayLayer(seu, assay = "SCT", which = "scale.data")
    }

    spatial_dir <- file.path(input_dir, "spatial")
    if (!dir.exists(spatial_dir)) stop("spatial/ directory not found under ", input_dir)

    pos_candidates <- c(
        file.path(spatial_dir, "tissue_positions_list.csv"),
        file.path(spatial_dir, "tissue_positions.csv")
    )
    pos_path <- pos_candidates[file.exists(pos_candidates)][1]
    if (is.na(pos_path)) stop("Could not find tissue_positions_list.csv or tissue_positions.csv under ", spatial_dir)

    pos_raw <- data.table::fread(pos_path)
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
    if (isTRUE(include_in_tissue_spots)) pos_dt <- pos_dt[in_tissue == 1L]
    keep <- pos_dt$barcode %in% colnames(raw_counts)
    if (!all(keep)) pos_dt <- pos_dt[keep]
    if (nrow(pos_dt) == 0L) stop("No overlap between positions and Seurat barcodes.")

    scalefactor_candidates <- c(
        file.path(spatial_dir, "scalefactors_json.json"),
        file.path(spatial_dir, "scalefactors_json_fullres.json")
    )
    sf_path <- scalefactor_candidates[file.exists(scalefactor_candidates)][1]
    if (is.na(sf_path)) stop("scalefactors_json.json not found under ", spatial_dir)
    sf <- jsonlite::fromJSON(sf_path)
    microns_per_pixel <- sf$microns_per_pixel
    spot_diameter_um <- sf$spot_diameter_microns
    spot_diameter_px <- sf$spot_diameter_fullres
    hires_scalef <- sf$tissue_hires_scalef %||% sf$tissue_hires_scale %||% NA_real_
    lowres_scalef <- sf$tissue_lowres_scalef %||% sf$tissue_lowres_scale %||% NA_real_
    if (is.null(spot_diameter_um) && !is.null(spot_diameter_px) && !is.null(microns_per_pixel)) {
        spot_diameter_um <- as.numeric(spot_diameter_px) * as.numeric(microns_per_pixel)
    }
    if (is.null(microns_per_pixel) && !is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
        microns_per_pixel <- as.numeric(spot_diameter_um) / as.numeric(spot_diameter_px)
    }
    if (is.null(microns_per_pixel) && is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
        default_spot_um <- getOption("geneSCOPE.default_visium_spot_um", 55)
        microns_per_pixel <- as.numeric(default_spot_um) / as.numeric(spot_diameter_px)
        spot_diameter_um <- as.numeric(default_spot_um)
        if (isTRUE(verbose)) message(log_prefix, " microns_per_pixel missing; assumed spot_diameter_microns = ", default_spot_um,
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

    pos_dt[, `:=`(
        x_um = pxl_col_in_fullres * microns_per_pixel,
        y_um = pxl_row_in_fullres * microns_per_pixel
    )]
    if (isTRUE(flip_y)) {
        y_max <- max(pos_dt$y_um, na.rm = TRUE)
        pos_dt[, y_um := y_max - y_um]
    }

    candidate_dirs <- unique(file.path(input_dir, c("filtered_feature_bc_matrix", "raw_feature_bc_matrix")))
    matrix_dir <- candidate_dirs[dir.exists(candidate_dirs)][1]
    feature_file <- if (!is.na(matrix_dir)) file.path(matrix_dir, "features.tsv.gz") else NA
    if (!is.na(feature_file) && !file.exists(feature_file)) feature_file <- file.path(matrix_dir, "features.tsv")
    if (!is.na(feature_file) && file.exists(feature_file)) {
        features_dt <- data.table::fread(feature_file, header = FALSE)
        gene_ids <- as.character(features_dt[[1]])
        gene_names <- make.unique(as.character(features_dt[[2]]))
        feature_type <- if (ncol(features_dt) >= 3) as.character(features_dt[[3]]) else rep(NA_character_, length(gene_names))
        idx <- match(gene_names, rownames(raw_counts))
        ok <- !is.na(idx)
        gene_meta <- data.frame(
            feature_id = gene_ids[ok],
            feature_type = feature_type[ok],
            platform = rep(data_type, sum(ok)),
            stringsAsFactors = FALSE,
            row.names = gene_names[ok]
        )
    } else {
        gene_meta <- data.frame(
            feature_id = rownames(raw_counts),
            feature_type = NA_character_,
            platform = rep(data_type, nrow(raw_counts)),
            stringsAsFactors = FALSE,
            row.names = rownames(raw_counts)
        )
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

    col_idx <- match(pos_dt$barcode, colnames(raw_counts))
    raw_counts <- raw_counts[, col_idx, drop = FALSE]
    colnames(raw_counts) <- pos_dt$barcode

    if (!is.null(sct_mat)) {
        col_idx2 <- match(pos_dt$barcode, colnames(sct_mat))
        sct_mat <- sct_mat[, col_idx2, drop = FALSE]
        colnames(sct_mat) <- pos_dt$barcode
    }

    counts_tbl <- Matrix::summary(raw_counts)
    counts_dt <- data.table::data.table(
        grid_id = colnames(raw_counts)[counts_tbl$j],
        gene = rownames(raw_counts)[counts_tbl$i],
        count = as.numeric(counts_tbl$x)
    )
    counts_dt <- counts_dt[count > 0]

    fmt_len <- formatC(grid_size_um, format = "fg", digits = 6)
    fmt_len <- trimws(fmt_len)
    fmt_len <- sub("\\.?0+$", "", fmt_len)
    grid_name <- paste0("grid", fmt_len)

    infer_hex_orientation <- function(cent, tol_deg = 15) {
        cn <- c("x", "y", "array_row", "array_col")
        if (!all(cn %in% names(cent))) return(list(orientation = NA_character_, horiz_frac = NA_real_, vert_frac = NA_real_, tol_deg = tol_deg))
        dt <- data.table::as.data.table(cent)[, .(x, y, array_row = as.integer(array_row), array_col = as.integer(array_col))]
        get_pairs <- function(dr, dc) {
            a <- dt
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
        vert <- abs(theta - 90) <= tol_deg
        hfrac <- mean(horiz, na.rm = TRUE)
        vfrac <- mean(vert, na.rm = TRUE)
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
        visium_orientation_vert_frac = ori$vert_frac,
        visium_orientation_tol_deg = ori$tol_deg,
        histology = list(
            lowres = NULL,
            hires = NULL
        )
    )
    roi_bbox <- c(
        xmin = min(grid_dt$xmin, na.rm = TRUE),
        xmax = max(grid_dt$xmax, na.rm = TRUE),
        ymin = min(grid_dt$ymin, na.rm = TRUE),
        ymax = max(grid_dt$ymax, na.rm = TRUE)
    )

    image_candidates_hires <- c(
        file.path(spatial_dir, "aligned_tissue_image.jpg"),
        file.path(spatial_dir, "tissue_hires_image.png"),
        file.path(spatial_dir, "tissue_hires_image.jpg")
    )
    image_candidates_lowres <- c(
        file.path(spatial_dir, "tissue_lowres_image.png"),
        file.path(spatial_dir, "detected_tissue_image.jpg")
    )
    hires_path <- image_candidates_hires[file.exists(image_candidates_hires)][1]
    lowres_path <- image_candidates_lowres[file.exists(image_candidates_lowres)][1]
    get_img_dim <- function(p) {
        if (is.na(p)) return(c(NA_integer_, NA_integer_))
        ext <- tolower(tools::file_ext(p))
        dims <- c(NA_integer_, NA_integer_)
        try({
            if (ext %in% c("png")) {
                if (requireNamespace("png", quietly = TRUE)) {
                    a <- png::readPNG(p)
                    dims <- c(dim(a)[2], dim(a)[1])
                }
            } else if (ext %in% c("jpg", "jpeg")) {
                if (requireNamespace("jpeg", quietly = TRUE)) {
                    a <- jpeg::readJPEG(p)
                    dims <- c(dim(a)[2], dim(a)[1])
                }
            }
        }, silent = TRUE)
        dims
    }
    hires_dim <- get_img_dim(hires_path)
    lowres_dim <- get_img_dim(lowres_path)
    scope_obj@grid[[grid_name]]$image_info <- list(
        hires_path = if (is.na(hires_path)) NULL else hires_path,
        lowres_path = if (is.na(lowres_path)) NULL else lowres_path,
        hires_dim_px = hires_dim,
        lowres_dim_px = lowres_dim,
        positions_path = pos_path,
        scalefactors_path = sf_path,
        tissue_hires_scalef = hires_scalef,
        tissue_lowres_scalef = lowres_scalef,
        spot_diameter_fullres = spot_diameter_px,
        y_origin = "top-left"
    )

    sf_lookup <- list(
        hires = suppressWarnings(as.numeric(hires_scalef)),
        lowres = suppressWarnings(as.numeric(lowres_scalef))
    )
    y_origin_img <- scope_obj@grid[[grid_name]]$image_info$y_origin
    if (is.null(y_origin_img)) {
        y_origin_img <- "top-left"
    }
    safe_attach_histology <- function(scope_obj,
                                      level,
                                      path) {
        if (is.null(path) || is.na(path)) return(scope_obj)
        tryCatch(
            attach_histology(
                scope_obj = scope_obj,
                grid_name = grid_name,
                png_path = path,
                json_path = sf_path,
                level = level,
                roi_bbox = roi_bbox,
                scalefactors = sf_lookup,
                y_origin = y_origin_img
            ),
            error = function(e) {
                if (isTRUE(verbose)) {
                    message(log_prefix, " Histology attach (", level, ") skipped: ", conditionMessage(e))
                }
                scope_obj
            }
        )
    }
    scope_obj <- safe_attach_histology(scope_obj, "hires", if (is.na(hires_path)) NULL else hires_path)
    scope_obj <- safe_attach_histology(scope_obj, "lowres", if (is.na(lowres_path)) NULL else lowres_path)

    if (!is.null(sct_mat)) {
        scope_obj@grid[[grid_name]]$SCT <- t(sct_mat)
    }

    scope_obj@stats$platform <- data_type
    scope_obj@meta.data$platform <- data_type
    if (isTRUE(verbose)) {
        message(log_prefix, " scope_object created with ", nrow(centroids_dt),
            " spots and grid_length = ", round(grid_size_um, 3), " um")
    }

    scope_obj
}

#' Create a scope_object from heterogeneous platforms.
#'
#' Auto-detects the dataset type (Xenium, CosMx, Visium) or honours the user
#' preference, forwarding the request to the appropriate builder.
#'
#' @param data_dir Path to the dataset root (outs directory for Xenium/Visium,
#'        project root for CosMx).
#' @param prefer Which backend to prefer (`"auto"`, `"xenium"`, `"cosmx"`, `"visium"`).
#' @param verbose Logical toggle for progress messages.
#' @param sctransform Logical enabling SCTransform when routing through the Visium Seurat helper.
#' @param return_both Logical; for Visium only, return `list(scope_obj, seurat)` when `TRUE`.
#' @param ... Additional parameters forwarded to the chosen builder. Xenium/CosMx
#'        require `grid_length`, Visium does not.
#' @return A `scope_object` (or a list when `return_both = TRUE` for Visium SCTransform workflows).
#' @export
createSCOPE <- function(data_dir = NULL,
                        prefer = c("auto", "xenium", "cosmx", "visium"),
                        verbose = TRUE,
                        sctransform = TRUE,
                        ...) {
    dots <- list(...)
    if (is.null(dots$verbose)) dots$verbose <- verbose
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

    auto_pick <- function(d) {
        if (is_visium(d)) return("visium")
        if (is_xenium(d)) return("xenium")
        if (is_cosmx(d)) return("cosmx")
        stop("Unable to detect data type under ", d,
            "; expected Xenium (parquet), Visium (spatial/ + matrix), or CosMx (flatFiles/).")
    }

    pick <- switch(prefer,
        xenium = "xenium",
        visium = "visium",
        cosmx = "cosmx",
        auto = auto_pick(data_dir)
    )

    dots$xenium_dir <- NULL; dots$visium_dir <- NULL; dots$cosmx_root <- NULL

    if (pick == "xenium") {
        if (isTRUE(verbose)) message("[geneSCOPE::createSCOPE/xenium] Dispatching to build_scope_from_xenium() …")
        if (is.null(dots$grid_length)) stop("grid_length is required for Xenium. Pass grid_length= ...")
        dots$xenium_dir <- data_dir
        return(do.call(createSCOPE_xenium, dots))
    }
    if (pick == "cosmx") {
        if (isTRUE(verbose)) message("[geneSCOPE::createSCOPE/cosmx] Dispatching to build_scope_from_cosmx() …")
        if (is.null(dots$grid_length)) stop("grid_length is required for CosMx. Pass grid_length= ...")
        dots$cosmx_root <- data_dir
        return(do.call(createSCOPE_cosmx, dots))
    }

    dots$visium_dir <- data_dir
    if (isTRUE(sctransform)) {
        if (is.null(dots$run_sctransform)) dots$run_sctransform <- TRUE
        if (isTRUE(verbose)) {
            message(sprintf("[geneSCOPE::createSCOPE visium] Dispatching to Seurat pipeline (SCTransform=%s) …",
                if (is.null(dots$run_sctransform)) "NULL" else as.character(dots$run_sctransform)))
        }
        rename_map <- list(
            visium_dir = "input_dir",
            include_in_tissue = "include_in_tissue_spots",
            grid_multiplier = "grid_multiplier",
            flip_y = "flip_y",
            SCTransform = "run_sctransform",
            sctransform = "run_sctransform"
        )
        for (old in names(rename_map)) {
            if (!is.null(dots[[old]]) && is.null(dots[[rename_map[[old]]]])) {
                dots[[rename_map[[old]]]] <- dots[[old]]
            }
            if (old != "grid_multiplier" && old != "flip_y") dots[[old]] <- NULL
        }
        if (!is.null(dots$visium_dir) && is.null(dots$input_dir)) dots$input_dir <- dots$visium_dir
        dots$visium_dir <- NULL
        res <- do.call(build_scope_from_visium_seurat, dots)
        res
    } else {
        if (isTRUE(verbose)) message("[geneSCOPE::createSCOPE visium] Dispatching to createSCOPE_visium (SCTransform=FALSE) …")
        obj <- do.call(createSCOPE_visium, dots)
        obj
    }
}
#' Xenium single-cell ingestion helper.
#'
#' Internal helper that powers `addSingleCells()` when the backend is Xenium.
#' @noRd
.add_singlecells_xenium <- function(scope_obj,
                                    xenium_dir,
                                    filter_genes = NULL,
                                    exclude_prefix = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword"),
                                    filter_barcodes = TRUE,
                                    verbose = TRUE) {
    cfg <- .build_singlecell_config(
        scope_obj = scope_obj,
        platform = "Xenium",
        xenium_dir = xenium_dir,
        cosmx_root = NULL,
        filter_genes = filter_genes,
        exclude_prefix = exclude_prefix,
        filter_barcodes = filter_barcodes,
        id_mode = NULL,
        force_all_genes = TRUE,
        verbose = verbose
    )
    .run_singlecell_pipeline(cfg)
}

#' CosMx single-cell ingestion helper.
#'
#' Internal helper that powers `addSingleCells()` when the backend is CosMx.
#' @noRd
.add_singlecells_cosmx <- function(scope_obj,
                                   cosmx_root,
                                   filter_genes = NULL,
                                   exclude_prefix = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword",
                                                      "SystemControl", "Negative"),
                                   id_mode = c("cell_id", "fov_cell", "auto"),
                                   filter_barcodes = TRUE,
                                   force_all_genes = FALSE,
                                   verbose = TRUE) {
    cfg <- .build_singlecell_config(
        scope_obj = scope_obj,
        platform = "CosMx",
        xenium_dir = NULL,
        cosmx_root = cosmx_root,
        filter_genes = filter_genes,
        exclude_prefix = exclude_prefix,
        filter_barcodes = filter_barcodes,
        id_mode = match.arg(id_mode),
        force_all_genes = force_all_genes,
        verbose = verbose
    )
    .run_singlecell_pipeline(cfg)
}

#' @title Add Single Cells (auto-dispatch)
#' @description
#'   Wrapper that dispatches to the Xenium or CosMx ingestion helper based on
#'   the recorded platform tag in `scope_obj` (`meta.data$platform` or
#'   `stats$platform`). You may also force the backend by providing `platform`,
#'   or by supplying a generic `data_dir`/`input_dir` that mirrors the
#'   the `createSCOPE()` entry point parameters.
#' @param scope_obj A `scope_object` with centroids loaded.
#' @param data_dir Optional generic directory pointing to single-cell outputs
#'        (alias for `xenium_dir`/`cosmx_root`; when provided the function will
#'        auto-detect the platform).
#' @param input_dir Alias for `data_dir`.
#' @param xenium_dir Optional. Xenium outs directory.
#' @param cosmx_root Optional. CosMx project root containing `flatFiles/`.
#' @param platform Character. One of `"auto"`, `"Xenium"`, `"CosMx"`.
#'        Default `"auto"`.
#' @param verbose Logical. Print progress messages.
#' @param ... Additional arguments forwarded to the chosen backend (e.g.
#'        `filter_genes`, `exclude_prefix`, `filter_barcodes`, `force_all_genes`).
#' @return The modified `scope_object` with `@cells$counts` added.
#' @importFrom Matrix sparseMatrix t
#' @importFrom data.table fread setDT
#' @export
addSingleCells <- function(scope_obj,
                           data_dir = NULL,
                           input_dir = NULL,
                           xenium_dir = NULL,
                           cosmx_root = NULL,
                           platform = c("auto", "Xenium", "CosMx"),
                           verbose = TRUE,
                           ...) {
    stopifnot(inherits(scope_obj, "scope_object"))
    platform <- match.arg(platform)
    dots <- list(...)

    if (!is.null(input_dir) && is.null(data_dir)) data_dir <- input_dir
    if (!is.null(dots$input_dir) && is.null(data_dir)) {
        data_dir <- dots$input_dir
    }
    dots$input_dir <- NULL

    if (!is.null(dots$data_dir) && is.null(data_dir)) {
        data_dir <- dots$data_dir
    }
    dots$data_dir <- NULL

    if (!is.null(data_dir) && !dir.exists(data_dir)) {
        stop("`data_dir` does not exist: ", data_dir)
    }

    if (is.null(xenium_dir) && !is.null(dots$xenium_dir)) {
        xenium_dir <- dots$xenium_dir
    }
    dots$xenium_dir <- NULL

    if (is.null(cosmx_root) && !is.null(dots$cosmx_root)) {
        cosmx_root <- dots$cosmx_root
    }
    dots$cosmx_root <- NULL

    if (!is.null(data_dir) && is.null(xenium_dir)) {
        guessed <- .guess_singlecell_platform(data_dir)
        if (identical(guessed$platform, "Xenium") && is.null(xenium_dir)) {
            xenium_dir <- guessed$path
        } else if (identical(guessed$platform, "CosMx") && is.null(cosmx_root)) {
            cosmx_root <- guessed$path
        }
        if (is.null(guessed$platform)) {
            # fall back to using the provided data_dir directly when platform is forced later
            if (is.null(xenium_dir)) xenium_dir <- data_dir
            if (is.null(cosmx_root)) cosmx_root <- data_dir
        }
    }

    chosen <- NULL
    if (!is.null(xenium_dir)) chosen <- "Xenium"
    if (!is.null(cosmx_root)) chosen <- "CosMx"

    if (is.null(chosen) && platform == "auto") {
        plat <- NA_character_
        if (!is.null(scope_obj@meta.data) && nrow(scope_obj@meta.data)) {
            if ("platform" %in% names(scope_obj@meta.data)) {
                vals <- as.character(scope_obj@meta.data$platform)
                vals <- vals[!is.na(vals)]
                if (length(vals)) {
                    freq <- sort(table(vals), decreasing = TRUE)
                    plat <- names(freq)[1]
                }
            }
        }
        if (is.na(plat) && !is.null(scope_obj@stats$platform)) {
            plat <- as.character(scope_obj@stats$platform)
        }
        if (!is.na(plat)) {
            if (grepl("cosmx", plat, ignore.case = TRUE)) chosen <- "CosMx"
            else if (grepl("xenium", plat, ignore.case = TRUE)) chosen <- "Xenium"
        }
    }
    if (!identical(platform, "auto")) {
        chosen <- platform
    }

    if (is.null(chosen)) {
        if (!is.null(xenium_dir)) chosen <- "Xenium"
        else if (!is.null(cosmx_root)) chosen <- "CosMx"
        else chosen <- "Xenium"
    }

    if (identical(chosen, "CosMx")) {
        if (is.null(cosmx_root)) {
            if (!is.null(data_dir)) cosmx_root <- data_dir
            if (is.null(cosmx_root)) stop("CosMx backend selected but `cosmx_root` (or data_dir) is missing.")
        }
        if (isTRUE(verbose)) message("[geneSCOPE::addSingleCells cosmx] Dispatching to addSingleCells_cosmx() …")
        root <- if (!is.null(cosmx_root)) cosmx_root else getOption("geneSCOPE.cosmx_root")
        helper_args <- dots
        helper_args$scope_obj <- scope_obj
        helper_args$cosmx_root <- root
        if (is.null(helper_args$verbose)) helper_args$verbose <- verbose
        return(do.call(.add_singlecells_cosmx, helper_args))
    }

    if (is.null(xenium_dir)) {
        if (!is.null(data_dir)) xenium_dir <- data_dir
        if (is.null(xenium_dir)) stop("Xenium backend selected but `xenium_dir` (or data_dir) is missing.")
    }
    if (isTRUE(verbose)) message("[geneSCOPE::addSingleCells xenium] Dispatching to addSingleCells_xenium() …")
    dir_path <- if (!is.null(xenium_dir)) xenium_dir else getOption("geneSCOPE.xenium_dir")
    helper_args <- dots
    helper_args$scope_obj <- scope_obj
    helper_args$xenium_dir <- dir_path
    if (is.null(helper_args$verbose)) helper_args$verbose <- verbose
    do.call(.add_singlecells_xenium, helper_args)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Dispatch the single-cell ingestion pipeline.
#'
#' Combined entrypoint that decides whether to run the Xenium or CosMx helper
#' according to the platform recorded in `cfg`.
#'
#' @param cfg Configuration list produced by `.build_singlecell_config()`.
#' @keywords internal
#' @noRd
.run_singlecell_pipeline <- function(cfg) {
    if (identical(cfg$platform, "Xenium")) {
        return(.run_xenium_pipeline(cfg))
    }
    if (identical(cfg$platform, "CosMx")) {
        return(.run_cosmx_pipeline(cfg))
    }
    stop("Unsupported platform: ", cfg$platform)
}

#' Construct a normalized configuration for single-cell ingestion.
#'
#' Collapses user-facing arguments (and aliases) into a single configuration
#' list consumed by downstream helpers.
#'
#' @param scope_obj A `scope_object`.
#' @param platform Character string (`"Xenium"` or `"CosMx"`).
#' @param xenium_dir Xenium outs directory (optional).
#' @param cosmx_root CosMx project root (optional).
#' @param filter_genes Optional gene allowlist.
#' @param exclude_prefix Vector of prefixes to remove.
#' @param filter_barcodes Logical toggling centroid filtering.
#' @param id_mode Identifier mode for CosMx ingestion.
#' @param force_all_genes Logical forcing CosMx to load the full panel.
#' @param verbose Logical flag controlling messages.
#' @keywords internal
#' @noRd
.build_singlecell_config <- function(scope_obj,
                                     platform,
                                     xenium_dir,
                                     cosmx_root,
                                     filter_genes,
                                     exclude_prefix,
                                     filter_barcodes,
                                     id_mode,
                                     force_all_genes,
                                     verbose) {
    stopifnot(inherits(scope_obj, "scope_object"))
    platform <- match.arg(platform, c("Xenium", "CosMx"))
    verbose <- isTRUE(verbose)
    filter_cells <- isTRUE(filter_barcodes)

    if (identical(platform, "Xenium")) {
        if (is.null(xenium_dir) || !dir.exists(xenium_dir)) {
            stop("`xenium_dir` must point to a valid Xenium outs/ directory.")
        }
    } else {
        if (is.null(cosmx_root) || !dir.exists(cosmx_root)) {
            stop("`cosmx_root` must point to a valid CosMx project directory.")
        }
        if (is.null(id_mode)) id_mode <- "auto"
        id_mode <- match.arg(id_mode, c("cell_id", "fov_cell", "auto"))
    }

    prefix <- switch(platform,
        Xenium = "[geneSCOPE::addSingleCells xenium]",
        CosMx = "[geneSCOPE::addSingleCells cosmx]"
    )

    centroid_cells <- NULL
    if (filter_cells) {
        centroid_cells <- scope_obj@coord$centroids$cell
        if (!length(centroid_cells)) {
            stop(prefix, " scope_obj@coord$centroids is empty, cannot determine cells to keep.")
        }
        if (verbose) {
            message(prefix, " Target cells: ", length(centroid_cells))
        }
    } else if (verbose) {
        message(prefix, " Barcode filtering disabled - will include all cells from source matrix")
    }

    list(
        scope_obj = scope_obj,
        platform = platform,
        verbose = verbose,
        prefix = prefix,
        paths = list(
            xenium = xenium_dir,
            cosmx = cosmx_root
        ),
        filters = list(
            gene_allowlist = filter_genes,
            gene_exclude_prefix = exclude_prefix,
            filter_cells = filter_cells
        ),
        cosmx = list(
            id_mode = if (identical(platform, "CosMx")) id_mode else NULL,
            force_all_genes = isTRUE(force_all_genes)
        ),
        targets = centroid_cells
    )
}

#' Guess the single-cell platform from a directory layout.
#'
#' @param path Directory path to inspect.
#' @return List with detected `platform` and normalised `path`.
#' @keywords internal
#' @noRd
.guess_singlecell_platform <- function(path) {
    if (is.null(path) || !dir.exists(path)) {
        return(list(platform = NULL, path = NULL))
    }
    candidates <- c(path, file.path(path, "outs"))
    for (cand in candidates) {
        if (dir.exists(cand) && file.exists(file.path(cand, "cell_feature_matrix.h5"))) {
            return(list(platform = "Xenium", path = cand))
        }
    }
    if (dir.exists(file.path(path, "flatFiles"))) {
        return(list(platform = "CosMx", path = path))
    }
    list(platform = NULL, path = NULL)
}

#' Execute the Xenium ingestion workflow.
#'
#' Loads the sparse Xenium matrix, applies filtering logic and writes data back
#' into the supplied `scope_object`.
#'
#' @param cfg Configuration list from `.build_singlecell_config()`.
#' @keywords internal
#' @noRd
.run_xenium_pipeline <- function(cfg) {
    prefix <- cfg$prefix
    verbose <- cfg$verbose
    if (verbose) message(prefix, " Loading cell-feature matrix from HDF5")

    h5_file <- .resolve_xenium_matrix_path(cfg$paths$xenium)
    components <- .read_xenium_sparse_components(h5_file)
    matrix_bundle <- .construct_xenium_sparse_matrix(components, prefix, verbose)

    filtered <- .filter_gene_panel(
        counts = matrix_bundle$counts,
        filters = cfg$filters,
        verbose = verbose,
        prefix = prefix
    )
    aligned <- .align_counts_to_targets(
        counts = filtered$counts,
        targets = cfg$targets,
        filter_cells = cfg$filters$filter_cells,
        verbose = verbose,
        prefix = prefix
    )
    if (!is.null(matrix_bundle$gene_map)) {
        attr(aligned$counts, "gene_map") <- matrix_bundle$gene_map
    }
    .store_scope_counts(cfg$scope_obj, aligned$counts, prefix, verbose)
}

#' Execute the CosMx ingestion workflow.
#'
#' Loads the expression matrix, aligns identifiers and stores the assembled
#' sparse matrix back into the `scope_object`.
#'
#' @param cfg Configuration list from `.build_singlecell_config()`.
#' @keywords internal
#' @noRd
.run_cosmx_pipeline <- function(cfg) {
    prefix <- cfg$prefix
    verbose <- cfg$verbose

    dataset_info <- .locate_cosmx_dataset(cfg$paths$cosmx)
    expr_mat <- dataset_info$expr_mat
    if (verbose) message(prefix, " Expression matrix: ", basename(expr_mat))

    header_cols <- .read_cosmx_header(expr_mat)
    col_info <- .identify_cosmx_gene_columns(header_cols, cfg$filters$gene_exclude_prefix)
    gene_panel <- .resolve_cosmx_gene_panel(
        filter_genes = cfg$filters$gene_allowlist,
        auto_genes = col_info$gene_cols,
        force_all_genes = cfg$cosmx$force_all_genes,
        verbose = verbose,
        prefix = prefix
    )

    selected_cols <- c(col_info$fov_col, col_info$cid_col, gene_panel)
    expr_dt <- .read_cosmx_expression_subset(expr_mat, selected_cols, verbose, prefix)

    keyed <- .assign_cosmx_cell_keys(
        expr_dt = expr_dt,
        cid_col = col_info$cid_col,
        fov_col = col_info$fov_col,
        id_mode = cfg$cosmx$id_mode,
        centroid_keys = cfg$targets,
        filter_cells = cfg$filters$filter_cells,
        verbose = verbose,
        prefix = prefix
    )

    counts <- .assemble_cosmx_sparse_matrix(
        expr_dt = keyed$expr_dt,
        gene_cols = gene_panel,
        column_order = keyed$column_order
    )

    .store_scope_counts(cfg$scope_obj, counts, prefix, verbose)
}

#' Locate the Xenium cell-feature matrix file.
#'
#' @param xenium_dir Base directory containing Xenium outputs.
#' @return Absolute path to `cell_feature_matrix.h5`.
#' @keywords internal
#' @noRd
.resolve_xenium_matrix_path <- function(xenium_dir) {
    h5_file <- file.path(xenium_dir, "cell_feature_matrix.h5")
    if (!file.exists(h5_file)) {
        stop("HDF5 file not found: ", h5_file)
    }
    h5_file
}

#' Read Xenium sparse components from HDF5.
#'
#' @param h5_file Path to the Xenium HDF5 matrix.
#' @return List containing barcodes, indices, data and gene mapping.
#' @keywords internal
#' @noRd
.read_xenium_sparse_components <- function(h5_file) {
    if (requireNamespace("rhdf5", quietly = TRUE)) {
        gene_id <- rhdf5::h5read(h5_file, "matrix/features/id")
        gene_name <- rhdf5::h5read(h5_file, "matrix/features/name")
        list(
            barcodes = rhdf5::h5read(h5_file, "matrix/barcodes"),
            gene_id = gene_id,
            gene_name = gene_name,
            data = rhdf5::h5read(h5_file, "matrix/data"),
            indices = rhdf5::h5read(h5_file, "matrix/indices"),
            indptr = rhdf5::h5read(h5_file, "matrix/indptr"),
            shape = as.integer(rhdf5::h5read(h5_file, "matrix/shape")),
            gene_map = data.frame(
                ensembl = gene_id,
                symbol = gene_name,
                stringsAsFactors = FALSE
            ),
            backend = "rhdf5"
        )
    } else if (requireNamespace("hdf5r", quietly = TRUE)) {
        h5file <- hdf5r::H5File$new(h5_file, mode = "r")
        on.exit(h5file$close(), add = TRUE)
        gene_id <- h5file[["matrix/features/id"]][]
        gene_name <- h5file[["matrix/features/name"]][]
        list(
            barcodes = h5file[["matrix/barcodes"]][],
            gene_id = gene_id,
            gene_name = gene_name,
            data = h5file[["matrix/data"]][],
            indices = h5file[["matrix/indices"]][],
            indptr = h5file[["matrix/indptr"]][],
            shape = as.integer(h5file[["matrix/shape"]][]),
            gene_map = data.frame(
                ensembl = gene_id,
                symbol = gene_name,
                stringsAsFactors = FALSE
            ),
            backend = "hdf5r"
        )
    } else {
        stop("Neither rhdf5 nor hdf5r package is available. Please install one of them.")
    }
}

#' Construct a `dgCMatrix` from Xenium components.
#'
#' @param components List returned by `.read_xenium_sparse_components()`.
#' @param prefix Log prefix for user messages.
#' @param verbose Logical controlling verbosity.
#' @return List with the counts matrix and optional gene map.
#' @keywords internal
#' @noRd
.construct_xenium_sparse_matrix <- function(components, prefix, verbose) {
    shape <- components$shape
    nrow_h5 <- shape[1]
    ncol_h5 <- shape[2]
    indptr <- as.integer(components$indptr)
    by_cols <- length(indptr) == (ncol_h5 + 1L)
    if (!by_cols) {
        stop("Current Xenium implementation only supports CSC-style storage (indptr length = ncol + 1).")
    }

    barcodes <- components$barcodes
    genes <- components$gene_name
    col_is_cell <- (ncol_h5 == length(barcodes))
    if (!col_is_cell && !(nrow_h5 == length(barcodes))) {
        stop("Shape doesn't match barcode/gene counts, cannot determine orientation.")
    }

    if (verbose) {
        message(
            prefix, "  Matrix: ", nrow_h5, "×", ncol_h5,
            " (", length(genes), " genes, ", length(barcodes), " cells)"
        )
    }

    counts <- NULL
    if (col_is_cell) {
        counts <- methods::new(
            "dgCMatrix",
            Dim = c(length(genes), length(barcodes)),
            x = as.numeric(components$data),
            i = as.integer(components$indices),
            p = indptr,
            Dimnames = list(make.unique(genes), barcodes)
        )
    } else {
        tmp <- methods::new(
            "dgCMatrix",
            Dim = c(length(barcodes), length(genes)),
            x = as.numeric(components$data),
            i = as.integer(components$indices),
            p = indptr,
            Dimnames = list(barcodes, make.unique(genes))
        )
        if (verbose) message(prefix, "  Transposing matrix layout")
        counts <- Matrix::t(tmp)
    }

    list(
        counts = counts,
        gene_map = components$gene_map
    )
}

#' Apply gene filtering rules to a sparse matrix.
#'
#' @param counts Sparse counts matrix.
#' @param filters List holding gene inclusion/exclusion information.
#' @param verbose Logical flag for logs.
#' @param prefix Log prefix string.
#' @return List with filtered matrix and retained gene IDs.
#' @keywords internal
#' @noRd
.filter_gene_panel <- function(counts, filters, verbose, prefix) {
    gene_names <- rownames(counts)
    keep_genes <- gene_names

    exclude_prefix <- filters$gene_exclude_prefix
    if (length(exclude_prefix)) {
        pattern <- paste0("^(?:", paste(exclude_prefix, collapse = "|"), ")")
        drop_mask <- grepl(pattern, keep_genes, perl = TRUE)
        if (any(drop_mask)) {
            keep_genes <- keep_genes[!drop_mask]
            if (verbose) {
                message(prefix, "  Excluded ", sum(drop_mask), " genes with prefixes: ",
                        paste(exclude_prefix, collapse = ", "))
            }
        }
    }

    allowlist <- filters$gene_allowlist
    if (!is.null(allowlist)) {
        keep_genes <- intersect(keep_genes, allowlist)
        if (verbose) {
            message(prefix, "  Filtered to ", length(keep_genes), " genes from target list")
        }
    }

    if (!length(keep_genes)) {
        stop("No genes remain after filtering. Please check exclude_prefix and filter_genes parameters.")
    }

    if (!identical(keep_genes, gene_names)) {
        counts <- counts[keep_genes, , drop = FALSE]
    }

    if (verbose) {
        message(prefix, "  Final gene count: ", nrow(counts), "/", length(gene_names))
    }

    list(counts = counts, genes = keep_genes)
}

#' Align counts matrix columns to centroid targets.
#'
#' @param counts Sparse counts matrix.
#' @param targets Character vector of centroid identifiers.
#' @param filter_cells Logical flag to enforce alignment.
#' @param verbose Logical logging flag.
#' @param prefix Log prefix string.
#' @return List containing the aligned matrix and cell vector.
#' @keywords internal
#' @noRd
.align_counts_to_targets <- function(counts,
                                     targets,
                                     filter_cells,
                                     verbose,
                                     prefix) {
    if (!isTRUE(filter_cells) || is.null(targets)) {
        if (verbose) {
            message(prefix, "  Keeping all ", ncol(counts), " cells from matrix")
        }
        return(list(counts = counts, cells = colnames(counts)))
    }

    keep_cells <- targets[targets %in% colnames(counts)]
    if (!length(keep_cells)) {
        stop("Centroid cell IDs do not overlap with matrix barcodes.")
    }

    if (verbose) {
        message(prefix, "  Cell overlap: ", length(keep_cells), "/", length(targets))
    }

    list(
        counts = counts[, keep_cells, drop = FALSE],
        cells = keep_cells
    )
}

#' Persist the counts matrix into a scope object.
#'
#' @param scope_obj A `scope_object`.
#' @param counts Sparse counts matrix to store.
#' @param prefix Log prefix string.
#' @param verbose Logical flag for logging.
#' @return Invisibly, the modified `scope_object`.
#' @keywords internal
#' @noRd
.store_scope_counts <- function(scope_obj, counts, prefix, verbose) {
    scope_obj@cells <- list(counts = counts)
    if (verbose) message(prefix, " Cell matrix integration completed")
    invisible(scope_obj)
}

#' Locate the CosMx dataset and expression matrix.
#'
#' @param cosmx_root Root directory of the CosMx project.
#' @return List containing dataset directory and expression matrix path.
#' @keywords internal
#' @noRd
.locate_cosmx_dataset <- function(cosmx_root) {
    ff_dir <- file.path(cosmx_root, "flatFiles")
    if (!dir.exists(ff_dir)) {
        stop("flatFiles/ not found: ", ff_dir)
    }
    ds_dir <- list.dirs(ff_dir, full.names = TRUE, recursive = FALSE)
    if (!length(ds_dir)) {
        stop("No dataset directory found under flatFiles/.")
    }
    dataset <- ds_dir[[1]]
    prefix <- basename(dataset)
    candidates <- c(
        file.path(dataset, paste0(prefix, "_exprMat_file.csv.gz")),
        file.path(dataset, paste0(prefix, "_exprMat_file.csv"))
    )
    expr_mat <- NULL
    for (path in candidates) {
        if (file.exists(path)) {
            expr_mat <- path
            break
        }
    }
    if (is.null(expr_mat)) {
        stop("No *_exprMat_file.csv(.gz) found under: ", dataset)
    }
    list(dataset = dataset, expr_mat = expr_mat)
}

#' Read only the header of a CosMx expression matrix.
#'
#' @param expr_mat Path to the expression matrix CSV/CSV.GZ.
#' @return Character vector of column names.
#' @keywords internal
#' @noRd
.read_cosmx_header <- function(expr_mat) {
    header_dt <- tryCatch(
        data.table::fread(expr_mat, nrows = 0L, showProgress = FALSE),
        error = function(e) NULL
    )
    if (!is.null(header_dt)) {
        return(names(header_dt))
    }
    con <- NULL
    on.exit(if (!is.null(con)) try(close(con), silent = TRUE), add = TRUE)
    con <- if (grepl("\\.gz$", expr_mat, ignore.case = TRUE)) {
        gzfile(expr_mat, open = "rt")
    } else {
        file(expr_mat, open = "rt")
    }
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line) || !nzchar(line)) {
        stop("Unable to read expression matrix header: ", expr_mat)
    }
    line <- sub("\\r$", "", line)
    strsplit(line, ",", fixed = TRUE)[[1]]
}

#' Identify the FOV/cell columns and gene columns in CosMx data.
#'
#' @param cols Character vector of column names.
#' @param exclude_prefix Prefixes used to discard control genes.
#' @return List with `fov_col`, `cid_col` and `gene_cols`.
#' @keywords internal
#' @noRd
.identify_cosmx_gene_columns <- function(cols, exclude_prefix) {
    base_cols <- tolower(cols)
    fov_col <- cols[which(base_cols == "fov")[1]]
    cid_col <- cols[which(base_cols %in% c("cell_id", "cellid"))[1]]
    if (is.na(fov_col) || is.na(cid_col)) {
        stop("Expression matrix missing fov/cell_ID column")
    }

    non_gene_exact <- c(
        "fov", "cell_id", "cellid", "cell", "barcode",
        "x", "y", "z",
        "umi", "reads", "total", "sum", "sizefactor",
        "feature_name", "gene", "target", "qv", "nucleus_distance",
        "cellcomp", "compartment", "compartmentlabel"
    )
    non_gene_substr <- c(
        "_px", "_um", "_mm", "global_", "local_", "centroid", "area",
        "x_", "y_", "_x", "_y"
    )
    is_non_gene <- base_cols %in% non_gene_exact
    for (ss in non_gene_substr) {
        is_non_gene <- is_non_gene | grepl(ss, base_cols, fixed = TRUE)
    }
    gene_cols <- cols[!is_non_gene]

    if (length(exclude_prefix)) {
        pattern <- paste0("^(?:", paste(exclude_prefix, collapse = "|"), ")")
        gene_cols <- gene_cols[!grepl(pattern, gene_cols, perl = TRUE)]
    }

    if (!length(gene_cols)) {
        stop("No gene columns identified. Check matrix header and prefix filters.")
    }

    list(
        fov_col = fov_col,
        cid_col = cid_col,
        gene_cols = gene_cols
    )
}

#' Determine the CosMx gene panel to read.
#'
#' @param filter_genes Optional user-provided gene list.
#' @param auto_genes Automatically inferred gene columns.
#' @param force_all_genes Logical forcing full gene panel usage.
#' @param verbose Logical flag for logging.
#' @param prefix Log prefix string.
#' @return Character vector of genes to retain.
#' @keywords internal
#' @noRd
.resolve_cosmx_gene_panel <- function(filter_genes,
                                      auto_genes,
                                      force_all_genes,
                                      verbose,
                                      prefix) {
    if (is.null(filter_genes)) {
        if (!isTRUE(force_all_genes) && verbose) {
            message(prefix, " No filter_genes supplied; reading auto-detected gene panel (", length(auto_genes), " genes)")
        }
        gene_cols <- auto_genes
    } else {
        gene_cols <- intersect(auto_genes, filter_genes)
        if (!length(gene_cols)) {
            stop("No overlap between filter_genes and available genes (after auto-detection).")
        }
    }
    if (verbose) {
        message(prefix, " Number of genes selected: ", length(gene_cols))
    }
    gene_cols
}

#' Read a subset of CosMx expression columns.
#'
#' @param expr_mat Path to the expression matrix file.
#' @param selected_cols Character vector of columns to load.
#' @param verbose Logical flag for progress reporting.
#' @param prefix Log prefix string.
#' @return `data.table` with the requested columns.
#' @keywords internal
#' @noRd
.read_cosmx_expression_subset <- function(expr_mat, selected_cols, verbose, prefix) {
    dt <- tryCatch(
        data.table::fread(expr_mat, select = selected_cols, showProgress = verbose),
        error = function(e) {
            if (verbose) {
                message(prefix, " fread failed; falling back to read.csv: ", e$message)
            }
            utils::read.csv(expr_mat, header = TRUE)[, selected_cols]
        }
    )
    data.table::setDT(dt)
    dt
}

#' Generate CosMx column keys compatible with centroid identifiers.
#'
#' @param expr_dt Expression data table.
#' @param cid_col Column name containing cell IDs.
#' @param fov_col Column name containing FOV indices.
#' @param id_mode Identifier mode (`"cell_id"`, `"fov_cell"`, `"auto"`).
#' @param centroid_keys Target centroid identifiers.
#' @param filter_cells Logical flag to restrict to centroid overlap.
#' @param verbose Logical for log messages.
#' @param prefix Log prefix string.
#' @return List with updated `expr_dt` and ordered cell key vector.
#' @keywords internal
#' @noRd
.assign_cosmx_cell_keys <- function(expr_dt,
                                    cid_col,
                                    fov_col,
                                    id_mode,
                                    centroid_keys,
                                    filter_cells,
                                    verbose,
                                    prefix) {
    id_mode_eff <- id_mode
    if (identical(id_mode_eff, "auto")) {
        if (isTRUE(filter_cells) && length(centroid_keys) &&
            (any(grepl("_", centroid_keys, fixed = TRUE)) || any(grepl("FOV", centroid_keys, ignore.case = TRUE)))) {
            id_mode_eff <- "fov_cell"
        } else {
            id_mode_eff <- "cell_id"
        }
    }

    make_key <- function(mode, fov_offset = 0L) {
        if (identical(mode, "cell_id") || is.na(fov_col)) {
            return(as.character(expr_dt[[cid_col]]))
        }
        fov_val <- suppressWarnings(as.integer(expr_dt[[fov_col]]))
        if (!is.na(fov_offset) && fov_offset != 0L) {
            fov_val <- fov_val + as.integer(fov_offset)
        }
        paste0(as.character(expr_dt[[cid_col]]), "_", as.character(fov_val))
    }

    expr_dt[, cell_key := make_key(id_mode_eff, 0L)]

    if (!isTRUE(filter_cells) || is.null(centroid_keys)) {
        column_order <- unique(expr_dt$cell_key)
        if (verbose) {
            message(prefix, " Keeping all ", length(column_order), " cells from matrix (no barcode filter)")
        }
        return(list(expr_dt = expr_dt, column_order = column_order))
    }

    target_unique <- unique(centroid_keys)
    candidates <- list(
        list(mode = id_mode_eff, offset = 0L),
        list(mode = if (identical(id_mode_eff, "cell_id")) "fov_cell" else "cell_id", offset = 0L),
        list(mode = "fov_cell", offset = +1L),
        list(mode = "fov_cell", offset = -1L)
    )

    best_idx <- 1L
    best_overlap <- -1L
    key_cache <- new.env(parent = emptyenv())
    for (i in seq_along(candidates)) {
        cand <- candidates[[i]]
        mode_i <- cand$mode
        offset_i <- cand$offset
        if (is.na(fov_col) && identical(mode_i, "fov_cell")) next
        cache_key <- paste(mode_i, offset_i, sep = ":")
        if (!exists(cache_key, envir = key_cache, inherits = FALSE)) {
            assign(cache_key, make_key(mode_i, offset_i), envir = key_cache)
        }
        keys_vec <- get(cache_key, envir = key_cache)
        overlap <- length(intersect(unique(keys_vec), target_unique))
        if (overlap > best_overlap) {
            best_overlap <- overlap
            best_idx <- i
        }
    }

    best_cand <- candidates[[best_idx]]
    expr_dt[, cell_key := make_key(best_cand$mode, best_cand$offset)]

    keep_cells <- intersect(unique(expr_dt$cell_key), target_unique)
    if (!length(keep_cells)) {
        stop("Cell identifiers in expression matrix do not overlap with centroids; check id_mode/FOV encoding.")
    }

    expr_dt <- expr_dt[cell_key %in% keep_cells]
    column_order <- centroid_keys[centroid_keys %in% keep_cells]
    if (verbose) {
        message(prefix, " Cell overlap: ", length(column_order), "/", length(centroid_keys))
    }
    list(expr_dt = expr_dt, column_order = column_order)
}

#' Assemble a CosMx sparse matrix.
#'
#' @param expr_dt `data.table` containing gene expression columns.
#' @param gene_cols Character vector of gene column names.
#' @param column_order Desired column order for cells.
#' @return A `dgCMatrix` of genes by cells.
#' @keywords internal
#' @noRd
.assemble_cosmx_sparse_matrix <- function(expr_dt, gene_cols, column_order) {
    if (!length(column_order)) {
        stop("No cells available after filtering; cannot build sparse matrix.")
    }

    col_index <- match(expr_dt$cell_key, column_order)
    gene_index <- seq_along(gene_cols)
    names(gene_index) <- gene_cols

    nz_entries <- lapply(gene_cols, function(gene) {
        values <- suppressWarnings(as.numeric(expr_dt[[gene]]))
        nz <- which(!is.na(values) & values != 0)
        if (!length(nz)) {
            return(NULL)
        }
        list(
            i = rep.int(gene_index[[gene]], length(nz)),
            j = col_index[nz],
            x = values[nz]
        )
    })
    nz_entries <- nz_entries[!vapply(nz_entries, is.null, logical(1))]
    if (!length(nz_entries)) {
        Matrix::sparseMatrix(
            i = integer(0),
            j = integer(0),
            x = numeric(0),
            dims = c(length(gene_cols), length(column_order)),
            dimnames = list(gene_cols, column_order)
        )
    } else {
        i_idx <- unlist(lapply(nz_entries, `[[`, "i"), use.names = FALSE)
        j_idx <- unlist(lapply(nz_entries, `[[`, "j"), use.names = FALSE)
        x_val <- unlist(lapply(nz_entries, `[[`, "x"), use.names = FALSE)
        Matrix::sparseMatrix(
            i = as.integer(i_idx),
            j = as.integer(j_idx),
            x = as.numeric(x_val),
            dims = c(length(gene_cols), length(column_order)),
            dimnames = list(gene_cols, column_order)
        )
    }
}
