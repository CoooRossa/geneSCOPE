#' Auto Histology Crop Bbox
#' @description
#' Internal helper for `.auto_histology_crop_bbox`.
#' @param roi_bbox Parameter value.
#' @param microns_per_pixel Parameter value.
#' @param level_scalefactor Parameter value.
#' @return Return value used internally.
#' @keywords internal
.auto_histology_crop_bbox <- function(roi_bbox,
                                      microns_per_pixel,
                                      level_scalefactor) {
    if (is.null(roi_bbox) || !is.finite(microns_per_pixel) || microns_per_pixel <= 0) {
        return(NULL)
    }
    req <- c("xmin", "xmax", "ymin", "ymax")
    if (!all(req %in% names(roi_bbox))) return(NULL)
    roi_vals <- setNames(as.numeric(roi_bbox[req]), req)
    if (any(!is.finite(roi_vals))) return(NULL)
    px_fullres <- roi_vals / microns_per_pixel
    scale_use <- .as_numeric_or_na(level_scalefactor)
    if (!is.finite(scale_use) || scale_use <= 0) scale_use <- 1
    px_level <- px_fullres * scale_use
    c(
        xmin = floor(px_level["xmin"]),
        xmax = ceiling(px_level["xmax"]),
        ymin = floor(px_level["ymin"]),
        ymax = ceiling(px_level["ymax"])
    )
}

#' Summarize coordinate ranges for debugging.
#' @description
#' Internal helper used for ingest diagnostics. Returns finite min/max for x/y
#' columns when present; otherwise returns NA values.
#' @param df A data.frame/data.table.
#' @param x_col Name of x column.
#' @param y_col Name of y column.
#' @return Named numeric vector with xmin/xmax/ymin/ymax/dx/dy.
#' @keywords internal
.summarize_xy_ranges <- function(df, x_col = "x", y_col = "y") {
    out <- c(xmin = NA_real_, xmax = NA_real_, ymin = NA_real_, ymax = NA_real_, dx = NA_real_, dy = NA_real_)
    if (is.null(df) || !is.data.frame(df)) return(out)
    if (!(x_col %in% names(df)) || !(y_col %in% names(df))) return(out)
    x <- suppressWarnings(as.numeric(df[[x_col]]))
    y <- suppressWarnings(as.numeric(df[[y_col]]))
    if (!length(x) || !length(y)) return(out)
    xr <- suppressWarnings(range(x, na.rm = TRUE))
    yr <- suppressWarnings(range(y, na.rm = TRUE))
    if (any(!is.finite(xr))) xr <- c(NA_real_, NA_real_)
    if (any(!is.finite(yr))) yr <- c(NA_real_, NA_real_)
    out["xmin"] <- xr[1]; out["xmax"] <- xr[2]
    out["ymin"] <- yr[1]; out["ymax"] <- yr[2]
    out["dx"] <- if (all(is.finite(xr))) diff(xr) else NA_real_
    out["dy"] <- if (all(is.finite(yr))) diff(yr) else NA_real_
    out
}

#' Read Xenium experiment metadata.
#' @description
#' Reads `experiment.xenium` and returns a small metadata list.
#' @param input_dir Xenium outs/ directory.
#' @return List (or NULL if missing/unreadable).
#' @keywords internal
.read_experiment_xenium <- function(input_dir) {
    exp_path <- file.path(input_dir, "experiment.xenium")
    if (!file.exists(exp_path)) return(NULL)
    if (!requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
    meta <- tryCatch(jsonlite::fromJSON(exp_path), error = function(e) NULL)
    if (is.null(meta) || !is.list(meta)) return(NULL)
    get1 <- function(k) if (!is.null(meta[[k]])) meta[[k]] else NULL
    list(
        major_version = suppressWarnings(as.integer(get1("major_version"))),
        minor_version = suppressWarnings(as.integer(get1("minor_version"))),
        patch_version = suppressWarnings(as.integer(get1("patch_version"))),
        chemistry_version = as.character(get1("chemistry_version")),
        pixel_size_um = suppressWarnings(as.numeric(get1("pixel_size"))),
        run_name = as.character(get1("run_name")),
        region_name = as.character(get1("region_name"))
    )
}

.format_xenium_version <- function(meta) {
    if (is.null(meta) || !is.list(meta)) return(NA_character_)
    mj <- suppressWarnings(as.integer(meta$major_version))
    mn <- suppressWarnings(as.integer(meta$minor_version))
    pt <- suppressWarnings(as.integer(meta$patch_version))
    if (!is.finite(mj)) return(NA_character_)
    if (!is.finite(mn)) mn <- 0L
    if (!is.finite(pt)) pt <- 0L
    paste0(mj, ".", mn, ".", pt)
}

.infer_xenium_generation <- function(meta) {
    mj <- suppressWarnings(as.integer(if (!is.null(meta)) meta$major_version else NA_integer_))
    if (!is.finite(mj)) return("unknown")
    if (mj >= 6L) return("modern")
    "legacy"
}

.ratio_close <- function(x, target, tol_rel = 0.25) {
    x <- suppressWarnings(as.numeric(x))
    target <- suppressWarnings(as.numeric(target))
    if (!is.finite(x) || !is.finite(target) || target == 0) return(FALSE)
    abs(x - target) / abs(target) <= tol_rel
}

.xenium_harmonize_coordinate_units <- function(state, ds, bounds_global, verbose) {
    parent <- "createSCOPE"
    if (is.null(state) || !is.list(state)) return(list(state = state, ds = ds, bounds_global = bounds_global))
    if (!is.null(state$coord_units) && isTRUE(state$coord_units$harmonized)) {
        return(list(state = state, ds = ds, bounds_global = bounds_global))
    }
    exp_meta <- state$experiment
    pix <- suppressWarnings(as.numeric(if (!is.null(exp_meta)) exp_meta$pixel_size_um else NA_real_))
    if (!is.finite(pix) || pix <= 0) {
        state$coord_units <- list(harmonized = TRUE, pixel_size_um = pix, action = "none", reason = "pixel_size_missing")
        return(list(state = state, ds = ds, bounds_global = bounds_global))
    }
    cent <- state$objects$centroids
    cen_rng <- .summarize_xy_ranges(cent, x_col = "x", y_col = "y")
    bdx <- suppressWarnings(as.numeric(bounds_global$xmax) - as.numeric(bounds_global$xmin))
    bdy <- suppressWarnings(as.numeric(bounds_global$ymax) - as.numeric(bounds_global$ymin))
    cdx <- cen_rng[["dx"]]
    cdy <- cen_rng[["dy"]]
    if (!is.finite(bdx) || !is.finite(bdy) || !is.finite(cdx) || !is.finite(cdy) || cdx <= 0 || cdy <= 0) {
        state$coord_units <- list(harmonized = TRUE, pixel_size_um = pix, action = "none", reason = "nonfinite_ranges")
        return(list(state = state, ds = ds, bounds_global = bounds_global))
    }

    rx <- bdx / cdx
    ry <- bdy / cdy
    r <- median(c(rx, ry), na.rm = TRUE)
    if (!is.finite(r) || r <= 0) {
        state$coord_units <- list(harmonized = TRUE, pixel_size_um = pix, action = "none", reason = "ratio_nonfinite")
        return(list(state = state, ds = ds, bounds_global = bounds_global))
    }

    # If transcripts look like pixels and centroids look like microns: bounds/centroids ~= 1/pixel_size
    if (.ratio_close(r, 1 / pix, tol_rel = 0.35)) {
        scale <- pix
        ds2 <- tryCatch(
            ds |> dplyr::mutate(
                x_location = x_location * scale,
                y_location = y_location * scale
            ),
            error = function(e) e
        )
        if (inherits(ds2, "error")) {
            warning("Auto unit harmonization detected transcript coordinates in pixels, but could not scale Arrow dataset: ", conditionMessage(ds2), call. = FALSE)
            state$coord_units <- list(harmonized = TRUE, pixel_size_um = pix, action = "none", reason = "arrow_mutate_failed")
            return(list(state = state, ds = ds, bounds_global = bounds_global))
        }
        ds <- ds2
        bounds_global$xmin <- bounds_global$xmin * scale
        bounds_global$xmax <- bounds_global$xmax * scale
        bounds_global$ymin <- bounds_global$ymin * scale
        bounds_global$ymax <- bounds_global$ymax * scale
        if (!is.null(state$datasets$transcripts_cache) && isTRUE(state$datasets$transcripts_cache_masked)) {
            cache_dt <- state$datasets$transcripts_cache
            if (is.data.frame(cache_dt) && all(c("x_location", "y_location") %in% names(cache_dt))) {
                cache_dt$x_location <- suppressWarnings(as.numeric(cache_dt$x_location)) * scale
                cache_dt$y_location <- suppressWarnings(as.numeric(cache_dt$y_location)) * scale
                state$datasets$transcripts_cache <- cache_dt
            }
        }
        .log_info(
            parent,
            "S08",
            paste0("auto unit harmonization: scaled transcripts by pixel_size_um=", format(pix, digits = 6, trim = TRUE), " (px->umm)"),
            verbose
        )
        state$coord_units <- list(harmonized = TRUE, pixel_size_um = pix, action = "scale_transcripts", scale = scale, ratio_before = r)
        return(list(state = state, ds = ds, bounds_global = bounds_global))
    }

    # If centroids look like pixels and transcripts look like microns: bounds/centroids ~= pixel_size
    if (.ratio_close(r, pix, tol_rel = 0.35)) {
        scale <- pix
        if (is.data.frame(cent) && all(c("x", "y") %in% names(cent))) {
            cent$x <- suppressWarnings(as.numeric(cent$x)) * scale
            cent$y <- suppressWarnings(as.numeric(cent$y)) * scale
            state$objects$centroids <- cent
        }
        .log_info(
            parent,
            "S08",
            paste0("auto unit harmonization: scaled centroids by pixel_size_um=", format(pix, digits = 6, trim = TRUE), " (px->umm)"),
            verbose
        )
        state$coord_units <- list(harmonized = TRUE, pixel_size_um = pix, action = "scale_centroids", scale = scale, ratio_before = r)
        return(list(state = state, ds = ds, bounds_global = bounds_global))
    }

    state$coord_units <- list(harmonized = TRUE, pixel_size_um = pix, action = "none", reason = "ratio_near_1", ratio = r)
    list(state = state, ds = ds, bounds_global = bounds_global)
}

#' Clip Points To Polygon
#' @description
#' Internal helper for `.clip_points_to_polygon`.
#' @param points Parameter value.
#' @param polygon Parameter value.
#' @param chunk_size Parameter value.
#' @param ncores Number of cores/threads to use.
#' @return Return value used internally.
#' @keywords internal
.clip_points_to_polygon <- function(points, polygon,
                                 chunk_size = 5e5, ncores = 1) {
    if (is.null(polygon) || nrow(points) == 0) {
        return(points)
    }

    ## -------- 0. Check sf package availability --------------------------------
    if (!requireNamespace("sf", quietly = TRUE)) {
        stop("Package 'sf' is required for polygon clipping. Please install it.")
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

    ## (1)-- If length 0, return empty result to avoid crash ------------------------------
    if (length(polygon) == 0) {
        return(points[0])
    } # empty data.table / data.frame

    ## (Optional) Only union if more than one feature; prevents EMPTY
    if (length(polygon) > 1) {
        polygon <- sf::st_union(polygon)
    }

    ## (2)-- Check again after union
    if (length(polygon) == 0) {
        return(points[0])
    }

    ## -------- 2. Preprocessing & parallel clipping (as before) -----------------------------
    stopifnot(all(c("x", "y") %in% names(points)))

    dt_pts <- as.data.table(points)

    idx_split <- split(
        seq_len(nrow(dt_pts)),
        ceiling(seq_len(nrow(dt_pts)) / chunk_size)
    )

    worker <- function(idx) {
        sub <- dt_pts[idx]
        inside <- lengths(sf::st_within(
            sf::st_as_sf(sub,
                coords = c("x", "y"),
            crs = sf::st_crs(polygon) %||% NA, remove = FALSE
            ),
            polygon
        )) > 0
        sub[inside]
    }

    res_list <- if (ncores > 1 && .Platform$OS.type != "windows") {
        mclapply(idx_split, worker, mc.cores = ncores)
    } else if (ncores > 1) {
        cl <- makeCluster(ncores)
        on.exit(stopCluster(cl))
        clusterEvalQ(cl, {
            if (!requireNamespace("sf", quietly = TRUE)) {
                stop("Package 'sf' is required on worker. Please install it.")
            }
            library(sf, quietly = TRUE)
        })
        clusterExport(cl, c("dt_pts", "polygon", "worker"),
            envir = environment()
        )
        parLapply(cl, idx_split, worker)
    } else {
        lapply(idx_split, worker)
    }

    rbindlist(res_list)
}

#' Map full-resolution Visium coordinates into a histology level.
#' @description
#' Internal helper for `.coords_fullres_to_level`.
#' @param df data.frame with spot columns.
#' @param x_col Column names storing Visium pxl_* values.
#' @param y_col Column names storing Visium pxl_* values.
#' @param histo Histology entry from `grid$histology[[level]]`.
#' @return data.frame with `x_img`/`y_img`.
#' @keywords internal
.coords_fullres_to_level <- function(df,
                                     x_col,
                                     y_col,
                                     histo) {
    stopifnot(all(c(x_col, y_col) %in% names(df)))
    if (is.null(histo$scalefactor) || !is.finite(histo$scalefactor)) {
        stop("Histology entry is missing a numeric scalefactor.")
    }
    if (is.null(histo$height) || !is.finite(histo$height)) {
        stop("Histology entry is missing `height` (pixels).")
    }
    sf <- as.numeric(histo$scalefactor)
    img_h <- as.numeric(histo$height)
    x_full <- as.numeric(df[[x_col]])
    y_full <- as.numeric(df[[y_col]])
    df$x_img <- x_full * sf
    df$y_img <- img_h - (y_full * sf)
    df
}

#' Map physical (micron) coordinates to histology pixel space using ROI bounds.
#' @description
#' Internal helper for `.coords_physical_to_level`.
#' @param df data.frame with numeric columns.
#' @param x_col Column names to transform.
#' @param y_col Column names to transform.
#' @param histo Histology entry from `grid$histology`.
#' @return data.frame with `x_img`/`y_img`.
#' @keywords internal
.coords_physical_to_level <- function(df,
                                      x_col,
                                      y_col,
                                      histo) {
    stopifnot(all(c(x_col, y_col) %in% names(df)))
    roi <- histo$roi_bbox
    if (is.null(roi)) stop("Histology entry missing `roi_bbox`; cannot map coordinates.")
    req <- c("xmin", "xmax", "ymin", "ymax")
    if (!all(req %in% names(roi))) {
        stop("roi_bbox must include xmin/xmax/ymin/ymax.")
    }
    width <- as.numeric(histo$width)
    height <- as.numeric(histo$height)
    if (!all(is.finite(c(width, height))) || any(c(width, height) <= 0)) {
        stop("Histology entry must define positive numeric width/height.")
    }

    rx <- roi["xmax"] - roi["xmin"]
    ry <- roi["ymax"] - roi["ymin"]
    if (rx == 0 || ry == 0) stop("roi_bbox has zero extent; cannot transform coordinates.")

    x_phys <- as.numeric(df[[x_col]])
    y_phys <- as.numeric(df[[y_col]])

    df$x_img <- (x_phys - roi["xmin"]) / rx * width
    df$y_img <- (1 - (y_phys - roi["ymin"]) / ry) * height
    df
}

#' Flip Coordinates
#' @description
#' Internal helper for `.flip_coordinates`.
#' @param data Parameter value.
#' @param y_max Numeric threshold.
#' @return Return value used internally.
#' @keywords internal
.flip_coordinates <- function(data, y_max) {
    stopifnot(is.numeric(y_max) && length(y_max) == 1)
    if (!"y" %in% names(data)) {
        stop("Input data must have a 'y' column for coordinates.")
    }
    # If data is a data.table, modify in place for efficiency
    if (is.data.table(data)) {
        data[, y := y_max - y]
    } else {
        data$y <- y_max - data$y
    }
    return(data)
}

#' Process Segmentation
#' @description
#' Internal helper for `.process_segmentation`.
#' @param seg_file Filesystem path.
#' @param tag Parameter value.
#' @param flip_y Parameter value.
#' @param y_max Numeric threshold.
#' @param keep_cells Logical flag.
#' @param ncores Number of cores/threads to use.
#' @return Return value used internally.
#' @keywords internal
.process_segmentation <- function(seg_file, tag = "cell",
                                 flip_y = FALSE, y_max = NULL,
                                 keep_cells = NULL, ncores = 1) {
    if (!file.exists(seg_file)) {
        stop("Segmentation file not found: ", seg_file)
    }

    seg_raw <- as.data.table(arrow::read_parquet(seg_file))
    # When an fov column exists, build a unique key paste(cell_id, fov) and update label_id accordingly
    fov_col <- NULL
    if ("fov" %in% names(seg_raw)) {
        fov_col <- "fov"
    } else if ("fov_name" %in% names(seg_raw)) {
        fov_col <- "fov_name"
    }
    if (!is.null(fov_col)) {
        seg_dt <- seg_raw[, .(
            cell     = paste0(as.character(cell_id), "_", as.character(get(fov_col))),
            x        = vertex_x,
            y        = vertex_y,
            label_id = paste0(as.character(cell_id), "_", as.character(get(fov_col)))
        )]
    } else {
        label_col <- if ("label_id" %in% names(seg_raw)) "label_id" else "cell_id"
        seg_dt <- seg_raw[, .(
            cell     = as.character(cell_id),
            x        = vertex_x,
            y        = vertex_y,
            label_id = as.character(get(label_col))
        )]
    }

    if (flip_y) {
        if (is.null(y_max)) {
            stop("y_max must be provided when flip_y = TRUE.")
        }
        seg_dt[, y := y_max - y]
    }

    if (!is.null(keep_cells)) {
        seg_dt <- seg_dt[cell %in% keep_cells]
    }

    ## Build polygons (ensure closure) - PDF suggestion
    split_idx <- split(seq_len(nrow(seg_dt)), seg_dt$label_id)
    build_poly <- function(idx) {
        sub <- seg_dt[idx, .(x, y)]
        sub <- sub[complete.cases(sub)] # (1) Remove NA
        if (nrow(sub) < 3) {
            return(NULL)
        }

        # (2) Ensure double storage
        coords <- as.matrix(sub)
        storage.mode(coords) <- "double"

        # Close polygon if not already closed
        if (!all(coords[1, ] == coords[nrow(coords), ])) {
            coords <- rbind(coords, coords[1, ])
        }

        # (3) try-catch, return NULL on failure
        tryCatch(sf::st_polygon(list(coords)), error = function(e) {
            message("label ", sub$label_id[1], " invalid: ", conditionMessage(e))
            NULL
        })
    }

    poly_list <- if (ncores > 1 && .Platform$OS.type != "windows") {
        mclapply(split_idx, build_poly, mc.cores = ncores)
    } else if (ncores > 1) {
        cl <- makeCluster(ncores)
        on.exit(stopCluster(cl))
        clusterEvalQ(cl, {
            if (!requireNamespace("sf", quietly = TRUE)) {
                stop("Package 'sf' is required on worker. Please install it.")
            }
            library(sf, quietly = TRUE)
        })
        clusterExport(cl, c("seg_dt", "build_poly"),
            envir = environment()
        )
        parLapply(cl, split_idx, build_poly)
    } else {
        lapply(split_idx, build_poly)
    }

    polygons <- sf::st_sfc(poly_list[!sapply(poly_list, is.null)])
    list(points = seg_dt, polygons = polygons)
}

#' Add Singlecells Cosmx
#' @description
#' Internal helper for `.add_singlecells_cosmx`.
#' @param scope_obj A `scope_object`.
#' @param cosmx_root Parameter value.
#' @param filter_genes Parameter value.
#' @param exclude_prefix Parameter value.
#' @param id_mode Parameter value.
#' @param filter_barcodes Parameter value.
#' @param force_all_genes Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
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

#' Add Singlecells Xenium
#' @description
#' Internal helper for `.add_singlecells_xenium`.
#' @param scope_obj A `scope_object`.
#' @param xenium_dir Filesystem path.
#' @param filter_genes Parameter value.
#' @param exclude_prefix Parameter value.
#' @param filter_barcodes Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
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

#' Align Counts To Targets
#' @description
#' Internal helper for `.align_counts_to_targets`.
#' @param counts Parameter value.
#' @param targets Parameter value.
#' @param filter_cells Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param prefix Parameter value.
#' @return Return value used internally.
#' @keywords internal
.align_counts_to_targets <- function(counts,
                                     targets,
                                     filter_cells,
                                     verbose,
                                     prefix) {
    parent <- "createSCOPE"
    if (!isTRUE(filter_cells) || is.null(targets)) {
        .log_info(parent, "S01", paste0("keeping all ", ncol(counts), " cells from matrix"), verbose)
        return(list(counts = counts, cells = colnames(counts)))
    }

    keep_cells <- targets[targets %in% colnames(counts)]
    if (!length(keep_cells)) {
        stop("Centroid cell IDs do not overlap with matrix barcodes.")
    }

    .log_info(parent, "S01", paste0("cell overlap: ", length(keep_cells), "/", length(targets)), verbose)

    list(
        counts = counts[, keep_cells, drop = FALSE],
        cells = keep_cells
    )
}

#' Resolve gridSelectQC output path for gene flags.
#' @description
#' Internal helper used by createSCOPE() to optionally exclude genes flagged as
#' unsuitable by `gridSelectQC()` / `grid_gene_bins.tsv.gz`.
#' @param qc_out_dir Path passed to gridSelectQC(out_dir=...).
#' @return Filesystem path to `grid_gene_bins.tsv(.gz)`.
#' @keywords internal
.gridselect_qc_resolve_gene_bins_path <- function(qc_out_dir) {
    if (is.null(qc_out_dir) || !nzchar(qc_out_dir)) return(NULL)
    qc_out_dir <- as.character(qc_out_dir)[1]
    if (is.na(qc_out_dir) || !nzchar(qc_out_dir)) return(NULL)
    candidates <- c(
        file.path(qc_out_dir, "plots", "grid_gene_bins.tsv.gz"),
        file.path(qc_out_dir, "grid_gene_bins.tsv.gz"),
        file.path(qc_out_dir, "plots", "grid_gene_bins.tsv"),
        file.path(qc_out_dir, "grid_gene_bins.tsv")
    )
    hit <- candidates[file.exists(candidates)]
    if (!length(hit)) return(NULL)
    hit[[1]]
}

.gridselect_qc_flag_to_logical <- function(x) {
    out <- rep(FALSE, length(x))
    if (!length(x)) return(out)
    if (is.logical(x)) {
        out[which(!is.na(x))] <- x[which(!is.na(x))]
        return(out)
    }
    if (is.numeric(x) || is.integer(x)) {
        out[which(!is.na(x))] <- x[which(!is.na(x))] != 0
        return(out)
    }
    idx <- which(!is.na(x))
    if (!length(idx)) return(out)
    v <- tolower(trimws(as.character(x[idx])))
    out[idx] <- v %in% c("true", "t", "1", "yes", "y")
    out
}

.gridselect_qc_read_tsv <- function(path) {
    if (is.null(path) || !nzchar(path)) stop("Path is empty.")
    path <- as.character(path)[1]
    if (is.na(path) || !nzchar(path)) stop("Path is empty.")
    is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)
    if (is_gz) {
        con <- gzfile(path, open = "rt")
        on.exit(try(close(con), silent = TRUE), add = TRUE)
        return(utils::read.delim(con,
            sep = "\t",
            stringsAsFactors = FALSE,
            check.names = FALSE
        ))
    }
    if (requireNamespace("data.table", quietly = TRUE)) {
        return(data.table::fread(path))
    }
    utils::read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

#' Load QC-flagged genes to exclude.
#' @description
#' Reads `grid_gene_bins.tsv(.gz)` produced by `gridSelectQC()` and returns a
#' vector of gene symbols to exclude.
#' @param qc_out_dir Path passed to gridSelectQC(out_dir=...).
#' @param qc_flag Optional character vector selecting flags to exclude
#' (`near_zero`, `scale_mismatch`, `loss`, `deg_01`, `any`). If NULL, excludes
#' all genes where `unsuitable_any==TRUE` (or any flag column when `unsuitable_any`
#' is unavailable).
#' @return Character vector of gene symbols.
#' @keywords internal
.gridselect_qc_genes_to_exclude <- function(qc_out_dir, qc_flag = NULL) {
    qc_path <- .gridselect_qc_resolve_gene_bins_path(qc_out_dir)
    if (is.null(qc_path)) {
        stop("QC gene bins file not found under qc_out_dir: ", qc_out_dir,
             "\nExpected one of: plots/grid_gene_bins.tsv.gz, grid_gene_bins.tsv.gz, plots/grid_gene_bins.tsv, grid_gene_bins.tsv")
    }

    dt <- tryCatch(.gridselect_qc_read_tsv(qc_path), error = function(e) {
        stop("Failed to read QC gene bins file: ", qc_path, " (", conditionMessage(e), ")")
    })

    if (is.null(dt) || !nrow(dt)) return(character())
    if (!"gene" %in% names(dt)) {
        stop("QC gene bins file missing required column `gene`: ", qc_path)
    }

    genes <- as.character(dt[["gene"]])
    ok_gene <- !is.na(genes) & nzchar(genes) & genes != "__SUMMARY__"
    if (!any(ok_gene)) return(character())

    requested <- qc_flag
    if (!is.null(requested)) {
        requested <- as.character(requested)
        requested <- requested[!is.na(requested) & nzchar(requested)]
        if (!length(requested)) requested <- NULL
    }

    cols <- character()
    if (is.null(requested)) {
        if ("unsuitable_any" %in% names(dt)) {
            cols <- "unsuitable_any"
        } else {
            cols <- intersect(c("flag_near_zero", "flag_scale_mismatch", "flag_loss", "flag_deg_01"), names(dt))
        }
    } else {
        norm <- tolower(requested)
        map <- list(
            near_zero = "flag_near_zero",
            scale_mismatch = "flag_scale_mismatch",
            loss = "flag_loss",
            deg_01 = "flag_deg_01",
            any = "unsuitable_any"
        )
        for (f in norm) {
            if (!f %in% names(map)) {
                stop("Unknown qc_flag: ", f, "; expected one of: near_zero, scale_mismatch, loss, deg_01, any")
            }
            cols <- c(cols, map[[f]])
        }
        cols <- unique(cols)
        missing_cols <- setdiff(cols, names(dt))
        if (length(missing_cols)) {
            stop("QC gene bins file missing required column(s): ", paste(missing_cols, collapse = ", "),
                 "\nFile: ", qc_path)
        }
    }

    if (!length(cols)) return(character())
    mask <- rep(FALSE, nrow(dt))
    for (col in cols) {
        mask <- mask | .gridselect_qc_flag_to_logical(dt[[col]])
    }
    unique(genes[ok_gene & mask])
}

# Resolve Xenium codebook/codeword mapping for newer transcript formats that
# emit `codeword_index` + `is_gene` instead of `feature_name`.
.resolve_xenium_codebook_path <- function(input_dir) {
    if (is.null(input_dir) || !length(input_dir)) return(NULL)
    input_dir <- as.character(input_dir)[1]
    if (is.na(input_dir) || !nzchar(input_dir) || !dir.exists(input_dir)) return(NULL)

    opt <- getOption("geneSCOPE.xenium_codebook_path", NULL)
    if (!is.null(opt)) {
        opt <- as.character(opt)[1]
        if (!is.na(opt) && nzchar(opt) && file.exists(opt)) {
            return(opt)
        }
        opt_rel <- file.path(input_dir, opt)
        if (!is.na(opt) && nzchar(opt) && file.exists(opt_rel)) {
            return(opt_rel)
        }
    }

    candidates <- c(
        file.path(input_dir, "codebook.csv.gz"),
        file.path(input_dir, "codebook.csv"),
        file.path(input_dir, "codebook.tsv.gz"),
        file.path(input_dir, "codebook.tsv"),
        file.path(input_dir, "analysis", "codebook.csv.gz"),
        file.path(input_dir, "analysis", "codebook.csv"),
        file.path(input_dir, "analysis", "codebook.tsv.gz"),
        file.path(input_dir, "analysis", "codebook.tsv")
    )
    hit <- candidates[file.exists(candidates)]
    if (length(hit)) return(hit[[1]])

    found <- tryCatch(
        list.files(input_dir,
            pattern = "codebook",
            recursive = TRUE,
            full.names = TRUE,
            ignore.case = TRUE
        ),
        error = function(e) character()
    )
    if (!length(found)) return(NULL)
    found <- found[grepl("\\.(csv|tsv)(\\.gz)?$", found, ignore.case = TRUE)]
    if (!length(found)) return(NULL)
    found[[1]]
}

.read_delim_any <- function(path, verbose = FALSE) {
    if (requireNamespace("data.table", quietly = TRUE)) {
        return(data.table::fread(path, showProgress = isTRUE(verbose)))
    }
    ext <- tolower(tools::file_ext(path))
    sep <- if (identical(ext, "csv")) "," else "\t"
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, open = "rt") else file(path, open = "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    utils::read.delim(con, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
}

.xenium_codebook_resolve_map <- function(input_dir, verbose = FALSE) {
    path <- .resolve_xenium_codebook_path(input_dir)
    if (is.null(path)) {
        return(list(map = NULL, path = NULL, index_column = NULL, name_column = NULL))
    }
    tb <- tryCatch(.read_delim_any(path, verbose = verbose), error = function(e) NULL)
    if (is.null(tb) || !nrow(tb)) {
        return(list(map = NULL, path = path, index_column = NULL, name_column = NULL))
    }

    cn <- tolower(names(tb))
    idx_candidates <- c("codeword_index", "codewordindex", "codeword_idx", "codeword", "codeword_id", "codewordid")
    name_candidates <- c("feature_name", "gene", "gene_name", "target", "target_name", "name", "display_name")
    idx_col <- names(tb)[match(TRUE, cn %in% idx_candidates)]
    name_col <- names(tb)[match(TRUE, cn %in% name_candidates)]
    if (!length(idx_col) || is.na(idx_col) || !length(name_col) || is.na(name_col)) {
        return(list(map = NULL, path = path, index_column = idx_col, name_column = name_col))
    }

    idx <- as.character(tb[[idx_col]])
    nm <- trimws(as.character(tb[[name_col]]))
    ok <- !is.na(idx) & nzchar(idx) & !is.na(nm) & nzchar(nm)
    if (!any(ok)) {
        return(list(map = NULL, path = path, index_column = idx_col, name_column = name_col))
    }
    map <- setNames(nm[ok], idx[ok])
    if (anyDuplicated(names(map))) {
        map <- map[!duplicated(names(map))]
    }
    list(map = map, path = path, index_column = idx_col, name_column = name_col)
}

.transcripts_spec_default <- function() {
    list(
        feature = list(
            mode = "column",
            column = "feature_name",
            is_gene_column = NULL,
            map = NULL,
            map_path = NULL
        )
    )
}

.resolve_transcripts_feature_spec <- function(ds, input_dir = NULL, verbose = FALSE) {
    cols <- tryCatch(names(ds), error = function(e) character())
    if (!length(cols)) cols <- character()
    spec <- .transcripts_spec_default()

    if ("is_gene" %in% cols) {
        spec$feature$is_gene_column <- "is_gene"
    }

    candidates <- c("feature_name", "gene", "gene_name", "target", "target_name", "codeword_name")
    hit <- candidates[candidates %in% cols]
    if (length(hit)) {
        spec$feature$mode <- "column"
        spec$feature$column <- hit[[1]]
        return(spec)
    }

    if ("codeword_index" %in% cols) {
        spec$feature$mode <- "codeword_index"
        spec$feature$column <- "codeword_index"
        cb <- .xenium_codebook_resolve_map(input_dir, verbose = verbose)
        spec$feature$map <- cb$map
        spec$feature$map_path <- cb$path
        return(spec)
    }

    spec$feature$mode <- "none"
    spec$feature$column <- NULL
    spec
}

.transcripts_feature_project_cols <- function(spec) {
    if (is.null(spec) || !is.list(spec)) spec <- .transcripts_spec_default()
    feat <- spec$feature
    cols <- character()
    if (!is.null(feat$column) && nzchar(feat$column)) cols <- c(cols, feat$column)
    if (!is.null(feat$is_gene_column) && nzchar(feat$is_gene_column)) cols <- c(cols, feat$is_gene_column)
    unique(cols)
}

.transcripts_extract_feature_name <- function(dt, spec) {
    if (is.null(spec) || !is.list(spec)) spec <- .transcripts_spec_default()
    feat <- spec$feature
    mode <- feat$mode
    col <- feat$column
    if (is.null(col) || !nzchar(col) || !col %in% names(dt)) return(NULL)

    if (identical(mode, "column")) {
        return(as.character(dt[[col]]))
    }
    if (identical(mode, "codeword_index")) {
        key <- as.character(dt[[col]])
        if (!is.null(feat$map)) {
            return(as.character(feat$map[key]))
        }
        return(as.character(key))
    }
    NULL
}

#' Ingest Molecule Filter Mask
#' @description
#' Internal helper to apply Xenium/CosMx molecule filters to an in-memory
#' `data.frame`/`data.table` chunk while preserving `dplyr::filter()` semantics
#' (i.e. NA predicates drop rows).
#' @param dt Chunk with transcript columns.
#' @param p Parameter list from `.initialize_dataset_builder_state()`.
#' @param qc_genes Optional character vector of QC-flagged genes to exclude.
#' @return Logical vector mask.
#' @keywords internal
.ingest_molecule_filter_mask <- function(dt, p, qc_genes = character()) {
    n <- nrow(dt)
    if (!is.finite(n) || n <= 0L) return(logical(0))
    mask <- rep(TRUE, n)

    spec <- p$transcripts_spec
    if (is.null(spec) || !is.list(spec)) spec <- .transcripts_spec_default()
    fn <- .transcripts_extract_feature_name(dt, spec)

    if (!is.null(p$filter_genes)) {
        if (is.null(fn)) stop("Missing transcript feature column; cannot apply filter_genes.")
        genes <- as.character(p$filter_genes)
        genes <- genes[!is.na(genes) & nzchar(genes)]
        if (!length(genes)) genes <- character()
        mask <- mask & !is.na(fn) & fn %in% genes
    }
    if (isTRUE(p$exclude_qc_unsuitable) && length(qc_genes)) {
        if (is.null(fn)) stop("Missing transcript feature column; cannot apply exclude_qc_unsuitable.")
        qc_genes <- as.character(qc_genes)
        qc_genes <- qc_genes[!is.na(qc_genes) & nzchar(qc_genes)]
        if (!length(qc_genes)) qc_genes <- character()
        mask <- mask & !is.na(fn) & !(fn %in% qc_genes)
    }
    if (isTRUE(p$filtermolecule)) {
        feat <- spec$feature
        has_name_map <- !identical(feat$mode, "codeword_index") || !is.null(feat$map)
        if (is.null(fn) || !isTRUE(has_name_map)) {
            ig_col <- feat$is_gene_column
            if (!is.null(ig_col) && nzchar(ig_col) && ig_col %in% names(dt)) {
                ig <- dt[[ig_col]]
                mask <- mask & !is.na(ig) & as.logical(ig)
            } else if (is.null(fn)) {
                stop("Missing transcript feature column; cannot apply exclude_prefix.")
            }
        } else {
            mask <- mask & !is.na(fn)
            prefixes <- p$exclude_prefix
            if (length(prefixes)) {
                prefixes <- as.character(prefixes)
                prefixes <- prefixes[!is.na(prefixes) & nzchar(prefixes)]
            }
            if (length(prefixes)) {
                for (pref in prefixes) {
                    mask <- mask & !startsWith(fn, pref)
                }
            }
        }
    }
    if (!is.null(p$max_dist_mol_nuc)) {
        if (!"nucleus_distance" %in% names(dt)) stop("Missing `nucleus_distance` column; cannot apply max_dist_mol_nuc.")
        nd <- dt[["nucleus_distance"]]
        mask <- mask & !is.na(nd) & nd <= p$max_dist_mol_nuc
    }
    if (!is.null(p$filterqv)) {
        if (!"qv" %in% names(dt)) stop("Missing `qv` column; cannot apply filterqv.")
        qv <- dt[["qv"]]
        mask <- mask & !is.na(qv) & qv >= p$filterqv
    }

    if (!is.null(fn) && length(fn) == n) {
        attr(mask, "feature_name") <- fn
    }
    mask
}

.resolve_transcripts_csv_path <- function(input_dir) {
    if (is.null(input_dir) || !length(input_dir)) return(NULL)
    input_dir <- as.character(input_dir)[1]
    if (is.na(input_dir) || !nzchar(input_dir)) return(NULL)
    candidates <- c(
        file.path(input_dir, "transcripts.csv.gz"),
        file.path(input_dir, "transcripts.csv")
    )
    hit <- candidates[file.exists(candidates)]
    if (!length(hit)) return(NULL)
    hit[[1]]
}

.resolve_transcripts_ingest_cache_dir <- function(input_dir) {
    if (is.null(input_dir) || !length(input_dir)) return(NULL)
    input_dir <- as.character(input_dir)[1]
    if (is.na(input_dir) || !nzchar(input_dir)) return(NULL)
    file.path(input_dir, "analysis", "genescope", "ingest_cache")
}

.resolve_pyarrow_python <- function(verbose = FALSE) {
    candidates <- c(
        Sys.getenv("GENESCOPE_PYTHON", unset = ""),
        Sys.getenv("RETICULATE_PYTHON", unset = ""),
        Sys.which("python3"),
        Sys.which("python")
    )
    if (requireNamespace("reticulate", quietly = TRUE)) {
        py_discover <- tryCatch(getFromNamespace("py_discover_config", "reticulate"), error = function(e) NULL)
        if (is.function(py_discover)) {
            cfg <- tryCatch(py_discover(required_module = "pyarrow"), error = function(e) NULL)
            cfg_py <- if (is.list(cfg) && !is.null(cfg$python)) as.character(cfg$python)[1] else ""
            if (!is.na(cfg_py) && nzchar(cfg_py)) {
                candidates <- c(candidates, cfg_py)
            }
        }
    }
    candidates <- unique(as.character(candidates))
    candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
    script_path <- tempfile("genescope_pyarrow_probe_", fileext = ".py")
    writeLines("import pyarrow", con = script_path, useBytes = TRUE)
    on.exit(unlink(script_path), add = TRUE)
    for (cand in candidates) {
        status <- tryCatch(
            suppressWarnings(system2(cand, script_path, stdout = FALSE, stderr = FALSE)),
            error = function(e) 1L
        )
        status <- suppressWarnings(as.integer(status)[1])
        if (is.finite(status) && identical(status, 0L)) {
            return(cand)
        }
    }
    if (isTRUE(verbose)) {
        message("No Python interpreter with pyarrow was found; skipping automatic Xenium transcript conversion.")
    }
    NULL
}

.transcripts_fallback_auth_key_id <- "genescope-transcripts-fallback-ed25519-20260424"

.transcripts_fallback_auth_public_key_pem <- paste(
    c(
        "-----BEGIN PUBLIC KEY-----",
        "MCowBQYDK2VwAyEAOkmCEEB0sxletJt/Az+dEFX1NAskbtlSosor6Y8ErsI=",
        "-----END PUBLIC KEY-----"
    ),
    collapse = "\n"
)

.transcripts_fallback_script_path <- function() {
    path <- tryCatch(
        system.file("python", "xenium_transcripts_fallback.py", package = "geneSCOPE"),
        error = function(e) ""
    )
    if (is.character(path) && length(path) && nzchar(path) && file.exists(path)) {
        return(path)
    }
    NULL
}

.transcripts_fallback_script_signature_path <- function() {
    path <- tryCatch(
        system.file("python", "xenium_transcripts_fallback.py.sig", package = "geneSCOPE"),
        error = function(e) ""
    )
    if (is.character(path) && length(path) && nzchar(path) && file.exists(path)) {
        return(path)
    }
    NULL
}

.normalize_transcripts_attestation_path <- function(path) {
    if (is.null(path) || !length(path)) return(NA_character_)
    path <- as.character(path)[1]
    if (is.na(path) || !nzchar(path)) return(NA_character_)
    normalizePath(path, winslash = "/", mustWork = FALSE)
}

.transcripts_attestation_scalar <- function(value, lowercase = FALSE) {
    if (is.null(value) || !length(value)) return(NA_character_)
    value <- as.character(value)[1]
    if (is.na(value)) return(NA_character_)
    value <- trimws(value)
    if (!nzchar(value)) return(NA_character_)
    if (isTRUE(lowercase)) value <- tolower(value)
    value
}

.transcripts_config_value <- function(option_name, env_name, default = NULL) {
    opt <- getOption(option_name, NULL)
    if (!is.null(opt) && length(opt)) {
        opt <- .transcripts_attestation_scalar(opt)
        if (!is.na(opt)) return(opt)
    }
    env <- Sys.getenv(env_name, unset = "")
    if (nzchar(env)) return(trimws(env))
    default
}

.transcripts_resolve_path_config <- function(option_name, env_name, default = NULL) {
    path <- .transcripts_config_value(option_name, env_name, default = default)
    if (is.null(path) || !length(path)) return(NULL)
    path <- as.character(path)[1]
    if (is.na(path) || !nzchar(path)) return(NULL)
    .normalize_transcripts_attestation_path(path)
}

.transcripts_integrity_sidecar_path <- function(path) {
    paste0(path, ".sha256")
}

.transcripts_fallback_auth_public_key_path <- function() {
    default_path <- tryCatch(
        system.file("trust", "transcripts-fallback-root.pem", package = "geneSCOPE"),
        error = function(e) ""
    )
    default_path <- if (is.character(default_path) && length(default_path) && nzchar(default_path) &&
        file.exists(default_path)) {
        default_path
    } else {
        NULL
    }
    .transcripts_resolve_path_config(
        "geneSCOPE.transcripts.fallback_trust_root_path",
        "GENESCOPE_TRANSCRIPTS_FALLBACK_TRUST_ROOT_PATH",
        default = default_path
    )
}

.transcripts_runtime_profile_path <- function() {
    .transcripts_resolve_path_config(
        "geneSCOPE.transcripts.runtime_profile_path",
        "GENESCOPE_TRANSCRIPTS_RUNTIME_PROFILE_PATH"
    )
}

.transcripts_runtime_profile_trust_root_path <- function() {
    .transcripts_resolve_path_config(
        "geneSCOPE.transcripts.runtime_profile_trust_root_path",
        "GENESCOPE_TRANSCRIPTS_RUNTIME_PROFILE_TRUST_ROOT_PATH",
        default = .transcripts_fallback_auth_public_key_path()
    )
}

.transcripts_runtime_profile_signature_path <- function(path) {
    if (is.null(path) || !length(path)) return(NULL)
    path <- as.character(path)[1]
    if (is.na(path) || !nzchar(path)) return(NULL)
    paste0(path, ".sig")
}

.transcripts_signature_trust_mode <- function() {
    .trust_chain_default_mode(
        .transcripts_config_value(
            "geneSCOPE.transcripts.trust_mode",
            "GENESCOPE_TRANSCRIPTS_TRUST_MODE",
            default = "ed25519"
        )
    )
}

.transcripts_runtime_trust_mode <- function(runtime_profile_path = .transcripts_runtime_profile_path()) {
    configured <- .transcripts_config_value(
        "geneSCOPE.transcripts.runtime_trust_mode",
        "GENESCOPE_TRANSCRIPTS_RUNTIME_TRUST_MODE",
        default = NULL
    )
    allowed <- c("observed_runtime", "strict_runtime_profile", "signed_runtime_profile")
    if (!is.null(configured)) {
        configured <- .transcripts_attestation_scalar(configured)
        if (!is.na(configured) && configured %in% allowed) {
            return(configured)
        }
        stop(
            "Unsupported transcript runtime trust mode '", configured,
            "'. Use one of: ", paste(allowed, collapse = ", "),
            ".",
            call. = FALSE
        )
    }
    if (!is.null(runtime_profile_path) && !is.na(runtime_profile_path) && nzchar(runtime_profile_path)) {
        return("signed_runtime_profile")
    }
    "observed_runtime"
}

.read_transcripts_runtime_profile <- function(path) {
    if (is.null(path) || !length(path) || !nzchar(path)) return(NULL)
    if (!file.exists(path)) {
        stop("Transcript runtime profile not found: ", path)
    }

    profile_to_list <- function(x) {
        if (is.null(x)) return(NULL)
        if (is.data.frame(x)) {
            if (nrow(x) < 1L) return(NULL)
            x <- as.list(x[1, , drop = TRUE])
        }
        if (!is.list(x) || is.null(names(x))) return(NULL)
        x
    }

    ext <- tolower(tools::file_ext(path))
    profile <- NULL
    if (identical(ext, "json")) {
        if (!requireNamespace("jsonlite", quietly = TRUE)) {
            stop("Package 'jsonlite' is required to read transcript runtime profiles.")
        }
        profile <- tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) e)
        if (inherits(profile, "error")) {
            stop("Failed to parse transcript runtime profile JSON: ", conditionMessage(profile))
        }
        profile <- profile_to_list(profile)
    } else {
        lines <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
        if (!length(lines)) {
            stop("Transcript runtime profile is empty: ", path)
        }
        lines <- trimws(lines)
        lines <- lines[nzchar(lines) & !grepl("^[[:space:]]*#", lines)]
        if (!length(lines)) {
            stop("Transcript runtime profile contains no usable entries: ", path)
        }
        keys <- character()
        vals <- character()
        for (line in lines) {
            parts <- strsplit(line, "=", fixed = TRUE)[[1]]
            if (length(parts) < 2L) next
            key <- trimws(parts[[1]])
            val <- trimws(paste(parts[-1], collapse = "="))
            if (!nzchar(key) || !nzchar(val)) next
            keys <- c(keys, key)
            vals <- c(vals, val)
        }
        if (!length(keys)) {
            stop("Transcript runtime profile contains no usable key=value entries: ", path)
        }
        profile <- as.list(vals)
        names(profile) <- keys
    }

    if (is.null(profile)) {
        stop("Transcript runtime profile contains no recognized data: ", path)
    }

    canon <- list()
    for (nm in names(profile)) {
        canon[[tolower(trimws(nm))]] <- .transcripts_attestation_scalar(profile[[nm]])
    }
    out <- list(
        python_sha256 = .transcripts_attestation_scalar(canon[["python_sha256"]], lowercase = TRUE),
        pyarrow_version = .transcripts_attestation_scalar(canon[["pyarrow_version"]]),
        openssl_sha256 = .transcripts_attestation_scalar(canon[["openssl_sha256"]], lowercase = TRUE),
        path = .normalize_transcripts_attestation_path(path),
        signature_path = .normalize_transcripts_attestation_path(.transcripts_runtime_profile_signature_path(path)),
        sha256 = .compute_file_sha256(path),
        trust_level = .transcripts_runtime_trust_mode(path)
    )
    if (is.na(out$python_sha256) && is.na(out$pyarrow_version) && is.na(out$openssl_sha256)) {
        stop(
            "Transcript runtime profile must define at least one of python_sha256, pyarrow_version, or openssl_sha256: ",
            path
        )
    }
    out
}

.probe_transcripts_runtime_openssl <- function() {
    openssl_bin <- Sys.which("openssl")
    if (!nzchar(openssl_bin)) {
        return(list(
            openssl_executable = NA_character_,
            openssl_version = NA_character_,
            openssl_sha256 = NA_character_
        ))
    }

    version_out <- tryCatch(
        suppressWarnings(system2(openssl_bin, "--version", stdout = TRUE, stderr = TRUE)),
        error = function(e) e
    )
    openssl_version <- NA_character_
    if (!inherits(version_out, "error") && length(version_out)) {
        openssl_version <- trimws(version_out[[1]])
    }

    list(
        openssl_executable = .normalize_transcripts_attestation_path(openssl_bin),
        openssl_version = openssl_version,
        openssl_sha256 = .compute_file_sha256(openssl_bin)
    )
}

.compare_transcripts_runtime_attestation <- function(runtime_attestation, runtime_profile = NULL) {
    if (!is.list(runtime_attestation) || !length(runtime_attestation)) {
        return(list(ok = FALSE, reason = "missing_runtime_attestation"))
    }

    rt <- function(name) {
        .transcripts_attestation_scalar(runtime_attestation[[name]])
    }
    rt_hex <- function(name) {
        .transcripts_attestation_scalar(runtime_attestation[[name]], lowercase = TRUE)
    }

    base_fields <- c(
        "runtime_trust_level",
        "runtime_python_version",
        "runtime_python_executable",
        "runtime_pyarrow_version"
    )
    missing_base <- base_fields[vapply(base_fields, function(field) is.na(rt(field)), logical(1))]
    if (length(missing_base)) {
        return(list(ok = FALSE, reason = paste0("missing_runtime_fields:", paste(missing_base, collapse = ","))))
    }

    trust_level <- rt("runtime_trust_level")
    if (!trust_level %in% c("observed_runtime", "strict_runtime_profile", "signed_runtime_profile")) {
        return(list(ok = FALSE, reason = paste0("invalid_runtime_trust_level:", trust_level)))
    }

    if (!is.null(runtime_profile) && is.list(runtime_profile) && length(runtime_profile)) {
        profile_trust_level <- .transcripts_attestation_scalar(runtime_profile$trust_level)
        if (!is.na(profile_trust_level) && !identical(trust_level, profile_trust_level)) {
            return(list(
                ok = FALSE,
                reason = paste0("runtime_trust_level_mismatch:", trust_level, "!=", profile_trust_level)
            ))
        }
        profile_mismatches <- character()
        profile_python_sha256 <- .transcripts_attestation_scalar(runtime_profile$python_sha256, lowercase = TRUE)
        if (!is.na(profile_python_sha256)) {
            if (!identical(rt_hex("runtime_python_executable_sha256"), profile_python_sha256)) {
                profile_mismatches <- c(profile_mismatches, "python_sha256")
            }
        }
        profile_pyarrow_version <- .transcripts_attestation_scalar(runtime_profile$pyarrow_version)
        if (!is.na(profile_pyarrow_version)) {
            if (!identical(rt("runtime_pyarrow_version"), profile_pyarrow_version)) {
                profile_mismatches <- c(profile_mismatches, "pyarrow_version")
            }
        }
        profile_openssl_sha256 <- .transcripts_attestation_scalar(runtime_profile$openssl_sha256, lowercase = TRUE)
        if (!is.na(profile_openssl_sha256)) {
            if (!identical(rt_hex("runtime_openssl_sha256"), profile_openssl_sha256)) {
                profile_mismatches <- c(profile_mismatches, "openssl_sha256")
            }
        }
        if (length(profile_mismatches)) {
            return(list(
                ok = FALSE,
                reason = paste0("strict_runtime_profile_mismatch:", paste(profile_mismatches, collapse = ","))
            ))
        }
        if (identical(profile_trust_level, "signed_runtime_profile")) {
            if (!identical(rt("runtime_profile_signature_verified"), "TRUE")) {
                return(list(ok = FALSE, reason = "missing_runtime_profile_signature_verification"))
            }
            profile_signature_sha256 <- .transcripts_attestation_scalar(runtime_profile$signature_sha256, lowercase = TRUE)
            if (!is.na(profile_signature_sha256) &&
                !identical(rt_hex("runtime_profile_signature_sha256"), profile_signature_sha256)) {
                return(list(ok = FALSE, reason = "runtime_profile_signature_sha256_mismatch"))
            }
        }
    }

    list(ok = TRUE, reason = "ok")
}

.attach_transcripts_runtime_attestation <- function(runtime_attestation, runtime_profile_path = NULL) {
    runtime_profile <- NULL
    runtime_profile_path <- .normalize_transcripts_attestation_path(runtime_profile_path)
    if (!is.na(runtime_profile_path) && nzchar(runtime_profile_path)) {
        runtime_profile <- .read_transcripts_runtime_profile(runtime_profile_path)
    }
    trust_mode <- .transcripts_runtime_trust_mode(runtime_profile_path)
    if (identical(trust_mode, "observed_runtime") && !is.null(runtime_profile)) {
        trust_mode <- runtime_profile$trust_level
    }
    profile_signature <- NULL
    if (!is.null(runtime_profile)) {
        profile_signature <- .verify_transcripts_runtime_profile_signature(
            runtime_profile$path,
            require_signature = identical(trust_mode, "signed_runtime_profile")
        )
        if (identical(trust_mode, "signed_runtime_profile") && !isTRUE(profile_signature$ok)) {
            stop("Transcript runtime profile signature verification failed: ", profile_signature$reason)
        }
    } else if (!identical(trust_mode, "observed_runtime")) {
        stop("Transcript runtime trust mode '", trust_mode, "' requires GENESCOPE_TRANSCRIPTS_RUNTIME_PROFILE_PATH.")
    }

    runtime <- list()
    runtime$runtime_trust_level <- trust_mode
    runtime$runtime_profile_path <- if (!is.null(runtime_profile)) runtime_profile$path else NA_character_
    runtime$runtime_profile_sha256 <- if (!is.null(runtime_profile)) runtime_profile$sha256 else NA_character_
    runtime$runtime_profile_signature_path <- if (!is.null(runtime_profile)) runtime_profile$signature_path else NA_character_
    runtime$runtime_profile_signature_sha256 <- if (is.list(profile_signature)) profile_signature$signature_sha256 else NA_character_
    runtime$runtime_profile_signature_verified <- if (is.list(profile_signature) && isTRUE(profile_signature$ok)) "TRUE" else "FALSE"
    runtime$runtime_profile_signature_key_id <- if (is.list(profile_signature)) profile_signature$signature_key_id else NA_character_
    runtime$runtime_profile_trust_root_sha256 <- if (is.list(profile_signature)) profile_signature$trust_root_sha256 else NA_character_
    runtime$runtime_python_version <- .transcripts_attestation_scalar(runtime_attestation$python_version)
    runtime$runtime_python_executable <- .normalize_transcripts_attestation_path(runtime_attestation$python_executable)
    runtime$runtime_python_executable_sha256 <- .transcripts_attestation_scalar(
        runtime_attestation$python_executable_sha256,
        lowercase = TRUE
    )
    runtime$runtime_pyarrow_version <- .transcripts_attestation_scalar(runtime_attestation$pyarrow_version)
    runtime$runtime_pyarrow_module_path <- .normalize_transcripts_attestation_path(runtime_attestation$pyarrow_module_path)
    runtime$runtime_openssl_executable <- .normalize_transcripts_attestation_path(runtime_attestation$openssl_executable)
    runtime$runtime_openssl_version <- .transcripts_attestation_scalar(runtime_attestation$openssl_version)
    runtime$runtime_openssl_sha256 <- .transcripts_attestation_scalar(runtime_attestation$openssl_sha256, lowercase = TRUE)

    if (!is.null(runtime_profile)) {
        checked <- .compare_transcripts_runtime_attestation(runtime, runtime_profile = runtime_profile)
        if (!isTRUE(checked$ok)) {
            stop("Transcript runtime profile mismatch: ", checked$reason)
        }
    }

    runtime
}

.run_transcripts_auth_openssl <- function(args) {
    openssl_bin <- Sys.which("openssl")
    if (!nzchar(openssl_bin)) {
        stop("OpenSSL is required for transcript fallback authenticity signing.")
    }
    out <- tryCatch(
        suppressWarnings(system2(openssl_bin, args, stdout = TRUE, stderr = TRUE)),
        error = function(e) e
    )
    if (inherits(out, "error")) {
        stop("OpenSSL invocation failed: ", conditionMessage(out))
    }
    status <- suppressWarnings(as.integer(attr(out, "status"))[1])
    if (!length(status) || !is.finite(status)) status <- 0L
    if (!identical(status, 0L)) {
        msg <- paste(out, collapse = "\n")
        if (!nzchar(msg)) msg <- paste0("openssl exited with status ", status)
        stop(msg)
    }
    out
}

.decode_transcripts_signature_file <- function(sig_path) {
    if (is.null(sig_path) || !nzchar(sig_path) || !file.exists(sig_path)) {
        stop("Transcript fallback signature file not found: ", sig_path)
    }
    sig_bin_path <- tempfile("genescope_transcripts_sig_", fileext = ".bin")
    .run_transcripts_auth_openssl(c(
        "base64", "-d", "-A",
        "-in", sig_path,
        "-out", sig_bin_path
    ))
    sig_info <- tryCatch(file.info(sig_bin_path), error = function(e) NULL)
    sig_size <- if (!is.null(sig_info)) sig_info$size[[1]] else NA_real_
    if (!file.exists(sig_bin_path) || !is.finite(sig_size) || sig_size <= 0) {
        stop("Transcript fallback signature payload could not be decoded.")
    }
    sig_bin_path
}

.verify_transcripts_fallback_script_signature <- function() {
    script_path <- .transcripts_fallback_script_path()
    if (is.null(script_path) || !file.exists(script_path)) {
        return(list(ok = FALSE, reason = "missing_script", script_path = script_path))
    }

    sig_path <- .transcripts_fallback_script_signature_path()
    if (is.null(sig_path) || !file.exists(sig_path)) {
        return(list(ok = FALSE, reason = "missing_script_signature", script_path = script_path))
    }

    verified <- .verify_ed25519_file_signature(
        data_path = script_path,
        sig_path = sig_path,
        trust_mode = .transcripts_signature_trust_mode(),
        trust_root_path = .transcripts_fallback_auth_public_key_path(),
        bundled_public_key_pem = .transcripts_fallback_auth_public_key_pem,
        key_id = .transcripts_fallback_auth_key_id
    )
    if (!isTRUE(verified$ok)) {
        return(list(ok = FALSE, reason = verified$reason, script_path = script_path))
    }

    list(
        ok = TRUE,
        reason = "ok",
        script_path = script_path,
        script_sha256 = .compute_file_sha256(script_path),
        signature_key_id = .transcripts_fallback_auth_key_id,
        script_trust_root_path = verified$trust_root_path,
        script_trust_root_sha256 = verified$trust_root_sha256,
        trust_mode = verified$trust_mode
    )
}

.verify_transcripts_runtime_profile_signature <- function(profile_path, require_signature = FALSE) {
    profile_path <- .normalize_transcripts_attestation_path(profile_path)
    if (is.null(profile_path) || !nzchar(profile_path) || !file.exists(profile_path)) {
        return(list(ok = FALSE, reason = "missing_runtime_profile", profile_path = profile_path))
    }

    sig_path <- .transcripts_runtime_profile_signature_path(profile_path)
    if (is.null(sig_path) || !file.exists(sig_path)) {
        if (isTRUE(require_signature)) {
            return(list(ok = FALSE, reason = "missing_runtime_profile_signature", profile_path = profile_path))
        }
        return(list(
            ok = FALSE,
            reason = "unsigned_runtime_profile",
            profile_path = profile_path,
            signature_path = .normalize_transcripts_attestation_path(sig_path)
        ))
    }

    verified <- .verify_ed25519_file_signature(
        data_path = profile_path,
        sig_path = sig_path,
        trust_mode = .transcripts_signature_trust_mode(),
        trust_root_path = .transcripts_runtime_profile_trust_root_path(),
        bundled_public_key_pem = .transcripts_fallback_auth_public_key_pem,
        key_id = .transcripts_fallback_auth_key_id
    )
    if (!isTRUE(verified$ok)) {
        return(list(ok = FALSE, reason = verified$reason, profile_path = profile_path))
    }

    list(
        ok = TRUE,
        reason = "ok",
        profile_path = profile_path,
        signature_path = .normalize_transcripts_attestation_path(sig_path),
        profile_sha256 = .compute_file_sha256(profile_path),
        signature_sha256 = verified$signature_sha256,
        signature_key_id = .transcripts_fallback_auth_key_id,
        trust_root_sha256 = verified$trust_root_sha256,
        trust_mode = verified$trust_mode
    )
}

.build_transcripts_integrity_report <- function(src_path,
                                               dst_path,
                                               mode,
                                               python_bin,
                                               src_sha256,
                                               src_size,
                                               dst_sha256,
                                               dst_size,
                                               timeout_secs,
                                               runtime_attestation = NULL) {
    script_check <- .verify_transcripts_fallback_script_signature()
    if (!isTRUE(script_check$ok)) {
        stop("Trusted transcript fallback script failed authenticity verification: ", script_check$reason)
    }

    if (is.null(runtime_attestation)) {
        runtime_attestation <- .probe_pyarrow_runtime(python_bin, verbose = FALSE)
    }
    if (is.list(runtime_attestation) && !is.null(runtime_attestation$runtime_trust_level)) {
        runtime <- runtime_attestation
    } else {
        runtime <- .attach_transcripts_runtime_attestation(
            runtime_attestation = runtime_attestation,
            runtime_profile_path = .transcripts_runtime_profile_path()
        )
    }

    list(
        attestation_type = "genescope_transcript_fallback_v2",
        attestation_version = "2",
        signature_key_id = script_check$signature_key_id,
        script_sha256 = script_check$script_sha256,
        script_trust_root_path = script_check$script_trust_root_path,
        script_trust_root_sha256 = script_check$script_trust_root_sha256,
        source_path = .normalize_transcripts_attestation_path(src_path),
        dest_path = .normalize_transcripts_attestation_path(dst_path),
        python_bin = as.character(python_bin)[1],
        mode = as.character(mode)[1],
        source_sha256 = src_sha256,
        source_size = src_size,
        dest_sha256 = dst_sha256,
        dest_size = dst_size,
        timeout_secs = timeout_secs,
        runtime_trust_level = runtime$runtime_trust_level,
        runtime_profile_path = runtime$runtime_profile_path,
        runtime_profile_sha256 = runtime$runtime_profile_sha256,
        runtime_profile_signature_path = runtime$runtime_profile_signature_path,
        runtime_profile_signature_sha256 = runtime$runtime_profile_signature_sha256,
        runtime_profile_signature_verified = runtime$runtime_profile_signature_verified,
        runtime_profile_signature_key_id = runtime$runtime_profile_signature_key_id,
        runtime_profile_trust_root_sha256 = runtime$runtime_profile_trust_root_sha256,
        runtime_python_version = runtime$runtime_python_version,
        runtime_python_executable = runtime$runtime_python_executable,
        runtime_python_executable_sha256 = runtime$runtime_python_executable_sha256,
        runtime_pyarrow_version = runtime$runtime_pyarrow_version,
        runtime_pyarrow_module_path = runtime$runtime_pyarrow_module_path,
        runtime_openssl_executable = runtime$runtime_openssl_executable,
        runtime_openssl_version = runtime$runtime_openssl_version,
        runtime_openssl_sha256 = runtime$runtime_openssl_sha256
    )
}

.check_transcripts_integrity_provenance <- function(report, src_path, dst_path) {
    if (!is.list(report) || !length(report)) {
        return(list(ok = FALSE, reason = "missing_attestation_report"))
    }

    script_check <- .verify_transcripts_fallback_script_signature()
    if (!isTRUE(script_check$ok)) {
        return(list(ok = FALSE, reason = paste0("script_authenticity_failed:", script_check$reason)))
    }

    report1 <- function(name) {
        val <- report[[name]]
        if (is.null(val) || !length(val)) return(NA_character_)
        as.character(val)[1]
    }

    checks <- c(
        attestation_type = identical(report1("attestation_type"), "genescope_transcript_fallback_v2"),
        attestation_version = identical(report1("attestation_version"), "2"),
        signature_key_id = identical(report1("signature_key_id"), script_check$signature_key_id),
        script_sha256 = identical(report1("script_sha256"), script_check$script_sha256),
        script_trust_root_sha256 = identical(report1("script_trust_root_sha256"), script_check$script_trust_root_sha256),
        source_path = identical(report1("source_path"), .normalize_transcripts_attestation_path(src_path)),
        dest_path = identical(report1("dest_path"), .normalize_transcripts_attestation_path(dst_path))
    )

    if (!all(checks)) {
        bad <- names(checks)[!checks]
        return(list(ok = FALSE, reason = paste0("provenance_mismatch:", paste(bad, collapse = ","))))
    }

    runtime_profile_path <- .transcripts_runtime_profile_path()
    runtime_profile <- NULL
    if (!is.null(runtime_profile_path) && nzchar(runtime_profile_path)) {
        runtime_profile <- tryCatch(
            .read_transcripts_runtime_profile(runtime_profile_path),
            error = function(e) e
        )
        if (inherits(runtime_profile, "error")) {
            return(list(ok = FALSE, reason = paste0("runtime_profile_error:", conditionMessage(runtime_profile))))
        }
    }

    runtime_check <- .compare_transcripts_runtime_attestation(
        report,
        runtime_profile = runtime_profile
    )
    if (!isTRUE(runtime_check$ok)) {
        return(list(ok = FALSE, reason = runtime_check$reason))
    }

    list(ok = TRUE, reason = "ok")
}
.compute_file_sha256 <- function(path) {
    if (is.null(path) || !length(path)) return(NA_character_)
    path <- as.character(path)[1]
    if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_character_)

    if (requireNamespace("openssl", quietly = TRUE)) {
        con <- file(path, open = "rb")
        on.exit(close(con), add = TRUE)
        hash_raw <- tryCatch(openssl::sha256(con), error = function(e) NULL)
        if (is.raw(hash_raw) && length(hash_raw)) {
            return(paste(sprintf("%02x", as.integer(hash_raw)), collapse = ""))
        }
    }

    if (requireNamespace("digest", quietly = TRUE)) {
        hash_val <- tryCatch(
            digest::digest(file = path, algo = "sha256", serialize = FALSE),
            error = function(e) NA_character_
        )
        hash_val <- unname(as.character(hash_val)[1])
        if (!is.na(hash_val) && nzchar(hash_val)) {
            return(hash_val)
        }
    }

    commands <- list(c("shasum", "-a", "256"), c("sha256sum"))
    for (cmd in commands) {
        exe <- Sys.which(cmd[[1]])
        if (!nzchar(exe)) next
        out <- tryCatch(
            system2(exe, c(cmd[-1], path), stdout = TRUE, stderr = FALSE),
            error = function(e) character()
        )
        status <- suppressWarnings(as.integer(attr(out, "status"))[1])
        if (length(status) && is.finite(status) && !identical(status, 0L)) next
        if (!length(out)) next
        fields <- strsplit(trimws(out[[1]]), "[[:space:]]+")[[1]]
        fields <- fields[nzchar(fields)]
        if (length(fields)) {
            return(unname(fields[[1]]))
        }
    }

    NA_character_
}

.read_transcripts_integrity_sidecar <- function(path) {
    report_path <- .transcripts_integrity_sidecar_path(path)
    if (!file.exists(report_path)) return(NULL)
    lines <- tryCatch(readLines(report_path, warn = FALSE), error = function(e) character())
    if (!length(lines)) return(NULL)

    parsed <- strsplit(lines, "=", fixed = TRUE)
    keys <- vapply(parsed, function(x) if (length(x)) x[[1]] else "", character(1))
    vals <- vapply(
        parsed,
        function(x) if (length(x) >= 2L) paste(x[-1], collapse = "=") else "",
        character(1)
    )
    keep <- nzchar(keys)
    if (!any(keep)) return(NULL)

    out <- as.list(vals[keep])
    names(out) <- keys[keep]
    out
}

.write_transcripts_integrity_sidecar <- function(path, report) {
    if (!is.list(report) || !length(report)) return(invisible(NULL))
    report_path <- .transcripts_integrity_sidecar_path(path)
    fields <- vapply(report, function(x) as.character(x)[1], character(1))
    keep <- !is.na(fields) & nzchar(fields)
    if (!any(keep)) return(invisible(NULL))
    writeLines(paste0(names(fields)[keep], "=", unname(fields[keep])), report_path, useBytes = TRUE)
    invisible(report_path)
}

.should_reuse_transcripts_conversion <- function(src_path, dst_path, verbose = FALSE) {
    if (!file.exists(src_path) || !file.exists(dst_path)) return(FALSE)

    src_info <- tryCatch(file.info(src_path), error = function(e) NULL)
    dst_info <- tryCatch(file.info(dst_path), error = function(e) NULL)
    if (is.null(src_info) || is.null(dst_info)) return(FALSE)

    src_mtime <- src_info$mtime[[1]]
    dst_mtime <- dst_info$mtime[[1]]
    dst_size <- dst_info$size[[1]]
    basic_ok <- is.finite(dst_size) && dst_size > 0 &&
        !is.na(src_mtime) && !is.na(dst_mtime) && dst_mtime >= src_mtime
    if (!basic_ok) return(FALSE)

    report <- .read_transcripts_integrity_sidecar(dst_path)
    if (!is.list(report) || is.null(report$source_sha256) || is.null(report$dest_sha256)) {
        if (isTRUE(verbose)) {
            message("Transcript conversion cache is missing trusted attestation metadata; regenerating ", basename(dst_path))
        }
        return(FALSE)
    }

    script_check <- .verify_transcripts_fallback_script_signature()
    if (!isTRUE(script_check$ok)) {
        if (isTRUE(verbose)) {
            message(
                "Transcript fallback script authenticity check failed (",
                script_check$reason,
                "); regenerating ",
                basename(dst_path)
            )
        }
        return(FALSE)
    }

    provenance_check <- .check_transcripts_integrity_provenance(report, src_path = src_path, dst_path = dst_path)
    if (!isTRUE(provenance_check$ok)) {
        if (isTRUE(verbose)) {
            message(
                "Transcript conversion provenance check failed (",
                provenance_check$reason,
                "); regenerating ",
                basename(dst_path)
            )
        }
        return(FALSE)
    }

    src_sha256 <- .compute_file_sha256(src_path)
    dst_sha256 <- .compute_file_sha256(dst_path)
    reuse_ok <- !is.na(src_sha256) && !is.na(dst_sha256) &&
        identical(unname(src_sha256), unname(report$source_sha256)) &&
        identical(unname(dst_sha256), unname(report$dest_sha256))
    if (!reuse_ok && isTRUE(verbose)) {
        message("Cached transcript conversion integrity mismatch; regenerating ", basename(dst_path))
    }
    reuse_ok
}

.probe_pyarrow_runtime <- function(python_bin, verbose = FALSE) {
    if (is.null(python_bin) || !nzchar(python_bin)) {
        stop("Python interpreter with pyarrow is required for transcript conversion.")
    }

    script_path <- tempfile("genescope_pyarrow_runtime_", fileext = ".py")
    writeLines(
        c(
            "import sys",
            "import pyarrow",
            "print(sys.version.split()[0])",
            "print(sys.executable)",
            "print(pyarrow.__version__)",
            "print(getattr(pyarrow, '__file__', ''))"
        ),
        con = script_path,
        useBytes = TRUE
    )
    on.exit(unlink(script_path), add = TRUE)

    out <- tryCatch(
        system2(python_bin, script_path, stdout = TRUE, stderr = TRUE),
        error = function(e) e
    )
    if (inherits(out, "error")) {
        stop("Failed to verify Python fallback runtime: ", conditionMessage(out))
    }
    status <- suppressWarnings(as.integer(attr(out, "status"))[1])
    if (!length(status) || !is.finite(status)) status <- 0L
    if (!identical(status, 0L)) {
        msg <- paste(out, collapse = "\n")
        if (!nzchar(msg)) msg <- paste0("python exited with status ", status)
        stop("Failed to verify Python fallback runtime: ", msg)
    }

    raw_runtime <- list(
        python_version = if (length(out) >= 1L) trimws(out[[1]]) else NA_character_,
        python_executable = if (length(out) >= 2L) trimws(out[[2]]) else NA_character_,
        pyarrow_version = if (length(out) >= 3L) trimws(out[[3]]) else NA_character_,
        pyarrow_module_path = if (length(out) >= 4L) trimws(out[[4]]) else NA_character_,
        python_executable_sha256 = NA_character_
    )
    raw_runtime$python_executable <- .normalize_transcripts_attestation_path(raw_runtime$python_executable)
    raw_runtime$python_executable_sha256 <- .compute_file_sha256(raw_runtime$python_executable)
    raw_runtime$pyarrow_module_path <- .normalize_transcripts_attestation_path(raw_runtime$pyarrow_module_path)
    openssl_runtime <- .probe_transcripts_runtime_openssl()
    raw_runtime$openssl_executable <- openssl_runtime$openssl_executable
    raw_runtime$openssl_version <- openssl_runtime$openssl_version
    raw_runtime$openssl_sha256 <- openssl_runtime$openssl_sha256

    runtime <- .attach_transcripts_runtime_attestation(
        runtime_attestation = raw_runtime,
        runtime_profile_path = .transcripts_runtime_profile_path()
    )

    if (!is.na(runtime$runtime_pyarrow_version) &&
        nzchar(runtime$runtime_pyarrow_version) &&
        utils::compareVersion(runtime$runtime_pyarrow_version, "10.0.0") < 0) {
        warning(
            "Python fallback detected pyarrow ", runtime$runtime_pyarrow_version,
            "; pyarrow >= 10.0.0 is recommended.",
            call. = FALSE,
            immediate. = TRUE
        )
    }
    if (isTRUE(verbose)) {
        message(
            "Using Python fallback runtime: python=",
            runtime$runtime_python_version,
            " pyarrow=",
            runtime$runtime_pyarrow_version,
            " trust=",
            runtime$runtime_trust_level
        )
    }
    runtime
}

.run_pyarrow_transcripts_conversion <- function(python_bin,
                                               src_path,
                                               dst_path,
                                               mode = c("parquet", "csv_gz"),
                                               verbose = FALSE,
                                               runtime_attestation = NULL) {
    mode <- match.arg(mode)
    if (is.null(python_bin) || !nzchar(python_bin)) {
        stop("Python interpreter with pyarrow is required for transcript conversion.")
    }
    if (is.null(src_path) || !nzchar(src_path) || !file.exists(src_path)) {
        stop("Source transcripts.parquet not found: ", src_path)
    }
    if (is.null(dst_path) || !nzchar(dst_path)) {
        stop("Destination path is required for transcript conversion.")
    }

    dst_dir <- dirname(dst_path)
    dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)

    tmp_path <- paste0(dst_path, ".tmp-", Sys.getpid(), "-", as.integer(Sys.time()))
    if (file.exists(tmp_path)) unlink(tmp_path)

    src_sha256 <- .compute_file_sha256(src_path)
    src_size <- tryCatch(file.info(src_path)$size[[1]], error = function(e) NA_real_)
    if (isTRUE(verbose)) {
        message(
            "Source fingerprint: sha256=",
            if (!is.na(src_sha256)) src_sha256 else "N/A",
            " size=",
            if (!is.na(src_size)) src_size else "N/A"
        )
    }

    script_check <- .verify_transcripts_fallback_script_signature()
    if (!isTRUE(script_check$ok)) {
        stop("Trusted transcript fallback script failed authenticity verification: ", script_check$reason)
    }
    script_path <- script_check$script_path

    if (is.null(runtime_attestation)) {
        runtime_attestation <- .probe_pyarrow_runtime(python_bin, verbose = verbose)
    }

    timeout_secs <- 300L
    warn_msgs <- character()
    timed_out <- FALSE
    sys2_args <- list(
        command = python_bin,
        args = c(script_path, src_path, tmp_path, mode),
        stdout = TRUE,
        stderr = TRUE
    )
    if ("timeout" %in% names(formals(system2))) {
        sys2_args$timeout <- timeout_secs
    }

    out <- tryCatch(
        withCallingHandlers(
            do.call(system2, sys2_args),
            warning = function(w) {
                warn_msgs <<- c(warn_msgs, conditionMessage(w))
                if (grepl("timed out", conditionMessage(w), ignore.case = TRUE)) {
                    timed_out <<- TRUE
                }
                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) e
    )
    if (inherits(out, "error")) {
        unlink(tmp_path)
        stop("Automatic transcript conversion failed: ", conditionMessage(out))
    }
    status <- suppressWarnings(as.integer(attr(out, "status"))[1])
    if (!length(status) || !is.finite(status)) status <- 0L
    if (isTRUE(timed_out) || identical(status, 124L)) {
        unlink(tmp_path)
        stop("Automatic transcript conversion timed out after ", timeout_secs, " seconds.")
    }
    if (!identical(status, 0L)) {
        unlink(tmp_path)
        msg <- paste(c(out, warn_msgs), collapse = "\n")
        if (!nzchar(msg)) msg <- paste0("python exited with status ", status)
        stop("Automatic transcript conversion failed: ", msg)
    }
    if (!file.exists(tmp_path)) {
        stop("Automatic transcript conversion did not create output: ", dst_path)
    }

    dst_size <- tryCatch(file.info(tmp_path)$size[[1]], error = function(e) NA_real_)
    if (is.na(dst_size) || dst_size <= 0) {
        unlink(tmp_path)
        stop("Automatic transcript conversion produced empty output: ", dst_path)
    }
    dst_sha256 <- .compute_file_sha256(tmp_path)
    if (isTRUE(verbose)) {
        message(
            "Output fingerprint: sha256=",
            if (!is.na(dst_sha256)) dst_sha256 else "N/A",
            " size=",
            dst_size
        )
    }

    if (file.exists(dst_path)) unlink(dst_path)
    ok_rename <- file.rename(tmp_path, dst_path)
    if (!isTRUE(ok_rename) || !file.exists(dst_path)) {
        unlink(tmp_path)
        stop("Automatic transcript conversion could not finalize output: ", dst_path)
    }

    report <- .build_transcripts_integrity_report(
        src_path = src_path,
        dst_path = dst_path,
        mode = mode,
        python_bin = python_bin,
        src_sha256 = src_sha256,
        src_size = src_size,
        dst_sha256 = dst_sha256,
        dst_size = dst_size,
        timeout_secs = timeout_secs,
        runtime_attestation = runtime_attestation
    )
    .write_transcripts_integrity_sidecar(dst_path, report)
    prov_check <- .check_transcripts_integrity_provenance(report, src_path = src_path, dst_path = dst_path)
    if (!isTRUE(prov_check$ok)) {
        unlink(c(dst_path, .transcripts_integrity_sidecar_path(dst_path)))
        stop(
            "Automatic transcript conversion authenticity verification failed: ",
            prov_check$reason
        )
    }

    if (isTRUE(verbose)) {
        message("Created transcript ", mode, " fallback: ", dst_path)
    }
    out_path <- dst_path
    attr(out_path, "genescope_conversion_report") <- report
    out_path
}

.ensure_transcripts_arrow_compat_parquet <- function(tf_parq, input_dir = NULL, verbose = FALSE) {
    if (is.null(tf_parq) || !nzchar(tf_parq) || !file.exists(tf_parq)) {
        stop("Parquet file not found: ", tf_parq)
    }
    cache_dir <- .resolve_transcripts_ingest_cache_dir(if (!is.null(input_dir) && nzchar(input_dir)) input_dir else dirname(tf_parq))
    if (is.null(cache_dir) || !nzchar(cache_dir)) {
        stop("Cannot determine ingest cache directory for Arrow-compatible transcripts.")
    }
    dst_path <- file.path(cache_dir, "transcripts.arrow_compat.parquet")
    if (.should_reuse_transcripts_conversion(tf_parq, dst_path, verbose = verbose)) {
        return(dst_path)
    }

    python_bin <- .resolve_pyarrow_python(verbose = verbose)
    if (is.null(python_bin)) {
        stop("No Python interpreter with pyarrow is available for Arrow-compatible transcript conversion.")
    }
    runtime_attestation <- .probe_pyarrow_runtime(python_bin, verbose = verbose)
    .run_pyarrow_transcripts_conversion(
        python_bin = python_bin,
        src_path = tf_parq,
        dst_path = dst_path,
        mode = "parquet",
        verbose = verbose,
        runtime_attestation = runtime_attestation
    )
}

.ensure_transcripts_csv_fallback <- function(tf_parq, input_dir = NULL, verbose = FALSE) {
    if (is.null(input_dir) || !length(input_dir) || is.na(input_dir) || !nzchar(input_dir)) {
        input_dir <- dirname(tf_parq)
    }
    dst_path <- file.path(input_dir, "transcripts.csv.gz")
    if (.should_reuse_transcripts_conversion(tf_parq, dst_path, verbose = verbose)) {
        return(dst_path)
    }

    python_bin <- .resolve_pyarrow_python(verbose = verbose)
    if (is.null(python_bin)) {
        stop("No Python interpreter with pyarrow is available for transcripts.csv.gz fallback generation.")
    }
    runtime_attestation <- .probe_pyarrow_runtime(python_bin, verbose = verbose)
    .run_pyarrow_transcripts_conversion(
        python_bin = python_bin,
        src_path = tf_parq,
        dst_path = dst_path,
        mode = "csv_gz",
        verbose = verbose,
        runtime_attestation = runtime_attestation
    )
}

.probe_transcripts_dataset_query <- function(ds) {
    cols_all <- tryCatch(names(ds), error = function(e) character())
    probe_cols <- intersect(c("x_location", "y_location", "feature_name", "codeword_index"), cols_all)
    probe_cols <- probe_cols[seq_len(min(length(probe_cols), 2L))]
    if (!length(probe_cols)) {
        return(list(ok = TRUE, cols = character(), error = NULL))
    }
    probe <- tryCatch(
        {
            if (inherits(ds, "Dataset") && "NewScan" %in% names(ds)) {
                scan_builder <- ds$NewScan()
                scan_builder$Project(probe_cols)
                scan_builder$UseThreads(TRUE)
                scan_builder$BatchSize(1024L)
                reader <- scan_builder$Finish()$ToRecordBatchReader()
            } else {
                reader <- arrow::as_record_batch_reader(
                    ds |>
                        dplyr::select(dplyr::all_of(probe_cols))
                )
            }
            reader$read_next_batch()
        },
        error = function(e) e
    )
    if (inherits(probe, "error")) {
        return(list(ok = FALSE, cols = probe_cols, error = probe))
    }
    list(ok = TRUE, cols = probe_cols, error = NULL)
}

.read_parquet_cols_as_dt <- function(path, cols, use_threads = TRUE) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("Package 'arrow' is required to read Parquet files.")
    }
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required to process transcripts.")
    }
    if (is.null(path) || !nzchar(path) || !file.exists(path)) {
        stop("Parquet file not found: ", path)
    }
    cols <- unique(as.character(cols))
    cols <- cols[!is.na(cols) & nzchar(cols)]
    if (!length(cols)) stop("No columns requested from Parquet.")

    rp <- getExportedValue("arrow", "read_parquet")
    fml <- names(formals(rp))
    args <- list(path)
    if ("as_data_frame" %in% fml) args$as_data_frame <- TRUE
    if ("col_select" %in% fml) {
        args$col_select <- cols
    } else if ("columns" %in% fml) {
        args$columns <- cols
    }
    if ("use_threads" %in% fml) args$use_threads <- isTRUE(use_threads)

    out <- tryCatch(do.call(rp, args), error = function(e) e)
    if (inherits(out, "error")) {
        args2 <- args
        if ("col_select" %in% names(args2)) args2$col_select <- NULL
        if ("columns" %in% names(args2)) args2$columns <- NULL
        out2 <- tryCatch(do.call(rp, args2), error = function(e) e)
        if (inherits(out2, "error")) {
            stop(conditionMessage(out))
        }
        out <- out2
    }

    dt <- data.table::as.data.table(out)
    missing <- setdiff(cols, names(dt))
    if (length(missing)) {
        stop("Parquet is missing required column(s): ", paste(missing, collapse = ", "))
    }
    dt[, cols, with = FALSE]
}

.read_transcripts_csv_cols_as_dt <- function(path, cols, verbose = FALSE) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required to read transcripts CSV.")
    }
    if (is.null(path) || !nzchar(path) || !file.exists(path)) {
        stop("Transcripts CSV not found: ", path)
    }
    cols <- unique(as.character(cols))
    cols <- cols[!is.na(cols) & nzchar(cols)]
    if (!length(cols)) stop("No columns requested from transcripts CSV.")

    is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)
    dt <- tryCatch(
        data.table::fread(path, select = cols, showProgress = isTRUE(verbose)),
        error = function(e) NULL
    )
    if (is.null(dt) && isTRUE(is_gz)) {
        has_zcat <- nzchar(Sys.which("zcat"))
        has_gzcat <- nzchar(Sys.which("gzcat"))
        has_pigz <- nzchar(Sys.which("pigz"))
        cmd <- NULL
        if (has_gzcat) {
            cmd <- paste(Sys.which("gzcat"), shQuote(path))
        } else if (has_zcat) {
            cmd <- paste(Sys.which("zcat"), shQuote(path))
        } else if (has_pigz) {
            cmd <- paste(Sys.which("pigz"), "-dc", shQuote(path))
        }
        if (!is.null(cmd) && nzchar(cmd)) {
            dt <- tryCatch(
                data.table::fread(cmd = cmd, select = cols, showProgress = isTRUE(verbose)),
                error = function(e) NULL
            )
        }
    }
    if (is.null(dt)) {
        stop("Failed to read transcripts CSV: ", path)
    }
    data.table::setDT(dt)
    missing <- setdiff(cols, names(dt))
    if (length(missing)) {
        stop("Transcripts CSV is missing required column(s): ", paste(missing, collapse = ", "))
    }
    dt
}

.read_transcripts_fallback_dt <- function(paths,
                                         input_dir,
                                         cols,
                                         verbose = FALSE,
                                         use_threads = TRUE,
                                         log_parent = "createSCOPE",
                                         log_step = "S08") {
    tf_parq <- NULL
    tf_csv <- NULL
    log_parent <- if (!is.null(log_parent) && nzchar(log_parent)) as.character(log_parent)[1] else "createSCOPE"
    log_step <- if (!is.null(log_step) && nzchar(log_step)) as.character(log_step)[1] else "S08"
    log_fallback <- function(msg) .log_info(log_parent, log_step, msg, verbose)
    format_size_gb <- function(path) {
        sz <- tryCatch(file.info(path)$size[[1]], error = function(e) NA_real_)
        if (!is.finite(sz) || sz <= 0) return(NA_character_)
        paste0(format(round(sz / 1024^3, 2), trim = TRUE, scientific = FALSE), "GB")
    }
    if (!is.null(paths) && is.list(paths)) {
        if (!is.null(paths$tf_parq)) tf_parq <- as.character(paths$tf_parq)[1]
        if (!is.null(paths$tf_csv)) tf_csv <- as.character(paths$tf_csv)[1]
    }
    if (is.null(tf_csv) || is.na(tf_csv) || !nzchar(tf_csv)) {
        tf_csv <- .resolve_transcripts_csv_path(input_dir)
    }

    err_parq <- NULL
    if (is.character(tf_parq) && !is.na(tf_parq) && nzchar(tf_parq) && file.exists(tf_parq)) {
        dt_parq <- tryCatch(.read_parquet_cols_as_dt(tf_parq, cols, use_threads = use_threads), error = function(e) e)
        if (!inherits(dt_parq, "error")) {
            log_fallback(paste0("fallback transcripts source=parquet path=", basename(tf_parq), " rows=", format(nrow(dt_parq), big.mark = ",")))
            return(list(dt = dt_parq, source = "parquet", path = tf_parq))
        }
        err_parq <- dt_parq
        log_fallback(paste0("fallback parquet read failed: ", conditionMessage(err_parq)))
    }

    err_csv_build <- NULL
    if ((is.null(tf_csv) || is.na(tf_csv) || !nzchar(tf_csv) || !file.exists(tf_csv)) &&
        is.character(tf_parq) && !is.na(tf_parq) && nzchar(tf_parq) && file.exists(tf_parq)) {
        log_fallback("no transcripts.csv(.gz) fallback found; attempting pyarrow conversion from transcripts.parquet")
        tf_csv_built <- tryCatch(
            .ensure_transcripts_csv_fallback(tf_parq, input_dir = input_dir, verbose = verbose),
            error = function(e) e
        )
        if (!inherits(tf_csv_built, "error")) {
            tf_csv <- tf_csv_built
            log_fallback(paste0("created transcripts CSV fallback: ", basename(tf_csv)))
        } else {
            err_csv_build <- tf_csv_built
            log_fallback(paste0("pyarrow CSV fallback conversion failed: ", conditionMessage(err_csv_build)))
        }
    }

    err_csv <- NULL
    if (is.character(tf_csv) && !is.na(tf_csv) && nzchar(tf_csv) && file.exists(tf_csv)) {
        tf_csv_size <- format_size_gb(tf_csv)
        log_fallback(
            paste0(
                "using transcripts CSV fallback: ", basename(tf_csv),
                if (!is.na(tf_csv_size)) paste0(" (", tf_csv_size, ")") else "",
                "; reading selected columns into memory"
            )
        )
        dt_csv <- tryCatch(.read_transcripts_csv_cols_as_dt(tf_csv, cols, verbose = verbose), error = function(e) e)
        if (!inherits(dt_csv, "error")) {
            log_fallback(paste0("fallback transcripts source=csv path=", basename(tf_csv), " rows=", format(nrow(dt_csv), big.mark = ",")))
            return(list(dt = dt_csv, source = "csv", path = tf_csv))
        }
        err_csv <- dt_csv
        log_fallback(paste0("fallback CSV read failed: ", conditionMessage(err_csv)))
    }

    msg <- character()
    if (!is.null(err_parq)) msg <- c(msg, paste0("parquet_error=", conditionMessage(err_parq)))
    if (!is.null(err_csv_build)) msg <- c(msg, paste0("csv_build_error=", conditionMessage(err_csv_build)))
    if (!is.null(err_csv)) msg <- c(msg, paste0("csv_error=", conditionMessage(err_csv)))
    if (!length(msg)) msg <- "no_transcripts_source_available"
    stop("Failed to read transcripts from Parquet/CSV (", paste(msg, collapse = "; "), ").")
}

#' Apply Dataset Filters
#' @description
#' Internal helper for `.apply_dataset_filters`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.apply_dataset_filters <- function(state) {
    p <- state$params
    ds <- state$datasets$ds
    parent <- "createSCOPE"
    verbose <- p$verbose
    ctx <- state$thread_context
    paths <- state$paths
    if (is.null(state$params$transcripts_spec) || !is.list(state$params$transcripts_spec)) {
        state$params$transcripts_spec <- .resolve_transcripts_feature_spec(
            ds = ds,
            input_dir = p$input_dir,
            verbose = verbose
        )
        p <- state$params
    }
    spec <- p$transcripts_spec

    ds_cols <- tryCatch(names(ds), error = function(e) character())
    if (!length(ds_cols)) ds_cols <- character()
    if (!is.null(p$max_dist_mol_nuc) && !"nucleus_distance" %in% ds_cols) {
        .log_info(parent, "S05", "missing nucleus_distance; skipping max_dist_mol_nuc filter", verbose)
        state$params$max_dist_mol_nuc <- NULL
        p$max_dist_mol_nuc <- NULL
    }
    if (!is.null(p$filterqv) && !"qv" %in% ds_cols) {
        .log_info(parent, "S07", "missing qv; skipping filterqv filter", verbose)
        state$params$filterqv <- NULL
        p$filterqv <- NULL
    }
    tf_path <- if (!is.null(paths$tf_parq)) paths$tf_parq else NA_character_
    size_gb <- tryCatch({
        info <- file.info(tf_path)
        if (is.finite(info$size)) info$size / 1024^3 else NA_real_
    }, error = function(e) NA_real_)
    mem_gb <- if (!is.null(ctx) && !is.null(ctx$sys_mem_gb) && is.finite(ctx$sys_mem_gb)) as.numeric(ctx$sys_mem_gb) else NA_real_
    mem_ratio <- if (is.finite(size_gb) && is.finite(mem_gb) && mem_gb > 0) size_gb / mem_gb else NA_real_
    cpu_target <- if (!is.null(ctx) && !is.null(ctx$ncores_safe)) as.integer(ctx$ncores_safe) else NA_integer_
    cpu_target <- ifelse(is.finite(cpu_target) && cpu_target > 0L, cpu_target, as.integer(getOption("arrow.num_threads", 1L)))
    cpu_target <- max(1L, cpu_target)

    ncores_arg <- suppressWarnings(as.integer(p$ncores)[1])
    user_requested <- length(ncores_arg) && is.finite(ncores_arg) && ncores_arg > 0L
    auto_throttle <- !user_requested && isTRUE(getOption("geneSCOPE.arrow_auto_throttle", TRUE))

    cpu_use <- cpu_target
    if (auto_throttle) {
        cap_factor <- NA_real_
        if (is.finite(mem_ratio)) {
            if (mem_ratio >= 0.50) {
                cap_factor <- 0.25
            } else if (mem_ratio >= 0.25) {
                cap_factor <- 0.50
            } else if (mem_ratio >= 0.10) {
                cap_factor <- 0.75
            }
        } else if (is.finite(size_gb)) {
            if (size_gb >= 40) {
                cap_factor <- 0.25
            } else if (size_gb >= 20) {
                cap_factor <- 0.50
            } else if (size_gb >= 5) {
                cap_factor <- 0.75
            }
        }
        if (!is.na(cap_factor)) {
            cpu_use <- suppressWarnings(as.integer(ceiling(cpu_target * cap_factor)))
            if (!is.finite(cpu_use) || cpu_use < 1L) cpu_use <- 1L
            if (cpu_target > 1L) cpu_use <- max(2L, cpu_use)
            cpu_use <- min(cpu_use, cpu_target)
        }
    }

    if (requireNamespace("arrow", quietly = TRUE) && cpu_use < cpu_target) {
        old_opt <- getOption("arrow.num_threads")
        old_cpu <- if (exists("cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
            tryCatch(arrow::cpu_count(), error = function(e) NA_integer_)
        } else {
            NA_integer_
        }
        old_io <- if (exists("io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
            tryCatch(arrow::io_thread_count(), error = function(e) NA_integer_)
        } else {
            NA_integer_
        }
        on.exit({
            try(options(arrow.num_threads = old_opt), silent = TRUE)

            restore_cpu <- suppressWarnings(as.integer(old_cpu))
            if (!is.finite(restore_cpu) || restore_cpu < 1L) restore_cpu <- cpu_target
            restore_cpu <- max(1L, restore_cpu)

            restore_io <- suppressWarnings(as.integer(old_io))
            if (!is.finite(restore_io) || restore_io < 1L) restore_io <- restore_cpu
            restore_io <- max(2L, restore_io)

            if (exists("set_cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
                try(arrow::set_cpu_count(restore_cpu), silent = TRUE)
            }
            if (exists("set_io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
                try(arrow::set_io_thread_count(restore_io), silent = TRUE)
            }
        }, add = TRUE)

        io_use <- max(2L, cpu_use)
        .log_info(
            parent,
            "S05",
            paste0(
                "throttling Arrow threads for dataset filters: cpu=", cpu_use,
                " io=", io_use,
                " (restore cpu=", cpu_target,
                if (is.finite(size_gb)) paste0(", transcripts_size_gb=", round(size_gb, 1)) else "",
                if (is.finite(mem_gb)) paste0(", sys_mem_gb=", round(mem_gb, 1)) else "",
                ")"
            ),
            verbose
        )
        try(options(arrow.num_threads = cpu_use), silent = TRUE)
        if (exists("set_cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
            try(arrow::set_cpu_count(cpu_use), silent = TRUE)
        }
        if (exists("set_io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
            try(arrow::set_io_thread_count(io_use), silent = TRUE)
        }
    }

    log_dataset_row_counts <- isTRUE(getOption("geneSCOPE.log_dataset_row_counts", FALSE))
    dataset_count <- function(ds_in) {
        if (!.log_enabled(verbose) || !isTRUE(log_dataset_row_counts)) return(NA_integer_)
        tryCatch({
            out <- ds_in |> summarise(n = n()) |> collect()
            as.integer(out$n[[1]])
        }, error = function(e) NA_integer_)
    }

    qc_genes <- character()

    step05 <- .log_step(parent, "S05", "filter by nucleus distance", verbose)
    step05$enter(paste0("max=", ifelse(is.null(p$max_dist_mol_nuc), "NA", p$max_dist_mol_nuc)))
    if (!is.null(p$filter_genes)) {
        .log_info(parent, "S05", paste0("filtering genes: ", length(p$filter_genes), " genes"), verbose)
        if ("feature_name" %in% names(ds)) {
            ds <- ds |> filter(feature_name %in% p$filter_genes)
        } else {
            .log_info(parent, "S05", "filter_genes deferred (missing feature_name column)", verbose)
        }
    }
    if (!is.null(p$max_dist_mol_nuc)) {
        ds <- ds |> filter(nucleus_distance <= p$max_dist_mol_nuc)
    }
    if (is.null(p$filter_genes) && is.null(p$max_dist_mol_nuc)) {
        step05$done("skipped")
    } else {
        kept <- dataset_count(ds)
        step05$done(paste0("retained=", ifelse(is.na(kept), "NA", kept)))
    }

    if (isTRUE(p$exclude_qc_unsuitable)) {
        step06a <- .log_step(parent, "S06a", "exclude QC-flagged genes", verbose)
        step06a$enter(paste0(
            "qc_out_dir=", ifelse(is.null(p$qc_out_dir), "NULL", p$qc_out_dir),
            " qc_flag=", ifelse(is.null(p$qc_flag), "NULL", paste(p$qc_flag, collapse = ","))
        ))
        qc_genes <- .gridselect_qc_genes_to_exclude(p$qc_out_dir, p$qc_flag)
        if (length(qc_genes)) {
            if ("feature_name" %in% names(ds)) {
                ds <- ds |> filter(!feature_name %in% qc_genes)
            } else {
                .log_info(parent, "S06a", "QC gene exclusion deferred (missing feature_name column)", verbose)
            }
        }
        kept <- dataset_count(ds)
        step06a$done(paste0(
            "excluded_genes=", length(qc_genes),
            " retained=", ifelse(is.na(kept), "NA", kept)
        ))
    }

    step06 <- .log_step(parent, "S06", "exclude molecule prefixes", verbose)
    if (isTRUE(p$filtermolecule)) {
        step06$enter(paste0("prefixes=", paste(p$exclude_prefix, collapse = ",")))
        step06$done("deferred_to_batch_scan")
    } else {
        step06$enter("skipped")
        step06$done("skipped")
    }

    step07 <- .log_step(parent, "S07", "filter by quality score", verbose)
    if (!is.null(p$filterqv)) {
        step07$enter(paste0("min=", p$filterqv))
        ds <- ds |> filter(qv >= p$filterqv)
        kept <- dataset_count(ds)
        step07$done(paste0("retained=", ifelse(is.na(kept), "NA", kept)))
    } else {
        step07$enter("skipped")
        step07$done("skipped")
    }

    step08 <- .log_step(parent, "S08", "compute bounds and transforms", verbose)
    step08$enter(paste0("flip_y=", p$flip_y))
    bounds_global <- tryCatch({
        cols_needed <- c("x_location", "y_location")
        if (!is.null(p$filter_genes) || isTRUE(p$exclude_qc_unsuitable) || isTRUE(p$filtermolecule)) {
            cols_needed <- c(cols_needed, .transcripts_feature_project_cols(spec))
        }
        if (!is.null(p$max_dist_mol_nuc)) cols_needed <- c(cols_needed, "nucleus_distance")
        if (!is.null(p$filterqv)) cols_needed <- c(cols_needed, "qv")
        cols_needed <- unique(cols_needed)

        bounds_scan <- function(ds_in, batch_size, report_every = 5000000L) {
            batch_size <- suppressWarnings(as.integer(batch_size))
            batch_size <- ifelse(is.finite(batch_size) && batch_size > 0L, batch_size, 500000L)
            report_every <- suppressWarnings(as.integer(report_every))
            report_every <- ifelse(is.finite(report_every) && report_every > 0L, report_every, 5000000L)

            reader <- tryCatch({
                if (inherits(ds_in, "Dataset")) {
                    scan_builder <- ds_in$NewScan()
                    scan_builder$Project(cols_needed)
                    scan_builder$UseThreads(TRUE)
                    scan_builder$BatchSize(batch_size)
                    scan_builder$Finish()$ToRecordBatchReader()
                } else {
                    arrow::as_record_batch_reader(ds_in |>
                        dplyr::select(dplyr::all_of(cols_needed)))
                }
            }, error = function(e) NULL)
            if (is.null(reader)) stop("Failed to create Arrow record batch reader for bounds scan.")

            xmin <- Inf
            xmax <- -Inf
            ymin <- Inf
            ymax <- -Inf
            scanned <- 0L
            next_report <- report_every

            repeat {
                batch <- tryCatch(reader$read_next_batch(), error = function(e) {
                    stop("Bounds scan failed: ", e$message)
                })
                if (is.null(batch)) break
                tbl <- tryCatch(as.data.frame(batch), error = function(e) NULL)
                if (is.null(tbl) || !nrow(tbl)) next
                scanned <- scanned + nrow(tbl)
                mask <- .ingest_molecule_filter_mask(tbl, p, qc_genes)
                if (length(mask) && any(mask)) {
                    x <- tbl$x_location[mask]
                    y <- tbl$y_location[mask]

                    xmin_b <- suppressWarnings(min(x, na.rm = TRUE))
                    xmax_b <- suppressWarnings(max(x, na.rm = TRUE))
                    ymin_b <- suppressWarnings(min(y, na.rm = TRUE))
                    ymax_b <- suppressWarnings(max(y, na.rm = TRUE))
                    if (is.finite(xmin_b)) xmin <- min(xmin, xmin_b)
                    if (is.finite(xmax_b)) xmax <- max(xmax, xmax_b)
                    if (is.finite(ymin_b)) ymin <- min(ymin, ymin_b)
                    if (is.finite(ymax_b)) ymax <- max(ymax, ymax_b)
                }
                if (.log_enabled(verbose) && scanned >= next_report) {
                    .log_info(parent, "S08", paste0("bounds scan progress rows=", format(scanned, big.mark = ",")), verbose)
                    next_report <- scanned + report_every
                }
            }

            if (!is.finite(xmin) || !is.finite(xmax) || !is.finite(ymin) || !is.finite(ymax)) {
                stop("No molecules available after filters; cannot compute bounds.")
            }
            data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
        }

        if (requireNamespace("arrow", quietly = TRUE)) {
            bounds_cpu_cap <- suppressWarnings(as.integer(getOption("geneSCOPE.arrow_bounds_cpu_cap", NA_integer_)))
            bounds_io_cap <- suppressWarnings(as.integer(getOption("geneSCOPE.arrow_bounds_io_cap", bounds_cpu_cap)))
            bounds_cpu_cap <- ifelse(is.finite(bounds_cpu_cap) && bounds_cpu_cap > 0L, bounds_cpu_cap, NA_integer_)
            bounds_io_cap <- ifelse(is.finite(bounds_io_cap) && bounds_io_cap > 0L, bounds_io_cap, NA_integer_)

            old_opt <- getOption("arrow.num_threads")
            old_cpu <- if (exists("cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
                tryCatch(arrow::cpu_count(), error = function(e) NA_integer_)
            } else {
                NA_integer_
            }
            old_io <- if (exists("io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
                tryCatch(arrow::io_thread_count(), error = function(e) NA_integer_)
            } else {
                NA_integer_
            }
	        on.exit({
	            try(options(arrow.num_threads = old_opt), silent = TRUE)

                restore_cpu <- suppressWarnings(as.integer(old_cpu))
                if (!is.finite(restore_cpu) || restore_cpu < 1L) restore_cpu <- cpu_target
                restore_cpu <- max(1L, restore_cpu)

                restore_io <- suppressWarnings(as.integer(old_io))
                if (!is.finite(restore_io) || restore_io < 1L) restore_io <- restore_cpu
                restore_io <- max(2L, restore_io)

	            if (exists("set_cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
	                try(arrow::set_cpu_count(restore_cpu), silent = TRUE)
	            }
	            if (exists("set_io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
	                try(arrow::set_io_thread_count(restore_io), silent = TRUE)
	            }
	        }, add = TRUE)

            cpu_now <- suppressWarnings(as.integer(old_cpu))
            if (!is.finite(cpu_now) || cpu_now < 1L) cpu_now <- cpu_target
	            io_now <- suppressWarnings(as.integer(old_io))
	            if (!is.finite(io_now) || io_now < 1L) io_now <- max(2L, cpu_now)
	            cpu_use_bounds <- cpu_now
	            io_use_bounds <- io_now
	            if (!is.na(bounds_cpu_cap)) cpu_use_bounds <- min(cpu_use_bounds, bounds_cpu_cap)
	            if (!is.na(bounds_io_cap)) io_use_bounds <- min(io_use_bounds, bounds_io_cap)
	            cpu_use_bounds <- max(1L, as.integer(cpu_use_bounds))
	            io_use_bounds <- max(2L, as.integer(io_use_bounds))

	            applied <- FALSE
            if (cpu_use_bounds < cpu_now) {
                applied <- TRUE
                try(options(arrow.num_threads = cpu_use_bounds), silent = TRUE)
                if (exists("set_cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
                    try(arrow::set_cpu_count(cpu_use_bounds), silent = TRUE)
                }
            }
            if (io_use_bounds < io_now) {
                applied <- TRUE
                if (exists("set_io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
                    try(arrow::set_io_thread_count(io_use_bounds), silent = TRUE)
                }
            }
            if (applied) {
                .log_info(
                    parent,
                    "S08",
                    paste0(
                        "Arrow bounds scan threads: cpu=", cpu_use_bounds,
                        " io=", io_use_bounds,
                        " (restore cpu=", ifelse(is.finite(old_cpu), old_cpu, "NA"),
                        " io=", ifelse(is.finite(old_io), old_io, "NA"),
                        ")"
                    ),
                    verbose
                )
            }

            set_arrow_threads <- function(cpu_threads, io_threads) {
                cpu_threads <- suppressWarnings(as.integer(cpu_threads))
                if (!is.finite(cpu_threads) || cpu_threads < 1L) cpu_threads <- 1L
                io_threads <- suppressWarnings(as.integer(io_threads))
                if (!is.finite(io_threads) || io_threads < 1L) io_threads <- cpu_threads
                io_threads <- max(2L, io_threads)
                try(options(arrow.num_threads = cpu_threads), silent = TRUE)
                if (exists("set_cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
                    try(arrow::set_cpu_count(cpu_threads), silent = TRUE)
                }
                if (exists("set_io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
                    try(arrow::set_io_thread_count(io_threads), silent = TRUE)
                }
                invisible(list(cpu = cpu_threads, io = io_threads))
            }

            retry_default <- c(32L, 16L, 8L, 4L, 2L, 1L)
            retry_opt <- getOption("geneSCOPE.arrow_bounds_retry_threads", retry_default)
            retry_vec <- suppressWarnings(as.integer(retry_opt))
            retry_vec <- retry_vec[is.finite(retry_vec) & retry_vec > 0L]
            cpu_trials <- unique(c(cpu_use_bounds, retry_vec))
            cpu_trials <- unique(pmax(1L, pmin(cpu_use_bounds, cpu_trials)))
            cpu_trials <- sort(cpu_trials, decreasing = TRUE)
            if (!length(cpu_trials)) cpu_trials <- cpu_use_bounds

            run_bounds_scan <- function(ds_in, dataset_label) {
                last_err <- NULL
                for (i in seq_along(cpu_trials)) {
                    cpu_i <- cpu_trials[[i]]
                    io_i <- max(2L, cpu_i)
                    if (i > 1L) {
                        .log_info(parent, "S08", paste0("bounds scan retry: dataset=", dataset_label, " cpu=", cpu_i, " io=", io_i), verbose)
                    }
                    set_arrow_threads(cpu_i, io_i)
                    res_i <- tryCatch(
                        bounds_scan(ds_in,
                            batch_size = p$chunk_pts,
                            report_every = getOption("geneSCOPE.bounds_progress_rows", 5000000L)
                        ),
                        error = function(e) e
                    )
                    if (!inherits(res_i, "error")) {
                        if (i > 1L) {
                            .log_info(parent, "S08", paste0("bounds scan recovered: dataset=", dataset_label, " cpu=", cpu_i, " io=", io_i), verbose)
                        }
                        return(list(ok = TRUE, result = res_i, cpu = cpu_i, io = io_i))
                    }
                    if (i == 1L && length(cpu_trials) > 1L) {
                        .log_info(
                            parent,
                            "S08",
                            paste0(
                                "bounds scan failed (dataset=", dataset_label,
                                " cpu=", cpu_i, " io=", io_i,
                                "); retrying: ", conditionMessage(res_i)
                            ),
                            verbose
                        )
                    }
                    last_err <- res_i
                }
                list(ok = FALSE, error = last_err)
            }

	            ds_queryable <- !isFALSE(state$datasets$ds_queryable)
	            probe_err <- if (!is.null(state$datasets$ds_query_probe_error) && nzchar(state$datasets$ds_query_probe_error)) {
	                state$datasets$ds_query_probe_error
	            } else {
	                "Arrow dataset query probe failed"
	            }
	            attempt <- if (isTRUE(ds_queryable)) {
	                run_bounds_scan(ds, "filtered")
	            } else {
	                list(ok = FALSE, error = simpleError(paste0("Arrow dataset query unavailable: ", probe_err)))
	            }
	            if (isTRUE(attempt$ok)) {
	                attempt$result
	            } else {
	                .log_info(
	                    parent,
	                    "S08",
	                    paste0("bounds scan failed on filtered dataset; retrying on raw dataset: ", conditionMessage(attempt$error)),
	                    verbose
	                )
	                raw_ds <- state$datasets$ds_raw
	                if (is.null(raw_ds)) raw_ds <- ds
	                attempt_raw <- if (isTRUE(ds_queryable)) {
	                    run_bounds_scan(raw_ds, "raw")
	                } else {
	                    list(ok = FALSE, error = attempt$error)
	                }
	                if (isTRUE(attempt_raw$ok)) {
	                    ds <- raw_ds
	                    attempt_raw$result
	                } else {
	                    can_try_compat <- is.character(paths$tf_parq) &&
	                        !is.na(paths$tf_parq) &&
	                        nzchar(paths$tf_parq) &&
	                        !grepl("transcripts\\.arrow_compat\\.parquet$", basename(paths$tf_parq), ignore.case = TRUE)
		                    if (isTRUE(can_try_compat)) {
		                        .log_info(parent, "S08", "attempting Arrow-compatible transcript rewrite via pyarrow before direct fallback read", verbose)
		                        compat_path <- tryCatch(
		                            .ensure_transcripts_arrow_compat_parquet(paths$tf_parq, input_dir = p$input_dir, verbose = verbose),
		                            error = function(e) e
		                        )
	                        if (!inherits(compat_path, "error")) {
	                            compat_ds <- tryCatch(arrow::open_dataset(compat_path), error = function(e) e)
	                            if (!inherits(compat_ds, "error")) {
	                                .log_info(parent, "S08", paste0("retrying bounds scan on Arrow-compatible transcripts: ", basename(compat_path)), verbose)
	                                attempt_compat <- run_bounds_scan(compat_ds, "compat")
	                                if (isTRUE(attempt_compat$ok)) {
	                                    state$paths$tf_parq_original <- paths$tf_parq
	                                    state$paths$tf_parq <- compat_path
	                                    state$datasets$ds <- compat_ds
	                                    state$datasets$ds_raw <- compat_ds
	                                    state$datasets$ds_queryable <- TRUE
	                                    state$datasets$ds_query_probe_error <- NULL
	                                    paths <- state$paths
	                                    ds <- compat_ds
	                                    attempt_raw <- attempt_compat
	                                } else {
	                                    .log_info(parent, "S08", paste0("Arrow-compatible bounds scan still failed: ", conditionMessage(attempt_compat$error)), verbose)
	                                }
	                            } else {
	                                .log_info(parent, "S08", paste0("Arrow-compatible transcripts open failed: ", conditionMessage(compat_ds)), verbose)
	                            }
	                        } else {
	                            .log_info(parent, "S08", paste0("Arrow-compatible transcript rewrite failed: ", conditionMessage(compat_path)), verbose)
	                        }
	                    }
	                }
	                if (isTRUE(attempt_raw$ok)) {
	                    ds <- if (!is.null(state$datasets$ds_raw)) state$datasets$ds_raw else raw_ds
	                    attempt_raw$result
	                } else {
	                    .log_info(
	                        parent,
	                        "S08",
	                        paste0(
	                            "bounds scan failed on raw Arrow dataset; falling back to direct transcripts read: ",
	                            conditionMessage(attempt_raw$error)
	                        ),
	                        verbose
	                    )
	                    ds <- raw_ds
	                    try(set_arrow_threads(cpu_target, max(2L, cpu_target)), silent = TRUE)

		                    cols_cache <- unique(c(cols_needed, .transcripts_feature_project_cols(spec)))
		                    .log_info(parent, "S08", "starting fallback transcripts load for bounds scan", verbose)
		                    fb <- tryCatch(
		                        .read_transcripts_fallback_dt(
		                            paths,
		                            p$input_dir,
		                            cols_cache,
		                            verbose = verbose,
		                            use_threads = TRUE,
		                            log_parent = parent,
		                            log_step = "S08"
		                        ),
		                        error = function(e) e
		                    )
		                    if (inherits(fb, "error")) {
		                        .log_info(parent, "S08", paste0("fallback transcripts read failed; retrying single-thread: ", conditionMessage(fb)), verbose)
		                        fb <- .read_transcripts_fallback_dt(
		                            paths,
		                            p$input_dir,
		                            cols_cache,
		                            verbose = verbose,
		                            use_threads = FALSE,
		                            log_parent = parent,
		                            log_step = "S08"
		                        )
		                    }

	                    dt_fb <- fb$dt
	                    mask_fb <- .ingest_molecule_filter_mask(dt_fb, p, qc_genes)
	                    if (!length(mask_fb) || !any(mask_fb)) {
	                        stop("No molecules available after filters; cannot compute bounds.")
	                    }

	                    fn_all <- attr(mask_fb, "feature_name")
	                    if (is.null(fn_all) || length(fn_all) != nrow(dt_fb)) {
	                        stop("Unable to resolve transcript feature names for fallback bounds scan; check Xenium codebook/feature columns.")
	                    }
	                    x_fb <- suppressWarnings(as.numeric(dt_fb$x_location[mask_fb]))
	                    y_fb <- suppressWarnings(as.numeric(dt_fb$y_location[mask_fb]))
	                    fn_fb <- as.character(fn_all[mask_fb])
	                    ok_fb <- is.finite(x_fb) & is.finite(y_fb) & !is.na(fn_fb) & nzchar(fn_fb)
	                    if (!any(ok_fb)) {
	                        stop("No molecules available after filters; cannot compute bounds.")
	                    }
	                    x_fb <- x_fb[ok_fb]
	                    y_fb <- y_fb[ok_fb]
	                    fn_fb <- fn_fb[ok_fb]

	                    bounds_fb <- data.frame(
	                        xmin = suppressWarnings(min(x_fb, na.rm = TRUE)),
	                        xmax = suppressWarnings(max(x_fb, na.rm = TRUE)),
	                        ymin = suppressWarnings(min(y_fb, na.rm = TRUE)),
	                        ymax = suppressWarnings(max(y_fb, na.rm = TRUE))
	                    )
	                    if (!all(is.finite(as.numeric(bounds_fb[1, ])))) {
	                        stop("No molecules available after filters; cannot compute bounds.")
	                    }

	                    cache_dt <- data.table::data.table(
	                        x_location = x_fb,
	                        y_location = y_fb,
	                        feature_name = fn_fb
	                    )
	                    state$datasets$transcripts_cache <- cache_dt
	                    state$datasets$transcripts_cache_masked <- TRUE
	                    state$datasets$transcripts_cache_source <- fb$source
	                    state$datasets$transcripts_cache_path <- fb$path
	                    .log_info(
	                        parent,
	                        "S08",
	                        paste0(
	                            "fallback bounds ok; cached transcripts rows=", format(nrow(cache_dt), big.mark = ","),
	                            " source=", fb$source
	                        ),
	                        verbose
	                    )

	                    bounds_fb
	                }
	            }
	        } else {
	            cols_needed <- c("x_location", "y_location")
	            if (!is.null(p$filter_genes) || isTRUE(p$exclude_qc_unsuitable) || isTRUE(p$filtermolecule)) {
	                cols_needed <- c(cols_needed, .transcripts_feature_project_cols(spec))
            }
            if (!is.null(p$max_dist_mol_nuc)) cols_needed <- c(cols_needed, "nucleus_distance")
            if (!is.null(p$filterqv)) cols_needed <- c(cols_needed, "qv")
            ds_bounds <- tryCatch(
                dplyr::select(ds, dplyr::all_of(unique(cols_needed))),
                error = function(e) ds
            )
            ds_bounds |>
                summarise(
                    xmin = min(x_location), xmax = max(x_location),
                    ymin = min(y_location), ymax = max(y_location)
                ) |>
                collect()
        }
    }, error = function(e) {
        stop("Error computing bounds: ", e$message)
    })

    # Auto-detect px/umm mismatches between centroids and transcripts using experiment.xenium pixel_size.
    if (requireNamespace("dplyr", quietly = TRUE)) {
        harm <- .xenium_harmonize_coordinate_units(state, ds, bounds_global, verbose)
        state <- harm$state
        ds <- harm$ds
        bounds_global <- harm$bounds_global
    } else if (isTRUE(getOption("geneSCOPE.debug_coords", FALSE))) {
        .log_info(parent, "S08", "note: package 'dplyr' missing; skipping auto unit harmonization", verbose)
    }

    y_max <- as.numeric(bounds_global$ymax)
    .log_info(
        parent,
        "S08",
        paste0(
            "bounds x=[", round(bounds_global$xmin), ",", round(bounds_global$xmax),
            "] y=[", round(bounds_global$ymin), ",", round(bounds_global$ymax), "]"
        ),
        verbose
    )

    if (isTRUE(getOption("geneSCOPE.debug_coords", FALSE))) {
        cen_rng_pre <- .summarize_xy_ranges(state$objects$centroids, x_col = "x", y_col = "y")
        bdx <- suppressWarnings(as.numeric(bounds_global$xmax) - as.numeric(bounds_global$xmin))
        bdy <- suppressWarnings(as.numeric(bounds_global$ymax) - as.numeric(bounds_global$ymin))
        cdx <- cen_rng_pre[["dx"]]
        cdy <- cen_rng_pre[["dy"]]
        sx <- if (is.finite(bdx) && is.finite(cdx) && cdx > 0) bdx / cdx else NA_real_
        sy <- if (is.finite(bdy) && is.finite(cdy) && cdy > 0) bdy / cdy else NA_real_
        .log_info(
            parent,
            "S08",
            paste0(
                "centroids(pre-flip) x=[", format(cen_rng_pre[["xmin"]], digits = 8, trim = TRUE),
                ",", format(cen_rng_pre[["xmax"]], digits = 8, trim = TRUE),
                "] y=[", format(cen_rng_pre[["ymin"]], digits = 8, trim = TRUE),
                ",", format(cen_rng_pre[["ymax"]], digits = 8, trim = TRUE),
                "] span_ratio(bounds/centroids) x=", format(sx, digits = 6, trim = TRUE),
                " y=", format(sy, digits = 6, trim = TRUE)
            ),
            verbose
        )
        if ((is.finite(sx) && (sx > 50 || sx < 0.02)) || (is.finite(sy) && (sy > 50 || sy < 0.02))) {
            warning(
                paste0(
                    "Centroid vs transcript coordinate spans differ greatly (bounds/centroids x=",
                    format(sx, digits = 6, trim = TRUE),
                    " y=", format(sy, digits = 6, trim = TRUE),
                    "). This often indicates a units mismatch (e.g., centroids in pixels but transcripts in microns, or vice versa) ",
                    "or a wrong flip_y setting."
                ),
                call. = FALSE
            )
        }
    }

    if (p$flip_y) {
        state$objects$centroids <- .flip_y_coordinates(state$objects$centroids, y_max)
        .log_info(parent, "S08", "applied y-axis flip", verbose)
    }

    if (isTRUE(getOption("geneSCOPE.debug_coords", FALSE))) {
        cen_rng_post <- .summarize_xy_ranges(state$objects$centroids, x_col = "x", y_col = "y")
        .log_info(
            parent,
            "S08",
            paste0(
                "centroids(post-flip) x=[", format(cen_rng_post[["xmin"]], digits = 8, trim = TRUE),
                ",", format(cen_rng_post[["xmax"]], digits = 8, trim = TRUE),
                "] y=[", format(cen_rng_post[["ymin"]], digits = 8, trim = TRUE),
                ",", format(cen_rng_post[["ymax"]], digits = 8, trim = TRUE),
                "] y_max_ref=", format(y_max, digits = 8, trim = TRUE)
            ),
            verbose
        )
    }
    step08$done(paste0("y_max=", round(y_max)))

    state$datasets$ds <- ds
    state$datasets$qc_genes <- qc_genes
    state$bounds <- list(bounds_global = bounds_global, y_max = y_max)
    state
}

#' Assemble Cosmx Sparse Matrix
#' @description
#' Internal helper for `.assemble_cosmx_sparse_matrix`.
#' @param expr_dt Parameter value.
#' @param gene_cols Parameter value.
#' @param column_order Parameter value.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.assemble_cosmx_sparse_matrix <- function(expr_dt, gene_cols, column_order) {
    if (!length(column_order)) {
        stop("No cells available after filtering; cannot build sparse matrix.")
    }

    nz_entries <- lapply(gene_cols, function(gene) {
        nz <- which(!is.na(suppressWarnings(as.numeric(expr_dt[[gene]]))) &
            suppressWarnings(as.numeric(expr_dt[[gene]])) != 0)
        if (!length(nz)) {
            return(NULL)
        }
        list(
            i = rep.int(setNames(seq_along(gene_cols), gene_cols)[[gene]], length(nz)),
            j = match(expr_dt$cell_key, column_order)[nz],
            x = suppressWarnings(as.numeric(expr_dt[[gene]]))[nz]
        )
    })
    nz_entries <- nz_entries[!vapply(nz_entries, is.null, logical(1))]
    if (!length(nz_entries)) {
        sparseMatrix(
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
        sparseMatrix(
            i = as.integer(i_idx),
            j = as.integer(j_idx),
            x = as.numeric(x_val),
            dims = c(length(gene_cols), length(column_order)),
            dimnames = list(gene_cols, column_order)
        )
    }
}

#' Assign Cosmx Cell Keys
#' @description
#' Internal helper for `.assign_cosmx_cell_keys`.
#' @param expr_dt Parameter value.
#' @param cid_col Parameter value.
#' @param fov_col Parameter value.
#' @param id_mode Parameter value.
#' @param centroid_keys Parameter value.
#' @param filter_cells Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param prefix Parameter value.
#' @return Return value used internally.
#' @keywords internal
.assign_cosmx_cell_keys <- function(expr_dt,
                                    cid_col,
                                    fov_col,
                                    id_mode,
                                    centroid_keys,
                                    filter_cells,
                                    verbose,
                                    prefix) {
    if (identical(id_mode, "auto")) {
        if (isTRUE(filter_cells) && length(centroid_keys) &&
            (any(grepl("_", centroid_keys, fixed = TRUE)) || any(grepl("FOV", centroid_keys, ignore.case = TRUE)))) {
            id_mode <- "fov_cell"
        } else {
            id_mode <- "cell_id"
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

    expr_dt[, cell_key := make_key(id_mode, 0L)]

    if (!isTRUE(filter_cells) || is.null(centroid_keys)) {
        column_order <- unique(expr_dt$cell_key)
        .log_info("createSCOPE", "S01", paste0("keeping all ", length(column_order), " cells from matrix (no barcode filter)"), verbose)
        return(list(expr_dt = expr_dt, column_order = column_order))
    }

    candidates <- list(
        list(mode = id_mode, offset = 0L),
        list(mode = if (identical(id_mode, "cell_id")) "fov_cell" else "cell_id", offset = 0L),
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
        overlap <- length(intersect(unique(keys_vec), unique(centroid_keys)))
        if (overlap > best_overlap) {
            best_overlap <- overlap
            best_idx <- i
        }
    }

    best_cand <- candidates[[best_idx]]
    expr_dt[, cell_key := make_key(best_cand$mode, best_cand$offset)]

    keep_cells <- intersect(unique(expr_dt$cell_key), unique(centroid_keys))
    if (!length(keep_cells)) {
        stop("Cell identifiers in expression matrix do not overlap with centroids; check id_mode/FOV encoding.")
    }

    expr_dt <- expr_dt[cell_key %in% keep_cells]
    column_order <- centroid_keys[centroid_keys %in% keep_cells]
    .log_info("createSCOPE", "S01", paste0("cell overlap: ", length(column_order), "/", length(centroid_keys)), verbose)
    list(expr_dt = expr_dt, column_order = column_order)
}

#' Build multi-resolution grid layers.
#' @description
#' Internal helper for `.build_grid_layers`.
#' Generates spatial grids across requested resolutions and writes them onto the
#' scope object.
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
    roi_geom <- if (!is.null(state$roi)) state$roi$geometry else NULL
    roi_source <- if (!is.null(state$roi) && !is.null(state$roi$source)) state$roi$source else if (is.null(roi_geom)) "none" else "file"
    use_counts_list <- is.null(roi_geom) || identical(roi_source, "auto_bbox")
    flip_y <- p$flip_y
    parent <- "createSCOPE"
    step13 <- .log_step(parent, "S13", "build grid layers", p$verbose)
    lg_unique <- sort(unique(p$grid_length))
    step13$enter(paste0("grid_length=", paste(lg_unique, collapse = ",")))
    .log_backend(parent, "S13", "grid_builder", "R", reason = "no_native_grid_builder", verbose = p$verbose)
    omp_info <- tryCatch(.native_openmp_info(), error = function(e) NULL)
    if (!is.null(omp_info)) {
        compiled <- omp_info$compiled_with_openmp
        .log_backend(
            parent,
            "S13",
            "openmp_compiled",
            ifelse(is.na(compiled), "NA", compiled),
            reason = "native_openmp_info",
            verbose = p$verbose
        )
        .log_backend(
            parent,
            "S13",
            "openmp_threads",
            ifelse(is.null(omp_info$omp_max_threads), "NA", omp_info$omp_max_threads),
            reason = "native_openmp_info",
            verbose = p$verbose
        )
    } else {
        .log_backend(parent, "S13", "openmp_compiled", "NA", reason = "native_openmp_info_unavailable", verbose = p$verbose)
    }

    debug_s13 <- isTRUE(getOption("geneSCOPE.debug_s13", FALSE))
    seg_filter_mode <- getOption("geneSCOPE.seg_filter_mode", "centroids")
    seg_filter_mode <- as.character(seg_filter_mode)[1]
    if (is.na(seg_filter_mode) || !nzchar(seg_filter_mode)) seg_filter_mode <- "centroids"
    seg_filter_mode <- match.arg(seg_filter_mode, c("centroids", "vertices"))

    seg_layers <- character()
    seg_dt_all <- NULL
    if ((p$min_seg_points > 0L && identical(seg_filter_mode, "vertices")) ||
        identical(p$grid_keep_mode, "segmentation_hybrid") ||
        identical(p$grid_keep_mode, "vertex_only")) {
        t0 <- if (debug_s13) proc.time()[["elapsed"]] else NULL
        seg_layers <- names(scope_obj@coord)[grepl("^segmentation_", names(scope_obj@coord)) &
            !grepl("_polys_", names(scope_obj@coord))]
        if (length(seg_layers)) {
            seg_dt_all <- rbindlist(scope_obj@coord[seg_layers], use.names = TRUE, fill = TRUE)
        } else {
            seg_dt_all <- data.table(x = numeric(), y = numeric())
        }
        if (debug_s13) {
            dt_ms <- as.integer(round((proc.time()[["elapsed"]] - t0) * 1000))
            .log_info(parent, "S13", paste0("debug: segmentation vertices prepared rows=", nrow(seg_dt_all), " elapsed=", dt_ms, "ms"), p$verbose)
        }
    }

    for (lg in lg_unique) {
        .log_info(parent, "S13", paste0("grid ", formatC(lg, format = "f", digits = 1), " um"), p$verbose)
        x0 <- floor(bounds$xmin / lg) * lg
        y0 <- if (flip_y) 0 else floor(bounds$ymin / lg) * lg

        cnt <- if (use_counts_list) {
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

        t_gene <- if (debug_s13) proc.time()[["elapsed"]] else NULL
        if (isTRUE(p$min_gene_types <= 1) && (is.infinite(p$max_gene_types) || isTRUE(p$max_gene_types >= .Machine$integer.max))) {
            keep_grid <- unique(cnt$grid_id)
        } else {
            gene_div <- cnt[, uniqueN(gene), by = grid_id]
            keep_grid <- gene_div[V1 >= p$min_gene_types & V1 <= p$max_gene_types, grid_id]
        }
        molecule_support <- cnt[, .(molecule_count = sum(count)), by = grid_id]
        cnt <- cnt[grid_id %in% keep_grid]
        if (debug_s13 && !is.null(t_gene)) {
            dt_ms <- as.integer(round((proc.time()[["elapsed"]] - t_gene) * 1000))
            .log_info(parent, "S13", paste0("debug: keep_grid=", length(keep_grid), " cnt_rows=", nrow(cnt), " gene_filter_elapsed=", dt_ms, "ms"), p$verbose)
        }

        if (p$min_seg_points > 0L) {
            t_seg <- if (debug_s13) proc.time()[["elapsed"]] else NULL
            valid_seg <- character()
            if (identical(seg_filter_mode, "centroids")) {
                centroids_dt <- scope_obj@coord$centroids
                if (!is.null(centroids_dt) && nrow(centroids_dt)) {
                    seg_cnt <- as.data.table(centroids_dt)[, .N, by = .(
                        grid_id = paste0(
                            "g",
                            (x - x0) %/% lg + 1L,
                            "_",
                            (y - y0) %/% lg + 1L
                        )
                    )]
                    valid_seg <- seg_cnt[N >= p$min_seg_points, grid_id]
                }
            } else if (!is.null(seg_dt_all) && nrow(seg_dt_all)) {
                seg_cnt <- seg_dt_all[, .N, by = .(
                    grid_id = paste0(
                        "g",
                        (x - x0) %/% lg + 1L,
                        "_",
                        (y - y0) %/% lg + 1L
                    )
                )]
                valid_seg <- seg_cnt[N >= p$min_seg_points, grid_id]
            }
            if (length(valid_seg)) {
                keep_grid <- intersect(keep_grid, valid_seg)
                cnt <- cnt[grid_id %in% keep_grid]
            }
            if (debug_s13 && !is.null(t_seg)) {
                dt_ms <- as.integer(round((proc.time()[["elapsed"]] - t_seg) * 1000))
                .log_info(parent, "S13", paste0("debug: seg_filter_mode=", seg_filter_mode, " valid_seg=", length(valid_seg), " seg_filter_elapsed=", dt_ms, "ms"), p$verbose)
            }
        }

        x_breaks <- seq(x0, bounds$xmax, by = lg)
        if (tail(x_breaks, 1L) < bounds$xmax) x_breaks <- c(x_breaks, bounds$xmax)
        max_y_use <- if (flip_y) y_max else bounds$ymax
        y_breaks <- seq(y0, max_y_use, by = lg)
        if (tail(y_breaks, 1L) < max_y_use) y_breaks <- c(y_breaks, max_y_use)

        xbins_eff <- length(x_breaks) - 1L
        ybins_eff <- length(y_breaks) - 1L
        if (!length(keep_grid)) {
            grid_dt <- data.table(
                gx = integer(),
                gy = integer(),
                grid_id = character(),
                xmin = numeric(),
                xmax = numeric(),
                ymin = numeric(),
                ymax = numeric(),
                center_x = numeric(),
                center_y = numeric(),
                width = numeric(),
                height = numeric(),
                idx = integer()
            )
        } else {
            keep_ids <- unique(as.character(keep_grid))
            parts <- data.table::tstrsplit(sub("^g", "", keep_ids), "_", fixed = TRUE)
            gx <- suppressWarnings(as.integer(parts[[1]]))
            gy <- suppressWarnings(as.integer(parts[[2]]))
            grid_dt <- data.table(gx = gx, gy = gy)
            grid_dt <- grid_dt[is.finite(gx) & is.finite(gy)]
            if (nrow(grid_dt)) {
                grid_dt <- grid_dt[gx >= 1L & gx <= xbins_eff & gy >= 1L & gy <= ybins_eff]
            }
            if (!nrow(grid_dt)) {
                grid_dt <- data.table(
                    gx = integer(),
                    gy = integer(),
                    grid_id = character(),
                    xmin = numeric(),
                    xmax = numeric(),
                    ymin = numeric(),
                    ymax = numeric(),
                    center_x = numeric(),
                    center_y = numeric(),
                    width = numeric(),
                    height = numeric(),
                    idx = integer()
                )
            } else {
                data.table::setorder(grid_dt, gx, gy)
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
                    idx = as.integer((gx - 1L) * ybins_eff + gy)
                )]
                if (!p$keep_partial_grid) {
                    grid_dt <- grid_dt[abs(width - lg) < 1e-6 & abs(height - lg) < 1e-6]
                }
            }
        }
        grid_keep_initial_n <- nrow(grid_dt)
        center_min_molecules_effective <- NA_integer_
        if (nrow(grid_dt) && !identical(p$grid_keep_mode, "molecule_driven")) {
            t_keep <- if (debug_s13) proc.time()[["elapsed"]] else NULL
            valid_centroid <- character()
            centroids_dt <- scope_obj@coord$centroids
            if (!is.null(centroids_dt) && nrow(centroids_dt)) {
                centroid_cnt <- as.data.table(centroids_dt)[, .N, by = .(
                    grid_id = paste0(
                        "g",
                        (x - x0) %/% lg + 1L,
                        "_",
                        (y - y0) %/% lg + 1L
                    )
                )]
                valid_centroid <- centroid_cnt[N >= 1L, grid_id]
            }
            valid_vertex <- character()
            if (!is.null(seg_dt_all) && nrow(seg_dt_all)) {
                vertex_cnt <- seg_dt_all[, .N, by = .(
                    grid_id = paste0(
                        "g",
                        (x - x0) %/% lg + 1L,
                        "_",
                        (y - y0) %/% lg + 1L
                    )
                )]
                valid_vertex <- vertex_cnt[N >= 1L, grid_id]
            }
            center_keep <- character()
            if (identical(p$grid_keep_mode, "segmentation_hybrid") && !is.null(roi_geom)) {
                center_dt <- grid_dt[, .(grid_id, x = center_x, y = center_y)]
                center_inside <- tryCatch(
                    .clip_points_to_region(
                        center_dt,
                        roi_geom,
                        chunk_size = max(10000L, as.integer(p$chunk_pts)),
                        workers = 1L,
                        parallel_backend = "serial"
                    ),
                    error = function(e) {
                        .log_info(parent, "S13", paste0("grid center ROI test failed; center fallback disabled: ", conditionMessage(e)), p$verbose)
                        center_dt[0]
                    }
                )
                if (nrow(center_inside)) {
                    center_candidate_ids <- setdiff(center_inside$grid_id, valid_vertex)
                    center_candidate_support <- molecule_support[grid_id %in% center_candidate_ids, molecule_count]
                    if (is.character(p$grid_center_min_molecules)) {
                        q_txt <- sub("^q", "", tolower(p$grid_center_min_molecules))
                        q_prob <- suppressWarnings(as.numeric(q_txt) / 100)
                        if (!is.finite(q_prob)) q_prob <- 0.80
                        q_prob <- min(max(q_prob, 0), 1)
                        if (length(center_candidate_support)) {
                            center_min_molecules_effective <- as.integer(ceiling(stats::quantile(
                                center_candidate_support,
                                probs = q_prob,
                                na.rm = TRUE,
                                names = FALSE,
                                type = 7
                            )))
                        } else {
                            center_min_molecules_effective <- .Machine$integer.max
                        }
                    } else {
                        center_min_molecules_effective <- suppressWarnings(as.integer(p$grid_center_min_molecules))
                    }
                    if (!is.finite(center_min_molecules_effective) || center_min_molecules_effective < 1L) {
                        center_min_molecules_effective <- 1L
                    }
                    center_keep <- intersect(
                        center_candidate_ids,
                        molecule_support[molecule_count >= center_min_molecules_effective, grid_id]
                    )
                }
            }
            keep_by_mode <- switch(
                p$grid_keep_mode,
                centroid_only = valid_centroid,
                vertex_only = valid_vertex,
                segmentation_hybrid = union(valid_vertex, center_keep),
                character()
            )
            grid_dt <- if (length(keep_by_mode)) grid_dt[grid_id %in% keep_by_mode] else grid_dt[0]
            .log_info(
                parent,
                "S13",
                paste0(
                    "grid keep mode=", p$grid_keep_mode,
                    " initial=", grid_keep_initial_n,
                    " final=", nrow(grid_dt),
                    " vertex=", length(valid_vertex),
                    " center=", length(center_keep),
                    " center_min_molecules=", p$grid_center_min_molecules,
                    " center_min_molecules_effective=", center_min_molecules_effective
                ),
                p$verbose
            )
            if (debug_s13 && !is.null(t_keep)) {
                dt_ms <- as.integer(round((proc.time()[["elapsed"]] - t_keep) * 1000))
                .log_info(parent, "S13", paste0("debug: grid_keep_elapsed=", dt_ms, "ms"), p$verbose)
            }
        }
        cnt <- cnt[grid_id %in% grid_dt$grid_id]

        scope_obj@grid[[paste0("grid", lg)]] <- list(
            grid_info = grid_dt,
            counts = cnt,
            grid_length = lg,
            xbins_eff = xbins_eff,
            ybins_eff = ybins_eff,
            grid_keep = list(
                mode = p$grid_keep_mode,
                center_min_molecules = p$grid_center_min_molecules,
                center_min_molecules_effective = center_min_molecules_effective,
                initial_bins = grid_keep_initial_n,
                final_bins = nrow(grid_dt)
            ),
            histology = list(
                lowres = NULL,
                hires = NULL
            )
        )
    }
    grid_sizes <- vapply(scope_obj@grid, function(g) nrow(g$grid_info), integer(1))
    step13$done(paste0("grid_layers=", length(lg_unique), " total_bins=", sum(grid_sizes)))
    state$objects$scope_obj <- scope_obj
    state
}

#' Build Singlecell Config
#' @description
#' Internal helper for `.build_singlecell_config`.
#' @param scope_obj A `scope_object`.
#' @param platform Parameter value.
#' @param xenium_dir Filesystem path.
#' @param cosmx_root Parameter value.
#' @param filter_genes Parameter value.
#' @param exclude_prefix Parameter value.
#' @param filter_barcodes Parameter value.
#' @param id_mode Parameter value.
#' @param force_all_genes Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
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

    parent <- "createSCOPE"

    centroid_cells <- NULL
    if (filter_cells) {
        centroid_cells <- scope_obj@coord$centroids$cell
        if (!length(centroid_cells)) {
            stop("scope_obj@coord$centroids is empty, cannot determine cells to keep.")
        }
        .log_info(parent, "S01", paste0("target cells: ", length(centroid_cells)), verbose)
    } else {
        .log_info(parent, "S01", "barcode filtering disabled - including all cells from source matrix", verbose)
    }

    list(
        scope_obj = scope_obj,
        platform = platform,
        verbose = verbose,
        prefix = parent,
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

#' Clip Points Within Roi
#' @description
#' Internal helper for `.clip_points_within_roi`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.clip_points_within_roi <- function(state) {
    p <- state$params
    roi <- state$roi$geometry
    roi_source <- if (!is.null(state$roi$source)) state$roi$source else if (is.null(roi)) "none" else "file"
    roi_bbox <- state$roi$bbox
    centroids_dt <- state$objects$centroids
    parent <- "createSCOPE"
    step10 <- .log_step(parent, "S10", "clip centroids to ROI", p$verbose)
    if (nrow(centroids_dt) &&
        identical(roi_source, "auto_bbox") &&
        !is.null(roi_bbox) &&
        all(c("xmin", "xmax", "ymin", "ymax") %in% names(roi_bbox))) {
        step10$enter(paste0("cells=", nrow(centroids_dt)))
        bb <- as.numeric(roi_bbox[c("xmin", "xmax", "ymin", "ymax")])
        if (all(is.finite(bb)) && bb[1] < bb[2] && bb[3] < bb[4]) {
            centroids_dt <- centroids_dt[
                x > bb[1] & x < bb[2] &
                    y > bb[3] & y < bb[4]
            ]
        }
        step10$done(paste0("retained=", nrow(centroids_dt)))
    } else if (!is.null(roi) && nrow(centroids_dt)) {
        step10$enter(paste0("cells=", nrow(centroids_dt)))
        centroids_dt <- tryCatch(
            .clip_points_to_region(
                centroids_dt,
                roi,
                chunk_size = p$chunk_pts,
                workers = state$thread_context$ncores_safe,
                parallel_backend = p$parallel_backend
            ),
            error = function(e) stop("Error clipping centroids: ", e$message)
        )
        step10$done(paste0("retained=", nrow(centroids_dt)))
    } else {
        step10$enter("skipped")
        step10$done("skipped")
    }
    state$objects$centroids <- centroids_dt
    state$objects$keep_cells <- if (nrow(centroids_dt)) unique(centroids_dt$cell) else NULL
    state
}

#' Configure Worker Threads
#' @description
#' Internal helper for `.configure_worker_threads`.
#' @param state Parameter value.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.configure_worker_threads <- function(state) {
    p <- state$params
    parent <- "createSCOPE"
    step02 <- .log_step(parent, "S02", "configure runtime resources and threads", p$verbose)

    os_type <- .detect_os()
    max_cores <- max(1L, detectCores())
    requested_cores <- if (is.null(p$ncores) || is.na(p$ncores)) max_cores else p$ncores
    target_cores <- min(requested_cores, max_cores)
    step02$enter(paste0("requested_ncores=", requested_cores, " max_cores=", max_cores))

    sys_mem_gb <- tryCatch({
        if (os_type == "linux") {
            mem_info <- readLines("/proc/meminfo")
            mem_total <- grep("MemTotal", mem_info, value = TRUE)
            mem_kb <- as.numeric(gsub("[^0-9]", "", mem_total))
            host_mem_gb <- mem_kb / 1024 / 1024
            if (host_mem_gb > 50000) warning("Detected unusually high memory: ", round(host_mem_gb, 1), "GB. This may indicate a parsing error.")

            cgroup_limit_bytes <- tryCatch({
                candidates <- c(
                    "/sys/fs/cgroup/memory.max", # cgroup v2
                    "/sys/fs/cgroup/memory/memory.limit_in_bytes", # cgroup v1
                    "/sys/fs/cgroup/memory.limit_in_bytes"
                )
                limit_bytes <- NA_real_
                for (cand in candidates) {
                    if (!file.exists(cand)) next
                    val <- trimws(readLines(cand, warn = FALSE)[1])
                    if (!length(val) || is.na(val) || !nzchar(val)) next
                    if (identical(val, "max")) next
                    num <- suppressWarnings(as.numeric(val))
                    if (is.finite(num) && num > 0 && num < 1e18) {
                        limit_bytes <- num
                        break
                    }
                }
                limit_bytes
            }, error = function(e) NA_real_)
            cgroup_mem_gb <- if (is.finite(cgroup_limit_bytes)) cgroup_limit_bytes / 1024^3 else NA_real_

            slurm_mem_mb <- suppressWarnings(as.numeric(Sys.getenv("SLURM_MEM_PER_NODE", NA_character_)))
            slurm_mem_gb <- if (is.finite(slurm_mem_mb) && slurm_mem_mb > 0) slurm_mem_mb / 1024 else NA_real_

            candidates_gb <- c(host_mem_gb)
            if (is.finite(cgroup_mem_gb) && cgroup_mem_gb > 0) candidates_gb <- c(candidates_gb, cgroup_mem_gb)
            if (is.finite(slurm_mem_gb) && slurm_mem_gb > 0) candidates_gb <- c(candidates_gb, slurm_mem_gb)
            min(candidates_gb)
        } else if (os_type == "macos") {
            mem_bytes <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
            mem_bytes / 1024^3
        } else {
            32
        }
    }, error = function(e) 32)

    # Aggressive mode: attempt the requested core count first, only backing off if cluster setup fails.
    test_cores <- target_cores
    ncores_safe <- 1L

    for (test_nc in seq(test_cores, 1, by = -1L)) {
        cl <- NULL
        ok <- tryCatch({
            if (test_nc > 1 && sys_mem_gb <= 100) {
                cl <- makeCluster(min(test_nc, 4))
                clusterEvalQ(cl, 1 + 1)
            }
            TRUE
        }, error = function(e) {
            .log_info(parent, "S02", paste0("testing ", test_nc, " cores failed, trying fewer: ", e$message), p$verbose)
            FALSE
        }, finally = {
            if (!is.null(cl)) try(stopCluster(cl), silent = TRUE)
        })
        if (isTRUE(ok)) {
            ncores_safe <- test_nc
            break
        }
    }

    mem_display <- if (sys_mem_gb >= 1024) paste0(round(sys_mem_gb / 1024, 1), "TB") else paste0(round(sys_mem_gb, 1), "GB")
    if (ncores_safe < target_cores) {
        .log_info(parent, "S02", paste0("core configuration: using ", ncores_safe, "/", target_cores, " cores (", mem_display, " RAM, ", os_type, ")"), p$verbose)
    } else {
        .log_info(parent, "S02", paste0("core configuration: using ", ncores_safe, " cores (", mem_display, " RAM, ", os_type, ")"), p$verbose)
    }

    omp_info <- tryCatch(.native_openmp_info(), error = function(e) NULL)
    if (!is.null(omp_info)) {
        .log_info(
            parent,
            "S02",
            paste0(
                "openmp compiled=",
                ifelse(is.na(omp_info$compiled_with_openmp), "NA", omp_info$compiled_with_openmp),
                " omp_max_threads=", ifelse(is.null(omp_info$omp_max_threads), "NA", omp_info$omp_max_threads),
                " omp_num_procs=", ifelse(is.null(omp_info$omp_num_procs), "NA", omp_info$omp_num_procs)
            ),
            p$verbose
        )
    } else {
        .log_info(parent, "S02", "openmp info unavailable", p$verbose)
    }

    thread_config <- .configure_threads_for("io_bound", ncores_safe, restore_after = TRUE)
    restore_fn <- attr(thread_config, "restore_function")
    ncores_io <- thread_config$r_threads
    ncores_cpp <- thread_config$openmp_threads

    state$thread_context <- list(
        thread_config = thread_config,
        restore_fn = restore_fn,
        ncores_safe = ncores_safe,
        ncores_io = ncores_io,
        ncores_cpp = ncores_cpp,
        os_type = os_type,
        sys_mem_gb = sys_mem_gb
    )
    step02$done(paste0("ncores_safe=", ncores_safe, " io_threads=", ncores_io, " omp_threads=", ncores_cpp))
    state
}

#' Construct a `dgCMatrix` from Xenium components.
#' @description
#' Internal helper for `.construct_xenium_sparse_matrix`.
#' @param components List returned by `.read_xenium_sparse_components()`.
#' @param prefix Log prefix for user messages.
#' @param verbose Logical controlling verbosity.
#' @return List with the counts matrix and optional gene map.
#' @keywords internal
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

    .log_info(
        "createSCOPE",
        "S01",
        paste0(
            "matrix: ", nrow_h5, "x", ncol_h5,
            " (", length(genes), " genes, ", length(barcodes), " cells)"
        ),
        verbose
    )

    counts <- NULL
    if (col_is_cell) {
        counts <- new(
            "dgCMatrix",
            Dim = c(length(genes), length(barcodes)),
            x = as.numeric(components$data),
            i = as.integer(components$indices),
            p = indptr,
            Dimnames = list(make.unique(genes), barcodes)
        )
    } else {
        tmp <- new(
            "dgCMatrix",
            Dim = c(length(barcodes), length(genes)),
            x = as.numeric(components$data),
            i = as.integer(components$indices),
            p = indptr,
            Dimnames = list(barcodes, make.unique(genes))
        )
        .log_info("createSCOPE", "S01", "transposing matrix layout", verbose)
        counts <- t(tmp)
    }

    list(
        counts = counts,
        gene_map = components$gene_map
    )
}

#' Filter Gene Panel
#' @description
#' Internal helper for `.filter_gene_panel`.
#' @param counts Parameter value.
#' @param filters Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param prefix Parameter value.
#' @return Return value used internally.
#' @keywords internal
.filter_gene_panel <- function(counts, filters, verbose, prefix) {
    gene_names <- rownames(counts)
    keep_genes <- gene_names

    exclude_prefix <- filters$gene_exclude_prefix
    if (length(exclude_prefix)) {
        pattern <- paste0("^(?:", paste(exclude_prefix, collapse = "|"), ")")
        drop_mask <- grepl(pattern, keep_genes, perl = TRUE)
        if (any(drop_mask)) {
            keep_genes <- keep_genes[!drop_mask]
            .log_info(
                "createSCOPE",
                "S01",
                paste0(
                    "excluded ", sum(drop_mask), " genes with prefixes: ",
                    paste(exclude_prefix, collapse = ", ")
                ),
                verbose
            )
        }
    }

    allowlist <- filters$gene_allowlist
    if (!is.null(allowlist)) {
        keep_genes <- intersect(keep_genes, allowlist)
        .log_info("createSCOPE", "S01", paste0("filtered to ", length(keep_genes), " genes from target list"), verbose)
    }

    if (!length(keep_genes)) {
        stop("No genes remain after filtering. Please check exclude_prefix and filter_genes parameters.")
    }

    if (!identical(keep_genes, gene_names)) {
        counts <- counts[keep_genes, , drop = FALSE]
    }

    .log_info("createSCOPE", "S01", paste0("final gene count: ", nrow(counts), "/", length(gene_names)), verbose)

    list(counts = counts, genes = keep_genes)
}

#' Finalize Output Object
#' @description
#' Internal helper for `.finalize_output_object`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.finalize_output_object <- function(state) {
    scope_obj <- state$objects$scope_obj
    parent <- "createSCOPE"
    step14 <- .log_step(parent, "S14", "finalize scope object", state$params$verbose)
    step14$enter()
    platform_label <- "geneSCOPE"
    if (length(scope_obj@grid) == 0 ||
        all(vapply(scope_obj@grid, function(g) nrow(g$grid_info) == 0L, logical(1)))) {
        warning("No effective grid layer generated - please check filters/parameters.")
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
    if (!is.null(state$experiment) && is.list(state$experiment)) {
        vstr <- .format_xenium_version(state$experiment)
        gen <- .infer_xenium_generation(state$experiment)
        pix <- suppressWarnings(as.numeric(state$experiment$pixel_size_um))
        scope_obj@stats$xenium_version <- vstr
        scope_obj@stats$xenium_generation <- gen
        scope_obj@stats$xenium_pixel_size_um <- pix
        scope_obj@stats$xenium_chemistry_version <- as.character(state$experiment$chemistry_version)
        scope_obj@stats$xenium_run_name <- as.character(state$experiment$run_name)
        scope_obj@stats$xenium_region_name <- as.character(state$experiment$region_name)
    }
    if (!is.null(state$coord_units) && is.list(state$coord_units)) {
        scope_obj@stats$xenium_coord_units_action <- as.character(state$coord_units$action)
        scope_obj@stats$xenium_coord_units_scale <- suppressWarnings(as.numeric(state$coord_units$scale))
        scope_obj@stats$xenium_coord_units_ratio_before <- suppressWarnings(as.numeric(state$coord_units$ratio_before))
    }
    cell_count <- if (!is.null(scope_obj@coord$centroids)) nrow(scope_obj@coord$centroids) else 0L
    grid_layers <- length(scope_obj@grid)
    step14$done(paste0("cells=", cell_count, " grid_layers=", grid_layers))
    state$objects$scope_obj <- scope_obj
    state
}

#' Guess Singlecell Platform
#' @description
#' Internal helper for `.guess_singlecell_platform`.
#' @param path Parameter value.
#' @return Return value used internally.
#' @keywords internal
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

#' Identify Cosmx Gene Columns
#' @description
#' Internal helper for `.identify_cosmx_gene_columns`.
#' @param cols Parameter value.
#' @param exclude_prefix Parameter value.
#' @return Return value used internally.
#' @keywords internal
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

#' Ingest Segmentation Geometries
#' @description
#' Internal helper for `.ingest_segmentation_geometries`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.ingest_segmentation_geometries <- function(state) {
    p <- state$params
    seg_files <- state$paths$seg_files
    scope_obj <- state$objects$scope_obj
    keep_cells <- state$objects$keep_cells
    y_max <- state$bounds$y_max
    parent <- "createSCOPE"
    step11 <- .log_step(parent, "S11", "ingest segmentation data", p$verbose)

    if (!length(seg_files)) {
        step11$enter("skipped")
        step11$done("skipped")
        state$objects$scope_obj <- scope_obj
        return(state)
    }

    step11$enter(paste0("seg_type=", p$seg_type, " files=", length(seg_files)))
    cell_key_mode <- if (!is.null(state$objects$cell_key_mode)) as.character(state$objects$cell_key_mode)[1] else "cell_id"
    cell_fov_col <- if (!is.null(state$objects$cell_fov_column)) as.character(state$objects$cell_fov_column)[1] else NULL
    for (tag in names(seg_files)) {
        .log_info(parent, "S11", paste0("processing ", tag, " segmentation"), p$verbose)
        seg_res <- tryCatch({
            .build_segmentation_geometries(
                path = seg_files[[tag]],
                label = tag,
                flip = p$flip_y,
                y_max = y_max,
                keep_ids = keep_cells,
                workers = state$thread_context$ncores_safe,
                cell_key_mode = cell_key_mode,
                fov_column = cell_fov_col
            )
        }, error = function(e) {
            warning("Error processing ", tag, " segmentation: ", e$message)
            list(points = data.table(), polygons = list())
        })
        scope_obj@coord[[paste0("segmentation_", tag)]] <- seg_res$points
        scope_obj@coord[[paste0("segmentation_polys_", tag)]] <- seg_res$polygons
    }
    step11$done("ok")
    state$objects$scope_obj <- scope_obj
    state
}

#' Initialize Dataset Builder State
#' @description
#' Internal helper for `.initialize_dataset_builder_state`.
#' @param input_dir Filesystem path.
#' @param grid_length Parameter value.
#' @param seg_type Parameter value.
#' @param filter_genes Parameter value.
#' @param max_dist_mol_nuc Numeric threshold.
#' @param filtermolecule Parameter value.
#' @param exclude_prefix Parameter value.
#' @param filterqv Parameter value.
#' @param coord_file Filesystem path.
#' @param ncores Number of cores/threads to use.
#' @param chunk_pts Parameter value.
#' @param min_gene_types Numeric threshold.
#' @param max_gene_types Numeric threshold.
#' @param min_seg_points Numeric threshold.
#' @param keep_partial_grid Logical flag.
#' @param verbose Logical; whether to emit progress messages.
#' @param flip_y Parameter value.
#' @param data_type Parameter value.
#' @return Return value used internally.
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
    parallel_backend = NULL,
    min_gene_types,
    max_gene_types,
    min_seg_points,
        grid_keep_mode = c("segmentation_hybrid", "molecule_driven", "vertex_only", "centroid_only"),
        grid_center_min_molecules = "q80",
    keep_partial_grid,
    verbose,
    flip_y,
    data_type = "xenium",
    exclude_qc_unsuitable = FALSE,
    qc_out_dir = NULL,
    qc_flag = NULL) {
    seg_type <- match.arg(seg_type, c("cell", "nucleus", "both", "none"))
    if (!dir.exists(input_dir)) {
        stop("`input_dir` does not exist: ", input_dir)
    }
    stopifnot(is.numeric(grid_length), length(grid_length) >= 1, all(grid_length > 0))
    stopifnot(is.logical(flip_y), length(flip_y) == 1)
    stopifnot(is.logical(keep_partial_grid), length(keep_partial_grid) == 1)
    if (is.null(parallel_backend)) parallel_backend <- "auto"
    parallel_backend <- as.character(parallel_backend)[1]
    if (is.na(parallel_backend) || !nzchar(parallel_backend)) parallel_backend <- "auto"
    grid_keep_mode <- match.arg(grid_keep_mode)
    grid_center_min_molecules <- grid_center_min_molecules[1]
    if (is.numeric(grid_center_min_molecules)) {
        grid_center_min_molecules <- suppressWarnings(as.integer(grid_center_min_molecules))
        if (!is.finite(grid_center_min_molecules) || grid_center_min_molecules < 1L) {
            stop("grid_center_min_molecules must be a positive integer or a quantile string like 'q80'.", call. = FALSE)
        }
    } else {
        grid_center_min_molecules <- tolower(as.character(grid_center_min_molecules))
        if (!grepl("^q(0|[1-9][0-9]?|100)(\\.[0-9]+)?$", grid_center_min_molecules)) {
            stop("grid_center_min_molecules must be a positive integer or a quantile string like 'q80'.", call. = FALSE)
        }
    }

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
            parallel_backend = parallel_backend,
            min_gene_types = min_gene_types,
            max_gene_types = max_gene_types,
            min_seg_points = min_seg_points,
            grid_keep_mode = grid_keep_mode,
            grid_center_min_molecules = grid_center_min_molecules,
            keep_partial_grid = keep_partial_grid,
            verbose = verbose,
            flip_y = flip_y,
            data_type = data_type,
            exclude_qc_unsuitable = exclude_qc_unsuitable,
            qc_out_dir = qc_out_dir,
            qc_flag = qc_flag
        ),
        thread_context = NULL,
        env_restore = NULL,
        paths = NULL,
        datasets = list(),
        roi = list(),
        objects = list()
    )
}

#' Initialize Output Object
#' @description
#' Internal helper for `.initialize_output_object`.
#' @param state Parameter value.
#' @return Return value used internally.
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

#' Limit Thread Environment
#' @description
#' Internal helper for `.limit_thread_environment`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.limit_thread_environment <- function(state) {
    p <- state$params
    ctx <- state$thread_context
    parent <- "createSCOPE"

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

    omp_threads <- suppressWarnings(as.integer(ctx$ncores_cpp))
    if (!is.finite(omp_threads) || omp_threads <= 0L) omp_threads <- 1L
    r_threads <- suppressWarnings(as.integer(ctx$ncores_io))
    if (!is.finite(r_threads) || r_threads <= 0L) r_threads <- suppressWarnings(as.integer(ctx$ncores_safe))
    if (!is.finite(r_threads) || r_threads <= 0L) r_threads <- 1L

    Sys.setenv(OMP_NUM_THREADS = as.character(omp_threads))
    Sys.setenv(OMP_THREAD_LIMIT = as.character(omp_threads))
    Sys.setenv(
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        BLAS_NUM_THREADS = "1",
        LAPACK_NUM_THREADS = "1"
    )
    if (!is.null(ctx$os_type) && identical(ctx$os_type, "linux")) Sys.setenv(OMP_SCHEDULE = "static")

    if (requireNamespace("data.table", quietly = TRUE)) {
        tryCatch(data.table::setDTthreads(r_threads), error = function(e) {
            .log_info(parent, "S02", "warning: could not set data.table threads", p$verbose)
        })
    }

    state$env_restore <- restore_env
    state
}

#' Load Centroid Table
#' @description
#' Internal helper for `.load_centroid_table`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.load_centroid_table <- function(state) {
    p <- state$params
    ctx <- state$thread_context
    paths <- state$paths
    parent <- "createSCOPE"
    .log_info(parent, "S04", "loading cell centroids", p$verbose)

	    .ensure_arrow_available("[geneSCOPE::.create_scope] Loading centroids")
	    centroids_dt <- .with_data_table_attached(function() {
	        arrow_old_opt <- getOption("arrow.num_threads")
	        arrow_old_cpu <- if (exists("cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
	            tryCatch(arrow::cpu_count(), error = function(e) NA_integer_)
	        } else {
	            NA_integer_
	        }
		        arrow_old_io <- if (exists("io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
		            tryCatch(arrow::io_thread_count(), error = function(e) NA_integer_)
		        } else {
		            NA_integer_
		        }
		        on.exit({
		            try(options(arrow.num_threads = arrow_old_opt), silent = TRUE)

                    restore_cpu <- suppressWarnings(as.integer(arrow_old_cpu))
                    if (!is.finite(restore_cpu) || restore_cpu < 1L) restore_cpu <- suppressWarnings(as.integer(ctx$ncores_io))
                    if (!is.finite(restore_cpu) || restore_cpu < 1L) restore_cpu <- suppressWarnings(as.integer(ctx$ncores_safe))
                    if (!is.finite(restore_cpu) || restore_cpu < 1L) restore_cpu <- 1L

                    restore_io <- suppressWarnings(as.integer(arrow_old_io))
                    if (!is.finite(restore_io) || restore_io < 1L) restore_io <- restore_cpu
                    restore_io <- max(2L, restore_io)

		            if (exists("set_cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
		                try(arrow::set_cpu_count(restore_cpu), silent = TRUE)
		            }
		            if (exists("set_io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
		                try(arrow::set_io_thread_count(restore_io), silent = TRUE)
		            }
		        }, add = TRUE)

	        extract_centroids <- function(cd) {
	            required <- c("cell_id", "x_centroid", "y_centroid")
	            missing <- setdiff(required, names(cd))
	            if (length(missing)) {
                stop("Centroids table missing required columns: ", paste(missing, collapse = ", "))
            }

            fov_col <- NA_character_
            if ("fov" %in% names(cd)) {
                fov_col <- "fov"
            } else if ("fov_name" %in% names(cd)) {
                fov_col <- "fov_name"
            }
            has_dups <- anyDuplicated(as.character(cd[["cell_id"]])) > 0L
            use_fov <- !is.na(fov_col) && isTRUE(has_dups)

            if (use_fov) {
                fov_val <- as.character(cd[[fov_col]])
                fov_val[is.na(fov_val) | !nzchar(fov_val)] <- "NA"
                cd[, cell := paste0(as.character(cell_id), "_", fov_val)]
            } else {
                cd[, cell := as.character(cell_id)]
            }
            out <- cd[, .(cell = cell, x = x_centroid, y = y_centroid)]
            attr(out, "genescope_cell_key_mode") <- if (use_fov) "cell_id_fov" else "cell_id"
            attr(out, "genescope_cell_fov_column") <- if (use_fov) fov_col else NULL
            out
        }

        read_cells_csv <- function(path) {
            if (!file.exists(path)) stop("Cells CSV not found: ", path)
            is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)
            has_zcat <- nzchar(Sys.which("zcat"))
            has_gzcat <- nzchar(Sys.which("gzcat"))
            has_pigz <- nzchar(Sys.which("pigz"))

            dt <- tryCatch(
                fread(path, showProgress = FALSE),
                error = function(e) NULL
            )
            if (is.null(dt) && isTRUE(is_gz)) {
                cmd <- NULL
                if (has_gzcat) {
                    cmd <- paste(Sys.which("gzcat"), shQuote(path))
                } else if (has_zcat) {
                    cmd <- paste(Sys.which("zcat"), shQuote(path))
                } else if (has_pigz) {
                    cmd <- paste(Sys.which("pigz"), "-dc", shQuote(path))
                }
                if (!is.null(cmd) && nzchar(cmd)) {
                    dt <- tryCatch(
                        fread(cmd = cmd, showProgress = FALSE),
                        error = function(e) NULL
                    )
                }
            }
            if (!is.null(dt)) return(dt)

            con <- if (isTRUE(is_gz)) gzfile(path, open = "rt") else file(path, open = "rt")
            on.exit(close(con), add = TRUE)
            df <- utils::read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
            as.data.table(df)
	        }

	        configure_arrow_threads <- function(n) {
	            n <- suppressWarnings(as.integer(n))
	            n <- ifelse(is.finite(n) && n > 0L, n, 1L)
	            io_n <- max(2L, n)
	            tryCatch({
	                options(arrow.num_threads = n)
	                if (exists("set_cpu_count", where = asNamespace("arrow"))) arrow::set_cpu_count(n)
	                if (exists("set_io_thread_count", where = asNamespace("arrow"))) arrow::set_io_thread_count(io_n)
	            }, error = function(e) {
	                .log_info(parent, "S04", paste0("warning: could not configure Arrow threading: ", e$message), p$verbose)
	            })
	        }

        read_cells_parquet <- function(nthreads) {
            configure_arrow_threads(nthreads)
            cell_data <- arrow::read_parquet(paths$cell_parq)
            cd <- as.data.table(cell_data)
            extract_centroids(cd)
        }

        tryCatch(
            read_cells_parquet(ctx$ncores_io),
            error = function(e_parq) {
                msg1 <- conditionMessage(e_parq)
                .log_info(parent, "S04", paste0("cells.parquet read failed: ", msg1), p$verbose)

                centroids_retry <- tryCatch(
                    read_cells_parquet(1),
                    error = function(e_retry) NULL
                )
                if (!is.null(centroids_retry)) {
                    .log_info(parent, "S04", "cells.parquet recovered with Arrow single-thread retry", p$verbose)
                    return(centroids_retry)
                }

                cells_csv_gz <- file.path(p$input_dir, "cells.csv.gz")
                cells_csv <- file.path(p$input_dir, "cells.csv")
                fallback <- if (file.exists(cells_csv_gz)) cells_csv_gz else if (file.exists(cells_csv)) cells_csv else NA_character_
                if (!is.character(fallback) || is.na(fallback)) {
                    stop(
                        "Error reading centroids from cells.parquet: ",
                        msg1,
                        "\nNo CSV fallback found (expected cells.csv.gz or cells.csv under: ",
                        p$input_dir,
                        ")."
                    )
                }
                .log_info(parent, "S04", paste0("falling back to ", basename(fallback)), p$verbose)
                cd <- tryCatch(
                    read_cells_csv(fallback),
                    error = function(e_csv) stop("Error reading centroids from ", basename(fallback), ": ", e_csv$message)
                )
                extract_centroids(cd)
            }
        )
    })

    if (isTRUE(getOption("geneSCOPE.debug_coords", FALSE))) {
        rng <- .summarize_xy_ranges(centroids_dt, x_col = "x", y_col = "y")
        .log_info(
            parent,
            "S04",
            paste0(
                "centroids range x=[", format(rng[["xmin"]], digits = 8, trim = TRUE),
                ",", format(rng[["xmax"]], digits = 8, trim = TRUE),
                "] y=[", format(rng[["ymin"]], digits = 8, trim = TRUE),
                ",", format(rng[["ymax"]], digits = 8, trim = TRUE),
                "]"
            ),
            p$verbose
        )
    }

    if (nrow(centroids_dt) == 0) stop("No centroids found in cells table")
    .log_info(parent, "S04", paste0("found ", nrow(centroids_dt), " cells"), p$verbose)

    cell_key_mode <- attr(centroids_dt, "genescope_cell_key_mode")
    if (is.null(cell_key_mode) || is.na(cell_key_mode) || !nzchar(cell_key_mode)) cell_key_mode <- "cell_id"
    state$objects$cell_key_mode <- cell_key_mode
    state$objects$cell_fov_column <- attr(centroids_dt, "genescope_cell_fov_column")
    attr(centroids_dt, "genescope_cell_key_mode") <- NULL
    attr(centroids_dt, "genescope_cell_fov_column") <- NULL

    state$objects$centroids <- centroids_dt
    state
}

#' Load Dataset Dependencies
#' @description
#' Internal helper for `.load_dataset_dependencies`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.load_dataset_dependencies <- function(state) {
    p <- state$params
    ctx <- state$thread_context
    parent <- "createSCOPE"

	    .log_info(parent, "S03", "initializing data processing environment", p$verbose)
	    .ensure_arrow_available("[geneSCOPE::.create_scope] Dataset dependencies")
	    tryCatch({
	        cpu_threads <- suppressWarnings(as.integer(ctx$ncores_safe))
	        cpu_threads <- ifelse(is.finite(cpu_threads) && cpu_threads > 0L, cpu_threads, 1L)
	        io_threads <- max(2L, cpu_threads)
	        options(arrow.num_threads = cpu_threads)
	        if (exists("set_cpu_count", where = asNamespace("arrow"))) arrow::set_cpu_count(cpu_threads)
	        if (exists("set_io_thread_count", where = asNamespace("arrow"))) arrow::set_io_thread_count(io_threads)
	    }, error = function(e) {
	        .log_info(parent, "S03", "warning: could not configure Arrow threading", p$verbose)
	    })
	    state
	}

#' Locate the CosMx dataset and expression matrix.
#' @description
#' Internal helper for `.locate_cosmx_dataset`.
#' @param cosmx_root Root directory of the CosMx project.
#' @return List containing dataset directory and expression matrix path.
#' @keywords internal
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

#' Open Transcript Dataset
#' @description
#' Internal helper for `.open_transcript_dataset`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.open_transcript_dataset <- function(state) {
    paths <- state$paths
    parent <- "createSCOPE"
    p <- state$params
    ctx <- state$thread_context
    cpu_target <- if (!is.null(ctx) && !is.null(ctx$ncores_safe)) as.integer(ctx$ncores_safe) else NA_integer_
    cpu_target <- ifelse(is.finite(cpu_target) && cpu_target > 0L, cpu_target, as.integer(getOption("arrow.num_threads", 1L)))
    cpu_target <- max(1L, cpu_target)
    cpu_open <- cpu_target
    size_gb <- tryCatch({
        info <- file.info(paths$tf_parq)
        if (is.finite(info$size)) info$size / 1024^3 else NA_real_
    }, error = function(e) NA_real_)
    mem_gb <- if (!is.null(ctx) && !is.null(ctx$sys_mem_gb) && is.finite(ctx$sys_mem_gb)) as.numeric(ctx$sys_mem_gb) else NA_real_
    mem_ratio <- if (is.finite(size_gb) && is.finite(mem_gb) && mem_gb > 0) size_gb / mem_gb else NA_real_

    ncores_arg <- suppressWarnings(as.integer(p$ncores)[1])
    user_requested <- length(ncores_arg) && is.finite(ncores_arg) && ncores_arg > 0L
    auto_throttle <- !user_requested && isTRUE(getOption("geneSCOPE.arrow_auto_throttle", TRUE))
    if (auto_throttle) {
        cap_factor <- NA_real_
        if (is.finite(mem_ratio)) {
            if (mem_ratio >= 0.50) {
                cap_factor <- 0.25
            } else if (mem_ratio >= 0.25) {
                cap_factor <- 0.50
            } else if (mem_ratio >= 0.10) {
                cap_factor <- 0.75
            }
        } else if (is.finite(size_gb)) {
            if (size_gb >= 40) {
                cap_factor <- 0.25
            } else if (size_gb >= 20) {
                cap_factor <- 0.50
            } else if (size_gb >= 5) {
                cap_factor <- 0.75
            }
        }
        if (!is.na(cap_factor)) {
            cpu_open <- suppressWarnings(as.integer(ceiling(cpu_target * cap_factor)))
            if (!is.finite(cpu_open) || cpu_open < 1L) cpu_open <- 1L
            if (cpu_target > 1L) cpu_open <- max(2L, cpu_open)
            cpu_open <- min(cpu_open, cpu_target)
        }
    }

    if (requireNamespace("arrow", quietly = TRUE) && cpu_open < cpu_target) {
        old_opt <- getOption("arrow.num_threads")
        old_cpu <- if (exists("cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
            tryCatch(arrow::cpu_count(), error = function(e) NA_integer_)
        } else {
            NA_integer_
        }
        old_io <- if (exists("io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
            tryCatch(arrow::io_thread_count(), error = function(e) NA_integer_)
        } else {
            NA_integer_
        }
        on.exit({
            try(options(arrow.num_threads = old_opt), silent = TRUE)

            restore_cpu <- suppressWarnings(as.integer(old_cpu))
            if (!is.finite(restore_cpu) || restore_cpu < 1L) restore_cpu <- cpu_target
            restore_cpu <- max(1L, restore_cpu)

            restore_io <- suppressWarnings(as.integer(old_io))
            if (!is.finite(restore_io) || restore_io < 1L) restore_io <- restore_cpu
            restore_io <- max(2L, restore_io)

            if (exists("set_cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
                try(arrow::set_cpu_count(restore_cpu), silent = TRUE)
            }
            if (exists("set_io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
                try(arrow::set_io_thread_count(restore_io), silent = TRUE)
            }
        }, add = TRUE)

        io_open <- max(2L, cpu_open)
        .log_info(
            parent,
            "S04",
            paste0(
                "throttling Arrow threads for open_dataset: cpu=", cpu_open,
                " io=", io_open,
                " (restore cpu=", cpu_target,
                if (is.finite(size_gb)) paste0(", transcripts_size_gb=", round(size_gb, 1)) else "",
                if (is.finite(mem_gb)) paste0(", sys_mem_gb=", round(mem_gb, 1)) else "",
                ")"
            ),
            state$params$verbose
        )
        try(options(arrow.num_threads = cpu_open), silent = TRUE)
        if (exists("set_cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
            try(arrow::set_cpu_count(cpu_open), silent = TRUE)
        }
        if (exists("set_io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
            try(arrow::set_io_thread_count(io_open), silent = TRUE)
        }
    }

	    ds <- tryCatch(arrow::open_dataset(paths$tf_parq), error = function(e) {
	        stop("Error opening transcripts dataset: ", e$message)
	    })
	    .log_info(parent, "S04", paste0("opened transcripts dataset: ", basename(paths$tf_parq)), state$params$verbose)
	    probe <- .probe_transcripts_dataset_query(ds)
	    state$datasets$ds_queryable <- isTRUE(probe$ok)
	    state$datasets$ds_query_probe_error <- if (inherits(probe$error, "error")) conditionMessage(probe$error) else NULL
	    if (!isTRUE(probe$ok)) {
	        .log_info(
	            parent,
	            "S04",
	            paste0(
	                "Arrow query probe failed for transcripts dataset; attempting compatibility rewrite: ",
	                conditionMessage(probe$error)
	            ),
	            state$params$verbose
	        )
	        compat_path <- tryCatch(
	            .ensure_transcripts_arrow_compat_parquet(paths$tf_parq, input_dir = p$input_dir, verbose = p$verbose),
	            error = function(e) e
	        )
	        if (!inherits(compat_path, "error")) {
	            ds_compat <- tryCatch(arrow::open_dataset(compat_path), error = function(e) e)
	            if (!inherits(ds_compat, "error")) {
	                probe_compat <- .probe_transcripts_dataset_query(ds_compat)
	                if (isTRUE(probe_compat$ok)) {
	                    .log_info(parent, "S04", paste0("using Arrow-compatible rewritten transcripts: ", basename(compat_path)), p$verbose)
	                    state$paths$tf_parq_original <- paths$tf_parq
	                    state$paths$tf_parq <- compat_path
	                    ds <- ds_compat
	                    state$datasets$ds_queryable <- TRUE
	                    state$datasets$ds_query_probe_error <- NULL
	                } else {
	                    .log_info(
	                        parent,
	                        "S04",
	                        paste0(
	                            "compatibility rewrite created dataset but Arrow query probe still failed: ",
	                            conditionMessage(probe_compat$error)
	                        ),
	                        p$verbose
	                    )
	                    state$datasets$ds_query_probe_error <- conditionMessage(probe_compat$error)
	                }
	            } else {
	                .log_info(
	                    parent,
	                    "S04",
	                    paste0("compatibility rewrite completed but open_dataset failed: ", conditionMessage(ds_compat)),
	                    p$verbose
	                )
	            }
	        } else {
	            .log_info(
	                parent,
	                "S04",
	                paste0("compatibility rewrite unavailable; will use direct fallback reads when needed: ", conditionMessage(compat_path)),
	                p$verbose
	            )
	        }
	    }
	    state$datasets$ds <- ds
	    state$datasets$ds_raw <- ds
	    state
	}

#' Prefetch Roi Molecules
#' @description
#' Internal helper for `.prefetch_roi_molecules`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.prefetch_roi_molecules <- function(state) {
    p <- state$params
    roi <- state$roi$geometry
    roi_source <- if (!is.null(state$roi$source)) state$roi$source else if (is.null(roi)) "none" else "file"
    ds <- state$datasets$ds
    if (is.null(p$transcripts_spec) || !is.list(p$transcripts_spec)) {
        state$params$transcripts_spec <- .resolve_transcripts_feature_spec(
            ds = ds,
            input_dir = p$input_dir,
            verbose = p$verbose
        )
        p <- state$params
    }
    spec <- p$transcripts_spec
    qc_genes <- state$datasets$qc_genes
    if (is.null(qc_genes)) qc_genes <- character()
    y_max <- state$bounds$y_max
    mol_small <- NULL
    parent <- "createSCOPE"
    step12 <- .log_step(parent, "S12", "prefetch ROI molecules", p$verbose)
	    if (!is.null(roi) && identical(roi_source, "file")) {
	        step12$enter("roi=provided")
	        cols <- c("x_location", "y_location", .transcripts_feature_project_cols(spec))
	        if (!is.null(p$max_dist_mol_nuc)) cols <- c(cols, "nucleus_distance")
	        if (!is.null(p$filterqv)) cols <- c(cols, "qv")
	        cols <- unique(cols)
	        mol_small <- NULL
	        use_direct_fallback <- isFALSE(state$datasets$ds_queryable)
	        if (!use_direct_fallback) {
	            mol_small <- tryCatch({
	                ds |>
	                    dplyr::select(dplyr::all_of(cols)) |>
	                    collect() |>
	                    as.data.table()
	            }, error = function(e) {
	                .log_info(parent, "S12", paste0("ROI prefetch Arrow collect failed; falling back to direct read: ", conditionMessage(e)), p$verbose)
	                NULL
	            })
	        }
	        if (is.null(mol_small)) {
	            fb <- tryCatch(
	                .read_transcripts_fallback_dt(state$paths, p$input_dir, cols, verbose = p$verbose, use_threads = TRUE),
	                error = function(e) e
	            )
	            if (inherits(fb, "error")) {
	                .log_info(parent, "S12", paste0("ROI prefetch fallback failed; retrying single-thread: ", conditionMessage(fb)), p$verbose)
	                fb <- tryCatch(
	                    .read_transcripts_fallback_dt(state$paths, p$input_dir, cols, verbose = p$verbose, use_threads = FALSE),
	                    error = function(e) e
	                )
	            }
	            if (inherits(fb, "error")) {
	                warning("Error pre-fetching molecules: ", conditionMessage(fb))
	                mol_small <- data.table(x = numeric(), y = numeric(), feature_name = character())
	            } else {
	                mol_small <- data.table::as.data.table(fb$dt)
	            }
	        }
	        if (nrow(mol_small) > 0) {
	            if ("x_location" %in% names(mol_small)) setnames(mol_small, "x_location", "x")
	            if ("y_location" %in% names(mol_small)) setnames(mol_small, "y_location", "y")
	            mask <- .ingest_molecule_filter_mask(mol_small, p, qc_genes)
            fn <- attr(mask, "feature_name")
            if (!is.null(fn) && length(fn) == nrow(mol_small)) {
                mol_small[, feature_name := fn]
            }
            if (length(mask)) {
                mol_small <- mol_small[mask, .(x, y, feature_name)]
            } else {
                mol_small <- mol_small[0, .(x, y, feature_name)]
            }
        }
        if (p$flip_y && nrow(mol_small) > 0) mol_small[, y := y_max - y]
        if (nrow(mol_small) > 0) {
            mol_small <- .clip_points_to_region(
                mol_small,
                roi,
                chunk_size = p$chunk_pts,
                workers = state$thread_context$ncores_safe,
                parallel_backend = p$parallel_backend
            )
	        }
        step12$done(paste0("retained=", nrow(mol_small)))
    } else {
        step12$enter("skipped")
        step12$done("skipped")
    }
    state$objects$mol_small <- mol_small
    state
}

#' Prepare Roi Geometry
#' @description
#' Internal helper for `.prepare_roi_geometry`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.prepare_roi_geometry <- function(state) {
    p <- state$params
    ds <- state$datasets$ds
    bounds <- state$bounds$bounds_global
    y_max <- state$bounds$y_max
    verbose <- p$verbose
    flip_y <- p$flip_y
    chunk_pts <- p$chunk_pts
    parent <- "createSCOPE"
    step09 <- .log_step(parent, "S09", "prepare ROI polygon", verbose)

    user_poly <- NULL
    roi_source <- NULL
    roi_bbox <- NULL
    if (!is.null(p$coord_file)) {
        roi_source <- "file"
        step09$enter("source=file")
        if (!file.exists(p$coord_file)) stop("Coordinate file not found: ", p$coord_file)

        ext <- tolower(file_ext(p$coord_file))
        sep_char <- switch(ext,
            csv = ",",
            tsv = "\t",
            txt = "\t",
            "\t"
        )
        poly_tb <- tryCatch(read.table(p$coord_file, header = TRUE, sep = sep_char, stringsAsFactors = FALSE), error = function(e) {
            stop("Error reading coordinate file: ", e$message)
        })
        names(poly_tb) <- tolower(names(poly_tb))
        if (!all(c("x", "y") %in% names(poly_tb))) {
            if (all(c("x_um", "y_um") %in% names(poly_tb))) {
                poly_tb$x <- poly_tb$x_um
                poly_tb$y <- poly_tb$y_um
            } else {
                stop("Coordinate file must contain 'x' and 'y' columns (or 'x_um' and 'y_um').")
            }
        }
        if (nrow(poly_tb) < 3) stop("ROI polygon must have at least 3 points")
        .log_info(parent, "S09", paste0("polygon vertices=", nrow(poly_tb)), verbose)

        poly_geom <- tryCatch({
            geom <- sf::st_polygon(list(as.matrix(poly_tb[, c("x", "y")])))
            sf::st_make_valid(geom)
        }, error = function(e) stop("Error creating polygon geometry: ", e$message))
        if (length(poly_geom) == 0) stop("ROI polygon became EMPTY after validation")

        if (flip_y) {
            poly_tb <- .flip_y_coordinates(poly_tb, y_max)
            poly_geom <- tryCatch(sf::st_make_valid(sf::st_polygon(list(as.matrix(poly_tb[, c("x", "y")])))), error = function(e) {
                stop("Error creating flipped polygon: ", e$message)
            })
        }
        user_poly <- sf::st_sfc(poly_geom)
        roi_bbox <- tryCatch({
            bb <- sf::st_bbox(user_poly)
            c(
                xmin = as.numeric(bb[["xmin"]]),
                xmax = as.numeric(bb[["xmax"]]),
                ymin = as.numeric(bb[["ymin"]]),
                ymax = as.numeric(bb[["ymax"]])
            )
        }, error = function(e) NULL)
    } else {
        roi_source <- "auto_bbox"
        step09$enter("source=auto")
        poly_coords <- data.frame(
            x = c(bounds$xmin, bounds$xmax, bounds$xmax, bounds$xmin, bounds$xmin),
            y = c(bounds$ymin, bounds$ymin, bounds$ymax, bounds$ymax, bounds$ymin)
        )
        if (flip_y) poly_coords <- .flip_y_coordinates(poly_coords, y_max)
        roi_bbox <- c(
            xmin = as.numeric(min(poly_coords$x, na.rm = TRUE)),
            xmax = as.numeric(max(poly_coords$x, na.rm = TRUE)),
            ymin = as.numeric(min(poly_coords$y, na.rm = TRUE)),
            ymax = as.numeric(max(poly_coords$y, na.rm = TRUE))
        )
        poly_geom <- tryCatch({
            geom <- sf::st_polygon(list(as.matrix(poly_coords[, c("x", "y")])))
            sf::st_make_valid(geom)
        }, error = function(e) stop("Error creating automatic ROI polygon: ", e$message))
        if (length(poly_geom) == 0) stop("Automatic ROI polygon became EMPTY after validation")
        user_poly <- sf::st_sfc(poly_geom)
        .log_info(
            parent,
            "S09",
            paste0(
                "auto ROI bounds x=[", round(min(poly_coords$x)), ",", round(max(poly_coords$x)),
                "] y=[", round(min(poly_coords$y)), ",", round(max(poly_coords$y)), "]"
            ),
            verbose
        )
    }

    state$roi <- list(
        geometry = user_poly,
        chunk_pts = chunk_pts,
        source = roi_source,
        bbox = roi_bbox
    )
    step09$done("ok")
    state
}

#' Read Cosmx Expression Subset
#' @description
#' Internal helper for `.read_cosmx_expression_subset`.
#' @param expr_mat Parameter value.
#' @param selected_cols Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param prefix Parameter value.
#' @return Return value used internally.
#' @keywords internal
.read_cosmx_expression_subset <- function(expr_mat, selected_cols, verbose, prefix) {
    dt <- tryCatch(
	        fread(expr_mat, select = selected_cols, showProgress = verbose),
        error = function(e) {
            .log_info("createSCOPE", "S01", paste0("fread failed; falling back to read.csv: ", e$message), verbose)
            read.csv(expr_mat, header = TRUE)[, selected_cols]
        }
    )
	    setDT(dt)
    dt
}

#' Read Cosmx Header
#' @description
#' Internal helper for `.read_cosmx_header`.
#' @param expr_mat Parameter value.
#' @return Return value used internally.
#' @keywords internal
.read_cosmx_header <- function(expr_mat) {
    header_dt <- tryCatch(
	        fread(expr_mat, nrows = 0L, showProgress = FALSE),
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

#' Read Histology Png
#' @description
#' Internal helper for `.read_histology_png`.
#' @param png_path Filesystem path.
#' @param crop_bbox_px Parameter value.
#' @return Return value used internally.
#' @keywords internal
.read_histology_png <- function(png_path,
                                crop_bbox_px = NULL) {
    if (is.null(png_path) || !file.exists(png_path)) {
        stop("Histology image not found: ", png_path)
    }
    ext <- tolower(file_ext(png_path))
    if (ext %in% c("png")) {
        if (!requireNamespace("png", quietly = TRUE)) {
            stop("Package 'png' is required to read histology PNG files.")
        }
        img <- png::readPNG(png_path)
    } else if (ext %in% c("jpg", "jpeg")) {
        if (!requireNamespace("jpeg", quietly = TRUE)) {
            stop("Package 'jpeg' is required to read histology JPEG files.")
        }
        img <- jpeg::readJPEG(png_path)
    } else {
        stop("Unsupported histology image extension: ", ext)
    }
    dims <- dim(img)
    if (length(dims) < 2L) stop("Unexpected PNG dimensions for ", png_path)
    h <- dims[1]
    w <- dims[2]

    adjust_bbox <- function(bbox, w_max, h_max) {
        if (is.null(bbox)) return(NULL)
        req <- c("xmin", "xmax", "ymin", "ymax")
        if (!all(req %in% names(bbox))) {
            stop("crop_bbox_px must be a named vector with xmin/xmax/ymin/ymax.")
        }
        bbox <- as.numeric(bbox[req])
        names(bbox) <- req
        bbox["xmin"] <- max(1, floor(bbox["xmin"]))
        bbox["xmax"] <- min(w_max, ceiling(bbox["xmax"]))
        bbox["ymin"] <- max(1, floor(bbox["ymin"]))
        bbox["ymax"] <- min(h_max, ceiling(bbox["ymax"]))
        if (bbox["xmin"] >= bbox["xmax"] || bbox["ymin"] >= bbox["ymax"]) {
            stop("Invalid crop_bbox_px; xmin/xmax or ymin/ymax collapsed.")
        }
        bbox
    }

    crop_box <- adjust_bbox(crop_bbox_px, w, h)
    if (!is.null(crop_box)) {
        img <- img[crop_box["ymin"]:crop_box["ymax"], crop_box["xmin"]:crop_box["xmax"], , drop = FALSE]
        h <- dim(img)[1]
        w <- dim(img)[2]
    }

    list(
        png = img,
        width = w,
        height = h,
        crop_bbox_px = crop_box,
        source_path = normalizePath(png_path)
    )
}

#' Read Scalefactors Json
#' @description
#' Internal helper for `.read_scalefactors_json`.
#' @param json_path Filesystem path.
#' @param warn_if_missing Parameter value.
#' @return Return value used internally.
#' @keywords internal
.read_scalefactors_json <- function(json_path,
                                    warn_if_missing = TRUE) {
    default <- list(
        hires = NA_real_,
        lowres = NA_real_,
        regist_target_img_scalef = NA_real_,
        microns_per_pixel = NA_real_,
        spot_diameter_fullres = NA_real_
    )
    if (is.null(json_path) || !nzchar(json_path)) {
        if (warn_if_missing) warning("scalefactors_json path is NULL; returning NA scalefactors.")
        return(default)
    }
    if (!file.exists(json_path)) {
        if (warn_if_missing) warning("scalefactors_json not found: ", json_path)
        return(default)
    }
    sf <- tryCatch(jsonlite::fromJSON(json_path),
        error = function(e) {
            stop("Failed to parse scalefactors JSON: ", conditionMessage(e))
        }
    )
    hires_raw <- if (!is.null(sf$tissue_hires_scalef)) {
        sf$tissue_hires_scalef
    } else if (!is.null(sf$tissue_hires_scale)) {
        sf$tissue_hires_scale
    } else {
        sf$hires_scalef
    }
    lowres_raw <- if (!is.null(sf$tissue_lowres_scalef)) {
        sf$tissue_lowres_scalef
    } else if (!is.null(sf$tissue_lowres_scale)) {
        sf$tissue_lowres_scale
    } else {
        sf$lowres_scalef
    }
    regist_raw <- if (!is.null(sf$regist_target_img_scalef)) {
        sf$regist_target_img_scalef
    } else if (!is.null(sf$regist_target_img_scale)) {
        sf$regist_target_img_scale
    } else {
        sf$regist_target_scalef
    }
    list(
        hires = suppressWarnings(as.numeric(hires_raw)),
        lowres = suppressWarnings(as.numeric(lowres_raw)),
        regist_target_img_scalef = suppressWarnings(as.numeric(regist_raw)),
        microns_per_pixel = suppressWarnings(as.numeric(sf$microns_per_pixel)),
        spot_diameter_fullres = suppressWarnings(as.numeric(sf$spot_diameter_fullres))
    )
}

#' Read Xenium Sparse Components
#' @description
#' Internal helper for `.read_xenium_sparse_components`.
#' @param h5_file Filesystem path.
#' @return Return value used internally.
#' @keywords internal
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

#' Resolve Cosmx Gene Panel
#' @description
#' Internal helper for `.resolve_cosmx_gene_panel`.
#' @param filter_genes Parameter value.
#' @param auto_genes Parameter value.
#' @param force_all_genes Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param prefix Parameter value.
#' @return Return value used internally.
#' @keywords internal
.resolve_cosmx_gene_panel <- function(filter_genes,
                                      auto_genes,
                                      force_all_genes,
                                      verbose,
                                      prefix) {
    if (is.null(filter_genes)) {
        if (!isTRUE(force_all_genes) && verbose) {
            .log_info("createSCOPE", "S01", paste0("no filter_genes supplied; using auto panel (", length(auto_genes), " genes)"), verbose)
        }
        gene_cols <- auto_genes
    } else {
        gene_cols <- intersect(auto_genes, filter_genes)
        if (!length(gene_cols)) {
            stop("No overlap between filter_genes and available genes (after auto-detection).")
        }
    }
    .log_info("createSCOPE", "S01", paste0("number of genes selected: ", length(gene_cols)), verbose)
    gene_cols
}

#' Normalize Histology Y Origin
#' @description
#' Internal helper for `.normalize_histology_y_origin`.
#' @param y_origin_choice Parameter value.
#' @return Return value used internally.
#' @keywords internal
.normalize_histology_y_origin <- function(y_origin_choice) {
    choices <- c("auto", "top-left", "bottom-left")
    if (is.null(y_origin_choice) || !length(y_origin_choice)) {
        return("auto")
    }
    y_origin_choice <- as.character(y_origin_choice[1])
    if (is.na(y_origin_choice) || !nzchar(y_origin_choice)) {
        return("auto")
    }
    y_origin_choice <- tolower(y_origin_choice)
    y_origin_choice <- gsub("_", "-", y_origin_choice, fixed = TRUE)
    if (!(y_origin_choice %in% choices)) {
        stop("y_origin must be one of: ", paste(choices, collapse = ", "))
    }
    y_origin_choice
}

#' Resolve Histology Y Origin
#' @description
#' Internal helper for `.resolve_histology_y_origin`.
#' @param grid_layer Layer name.
#' @param y_origin_choice Parameter value.
#' @return Return value used internally.
#' @keywords internal
.resolve_histology_y_origin <- function(grid_layer, y_origin_choice) {
    y_origin_choice <- .normalize_histology_y_origin(y_origin_choice)
    if (!identical(y_origin_choice, "auto")) {
        return(y_origin_choice)
    }
    candidates <- character()
    if (!is.null(grid_layer$image_info$y_origin)) {
        candidates <- c(candidates, grid_layer$image_info$y_origin)
    }
    if (!is.null(grid_layer$histology)) {
        if (length(Filter(function(x) !is.null(x) && !is.null(x$y_origin), grid_layer$histology))) {
            candidates <- c(
                candidates,
                vapply(
                    Filter(function(x) !is.null(x) && !is.null(x$y_origin), grid_layer$histology),
                    `[[`,
                    character(1),
                    "y_origin"
                )
            )
        }
    }
    candidates <- candidates[candidates %in% c("top-left", "bottom-left")]
    if (length(candidates)) {
        return(candidates[1])
    }
    "top-left"
}

#' Resolve Input Paths
#' @description
#' Internal helper for `.resolve_input_paths`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.resolve_input_paths <- function(state) {
    p <- state$params
    parent <- "createSCOPE"
    exp_path <- file.path(p$input_dir, "experiment.xenium")
    exp_meta <- .read_experiment_xenium(p$input_dir)
    if (!is.null(exp_meta) && .log_enabled(p$verbose)) {
        vstr <- .format_xenium_version(exp_meta)
        gen <- .infer_xenium_generation(exp_meta)
        pix <- suppressWarnings(as.numeric(exp_meta$pixel_size_um))
        .log_info(
            parent,
            "S03",
            paste0(
                "Xenium experiment: version=", ifelse(is.na(vstr), "NA", vstr),
                " generation=", gen,
                if (is.finite(pix)) paste0(" pixel_size_um=", format(pix, digits = 6, trim = TRUE)) else ""
            ),
            p$verbose
        )
    }
    tf_parq <- file.path(p$input_dir, "transcripts.parquet")
    tf_csv <- .resolve_transcripts_csv_path(p$input_dir)
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

    .log_info(
        parent,
        "S03",
        paste0(
            "resolved inputs: transcripts=", basename(tf_parq),
            " transcripts_csv=", ifelse(is.null(tf_csv), "none", basename(tf_csv)),
            " cells=", basename(cell_parq),
            " seg_files=", if (length(seg_files)) paste(names(seg_files), collapse = ",") else "none"
        ),
        p$verbose
    )

    state$paths <- list(
        tf_parq = tf_parq,
        tf_csv = tf_csv,
        cell_parq = cell_parq,
        seg_files = seg_files,
        experiment_xenium = if (file.exists(exp_path)) exp_path else NULL
    )
    state$experiment <- exp_meta
    state
}

#' Resolve Xenium Matrix Path
#' @description
#' Internal helper for `.resolve_xenium_matrix_path`.
#' @param xenium_dir Filesystem path.
#' @return Return value used internally.
#' @keywords internal
.resolve_xenium_matrix_path <- function(xenium_dir) {
    h5_file <- file.path(xenium_dir, "cell_feature_matrix.h5")
    if (!file.exists(h5_file)) {
        stop("HDF5 file not found: ", h5_file)
    }
    h5_file
}

#' Execute the CosMx ingestion workflow.
#' @description
#' Internal helper for `.run_cosmx_pipeline`.
#' Loads the expression matrix, aligns identifiers and stores the assembled
#' sparse matrix back into the `scope_object`.
#' @param cfg Configuration list from `.build_singlecell_config()`.
#' @return Return value used internally.
#' @keywords internal
.run_cosmx_pipeline <- function(cfg) {
    parent <- "createSCOPE"
    verbose <- cfg$verbose

    dataset_info <- .locate_cosmx_dataset(cfg$paths$cosmx)
    expr_mat <- dataset_info$expr_mat
    .log_info(parent, "S01", paste0("expression matrix: ", basename(expr_mat)), verbose)

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

    .store_scope_counts(cfg$scope_obj, counts, parent, verbose)
}

#' Run Singlecell Pipeline
#' @description
#' Internal helper for `.run_singlecell_pipeline`.
#' @param cfg Parameter value.
#' @return Return value used internally.
#' @keywords internal
.run_singlecell_pipeline <- function(cfg) {
    if (identical(cfg$platform, "Xenium")) {
        return(.run_xenium_pipeline(cfg))
    }
    if (identical(cfg$platform, "CosMx")) {
        return(.run_cosmx_pipeline(cfg))
    }
    stop("Unsupported platform: ", cfg$platform)
}

#' Run Xenium Pipeline
#' @description
#' Internal helper for `.run_xenium_pipeline`.
#' @param cfg Parameter value.
#' @return Return value used internally.
#' @keywords internal
.run_xenium_pipeline <- function(cfg) {
    parent <- "createSCOPE"
    verbose <- cfg$verbose
    .log_info(parent, "S01", "loading cell-feature matrix from HDF5", verbose)

    h5_file <- .resolve_xenium_matrix_path(cfg$paths$xenium)
    components <- .read_xenium_sparse_components(h5_file)
    matrix_bundle <- .construct_xenium_sparse_matrix(components, parent, verbose)

    filtered <- .filter_gene_panel(
        counts = matrix_bundle$counts,
        filters = cfg$filters,
        verbose = verbose,
        prefix = parent
    )
    aligned <- .align_counts_to_targets(
        counts = filtered$counts,
        targets = cfg$targets,
        filter_cells = cfg$filters$filter_cells,
        verbose = verbose,
        prefix = parent
    )
    if (!is.null(matrix_bundle$gene_map)) {
        attr(aligned$counts, "gene_map") <- matrix_bundle$gene_map
    }
    .store_scope_counts(cfg$scope_obj, aligned$counts, parent, verbose)
}

#' Scan Molecule Batches
#' @description
#' Internal helper for `.scan_molecule_batches`.
#' @param state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.scan_molecule_batches <- function(state) {
    p <- state$params
    ds <- state$datasets$ds
    parent <- "createSCOPE"
    step12b <- .log_step(parent, "S12b", "scan molecules into grid bins", p$verbose)
    step12b$enter()
    scanned_total <- NA_integer_
    if (is.null(p$transcripts_spec) || !is.list(p$transcripts_spec)) {
        state$params$transcripts_spec <- .resolve_transcripts_feature_spec(
            ds = ds,
            input_dir = p$input_dir,
            verbose = p$verbose
        )
        p <- state$params
    }
    spec <- p$transcripts_spec
    qc_genes <- state$datasets$qc_genes
    if (is.null(qc_genes)) qc_genes <- character()
    bounds <- state$bounds$bounds_global
    y_max <- state$bounds$y_max
    ctx <- state$thread_context
    paths <- state$paths
    tf_path <- if (!is.null(paths$tf_parq)) paths$tf_parq else NA_character_
    size_gb <- tryCatch({
        info <- file.info(tf_path)
        if (is.finite(info$size)) info$size / 1024^3 else NA_real_
    }, error = function(e) NA_real_)
    mem_gb <- if (!is.null(ctx) && !is.null(ctx$sys_mem_gb) && is.finite(ctx$sys_mem_gb)) as.numeric(ctx$sys_mem_gb) else NA_real_
    mem_ratio <- if (is.finite(size_gb) && is.finite(mem_gb) && mem_gb > 0) size_gb / mem_gb else NA_real_
    cpu_target <- if (!is.null(ctx) && !is.null(ctx$ncores_safe)) as.integer(ctx$ncores_safe) else NA_integer_
    cpu_target <- ifelse(is.finite(cpu_target) && cpu_target > 0L, cpu_target, as.integer(getOption("arrow.num_threads", 1L)))
    cpu_target <- max(1L, cpu_target)
    cpu_use <- cpu_target

    ncores_arg <- suppressWarnings(as.integer(p$ncores)[1])
    user_requested <- length(ncores_arg) && is.finite(ncores_arg) && ncores_arg > 0L
    auto_throttle <- !user_requested && isTRUE(getOption("geneSCOPE.arrow_auto_throttle", TRUE))
    if (auto_throttle) {
        cap_factor <- NA_real_
        if (is.finite(mem_ratio)) {
            if (mem_ratio >= 0.50) {
                cap_factor <- 0.25
            } else if (mem_ratio >= 0.25) {
                cap_factor <- 0.50
            } else if (mem_ratio >= 0.10) {
                cap_factor <- 0.75
            }
        } else if (is.finite(size_gb)) {
            if (size_gb >= 40) {
                cap_factor <- 0.25
            } else if (size_gb >= 20) {
                cap_factor <- 0.50
            } else if (size_gb >= 5) {
                cap_factor <- 0.75
            }
        }
        if (!is.na(cap_factor)) {
            cpu_use <- suppressWarnings(as.integer(ceiling(cpu_target * cap_factor)))
            if (!is.finite(cpu_use) || cpu_use < 1L) cpu_use <- 1L
            if (cpu_target > 1L) cpu_use <- max(2L, cpu_use)
            cpu_use <- min(cpu_use, cpu_target)
        }
    }
    if (requireNamespace("arrow", quietly = TRUE) && cpu_use < cpu_target) {
        old_opt <- getOption("arrow.num_threads")
        old_cpu <- if (exists("cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
            tryCatch(arrow::cpu_count(), error = function(e) NA_integer_)
        } else {
            NA_integer_
        }
        old_io <- if (exists("io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
            tryCatch(arrow::io_thread_count(), error = function(e) NA_integer_)
        } else {
            NA_integer_
        }
        on.exit({
            try(options(arrow.num_threads = old_opt), silent = TRUE)

            restore_cpu <- suppressWarnings(as.integer(old_cpu))
            if (!is.finite(restore_cpu) || restore_cpu < 1L) restore_cpu <- cpu_target
            restore_cpu <- max(1L, restore_cpu)

            restore_io <- suppressWarnings(as.integer(old_io))
            if (!is.finite(restore_io) || restore_io < 1L) restore_io <- restore_cpu
            restore_io <- max(2L, restore_io)

            if (exists("set_cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
                try(arrow::set_cpu_count(restore_cpu), silent = TRUE)
            }
            if (exists("set_io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
                try(arrow::set_io_thread_count(restore_io), silent = TRUE)
            }
        }, add = TRUE)

        io_use <- max(2L, cpu_use)
        .log_info(
            parent,
            "S12b",
            paste0(
                "throttling Arrow threads for scan: cpu=", cpu_use,
                " io=", io_use,
                " (restore cpu=", cpu_target,
                if (is.finite(size_gb)) paste0(", transcripts_size_gb=", round(size_gb, 1)) else "",
                if (is.finite(mem_gb)) paste0(", sys_mem_gb=", round(mem_gb, 1)) else "",
                ")"
            ),
            p$verbose
        )
        try(options(arrow.num_threads = cpu_use), silent = TRUE)
        if (exists("set_cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
            try(arrow::set_cpu_count(cpu_use), silent = TRUE)
        }
        if (exists("set_io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
            try(arrow::set_io_thread_count(io_use), silent = TRUE)
        }
    }
    counts_list <- setNames(vector("list", length(p$grid_length)), paste0("lg", p$grid_length))
			    for (i in seq_along(counts_list)) {
			        counts_list[[i]] <- data.table(grid_id = character(), gene = character(), count = integer())
			    }
    parent <- "createSCOPE"
    roi_geom <- if (!is.null(state$roi)) state$roi$geometry else NULL
    roi_source <- if (!is.null(state$roi) && !is.null(state$roi$source)) state$roi$source else if (is.null(roi_geom)) "none" else "file"

    compute_counts_from_filtered_dt <- function(dt_in) {
        out <- setNames(vector("list", length(p$grid_length)), paste0("lg", p$grid_length))
        for (i in seq_along(out)) {
            out[[i]] <- data.table::data.table(grid_id = character(), gene = character(), count = integer())
        }
        if (is.null(dt_in)) return(out)

        dt_scan <- data.table::as.data.table(dt_in)
        needed <- c("x_location", "y_location", "feature_name")
        missing <- setdiff(needed, names(dt_scan))
        if (length(missing)) {
            stop("Transcripts table missing required column(s) for scan: ", paste(missing, collapse = ", "))
        }

        dt_scan <- dt_scan[is.finite(x_location) & is.finite(y_location) & !is.na(feature_name) & nzchar(feature_name)]
        if (!nrow(dt_scan)) return(out)

        if (identical(roi_source, "auto_bbox")) {
            bx <- suppressWarnings(as.numeric(bounds$xmin))
            bx2 <- suppressWarnings(as.numeric(bounds$xmax))
            by <- suppressWarnings(as.numeric(bounds$ymin))
            by2 <- suppressWarnings(as.numeric(bounds$ymax))
            if (all(is.finite(c(bx, bx2, by, by2))) && bx < bx2 && by < by2) {
                dt_scan <- dt_scan[
                    x_location > bx & x_location < bx2 &
                        y_location > by & y_location < by2
                ]
                if (!nrow(dt_scan)) return(out)
            }
        }

        chunk_size <- suppressWarnings(as.integer(p$chunk_pts))
        chunk_size <- ifelse(is.finite(chunk_size) && chunk_size > 0L, chunk_size, 500000L)

        x0_lg_all <- vapply(p$grid_length, function(lg) floor(bounds$xmin / lg) * lg, numeric(1))
        y0_lg_all <- vapply(
            p$grid_length,
            function(lg) if (p$flip_y) 0 else floor(bounds$ymin / lg) * lg,
            numeric(1)
        )

        accum_parts <- vector("list", length(p$grid_length))
        for (k in seq_along(accum_parts)) accum_parts[[k]] <- list()

        n_rows <- nrow(dt_scan)
        for (start in seq.int(1L, n_rows, by = chunk_size)) {
            end <- min(n_rows, start + chunk_size - 1L)
            dt <- dt_scan[start:end]
            if (p$flip_y) dt[, y_location := y_max - y_location]

            for (k in seq_along(p$grid_length)) {
                lg <- p$grid_length[k]
                accum <- dt[, .N, by = .(
                    gx = (x_location - x0_lg_all[k]) %/% lg + 1L,
                    gy = (y_location - y0_lg_all[k]) %/% lg + 1L,
                    gene = feature_name
                )]
                accum[, grid_id := paste0("g", gx, "_", gy)]
                accum[, c("gx", "gy") := NULL]
                data.table::setnames(accum, "N", "count")
                data.table::setcolorder(accum, c("grid_id", "gene", "count"))
                accum_parts[[k]][[length(accum_parts[[k]]) + 1L]] <- accum
            }
        }

        for (k in seq_along(out)) {
            parts_k <- accum_parts[[k]]
            if (!length(parts_k)) next
            merged <- data.table::rbindlist(parts_k, use.names = TRUE)
            out[[k]] <- merged[, .(count = sum(count)), by = .(grid_id, gene)]
        }

        out
    }

    do_scan <- is.null(roi_geom) || identical(roi_source, "auto_bbox")
    if (do_scan) {
        cache_dt <- state$datasets$transcripts_cache
        cache_masked <- isTRUE(state$datasets$transcripts_cache_masked)
        cache_src <- if (!is.null(state$datasets$transcripts_cache_source)) as.character(state$datasets$transcripts_cache_source)[1] else NA_character_
        if (!is.null(cache_dt) && isTRUE(cache_masked) && nrow(cache_dt) > 0L) {
            .log_info(
                parent,
                "S12b",
                paste0("using cached transcripts backend source=", ifelse(is.na(cache_src), "unknown", cache_src)),
                p$verbose
            )
            counts_list <- compute_counts_from_filtered_dt(cache_dt)
        } else {
        ds_scan <- ds
        if (identical(roi_source, "auto_bbox")) {
            bx <- suppressWarnings(as.numeric(bounds$xmin))
            bx2 <- suppressWarnings(as.numeric(bounds$xmax))
            by <- suppressWarnings(as.numeric(bounds$ymin))
            by2 <- suppressWarnings(as.numeric(bounds$ymax))
            if (all(is.finite(c(bx, bx2, by, by2))) && bx < bx2 && by < by2) {
                ds_scan <- ds_scan |> filter(
                    x_location > bx, x_location < bx2,
                    y_location > by, y_location < by2
                )
            }
        }
        .log_info(parent, "S12b", "starting multi-resolution scan", p$verbose)
        x0_lg_all <- vapply(p$grid_length, function(lg) floor(bounds$xmin / lg) * lg, numeric(1))
        y0_lg_all <- vapply(
            p$grid_length,
            function(lg) if (p$flip_y) 0 else floor(bounds$ymin / lg) * lg,
            numeric(1)
        )
        process_tbl <- function(tbl) {
            dt <- data.table::setDT(tbl)
            mask <- .ingest_molecule_filter_mask(dt, p, qc_genes)
            fn <- attr(mask, "feature_name")
            if (!is.null(fn) && length(fn) == nrow(dt)) {
                dt[, feature_name := fn]
            }
            if (length(mask)) {
                dt <- dt[mask, .(x_location, y_location, feature_name)]
            } else {
                dt <- dt[0, .(x_location, y_location, feature_name)]
            }
            if (!nrow(dt)) {
                empty <- data.table(grid_id = character(), gene = character(), count = integer())
                return(rep.int(list(empty), length(p$grid_length)))
            }
            if (p$flip_y) dt[, y_location := y_max - y_location]
            out <- vector("list", length(p$grid_length))
            for (k in seq_along(p$grid_length)) {
                lg <- p$grid_length[k]
                accum <- dt[, .N, by = .(
                    gx = (x_location - x0_lg_all[k]) %/% lg + 1L,
                    gy = (y_location - y0_lg_all[k]) %/% lg + 1L,
                    gene = feature_name
                )]
                accum[, grid_id := paste0("g", gx, "_", gy)]
                accum[, c("gx", "gy") := NULL]
                setnames(accum, "N", "count")
                data.table::setcolorder(accum, c("grid_id", "gene", "count"))
                out[[k]] <- accum
            }
            out
        }

        scan_cols <- c("x_location", "y_location", .transcripts_feature_project_cols(spec))
        if (!is.null(p$max_dist_mol_nuc)) scan_cols <- c(scan_cols, "nucleus_distance")
        if (!is.null(p$filterqv)) scan_cols <- c(scan_cols, "qv")
        scan_cols <- unique(scan_cols)

        set_arrow_threads <- function(cpu_threads, io_threads) {
            cpu_threads <- suppressWarnings(as.integer(cpu_threads))
            if (!is.finite(cpu_threads) || cpu_threads < 1L) cpu_threads <- 1L
            io_threads <- suppressWarnings(as.integer(io_threads))
            if (!is.finite(io_threads) || io_threads < 1L) io_threads <- cpu_threads
            io_threads <- max(2L, io_threads)
            try(options(arrow.num_threads = cpu_threads), silent = TRUE)
            if (exists("set_cpu_count", where = asNamespace("arrow"), inherits = FALSE)) {
                try(arrow::set_cpu_count(cpu_threads), silent = TRUE)
            }
            if (exists("set_io_thread_count", where = asNamespace("arrow"), inherits = FALSE)) {
                try(arrow::set_io_thread_count(io_threads), silent = TRUE)
            }
            invisible(list(cpu = cpu_threads, io = io_threads))
        }

        retry_base <- getOption("geneSCOPE.arrow_bounds_retry_threads", c(32L, 16L, 8L, 4L, 2L, 1L))
        retry_opt <- getOption("geneSCOPE.arrow_scan_retry_threads", retry_base)
        retry_vec <- suppressWarnings(as.integer(retry_opt))
        retry_vec <- retry_vec[is.finite(retry_vec) & retry_vec > 0L]
        cpu_start <- suppressWarnings(as.integer(cpu_use))
        if (!is.finite(cpu_start) || cpu_start < 1L) cpu_start <- cpu_target
        cpu_trials <- unique(c(cpu_start, retry_vec))
        cpu_trials <- unique(pmax(1L, pmin(cpu_start, cpu_trials)))
        cpu_trials <- sort(cpu_trials, decreasing = TRUE)
        if (!length(cpu_trials)) cpu_trials <- cpu_start

        scan_once <- function() {
            scanned_local <- 0L
            accum_parts <- vector("list", length(p$grid_length))
            for (k in seq_along(accum_parts)) accum_parts[[k]] <- list()

            reader <- NULL
            reader <- tryCatch({
                if (inherits(ds_scan, "Dataset") && "NewScan" %in% names(ds_scan)) {
                    scan_builder <- ds_scan$NewScan()
                    scan_builder$Project(scan_cols)
                    scan_builder$UseThreads(TRUE)
                    scan_builder$BatchSize(as.integer(p$chunk_pts))
                    scan_builder$Finish()$ToRecordBatchReader()
                } else if (inherits(ds_scan, "Dataset") &&
                    exists("Scanner", where = asNamespace("arrow"), inherits = FALSE) &&
                    "create" %in% names(arrow::Scanner)) {
                    scanner <- arrow::Scanner$create(ds_scan)
                    if ("Project" %in% names(scanner)) scanner$Project(scan_cols)
                    if ("UseThreads" %in% names(scanner)) scanner$UseThreads(TRUE)
                    if ("BatchSize" %in% names(scanner)) scanner$BatchSize(as.integer(p$chunk_pts))
                    if ("ToRecordBatchReader" %in% names(scanner)) scanner$ToRecordBatchReader() else NULL
                } else {
                    arrow::as_record_batch_reader(ds_scan |>
                        dplyr::select(dplyr::all_of(scan_cols)))
                }
            }, error = function(e) NULL)

            if (!is.null(reader)) {
                pending <- list()
                pending_n <- 0L
                repeat {
                    batch <- tryCatch(reader$read_next_batch(), error = function(e) e)
                    if (inherits(batch, "error")) {
                        stop("Arrow scan failed: ", conditionMessage(batch))
                    }
                    if (is.null(batch)) break
                    tbl <- tryCatch(as.data.frame(batch), error = function(e) e)
                    if (inherits(tbl, "error")) {
                        stop("Arrow batch conversion failed: ", conditionMessage(tbl))
                    }
                    # Defensive: some Arrow builds can yield occasional empty batches; skip rather than stop early.
                    # End-of-stream is signaled by `batch` being NULL above.
                    if (is.null(tbl) || nrow(tbl) == 0) next

                    scanned_local <- scanned_local + nrow(tbl)
                    scanned_total <- scanned_local
                    pending[[length(pending) + 1L]] <- tbl
                    pending_n <- pending_n + nrow(tbl)
                    if (pending_n >= p$chunk_pts) {
                        chunk_tbl <- rbindlist(pending, use.names = TRUE)
                        pending <- list()
                        pending_n <- 0L
                        accums <- process_tbl(chunk_tbl)
                        for (k in seq_along(accums)) {
                            accum_parts[[k]][[length(accum_parts[[k]]) + 1L]] <- accums[[k]]
                        }
                    }
                }
                if (length(pending) && pending_n > 0L) {
                    chunk_tbl <- rbindlist(pending, use.names = TRUE)
                    pending <- list()
                    pending_n <- 0L
                    accums <- process_tbl(chunk_tbl)
                    for (k in seq_along(accums)) {
                        accum_parts[[k]][[length(accum_parts[[k]]) + 1L]] <- accums[[k]]
                    }
                }
            } else {
                .log_info(
                    parent,
                    "S12b",
                    "warning: streaming Arrow scan unavailable; using slow slice/collect fallback (consider upgrading Arrow)",
                    p$verbose
                )
                repeat {
                    tbl <- tryCatch({
                        if (scanned_local == 0) {
                            ds_scan |>
                                dplyr::select(dplyr::all_of(scan_cols)) |>
                                slice_head(n = p$chunk_pts) |>
                                collect()
                        } else {
                            batch <- ds_scan |>
                                dplyr::select(dplyr::all_of(scan_cols)) |>
                                slice(scanned_local + 1, scanned_local + p$chunk_pts) |>
                                collect()
                            if (nrow(batch) == 0) NULL else batch
                        }
                    }, error = function(e) e)
                    if (inherits(tbl, "error")) {
                        stop("Arrow scan collect failed: ", conditionMessage(tbl))
                    }
                    if (is.null(tbl) || nrow(tbl) == 0) break

                    scanned_local <- scanned_local + nrow(tbl)
                    scanned_total <- scanned_local
                    accums <- process_tbl(tbl)
                    for (k in seq_along(accums)) {
                        accum_parts[[k]][[length(accum_parts[[k]]) + 1L]] <- accums[[k]]
                    }
                }
            }

            list(accum_parts = accum_parts, scanned = scanned_local)
        }

	        ds_queryable <- !isFALSE(state$datasets$ds_queryable)
	        probe_err <- if (!is.null(state$datasets$ds_query_probe_error) && nzchar(state$datasets$ds_query_probe_error)) {
	            state$datasets$ds_query_probe_error
	        } else {
	            "Arrow dataset query probe failed"
	        }
	        scan_out <- NULL
	        last_err <- NULL
	        if (isTRUE(ds_queryable)) {
	            for (i in seq_along(cpu_trials)) {
	                cpu_i <- cpu_trials[[i]]
	                io_i <- max(2L, cpu_i)
	                if (i > 1L) {
	                    .log_info(parent, "S12b", paste0("scan retry: cpu=", cpu_i, " io=", io_i), p$verbose)
	                }
	                set_arrow_threads(cpu_i, io_i)
	                res_i <- tryCatch(scan_once(), error = function(e) e)
	                if (!inherits(res_i, "error")) {
	                    scan_out <- res_i
	                    if (i > 1L) {
	                        .log_info(parent, "S12b", paste0("scan recovered: cpu=", cpu_i, " io=", io_i), p$verbose)
	                    }
	                    break
	                }
	                if (i == 1L && length(cpu_trials) > 1L) {
	                    .log_info(parent, "S12b", paste0("scan failed; retrying: ", conditionMessage(res_i)), p$verbose)
	                }
	                last_err <- res_i
	            }
	        } else {
	            last_err <- simpleError(paste0("Arrow dataset query unavailable: ", probe_err))
	        }
	        if (is.null(scan_out)) {
	            can_try_compat <- is.character(paths$tf_parq) &&
	                !is.na(paths$tf_parq) &&
	                nzchar(paths$tf_parq) &&
	                !grepl("transcripts\\.arrow_compat\\.parquet$", basename(paths$tf_parq), ignore.case = TRUE)
	            if (isTRUE(can_try_compat)) {
	                .log_info(parent, "S12b", "attempting Arrow-compatible transcript rewrite via pyarrow before direct fallback scan", p$verbose)
	                compat_path <- tryCatch(
	                    .ensure_transcripts_arrow_compat_parquet(paths$tf_parq, input_dir = p$input_dir, verbose = p$verbose),
	                    error = function(e) e
	                )
	                if (!inherits(compat_path, "error")) {
	                    compat_ds <- tryCatch(arrow::open_dataset(compat_path), error = function(e) e)
	                    if (!inherits(compat_ds, "error")) {
	                        .log_info(parent, "S12b", paste0("retrying Arrow scan on Arrow-compatible transcripts: ", basename(compat_path)), p$verbose)
	                        state$paths$tf_parq_original <- paths$tf_parq
	                        state$paths$tf_parq <- compat_path
	                        state$datasets$ds <- compat_ds
	                        state$datasets$ds_raw <- compat_ds
	                        state$datasets$ds_queryable <- TRUE
	                        state$datasets$ds_query_probe_error <- NULL
	                        paths <- state$paths
	                        ds_scan <- compat_ds
	                        last_err <- NULL
	                        for (i in seq_along(cpu_trials)) {
	                            cpu_i <- cpu_trials[[i]]
	                            io_i <- max(2L, cpu_i)
	                            if (i > 1L) {
	                                .log_info(parent, "S12b", paste0("compat scan retry: cpu=", cpu_i, " io=", io_i), p$verbose)
	                            }
	                            set_arrow_threads(cpu_i, io_i)
	                            res_i <- tryCatch(scan_once(), error = function(e) e)
	                            if (!inherits(res_i, "error")) {
	                                scan_out <- res_i
	                                break
	                            }
	                            last_err <- res_i
	                        }
	                    } else {
	                        .log_info(parent, "S12b", paste0("Arrow-compatible transcripts open failed: ", conditionMessage(compat_ds)), p$verbose)
	                    }
	                } else {
	                    .log_info(parent, "S12b", paste0("Arrow-compatible transcript rewrite failed: ", conditionMessage(compat_path)), p$verbose)
	                }
	            }
	        }
	        if (is.null(scan_out)) {
	            .log_info(
	                parent,
	                "S12b",
	                paste0("Arrow scan failed; falling back to direct transcripts read: ", conditionMessage(last_err)),
                p$verbose
            )
            try(set_arrow_threads(cpu_target, max(2L, cpu_target)), silent = TRUE)

            .log_info(parent, "S12b", "starting fallback transcripts load for grid scan", p$verbose)
            fb <- tryCatch(
                .read_transcripts_fallback_dt(
                    paths,
                    p$input_dir,
                    scan_cols,
                    verbose = p$verbose,
                    use_threads = TRUE,
                    log_parent = parent,
                    log_step = "S12b"
                ),
                error = function(e) e
            )
            if (inherits(fb, "error")) {
                .log_info(parent, "S13", paste0("fallback transcripts read failed; retrying single-thread: ", conditionMessage(fb)), p$verbose)
                fb <- .read_transcripts_fallback_dt(
                    paths,
                    p$input_dir,
                    scan_cols,
                    verbose = p$verbose,
                    use_threads = FALSE,
                    log_parent = parent,
                    log_step = "S12b"
                )
            }
            dt_fb <- fb$dt
            mask_fb <- .ingest_molecule_filter_mask(dt_fb, p, qc_genes)
            fn_fb <- attr(mask_fb, "feature_name")
            if (!is.null(fn_fb) && length(fn_fb) == nrow(dt_fb)) {
                dt_fb[, feature_name := fn_fb]
            }
            if (!length(mask_fb) || !any(mask_fb)) {
                stop("No molecules available after filters; cannot scan molecules.")
            }
            dt_filtered <- dt_fb[mask_fb, .(x_location, y_location, feature_name)]
            dt_filtered <- dt_filtered[is.finite(x_location) & is.finite(y_location) & !is.na(feature_name) & nzchar(feature_name)]
            if (!nrow(dt_filtered)) {
                stop("No molecules available after filters; cannot scan molecules.")
            }

            state$datasets$transcripts_cache <- dt_filtered
            state$datasets$transcripts_cache_masked <- TRUE
            state$datasets$transcripts_cache_source <- fb$source
            state$datasets$transcripts_cache_path <- fb$path

            counts_list <- compute_counts_from_filtered_dt(dt_filtered)
        } else {
            accum_parts <- scan_out$accum_parts
            for (k in seq_along(counts_list)) {
                parts_k <- accum_parts[[k]]
                if (!length(parts_k)) next
                merged <- rbindlist(parts_k, use.names = TRUE)
                counts_list[[k]] <- merged[, .(count = sum(count)), by = .(grid_id, gene)]
            }
        }
    }
    }
    state$objects$counts_list <- counts_list
    step12b$done(if (!is.na(scanned_total)) paste0("scanned=", format(scanned_total, big.mark = ",")) else "ok")
    state
}

#' Store Scope Counts
#' @description
#' Internal helper for `.store_scope_counts`.
#' @param scope_obj A `scope_object`.
#' @param counts Parameter value.
#' @param prefix Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.store_scope_counts <- function(scope_obj, counts, prefix, verbose) {
    scope_obj@cells <- list(counts = counts)
    .log_info("createSCOPE", "S01", "cell matrix integration completed", verbose)
    invisible(scope_obj)
}

#' Build Scope From Visium Seurat
#' @description
#' Internal helper for `.build_scope_from_visium_seurat`.
#' @param input_dir Filesystem path.
#' @param run_sctransform Parameter value.
#' @param sct_return_only_var_genes Parameter value.
#' @param vars_to_regress Parameter value.
#' @param glmGamPoi Parameter value.
#' @param include_in_tissue_spots Parameter value.
#' @param grid_multiplier Parameter value.
#' @param flip_y Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param data_type Parameter value.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.build_scope_from_visium_seurat <- function(input_dir,
                                          run_sctransform = TRUE,
                                          sct_return_only_var_genes = FALSE,
                                          vars_to_regress = NULL,
                                          glmGamPoi = TRUE,
                                          include_in_tissue_spots = TRUE,
                                          grid_multiplier = 1,
                                          flip_y = FALSE,
                                          verbose = TRUE,
                                          data_type = "visium") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' is required for .build_scope_from_visium_seurat(). Please install it.")
    }
    if (file.exists(file.path(input_dir, "filtered_feature_bc_matrix.h5")) &&
        !requireNamespace("hdf5r", quietly = TRUE)) {
        stop(
            "Visium HDF5 ingest requires package 'hdf5r' before calling ",
            "Seurat::Load10X_Spatial(). Install hdf5r in the active R library; ",
            "this is an ingest dependency issue, not a Lee's L or clustering failure.",
            call. = FALSE
        )
    }
    .with_data_table_attached(function() {
        inputs_validated <- .visium_scope_validate_inputs(
            input_dir = input_dir,
            grid_multiplier = grid_multiplier,
            include_in_tissue_spots = include_in_tissue_spots,
            flip_y = flip_y,
            data_type = data_type
        )
        spec <- .visium_scope_spec_build(inputs_validated)
        payload <- .visium_scope_payload_build(spec)
        include_in_tissue_spots <- spec$include_in_tissue_spots
        flip_y <- spec$flip_y
        grid_multiplier <- spec$grid_multiplier
        data_type <- spec$data_type

        parent <- "createSCOPE"

    .log_info(parent, "S04", "calling Seurat::Load10X_Spatial()", verbose)
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
        .log_info(parent, "S04", "calling Seurat::SCTransform()", verbose)
        seu <- Seurat::SCTransform(
            object = seu,
            assay = "Spatial",
            new.assay.name = "SCT",
            return.only.var.genes = sct_return_only_var_genes,
            vars.to.regress = vars_to_regress,
            method = if (isTRUE(glmGamPoi) && requireNamespace("glmGamPoi", quietly = TRUE)) "glmGamPoi" else "poisson",
            verbose = isTRUE(verbose)
        )
        Seurat::DefaultAssay(seu) <- "SCT"
    }

    getAssayLayer <- function(obj, assay, which) {
        ga <- if (requireNamespace("SeuratObject", quietly = TRUE)) {
            getExportedValue("SeuratObject", "GetAssayData")
        } else {
            Seurat::GetAssayData
        }
        res <- tryCatch(
            ga(object = obj, assay = assay, layer = which),
            error = function(e) e
        )
        if (inherits(res, "error")) {
            msg <- conditionMessage(res)
            if (grepl("unused argument.*layer", msg)) {
                return(ga(object = obj, assay = assay, slot = which))
            }
            stop(msg)
        }
        res
    }

    raw_counts <- getAssayLayer(seu, assay = "Spatial", which = "counts")
    if (!inherits(raw_counts, "dgCMatrix")) raw_counts <- as(raw_counts, "dgCMatrix")
    sct_mat <- NULL
    if (isTRUE(run_sctransform)) {
        sct_mat <- getAssayLayer(seu, assay = "SCT", which = "scale.data")
    }

    spatial_dir <- spec$spatial_dir
    pos_path <- payload$pos_path

	    pos_raw <- fread(pos_path)
    standardize_positions <- function(dt) {
	        if (!is.data.table(dt)) dt <- as.data.table(dt)
        if (ncol(dt) < 6) stop("Visium position file must contain at least six columns.")
        if (all(grepl("^V[0-9]+$", names(dt)))) {
            cols <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
	            setnames(dt, cols[seq_len(ncol(dt))])
        } else {
            lower <- tolower(names(dt))
            rename_first <- function(target, aliases) {
                idx <- which(lower %in% aliases)
                if (length(idx)) {
	                    setnames(dt, idx[1], target)
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
        if (length(setdiff(required, names(dt)))) {
            stop("Position file missing required columns: ", paste(setdiff(required, names(dt)), collapse = ", "))
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
    if (isTRUE(include_in_tissue_spots)) pos_dt <- pos_dt[in_tissue == 1L]
    keep <- pos_dt$barcode %in% colnames(raw_counts)
    if (!all(keep)) pos_dt <- pos_dt[keep]
    if (nrow(pos_dt) == 0L) stop("No overlap between positions and Seurat barcodes.")

    sf_path <- payload$scalefactor_path
    if (is.na(sf_path)) stop("scalefactors_json.json not found under ", spatial_dir)
    sf <- jsonlite::fromJSON(sf_path)
    microns_per_pixel <- sf$microns_per_pixel
    spot_diameter_um <- sf$spot_diameter_microns
    spot_diameter_px <- sf$spot_diameter_fullres
    hires_scalef <- sf$tissue_hires_scalef %||% sf$tissue_hires_scale %||% NA_real_
    lowres_scalef <- sf$tissue_lowres_scalef %||% sf$tissue_lowres_scale %||% NA_real_
    regist_target_img_scalef <- sf$regist_target_img_scalef %||% sf$regist_target_img_scale %||% NA_real_
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
        .log_info(
            parent,
            "S04",
            paste0(
                "microns_per_pixel missing; assumed spot_diameter_microns=",
                default_spot_um,
                ", computed microns_per_pixel=",
                format(round(microns_per_pixel, 9), scientific = FALSE)
            ),
            verbose
        )
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
    y_flip_reference_um <- NA_real_
    if (isTRUE(flip_y)) {
        y_max <- max(pos_dt$y_um, na.rm = TRUE)
        y_flip_reference_um <- y_max
        pos_dt[, y_um := y_max - y_um]
    }

    feature_file <- payload$feature_file
    if (!is.na(feature_file) && file.exists(feature_file)) {
	        features_dt <- fread(feature_file, header = FALSE)
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

	    centroids_dt <- data.table(
	        cell = pos_dt$barcode,
	        x = pos_dt$x_um,
	        y = pos_dt$y_um,
	        array_row = pos_dt$array_row,
	        array_col = pos_dt$array_col,
	        in_tissue = pos_dt$in_tissue
	    )

    if (anyNA(pos_dt$array_row) || anyNA(pos_dt$array_col)) {
        stop("Visium positions missing array_row/array_col values; cannot build grid indices.")
    }
    col_min <- min(pos_dt$array_col, na.rm = TRUE)
    row_min <- min(pos_dt$array_row, na.rm = TRUE)
    if (!is.finite(col_min) || !is.finite(row_min)) {
        stop("Visium positions contain non-finite array_row/array_col values.")
    }
    col_shift <- if (col_min < 1L) 1L - col_min else 0L
    row_shift <- if (row_min < 1L) 1L - row_min else 0L
    gx_vals <- as.integer(pos_dt$array_col + col_shift)
    gy_vals <- as.integer(pos_dt$array_row + row_shift)

	    grid_dt <- data.table(
	        grid_id = pos_dt$barcode,
	        gx = gx_vals,
	        gy = gy_vals,
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
    xbins_eff <- max(grid_dt$gx, na.rm = TRUE)
    ybins_eff <- max(grid_dt$gy, na.rm = TRUE)

    col_idx <- match(pos_dt$barcode, colnames(raw_counts))
    raw_counts <- raw_counts[, col_idx, drop = FALSE]
    colnames(raw_counts) <- pos_dt$barcode

    if (!is.null(sct_mat)) {
        col_idx2 <- match(pos_dt$barcode, colnames(sct_mat))
        sct_mat <- sct_mat[, col_idx2, drop = FALSE]
        colnames(sct_mat) <- pos_dt$barcode
    }

    counts_tbl <- summary(raw_counts)
	    counts_dt <- data.table(
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
	        dt <- as.data.table(cent)[, .(x, y, array_row = as.integer(array_row), array_col = as.integer(array_col))]
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
	        edges <- rbindlist(list(
	            get_pairs(0, 1), get_pairs(1, 0), get_pairs(1, 1), get_pairs(1, -1)
	        ), use.names = TRUE, fill = TRUE)
        if (is.null(edges) || !nrow(edges)) return(list(orientation = NA_character_, horiz_frac = NA_real_, vert_frac = NA_real_, tol_deg = tol_deg))
        theta <- abs(atan2(edges$dy, edges$dx) * 180 / pi)
        hfrac <- mean(pmin(theta, 180 - theta) <= tol_deg, na.rm = TRUE)
        vfrac <- mean(abs(theta - 90) <= tol_deg, na.rm = TRUE)
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
        ext <- tolower(file_ext(p))
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
        regist_target_img_scalef = regist_target_img_scalef,
        spot_diameter_fullres = spot_diameter_px,
        y_origin = if (isTRUE(flip_y)) "bottom-left" else "top-left",
        flip_y_applied = isTRUE(flip_y),
        y_flip_reference_um = if (isTRUE(flip_y)) y_flip_reference_um else NA_real_
    )

    sf_lookup <- list(
        hires = suppressWarnings(as.numeric(hires_scalef)),
        lowres = suppressWarnings(as.numeric(lowres_scalef)),
        regist_target_img_scalef = suppressWarnings(as.numeric(regist_target_img_scalef))
    )
    y_origin_img <- scope_obj@grid[[grid_name]]$image_info$y_origin
    if (is.null(y_origin_img)) {
        y_origin_img <- if (isTRUE(flip_y)) "bottom-left" else "top-left"
    }
    safe_.attach_histology <- function(scope_obj,
                                      level,
                                      path) {
        if (is.null(path) || is.na(path)) return(scope_obj)
        tryCatch(
            .attach_histology(
                scope_obj = scope_obj,
                grid_name = grid_name,
                png_path = path,
                json_path = sf_path,
                level = level,
                roi_bbox = roi_bbox,
                scalefactors = sf_lookup,
                y_origin = y_origin_img,
                y_flip_reference_um = if (isTRUE(flip_y)) y_flip_reference_um else NULL
            ),
            error = function(e) {
                .log_info(parent, "S14", paste0("histology attach (", level, ") skipped: ", conditionMessage(e)), verbose)
                scope_obj
            }
        )
    }
    scope_obj <- safe_.attach_histology(scope_obj, "hires", if (is.na(hires_path)) NULL else hires_path)
    scope_obj <- safe_.attach_histology(scope_obj, "lowres", if (is.na(lowres_path)) NULL else lowres_path)

        if (!is.null(sct_mat)) {
            scope_obj@grid[[grid_name]]$SCT <- t(sct_mat)
        }

        scope_obj@stats$platform <- data_type
        scope_obj@meta.data$platform <- data_type
        .log_info(
            parent,
            "S14",
            paste0("scope_object created with ", nrow(centroids_dt), " spots and grid_length=", round(grid_size_um, 3), " um"),
            verbose
        )

        scope_obj
    }, context = "[geneSCOPE::.create_scope] Visium/Seurat builder")
}
#' Visium Validate Inputs
#' @description
#' Internal helper for `.visium_validate_inputs`.
#' @param input_dir Filesystem path.
#' @param grid_multiplier Parameter value.
#' @param data_type Parameter value.
#' @return Return value used internally.
#' @keywords internal
.visium_validate_inputs <- function(input_dir, grid_multiplier, data_type) {
    stopifnot(dir.exists(input_dir))
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' is required for .build_scope_from_visium_seurat(). Please install it.")
    }
    if (!is.numeric(grid_multiplier) || length(grid_multiplier) != 1L || !is.finite(grid_multiplier) || grid_multiplier <= 0) {
        stop("grid_multiplier must be a single positive finite numeric value.")
    }
    sprintf("[geneSCOPE::build_scope/%s]", data_type)
}

#' Visium Load Seurat Data
#' @description
#' Internal helper for `.visium_load_seurat_data`.
#' @param input_dir Filesystem path.
#' @param run_sctransform Parameter value.
#' @param sct_return_only_var_genes Parameter value.
#' @param vars_to_regress Parameter value.
#' @param glmGamPoi Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param data_type Parameter value.
#' @return Return value used internally.
#' @keywords internal
.visium_load_seurat_data <- function(input_dir,
                                     run_sctransform,
                                     sct_return_only_var_genes,
                                     vars_to_regress,
                                     glmGamPoi,
                                     verbose,
                                     data_type) {
    if (file.exists(file.path(input_dir, "filtered_feature_bc_matrix.h5")) &&
        !requireNamespace("hdf5r", quietly = TRUE)) {
        stop(
            "Visium HDF5 ingest requires package 'hdf5r' before calling ",
            "Seurat::Load10X_Spatial(). Install hdf5r in the active R library; ",
            "this is an ingest dependency issue, not a Lee's L or clustering failure.",
            call. = FALSE
        )
    }
    .with_data_table_attached(function() {
        .log_info("createSCOPE", "S04", "calling Seurat::Load10X_Spatial()", verbose)
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
            .log_info("createSCOPE", "S04", "calling Seurat::SCTransform()", verbose)
            seu <- Seurat::SCTransform(
                object = seu,
                assay = "Spatial",
                new.assay.name = "SCT",
                return.only.var.genes = sct_return_only_var_genes,
                vars.to.regress = vars_to_regress,
                method = if (isTRUE(glmGamPoi) && requireNamespace("glmGamPoi", quietly = TRUE)) "glmGamPoi" else "poisson",
                verbose = isTRUE(verbose)
            )
            Seurat::DefaultAssay(seu) <- "SCT"
        }

            getAssayLayer <- function(obj, assay, which) {
                ga <- if (requireNamespace("SeuratObject", quietly = TRUE)) {
                    getExportedValue("SeuratObject", "GetAssayData")
                } else {
                    Seurat::GetAssayData
                }
                res <- tryCatch(
                    ga(object = obj, assay = assay, layer = which),
                    error = function(e) e
                )
                if (inherits(res, "error")) {
                    msg <- conditionMessage(res)
                    if (grepl("unused argument.*layer", msg)) {
                        return(ga(object = obj, assay = assay, slot = which))
                    }
                    stop(msg)
                }
                res
            }

        raw_counts <- getAssayLayer(seu, assay = "Spatial", which = "counts")
        if (!inherits(raw_counts, "dgCMatrix")) raw_counts <- as(raw_counts, "dgCMatrix")
        sct_mat <- NULL
        if (isTRUE(run_sctransform)) {
            sct_mat <- getAssayLayer(seu, assay = "SCT", which = "scale.data")
        }

        list(raw_counts = raw_counts, sct_mat = sct_mat, data_type = data_type)
    }, context = "[geneSCOPE::.create_scope] Visium/Seurat load")
}

#' Visium Load Scalefactors
#' @description
#' Internal helper for `.visium_load_scalefactors`.
#' @param spatial_dir Filesystem path.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.visium_load_scalefactors <- function(spatial_dir, verbose) {
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
    regist_target_img_scalef <- sf$regist_target_img_scalef %||% sf$regist_target_img_scale %||% NA_real_
    list(
        microns_per_pixel = microns_per_pixel,
        spot_diameter_um = spot_diameter_um,
        spot_diameter_px = spot_diameter_px,
        hires_scalef = hires_scalef,
        lowres_scalef = lowres_scalef,
        regist_target_img_scalef = regist_target_img_scalef,
        sf = sf
    )
}

#' Visium Resolve Scaling
#' @description
#' Internal helper for `.visium_resolve_scaling`.
#' @param sf_info Parameter value.
#' @param grid_multiplier Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.visium_resolve_scaling <- function(sf_info, grid_multiplier, verbose) {
    microns_per_pixel <- sf_info$microns_per_pixel
    spot_diameter_um <- sf_info$spot_diameter_um
    spot_diameter_px <- sf_info$spot_diameter_px
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
        .log_info(
            "createSCOPE",
            "S04",
            paste0(
                "microns_per_pixel missing; assumed spot_diameter_microns=",
                default_spot_um,
                ", computed microns_per_pixel=",
                format(round(microns_per_pixel, 9), scientific = FALSE)
            ),
            verbose
        )
    }
    if (is.null(microns_per_pixel)) stop("scalefactors missing microns_per_pixel; cannot convert coordinates.")
    if (is.null(spot_diameter_um) && !is.null(spot_diameter_px)) {
        spot_diameter_um <- as.numeric(spot_diameter_px) * as.numeric(microns_per_pixel)
    }
    if (is.null(spot_diameter_um)) stop("scalefactors missing spot diameter information.")

    microns_per_pixel <- as.numeric(microns_per_pixel)
    base_spot_um <- as.numeric(spot_diameter_um)
    grid_size_um <- base_spot_um * as.numeric(grid_multiplier)

    list(
        microns_per_pixel = microns_per_pixel,
        base_spot_um = base_spot_um,
        grid_size_um = grid_size_um
    )
}

#' Visium Load Positions
#' @description
#' Internal helper for `.visium_load_positions`.
#' @param spatial_dir Filesystem path.
#' @param include_in_tissue_spots Parameter value.
#' @param raw_counts Parameter value.
#' @param microns_per_pixel Parameter value.
#' @param flip_y Parameter value.
#' @return Return value used internally.
#' @keywords internal
.visium_load_positions <- function(spatial_dir, include_in_tissue_spots, raw_counts, microns_per_pixel, flip_y) {
    pos_candidates <- c(
        file.path(spatial_dir, "tissue_positions_list.csv"),
        file.path(spatial_dir, "tissue_positions.csv")
    )
    pos_path <- pos_candidates[file.exists(pos_candidates)][1]
    if (is.na(pos_path)) stop("Could not find tissue_positions_list.csv or tissue_positions.csv under ", spatial_dir)

	    pos_raw <- fread(pos_path)
    standardize_positions <- function(dt) {
	        if (!is.data.table(dt)) dt <- as.data.table(dt)
        if (ncol(dt) < 6) stop("Visium position file must contain at least six columns.")
        if (all(grepl("^V[0-9]+$", names(dt)))) {
            cols <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
	            setnames(dt, cols[seq_len(ncol(dt))])
        } else {
            lower <- tolower(names(dt))
            rename_first <- function(target, aliases) {
                idx <- which(lower %in% aliases)
                if (length(idx)) {
	                    setnames(dt, idx[1], target)
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
        if (length(setdiff(required, names(dt)))) {
            stop("Position file missing required columns: ", paste(setdiff(required, names(dt)), collapse = ", "))
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
    if (isTRUE(include_in_tissue_spots)) pos_dt <- pos_dt[in_tissue == 1L]
    keep <- pos_dt$barcode %in% colnames(raw_counts)
    if (!all(keep)) pos_dt <- pos_dt[keep]
    if (nrow(pos_dt) == 0L) stop("No overlap between positions and Seurat barcodes.")

    pos_dt[, `:=`(
        x_um = pxl_col_in_fullres * microns_per_pixel,
        y_um = pxl_row_in_fullres * microns_per_pixel
    )]
    y_flip_reference_um <- NA_real_
    if (isTRUE(flip_y)) {
        y_max <- max(pos_dt$y_um, na.rm = TRUE)
        y_flip_reference_um <- y_max
        pos_dt[, y_um := y_max - y_um]
    }
    list(pos_dt = pos_dt, y_flip_reference_um = y_flip_reference_um)
}

#' Visium Build Grid Components
#' @description
#' Internal helper for `.visium_build_grid_components`.
#' @param pos_dt Parameter value.
#' @param grid_multiplier Parameter value.
#' @param microns_per_pixel Parameter value.
#' @param raw_counts Parameter value.
#' @param sct_mat Parameter value.
#' @param data_type Parameter value.
#' @return Return value used internally.
#' @keywords internal
.visium_build_grid_components <- function(pos_dt, grid_multiplier, microns_per_pixel, raw_counts, sct_mat, data_type) {
    base_spot_um <- max(1, as.numeric(grid_multiplier) * 0 + 1) # placeholder to satisfy R CMD check for global vars
    grid_size_um <- unique(pos_dt$x_um * 0 + grid_multiplier * 0 + 1) # overwritten by caller
	    setDT(pos_dt)

    col_map <- setNames(
        seq_along(sort(unique(pos_dt$array_col))),
        sort(unique(pos_dt$array_col))
    )
    row_map <- setNames(
        seq_along(sort(unique(pos_dt$array_row))),
        sort(unique(pos_dt$array_row))
    )

    list(
        pos_dt = pos_dt,
        col_map = col_map,
        row_map = row_map,
        raw_counts = raw_counts,
        sct_mat = sct_mat,
        data_type = data_type,
        microns_per_pixel = microns_per_pixel,
        base_spot_um = base_spot_um,
        grid_size_um = grid_size_um
    )
}

#' Internal helper for correlation workflows
#' @description
#' Internal helper for `.compute_correlation_assemble_outputs`.
#' @param scope_obj Internal parameter
#' @param level Parameter value.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param store_layer Layer name.
#' @param cor_matrix Parameter value.
#' @param fdr_matrix Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.compute_correlation_assemble_outputs <- function(scope_obj,
                                                level,
                                                grid_name,
                                                store_layer,
                                                cor_matrix,
                                                fdr_matrix,
                                                verbose) {
  if (is.null(scope_obj@stats) || !is.list(scope_obj@stats)) {
    scope_obj@stats <- list()
  }
  if (level == "cell") {
    if (is.null(scope_obj@stats[["cell"]]) || !is.list(scope_obj@stats[["cell"]])) {
      scope_obj@stats[["cell"]] <- list()
    }
    scope_obj@stats[["cell"]][[store_layer]] <- cor_matrix
  } else {
    if (is.null(scope_obj@stats[[grid_name]]) || !is.list(scope_obj@stats[[grid_name]])) {
      scope_obj@stats[[grid_name]] <- list()
    }
    scope_obj@stats[[grid_name]][[store_layer]] <- cor_matrix
  }

  if (level == "cell") {
    scope_obj@cells[[store_layer]] <- cor_matrix
  } else {
    scope_obj@grid[[grid_name]][[store_layer]] <- cor_matrix
  }

  if (!is.null(fdr_matrix)) {
    fdr_layer_name <- paste0(store_layer, "_FDR")
    if (level == "cell") {
      scope_obj@stats[["cell"]][[fdr_layer_name]] <- fdr_matrix
      scope_obj@cells[[fdr_layer_name]] <- fdr_matrix
      .log_info("computeCorrelation", "S06",
        paste0("FDR matrix stored in @cells$", fdr_layer_name), verbose
      )
    } else {
      scope_obj@stats[[grid_name]][[fdr_layer_name]] <- fdr_matrix
      scope_obj@grid[[grid_name]][[fdr_layer_name]] <- fdr_matrix
      .log_info("computeCorrelation", "S06",
        paste0("FDR matrix stored in @grid[['", grid_name, "']]$", fdr_layer_name), verbose
      )
    }
  }

  if (level == "cell") {
    .log_info("computeCorrelation", "S06",
      paste0("Correlation matrix stored in @cells$", store_layer), verbose
    )
    .log_info("computeCorrelation", "S06",
      paste0("Correlation matrix stored in @stats[['cell']]$", store_layer), verbose
    )
  } else {
    .log_info("computeCorrelation", "S06",
      paste0("Correlation matrix stored in @grid[['", grid_name, "']]$", store_layer), verbose
    )
    .log_info("computeCorrelation", "S06",
      paste0("Correlation matrix stored in @stats[['", grid_name, "']]$", store_layer), verbose
    )
  }

  scope_obj
}

#' Internal helper for correlation workflows
#' @description
#' Internal helper for `.compute_correlation_validate_inputs`.
#' @param scope_obj Internal parameter
#' @param level Parameter value.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param layer Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @return Internal helper result
#' @keywords internal
.compute_correlation_validate_inputs <- function(scope_obj,
                                               level,
                                               grid_name,
                                               layer,
                                               verbose) {
  .log_info("computeCorrelation", "S01",
    paste0("Loading expression matrix from level: ", level, " layer=", layer), verbose
  )
  if (level == "cell") {
    if (is.null(scope_obj@cells) || is.null(scope_obj@cells[[layer]])) {
      stop("Cell layer '", layer, "' not found in @cells")
    }
    expr_mat <- scope_obj@cells[[layer]]
  } else {
    grid_name <- if (is.null(grid_name)) {
      names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), .select_grid_layer(scope_obj, grid_name))]
    } else grid_name
    .log_info("computeCorrelation", "S01",
      paste0("Resolved grid_name=", grid_name, " layer=", layer), verbose
    )
    if (is.null(.select_grid_layer(scope_obj, grid_name)[[layer]])) {
      stop("Grid layer '", layer, "' not found in @grid[['", grid_name, "']]")
    }
    expr_mat <- t(
      if (!inherits(.select_grid_layer(scope_obj, grid_name)[[layer]], "dgCMatrix")) {
        as(.select_grid_layer(scope_obj, grid_name)[[layer]], "dgCMatrix")
      } else {
        .select_grid_layer(scope_obj, grid_name)[[layer]]
      }
    )
  }

  if (!inherits(expr_mat, "dgCMatrix")) {
    .log_info("computeCorrelation", "S01", "Converting to sparse matrix...", verbose)
    expr_mat <- as(expr_mat, "dgCMatrix")
  }
  .validate_core_numeric_input(
    expr_mat,
    caller = "computeCorrelation",
    observation_axis = "columns"
  )

  sample_label <- if (level == "grid") "spots" else "cells"
  list(expr_mat = expr_mat, grid_name = grid_name, sample_label = sample_label)
}

#' Internal helper for correlation workflows
#' @description
#' Internal helper for `.compute_correlation_restore_env`.
#' @param old_blas_env Internal parameter
#' @param old_mkl_env Internal parameter
#' @return Internal helper result
#' @keywords internal
.compute_correlation_restore_env <- function(old_blas_env, old_mkl_env) {
  if (is.na(old_blas_env)) {
    Sys.unsetenv("OPENBLAS_NUM_THREADS")
  } else {
    Sys.setenv(OPENBLAS_NUM_THREADS = old_blas_env)
  }
  if (is.na(old_mkl_env)) {
    Sys.unsetenv("MKL_NUM_THREADS")
  } else {
    Sys.setenv(MKL_NUM_THREADS = old_mkl_env)
  }
}

#' Internal helper for correlation workflows
#' @description
#' Internal helper for `.compute_correlation_resolve_runtime_config`.
#' @param ncores Internal parameter
#' @param verbose Internal parameter
#' @return Internal helper result
#' @keywords internal
.compute_correlation_resolve_runtime_config <- function(ncores, verbose) {
  .log_info("computeCorrelation", "S01", "Configuring thread settings...", verbose)
  max_cores <- detectCores()
  if (is.na(max_cores) || max_cores <= 0) {
    max_cores <- 1
    .log_info("computeCorrelation", "S01",
      "!!! Warning: Could not detect cores, using single core !!!", verbose
    )
  }
  ncores_safe <- max(1L, min(ncores, max_cores))
  thread_source <- if (ncores_safe < ncores) "clamped" else "requested"
  .log_backend("computeCorrelation", "S01", "threads",
    paste0(ncores_safe, " source=", thread_source, " requested=", ncores),
    reason = if (thread_source == "clamped") "clamped" else NULL,
    verbose = verbose
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
  }
  old_blas_env <- Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA)
  old_mkl_env <- Sys.getenv("MKL_NUM_THREADS", unset = NA)
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1")
  .log_backend("computeCorrelation", "S01", "blas_threads",
    "OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1", verbose = verbose
  )
  list(ncores_safe = ncores_safe, old_blas_env = old_blas_env, old_mkl_env = old_mkl_env)
}

#' Internal helper for data construction helpers
#' @description
#' Internal helper for `.build_scope_from_standardized_input_payload`.
#' @param payload Internal parameter
#' @param filter_genes Parameter value.
#' @param max_distance_to_nucleus Numeric threshold.
#' @param filter_molecules Parameter value.
#' @param exclude_prefixes Parameter value.
#' @param min_quality Numeric threshold.
#' @param roi_path Filesystem path.
#' @param num_workers Parameter value.
#' @param chunk_size Parameter value.
#' @param min_gene_types Numeric threshold.
#' @param max_gene_types Numeric threshold.
#' @param min_segmentation_points Numeric threshold.
#' @param keep_partial_tiles Logical flag.
#' @param verbose Logical; whether to emit progress messages.
#' @param flip_y Parameter value.
#' @param data_type Parameter value.
#' @return Internal helper result
#' @keywords internal
.build_scope_from_standardized_input_payload <- function(payload,
                                                        filter_genes = NULL,
                                                        max_distance_to_nucleus = 25,
                                                        filter_molecules = TRUE,
                                                        exclude_prefixes = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword",
                                                                             "SystemControl", "Negative"),
                                                        min_quality = 20,
                                                        roi_path = NULL,
                                                        exclude_qc_unsuitable = FALSE,
                                                        qc_out_dir = NULL,
                                                        qc_flag = NULL,
                                                        num_workers = 1L,
                                                        chunk_size = 5e5,
                                                        min_gene_types = 1,
                                                        max_gene_types = Inf,
                                                        min_segmentation_points = 0,
                                                        grid_keep_mode = c("segmentation_hybrid", "molecule_driven", "vertex_only", "centroid_only"),
                                                        grid_center_min_molecules = "q80",
                                                        keep_partial_tiles = FALSE,
                                                        verbose = TRUE,
                                                        flip_y = TRUE,
                                                        data_type = NULL) {
    if (is.null(data_type)) {
        data_type <- if (!is.null(payload$data_type)) payload$data_type else "xenium"
    }
    grid_keep_mode <- match.arg(grid_keep_mode)
    .trace_contract_field(payload, "grid_length", ".build_scope_from_standardized_input_payload.input")
    state <- .initialize_dataset_builder_state(
        input_dir = payload$work_dir,
        grid_length = payload$grid_length,
        seg_type = payload$segmentation_strategy,
        filter_genes = filter_genes,
        max_dist_mol_nuc = max_distance_to_nucleus,
        filtermolecule = filter_molecules,
        exclude_prefix = exclude_prefixes,
        filterqv = min_quality,
        coord_file = roi_path,
        ncores = num_workers,
        chunk_pts = chunk_size,
        exclude_qc_unsuitable = exclude_qc_unsuitable,
        qc_out_dir = qc_out_dir,
        qc_flag = qc_flag,
        min_gene_types = min_gene_types,
        max_gene_types = max_gene_types,
        min_seg_points = min_segmentation_points,
        grid_keep_mode = grid_keep_mode,
        grid_center_min_molecules = grid_center_min_molecules,
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

    .with_data_table_attached(function() {
        state_inner <- state
        state_inner <- .load_dataset_dependencies(state_inner)
        state_inner <- .resolve_input_paths(state_inner)
        state_inner <- .load_centroid_table(state_inner)
        state_inner <- .open_transcript_dataset(state_inner)
        state_inner <- .apply_dataset_filters(state_inner)
        state_inner <- .prepare_roi_geometry(state_inner)
        state_inner <- .clip_points_within_roi(state_inner)
        state_inner <- .initialize_output_object(state_inner)
        state_inner <- .ingest_segmentation_geometries(state_inner)
        state_inner <- .prefetch_roi_molecules(state_inner)
        state_inner <- .scan_molecule_batches(state_inner)
        state_inner <- .build_grid_layers(state_inner)
        state_inner <- .finalize_output_object(state_inner)
        state_inner$objects$scope_obj
    }, context = "[geneSCOPE::.create_scope] ingest pipeline")
}

#' Attach histology assets to a scope object
#' @description
#' Internal helper for `.add_scope_histology`.
#' Adds Visium-style histology images and scalefactors to an existing grid layer.
#' @param scope_obj A `scope_object` to annotate.
#' @param grid_name Grid layer name to update.
#' @param png_path Path to the histology PNG.
#' @param json_path Optional Visium scalefactors JSON; used when `scalefactors` is not provided.
#' @param level Image resolution level (`lowres` or `hires`).
#' @param coord_type Coordinate system for scaling (`visium` or `manual`).
#' @param scalefactors Optional list of scalefactors to use instead of the JSON file.
#' @param y_origin Coordinate origin policy (`auto`, `top-left`, `bottom-left`).
#' @param crop_bbox_px Parameter value.
#' @param roi_bbox Parameter value.
#' @return Updated `scope_object` with histology metadata.
#' @keywords internal
.add_scope_histology <- function(scope_obj,
                              grid_name,
                              png_path,
                              json_path = NULL,
                              level = c("lowres", "hires"),
                              crop_bbox_px = NULL,
                              roi_bbox = NULL,
                              coord_type = c("visium", "manual"),
                              scalefactors = NULL,
                              y_origin = c("auto", "top-left", "bottom-left")) {
    if (missing(level) || is.null(level) || !length(level)) {
        base <- tolower(basename(png_path))
        if (grepl("aligned_tissue_image", base, fixed = TRUE) ||
            grepl("tissue_hires_image", base, fixed = TRUE)) {
            level <- "hires"
        } else if (grepl("tissue_lowres_image", base, fixed = TRUE) ||
            grepl("detected_tissue_image", base, fixed = TRUE)) {
            level <- "lowres"
        } else {
            level <- "lowres"
        }
    } else {
        level <- match.arg(level)
    }
    coord_type <- match.arg(coord_type)
    y_origin <- .normalize_histology_y_origin(y_origin)
    if (!is(scope_obj, "scope_object")) {
        stop("`scope_obj` must be a scope_object.")
    }
    if (is.null(grid_name) || !nzchar(grid_name) || !(grid_name %in% names(scope_obj@grid))) {
        stop("grid '", grid_name, "' does not exist in scope_obj@grid.")
    }
    if (is.null(png_path) || !file.exists(png_path)) {
        stop("Image file not found: ", png_path)
    }
    sf_use <- scalefactors
    if (is.null(sf_use)) {
        if (is.null(json_path)) {
            stop("Missing json_path or scalefactors; at least one source is required.")
        }
        sf_use <- .read_scalefactors_json(json_path, warn_if_missing = TRUE)
    }
    grid_layer <- scope_obj@grid[[grid_name]]
    y_flip_ref_arg <- NULL
    if (!is.null(grid_layer$image_info) && !is.null(grid_layer$image_info$y_flip_reference_um)) {
        cand <- suppressWarnings(as.numeric(grid_layer$image_info$y_flip_reference_um))
        if (length(cand) && is.finite(cand[1])) y_flip_ref_arg <- cand[1]
    }
    .attach_histology(
        scope_obj = scope_obj,
        grid_name = grid_name,
        png_path = png_path,
        json_path = json_path,
        level = level,
        crop_bbox_px = crop_bbox_px,
        roi_bbox = roi_bbox,
        coord_type = coord_type,
        scalefactors = sf_use,
        y_origin = y_origin,
        y_flip_reference_um = y_flip_ref_arg
    )
}

#' Add single-cell spatial data
#' @description
#' Internal helper for `.add_single_cells`.
#' Dispatches to platform-specific loaders for CosMx or Xenium inputs and attaches
#' single-cell layers onto a scope object.
#' @param scope_obj A `scope_object` to extend.
#' @param xenium_dir Optional Xenium run directory.
#' @param cosmx_root Optional CosMx project root.
#' @param platform Deprecated platform hint (`auto`, `Xenium`, or `CosMx`).
#' @param prefer Preferred single-cell loader (`auto`, `xenium`, or `cosmx`).
#' @param verbose Emit progress messages when TRUE.
#' @param ... Additional arguments passed to the platform-specific loader.
#'   `data_dir` is accepted as a legacy alias for `xenium_dir` / `cosmx_root`.
#' @return Updated `scope_object` with single-cell layers.
#' @keywords internal
.add_single_cells <- function(scope_obj,
                           xenium_dir = NULL,
                           cosmx_root = NULL,
                           platform = NULL,
                           prefer = c("auto", "xenium", "cosmx"),
                           verbose = TRUE,
                           ...) {
  stopifnot(inherits(scope_obj, "scope_object"))
  parent <- "createSCOPE"
  step01 <- .log_step(parent, "S01", "dispatch single-cell loader", verbose)
  dots <- list(...)
  prefer_missing <- missing(prefer)
  prefer <- match.arg(prefer)
  if (!is.null(platform)) {
    platform_value <- match.arg(platform, c("auto", "Xenium", "CosMx"))
    platform_prefer <- switch(tolower(platform_value),
      auto = "auto",
      xenium = "xenium",
      cosmx = "cosmx"
    )
    .Deprecated(
      "prefer",
      package = "geneSCOPE",
      msg = "addSingleCells(platform=) is deprecated; use addSingleCells(prefer=) instead."
    )
    if (prefer_missing) {
      prefer <- platform_prefer
    } else if (!identical(platform_prefer, prefer)) {
      stop("addSingleCells(platform=) conflicts with prefer=. Use only prefer= to select the single-cell loader.", call. = FALSE)
    }
  }
  step01$enter(paste0("prefer=", prefer))

  # Legacy alias: allow `data_dir` to specify platform directory without
  # changing exported addSingleCells() formals.
  if (!is.null(dots$data_dir)) {
    data_dir <- dots$data_dir
    dots$data_dir <- NULL
    if (is.null(xenium_dir) && is.null(cosmx_root)) {
      guess <- .guess_singlecell_platform(data_dir)
      if (!is.null(guess$platform)) {
        if (identical(guess$platform, "Xenium")) {
          xenium_dir <- guess$path
        } else if (identical(guess$platform, "CosMx")) {
          cosmx_root <- guess$path
        }
      } else if (!is.null(data_dir) && nzchar(data_dir)) {
        if (identical(prefer, "cosmx")) {
          cosmx_root <- data_dir
        } else if (identical(prefer, "xenium")) {
          xenium_dir <- data_dir
        }
      }
    }
  }

  # 1) If user passed an explicit directory, prefer that
  pick <- NULL
  if (!is.null(xenium_dir)) pick <- "Xenium"
  if (!is.null(cosmx_root)) pick <- "CosMx"

  # 2) Auto-detect from scope_obj if not forced
  if (is.null(pick) && prefer == "auto") {
    plat <- NA_character_
    if (!is.null(scope_obj@meta.data) && nrow(scope_obj@meta.data)) {
      if ("platform" %in% names(scope_obj@meta.data)) {
        v <- scope_obj@meta.data$platform
        if (length(v)) plat <- as.character(v[which.max(tabulate(match(v, unique(v))))])
      }
    }
    if (is.na(plat) && !is.null(scope_obj@stats$platform)) plat <- as.character(scope_obj@stats$platform)
    if (!is.na(plat)) {
      if (grepl("cosmx", plat, ignore.case = TRUE)) pick <- "CosMx"
      else if (grepl("xenium", plat, ignore.case = TRUE)) pick <- "Xenium"
    }
  }

  # 3) Fallback default
  if (is.null(pick) && !identical(prefer, "auto")) {
    pick <- if (identical(prefer, "cosmx")) "CosMx" else "Xenium"
  }
  if (is.null(pick)) pick <- if (!is.null(xenium_dir)) "Xenium" else if (!is.null(cosmx_root)) "CosMx" else "Xenium"

  if (identical(pick, "CosMx")) {
    .log_info(parent, "S01", "dispatching to .add_single_cells_cosmx()", verbose)
    step01$done("builder=cosmx")
    cr <- if (!is.null(cosmx_root)) cosmx_root else getOption("geneSCOPE.cosmx_root")
    return(do.call(.add_single_cells_cosmx, c(list(scope_obj = scope_obj, cosmx_root = cr, verbose = verbose), dots)))
  } else {
    .log_info(parent, "S01", "dispatching to .add_single_cells_xenium()", verbose)
    step01$done("builder=xenium")
    xd <- if (!is.null(xenium_dir)) xenium_dir else getOption("geneSCOPE.xenium_dir")
    return(do.call(.add_single_cells_xenium, c(list(scope_obj = scope_obj, xenium_dir = xd, verbose = verbose), dots)))
  }
}

#' Create a scope object from platform outputs
#' @description
#' Internal helper for `.create_scope`.
#' Auto-detects Xenium/CosMx/Visium inputs under `data_dir` (or explicit
#' platform-specific arguments) and dispatches to the matching builder.
#' @param data_dir Root directory containing Xenium/CosMx/Visium outputs.
#' @param prefer Override platform detection (`auto`, `xenium`, `cosmx`, `visium`).
#' @param verbose Emit progress messages when TRUE.
#' @param sctransform Run SCTransform when using the Visium Seurat builder.
#' @param ... Additional arguments (currently unused).
#' @return A constructed `scope_object`.
#' @keywords internal
.create_scope <- function(data_dir = NULL,
                        prefer = c("auto", "xenium", "cosmx", "visium"),
                        verbose = TRUE,
                        sctransform = TRUE,
                        ...) {
    dots <- list(...)
    parent <- "createSCOPE"
    step01 <- .log_step(parent, "S01", "parse inputs and dispatch builder", verbose)
    prefer_missing <- missing(prefer)
    if (!is.null(dots$platform)) {
        platform <- tolower(as.character(dots$platform)[1])
        if (!platform %in% c("auto", "xenium", "cosmx", "visium")) {
            stop("createSCOPE(platform=) must be one of 'auto', 'xenium', 'cosmx', or 'visium'. Use prefer= instead.", call. = FALSE)
        }
        .Deprecated(
            "prefer",
            package = "geneSCOPE",
            msg = "createSCOPE(platform=) is deprecated and no longer silently ignored; use createSCOPE(prefer=) instead."
        )
        if (prefer_missing) {
            prefer <- platform
        } else {
            prefer_requested <- match.arg(prefer)
            if (!identical(platform, prefer_requested)) {
                stop("createSCOPE(platform=) conflicts with prefer=. Use only prefer= to select the platform.", call. = FALSE)
            }
        }
        dots$platform <- NULL
    }
    if (is.null(dots$verbose)) dots$verbose <- verbose
    if (is.null(data_dir)) {
        if (!is.null(dots$xenium_dir)) data_dir <- dots$xenium_dir
        else if (!is.null(dots$visium_dir)) data_dir <- dots$visium_dir
        else if (!is.null(dots$cosmx_root)) data_dir <- dots$cosmx_root
    }
    if (is.null(data_dir)) stop("Please provide data_dir (or xenium_dir/visium_dir/cosmx_root).")
    if (!dir.exists(data_dir)) stop("Directory not found: ", data_dir)

    prefer <- match.arg(prefer)
    step01$enter(paste0("prefer=", prefer, " data_dir=", basename(data_dir)))

    is_xenium <- function(d) {
        file.exists(file.path(d, "transcripts.parquet")) ||
            file.exists(file.path(d, "cells.parquet"))
    }
    is_visium <- function(d) {
        dir.exists(file.path(d, "spatial")) && (
            file.exists(file.path(d, "filtered_feature_bc_matrix.h5")) ||
            dir.exists(file.path(d, "filtered_feature_bc_matrix")) ||
            dir.exists(file.path(d, "raw_feature_bc_matrix")) ||
            file.exists(file.path(d, "spatial", "scalefactors_json.json"))
        )
    }
    is_cosmx <- function(d) {
        dir.exists(file.path(d, "flatFiles")) ||
            dir.exists(file.path(d, "DecodedFiles"))
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
    .log_info(parent, "S01", paste0("selected builder=", pick), verbose)

    dots$xenium_dir <- NULL; dots$visium_dir <- NULL; dots$cosmx_root <- NULL

    if (pick == "xenium") {
        .log_info(parent, "S01", "dispatching to builder", verbose)
        if (is.null(dots$grid_length)) stop("grid_length is required for Xenium. Pass grid_length= ...")
        step01$done("builder=xenium")
        dots$xenium_dir <- data_dir
        return(do.call(createSCOPE_xenium, dots))
    }
    if (pick == "cosmx") {
        .log_info(parent, "S01", "dispatching to builder", verbose)
        if (is.null(dots$grid_length)) stop("grid_length is required for CosMx. Pass grid_length= ...")
        step01$done("builder=cosmx")
        dots$cosmx_root <- data_dir
        return(do.call(createSCOPE_cosmx, dots))
    }

    dots$visium_dir <- data_dir
    if (is.null(dots$run_sctransform)) dots$run_sctransform <- isTRUE(sctransform)
    if (!isTRUE(sctransform) || !isTRUE(dots$run_sctransform)) {
        stop(
            "Visium construction requires sctransform=TRUE in this release. ",
            "The non-SCT Visium builder is intentionally unsupported because ",
            "previous non-SCT closures produced q-guard fallback and unstable ",
            "biological interpretation. Re-run createSCOPE(..., prefer='visium', ",
            "sctransform=TRUE).",
            call. = FALSE
        )
    }
    if (isTRUE(sctransform)) {
        .log_info(
            parent,
            "S01",
            paste0("dispatching to Seurat pipeline (SCTransform=", if (is.null(dots$run_sctransform)) "NULL" else as.character(dots$run_sctransform), ")"),
            verbose
        )
        step01$done("builder=visium_seurat")
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
        res <- do.call(.build_scope_from_visium_seurat, dots)
        res
    }
}

# ---- Underscore-named helpers (moved from api_data_construction.R) ----
#' Build a scope object from CosMx outputs
#' @description
#' Consumes CosMx parquet/flatFile exports and constructs a `scope_object`
#' with grid layers, metadata, and segmentation.
#' @param input_dir Root directory containing CosMx outputs.
#' @param grid_sizes Grid side length(s) to use for tiling.
#' @param segmentation_strategy Cell/nucleus segmentation strategy.
#' @param dataset_id Optional dataset identifier for naming outputs.
#' @param pixel_size_um Pixel size (microns) for coordinate scaling.
#' @param verbose Emit progress messages when TRUE.
#' @param allow_flatfile_generation Generate parquet from flatFiles when parquet is missing.
#' @param derive_cells_from_polygons Derive cell centroids from polygon CSV when needed.
#' @param data_type Data type label (default `cosmx`).
#' @param ... Additional arguments (currently unused).
#' @return A constructed `scope_object`.
#' @keywords internal
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
    parent <- "createSCOPE"

    .log_info(parent, "S03", paste0("input_dir=", input_dir), verbose)

    transcripts_path <- file.path(input_dir, "transcripts.parquet")
    cells_path <- file.path(input_dir, "cells.parquet")
    primary_cell_boundary <- file.path(input_dir, "cell_boundaries.parquet")
    fallback_cell_boundary <- file.path(input_dir, "segmentation_boundaries.parquet")
    nucleus_boundary <- file.path(input_dir, "nucleus_boundaries.parquet")
    .log_key_values(parent, list(
        platform = "cosmx",
        layout_flatFiles = dir.exists(file.path(input_dir, "flatFiles")),
        layout_DecodedFiles = dir.exists(file.path(input_dir, "DecodedFiles")),
        transcripts_parquet = file.exists(transcripts_path),
        cells_parquet = file.exists(cells_path),
        tx_file = file.exists(file.path(input_dir, "tx_file.csv")) || file.exists(file.path(input_dir, "tx_file.csv.gz")),
        exprMat = file.exists(file.path(input_dir, "exprMat_file.csv")) || file.exists(file.path(input_dir, "exprMat_file.csv.gz"))
    ))

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
        .log_info(parent, "S04", "transcripts.parquet not found; generating from flatFiles", verbose)
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
        # data.table loaded via Imports; arrow loaded via Suggests
        .log_info(parent, "S04", paste0("reading ", tx_path), verbose)
        tx <- tryCatch({ fread(tx_path) }, error = function(e) stop("Failed to read ", tx_path, ": ", e$message))
        cn <- tolower(names(tx)); names(tx) <- cn
        headerless <- all(grepl("^v[0-9]+$", cn))
        if (!headerless && all(c("x_global_px", "y_global_px", "target") %in% cn)) {
            xg <- tx[["x_global_px"]]; yg <- tx[["y_global_px"]]; tgt <- tx[["target"]]
            qv <- if ("qv" %in% cn) tx[["qv"]] else Inf
            nd <- if ("nucleus_distance" %in% cn) tx[["nucleus_distance"]] else 0
        } else {
            if (ncol(tx) < 9L) stop("headerless tx_file has insufficient columns (<9)")
            .log_info(parent, "S04", "tx_file headerless; using positional mapping col6/col7/col9 (col8 nucleus_distance)", verbose)
            xg <- tx[[6]]; yg <- tx[[7]]; tgt <- tx[[9]]
            qv <- if (ncol(tx) >= 10 && is.numeric(tx[[10]])) tx[[10]] else Inf
            nd <- if (ncol(tx) >= 8 && is.numeric(tx[[8]])) tx[[8]] else 0
        }
        tx_out <- data.table(
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
        .log_info(parent, "S04", paste0("writing ", tf_parq), verbose)
        arrow::write_parquet(tx_out, tf_parq)
    } else if (verbose) {
        .log_info(parent, "S04", "using existing transcripts.parquet", verbose)
    }

    have_cell_seg <- file.exists(seg_cell) || file.exists(seg_any)
    have_nuc_seg <- file.exists(nuc_parq)

    if (seg_type %in% c("cell", "both") && !have_cell_seg) {
        .log_info(parent, "S11", "cell segmentation parquet not found; downgrade segmentation strategy", verbose)
        seg_type <- if (have_nuc_seg) "nucleus" else "none"
    }
    if (seg_type %in% c("nucleus", "both") && !have_nuc_seg) {
        .log_info(parent, "S11", "nucleus segmentation parquet not found; downgrade segmentation strategy", verbose)
        seg_type <- if (have_cell_seg) "cell" else "none"
    }

    if (!file.exists(cell_parq) && build_from_flatfiles && build_cells_from_polygons) {
        poly_path <- NULL
        ff_dir <- file.path(input_dir, "flatFiles")
        cand <- c(Sys.glob(file.path(ff_dir, "*/*-polygons.csv.gz")), Sys.glob(file.path(ff_dir, "*/*_polygons.csv.gz")))
        if (length(cand)) poly_path <- cand[[1]]
        if (!is.null(poly_path)) {
            .log_info(parent, "S04", paste0("cells.parquet not found; deriving centroids from polygons: ", basename(poly_path)), verbose)
            pol <- fread(poly_path)
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
            dt <- data.table(fov = as.integer(fov), cell_id = as.character(cid),
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
            .log_info(parent, "S04", paste0("writing ", cell_parq), verbose)
            arrow::write_parquet(cent, cell_parq)
        } else if (verbose) {
            .log_info(parent, "S04", "polygons not found; cannot derive cells.parquet", verbose)
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

    .log_info(parent, "S01", "delegating to parquet builder", verbose)
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

#' Build a scope object from 10x Visium data
#' @param input_dir Visium output directory containing feature matrices.
#' @param use_filtered_matrix Use the filtered matrix when available.
#' @param include_in_tissue_spots Restrict to spots flagged as in-tissue.
#' @param grid_multiplier Optional multiplier for grid resolution.
#' @param flip_y Whether to flip the y-axis when loading spot coordinates.
#' @param verbose Emit progress messages when TRUE.
#' @return A constructed `scope_object`.
#' @keywords internal
build_scope_from_visium <- function(input_dir,
                                    use_filtered_matrix = TRUE,
                                    include_in_tissue_spots = TRUE,
                                    grid_multiplier = 1,
                                    flip_y = FALSE,
                                    verbose = TRUE) {
    stopifnot(dir.exists(input_dir))
    parent <- "createSCOPE"
    step15 <- .log_step(parent, "S15", "done summary", verbose)

    step02 <- .log_step(parent, "S02", "configure runtime resources and threads", verbose)
    step02$enter("visium builder")
    omp_info <- tryCatch(.native_openmp_info(), error = function(e) NULL)
    if (!is.null(omp_info)) {
        .log_info(
            parent,
            "S02",
            paste0(
                "openmp compiled=",
                ifelse(is.na(omp_info$compiled_with_openmp), "NA", omp_info$compiled_with_openmp),
                " omp_max_threads=", ifelse(is.null(omp_info$omp_max_threads), "NA", omp_info$omp_max_threads),
                " omp_num_procs=", ifelse(is.null(omp_info$omp_num_procs), "NA", omp_info$omp_num_procs)
            ),
            verbose
        )
    } else {
        .log_info(parent, "S02", "openmp info unavailable", verbose)
    }
    step02$done("skipped")

    step03 <- .log_step(parent, "S03", "initialize ingest environment", verbose)
    step03$enter(paste0("input_dir=", basename(input_dir)))
    step03$done("ok")

    # use input_dir directly (no alias)
    if (!is.numeric(grid_multiplier) || length(grid_multiplier) != 1L || !is.finite(grid_multiplier)) {
        stop("grid_multiplier must be a single finite numeric value.")
    }
    if (grid_multiplier <= 0) {
        stop("grid_multiplier must be positive.")
    }

    .with_data_table_attached(function() {
        step04 <- .log_step(parent, "S04", "load centroids and transcripts", verbose)
        step04$enter("start")
        candidate_dirs <- unique(file.path(input_dir, c(
            if (isTRUE(use_filtered_matrix)) "filtered_feature_bc_matrix",
            if (!isTRUE(use_filtered_matrix)) "raw_feature_bc_matrix",
            "filtered_feature_bc_matrix",
            "raw_feature_bc_matrix"
        )))
        matrix_dir <- candidate_dirs[dir.exists(candidate_dirs)][1]

        if (is.na(matrix_dir)) {
            stop("Could not find filtered_feature_bc_matrix/ or raw_feature_bc_matrix/ under ", input_dir)
        }
        .log_info(parent, "S04", paste0("using matrix directory: ", basename(matrix_dir)), verbose)

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

        counts_mt <- readMM(if (grepl("\\.gz$", matrix_file, ignore.case = TRUE)) gzfile(matrix_file) else matrix_file)
        counts_mat <- as(counts_mt, "CsparseMatrix")
        rm(counts_mt)

        features_dt <- fread(feature_file, header = FALSE)
        barcodes_dt <- fread(barcode_file, header = FALSE)
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

        spatial_dir <- file.path(input_dir, "spatial")
        if (!dir.exists(spatial_dir)) {
            stop("spatial/ directory not found under ", input_dir)
        }

        pos_candidates <- c(
            file.path(spatial_dir, "tissue_positions_list.csv"),
            file.path(spatial_dir, "tissue_positions.csv")
        )
        pos_path <- pos_candidates[file.exists(pos_candidates)][1]
        if (is.na(pos_path)) {
            stop("Could not find tissue_positions_list.csv or tissue_positions.csv under ", spatial_dir)
        }
        pos_raw <- fread(pos_path)
    .log_info(parent, "S04", paste0("loaded ", nrow(pos_raw), " spot positions"), verbose)

    standardize_positions <- function(dt) {
        if (!is.data.table(dt)) dt <- as.data.table(dt)
        if (ncol(dt) < 6) {
            stop("Visium position file must contain at least six columns.")
        }
        if (all(grepl("^V[0-9]+$", names(dt)))) {
            cols <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
            setnames(dt, cols[seq_len(ncol(dt))])
        } else {
            rename_first <- function(target, aliases) {
                curr_lower <- tolower(names(dt))
                if (length(which(curr_lower %in% aliases))) {
                    setnames(dt, which(curr_lower %in% aliases)[1], target)
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
        if (length(setdiff(required, names(dt)))) {
            stop("Position file missing required columns: ", paste(setdiff(required, names(dt)), collapse = ", "))
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
    if (isTRUE(include_in_tissue_spots)) {
        pos_dt <- pos_dt[in_tissue == 1L]
        .log_info(parent, "S04", paste0("retained ", nrow(pos_dt), " in-tissue spots"), verbose)
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

    nonzero_genes <- rowSums(counts_mat) > 0
    if (!any(nonzero_genes)) {
        stop("All genes have zero counts across the selected spots.")
    }
    if (any(!nonzero_genes)) {
        counts_mat <- counts_mat[nonzero_genes, , drop = FALSE]
    }
    gene_meta <- gene_meta[rownames(counts_mat), , drop = FALSE]
    step04$done(paste0("spots=", nrow(pos_dt), " genes=", nrow(counts_mat)))

    step05 <- .log_step(parent, "S05", "filter by nucleus distance", verbose)
    step05$enter("skipped")
    step05$done("skipped")

    step06 <- .log_step(parent, "S06", "exclude molecule prefixes", verbose)
    step06$enter("skipped")
    step06$done("skipped")

    step07 <- .log_step(parent, "S07", "filter by quality score", verbose)
    step07$enter("skipped")
    step07$done("skipped")

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

    step08 <- .log_step(parent, "S08", "compute bounds and transforms", verbose)
    step08$enter(paste0("flip_y=", flip_y))
    pos_dt[, `:=`(
        x = pxl_col_in_fullres * microns_per_pixel,
        y = pxl_row_in_fullres * microns_per_pixel
    )]
    if (isTRUE(flip_y)) {
        ymax <- max(pos_dt$y, na.rm = TRUE)
        pos_dt[, y := ymax - y]
    }
    bounds_x <- range(pos_dt$x, na.rm = TRUE)
    bounds_y <- range(pos_dt$y, na.rm = TRUE)
    .log_info(parent, "S08", paste0("bounds x=[", round(bounds_x[1]), ",", round(bounds_x[2]), "] y=[", round(bounds_y[1]), ",", round(bounds_y[2]), "]"), verbose)
    step08$done("ok")

    step09 <- .log_step(parent, "S09", "prepare ROI polygon", verbose)
    step09$enter("skipped")
    step09$done("skipped")

    step10 <- .log_step(parent, "S10", "clip centroids to ROI", verbose)
    step10$enter("skipped")
    step10$done("skipped")

    step11 <- .log_step(parent, "S11", "ingest segmentation data", verbose)
    step11$enter("skipped")
    step11$done("skipped")

    step12 <- .log_step(parent, "S12", "prefetch ROI molecules", verbose)
    step12$enter("skipped")
    step12$done("skipped")

    step13 <- .log_step(parent, "S13", "build grid layers", verbose)
    step13$enter("skipped")
    step13$done("skipped")

    step14 <- .log_step(parent, "S14", "finalize scope object", verbose)
    step14$enter()

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
        include_in_tissue = include_in_tissue_spots,
        grid_multiplier = grid_multiplier
    )
    step14$done(paste0("spots=", nrow(pos_dt), " grid_multiplier=", grid_multiplier))
    step15$enter("summary")
    step15$done(paste0("spots=", nrow(pos_dt), " genes=", nrow(counts_mat)))
    scope_obj
    }, context = "[geneSCOPE::.create_scope] Visium builder")
}

# Internal helper for the Xenium builder that consumes the payload emitted by scope_input_module.

#' Build a scope object from Xenium outputs
#' @param input_dir Xenium run directory containing molecule and cell outputs.
#' @param grid_sizes Grid side length(s) for spatial binning.
#' @param segmentation_strategy Cell/nucleus segmentation strategy.
#' @param filter_genes Optional gene whitelist for filtering.
#' @param max_distance_to_nucleus Maximum nucleus distance when deriving cells.
#' @param filter_molecules Whether to drop molecules failing quality checks.
#' @param exclude_prefixes Transcript prefixes to exclude.
#' @param min_quality Minimum molecule quality score.
#' @param roi_path Optional ROI to subset the dataset.
#' @param num_workers Number of workers for preprocessing.
#' @param chunk_size Chunk size for streaming molecule data.
#' @param parallel_backend Parallel backend for ROI clipping (`auto`, `fork`, `psock`, `serial`).
#' @param min_gene_types Minimum gene types per cell to retain.
#' @param max_gene_types Maximum gene types per cell to retain.
#' @param min_segmentation_points Minimum segmentation vertices per cell.
#' @param keep_partial_tiles Keep partially filled tiles when TRUE.
#' @param verbose Emit progress messages when TRUE.
#' @param flip_y Whether to flip the Y axis on import.
#' @param data_type Data type label (default `xenium`).
#' @return A constructed `scope_object`.
#' @keywords internal
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
                                    parallel_backend = c("auto", "fork", "psock", "serial"),
                                    min_gene_types = 1,
                                    max_gene_types = Inf,
                                    min_segmentation_points = 0,
                                    grid_keep_mode = c("segmentation_hybrid", "molecule_driven", "vertex_only", "centroid_only"),
                                    grid_center_min_molecules = "q80",
                                    keep_partial_tiles = FALSE,
                                    verbose = TRUE,
                                    flip_y = TRUE,
                                    data_type = "xenium",
                                    exclude_qc_unsuitable = FALSE,
                                    qc_out_dir = NULL,
                                    qc_flag = NULL) {
    parallel_backend <- match.arg(parallel_backend)
    grid_keep_mode <- match.arg(grid_keep_mode)
    parent <- "createSCOPE"
    step15 <- .log_step(parent, "S15", "done summary", verbose)
    tx_total <- NA_integer_

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
        parallel_backend = parallel_backend,
        min_gene_types = min_gene_types,
        max_gene_types = max_gene_types,
        min_seg_points = min_segmentation_points,
        grid_keep_mode = grid_keep_mode,
        grid_center_min_molecules = grid_center_min_molecules,
        keep_partial_grid = keep_partial_tiles,
        verbose = verbose,
        flip_y = flip_y,
        data_type = data_type,
        exclude_qc_unsuitable = exclude_qc_unsuitable,
        qc_out_dir = qc_out_dir,
        qc_flag = qc_flag
    )

    state <- .configure_worker_threads(state)
    if (!is.null(state$thread_context$restore_fn)) {
        on.exit(state$thread_context$restore_fn(), add = TRUE)
    }
    state <- .limit_thread_environment(state)
    if (!is.null(state$env_restore)) on.exit(state$env_restore(), add = TRUE)

    step03 <- .log_step(parent, "S03", "initialize ingest environment", verbose)
    step03$enter(paste0("input_dir=", basename(input_dir)))
    state <- .load_dataset_dependencies(state)
    state <- .resolve_input_paths(state)
    step03$done("ok")

    step04 <- .log_step(parent, "S04", "load centroids and transcripts", verbose)
    step04$enter("start")
    state <- .load_centroid_table(state)
    state <- .open_transcript_dataset(state)
    if (.log_enabled(verbose)) {
        tx_total <- tryCatch({
            ds_obj <- state$datasets$ds
            if (inherits(ds_obj, "Dataset")) {
                as.integer(ds_obj$num_rows)
            } else {
                ds_obj |> summarise(n = n()) |> collect() |> pull(n)
            }
        }, error = function(e) NA_integer_)
    }
    step04$done(paste0("cells=", nrow(state$objects$centroids), " molecules=", ifelse(is.na(tx_total), "NA", tx_total)))
    state <- .apply_dataset_filters(state)
    state <- .prepare_roi_geometry(state)
    state <- .clip_points_within_roi(state)
    state <- .initialize_output_object(state)
    state <- .ingest_segmentation_geometries(state)
    state <- .prefetch_roi_molecules(state)
    state <- .scan_molecule_batches(state)
    state <- .build_grid_layers(state)
    state <- .finalize_output_object(state)

    scope_obj <- state$objects$scope_obj
    cell_count <- if (!is.null(scope_obj@coord$centroids)) nrow(scope_obj@coord$centroids) else 0L
    grid_layers <- length(scope_obj@grid)
    step15$enter("summary")
    step15$done(paste0("cells=", cell_count, " molecules=", ifelse(is.na(tx_total), "NA", tx_total), " grids=", grid_layers))
    scope_obj
}

#' Create a scope object from CosMx inputs
#' @description
#' Thin wrapper over `build_scope_from_cosmx` that normalizes parameter names
#' and provides a user-facing entrypoint.
#' @param cosmx_root Root directory containing CosMx outputs.
#' @param grid_length Grid side length(s) to build.
#' @param seg_type Segmentation strategy (`cell`, `nucleus`, `both`).
#' @param dataset_id Optional dataset identifier.
#' @param pixel_size_um Pixel size (microns) for coordinate scaling.
#' @param verbose Emit progress messages when TRUE.
#' @param build_from_flatfiles Allow generation from flatFiles when parquet is missing.
#' @param build_cells_from_polygons Derive cell centroids from polygon CSVs when needed.
#' @param data_type Data type label.
#' @param ... Additional arguments (currently unused).
#' @return A constructed `scope_object`.
#' @keywords internal
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

#' Create a scope object from Visium outputs
#' @description
#' Legacy Visium constructor. Visium construction is SCT-only in this release;
#' call `createSCOPE(..., prefer = "visium", sctransform = TRUE)` instead.
#' @param visium_dir Visium output directory (with `spatial/`).
#' @param use_filtered Use the filtered matrix when available.
#' @param include_in_tissue Restrict to spots flagged as in-tissue.
#' @param grid_multiplier Optional multiplier for grid resolution.
#' @param flip_y Whether to flip the y-axis when loading positions.
#' @param verbose Emit progress messages when TRUE.
#' @return A constructed `scope_object`.
#' @keywords internal
createSCOPE_visium <- function(visium_dir,
                               use_filtered = TRUE,
                               include_in_tissue = TRUE,
                               grid_multiplier = 1,
                               flip_y = FALSE,
                               verbose = TRUE) {
    stop(
        "createSCOPE_visium() uses the legacy non-SCT Visium builder and is ",
        "unsupported in this release. Use createSCOPE(..., prefer = 'visium', ",
        "sctransform = TRUE) so the Seurat/SCTransform Visium path constructs ",
        "the grid layer required by downstream geneSCOPE statistics.",
        call. = FALSE
    )
}

#' Create a scope object from Xenium outputs
#' @description
#' Wrapper over `build_scope_from_xenium` with legacy parameter names.
#' @param xenium_dir Xenium run directory.
#' @param grid_length Grid side length(s) to build.
#' @param seg_type Segmentation strategy (`cell`, `nucleus`, `both`).
#' @param data_type Data type label for the resulting scope object.
#' @param parallel_backend Parallel backend for ROI clipping (`auto`, `fork`, `psock`, `serial`).
#' @param ... Additional arguments (currently unused).
#' @return A constructed `scope_object`.
#' @keywords internal
createSCOPE_xenium <- function(xenium_dir,
                               grid_length,
                               seg_type = c("cell", "nucleus", "both"),
                               data_type = "xenium",
                               parallel_backend = c("auto", "fork", "psock", "serial"),
                               ...) {
    seg_type <- match.arg(seg_type)
    parallel_backend <- match.arg(parallel_backend)
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
            data_type = data_type,
            parallel_backend = parallel_backend
        ),
        dots
    )
    do.call(build_scope_from_xenium, build_args)
}
