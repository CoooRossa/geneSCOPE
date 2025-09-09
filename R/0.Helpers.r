#' @title Helper Functions for geneSCOPE Package
#' @description Internal helper functions used throughout the package

#' @noRd
.selectGridLayer <- function(coordObj, grid_name = NULL, verbose = FALSE) {
    # Removed verbose message to avoid redundancy with main functions

    if (is.null(coordObj@grid) || length(coordObj@grid) == 0) {
        stop("No grid layers found in coordObj")
    }

    if (is.null(grid_name)) {
        if (length(coordObj@grid) == 1) {
            # Removed verbose message to avoid redundancy
            return(coordObj@grid[[1]])
        } else {
            stop("Multiple grid layers found. Please specify grid_name.")
        }
    }

    if (!grid_name %in% names(coordObj@grid)) {
        stop("Grid layer '", grid_name, "' not found.")
    }

    # Removed verbose message to avoid redundancy with main functions
    return(coordObj@grid[[grid_name]])
}

#' @noRd
.checkGridContent <- function(coordObj, grid_name, verbose = FALSE) {
    # Removed verbose message to avoid redundancy with main functions

    g_layer <- .selectGridLayer(coordObj, grid_name)

    required_elements <- c("grid_info", "counts")
    missing <- setdiff(required_elements, names(g_layer))
    if (length(missing) > 0) {
        stop(
            "Grid layer '", grid_name, "' missing required elements: ",
            paste(missing, collapse = ", ")
        )
    }

    # Removed verbose message to avoid redundancy with main functions
    invisible(TRUE)
}

#' @noRd
.getGeneSubset <- function(coordObj, genes = NULL, cluster_col = NULL, cluster_num = NULL, verbose = FALSE) {
    # Removed verbose message to avoid redundancy with main functions

    if (!is.null(genes)) {
        return(genes)
    }

    if (!is.null(cluster_col) && !is.null(cluster_num)) {
        if (is.null(coordObj@meta.data) || !cluster_col %in% colnames(coordObj@meta.data)) {
            stop("Cluster column '", cluster_col, "' not found in meta.data")
        }

        cluster_values <- coordObj@meta.data[[cluster_col]]
        selected_genes <- rownames(coordObj@meta.data)[cluster_values == cluster_num]
        selected_genes <- selected_genes[!is.na(selected_genes)]

        if (length(selected_genes) == 0) {
            stop("No genes found for cluster ", cluster_num, " in column ", cluster_col)
        }

        return(selected_genes)
    }

    # If no subset specified, return all genes
    if (!is.null(coordObj@meta.data)) {
        all_genes <- rownames(coordObj@meta.data)
        return(all_genes)
    }

    stop("No gene subset specified and no meta.data available")
}

#' @noRd
.getLeeMatrix <- function(coordObj, grid_name = NULL, lee_layer = NULL, verbose = FALSE) {

    ## ---- 0. Select grid sub-layer ------------------------------------------------
    g_layer <- .selectGridLayer(coordObj, grid_name, verbose = verbose)
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
            "' does not contain a valid Lee's L matrix."
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

#' @noRd
.parse_q <- function(qstr) {
    if (!is.character(qstr) || length(qstr) != 1 ||
        !grepl("^[qQ][0-9]+\\.?[0-9]*$", qstr)) {
        stop("pct string must look like 'q90' / 'q99.5' …")
    }
    val <- as.numeric(sub("^[qQ]", "", qstr)) / 100
    if (is.na(val) || val < 0 || val > 1) {
        stop("pct out of range 0–100")
    }
    val
}

#' @noRd
.filter_matrix_by_quantile <- function(mat, pct_min = "q0", pct_max = "q100") {
    stopifnot(is.matrix(mat))
    pmin <- .parse_q(pct_min)
    pmax <- .parse_q(pct_max)
    if (pmin > pmax) stop("pct_min > pct_max")
    vec <- as.vector(mat)
    thrL <- as.numeric(stats::quantile(vec, pmin, na.rm = TRUE))
    thrU <- as.numeric(stats::quantile(vec, pmax, na.rm = TRUE))
    mat[vec < thrL | vec > thrU] <- 0
    Matrix::drop0(mat)
}
