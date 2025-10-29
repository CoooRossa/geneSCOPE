#' @title Helper Functions for geneSCOPE Package
#' @description Internal helper functions used throughout the package

#' @noRd
.selectGridLayer <- function(scope_obj, grid_name = NULL, verbose = FALSE) {
    # Removed verbose message to avoid redundancy with main functions

    if (is.null(scope_obj@grid) || length(scope_obj@grid) == 0) {
        stop("No grid layers found in scope_obj")
    }

    if (is.null(grid_name)) {
        if (length(scope_obj@grid) == 1) {
            # Removed verbose message to avoid redundancy
            return(scope_obj@grid[[1]])
        } else {
            stop("Multiple grid layers found. Please specify grid_name.")
        }
    }

    if (!grid_name %in% names(scope_obj@grid)) {
        stop("Grid layer '", grid_name, "' not found.")
    }

    # Removed verbose message to avoid redundancy with main functions
    return(scope_obj@grid[[grid_name]])
}

#' @noRd
.checkGridContent <- function(scope_obj, grid_name, verbose = FALSE) {
    # Removed verbose message to avoid redundancy with main functions

    g_layer <- .selectGridLayer(scope_obj, grid_name)

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
.getGeneSubset <- function(scope_obj, genes = NULL, cluster_col = NULL, cluster_num = NULL, verbose = FALSE) {
    # Removed verbose message to avoid redundancy with main functions

    if (!is.null(genes)) {
        return(genes)
    }

    if (!is.null(cluster_col) && !is.null(cluster_num)) {
        if (is.null(scope_obj@meta.data) || !cluster_col %in% colnames(scope_obj@meta.data)) {
            stop("Cluster column '", cluster_col, "' not found in meta.data")
        }

        cluster_values <- scope_obj@meta.data[[cluster_col]]
        selected_genes <- rownames(scope_obj@meta.data)[cluster_values == cluster_num]
        selected_genes <- selected_genes[!is.na(selected_genes)]

        if (length(selected_genes) == 0) {
            stop("No genes found for cluster ", cluster_num, " in column ", cluster_col)
        }

        return(selected_genes)
    }

    # If no subset specified, return all genes
    if (!is.null(scope_obj@meta.data)) {
        all_genes <- rownames(scope_obj@meta.data)
        return(all_genes)
    }

    stop("No gene subset specified and no meta.data available")
}

#' @noRd
.getLeeMatrix <- function(scope_obj, grid_name = NULL, lee_layer = NULL, verbose = FALSE) {
    ## ---- 0. Select grid sub-layer ------------------------------------------------
    g_layer <- .selectGridLayer(scope_obj, grid_name, verbose = verbose)
    if (is.null(grid_name)) { # Write back the actual name
        grid_name <- names(scope_obj@grid)[
            vapply(scope_obj@grid, identical, logical(1), g_layer)
        ]
    }

    ## ---- 1. Auto-detect LeeStats layer name ----------------------------------------
    if (is.null(lee_layer)) {
        cand <- character(0)
        if (!is.null(scope_obj@stats[[grid_name]])) {
            cand <- names(scope_obj@stats[[grid_name]])
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
    if (!is.null(scope_obj@stats[[grid_name]]) &&
        !is.null(scope_obj@stats[[grid_name]][[lee_layer]])) {
        leeStat <- scope_obj@stats[[grid_name]][[lee_layer]]
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
.getPearsonMatrix <- function(scope_obj,
                              grid_name = NULL,
                              level = c("grid", "cell")) {
    level <- match.arg(level)

    ## ---------- 1. Set layer name & final target ------------------------------------
    corr_name <- "pearson_cor"
    f_cell_suf <- "_cell" # single-cell layer suffix

    ## ---------- 2. Get matrix --------------------------------------------------
    if (level == "grid") {
        g_layer <- .selectGridLayer(scope_obj, grid_name)
        if (is.null(grid_name)) {
            grid_name <- names(scope_obj@grid)[
                vapply(scope_obj@grid, identical, logical(1), g_layer)
            ]
        }

        ##   2a. New version: @stats[[grid_name]]
        rmat <- if (!is.null(scope_obj@stats[[grid_name]]) &&
            !is.null(scope_obj@stats[[grid_name]][[corr_name]])) {
            scope_obj@stats[[grid_name]][[corr_name]]
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

        ## Lee's L alignment
        Lmat <- tryCatch(
            .getLeeMatrix(scope_obj, grid_name = grid_name),
            error = function(e) NULL
        )
        if (!is.null(Lmat)) {
            common <- intersect(rownames(rmat), rownames(Lmat))
            rmat <- rmat[common, common, drop = FALSE]
        }
    } else { # ---------- single-cell ----------

        ##   2a. New version: @stats[["cell"]]
        rmat <- if (!is.null(scope_obj@stats[["cell"]]) &&
            !is.null(scope_obj@stats[["cell"]][[paste0(corr_name, f_cell_suf)]])) {
            scope_obj@stats[["cell"]][[paste0(corr_name, f_cell_suf)]]
        } else {
            NULL
        }

        ##   2b. Fallback: @cells$pearson_cor
        if (is.null(rmat) && !is.null(scope_obj@cells[[corr_name]])) {
            rmat <- scope_obj@cells[[corr_name]]
        }

        if (is.null(rmat)) {
            stop("Pearson matrix at cell level not found in scope_obj.")
        }

        ## Lee's L alignment (if single-cell layer exists)
        Lmat <- NULL
        if (!is.null(scope_obj@stats[["cell"]])) {
            layer_L <- intersect(
                grep("^LeeStats_", names(scope_obj@stats[["cell"]]),
                    value = TRUE
                ),
                paste0("LeeStats_Xz", f_cell_suf)
            )
            if (length(layer_L)) {
                Lmat <- scope_obj@stats[["cell"]][[layer_L[1]]]$L
            }
        }
        if (is.null(Lmat) && !is.null(scope_obj@cells$LeeStats_Xz)) {
            Lmat <- scope_obj@cells$LeeStats_Xz$L
        }

        if (!is.null(Lmat)) {
            common <- intersect(rownames(rmat), rownames(Lmat))
            rmat <- rmat[common, common, drop = FALSE]
        }
    }

    return(rmat)
}

#' @noRd
#' Compute consensus weights on a given edge list without allocating a dense N×N matrix.
#' @param memb_mat Integer matrix (genes × runs), each column is a community assignment.
#' @param edge_i Integer vector (1-based) indices of edge endpoints (length m).
#' @param edge_j Integer vector (1-based) indices of edge endpoints (length m).
#' @return Numeric vector of length m: fraction of runs in which the two endpoints co-cluster.
.consensus_on_edges <- function(memb_mat, edge_i, edge_j, n_threads = NULL) {
    stopifnot(is.matrix(memb_mat), length(edge_i) == length(edge_j))
    m <- length(edge_i)
    if (m == 0L) return(numeric(0))
    runs <- ncol(memb_mat)
    if (is.null(n_threads)) n_threads <- getSafeThreadCount(8L)
    # Fast path via C++ when available
    if (exists("consensus_on_edges_omp", mode = "function")) {
        counts <- consensus_on_edges_omp(as.integer(edge_i), as.integer(edge_j),
                                         memb = matrix(as.integer(memb_mat), nrow(memb_mat), ncol(memb_mat)),
                                         n_threads = as.integer(n_threads))
        return(as.numeric(counts) / runs)
    }
    # R fallback (block-wise compare to limit memory)
    mm <- memb_mat
    storage.mode(mm) <- "integer"
    counts <- integer(m)
    block <- max(1L, min(64L, runs))
    for (start in seq.int(1L, runs, by = block)) {
        end <- min(runs, start + block - 1L)
        comp <- mm[edge_i, start:end, drop = FALSE] == mm[edge_j, start:end, drop = FALSE]
        counts <- counts + as.integer(rowSums(comp))
    }
    counts / runs
}

#' @noRd
#' Map CMH/weight from a similarity graph to (ei,ej) edges among kept_genes.
#' Falls back to data.table when C++ helper is unavailable.
.cmh_lookup_pairs <- function(kept_genes, ei, ej, g_sim, weight_col = c("CMH", "weight"), n_threads = NULL) {
    weight_col <- match.arg(weight_col)
    if (is.null(n_threads)) n_threads <- getSafeThreadCount(8L)
    ed_sim <- igraph::as_data_frame(g_sim, what = "edges")
    wcol <- if (!is.null(ed_sim[[weight_col]])) weight_col else if (!is.null(ed_sim$weight)) "weight" else stop("CMH graph lacks weight column")
    si <- match(ed_sim$from, kept_genes)
    sj <- match(ed_sim$to, kept_genes)
    ok <- which(!is.na(si) & !is.na(sj))
    si <- as.integer(si[ok]); sj <- as.integer(sj[ok]); sw <- as.numeric(ed_sim[[wcol]][ok])
    fallback <- stats::median(sw, na.rm = TRUE)
    if (exists("cmh_lookup_rcpp", mode = "function")) {
        return(cmh_lookup_rcpp(as.integer(ei), as.integer(ej), si, sj, sw, fallback, as.integer(n_threads)))
    }
    # R fallback: keyed join via data.table
    dt_sim <- data.table::data.table(from = pmin(ed_sim$from, ed_sim$to), to = pmax(ed_sim$from, ed_sim$to), w = ed_sim[[wcol]])
    data.table::setkey(dt_sim, from, to)
    g1 <- kept_genes[ei]; g2 <- kept_genes[ej]
    dt_lp <- data.table::data.table(from = pmin(g1, g2), to = pmax(g1, g2), idx = seq_along(ei))
    data.table::setkey(dt_lp, from, to)
    ans <- dt_lp[dt_sim]
    ww <- ans$w; ww[is.na(ww)] <- fallback
    ww
}

#' @noRd
#' Compute per-gene intra-cluster degree using only existing edges.
#' @param edge_i,j Integer endpoint indices (1-based) within the kept gene set.
#' @param memb Integer cluster labels (length = number of kept genes), may contain NA.
#' @return Integer vector of intra-cluster degrees per gene.
.intra_degree_from_edges <- function(edge_i, edge_j, memb) {
    n <- length(memb)
    if (length(edge_i) == 0L) return(integer(n))
    same <- (memb[edge_i] == memb[edge_j]) & !is.na(memb[edge_i]) & !is.na(memb[edge_j])
    if (!any(same)) return(integer(n))
    i_tab <- tabulate(edge_i[same], nbins = n)
    j_tab <- tabulate(edge_j[same], nbins = n)
    i_tab + j_tab
}

#' @noRd
#' Extract Pearson correlations only for the requested gene pairs without
#' constructing a full kept_genes × kept_genes dense submatrix.
#' @param scope_obj scope_object
#' @param grid_name character grid layer name
#' @param level "grid" or "cell" (default used by cluster code is "cell")
#' @param kept_genes character vector of genes (in the same order used by edge_i/j)
#' @param edge_i,j integer vectors (1-based) indexing into kept_genes
#' @return numeric vector of correlations for each pair (edge_i[k], edge_j[k]).
.get_pairwise_cor_for_edges <- function(scope_obj, grid_name, level = c("cell", "grid"),
                                        kept_genes, edge_i, edge_j) {
    level <- match.arg(level)
    rMat <- .getPearsonMatrix(scope_obj, grid_name = grid_name, level = level)
    if (!is.matrix(rMat)) {
        # Best-effort: try to coerce only if not file-backed; avoid huge subsetting
        rMat <- try(as.matrix(rMat), silent = TRUE)
        if (inherits(rMat, "try-error") || !is.matrix(rMat)) {
            stop("Pearson correlation matrix is not a base matrix; cannot extract pairs efficiently.")
        }
    }
    ridx <- match(kept_genes, rownames(rMat))
    if (anyNA(ridx)) {
        stop("Some kept_genes not present in Pearson matrix.")
    }
    rMat[cbind(ridx[edge_i], ridx[edge_j])]
}

#' @noRd
#' Sparse-safe symmetric pmax for matrices. Keeps sparse structure when possible.
.symmetric_pmax <- function(M) {
    if (inherits(M, "sparseMatrix")) {
        M <- methods::as(M, "dgCMatrix")
        M <- Matrix::pmax(M, Matrix::t(M))
        Matrix::diag(M) <- 0
        return(M)
    } else {
        M <- pmax(M, t(M))
        diag(M) <- 0
        return(M)
    }
}

#' @noRd
#' Cluster-level metrics without densifying matrices.
#' within_cons: mean consensus weight on intra-cluster edges.
#' conductance: Wcut/(Wcut+Win) using weighted adjacency W (sparse).
.cluster_metrics_sparse <- function(members, cons_mat, W_mat, genes) {
    cl_ids <- sort(na.omit(unique(members)))
    if (!inherits(cons_mat, "sparseMatrix")) cons_mat <- methods::as(cons_mat, "dgCMatrix")
    if (!inherits(W_mat, "sparseMatrix")) W_mat <- methods::as(W_mat, "dgCMatrix")
    TTc <- methods::as(cons_mat, "TsparseMatrix")
    maskU_c <- (TTc@i < TTc@j)
    ei_c <- TTc@i[maskU_c] + 1L; ej_c <- TTc@j[maskU_c] + 1L; ex_c <- TTc@x[maskU_c]
    TTw <- methods::as(W_mat, "TsparseMatrix")
    maskU_w <- (TTw@i < TTw@j)
    ei_w <- TTw@i[maskU_w] + 1L; ej_w <- TTw@j[maskU_w] + 1L; ex_w <- TTw@x[maskU_w]
    res <- lapply(cl_ids, function(cid) {
        idx <- which(members == cid)
        n <- length(idx)
        if (n <= 1) return(data.frame(cluster = cid, size = n, within_cons = NA_real_, conductance = NA_real_))
        in_cl <- logical(length(members)); in_cl[idx] <- TRUE
        m_in_c <- in_cl[ei_c] & in_cl[ej_c]
        within_cons <- if (any(m_in_c)) mean(ex_c[m_in_c], na.rm = TRUE) else NA_real_
        m_in_w <- in_cl[ei_w] & in_cl[ej_w]
        m_cut_w <- xor(in_cl[ei_w], in_cl[ej_w])
        Win <- if (any(m_in_w)) sum(ex_w[m_in_w]) else 0
        Wcut <- if (any(m_cut_w)) sum(ex_w[m_cut_w]) else 0
        conductance <- Wcut / (Wcut + Win + 1e-12)
        data.frame(cluster = cid, size = n, within_cons = within_cons, conductance = conductance)
    })
    do.call(rbind, res)
}

#' @noRd
#' Gene-level metrics without densifying: p_in, p_best_out, w_in
.gene_metrics_sparse <- function(members, cons_mat, W_mat, genes) {
    if (!inherits(cons_mat, "sparseMatrix")) cons_mat <- methods::as(cons_mat, "dgCMatrix")
    if (!inherits(W_mat, "sparseMatrix")) W_mat <- methods::as(W_mat, "dgCMatrix")
    cl_ids <- sort(na.omit(unique(members)))
    n <- length(members)
    p_in <- rep(NA_real_, n); p_best_out <- rep(NA_real_, n); w_in <- rep(NA_real_, n)
    for (cid in cl_ids) {
        idx <- which(members == cid)
        if (!length(idx)) next
        # p_in: row means on intra-cluster consensus (excluding diag)
        subC <- cons_mat[idx, idx, drop = FALSE]
        # diag is zero; rowSums(subC) equals sum of weights to same cluster
        denom <- pmax(length(idx) - 1L, 1L)
        p_in[idx] <- Matrix::rowSums(subC) / denom
        # w_in: weighted adjacency to same cluster
        subW <- W_mat[idx, idx, drop = FALSE]
        # subtract zero diag implicitly; rowSums suffices
        w_in[idx] <- Matrix::rowSums(subW)
        # p_best_out: for each other cluster, compute row mean and take max
        if (length(cl_ids) > 1) {
            best <- rep(0, length(idx))
            for (cj in setdiff(cl_ids, cid)) {
                jdx <- which(members == cj)
                if (!length(jdx)) next
                block <- cons_mat[idx, jdx, drop = FALSE]
                denom_j <- pmax(length(jdx), 1L)
                # row mean to cluster cj
                mean_vec <- Matrix::rowSums(block) / denom_j
                # track maximum
                best <- pmax(best, mean_vec)
            }
            p_best_out[idx] <- best
        } else {
            p_best_out[idx] <- 0
        }
    }
    data.frame(gene = genes, cluster = members, p_in = p_in, p_best_out = p_best_out, w_in = w_in, stringsAsFactors = FALSE)
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

    seg_raw <- data.table::as.data.table(arrow::read_parquet(seg_file))
    # 若存在 fov 列：构造唯一键 key = paste(cell_id, fov)，并把 label_id 也更新为该唯一键
    if ("fov" %in% names(seg_raw)) {
        seg_dt <- seg_raw[, .(
            cell     = paste0(as.character(cell_id), "_", as.character(fov)),
            x        = vertex_x,
            y        = vertex_y,
            label_id = paste0(as.character(cell_id), "_", as.character(fov))
        )]
    } else {
        seg_dt <- seg_raw[, .(
            cell     = as.character(cell_id),
            x        = vertex_x,
            y        = vertex_y,
            label_id = as.character(label_id)
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
    stopifnot(is.matrix(mat) || inherits(mat, "Matrix"))
    pmin <- .parse_q(pct_min)
    pmax <- .parse_q(pct_max)
    if (pmin > pmax) stop("pct_min > pct_max")

    # Fast path for sparse matrices: operate only on existing (upper-tri) entries
    if (inherits(mat, "sparseMatrix")) {
        TT <- methods::as(mat, "TsparseMatrix")
        # Use only upper triangle to avoid double work; keep sign information
        mU <- (TT@i < TT@j)
        if (!any(mU)) return(Matrix::drop0(mat * 0))
        xi <- TT@i[mU] + 1L
        xj <- TT@j[mU] + 1L
        xv <- TT@x[mU]
        # Split by sign; negatives use absolute value for quantiles (consistent with previous behavior)
        pos_m <- xv >= 0
        neg_m <- !pos_m
        keep_idx <- logical(length(xv))
        if (any(pos_m)) {
            pos_vec <- xv[pos_m]
            # Guard against all-NA or length-1 edge cases
            if (length(pos_vec) == 1L) {
                pos_thrL <- pos_vec; pos_thrU <- pos_vec
            } else {
                pos_thrL <- as.numeric(stats::quantile(pos_vec, pmin, na.rm = TRUE, type = 7))
                pos_thrU <- as.numeric(stats::quantile(pos_vec, pmax, na.rm = TRUE, type = 7))
            }
            keep_idx[pos_m] <- (pos_vec >= pos_thrL) & (pos_vec <= pos_thrU)
        }
        if (any(neg_m)) {
            neg_abs <- abs(xv[neg_m])
            if (length(neg_abs) == 1L) {
                neg_thrL <- neg_abs; neg_thrU <- neg_abs
            } else {
                neg_thrL <- as.numeric(stats::quantile(neg_abs, pmin, na.rm = TRUE, type = 7))
                neg_thrU <- as.numeric(stats::quantile(neg_abs, pmax, na.rm = TRUE, type = 7))
            }
            keep_idx[neg_m] <- (neg_abs >= neg_thrL) & (neg_abs <= neg_thrU)
        }
        if (!any(keep_idx)) return(Matrix::drop0(mat * 0))
        ii <- xi[keep_idx]; jj <- xj[keep_idx]; vv <- xv[keep_idx]
        # Materialize symmetric sparse matrix of kept edges
        res <- Matrix::sparseMatrix(
            i = c(ii, jj), j = c(jj, ii), x = c(vv, vv),
            dims = dim(mat), dimnames = dimnames(mat)
        )
        # Zero diagonal by construction; drop explicit zeros
        return(Matrix::drop0(res))
    }

    # Dense path: operate on upper-tri only; for very large matrices, use sampling to estimate quantiles
    n <- nrow(mat)
    if (n == 0L) return(mat)
    ut_mask <- upper.tri(mat, diag = FALSE)
    vals <- mat[ut_mask]
    # Split by sign
    pos_vals <- vals[vals >= 0]
    neg_vals_abs <- abs(vals[vals < 0])
    # Approximate quantiles for very large vectors to avoid O(N log N) on tens of millions of entries
    max_samples <- getOption("geneSCOPE.quantile_max_samples", 5e6L)
    qfun <- function(x, p) {
        if (!length(x)) return(NA_real_)
        if (length(x) > max_samples) {
            # Sample without replacement for a robust approximation
            idx <- sample.int(length(x), max_samples)
            x <- x[idx]
        }
        as.numeric(stats::quantile(x, p, na.rm = TRUE, type = 7))
    }
    pos_thrL <- qfun(pos_vals, pmin); pos_thrU <- qfun(pos_vals, pmax)
    neg_thrL <- qfun(neg_vals_abs, pmin); neg_thrU <- qfun(neg_vals_abs, pmax)

    # Build a logical keep mask for upper-tri, then mirror to full matrix
    keep_ut <- logical(length(vals))
    if (length(pos_vals)) {
        kpos <- (vals >= 0) & (vals >= pos_thrL) & (vals <= pos_thrU)
        keep_ut <- keep_ut | kpos
    }
    if (length(neg_vals_abs)) {
        kneg <- (vals < 0) & (abs(vals) >= neg_thrL) & (abs(vals) <= neg_thrU)
        keep_ut <- keep_ut | kneg
    }
    if (!any(keep_ut)) return(mat * 0)

    # Initialize zero matrix and fill symmetric kept entries
    res <- mat * 0
    res[ut_mask][keep_ut] <- vals[keep_ut]
    res <- res + t(res)
    return(res)
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
.cluster_metrics <- function(members, cons_mat, W_mat) {
    cl_ids <- sort(na.omit(unique(members)))
    res <- lapply(cl_ids, function(cid) {
        idx <- which(members == cid)
        n <- length(idx)
        if (n <= 1) return(data.frame(cluster = cid, size = n, within_cons = NA_real_, conductance = NA_real_))
        C <- as.matrix(cons_mat[idx, idx])
        within_cons <- if (n > 1) mean(C[upper.tri(C)], na.rm = TRUE) else NA_real_
        Win <- sum(W_mat[idx, idx]) / 2
        Wcut <- sum(W_mat[idx, -idx, drop = FALSE])
        conductance <- Wcut / (Wcut + Win + 1e-12)
        data.frame(cluster = cid, size = n, within_cons = within_cons, conductance = conductance)
    })
    do.call(rbind, res)
}
#' @noRd
.gene_metrics <- function(members, cons_mat, W_mat) {
    cl_ids <- sort(na.omit(unique(members)))
    # precompute cluster membership indices
    idx_list <- lapply(cl_ids, function(cid) which(members == cid))
    names(idx_list) <- cl_ids
    n <- length(members)
    p_in <- numeric(n); p_best_out <- numeric(n); w_in <- numeric(n)
    for (i in seq_len(n)) {
        cid <- members[i]
        if (is.na(cid)) { p_in[i] <- NA; p_best_out[i] <- NA; w_in[i] <- NA; next }
        idx <- idx_list[[as.character(cid)]]
        idx_noi <- setdiff(idx, i)
        if (length(idx_noi)) {
            p_in[i] <- mean(cons_mat[i, idx_noi], na.rm = TRUE)
            w_in[i] <- sum(W_mat[i, idx_noi])
        } else {
            p_in[i] <- 0; w_in[i] <- 0
        }
        other <- setdiff(cl_ids, cid)
        if (length(other)) {
            means <- vapply(other, function(cj) {
                jj <- idx_list[[as.character(cj)]]
                if (length(jj)) mean(cons_mat[i, jj], na.rm = TRUE) else 0
            }, numeric(1))
            p_best_out[i] <- max(means, na.rm = TRUE)
        } else {
            p_best_out[i] <- 0
        }
    }
    data.frame(gene = kept_genes, cluster = members, p_in = p_in, p_best_out = p_best_out, w_in = w_in, stringsAsFactors = FALSE)
}
