#' @title Shared spatial utilities for scope creation pipelines
#' @description
#'   A lightweight collection of geometry helpers extracted from legacy
#'   `0.Helpers.r` so that Xenium, CosMx, and Visium loaders can share a
#'   consistent toolkit without duplicating logic.

#' Flip Y coordinates around a maximum bound.
#' Useful for mirroring image-based coordinate systems to mathematical axes.
#'
#' @param coords data.frame/data.table with a numeric column `y`.
#' @param y_max Numeric scalar; the reference maximum used for flipping.
#' @return The same object with `y` replaced by `y_max - y`.
#' @noRd
flip_y_coordinates <- function(coords, y_max) {
    stopifnot(is.numeric(y_max), length(y_max) == 1L)
    if (!"y" %in% names(coords)) {
        stop("Input data must expose a 'y' column.")
    }
    if (data.table::is.data.table(coords)) {
        coords[, y := y_max - y]
        coords
    } else {
        coords$y <- y_max - coords$y
        coords
    }
}

#' Clip point cloud to an sf/sfc region.
#' Keeps only spatial points that fall within the specified polygonal footprint.
#'
#' @param points data.frame/data.table with columns `x` and `y`.
#' @param region sf/sfc geometry describing the region of interest.
#' @param chunk_size Integer; number of rows per spatial join chunk.
#' @param workers Integer; parallel workers (mclapply on POSIX, snow on Windows).
#' @return A data.table containing only points inside the region.
#' @noRd
clip_points_to_region <- function(points,
                                  region,
                                  chunk_size = 5e5,
                                  workers = 1L) {
    if (is.null(region) || nrow(points) == 0L) {
        return(points)
    }
    if (inherits(region, "sf")) {
        region <- sf::st_geometry(region)
    }
    if (!inherits(region, "sfc")) {
        if (all(vapply(region, inherits, logical(1), "sfg"))) {
            region <- sf::st_sfc(region)
        } else {
            stop("`region` must be an sf/sfc object (or list of sfg).")
        }
    }
    if (length(region) == 0L) {
        return(points[0])
    }
    if (length(region) > 1L) {
        region <- sf::st_union(region)
    }
    if (length(region) == 0L) {
        return(points[0])
    }

    stopifnot(all(c("x", "y") %in% names(points)))
    dt_pts <- data.table::as.data.table(points)
    crs_use <- sf::st_crs(region) %||% NA
    chunk_ids <- split(seq_len(nrow(dt_pts)),
        ceiling(seq_len(nrow(dt_pts)) / chunk_size)
    )

    clip_worker <- function(idx) {
        subset_dt <- dt_pts[idx]
        inside_flag <- lengths(sf::st_within(
            sf::st_as_sf(subset_dt,
                coords = c("x", "y"),
                crs = crs_use, remove = FALSE
            ),
            region
        )) > 0
        subset_dt[inside_flag]
    }

    if (workers > 1L && .Platform$OS.type != "windows") {
        res_list <- parallel::mclapply(chunk_ids, clip_worker, mc.cores = workers)
    } else if (workers > 1L) {
        cluster <- parallel::makeCluster(workers)
        on.exit(parallel::stopCluster(cluster), add = TRUE)
        parallel::clusterEvalQ(cluster, library(sf))
        parallel::clusterExport(cluster,
            c("dt_pts", "region", "crs_use", "clip_worker"),
            envir = environment()
        )
        res_list <- parallel::parLapply(cluster, chunk_ids, clip_worker)
    } else {
        res_list <- lapply(chunk_ids, clip_worker)
    }
    data.table::rbindlist(res_list)
}

#' Read segmentation vertices and rebuild polygons.
#' Converts raw vertex listings into usable point clouds and boundary polygons.
#'
#' @param path Character; Parquet file containing segmentation vertices.
#' @param label Character tag written to output list names.
#' @param flip Logical; whether to flip Y coordinates using `y_max`.
#' @param y_max Numeric scalar; mandatory when `flip = TRUE`.
#' @param keep_ids Optional character vector of cell IDs to retain.
#' @param workers Integer; parallel workers during polygon construction.
#' @return A list with `points` (data.table) and `polygons` (sfc).
#' @noRd
build_segmentation_geometries <- function(path,
                                          label = "cell",
                                          flip = FALSE,
                                          y_max = NULL,
                                          keep_ids = NULL,
                                          workers = 1L) {
    if (!file.exists(path)) {
        stop("Segmentation file not found: ", path)
    }
    raw_dt <- data.table::as.data.table(arrow::read_parquet(path))
    if ("fov" %in% names(raw_dt)) {
        seg_dt <- raw_dt[, .(
            cell = paste0(as.character(cell_id), "_", as.character(fov)),
            x = vertex_x,
            y = vertex_y,
            label_id = paste0(as.character(cell_id), "_", as.character(fov))
        )]
    } else {
        seg_dt <- raw_dt[, .(
            cell = as.character(cell_id),
            x = vertex_x,
            y = vertex_y,
            label_id = as.character(label_id)
        )]
    }

    if (flip) {
        if (is.null(y_max)) stop("`y_max` must be supplied when `flip = TRUE`.")
        seg_dt[, y := y_max - y]
    }
    if (!is.null(keep_ids)) {
        seg_dt <- seg_dt[cell %in% keep_ids]
    }

    index_map <- split(seq_len(nrow(seg_dt)), seg_dt$label_id)
    polygon_builder <- function(idx) {
        sub_dt <- seg_dt[idx, .(x, y)]
        sub_dt <- sub_dt[complete.cases(sub_dt)]
        if (nrow(sub_dt) < 3L) return(NULL)
        coords <- as.matrix(sub_dt)
        storage.mode(coords) <- "double"
        if (!all(coords[1, ] == coords[nrow(coords), ])) {
            coords <- rbind(coords, coords[1, ])
        }
        tryCatch(sf::st_polygon(list(coords)), error = function(e) {
            message("[build_segmentation_geometries] Invalid polygon in ", label, ": ", conditionMessage(e))
            NULL
        })
    }

    if (workers > 1L && .Platform$OS.type != "windows") {
        poly_list <- parallel::mclapply(index_map, polygon_builder, mc.cores = workers)
    } else if (workers > 1L) {
        cluster <- parallel::makeCluster(workers)
        on.exit(parallel::stopCluster(cluster), add = TRUE)
        parallel::clusterEvalQ(cluster, library(sf))
        parallel::clusterExport(cluster, c("seg_dt", "polygon_builder"),
            envir = environment()
        )
        poly_list <- parallel::parLapply(cluster, index_map, polygon_builder)
    } else {
        poly_list <- lapply(index_map, polygon_builder)
    }
    list(points = seg_dt, polygons = sf::st_sfc(poly_list[!vapply(poly_list, is.null, logical(1L))]))
}
#' @title Helper Functions for geneSCOPE Package
#' @description Internal helper functions used throughout the package

#' Select a grid layer from the scope object, ensuring valid input.
#' Centralizes layer-naming checks so downstream helpers can assume consistent inputs.
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

#' Confirm that a grid layer contains the required data slots.
#' Guards against incomplete preprocessing before expensive computations run.
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

#' Decide which genes to use based on explicit lists or cluster labels.
#' Provides a single entry-point for deriving analysis subsets from flexible inputs.
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

#' Retrieve the Lee's L matrix for a grid layer, resolving default names.
#' Normalizes the various storage locations used across versions of geneSCOPE objects.
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
    # Note: do NOT intersect Lee's L with Pearson gene set here.
    # This keeps the full L matrix as provided by the LeeStats layer
    # and avoids unintended shrinkage of the gene set during Stage-1.

    Lmat
}

#' Fetch the Pearson correlation matrix and align it with Lee's L when possible.
#' Handles both grid-level and cell-level storage patterns with sensible fallbacks.
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

#' Summarize how often each edge co-clusters across runs.
#' Computes consensus weights on the provided edge list without allocating a dense N×N matrix.
#' @noRd
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

#' Transfer CMH or edge weights from a similarity graph to requested pairs.
#' Maps precomputed weights onto requested edges and falls back to pure R when C++ helpers are absent.
#' @noRd
.cmh_lookup_pairs <- function(kept_genes, ei, ej, g_sim, weight_col = c("CMH", "weight"), n_threads = NULL) {
    weight_col <- match.arg(weight_col)
    if (is.null(n_threads)) {
        # Robust threads selection across versions
        n_threads <- tryCatch({
            if (exists("safe_thread_count", mode = "function", inherits = TRUE)) {
                safe_thread_count()
            } else if (exists("getSafeThreadCount", mode = "function", inherits = TRUE)) {
                f <- get("getSafeThreadCount", mode = "function")
                an <- tryCatch(names(formals(f)), error = function(e) character(0))
                if ("max_requested" %in% an) f(max_requested = 8L)
                else if ("requested_cores" %in% an) f(requested_cores = 8L)
                else f(8L)
            } else {
                max(1L, min(8L, as.integer(parallel::detectCores()) - 1L))
            }
        }, error = function(e) 1L)
    }
    ed_sim <- igraph::as_data_frame(g_sim, what = "edges")
    wcol <- if (!is.null(ed_sim[[weight_col]])) weight_col else if (!is.null(ed_sim$weight)) "weight" else stop("CMH graph lacks weight column")
    si <- match(ed_sim$from, kept_genes)
    sj <- match(ed_sim$to, kept_genes)
    ok <- which(!is.na(si) & !is.na(sj))
    if (!length(ok)) {
        return(rep_len(1, length(ei)))
    }
    si <- as.integer(si[ok])
    sj <- as.integer(sj[ok])
    sw <- as.numeric(ed_sim[[wcol]][ok])
    fallback <- stats::median(sw, na.rm = TRUE)
    if (!is.finite(fallback)) fallback <- 1
    if (exists("cmh_lookup_rcpp", mode = "function")) {
        return(cmh_lookup_rcpp(as.integer(ei), as.integer(ej), si, sj, sw, fallback, as.integer(n_threads)))
    }
    # R fallback: direct keyed lookup to maintain length(ei)
    sim_names <- paste(
        pmin(kept_genes[si], kept_genes[sj]),
        pmax(kept_genes[si], kept_genes[sj]),
        sep = "|"
    )
    w_map <- sw
    names(w_map) <- sim_names
    g1 <- kept_genes[ei]
    g2 <- kept_genes[ej]
    edge_keys <- paste(pmin(g1, g2), pmax(g1, g2), sep = "|")
    ww <- w_map[edge_keys]
    ww[is.na(ww)] <- fallback
    as.numeric(ww)
}

#' Count intra-cluster degrees using existing edge endpoints only.
#' Tallies within-cluster degrees without materializing a full adjacency structure.
#' @noRd
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

#' Extract Pearson correlations for specific edge pairs without densifying.
#' Computes only the requested pairwise correlations, keeping large Pearson matrices in sparse form.
#' @noRd
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

#' Symmetrize a matrix via pmax while preserving sparse representations.
#' Provides a sparse-friendly alternative when `pmax` is unavailable or dense conversion is costly.
#' @noRd
.symmetric_pmax <- function(M) {
    if (inherits(M, "sparseMatrix")) {
        # Robust to Matrix versions without exported pmax(): use algebraic identity
        # max(A,B) = (A + B + |A - B|) / 2, then zero the diagonal
        M <- methods::as(M, "dgCMatrix")
        Mt <- Matrix::t(M)
        M <- (M + Mt + abs(M - Mt)) * 0.5
        Matrix::diag(M) <- 0
        return(Matrix::drop0(M))
    } else {
        M <- pmax(M, t(M))
        diag(M) <- 0
        return(M)
    }
}

#' Compute cluster-level quality metrics while keeping matrices sparse.
#' Summarizes within-cluster consensus and conductance without converting matrices to dense form.
#' @noRd
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

#' Derive per-gene metrics (p_in, w_in, etc.) using sparse matrices.
#' Quantifies intra- and inter-cluster affinities while avoiding dense memory usage.
#' @noRd
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

#' Run Leiden community detection through NetworKit with reticulate.
#' Provides a high-performance parallel backend while keeping a pure-R fallback.
#' @noRd
.community_detect_networkit <- function(graph, resolution = 1, threads = 1, debug = FALSE) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("reticulate package required for networkit backend", call. = FALSE)
    }

    edges <- igraph::as_data_frame(graph, what = "edges")
    n <- igraph::vcount(graph)
    node_names <- igraph::as_ids(igraph::V(graph))
    if (is.null(node_names)) node_names <- as.character(seq_len(n))

    if (!nrow(edges)) {
        memb <- rep(1L, n)
        names(memb) <- node_names
        return(memb)
    }

    if (!"weight" %in% names(edges)) edges$weight <- 1
    edges$weight[is.na(edges$weight)] <- 1

    idx_map <- setNames(seq_len(n) - 1L, node_names)
    edges$from_idx <- idx_map[as.character(edges$from)]
    edges$to_idx   <- idx_map[as.character(edges$to)]

    tmp <- tempfile(fileext = ".tsv")
    on.exit(unlink(tmp), add = TRUE)
    cols <- edges[, c("from_idx", "to_idx", "weight")]
    if (requireNamespace("data.table", quietly = TRUE)) {
        data.table::fwrite(cols, tmp, sep = " ", col.names = FALSE)
    } else {
        write.table(cols, file = tmp, sep = " ", col.names = FALSE,
                    row.names = FALSE, quote = FALSE)
    }

    nk <- reticulate::import("networkit", delay_load = TRUE)
    threads_use <- suppressWarnings(as.integer(threads))
    if (length(threads_use) == 0L || is.na(threads_use) || threads_use < 1L) {
        threads_use <- 1L
    }
    try(nk$setNumberOfThreads(threads_use), silent = TRUE)
    reader <- tryCatch(
        nk$graphio$EdgeListReader(" ", 0L, "#", TRUE, FALSE),
        error = function(e_new) {
            tryCatch(
                nk$graphio$EdgeListReader(weighted = TRUE, directed = FALSE),
                error = function(e_old) stop(e_old)
            )
        }
    )
    G <- reader$read(tmp)
    if (debug) message(sprintf("[networkit] ParallelLeiden threads=%d, resolution=%.4f",
                               threads_use, resolution))
    # NetworKit ≥11 uses `gamma` argument; older releases still expect `resolution`.
    create_algo <- function(use_gamma = TRUE) {
        res <- as.numeric(resolution)
        if (use_gamma) {
            nk$community$ParallelLeiden(G, gamma = res)
        } else {
            nk$community$ParallelLeiden(G, resolution = res)
        }
    }
    algo <- tryCatch(
        create_algo(TRUE),
        error = function(e_gamma) {
            msg <- conditionMessage(e_gamma)
            is_py_type_error <- inherits(e_gamma, "python.builtin.TypeError") ||
                inherits(e_gamma, "reticulate.python.builtin.type_error")
            has_gamma_kw <- grepl("unexpected keyword argument 'gamma'", msg, fixed = TRUE)
            has_cinit_arity <- grepl("__cinit__", msg, fixed = TRUE)
            if (is_py_type_error || has_gamma_kw || has_cinit_arity) {
                tryCatch(
                    create_algo(FALSE),
                    error = function(e_resolution) {
                        stop(e_resolution)
                    }
                )
            } else {
                stop(e_gamma)
            }
        }
    )
    algo$run()
    comm <- algo$getPartition()
    membership <- reticulate::py_to_r(comm$getVector())
    membership <- as.integer(membership) + 1L
    names(membership) <- node_names
    membership
}

#' Provide a backend-agnostic wrapper around Leiden/Louvain clustering.
#' Switches between igraph and NetworKit implementations with graceful fallbacks.
#' @noRd
.community_detect <- function(graph, algo = c("leiden", "louvain"), resolution = 1,
                               objective = c("CPM", "modularity"),
                               backend = c("igraph", "networkit"),
                               threads = 1, debug = FALSE) {
    algo <- match.arg(algo)
    objective <- match.arg(objective)
    backend <- match.arg(backend)

    if (algo == "leiden" && identical(backend, "networkit")) {
        memb <- try(.community_detect_networkit(graph, resolution = resolution,
                                               threads = threads, debug = debug), silent = TRUE)
        if (!inherits(memb, "try-error")) return(memb)
        if (debug) message("[networkit] Leiden failed, falling back to igraph implementation")
        backend <- "igraph"
    }

    if (identical(backend, "igraph")) {
        if (algo == "leiden") {
            comm <- try(igraph::cluster_leiden(graph, weights = igraph::E(graph)$weight,
                                               resolution = resolution, objective_function = objective), silent = TRUE)
            if (inherits(comm, "try-error")) {
                comm <- try(igraph::cluster_leiden(graph, weights = igraph::E(graph)$weight,
                                                   resolution_parameter = resolution, objective_function = objective), silent = TRUE)
            }
            if (inherits(comm, "try-error")) {
                comm <- try(igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight,
                                                     resolution = resolution, objective_function = objective), silent = TRUE)
                if (inherits(comm, "try-error")) {
                    comm <- try(igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight,
                                                         resolution_parameter = resolution, objective_function = objective), silent = TRUE)
                    if (inherits(comm, "try-error")) comm <- igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight)
                }
            }
            memb <- igraph::membership(comm)
            if (is.null(names(memb))) names(memb) <- igraph::as_ids(igraph::V(graph))
            return(memb)
        } else {
            comm <- try(igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight,
                                                 resolution = resolution, objective_function = objective), silent = TRUE)
            if (inherits(comm, "try-error")) {
                comm <- try(igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight,
                                                     resolution_parameter = resolution, objective_function = objective), silent = TRUE)
                if (inherits(comm, "try-error")) comm <- igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight)
            }
            memb <- igraph::membership(comm)
            if (is.null(names(memb))) names(memb) <- igraph::as_ids(igraph::V(graph))
            return(memb)
        }
    }

    stop("Unsupported backend", call. = FALSE)
}

#' Coerce diverse matrix-like inputs into base numeric matrices.
#' Standardizes sparse, bigmemory, and other exotic structures into plain numeric arrays.
#' @noRd
.ensure_numeric_matrix <- function(mat,
                                   nrow = NULL,
                                   ncol = NULL,
                                   rownames_out = NULL,
                                   colnames_out = NULL,
                                   as_sparse = FALSE) {
    if (is.null(mat)) return(NULL)

    coerce_sparse_pattern <- function(x) {
        idx <- Matrix::summary(x)
        rr <- idx$i
        cc <- idx$j
        vv <- idx$x
        if (is.null(vv)) vv <- rep(1, length(rr))
        nr <- nrow(x)
        nc <- ncol(x)
        dense <- matrix(0, nrow = nr, ncol = nc)
        dense[cbind(rr, cc)] <- vv
        dimnames(dense) <- dimnames(x)
        dense
    }

    out <- mat
    orig <- mat
    if (inherits(out, "big.matrix")) {
        if (requireNamespace("bigmemory", quietly = TRUE)) {
            out <- bigmemory::as.matrix(out)
        } else {
            stop("bigmemory package required to coerce big.matrix objects", call. = FALSE)
        }
    }
    if (inherits(out, "big.matrix")) {
        if (!requireNamespace("bigmemory", quietly = TRUE)) {
            stop("bigmemory package required to coerce big.matrix objects", call. = FALSE)
        }
        out <- bigmemory::as.matrix(out)
    }

    if (inherits(out, "Matrix")) {
        if (inherits(out, "nMatrix")) {
            out <- coerce_sparse_pattern(out)
        } else {
            out <- tryCatch(as.matrix(out), error = function(e) NULL)
            if (is.null(out)) {
                out <- coerce_sparse_pattern(methods::as(orig, "ngCMatrix"))
            }
        }
    } else if (!is.matrix(out)) {
        out <- tryCatch(as.matrix(out), error = function(e) NULL)
        if (is.null(out) && is.atomic(mat)) {
            if (is.null(nrow) || is.null(ncol)) stop("Cannot reshape atomic vector without target dimensions")
            out <- matrix(as.numeric(mat), nrow = nrow, ncol = ncol)
        } else if (is.null(out)) {
            dm <- dim(mat)
            if (!is.null(dm) && length(dm) == 2L) {
                out <- matrix(as.numeric(mat), nrow = dm[1], ncol = dm[2])
            }
        }
    }

    if (is.null(out)) {
        cls <- paste(class(orig), collapse = ", ")
        stop(sprintf("Failed to coerce object of class [%s] to numeric matrix", cls), call. = FALSE)
    }
    storage.mode(out) <- "double"
    if (!is.null(rownames_out)) rownames(out) <- rownames_out
    if (!is.null(colnames_out)) colnames(out) <- colnames_out
    if (as_sparse) {
        if (inherits(out, "sparseMatrix")) {
            out <- methods::as(out, "CsparseMatrix")
        } else {
            out <- methods::as(Matrix::Matrix(out, sparse = TRUE), "CsparseMatrix")
        }
        out <- Matrix::drop0(out)
    }
    out
}

#' Attempt to reshape an object safely to target dimensions.
#' Performs guarded `dim<-` assignment and falls back to matrix coercion if needed.
#' @noRd
.safe_set_dim <- function(M, target_dim) {
    if (is.null(target_dim)) return(M)
    ok <- tryCatch({
        dim(M) <- target_dim
        TRUE
    }, error = function(e) FALSE)
    if (ok) return(M)
    M_coerced <- try(as.matrix(M), silent = TRUE)
    if (inherits(M_coerced, "try-error")) {
        stop("Failed to coerce object to requested dimensions when aligning FDR matrix", call. = FALSE)
    }
    dim(M_coerced) <- target_dim
    M_coerced
}

#' Assign dimension names defensively, honoring debug options.
#' Applies `dimnames<-` but tolerates failures when debugging overrides are set.
#' @noRd
.safe_set_dimnames <- function(M, target_dn) {
    if (is.null(target_dn)) return(M)
    ok <- tryCatch({
        dimnames(M) <- target_dn
        TRUE
    }, error = function(e) FALSE)
    if (!ok && isTRUE(getOption("geneSCOPE.debug_dimnames", FALSE))) {
        message("[align_and_filter_FDR] dimnames assignment skipped; proceeding without named FDR matrix.")
    }
    M
}

#' Align an FDR matrix with L and zero out edges above threshold.
#' Keeps Lee statistics and FDR control matrices synchronized before masking edges.
#' @noRd
.align_and_filter_FDR <- function(L, L_raw, FDRmat, FDR_max) {
    if (is.null(FDRmat)) return(L)
    FM <- FDRmat

    # 1) Robust coercion to base matrix
    if (inherits(FM, "big.matrix")) {
        FM <- try({
            if (requireNamespace("bigmemory", quietly = TRUE)) bigmemory::as.matrix(FM) else FM[, ]
        }, silent = TRUE)
        if (inherits(FM, "try-error")) {
            warning("[align_and_filter_FDR] FDR is big.matrix but cannot coerce to base matrix; disabling FDR filter.")
            return(L)
        }
    } else if (!is.matrix(FM)) {
        FM_try <- try(as.matrix(FM), silent = TRUE)
        if (!inherits(FM_try, "try-error") && is.matrix(FM_try)) {
            FM <- FM_try
        } else {
            warning("[align_and_filter_FDR] FDR object cannot be coerced to matrix (class=", paste(class(FDRmat), collapse=","), "); disabling FDR filter.")
            return(L)
        }
    }

    # 2) Align by names (preferred) or by shape
    policy <- getOption("geneSCOPE.fdr_align_policy", "by_name")
    if (is.null(dim(FM))) {
        FM <- matrix(FM, nrow = nrow(L_raw), ncol = ncol(L_raw))
        FM <- .safe_set_dimnames(FM, dimnames(L_raw))
    } else if (!(identical(dim(FM), dim(L_raw)) &&
                 identical(rownames(FM), rownames(L_raw)) &&
                 identical(colnames(FM), colnames(L_raw)))) {
        if (identical(policy, "by_shape")) {
            FM <- .safe_set_dim(FM, dim(L_raw))
            FM <- .safe_set_dimnames(FM, dimnames(L_raw))
        } else {
            if (!is.null(rownames(FM)) && !is.null(colnames(FM)) &&
                all(rownames(L_raw) %in% rownames(FM)) && all(colnames(L_raw) %in% colnames(FM))) {
                FM <- FM[rownames(L_raw), colnames(L_raw), drop = FALSE]
            } else {
                # fall back to shape-only alignment, but keep numeric matrix where possible
                FM <- .safe_set_dim(FM, dim(L_raw))
                FM <- .safe_set_dimnames(FM, dimnames(L_raw))
            }
        }
    }
    FM <- .safe_set_dimnames(FM, dimnames(L_raw))

    # 3) Subset FDR to current kept rows of L
    target_rows <- rownames(L)
    if (!is.null(target_rows) && !is.null(rownames(FM))) {
        idx <- match(target_rows, rownames(FM))
        if (anyNA(idx)) {
            stop("FDR matrix cannot be aligned to current gene subset", call. = FALSE)
        }
        FM <- FM[idx, idx, drop = FALSE]
    } else if (nrow(FM) != nrow(L) || ncol(FM) != ncol(L)) {
        stop("FDR matrix dimensions do not match current subset", call. = FALSE)
    }

    # 4) Apply mask
    if (inherits(L, "sparseMatrix")) {
        LT <- methods::as(L, "TsparseMatrix")
        if (length(LT@x)) {
            rows <- LT@i + 1L
            cols <- LT@j + 1L
            mask <- (FM[cbind(rows, cols)] > FDR_max)
            if (any(mask)) LT@x[mask] <- 0
        }
        return(Matrix::drop0(methods::as(LT, "CsparseMatrix")))
    }
    L[FM > FDR_max] <- 0
    L
}

#' Flip y-coordinates in-place for data frames or data.tables.
#' Provides a lightweight mirror transform while preserving data.table semantics.
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

#' Clip point coordinates to an sf polygon in chunks for memory safety.
#' Processes large point clouds in parallel-friendly batches to avoid memory spikes.
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

#' Convert segmentation vertices into point tables and polygons.
#' Handles optional axis flips and ID filtering so downstream workflows see clean geometry.
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

#' Parse a quantile shorthand (e.g., q95) into a probability.
#' Supports the qXX notation used throughout geneSCOPE configuration parameters.
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

#' Retain matrix entries that fall within chosen quantile bounds.
#' Applies percentile-based masking to both dense and sparse matrices with zero-awareness.
#' @noRd
.filter_matrix_by_quantile <- function(mat, pct_min = "q0", pct_max = "q100") {
    stopifnot(is.matrix(mat) || inherits(mat, "Matrix"))
    pmin <- .parse_q(pct_min)
    pmax <- .parse_q(pct_max)
    if (pmin > pmax) stop("pct_min > pct_max")

    legacy_mode <- !identical(getOption("geneSCOPE.filter_sparse_mode", "legacy"), "sparse")

    include_zero_mass <- !identical(getOption("geneSCOPE.filter_sparse_ignore_zero", FALSE), TRUE)

    adjust_quantile_with_zeros <- function(vec, p, zero_frac) {
        if (!length(vec)) return(NA_real_)
        if (!is.finite(zero_frac) || zero_frac <= 0) {
            return(as.numeric(stats::quantile(vec, p, na.rm = TRUE, type = 7)))
        }
        if (zero_frac >= 1) return(0)
        if (p <= zero_frac) return(0)
        adj_p <- (p - zero_frac) / (1 - zero_frac)
        as.numeric(stats::quantile(vec, adj_p, na.rm = TRUE, type = 7))
    }

    if (legacy_mode) {
        dense_mat <- as.matrix(mat)
        vec <- as.vector(dense_mat)
        pos_mask <- vec >= 0
        neg_mask <- vec < 0

        res <- dense_mat * 0
        if (any(pos_mask, na.rm = TRUE)) {
            pos_vec <- vec[pos_mask]
            thrL <- stats::quantile(pos_vec, pmin, na.rm = TRUE, type = 7)
            thrU <- stats::quantile(pos_vec, pmax, na.rm = TRUE, type = 7)
            keep <- pos_mask & vec >= thrL & vec <= thrU
            res[keep] <- dense_mat[keep]
        }
        if (any(neg_mask, na.rm = TRUE)) {
            neg_vec <- abs(vec[neg_mask])
            thrL <- stats::quantile(neg_vec, pmin, na.rm = TRUE, type = 7)
            thrU <- stats::quantile(neg_vec, pmax, na.rm = TRUE, type = 7)
            keep <- neg_mask & abs(vec) >= thrL & abs(vec) <= thrU
            res[keep] <- dense_mat[keep]
        }
        return(Matrix::drop0(Matrix::Matrix(res, sparse = TRUE)))
    }

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
                if (include_zero_mass) {
                    n_nodes <- nrow(mat)
                    total_pairs <- (as.double(n_nodes) * (n_nodes - 1)) / 2
                    zero_count <- max(total_pairs - length(pos_vec), 0)
                    zero_frac <- if (total_pairs > 0) zero_count / (zero_count + length(pos_vec)) else 0
                    pos_thrL <- adjust_quantile_with_zeros(pos_vec, pmin, zero_frac)
                    pos_thrU <- adjust_quantile_with_zeros(pos_vec, pmax, zero_frac)
                } else {
                    pos_thrL <- as.numeric(stats::quantile(pos_vec, pmin, na.rm = TRUE, type = 7))
                    pos_thrU <- as.numeric(stats::quantile(pos_vec, pmax, na.rm = TRUE, type = 7))
                }
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

#' Compute block identifiers from grid coordinates and block size.
#' Groups grid tiles into coarse blocks so block-aware statistics can be aggregated quickly.
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

#' Summarize cluster metrics using dense matrices for smaller cases.
#' Provides a simple dense fallback when sparse machinery is unnecessary.
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
#' Produce per-gene metrics when dense consensus/weight matrices are available.
#' Reports intra-cluster strength and the strongest external attraction per gene.
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
