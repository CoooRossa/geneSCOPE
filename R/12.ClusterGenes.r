#' @title Gene Clustering (results written to scope_obj@meta.data)
#' @description Based on Lee's L statistics for filtering, graph construction, multiple community detection and optional consensus, output cluster numbers to scope_obj@meta.data.
#' @export
#' @noRd 
# Lightweight worker: do NOT capture parent environment
# Accepts 'ii' and '...' to swallow any future-specific extras
.geneSCOPE_cluster_once <- function(ii, graph, algo = c("leiden", "louvain"), resolution = 1, objective = c("CPM", "modularity"), ...) {
    # Guard against thread oversubscription inside workers
    try({ if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1) }, silent = TRUE)
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
    try({ if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1) }, silent = TRUE)
    algo <- match.arg(algo)
    objective <- match.arg(objective)
    if (algo == "leiden") {
        comm <- try(igraph::cluster_leiden(graph, weights = igraph::E(graph)$weight, resolution = resolution, objective_function = objective), silent = TRUE)
        if (inherits(comm, "try-error")) comm <- igraph::cluster_leiden(graph, weights = igraph::E(graph)$weight, resolution_parameter = resolution, objective_function = objective)
    } else {
        comm <- try(igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight, resolution = resolution, objective_function = objective), silent = TRUE)
        if (inherits(comm, "try-error")) {
            comm <- try(igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight, resolution_parameter = resolution, objective_function = objective), silent = TRUE)
            if (inherits(comm, "try-error")) comm <- igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight)
        }
    }
    as.integer(igraph::membership(comm))
}
#' @noRd 
# Batch worker to avoid capturing large environments; only depends on arguments.
.geneSCOPE_cluster_batch <- function(batch_ids, graph, algo = c("leiden", "louvain"), resolution = 1, objective = c("CPM", "modularity")) {
    res <- vector("list", length(batch_ids))
    for (k in seq_along(batch_ids)) {
        # Use the single-run worker to keep identical code paths and RNG behavior under future.seed
        res[[k]] <- .geneSCOPE_cluster_once(batch_ids[[k]], graph = graph, algo = algo, resolution = resolution, objective = objective)
    }
    res
}
#' @export
clusterGenes <- function(
    scope_obj,
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    L_min = 0,
    use_FDR = TRUE,
    FDR_max = 0.05,
    pct_min = "q0",
    drop_isolated = TRUE,
    algo = c("leiden", "louvain"),
    resolution = 1,
    objective = c("CPM", "modularity"),
    cluster_name = NULL,
    use_log1p_weight = TRUE,
    use_consensus = TRUE,
    graph_slot_name = "g_consensus",
    n_restart = 100,
    consensus_thr = 0,
    keep_cross_stable = TRUE,
    CI95_filter = FALSE,
    curve_layer = "LR_curve",
    CI_rule = c("remove_within", "remove_outside"),
    use_cmh_weight = FALSE,
    cmh_slot = "g_morisita",
    # smoothing for corrected weights
    post_smooth = TRUE,
    post_smooth_quant = c(0.05, 0.95),
    post_smooth_power = 0.5,
    # conservative subclustering and QC controls
    enable_subcluster = TRUE,
    sub_min_size = 10,
    sub_min_child_size = 5,
    sub_resolution_factor = 1.3,
    sub_within_cons_max = 0.5,
    sub_conductance_min = 0.6,
    sub_improve_within_cons_min = 0.07,
    sub_max_groups = 3,
    enable_qc_filter = TRUE,
    qc_gene_intra_cons_min = 0.25,
    qc_gene_best_out_cons_min = 0.6,
    qc_gene_intra_weight_q = 0.05,
    # preserve connectivity inside stage-1 clusters after corrections
    keep_stage1_backbone = TRUE,
    backbone_floor_q = 0.02,
    # optional report output
    return_report = FALSE,
    verbose = TRUE,
    debug = FALSE,
    ncores = 16,
    future_strategy = c("auto", "multisession", "multicore", "inherit", "sequential"),
    # parallel controls
    future_scheduling_stage1 = 1,
    future_scheduling_stage2 = 0.5,
    restart_batch_size = NULL,
    aggressive_cpu = FALSE,
    stage2_cluster_parallel = FALSE,
    stage2_max_workers = NULL,
    mode = c("fast", "robust"),            # speed/memory profile
    adaptive_restarts = FALSE,               # do not cap restarts by default
    target_pair_ops = 2.5e8,                # total pairwise ops budget for consensus (fast mode)
    aggressive_gc = TRUE                    # free big objects between stages
    ) {

    algo <- match.arg(algo)
    objective <- match.arg(objective)
    CI_rule <- match.arg(CI_rule)
    mode <- match.arg(mode)
    fast_mode <- identical(mode, "fast")

    ## 1. Select grid layer
    g_layer <- .selectGridLayer(scope_obj, grid_name)
    gname <- names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    .checkGridContent(scope_obj, gname)

    ## 1.1 Configure future plan locally (multisession on Windows/RStudio, multicore otherwise)
    old_plan <- NULL
    if (requireNamespace("future", quietly = TRUE)) {
        strategy <- match.arg(future_strategy)
        os_type <- tryCatch(detectOS(), error = function(e) if (.Platform$OS.type == "windows") "windows" else "linux")
        # Save and restore current plan only if we change it
        set_plan <- function(kind) {
            if (is.null(old_plan)) old_plan <<- tryCatch(future::plan(), error = function(e) NULL)
            on.exit({ if (!is.null(old_plan)) try(future::plan(old_plan), silent = TRUE) }, add = TRUE)
            switch(kind,
                multisession = future::plan(future::multisession, workers = ncores),
                multicore    = future::plan(future::multicore,    workers = ncores),
                sequential   = future::plan(future::sequential),
                NULL
            )
        }

        if (identical(strategy, "inherit")) {
            # do not touch user's current plan
        } else if (identical(strategy, "multisession")) {
            set_plan("multisession")
        } else if (identical(strategy, "multicore")) {
            set_plan("multicore")
        } else if (identical(strategy, "sequential")) {
            set_plan("sequential")
        } else { # auto
            # Linux/macOS -> multicore (CoW), Windows/RStudio -> multisession
            is_rstudio <- tryCatch(nzchar(Sys.getenv("RSTUDIO")), error = function(e) FALSE)
            if (.Platform$OS.type == "windows" || is_rstudio) {
                set_plan("multisession")
            } else {
                set_plan("multicore")
            }
        }
    }

    ## 1.2 Configure threads using package helpers (zzz.R)
    # Follow computeL (6.LeesL.r) pattern: mixed task with restore
    thread_cfg <- try(configureThreadsFor(task_type = "mixed", ncores_requested = ncores, restore_after = TRUE), silent = TRUE)
    if (!inherits(thread_cfg, "try-error")) {
        on.exit({
            restore_fn <- attr(thread_cfg, "restore_function")
            if (!is.null(restore_fn)) try(restore_fn(), silent = TRUE)
        }, add = TRUE)
    }

    ## 2. Read statistics (robust: use helper for L; avoid `$` on possibly atomic)
    if (is.null(scope_obj@stats[[gname]]) || is.null(scope_obj@stats[[gname]][[lee_stats_layer]])) {
        stop("Statistics layer missing: ", lee_stats_layer)
    }
    leeStat <- scope_obj@stats[[gname]][[lee_stats_layer]]
    # Use helper to safely fetch L; it tolerates @stats/@grid placement and checks validity
    L_raw <- .getLeeMatrix(scope_obj, grid_name = gname, lee_layer = lee_stats_layer)
    if (is.null(L_raw)) stop("LeeStats layer missing L matrix")
    genes_all <- rownames(L_raw)

    FDRmat <- if (use_FDR) {
        FDR_candidate <- NULL
        if (!is.null(leeStat) && !is.null(leeStat[["FDR"]])) {
            FDR_candidate <- leeStat[["FDR"]]
        }
        if (is.null(FDR_candidate)) stop("use_FDR=TRUE but FDR matrix not found")
        FDR_candidate
    } else {
        matrix(0, nrow(L_raw), ncol(L_raw))
    }

    ## 3. Threshold filtering (base, before any correction)
    A <- L_raw
    A[abs(A) < L_min] <- 0
    if (use_FDR) {
        FM2 <- FDRmat
        # Coerce to matrix if needed; ensure dimensions align with L_raw (avoid dimnames<- on non-arrays)
        if (!is.matrix(FM2)) {
            FM2_try <- try(as.matrix(FM2), silent = TRUE)
            if (!inherits(FM2_try, "try-error")) FM2 <- FM2_try
        }
        if (is.null(dim(FM2))) {
            FM2 <- matrix(FM2, nrow = nrow(L_raw), ncol = ncol(L_raw))
        } else if (!identical(dim(FM2), dim(L_raw))) {
            # Try name-based reindex; otherwise reshape to L_raw dims
            if (!is.null(rownames(FM2)) && !is.null(colnames(FM2)) &&
                all(rownames(L_raw) %in% rownames(FM2)) &&
                all(colnames(L_raw) %in% colnames(FM2))) {
                FM2 <- FM2[rownames(L_raw), colnames(L_raw), drop = FALSE]
            } else {
                FM2 <- matrix(as.vector(FM2), nrow = nrow(L_raw), ncol = ncol(L_raw))
            }
        }
        A[FM2 > FDR_max] <- 0
    }
    A[A < 0] <- 0
    A <- .symmetric_pmax(A)
    A <- .filter_matrix_by_quantile(A, pct_min, "q100")
    # enforce sparse to avoid dense O(n^2) scans downstream
    if (!inherits(A, "sparseMatrix")) A <- Matrix::drop0(Matrix::Matrix(A, sparse = TRUE)) else A <- Matrix::drop0(A)

    ## 4. Base graph, without CI95/CMH corrections
    keep <- if (drop_isolated) Matrix::rowSums(A != 0) > 0 else rep(TRUE, nrow(A))
    kept_genes <- genes_all[keep]
    if (length(kept_genes) < 2) stop("Fewer than 2 genes after filtering")

    A_sub <- A[kept_genes, kept_genes, drop = FALSE]
    # Robust triplet extraction (avoid summary/as.data.frame dimnames pitfalls)
    if (inherits(A_sub, "sparseMatrix")) {
        TT <- methods::as(A_sub, "TsparseMatrix")
        xvals <- tryCatch(TT@x, error = function(e) rep(1, length(TT@i)))
        ed0 <- data.frame(i = TT@i + 1L, j = TT@j + 1L, x = xvals, stringsAsFactors = FALSE)
    } else {
        idx <- which(A_sub > 0, arr.ind = TRUE)
        ed0 <- data.frame(i = idx[,1], j = idx[,2], x = A_sub[idx], stringsAsFactors = FALSE)
    }
    ed0 <- ed0[ed0$i < ed0$j & ed0$x > 0, ]
    if (!nrow(ed0)) stop("All edges filtered out")

    edges_df <- data.frame(
        from = kept_genes[ed0$i],
        to = kept_genes[ed0$j],
        L_raw_edge = ed0$x,
        stringsAsFactors = FALSE
    )
    edges_df$weight <- if (use_log1p_weight) log1p(edges_df$L_raw_edge) else edges_df$L_raw_edge
    g0 <- igraph::graph_from_data_frame(edges_df[, c("from", "to", "weight")],
        directed = FALSE, vertices = kept_genes
    )

    ## 5. First-round multi-run community detection / consensus
    # Note: avoid local wrapper closures that capture large environments; call igraph API inline with try/ fallback

    old_max_size <- getOption("future.globals.maxSize")
    # In forked (multicore) futures, globals are CoW; avoid false positives by lifting the cap.
    is_rstudio <- tryCatch(nzchar(Sys.getenv("RSTUDIO")), error = function(e) FALSE)
    is_windows <- identical(.Platform$OS.type, "windows")
    if (!is_windows && !is_rstudio) {
        options(future.globals.maxSize = +Inf)
    } else {
        # multisession: keep a conservative floor to reduce accidental huge transfers
        min_max_size <- 2 * 1024^3
        new_max_size <- old_max_size
        if (is.null(new_max_size) || (!is.infinite(new_max_size) && new_max_size < min_max_size)) {
            new_max_size <- min_max_size
        }
        options(future.globals.maxSize = new_max_size)
    }
    on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)

    # graph is constructed once and reused across restarts; avoid per-iteration rebuilds

    # Adaptive restart capping in fast mode to bound n^2 * restarts cost
    n_restart_eff <- n_restart
    if (fast_mode && adaptive_restarts && use_consensus) {
        ng <- length(kept_genes)
        pair_per_run <- (as.double(ng) * (ng - 1)) / 2
        if (is.finite(pair_per_run) && pair_per_run > 0) {
            cap <- max(50L, as.integer(floor(target_pair_ops / pair_per_run)))
            if (cap < n_restart_eff) {
                if (verbose) message("[geneSCOPE::clusterGenes] Adaptive restart cap: ", n_restart_eff, " -> ", cap)
                n_restart_eff <- cap
            }
        }
    }

    # Stage-1: restart batching to reduce scheduling/serialization overhead
    rvec <- seq_len(n_restart_eff)
    workers <- tryCatch(future::nbrOfWorkers(), error = function(e) NA_integer_)
    if (is.na(workers) || workers < 1) workers <- ncores
    if (is.null(restart_batch_size)) {
        bs <- if (n_restart_eff >= 2 * workers) max(3L, min(20L, as.integer(ceiling(n_restart_eff / workers)))) else 1L
        if (isTRUE(aggressive_cpu)) bs <- max(3L, min(10L, bs))
    } else bs <- as.integer(max(1L, restart_batch_size))
    batches_stage1 <- split(rvec, ceiling(rvec / bs))
    if (isTRUE(verbose) || isTRUE(debug)) {
        message(sprintf("[Stage1] workers=%d, n_restart=%d, batch_size=%d, n_batches=%d, scheduling=%.2f",
                        workers, n_restart_eff, bs, length(batches_stage1), as.numeric(future_scheduling_stage1)))
    }
    memb_list1 <- future.apply::future_lapply(
        X = batches_stage1,
        FUN = .geneSCOPE_cluster_batch,
        graph = g0,
        algo = algo,
        resolution = resolution,
        objective = objective,
        future.seed = TRUE,
        future.globals = FALSE,
        future.packages = c("igraph", "Matrix", "data.table"),
        future.scheduling = future_scheduling_stage1
    )
    memb_list1 <- unlist(memb_list1, recursive = FALSE, use.names = FALSE)
    # Robust bind to matrix
    memb_mat1 <- try(do.call(cbind, memb_list1), silent = TRUE)
    if (inherits(memb_mat1, "try-error") || is.null(dim(memb_mat1))) {
        memb_mat1 <- matrix(unlist(memb_list1, use.names = FALSE), nrow = length(kept_genes))
    }
    rownames(memb_mat1) <- kept_genes

    if (!use_consensus) {
        memb_final1 <- memb_mat1[, 1]
        cons <- Matrix::sparseMatrix(dims = c(length(kept_genes), length(kept_genes)))
    } else {
        # Build consensus only on base edges (vectorized over restart blocks)
        E0 <- nrow(ed0)
        if (isTRUE(verbose) || isTRUE(debug)) {
            message(sprintf("[Stage1] consensus inputs: edges(E0)=%d, restarts(R)=%d, threads=%d, cpp_available=%s",
                            E0, ncol(memb_mat1), ncores, as.character(exists("consensus_on_edges_omp"))))
        }
        cons_w <- .consensus_on_edges(memb_mat1, ed0$i, ed0$j, n_threads = ncores)
        keep_e <- cons_w >= consensus_thr & is.finite(cons_w)
        if (any(keep_e)) {
            g_cons_tmp1 <- igraph::graph_from_data_frame(
                data.frame(from = kept_genes[ed0$i[keep_e]], to = kept_genes[ed0$j[keep_e]], weight = cons_w[keep_e]),
                directed = FALSE, vertices = kept_genes
            )
        } else {
            g_cons_tmp1 <- g0
        }
        # community on consensus graph with API fallback
        comm1 <- try(igraph::cluster_leiden(g_cons_tmp1, weights = igraph::E(g_cons_tmp1)$weight, resolution = resolution, objective_function = objective), silent = TRUE)
        if (inherits(comm1, "try-error")) {
            comm1 <- try(igraph::cluster_leiden(g_cons_tmp1, weights = igraph::E(g_cons_tmp1)$weight, resolution_parameter = resolution, objective_function = objective), silent = TRUE)
        }
        if (inherits(comm1, "try-error")) {
            comm1 <- try(igraph::cluster_louvain(g_cons_tmp1, weights = igraph::E(g_cons_tmp1)$weight, resolution = resolution, objective_function = objective), silent = TRUE)
            if (inherits(comm1, "try-error")) {
                comm1 <- try(igraph::cluster_louvain(g_cons_tmp1, weights = igraph::E(g_cons_tmp1)$weight, resolution_parameter = resolution, objective_function = objective), silent = TRUE)
                if (inherits(comm1, "try-error")) comm1 <- igraph::cluster_louvain(g_cons_tmp1, weights = igraph::E(g_cons_tmp1)$weight)
            }
        }
        memb_final1 <- igraph::membership(comm1)

        # sparse symmetric consensus matrix (only kept edges)
        if (any(keep_e)) {
            cons <- Matrix::sparseMatrix(i = ed0$i[keep_e], j = ed0$j[keep_e], x = cons_w[keep_e],
                                         dims = c(length(kept_genes), length(kept_genes)))
            cons <- cons + Matrix::t(cons)
        } else {
            cons <- Matrix::sparseMatrix(dims = c(length(kept_genes), length(kept_genes)))
        }
    }

    ## Early exit: if no post-correction requested, return stage-1 result
    if (!use_cmh_weight && !CI95_filter) {
        memb_final <- memb_final1
        # use sparse consensus built above (may be empty)

        # Clean isolated / micro clusters using base adjacency edges only (O(E))
        if (inherits(A_sub, "sparseMatrix")) {
            TT_tmp <- methods::as(A_sub[kept_genes, kept_genes, drop = FALSE], "TsparseMatrix")
            msk <- (TT_tmp@x > 0) & (TT_tmp@i < TT_tmp@j)
            if (any(msk)) { ei <- TT_tmp@i[msk] + 1L; ej <- TT_tmp@j[msk] + 1L } else { ei <- ej <- integer(0) }
        } else {
            idx_e <- which(upper.tri(A_sub[kept_genes, kept_genes, drop = FALSE]) &
                               A_sub[kept_genes, kept_genes, drop = FALSE] > 0, arr.ind = TRUE)
            ei <- idx_e[, 1]; ej <- idx_e[, 2]
        }
        di <- .intra_degree_from_edges(ei, ej, memb_final)
        memb_final[di == 0] <- NA
        tbl <- table(memb_final, useNA = "no")
        small <- names(tbl)[tbl < 2]
        memb_final[memb_final %in% small] <- NA
        stable_cl <- names(tbl)[tbl >= 2]
        names(memb_final) <- kept_genes

        # Build consensus graph (sparse upper-tri only)
        if (inherits(cons, "sparseMatrix")) {
            TTc <- methods::as(cons, "TsparseMatrix")
            maskc <- (TTc@x > 0) & (TTc@i < TTc@j)
            if (any(maskc)) {
                edges_keep <- cbind(TTc@i[maskc] + 1L, TTc@j[maskc] + 1L)
            } else {
                edges_keep <- matrix(numeric(0), ncol = 2)
            }
        } else {
            edges_keep <- which(cons > 0, arr.ind = TRUE)
            edges_keep <- edges_keep[edges_keep[, 1] < edges_keep[, 2], ]
        }
        intra_df <- data.frame(
            from = kept_genes[edges_keep[, 1]],
            to = kept_genes[edges_keep[, 2]],
            weight = A_sub[cbind(edges_keep[, 1], edges_keep[, 2])],
            sign = ifelse(L_raw[cbind(
                match(kept_genes[edges_keep[, 1]], rownames(L_raw)),
                match(kept_genes[edges_keep[, 2]], colnames(L_raw))
            )] < 0, "neg", "pos"),
            stringsAsFactors = FALSE
        )

        cross_df <- NULL
        if (keep_cross_stable && length(stable_cl) > 1) {
            idx_stable <- memb_final %in% stable_cl
            genes_stable <- kept_genes[idx_stable]
            L_sub_raw <- L_raw[genes_stable, genes_stable]
            if (inherits(L_sub_raw, "sparseMatrix")) {
                TTls <- methods::as(L_sub_raw, "TsparseMatrix")
                maskls <- (TTls@x >= L_min) & (TTls@i < TTls@j)
                if (any(maskls)) {
                    ij <- cbind(TTls@i[maskls] + 1L, TTls@j[maskls] + 1L)
                } else {
                    ij <- matrix(numeric(0), ncol = 2)
                }
            } else {
                ij <- which(L_sub_raw >= L_min & upper.tri(L_sub_raw), arr.ind = TRUE)
            }
            if (nrow(ij)) {
                g1 <- genes_stable[ij[, 1]]
                g2 <- genes_stable[ij[, 2]]
                cl1 <- memb_final[g1]
                cl2 <- memb_final[g2]
                keep_cross <- cl1 != cl2 & !is.na(cl1) & !is.na(cl2)
                if (any(keep_cross)) {
                    g1 <- g1[keep_cross]
                    g2 <- g2[keep_cross]
                    Lval <- L_sub_raw[cbind(match(g1, genes_stable), match(g2, genes_stable))]
                    use_log1p_local <- use_log1p_weight
                    wL <- if (use_log1p_local) log1p(Lval) else Lval
                    wgt <- wL
                    mask <- wgt > 0
                    if (any(mask)) {
                        cross_df <- data.frame(
                            from = g1[mask], to = g2[mask],
                            weight = wgt[mask],
                            sign = ifelse(Lval[mask] < 0, "neg", "pos")
                        )
                    }
                }
            }
        }

        edges_df_final <- if (!is.null(cross_df)) unique(rbind(intra_df, cross_df)) else intra_df
        g_cons <- igraph::graph_from_data_frame(edges_df_final,
            directed = FALSE,
            vertices = kept_genes
        )

        # Continuous numbering
        if (length(memb_final)) {
            valid_idx <- !is.na(memb_final)
            old_ids <- sort(unique(memb_final[valid_idx]))
            id_map <- setNames(seq_along(old_ids), old_ids)
            memb_final[valid_idx] <- id_map[as.character(memb_final[valid_idx])]
        }

        # Write to meta.data
        if (is.null(cluster_name) || !nzchar(cluster_name)) {
            cluster_name <- sprintf("modL%.2f", L_min)
        }
        if (is.null(scope_obj@meta.data)) {
            scope_obj@meta.data <- data.frame(row.names = genes_all)
        } else {
            miss_rows <- setdiff(genes_all, rownames(scope_obj@meta.data))
            if (length(miss_rows)) {
                scope_obj@meta.data <- rbind(scope_obj@meta.data, data.frame(row.names = miss_rows))
            }
        }
        scope_obj@meta.data[genes_all, cluster_name] <- NA_integer_
        scope_obj@meta.data[kept_genes, cluster_name] <- as.integer(memb_final[kept_genes])

        if (verbose) {
            nz <- sum(!is.na(scope_obj@meta.data[genes_all, cluster_name]))
            message("[geneSCOPE::clusterGenes] (stage1 only) Written clustering column '", cluster_name, "': ",
                    nz, "/", length(genes_all), " genes (",
                    sprintf("%.1f%%", 100 * nz / length(genes_all)), ")")
        }

        # Save consensus graph
        if (is.null(scope_obj@stats[[gname]][[lee_stats_layer]])) {
            scope_obj@stats[[gname]][[lee_stats_layer]] <- list()
        }
        scope_obj@stats[[gname]][[lee_stats_layer]][[graph_slot_name]] <- g_cons

        # Return
        if (isTRUE(return_report)) {
            report <- new.env(parent = emptyenv())
            report$note <- "clusterGenes early-exit: no CI95/CMH correction requested; returned stage1 result."
            return(list(scope_obj = scope_obj, report = report))
        } else {
            return(scope_obj)
        }
    }

    ## Optional: free large Stage-1 objects before post-corrections
    if (isTRUE(aggressive_gc)) {
        to_drop <- intersect(ls(), c("memb_list1", "memb_mat1", "g0", "edges_df", "edge_list", "edge_weights", "ed0", "cons"))
        if (length(to_drop)) rm(list = to_drop)
        gc(FALSE)
    }

    ## 6. Correct L by CI95 and/or CMH (post step)
    # start from raw L among kept genes
    L_post <- L_raw[kept_genes, kept_genes, drop = FALSE]
    # apply baseline thresholding again to keep consistency
    L_post[abs(L_post) < L_min] <- 0
    if (use_FDR) {
        # robust subsetting regardless of FDR dimnames presence; avoid dimnames<- on non-arrays
        FM <- FDRmat
        if (!is.matrix(FM)) {
            FM_try <- try(as.matrix(FM), silent = TRUE)
            if (!inherits(FM_try, "try-error")) FM <- FM_try
        }
        if (is.null(dim(FM))) {
            FM <- matrix(FM, nrow = nrow(L_raw), ncol = ncol(L_raw))
        } else if (!identical(dim(FM), dim(L_raw))) {
            if (!is.null(rownames(FM)) && !is.null(colnames(FM)) &&
                all(rownames(L_raw) %in% rownames(FM)) && all(colnames(L_raw) %in% colnames(FM))) {
                FM <- FM[rownames(L_raw), colnames(L_raw), drop = FALSE]
            } else {
                FM <- matrix(as.vector(FM), nrow = nrow(L_raw), ncol = ncol(L_raw))
            }
        }
        ridx <- match(kept_genes, rownames(L_raw))
        FDR_sub <- FM[ridx, ridx, drop = FALSE]
        L_post[FDR_sub > FDR_max] <- 0
    }
    L_post[L_post < 0] <- 0
    L_post <- .symmetric_pmax(L_post)
    L_post <- .filter_matrix_by_quantile(L_post, pct_min, "q100")
    if (!inherits(L_post, "sparseMatrix")) L_post <- Matrix::drop0(Matrix::Matrix(L_post, sparse = TRUE)) else L_post <- Matrix::drop0(L_post)

    # CI95 masking on L (if requested)
    if (CI95_filter) {
        curve_obj <- leeStat[[curve_layer]]
        if (is.null(curve_obj)) stop("CI95_filter=TRUE but not found: ", curve_layer)
        curve_mat <- if (inherits(curve_obj, "big.matrix")) {
            bigmemory::as.matrix(curve_obj)
        } else if (is.matrix(curve_obj)) curve_obj else as.matrix(as.data.frame(curve_obj))
        xp <- as.numeric(curve_mat[, "Pear"])
        lo <- as.numeric(curve_mat[, "lo95"])
        hi <- as.numeric(curve_mat[, "hi95"])
        if (anyNA(c(xp, lo, hi))) stop("LR curve contains NA")
        f_lo <- approxfun(xp, lo, rule = 2)
        f_hi <- approxfun(xp, hi, rule = 2)

        # Only evaluate Pearson on existing positive edges (upper triangle)
        if (inherits(L_post, "sparseMatrix")) {
            TTci <- methods::as(L_post, "TsparseMatrix")
            msk <- (TTci@x > 0) & (TTci@i < TTci@j)
            if (any(msk)) {
                ei <- TTci@i[msk] + 1L; ej <- TTci@j[msk] + 1L
            } else { ei <- ej <- integer(0) }
        } else {
            idx_ci <- which(upper.tri(L_post) & L_post > 0, arr.ind = TRUE)
            ei <- idx_ci[, 1]; ej <- idx_ci[, 2]
        }
        if (length(ei)) {
            # Try C++ fast path using full correlation matrix; fallback to R helper if needed
            drop_applied <- FALSE
            rMat_try <- try(.getPearsonMatrix(scope_obj, grid_name = gname, level = "cell"), silent = TRUE)
            if (!inherits(rMat_try, "try-error")) {
                rMat <- rMat_try
                if (!is.matrix(rMat)) {
                    rMat2 <- try(as.matrix(rMat), silent = TRUE)
                    if (!inherits(rMat2, "try-error") && is.matrix(rMat2)) rMat <- rMat2 else rMat <- NULL
                }
                if (!is.null(rMat)) {
                    ridx <- match(kept_genes, rownames(rMat))
                    L_vals <- as.numeric(L_post[cbind(ei, ej)])
                    rule_int <- if (CI_rule == "remove_within") 0L else 1L
                    if (isTRUE(verbose) || isTRUE(debug)) {
                        message(sprintf("[Stage2.CI95] using C++ ci95_drop_mask_edges_omp: edges=%d, threads=%d", length(ei), ncores))
                    }
                    drop_mask <- ci95_drop_mask_edges_omp(rMat, as.integer(ridx), as.integer(ei), as.integer(ej),
                                                          L_vals, xp, lo, hi, rule = rule_int, n_threads = ncores)
                    drop_mask <- as.logical(drop_mask)
                    if (any(drop_mask, na.rm = TRUE)) {
                        ii <- ei[drop_mask]; jj <- ej[drop_mask]
                        L_post[cbind(ii, jj)] <- 0
                        L_post[cbind(jj, ii)] <- 0
                    }
                    drop_applied <- TRUE
                }
            }
            if (!drop_applied) {
                if (isTRUE(verbose) || isTRUE(debug)) {
                    message(sprintf("[Stage2.CI95] using R fallback for edges=%d", length(ei)))
                }
                r_vals <- .get_pairwise_cor_for_edges(scope_obj, grid_name = gname, level = "cell",
                                                     kept_genes = kept_genes, edge_i = ei, edge_j = ej)
                L_vals <- L_post[cbind(ei, ej)]
                if (CI_rule == "remove_within") {
                    drop_mask <- (L_vals <= f_hi(r_vals))
                } else {
                    drop_mask <- (L_vals < f_lo(r_vals) | L_vals > f_hi(r_vals))
                }
                if (any(drop_mask)) {
                    ii <- ei[drop_mask]; jj <- ej[drop_mask]
                    L_post[cbind(ii, jj)] <- 0
                    L_post[cbind(jj, ii)] <- 0
                }
            }
        }
    }

    # CMH correction on L (if requested) — scale per-edge L by CMH similarity
    if (use_cmh_weight && !is.null(leeStat[[cmh_slot]])) {
        g_sim <- leeStat[[cmh_slot]]
        ed_sim <- igraph::as_data_frame(g_sim, what = "edges")
        cmh_col <- if (!is.null(ed_sim$CMH)) "CMH" else if (!is.null(ed_sim$weight)) "weight" else NA_character_
        if (is.na(cmh_col)) stop("CMH graph lacks 'CMH' (or 'weight') attribute")

        # extract upper-tri positive edges once
        if (inherits(L_post, "sparseMatrix")) {
            TTu <- methods::as(L_post, "TsparseMatrix")
            masku <- (TTu@x > 0) & (TTu@i < TTu@j)
            if (any(masku)) { ei <- TTu@i[masku] + 1L; ej <- TTu@j[masku] + 1L } else { ei <- ej <- integer(0) }
        } else {
            idxu <- which(upper.tri(L_post) & L_post > 0, arr.ind = TRUE)
            ei <- idxu[, 1]; ej <- idxu[, 2]
        }
        if (length(ei)) {
            w_vec <- .cmh_lookup_pairs(kept_genes = kept_genes, ei = ei, ej = ej, g_sim = g_sim, weight_col = cmh_col, n_threads = ncores)
            # scale L symmetrically
            L_post[cbind(ei, ej)] <- L_post[cbind(ei, ej)] * w_vec
            L_post[cbind(ej, ei)] <- L_post[cbind(ej, ei)] * w_vec
        }
    }

    # ensure symmetry and zero diagonal after corrections
    L_post <- .symmetric_pmax(L_post)
    if (!inherits(L_post, "sparseMatrix")) L_post <- Matrix::drop0(Matrix::Matrix(L_post, sparse = TRUE)) else L_post <- Matrix::drop0(L_post)

    ## Apply FDR filter on corrected L if requested (robust to missing dimnames)
    if (use_FDR && is.matrix(FDRmat)) {
        # ensure FDRmat aligned to L_raw order and subset by kept_genes
        FM <- FDRmat
        if (is.null(dimnames(FM))) {
            dimnames(FM) <- dimnames(L_raw)
        } else {
            if (!(identical(rownames(FM), rownames(L_raw)) && identical(colnames(FM), colnames(L_raw)))) {
                if (!is.null(rownames(FM)) && !is.null(colnames(FM)) &&
                    all(rownames(L_raw) %in% rownames(FM)) && all(colnames(L_raw) %in% colnames(FM))) {
                    FM <- FM[rownames(L_raw), colnames(L_raw), drop = FALSE]
                } else {
                    dimnames(FM) <- dimnames(L_raw)
                }
            }
        }
        # Use numeric indices as a final fallback to avoid name issues
        ridx <- match(kept_genes, rownames(L_raw))
        cidx <- ridx
        FDR_sub <- FM[ridx, cidx, drop = FALSE]
        L_post[FDR_sub > FDR_max] <- 0
    }

    ## 7. Preserve an intra-cluster backbone to avoid breaking small clusters
    if (isTRUE(keep_stage1_backbone)) {
        # determine a small positive floor from current corrected L
        pos_vals <- as.numeric(L_post[L_post > 0])
        if (length(pos_vals) == 0) {
            w_floor <- 1e-6
        } else {
            qv <- stats::quantile(pos_vals, probs = min(max(backbone_floor_q, 0), 0.25), na.rm = TRUE)
            w_floor <- max(1e-8, as.numeric(qv))
        }

        # fallback similarity from raw L if corrected L has no edges
        L_raw_pos <- pmax(L_raw, 0)

        cl_ids_stage1 <- sort(na.omit(unique(memb_final1)))
        for (cid in cl_ids_stage1) {
            genes_c <- kept_genes[memb_final1[kept_genes] == cid]
            if (length(genes_c) < 2) next
            # prefer corrected weights within the cluster; if empty, fall back to raw L
            M <- L_post[genes_c, genes_c, drop = FALSE]
            if (!any(M > 0, na.rm = TRUE)) {
                M <- L_raw_pos[genes_c, genes_c, drop = FALSE]
            }
            # build graph and compute MST on available positive edges
            if (inherits(M, "sparseMatrix")) {
                TTm <- methods::as(M, "TsparseMatrix")
                maskm <- (TTm@x > 0) & (TTm@i < TTm@j)
                if (any(maskm)) {
                    idx <- cbind(TTm@i[maskm] + 1L, TTm@j[maskm] + 1L)
                } else {
                    idx <- matrix(numeric(0), ncol = 2)
                }
            } else {
                idx <- which(upper.tri(M) & M > 0, arr.ind = TRUE)
            }
            if (nrow(idx) == 0) next
            edf <- data.frame(from = genes_c[idx[,1]], to = genes_c[idx[,2]], w = M[idx], stringsAsFactors = FALSE)
            gtmp <- igraph::graph_from_data_frame(edf, directed = FALSE, vertices = genes_c)
            if (igraph::ecount(gtmp) == 0) next
            mst_c <- igraph::mst(gtmp, weights = 1 / (igraph::E(gtmp)$w + 1e-9))
            if (igraph::ecount(mst_c) == 0) next
            ep <- igraph::as_data_frame(mst_c, what = "edges")
            # raise L_post on these edges to at least w_floor
            for (k in seq_len(nrow(ep))) {
                i <- ep$from[k]; j <- ep$to[k]
                L_post[i, j] <- max(L_post[i, j], w_floor)
                L_post[j, i] <- L_post[i, j]
            }
        }
    }

    ## 8. Build corrected weighted graph with smoothing
    if (inherits(L_post, "sparseMatrix")) {
        TT1 <- methods::as(L_post, "TsparseMatrix")
        xvals1 <- tryCatch(TT1@x, error = function(e) rep(1, length(TT1@i)))
        ed1 <- data.frame(i = TT1@i + 1L, j = TT1@j + 1L, x = xvals1, stringsAsFactors = FALSE)
    } else {
        idx1 <- which(L_post > 0, arr.ind = TRUE)
        ed1 <- data.frame(i = idx1[,1], j = idx1[,2], x = L_post[idx1], stringsAsFactors = FALSE)
    }
    ed1 <- ed1[ed1$i < ed1$j & ed1$x > 0, ]
    if (!nrow(ed1)) stop("No edges remain after post-correction")

    edges_corr <- data.frame(
        from = kept_genes[ed1$i],
        to = kept_genes[ed1$j],
        L_corr = ed1$x,
        stringsAsFactors = FALSE
    )
    edges_corr$w_raw <- if (use_log1p_weight) log1p(edges_corr$L_corr) else edges_corr$L_corr

    if (post_smooth) {
        q <- stats::quantile(edges_corr$w_raw, probs = post_smooth_quant, na.rm = TRUE)
        w_clipped <- pmin(pmax(edges_corr$w_raw, q[1]), q[2])
        rng <- (q[2] - q[1]); if (rng <= 0) rng <- max(1e-8, max(w_clipped) - min(w_clipped))
        w01 <- (w_clipped - min(w_clipped)) / (rng + 1e-12)
        if (!is.null(post_smooth_power) && is.finite(post_smooth_power) && post_smooth_power != 1) {
            w01 <- w01 ^ post_smooth_power
        }
        edges_corr$weight <- w01
    } else {
        edges_corr$weight <- edges_corr$w_raw
    }

    g1 <- igraph::graph_from_data_frame(edges_corr[, c("from", "to", "weight")],
        directed = FALSE, vertices = kept_genes
    )

    ## 8. Second-round multi-run community detection / consensus
    ## Constraint: stage-2 can only split stage-1 clusters (no merges)
    ng <- length(kept_genes)
    memb_final <- rep(NA_integer_, ng)
    names(memb_final) <- kept_genes
    # accumulate sparse consensus triplets; materialize after loop
    cons_I <- integer(0); cons_J <- integer(0); cons_X <- numeric(0)

    cl_ids_stage1 <- sort(na.omit(unique(memb_final1)))
    next_global_id <- 1L
    for (cid in cl_ids_stage1) {
        idx <- which(memb_final1 == cid)
        genes_c <- kept_genes[idx]

        # subgraph of corrected-weight graph restricted to stage-1 cluster
        g1_sub <- igraph::induced_subgraph(g1, vids = genes_c)
        if (length(idx) <= 1 || igraph::ecount(g1_sub) == 0) {
            memb_final[genes_c] <- next_global_id
            next_global_id <- next_global_id + 1L
            next
        }

        ew <- igraph::E(g1_sub)$weight
        el <- igraph::as_data_frame(g1_sub, what = "edges")[, c("from", "to")]

        rvec <- seq_len(n_restart_eff)
        workers <- tryCatch(future::nbrOfWorkers(), error = function(e) NA_integer_)
        if (is.na(workers) || workers < 1) workers <- ncores
        if (is.null(restart_batch_size)) {
            bs <- if (n_restart_eff >= 2 * workers) max(3L, min(20L, as.integer(ceiling(n_restart_eff / workers)))) else 1L
            if (isTRUE(aggressive_cpu)) bs <- max(3L, min(10L, bs))
        } else bs <- as.integer(max(1L, restart_batch_size))
        batches_stage2 <- split(rvec, ceiling(rvec / bs))
        if (isTRUE(verbose) || isTRUE(debug)) {
            message(sprintf("[Stage2] workers=%d, n_restart=%d, batch_size=%d, n_batches=%d, scheduling=%.2f",
                            workers, n_restart_eff, bs, length(batches_stage2), as.numeric(future_scheduling_stage1)))
        }
        memb_list_sub <- future.apply::future_lapply(
            X = batches_stage2,
            FUN = .geneSCOPE_cluster_batch,
            graph = g1_sub,
            algo = algo,
            resolution = resolution,
            objective = objective,
            future.seed = TRUE,
            future.globals = FALSE,
            future.packages = c("igraph", "Matrix", "data.table"),
            future.scheduling = future_scheduling_stage1
        )
        memb_list_sub <- unlist(memb_list_sub, recursive = FALSE, use.names = FALSE)
        memb_mat_sub <- try(do.call(cbind, memb_list_sub), silent = TRUE)
        if (inherits(memb_mat_sub, "try-error") || is.null(dim(memb_mat_sub))) {
            memb_mat_sub <- matrix(unlist(memb_list_sub, use.names = FALSE), nrow = length(genes_c))
        }
        rownames(memb_mat_sub) <- igraph::V(g1_sub)$name

        # per-cluster consensus or single-run assignment (edge-based)
        if (!use_consensus) {
            memb_sub <- memb_mat_sub[, 1]
            cons_w <- as.integer(memb_sub[match(el$from, genes_c)] == memb_sub[match(el$to, genes_c)])
        } else {
            ii <- match(el$from, genes_c); jj <- match(el$to, genes_c)
            if (isTRUE(verbose) || isTRUE(debug)) {
                message(sprintf("[Stage2] consensus on subgraph: |E|=%d, R=%d, threads=%d, cpp_available=%s",
                                nrow(el), ncol(memb_mat_sub), ncores, as.character(exists("consensus_on_edges_omp"))))
            }
            cons_w <- .consensus_on_edges(memb_mat_sub, ii, jj, n_threads = ncores)
        }
        keep_e <- (cons_w >= consensus_thr) & is.finite(cons_w)
        if (any(keep_e)) {
            g_cons_sub <- igraph::graph_from_data_frame(
                data.frame(from = el$from[keep_e], to = el$to[keep_e], weight = cons_w[keep_e]),
                directed = FALSE, vertices = genes_c
            )
        } else {
            # fall back to unweighted subgraph if consensus empty
            g_cons_sub <- igraph::graph_from_data_frame(el[FALSE, ], directed = FALSE, vertices = genes_c)
        }
        if (algo == "leiden") {
            comm_sub <- try(igraph::cluster_leiden(g_cons_sub, weights = igraph::E(g_cons_sub)$weight, resolution = resolution, objective_function = objective), silent = TRUE)
            if (inherits(comm_sub, "try-error")) comm_sub <- igraph::cluster_leiden(g_cons_sub, weights = igraph::E(g_cons_sub)$weight, resolution_parameter = resolution, objective_function = objective)
        } else {
            comm_sub <- try(igraph::cluster_louvain(g_cons_sub, weights = igraph::E(g_cons_sub)$weight, resolution = resolution, objective_function = objective), silent = TRUE)
            if (inherits(comm_sub, "try-error")) {
                comm_sub <- try(igraph::cluster_louvain(g_cons_sub, weights = igraph::E(g_cons_sub)$weight, resolution_parameter = resolution, objective_function = objective), silent = TRUE)
                if (inherits(comm_sub, "try-error")) comm_sub <- igraph::cluster_louvain(g_cons_sub, weights = igraph::E(g_cons_sub)$weight)
            }
        }
        memb_sub <- igraph::membership(comm_sub)

        # Only splits inside this parent cluster: assign new global ids per subcluster
        sub_ids <- sort(unique(memb_sub))
        if (length(sub_ids) <= 1) {
            memb_final[genes_c] <- next_global_id
            next_global_id <- next_global_id + 1L
        } else {
            for (sid in sub_ids) {
                gs <- names(memb_sub)[memb_sub == sid]
                memb_final[gs] <- next_global_id
                next_global_id <- next_global_id + 1L
            }
        }

        # accumulate sparse consensus edges to global coordinates
        if (any(keep_e)) {
            gi <- match(el$from[keep_e], kept_genes)
            gj <- match(el$to[keep_e], kept_genes)
            cons_I <- c(cons_I, gi); cons_J <- c(cons_J, gj); cons_X <- c(cons_X, cons_w[keep_e])
        }
    }

    # materialize sparse global consensus matrix (symmetric)
    cons <- if (length(cons_I)) {
        Matrix::sparseMatrix(i = c(cons_I, cons_J), j = c(cons_J, cons_I), x = c(cons_X, cons_X),
                             dims = c(ng, ng), dimnames = list(kept_genes, kept_genes))
    } else {
        Matrix::sparseMatrix(dims = c(ng, ng), dimnames = list(kept_genes, kept_genes))
    }

    ## 8b. Stability metrics, optional conservative subclustering and QC
    # Build weighted adjacency for post-corrected graph
    W <- igraph::as_adjacency_matrix(g1, attr = "weight", sparse = TRUE)
    rownames(W) <- colnames(W) <- kept_genes

    # helper: compute cluster-level metrics
    metrics_before <- .cluster_metrics_sparse(memb_final, cons, W, kept_genes)
    genes_before <- .gene_metrics_sparse(memb_final, cons, W, kept_genes)

    # Conservative QC: mark very unstable/weakly connected genes as NA
    if (enable_qc_filter) {
        q_w <- stats::quantile(genes_before$w_in, probs = qc_gene_intra_weight_q, na.rm = TRUE)
        bad <- which(!is.na(genes_before$p_in) &
                     genes_before$p_in < qc_gene_intra_cons_min &
                     genes_before$p_best_out > qc_gene_best_out_cons_min &
                     genes_before$w_in <= q_w)
        if (length(bad)) memb_final[match(genes_before$gene[bad], kept_genes)] <- NA
    }

    # Conservative subclustering: attempt only for clearly weak clusters
    if (enable_subcluster) {
        next_id <- if (length(na.omit(memb_final))) max(na.omit(memb_final)) + 1L else 1L
        for (row in seq_len(nrow(metrics_before))) {
            cid <- metrics_before$cluster[row]
            size <- metrics_before$size[row]
            wcons <- metrics_before$within_cons[row]
            condc <- metrics_before$conductance[row]
            if (is.na(size) || size < sub_min_size || is.na(wcons) || is.na(condc)) next
            if (!(wcons < sub_within_cons_max && condc > sub_conductance_min)) next

            idx <- which(memb_final == cid)
            # Build subgraph and re-run consensus at higher resolution
            g_sub <- igraph::induced_subgraph(g1, vids = kept_genes[idx])
            if (igraph::ecount(g_sub) == 0) next
            ew <- igraph::E(g_sub)$weight
            el <- igraph::as_data_frame(g_sub, what = "edges")[, c("from", "to")]
            memb_list2 <- future.apply::future_lapply(
                X = {
                    rvec <- seq_len(n_restart_eff)
                    if (is.null(restart_batch_size)) {
                        workers <- tryCatch(future::nbrOfWorkers(), error = function(e) NA_integer_)
                        if (is.na(workers) || workers < 1) workers <- ncores
                        bs <- if (n_restart_eff >= 2 * workers) max(3L, min(20L, as.integer(ceiling(n_restart_eff / workers)))) else 1L
                    } else bs <- as.integer(max(1L, restart_batch_size))
                    split(rvec, ceiling(rvec / bs))
                },
                FUN = .geneSCOPE_cluster_batch,
                graph = g_sub,
                algo = algo,
                resolution = resolution * sub_resolution_factor,
                objective = objective,
                future.seed = TRUE,
                future.globals = FALSE,
                future.packages = c("igraph", "Matrix", "data.table"),
                future.scheduling = future_scheduling_stage1
            )
            memb_mat_sub <- try(do.call(cbind, memb_list2), silent = TRUE)
            if (inherits(memb_mat_sub, "try-error") || is.null(dim(memb_mat_sub))) memb_mat_sub <- matrix(unlist(memb_list2, use.names = FALSE), ncol = 1)
            rownames(memb_mat_sub) <- igraph::V(g_sub)$name
            # Edge-based consensus inside parent cluster
            ns <- length(idx)
            el2 <- igraph::as_data_frame(g_sub, what = "edges")[, c("from", "to")]
            ii2 <- match(el2$from, igraph::V(g_sub)$name)
            jj2 <- match(el2$to,   igraph::V(g_sub)$name)
            cons_w2 <- .consensus_on_edges(memb_mat_sub, ii2, jj2)
            g_cons_sub <- igraph::graph_from_data_frame(
                data.frame(from = el2$from, to = el2$to, weight = cons_w2),
                directed = FALSE, vertices = igraph::V(g_sub)$name
            )
            if (algo == "leiden") {
                comm_sub <- try(igraph::cluster_leiden(g_cons_sub, weights = igraph::E(g_cons_sub)$weight, resolution = resolution * sub_resolution_factor, objective_function = objective), silent = TRUE)
                if (inherits(comm_sub, "try-error")) comm_sub <- igraph::cluster_leiden(g_cons_sub, weights = igraph::E(g_cons_sub)$weight, resolution_parameter = resolution * sub_resolution_factor, objective_function = objective)
            } else {
                comm_sub <- try(igraph::cluster_louvain(g_cons_sub, weights = igraph::E(g_cons_sub)$weight, resolution = resolution * sub_resolution_factor, objective_function = objective), silent = TRUE)
                if (inherits(comm_sub, "try-error")) {
                    comm_sub <- try(igraph::cluster_louvain(g_cons_sub, weights = igraph::E(g_cons_sub)$weight, resolution_parameter = resolution * sub_resolution_factor, objective_function = objective), silent = TRUE)
                    if (inherits(comm_sub, "try-error")) comm_sub <- igraph::cluster_louvain(g_cons_sub, weights = igraph::E(g_cons_sub)$weight)
                }
            }
            memb_sub <- igraph::membership(comm_sub)

            # evaluate split
            sub_ids <- sort(unique(memb_sub))
            if (length(sub_ids) <= 1 || length(sub_ids) > sub_max_groups) next
            sizes <- vapply(sub_ids, function(sid) sum(memb_sub == sid), integer(1))
            if (any(sizes < sub_min_child_size)) next
            # within-consensus improvement (weighted) – approximate by averaging
            wcons_sub <- vapply(sub_ids, function(sid) {
                genes_sid <- names(memb_sub)[memb_sub == sid]
                if (length(genes_sid) <= 1) return(0)
                msk <- el2$from %in% genes_sid & el2$to %in% genes_sid
                if (any(msk)) mean(cons_w2[msk], na.rm = TRUE) else 0
            }, numeric(1))
            wcons_old <- wcons
            wcons_new <- sum(wcons_sub * sizes) / sum(sizes)
            if (!(is.finite(wcons_new) && (wcons_new >= (wcons_old + sub_improve_within_cons_min)))) next

            # accept split: assign new ids
            # map each subcluster to a new global id
            for (sid in sub_ids) {
                genes_sid <- names(memb_sub)[memb_sub == sid]
                memb_final[match(genes_sid, kept_genes)] <- next_id
                next_id <- next_id + 1L
            }
        }
    }

    ## 9. Clean isolated / micro clusters (based on corrected L) — O(E)
    if (inherits(L_post, "sparseMatrix")) {
        TTi <- methods::as(L_post[kept_genes, kept_genes, drop = FALSE], "TsparseMatrix")
        msk <- (TTi@x > 0) & (TTi@i < TTi@j)
        if (any(msk)) { ei <- TTi@i[msk] + 1L; ej <- TTi@j[msk] + 1L } else { ei <- ej <- integer(0) }
    } else {
        idx_in <- which(upper.tri(L_post[kept_genes, kept_genes, drop = FALSE]) & L_post[kept_genes, kept_genes, drop = FALSE] > 0, arr.ind = TRUE)
        ei <- idx_in[, 1]; ej <- idx_in[, 2]
    }
    di <- .intra_degree_from_edges(ei, ej, memb_final)
    memb_final[di == 0] <- NA
    tbl <- table(memb_final, useNA = "no")
    small <- names(tbl)[tbl < 2]
    memb_final[memb_final %in% small] <- NA
    stable_cl <- names(tbl)[tbl >= 2]
    names(memb_final) <- kept_genes

    # prepare optional report
    if (return_report) {
        report <- new.env(parent = emptyenv())
        report$thresholds <- list(
            sub_min_size = sub_min_size,
            sub_min_child_size = sub_min_child_size,
            sub_resolution_factor = sub_resolution_factor,
            sub_within_cons_max = sub_within_cons_max,
            sub_conductance_min = sub_conductance_min,
            sub_improve_within_cons_min = sub_improve_within_cons_min,
            qc_gene_intra_cons_min = qc_gene_intra_cons_min,
            qc_gene_best_out_cons_min = qc_gene_best_out_cons_min,
            qc_gene_intra_weight_q = qc_gene_intra_weight_q
        )
        report$metrics_before <- metrics_before
        report$genes_before <- genes_before
        # recompute after- metrics using updated memb_final
        report$metrics_after <- .cluster_metrics(memb_final, cons, W)
        report$genes_after <- .gene_metrics(memb_final, cons, W)
    }

    ## 10. Build final consensus graph using corrected (and smoothed) weights
    if (inherits(cons, "sparseMatrix")) {
        TTc <- methods::as(cons, "TsparseMatrix")
        maskc <- (TTc@x > 0) & (TTc@i < TTc@j)
        if (any(maskc)) {
            edges_keep <- cbind(TTc@i[maskc] + 1L, TTc@j[maskc] + 1L)
        } else {
            edges_keep <- matrix(numeric(0), ncol = 2)
        }
    } else {
        edges_keep <- which(cons > 0, arr.ind = TRUE)
        edges_keep <- edges_keep[edges_keep[, 1] < edges_keep[, 2], ]
    }

    # mapping from edge key to smoothed weight
    key_corr <- paste(pmin(edges_corr$from, edges_corr$to), pmax(edges_corr$from, edges_corr$to), sep = "|")
    w_map <- setNames(edges_corr$weight, key_corr)

    intra_keys <- paste(
        pmin(kept_genes[edges_keep[, 1]], kept_genes[edges_keep[, 2]]),
        pmax(kept_genes[edges_keep[, 1]], kept_genes[edges_keep[, 2]]),
        sep = "|"
    )
    w_vec <- as.numeric(w_map[intra_keys])
    keep_mask <- !is.na(w_vec) & w_vec > 0
    intra_df <- data.frame(
        from = kept_genes[edges_keep[keep_mask, 1]],
        to   = kept_genes[edges_keep[keep_mask, 2]],
        weight = w_vec[keep_mask],
        sign = ifelse(L_raw[cbind(
            match(kept_genes[edges_keep[keep_mask, 1]], rownames(L_raw)),
            match(kept_genes[edges_keep[keep_mask, 2]], colnames(L_raw))
        )] < 0, "neg", "pos"),
        stringsAsFactors = FALSE
    )

    cross_df <- NULL
    if (keep_cross_stable && length(stable_cl) > 1) {
        idx_stable <- memb_final %in% stable_cl
        genes_stable <- kept_genes[idx_stable]

        # consider cross-cluster edges among stable genes
        L_sub_corr <- L_post[genes_stable, genes_stable, drop = FALSE]
        if (inherits(L_sub_corr, "sparseMatrix")) {
            TTs <- methods::as(L_sub_corr, "TsparseMatrix")
            masksc <- (TTs@x >= L_min) & (TTs@i < TTs@j)
            if (any(masksc)) {
                ij <- cbind(TTs@i[masksc] + 1L, TTs@j[masksc] + 1L)
            } else {
                ij <- matrix(numeric(0), ncol = 2)
            }
        } else {
            ij <- which(L_sub_corr >= L_min & upper.tri(L_sub_corr), arr.ind = TRUE)
        }
        if (nrow(ij)) {
            g1 <- genes_stable[ij[, 1]]
            g2 <- genes_stable[ij[, 2]]
            cl1 <- memb_final[g1]; cl2 <- memb_final[g2]
            keep_cross <- cl1 != cl2 & !is.na(cl1) & !is.na(cl2)
            if (any(keep_cross)) {
                g1 <- g1[keep_cross]; g2 <- g2[keep_cross]
                kx <- paste(pmin(g1, g2), pmax(g1, g2), sep = "|")
                w2 <- as.numeric(w_map[kx])
                mask <- !is.na(w2) & w2 > 0
                if (any(mask)) {
                    # sign from raw L
                    Lval_raw <- L_raw[cbind(match(g1[mask], rownames(L_raw)), match(g2[mask], colnames(L_raw)))]
                    cross_df <- data.frame(
                        from = g1[mask], to = g2[mask],
                        weight = w2[mask],
                        sign = ifelse(Lval_raw < 0, "neg", "pos")
                    )
                }
            }
        }
    }

    edges_df_final <- if (!is.null(cross_df)) unique(rbind(intra_df, cross_df)) else intra_df
    g_cons <- igraph::graph_from_data_frame(edges_df_final,
        directed = FALSE,
        vertices = kept_genes
    )

    ## 11. Continuous numbering
    if (length(memb_final)) {
        valid_idx <- !is.na(memb_final)
        old_ids <- sort(unique(memb_final[valid_idx]))
        id_map <- setNames(seq_along(old_ids), old_ids)
        memb_final[valid_idx] <- id_map[as.character(memb_final[valid_idx])]
    }

    ## 12. Write to meta.data
    if (is.null(cluster_name) || !nzchar(cluster_name)) {
        cluster_name <- sprintf("modL%.2f_post", L_min)
    }
    if (is.null(scope_obj@meta.data)) {
        scope_obj@meta.data <- data.frame(row.names = genes_all)
    } else {
        miss_rows <- setdiff(genes_all, rownames(scope_obj@meta.data))
        if (length(miss_rows)) {
            scope_obj@meta.data <- rbind(scope_obj@meta.data, data.frame(row.names = miss_rows))
        }
    }
    scope_obj@meta.data[genes_all, cluster_name] <- NA_integer_
    scope_obj@meta.data[kept_genes, cluster_name] <- as.integer(memb_final[kept_genes])

    if (verbose) {
        nz <- sum(!is.na(scope_obj@meta.data[genes_all, cluster_name]))
        message("[geneSCOPE::clusterGenes] Written clustering column '", cluster_name, "': ",
                nz, "/", length(genes_all), " genes (",
                sprintf("%.1f%%", 100 * nz / length(genes_all)), ")")
    }

    ## 13. Save consensus graph
    if (is.null(scope_obj@stats[[gname]][[lee_stats_layer]])) {
        scope_obj@stats[[gname]][[lee_stats_layer]] <- list()
    }
    scope_obj@stats[[gname]][[lee_stats_layer]][[graph_slot_name]] <- g_cons

    ## 14. Return
    if (isTRUE(return_report)) {
        return(list(scope_obj = scope_obj, report = report))
    } else {
        return(scope_obj)
    }
}
