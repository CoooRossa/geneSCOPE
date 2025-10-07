#' @title Gene Clustering (results written to scope_obj@meta.data)
#' @description Based on Lee's L statistics for filtering, graph construction, multiple community detection and optional consensus, output cluster numbers to scope_obj@meta.data.
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
    verbose = TRUE) {

    algo <- match.arg(algo)
    objective <- match.arg(objective)
    CI_rule <- match.arg(CI_rule)

    ## 1. Select grid layer
    g_layer <- .selectGridLayer(scope_obj, grid_name)
    gname <- names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    .checkGridContent(scope_obj, gname)

    ## 2. Read statistics
    if (is.null(scope_obj@stats[[gname]]) ||
        is.null(scope_obj@stats[[gname]][[lee_stats_layer]])) {
        stop("Statistics layer missing: ", lee_stats_layer)
    }
    leeStat <- scope_obj@stats[[gname]][[lee_stats_layer]]
    L_raw <- leeStat$L
    if (is.null(L_raw)) stop("LeeStats layer missing L matrix")
    genes_all <- rownames(L_raw)

    FDRmat <- if (use_FDR) {
        if (is.null(leeStat$FDR)) stop("use_FDR=TRUE but FDR matrix not found")
        leeStat$FDR
    } else {
        matrix(0, nrow(L_raw), ncol(L_raw), dimnames = dimnames(L_raw))
    }

    ## 3. Threshold filtering (base, before any correction)
    A <- L_raw
    A[abs(A) < L_min] <- 0
    if (use_FDR) A[FDRmat > FDR_max] <- 0
    A[A < 0] <- 0
    A <- pmax(A, t(A))
    diag(A) <- 0
    A <- .filter_matrix_by_quantile(A, pct_min, "q100")

    ## 4. Base graph, without CI95/CMH corrections
    keep <- if (drop_isolated) Matrix::rowSums(A != 0) > 0 else rep(TRUE, nrow(A))
    kept_genes <- genes_all[keep]
    if (length(kept_genes) < 2) stop("Fewer than 2 genes after filtering")

    A_sub <- A[kept_genes, kept_genes, drop = FALSE]
    ed0 <- summary(A_sub)
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
    algo_fun <- if (algo == "leiden") igraph::cluster_leiden else igraph::cluster_louvain

    old_max_size <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 2 * 1024^3)
    on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)

    edge_weights <- igraph::E(g0)$weight
    vertex_names <- igraph::V(g0)$name
    edge_list <- igraph::as_data_frame(g0, what = "edges")[, c("from", "to")]

    memb_mat1 <- future.apply::future_sapply(
        seq_len(n_restart),
        function(ii) {
            g_local <- igraph::graph_from_data_frame(
                cbind(edge_list, weight = edge_weights),
                directed = FALSE,
                vertices = vertex_names
            )
            if (algo == "leiden") {
                igraph::cluster_leiden(
                    g_local,
                    weights = igraph::E(g_local)$weight,
                    resolution_parameter = resolution,
                    objective_function = objective
                )$membership
            } else {
                igraph::cluster_louvain(
                    g_local,
                    weights = igraph::E(g_local)$weight,
                    resolution_parameter = resolution,
                    objective_function = objective
                )$membership
            }
        },
        future.seed = TRUE
    )
    rownames(memb_mat1) <- kept_genes

    if (!use_consensus) {
        memb_final1 <- memb_mat1[, 1]
        cons1 <- outer(memb_final1, memb_final1, "==") * 1
        diag(cons1) <- 0
    } else {
        ng <- length(kept_genes)
        cons1 <- matrix(0, ng, ng)
        for (k in seq_len(ncol(memb_mat1))) cons1 <- cons1 + (outer(memb_mat1[, k], memb_mat1[, k], "=="))
        cons1 <- cons1 / ncol(memb_mat1)
        diag(cons1) <- 0
        cons1[cons1 < consensus_thr] <- 0
        cons1 <- Matrix::drop0(cons1)
        g_cons_tmp1 <- igraph::graph_from_adjacency_matrix(cons1, mode = "undirected", weighted = TRUE, diag = FALSE)
        comm1 <- algo_fun(
            g_cons_tmp1,
            weights = igraph::E(g_cons_tmp1)$weight,
            resolution_parameter = resolution,
            objective_function = objective
        )
        memb_final1 <- igraph::membership(comm1)
    }

    ## Early exit: if no post-correction requested, return stage-1 result
    if (!use_cmh_weight && !CI95_filter) {
        memb_final <- memb_final1
        cons <- cons1

        # Clean isolated / micro clusters using base adjacency A_sub
        di <- Matrix::rowSums(A_sub[kept_genes, kept_genes] > 0 &
            outer(memb_final, memb_final, "=="))
        memb_final[di == 0] <- NA
        tbl <- table(memb_final, useNA = "no")
        small <- names(tbl)[tbl < 2]
        memb_final[memb_final %in% small] <- NA
        stable_cl <- names(tbl)[tbl >= 2]
        names(memb_final) <- kept_genes

        # Build consensus graph (same as clusterGenes base path)
        edges_keep <- which(cons > 0, arr.ind = TRUE)
        edges_keep <- edges_keep[edges_keep[, 1] < edges_keep[, 2], ]
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
            ij <- which(L_sub_raw >= L_min & upper.tri(L_sub_raw), arr.ind = TRUE)
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
            message("[geneSCOPE::clusterGenes] (step1 only) Written clustering column '", cluster_name, "': ",
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

    ## 6. Correct L by CI95 and/or CMH (post step)
    # start from raw L among kept genes
    L_post <- L_raw[kept_genes, kept_genes, drop = FALSE]
    # apply baseline thresholding again to keep consistency
    L_post[abs(L_post) < L_min] <- 0
    if (use_FDR) {
        # robust subsetting regardless of FDR dimnames presence
        FM <- FDRmat
        if (!is.matrix(FM)) FM <- as.matrix(FM)
        if (is.null(dimnames(FM))) {
            dimnames(FM) <- dimnames(L_raw)
        } else if (!(identical(rownames(FM), rownames(L_raw)) && identical(colnames(FM), colnames(L_raw)))) {
            if (!is.null(rownames(FM)) && !is.null(colnames(FM)) &&
                all(rownames(L_raw) %in% rownames(FM)) && all(colnames(L_raw) %in% colnames(FM))) {
                FM <- FM[rownames(L_raw), colnames(L_raw), drop = FALSE]
            } else {
                dimnames(FM) <- dimnames(L_raw)
            }
        }
        ridx <- match(kept_genes, rownames(L_raw))
        FDR_sub <- FM[ridx, ridx, drop = FALSE]
        L_post[FDR_sub > FDR_max] <- 0
    }
    L_post[L_post < 0] <- 0
    L_post <- pmax(L_post, t(L_post)); diag(L_post) <- 0
    L_post <- .filter_matrix_by_quantile(L_post, pct_min, "q100")

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
        rMat <- .getPearsonMatrix(scope_obj, grid_name = gname, level = "cell")
        rMat <- as.matrix(rMat[kept_genes, kept_genes, drop = FALSE])
        mask <- if (CI_rule == "remove_within") (L_post <= f_hi(rMat)) else (L_post < f_lo(rMat) | L_post > f_hi(rMat))
        L_post[mask] <- 0
    }

    # CMH correction on L (if requested) — scale per-edge L by CMH similarity
    if (use_cmh_weight && !is.null(leeStat[[cmh_slot]])) {
        g_sim <- leeStat[[cmh_slot]]
        ed_sim <- igraph::as_data_frame(g_sim, what = "edges")
        key_sim <- paste(pmin(ed_sim$from, ed_sim$to), pmax(ed_sim$from, ed_sim$to), sep = "|")
        sim_map <- setNames(ed_sim$CMH, key_sim)

        # build key map for kept genes
        idx <- which(upper.tri(L_post) & L_post > 0, arr.ind = TRUE)
        if (nrow(idx)) {
            g1 <- kept_genes[idx[, 1]]
            g2 <- kept_genes[idx[, 2]]
            k <- paste(pmin(g1, g2), pmax(g1, g2), sep = "|")
            s <- sim_map[k]
            s[is.na(s)] <- median(sim_map, na.rm = TRUE)
            # scale L by CMH
            L_post[cbind(idx[, 1], idx[, 2])] <- L_post[cbind(idx[, 1], idx[, 2])] * s
            L_post[cbind(idx[, 2], idx[, 1])] <- L_post[cbind(idx[, 2], idx[, 1])] * s
        }
    }

    # ensure symmetry and zero diagonal after corrections
    L_post <- pmax(L_post, t(L_post)); diag(L_post) <- 0

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
            idx <- which(upper.tri(M) & M > 0, arr.ind = TRUE)
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
    ed1 <- summary(L_post)
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
    cons <- matrix(0, ng, ng, dimnames = list(kept_genes, kept_genes))

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

        memb_mat_sub <- future.apply::future_sapply(
            seq_len(n_restart),
            function(ii) {
                g_local <- igraph::graph_from_data_frame(cbind(el, weight = ew), directed = FALSE, vertices = igraph::V(g1_sub)$name)
                if (algo == "leiden") {
                    igraph::cluster_leiden(
                        g_local,
                        weights = igraph::E(g_local)$weight,
                        resolution_parameter = resolution,
                        objective_function = objective
                    )$membership
                } else {
                    igraph::cluster_louvain(
                        g_local,
                        weights = igraph::E(g_local)$weight,
                        resolution_parameter = resolution,
                        objective_function = objective
                    )$membership
                }
            },
            future.seed = TRUE
        )
        if (is.null(dim(memb_mat_sub))) memb_mat_sub <- matrix(memb_mat_sub, ncol = 1)
        rownames(memb_mat_sub) <- igraph::V(g1_sub)$name

        # per-cluster consensus or single-run assignment
        if (!use_consensus) {
            memb_sub <- memb_mat_sub[, 1]
            cons_sub <- outer(memb_sub, memb_sub, "==") * 1
            diag(cons_sub) <- 0
        } else {
            ns <- length(genes_c)
            cons_sub <- matrix(0, ns, ns, dimnames = list(genes_c, genes_c))
            for (k in seq_len(ncol(memb_mat_sub))) cons_sub <- cons_sub + (outer(memb_mat_sub[, k], memb_mat_sub[, k], "=="))
            cons_sub <- cons_sub / ncol(memb_mat_sub)
            diag(cons_sub) <- 0
            cons_sub[cons_sub < consensus_thr] <- 0
            cons_sub <- Matrix::drop0(cons_sub)
            g_cons_sub <- igraph::graph_from_adjacency_matrix(cons_sub, mode = "undirected", weighted = TRUE, diag = FALSE)
            comm_sub <- if (algo == "leiden") igraph::cluster_leiden(g_cons_sub, weights = igraph::E(g_cons_sub)$weight, resolution_parameter = resolution, objective_function = objective) else igraph::cluster_louvain(g_cons_sub, weights = igraph::E(g_cons_sub)$weight, resolution_parameter = resolution, objective_function = objective)
            memb_sub <- igraph::membership(comm_sub)
        }

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

        # fill block of consensus matrix
        cons[genes_c, genes_c] <- as.matrix(cons_sub[genes_c, genes_c])
    }

    ## 8b. Stability metrics, optional conservative subclustering and QC
    # Build weighted adjacency for post-corrected graph
    W <- as.matrix(igraph::as_adjacency_matrix(g1, attr = "weight", sparse = TRUE))
    rownames(W) <- colnames(W) <- kept_genes

    metrics_before <- .cluster_metrics(memb_final, cons, W)
    genes_before <- .gene_metrics(memb_final, cons, W)

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
            memb_mat_sub <- future.apply::future_sapply(
                seq_len(min(n_restart, 50)),
                function(ii) {
                    g_local <- igraph::graph_from_data_frame(cbind(el, weight = ew), directed = FALSE, vertices = igraph::V(g_sub)$name)
                    if (algo == "leiden") {
                        igraph::cluster_leiden(
                            g_local,
                            weights = igraph::E(g_local)$weight,
                            resolution_parameter = resolution * sub_resolution_factor,
                            objective_function = objective
                        )$membership
                    } else {
                        igraph::cluster_louvain(
                            g_local,
                            weights = igraph::E(g_local)$weight,
                            resolution_parameter = resolution * sub_resolution_factor,
                            objective_function = objective
                        )$membership
                    }
                },
                future.seed = TRUE
            )
            if (is.null(dim(memb_mat_sub))) memb_mat_sub <- matrix(memb_mat_sub, ncol = 1)
            rownames(memb_mat_sub) <- igraph::V(g_sub)$name
            # consensus inside parent cluster
            ns <- length(idx)
            cons_sub <- matrix(0, ns, ns, dimnames = list(igraph::V(g_sub)$name, igraph::V(g_sub)$name))
            for (k in seq_len(ncol(memb_mat_sub))) cons_sub <- cons_sub + (outer(memb_mat_sub[, k], memb_mat_sub[, k], "=="))
            cons_sub <- cons_sub / ncol(memb_mat_sub)
            diag(cons_sub) <- 0
            g_cons_sub <- igraph::graph_from_adjacency_matrix(cons_sub, mode = "undirected", weighted = TRUE, diag = FALSE)
            comm_sub <- if (algo == "leiden") igraph::cluster_leiden(g_cons_sub, weights = igraph::E(g_cons_sub)$weight, resolution_parameter = resolution * sub_resolution_factor, objective_function = objective) else igraph::cluster_louvain(g_cons_sub, weights = igraph::E(g_cons_sub)$weight, resolution_parameter = resolution * sub_resolution_factor, objective_function = objective)
            memb_sub <- igraph::membership(comm_sub)

            # evaluate split
            sub_ids <- sort(unique(memb_sub))
            if (length(sub_ids) <= 1 || length(sub_ids) > sub_max_groups) next
            sizes <- vapply(sub_ids, function(sid) sum(memb_sub == sid), integer(1))
            if (any(sizes < sub_min_child_size)) next
            # within-consensus improvement (weighted)
            wcons_sub <- vapply(sub_ids, function(sid) {
                jj <- which(names(memb_sub) %in% names(memb_sub)[memb_sub == sid])
                if (length(jj) > 1) mean(cons_sub[jj, jj][upper.tri(cons_sub[jj, jj])], na.rm = TRUE) else 0
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

    ## 9. Clean isolated / micro clusters (based on corrected L)
    di <- Matrix::rowSums((L_post[kept_genes, kept_genes] > 0) & outer(memb_final, memb_final, "=="))
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
    edges_keep <- which(cons > 0, arr.ind = TRUE)
    edges_keep <- edges_keep[edges_keep[, 1] < edges_keep[, 2], ]

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
        ij <- which(L_sub_corr >= L_min & upper.tri(L_sub_corr), arr.ind = TRUE)
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

