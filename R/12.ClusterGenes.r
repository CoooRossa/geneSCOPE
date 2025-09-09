#' @title Gene Clustering (results written to coordObj@meta.data)
#' @description Based on Lee's L statistics for filtering, graph construction, multiple community detection and optional consensus, output cluster numbers to coordObj@meta.data.
#' @export
clusterGenes <- function(
    coordObj,
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
    verbose = TRUE) {
    algo <- match.arg(algo)
    objective <- match.arg(objective)
    CI_rule <- match.arg(CI_rule)

    ## 1. Select grid layer
    g_layer <- .selectGridLayer(coordObj, grid_name)
    gname <- names(coordObj@grid)[vapply(coordObj@grid, identical, logical(1), g_layer)]
    .checkGridContent(coordObj, gname)

    ## 2. Read statistics
    if (is.null(coordObj@stats[[gname]]) ||
        is.null(coordObj@stats[[gname]][[lee_stats_layer]])) {
        stop("Statistics layer missing: ", lee_stats_layer)
    }
    leeStat <- coordObj@stats[[gname]][[lee_stats_layer]]
    L_raw <- leeStat$L
    if (is.null(L_raw)) stop("LeeStats layer missing L matrix")
    genes_all <- rownames(L_raw)

    FDRmat <- if (use_FDR) {
        if (is.null(leeStat$FDR)) stop("use_FDR=TRUE but FDR matrix not found")
        leeStat$FDR
    } else {
        matrix(0, nrow(L_raw), ncol(L_raw), dimnames = dimnames(L_raw))
    }

    ## 3. Threshold filtering
    A <- L_raw
    A[abs(A) < L_min] <- 0
    if (use_FDR) A[FDRmat > FDR_max] <- 0
    A[A < 0] <- 0
    A <- pmax(A, t(A))
    diag(A) <- 0
    A <- .filter_matrix_by_quantile(A, pct_min, "q100")

    ## 4. CI filtering (optional)
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
        rMat <- .getPearsonMatrix(coordObj, grid_name = gname, level = "cell")
        rMat <- as.matrix(rMat[genes_all, genes_all, drop = FALSE])
        mask <- if (CI_rule == "remove_within") (L_raw <= f_hi(rMat)) else (L_raw < f_lo(rMat) | L_raw > f_hi(rMat))
        A[mask] <- 0
    }

    ## 5. Remove isolated
    keep <- if (drop_isolated) Matrix::rowSums(A != 0) > 0 else rep(TRUE, nrow(A))
    kept_genes <- genes_all[keep]
    if (length(kept_genes) < 2) stop("Fewer than 2 genes after filtering")

    ## 6. CMH weight (optional)
    sim_map <- NULL
    if (use_cmh_weight && !is.null(leeStat[[cmh_slot]])) {
        g_sim <- leeStat[[cmh_slot]]
        ed_sim <- igraph::as_data_frame(g_sim, what = "edges")
        key_sim <- paste(pmin(ed_sim$from, ed_sim$to), pmax(ed_sim$from, ed_sim$to), sep = "|")
        sim_map <- setNames(ed_sim$CMH, key_sim)
    }

    ## 7. Edge data frame
    A_sub <- A[kept_genes, kept_genes, drop = FALSE]
    ed0 <- summary(A_sub)
    ed0 <- ed0[ed0$i < ed0$j & ed0$x > 0, ]
    if (!nrow(ed0)) stop("All edges filtered out")
    edges_df <- data.frame(
        from = kept_genes[ed0$i],
        to = kept_genes[ed0$j],
        L_raw = ed0$x,
        stringsAsFactors = FALSE
    )
    use_log1p_local <- if (use_cmh_weight) TRUE else use_log1p_weight
    edges_df$w_L <- if (use_log1p_local) log1p(edges_df$L_raw) else edges_df$L_raw
    if (!is.null(sim_map)) {
        k <- paste(pmin(edges_df$from, edges_df$to), pmax(edges_df$from, edges_df$to), sep = "|")
        s <- sim_map[k]
        s[is.na(s)] <- median(sim_map, na.rm = TRUE)
    } else {
        s <- 1
    }
    edges_df$weight <- edges_df$w_L * s
    edges_df <- edges_df[edges_df$weight > 0, ]
    if (!nrow(edges_df)) stop("No valid edges after weighting")
    g0 <- igraph::graph_from_data_frame(edges_df[, c("from", "to", "weight")],
        directed = FALSE, vertices = kept_genes
    )

    ## 8. Multiple communities / consensus
    algo_fun <- if (algo == "leiden") igraph::cluster_leiden else igraph::cluster_louvain

    # Increase future.globals.maxSize temporarily to handle large objects
    old_max_size <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 2 * 1024^3) # 2GB
    on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)

    # Extract necessary components to minimize data transfer
    edge_weights <- igraph::E(g0)$weight
    vertex_names <- igraph::V(g0)$name
    edge_list <- igraph::as_data_frame(g0, what = "edges")[, c("from", "to")]

    memb_mat <- future.apply::future_sapply(
        seq_len(n_restart),
        function(ii) {
            # Reconstruct minimal graph in worker
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
    rownames(memb_mat) <- kept_genes
    ng <- length(kept_genes)

    if (!use_consensus) {
        memb_final <- memb_mat[, 1]
        cons <- outer(memb_final, memb_final, "==") * 1
        diag(cons) <- 0
    } else {
        cons <- matrix(0, ng, ng)
        for (k in seq_len(ncol(memb_mat))) {
            cons <- cons + (outer(memb_mat[, k], memb_mat[, k], "=="))
        }
        cons <- cons / ncol(memb_mat)
        diag(cons) <- 0
        cons[cons < consensus_thr] <- 0
        cons <- Matrix::drop0(cons)
        g_cons_tmp <- igraph::graph_from_adjacency_matrix(
            cons,
            mode = "undirected", weighted = TRUE, diag = FALSE
        )
        comm_fin <- algo_fun(
            g_cons_tmp,
            weights = igraph::E(g_cons_tmp)$weight,
            resolution_parameter = resolution,
            objective_function = objective
        )
        memb_final <- igraph::membership(comm_fin)
    }

    ## 9. Clean isolated / micro clusters
    di <- Matrix::rowSums(A_sub[kept_genes, kept_genes] > 0 &
        outer(memb_final, memb_final, "=="))
    memb_final[di == 0] <- NA
    tbl <- table(memb_final, useNA = "no")
    small <- names(tbl)[tbl < 2]
    memb_final[memb_final %in% small] <- NA
    stable_cl <- names(tbl)[tbl >= 2]
    names(memb_final) <- kept_genes

    ## 10. Consensus graph (including cross-cluster stable edges optional)
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
        g1 <- genes_stable[ij[, 1]]
        g2 <- genes_stable[ij[, 2]]
        cl1 <- memb_final[g1]
        cl2 <- memb_final[g2]
        keep_cross <- cl1 != cl2 & !is.na(cl1) & !is.na(cl2)
        if (any(keep_cross)) {
            g1 <- g1[keep_cross]
            g2 <- g2[keep_cross]
            Lval <- L_sub_raw[cbind(match(g1, genes_stable), match(g2, genes_stable))]
            wL <- if (use_log1p_local) log1p(Lval) else Lval
            if (!is.null(sim_map)) {
                kx <- paste(pmin(g1, g2), pmax(g1, g2), sep = "|")
                sv <- sim_map[kx]
                sv[is.na(sv)] <- median(sim_map, na.rm = TRUE)
            } else {
                sv <- 1
            }
            wgt <- wL * sv
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
    edges_df_final <- unique(rbind(intra_df, cross_df))
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

    ## 12. Write to meta.data (key fix point)
    if (is.null(cluster_name) || !nzchar(cluster_name)) {
        cluster_name <- sprintf("modL%.2f", L_min)
    }
    # Ensure meta.data exists and contains all gene rows
    if (is.null(coordObj@meta.data)) {
        coordObj@meta.data <- data.frame(row.names = genes_all)
    } else {
        miss_rows <- setdiff(genes_all, rownames(coordObj@meta.data))
        if (length(miss_rows)) {
            coordObj@meta.data <- rbind(
                coordObj@meta.data,
                data.frame(row.names = miss_rows)
            )
        }
    }
    # First set entire column to NA, ensuring column exists
    coordObj@meta.data[genes_all, cluster_name] <- NA_integer_
    coordObj@meta.data[kept_genes, cluster_name] <- as.integer(memb_final[kept_genes])

    if (verbose) {
        nz <- sum(!is.na(coordObj@meta.data[genes_all, cluster_name]))
        message(
            "[geneSCOPE] Written clustering column '", cluster_name, "': ",
            nz, "/", length(genes_all), " genes (",
            sprintf("%.1f%%", 100 * nz / length(genes_all)), ")"
        )
    }

    ## 13. Save consensus graph
    if (is.null(coordObj@stats[[gname]][[lee_stats_layer]])) {
        coordObj@stats[[gname]][[lee_stats_layer]] <- list()
    }
    coordObj@stats[[gname]][[lee_stats_layer]][[graph_slot_name]] <- g_cons

    ## 14. Return (explicit return to avoid upper layer loss)
    coordObj
}
