#' Modules:
#'   * .extract_scope_layers()        – extract primary matrix / aux slots
#'   * .threshold_similarity_graph()  – base thresholding + adjacency build
#'   * .stage1_consensus_workflow()   – multi-run community detection + consensus
#'   * .stage2_refine_workflow()      – CI95 / MH corrections + Stage-2
#'   * clusterGenes()                 – user facing wrapper (scope_obj in/out)
# ---------------------------------------------------------------------------
# Utility: build default configuration list
# ---------------------------------------------------------------------------
#' Build default clustering configuration.
#'
#' @param min_cutoff Numeric Lee's L cutoff.
#' @param use_significance Logical; apply FDR alignment.
#' @param significance_max Maximum allowed FDR value.
#' @param pct_min Quantile filter specification.
#' @param drop_isolated Logical; drop zero-degree nodes.
#' @param algo Stage-1 algorithm name.
#' @param stage2_algo Optional Stage-2 algorithm name.
#' @param resolution CPM resolution parameter.
#' @param gamma Modularity resolution parameter.
#' @param objective Objective name (`"CPM"` or `"modularity"`).
#' @param use_log1p_weight Logical; log1p transform weights.
#' @param use_consensus Logical; enable consensus stage.
#' @param consensus_thr Numeric consensus threshold.
#' @param n_restart Integer number of restarts.
#' @param n_threads Integer thread count.
#' @param mode Character mode flag.
#' @param large_n_threshold Integer node threshold for aggressive mode.
#' @param nk_leiden_iterations Integer NetworKit iteration count.
#' @param nk_leiden_randomize Logical; randomize NetworKit initialisation.
#' @param aggr_future_workers Integer future workers.
#' @param aggr_batch_size Integer restart batch size.
#' @param enable_subcluster Logical; attempt Stage-2 subclustering.
#' @param sub_min_size Minimum Stage-2 cluster size.
#' @param sub_min_child_size Minimum Stage-2 child size.
#' @param sub_resolution_factor Multiplicative Stage-2 resolution factor.
#' @param sub_within_cons_max Maximum within-consensus for subclustering.
#' @param sub_conductance_min Minimum conductance for subclustering.
#' @param sub_improve_within_cons_min Minimum within-consensus gain.
#' @param sub_max_groups Maximum Stage-2 splits.
#' @param enable_qc_filter Logical; enable QC pruning.
#' @param qc_gene_intra_cons_min Minimum within-consensus per gene.
#' @param qc_gene_best_out_cons_min Maximum best-out consensus per gene.
#' @param qc_gene_intra_weight_q Quantile threshold for intra weights.
#' @param keep_cross_stable Logical; retain cross-cluster stable edges.
#' @param min_cluster_size Minimum cluster size to keep.
#' @param CI95_filter Logical; enable CI95 filtering.
#' @param CI_rule Character CI rule.
#' @param curve_layer Name of curve layer.
#' @param similarity_slot Similarity slot name.
#' @param significance_slot Significance slot name.
#' @param use_mh_weight Logical; apply Morisita-Horn weights.
#' @param mh_slot MH slot name.
#' @param post_smooth Logical; smooth post Stage-2 weights.
#' @param post_smooth_quant Numeric vector of smoothing quantiles.
#' @param post_smooth_power Numeric smoothing power.
#' @param keep_stage1_backbone Logical; keep Stage-1 backbone.
#' @param backbone_floor_q Quantile floor for backbone.
#'
#' @return A named configuration list consumed by downstream helpers.
#' @noRd
.build_pipeline_config <- function(
    min_cutoff = 0,
    use_significance = TRUE,
    significance_max = 0.05,
    pct_min = "q0",
    drop_isolated = TRUE,
    algo = "leiden",
    stage2_algo = NULL,
    resolution = 1,
    gamma = 1,
    objective = "CPM",
    use_log1p_weight = TRUE,
    use_consensus = TRUE,
    consensus_thr = 0.95,
    n_restart = 100,
    n_threads = NULL,
    mode = "auto",
    large_n_threshold = 1000L,
    nk_leiden_iterations = 10L,
    nk_leiden_randomize = TRUE,
    aggr_future_workers = 2L,
    aggr_batch_size = NULL,
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
    keep_cross_stable = TRUE,
    min_cluster_size = 2L,
    CI95_filter = FALSE,
    CI_rule = "remove_within",
    curve_layer = "LR_curve",
    similarity_slot = "L",
    significance_slot = "FDR",
    use_mh_weight = FALSE,
    mh_slot = "g_mh",
    post_smooth = TRUE,
    post_smooth_quant = c(0.05, 0.95),
    post_smooth_power = 0.5,
    keep_stage1_backbone = TRUE,
    backbone_floor_q = 0.02,
    hotspot_k = NULL,
    hotspot_min_module_size = 5L) {

    algo <- match.arg(algo, c("leiden", "louvain", "hotspot-like", "hotspot"))
    if (identical(algo, "hotspot")) algo <- "hotspot-like"
    if (!is.null(stage2_algo)) {
        stage2_algo <- match.arg(stage2_algo, c("leiden", "louvain", "hotspot-like", "hotspot"))
        if (identical(stage2_algo, "hotspot")) stage2_algo <- "hotspot-like"
    }

    list(
        min_cutoff = min_cutoff,
        use_significance = use_significance,
        significance_max = significance_max,
        pct_min = pct_min,
        drop_isolated = drop_isolated,
        algo = algo,
        stage2_algo = stage2_algo,
        resolution = resolution,
        gamma = gamma,
        objective = objective,
        use_log1p_weight = use_log1p_weight,
        use_consensus = use_consensus,
        consensus_thr = consensus_thr,
        n_restart = n_restart,
        n_threads = n_threads,
        mode = mode,
        large_n_threshold = as.integer(large_n_threshold),
        nk_leiden_iterations = nk_leiden_iterations,
        nk_leiden_randomize = nk_leiden_randomize,
        aggr_future_workers = aggr_future_workers,
        aggr_batch_size = aggr_batch_size,
        enable_subcluster = enable_subcluster,
        sub_min_size = sub_min_size,
        sub_min_child_size = sub_min_child_size,
        sub_resolution_factor = sub_resolution_factor,
        sub_within_cons_max = sub_within_cons_max,
        sub_conductance_min = sub_conductance_min,
        sub_improve_within_cons_min = sub_improve_within_cons_min,
        sub_max_groups = sub_max_groups,
        enable_qc_filter = enable_qc_filter,
        qc_gene_intra_cons_min = qc_gene_intra_cons_min,
        qc_gene_best_out_cons_min = qc_gene_best_out_cons_min,
        qc_gene_intra_weight_q = qc_gene_intra_weight_q,
        keep_cross_stable = keep_cross_stable,
        min_cluster_size = as.integer(min_cluster_size),
        CI95_filter = CI95_filter,
        CI_rule = match.arg(CI_rule, c("remove_within", "remove_outside")),
        curve_layer = curve_layer,
        similarity_slot = similarity_slot,
        significance_slot = significance_slot,
        use_mh_weight = use_mh_weight,
        mh_slot = mh_slot,
        post_smooth = post_smooth,
        post_smooth_quant = post_smooth_quant,
        post_smooth_power = post_smooth_power,
        keep_stage1_backbone = keep_stage1_backbone,
        backbone_floor_q = backbone_floor_q,
        hotspot_k = if (is.null(hotspot_k)) NULL else as.integer(hotspot_k),
        hotspot_min_module_size = as.integer(hotspot_min_module_size)
    )
}

# ---------------------------------------------------------------------------
# Preparation: extract matrices/aux data from container object
# ---------------------------------------------------------------------------
#' Extract primary matrices from a CosMx scope object.
#'
#' @param scope_obj CosMx scope object.
#' @param grid_name Character grid identifier.
#' @param stats_layer Statistics layer name.
#' @param similarity_slot Slot holding the similarity matrix.
#' @param significance_slot Slot holding the significance matrix.
#' @param use_significance Logical; whether to return significance matrix.
#' @param verbose Logical; emit progress messages.
#'
#' @return A list containing similarity and significance matrices plus
#'   metadata handles required by downstream stages.
#' @noRd
.extract_scope_layers <- function(scope_obj,
                                        grid_name,
                                        stats_layer = "LeeStats_Xz",
                                        similarity_slot = "L",
                                        significance_slot = "FDR",
                                        use_significance = TRUE,
                                        verbose = TRUE) {
    if (isTRUE(verbose)) message("[cluster] Prep: extracting container layers…")

    g_layer <- .selectGridLayer(scope_obj, grid_name)
    gname <- names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    .checkGridContent(scope_obj, gname)

    if (is.null(scope_obj@stats[[gname]]) ||
        is.null(scope_obj@stats[[gname]][[stats_layer]])) {
        stop("Statistics layer missing: ", grid_name, "/", stats_layer)
    }

    stats_obj <- scope_obj@stats[[gname]][[stats_layer]]
    sim_mat <- stats_obj[[similarity_slot]]
    if (is.null(sim_mat)) stop("Similarity slot missing in stats layer: ", similarity_slot)
    sig_mat <- if (use_significance) stats_obj[[significance_slot]] else NULL

    list(
        grid_name = gname,
        stats_layer = stats_layer,
        similarity = sim_mat,
        significance = sig_mat,
        stats_object = stats_obj,
        scope_stats = scope_obj@stats[[gname]],
        L_raw = sim_mat,          # backward compatibility
        FDR = sig_mat,            # backward compatibility
        lee_stats = stats_obj     # backward compatibility
    )
}

# ---------------------------------------------------------------------------
# Shared helpers reused across stages
# ---------------------------------------------------------------------------
#' Run community detection with graceful fallbacks.
#'
#' @param g igraph object.
#' @param algo_name Character algorithm name.
#' @param res_param Numeric resolution parameter.
#' @param objective Objective string.
#'
#' @return An igraph clustering object.
#' @noRd
.run_graph_algorithm <- function(g, algo_name, res_param, objective) {
    w <- igraph::E(g)$weight
    if (identical(algo_name, "leiden")) {
        res <- try(igraph::cluster_leiden(
            g,
            weights = w,
            resolution_parameter = res_param,
            objective_function = objective
        ), silent = TRUE)
        if (inherits(res, "try-error")) {
            res <- try(igraph::cluster_leiden(
                g,
                weights = w,
                resolution = res_param
            ), silent = TRUE)
            if (inherits(res, "try-error")) res <- igraph::cluster_leiden(g, weights = w)
        }
        return(res)
    }
    formn <- try(names(formals(igraph::cluster_louvain)), silent = TRUE)
    if (!inherits(formn, "try-error") && "resolution" %in% formn) {
        return(igraph::cluster_louvain(g, weights = w, resolution = res_param))
    }
    igraph::cluster_louvain(g, weights = w)
}

`%||%` <- function(x, y) {
    if (!is.null(x)) x else y
}

.run_single_partition <- function(graph,
                                  backend = c("igraph", "networkit"),
                                  algo_name = c("leiden", "louvain"),
                                  res_param,
                                  objective,
                                  threads = 1L,
                                  nk_opts = list()) {
    backend <- match.arg(backend)
    algo_name <- match.arg(algo_name)
    vertex_names <- igraph::V(graph)$name

    if (backend == "networkit") {
        edge_df <- igraph::as_data_frame(graph, what = "edges")[, c("from", "to"), drop = FALSE]
        weights <- igraph::E(graph)$weight
        seeds <- nk_opts$seeds %||% 1L
        if (identical(algo_name, "leiden")) {
            memb_mat <- .networkit_parallel_leiden_runs(
                vertex_names = vertex_names,
                edge_table = edge_df,
                edge_weights = weights,
                gamma = as.numeric(res_param),
                iterations = nk_opts$iterations %||% 10L,
                randomize = nk_opts$randomize %||% TRUE,
                threads = threads,
                seeds = as.integer(seeds)
            )
        } else {
            memb_mat <- .networkit_plm_runs(
                vertex_names = vertex_names,
                edge_table = edge_df,
                edge_weights = weights,
                gamma = as.numeric(res_param),
                refine = nk_opts$refine %||% TRUE,
                threads = threads,
                seeds = as.integer(seeds)
            )
        }
        membership <- as.integer(memb_mat[, 1])
        names(membership) <- vertex_names
        return(membership)
    }

    clustering <- .run_graph_algorithm(graph, algo_name, res_param, objective)
    membership <- igraph::membership(clustering)
    if (is.null(names(membership))) names(membership) <- vertex_names
    membership
}

#' Compute consensus matrix in pure R (sparse).
#'
#' @param memb Integer matrix of memberships (genes x runs).
#' @param thr Numeric consensus threshold.
#'
#' @return Symmetric sparseMatrix of mean co-occurrence scores.
#' @noRd
.compute_consensus_sparse_R <- function(memb, thr = 0) {
    if (is.null(dim(memb))) memb <- matrix(as.integer(memb), nrow = length(memb))
    if (!is.matrix(memb)) memb <- as.matrix(memb)
    storage.mode(memb) <- "integer"
    n <- nrow(memb); r <- ncol(memb)
    if (!is.finite(r) || r < 1L) return(Matrix::Matrix(0, n, n, sparse = TRUE))
    acc <- Matrix::Matrix(0, n, n, sparse = FALSE)
    for (jj in seq_len(r)) {
        labs <- as.integer(memb[, jj])
        keep <- !is.na(labs)
        if (!any(keep)) next
        labs2 <- labs[keep]
        K <- suppressWarnings(max(labs2))
        if (!is.finite(K) || K < 1L) next
        ind <- Matrix::sparseMatrix(i = which(keep), j = labs2, x = 1, dims = c(n, K))
        acc <- acc + as.matrix(ind %*% Matrix::t(ind))
    }
    acc <- acc / max(1L, r)
    diag(acc) <- 0
    M <- Matrix::drop0(Matrix::Matrix(acc, sparse = TRUE))
    rn <- rownames(memb)
    if (!is.null(rn) && length(rn) == nrow(M)) {
        try({ dimnames(M) <- list(rn, rn) }, silent = TRUE)
    }
    if (thr > 0) {
        TT <- methods::as(M, "TsparseMatrix")
        if (length(TT@x)) {
            keep <- TT@x >= thr
            M <- Matrix::sparseMatrix(i = TT@i[keep] + 1L, j = TT@j[keep] + 1L, x = TT@x[keep], dims = dim(M), symmetric = TRUE)
            try({ dimnames(M) <- list(rn, rn) }, silent = TRUE)
        }
    }
    M
}

# ---------------------------------------------------------------------------
# Metrics helpers reused across Stage-2 and downstream diagnostics
# ---------------------------------------------------------------------------
#' Compute cluster-level diagnostics from consensus and weight matrices.
#'
#' @param memberships Integer vector of cluster assignments (NA allowed).
#' @param consensus_matrix Square matrix of consensus weights (gene x gene).
#' @param weight_matrix Square matrix of edge weights (gene x gene).
#'
#' @return A data.frame with one row per cluster containing columns
#'   `cluster`, `size`, `within_cons`, and `conductance`.
#' @noRd
.summarise_cluster_metrics <- function(memberships,
                                    consensus_matrix,
                                    weight_matrix) {
    if (!length(memberships)) {
        return(data.frame(
            cluster = integer(0),
            size = integer(0),
            within_cons = numeric(0),
            conductance = numeric(0)
        ))
    }
    cl_ids <- sort(na.omit(unique(memberships)))
    if (!length(cl_ids)) {
        return(data.frame(
            cluster = integer(0),
            size = integer(0),
            within_cons = numeric(0),
            conductance = numeric(0)
        ))
    }
    res <- lapply(cl_ids, function(cid) {
        idx <- which(memberships == cid)
        n <- length(idx)
        if (n <= 1) {
            return(data.frame(
                cluster = cid,
                size = n,
                within_cons = NA_real_,
                conductance = NA_real_
            ))
        }
        cons_block <- as.matrix(consensus_matrix[idx, idx])
        within_cons <- mean(cons_block[upper.tri(cons_block)], na.rm = TRUE)
        W_block <- as.matrix(weight_matrix[idx, idx])
        Win <- sum(W_block) / 2
        Wcut <- sum(weight_matrix[idx, -idx, drop = FALSE])
        conductance <- Wcut / (Wcut + Win + 1e-12)
        data.frame(
            cluster = cid,
            size = n,
            within_cons = within_cons,
            conductance = conductance
        )
    })
    do.call(rbind, res)
}

#' Compute per-gene membership diagnostics within clusters.
#'
#' @param memberships Integer vector of cluster assignments (NA allowed).
#' @param consensus_matrix Square consensus matrix (gene x gene).
#' @param weight_matrix Square weight matrix (gene x gene).
#' @param genes Character vector of gene identifiers corresponding to rows.
#'
#' @return A data.frame with columns `gene`, `cluster`, `p_in`, `p_best_out`,
#'   and `w_in`.
#' @noRd
.summarise_gene_membership <- function(memberships,
                                         consensus_matrix,
                                         weight_matrix,
                                         genes = NULL) {
    n <- length(memberships)
    if (is.null(genes)) {
        genes <- rownames(consensus_matrix)
        if (is.null(genes)) genes <- as.character(seq_len(n))
    }
    if (length(genes) != n) {
        stop("Length of `genes` must match length of `memberships`.", call. = FALSE)
    }
    cl_ids <- sort(na.omit(unique(memberships)))
    idx_list <- lapply(cl_ids, function(cid) which(memberships == cid))
    names(idx_list) <- as.character(cl_ids)
    p_in <- numeric(n); p_best_out <- numeric(n); w_in <- numeric(n)
    for (i in seq_len(n)) {
        cid <- memberships[i]
        if (is.na(cid)) {
            p_in[i] <- NA_real_
            p_best_out[i] <- NA_real_
            w_in[i] <- NA_real_
            next
        }
        idx <- idx_list[[as.character(cid)]]
        idx_noi <- setdiff(idx, i)
        if (length(idx_noi)) {
            p_in[i] <- mean(consensus_matrix[i, idx_noi], na.rm = TRUE)
            w_in[i] <- sum(weight_matrix[i, idx_noi])
        } else {
            p_in[i] <- 0
            w_in[i] <- 0
        }
        other <- setdiff(cl_ids, cid)
        if (length(other)) {
            means <- vapply(other, function(cj) {
                jj <- idx_list[[as.character(cj)]]
                if (length(jj)) mean(consensus_matrix[i, jj], na.rm = TRUE) else 0
            }, numeric(1))
            p_best_out[i] <- max(means, na.rm = TRUE)
        } else {
            p_best_out[i] <- 0
        }
    }
    data.frame(
        gene = genes,
        cluster = memberships,
        p_in = p_in,
        p_best_out = p_best_out,
        w_in = w_in,
        stringsAsFactors = FALSE
    )
}

#' Prune unsupported or too-small clusters after Stage-2 correction.
#'
#' @param memberships Integer vector of cluster assignments (NA allowed).
#' @param L_matrix Square matrix of corrected Lee's L weights (gene x gene).
#' @param min_cluster_size Minimum cluster size to retain (default 2).
#'
#' @return A list with elements `membership` (pruned vector) and
#'   `stable_clusters` (character vector of surviving cluster ids).
#' @noRd
.prune_unstable_memberships <- function(memberships,
                                      L_matrix,
                                      min_cluster_size = 2L) {
    if (!length(memberships)) {
        return(list(membership = memberships, stable_clusters = character(0)))
    }
    if (is.null(dim(L_matrix)) ||
        nrow(L_matrix) != length(memberships) ||
        ncol(L_matrix) != length(memberships)) {
        stop("`L_matrix` must be square and aligned with `memberships`.", call. = FALSE)
    }
    min_cluster_size <- max(1L, as.integer(min_cluster_size))
    logical_L <- L_matrix > 0
    if (!is.matrix(logical_L)) logical_L <- as.matrix(logical_L)
    same_cluster <- outer(memberships, memberships, "==")
    di <- rowSums(logical_L & same_cluster, na.rm = TRUE)
    memberships[di == 0] <- NA
    tbl <- table(memberships, useNA = "no")
    stable_clusters <- character(0)
    if (length(tbl)) {
        if (min_cluster_size > 1L) {
            small <- names(tbl)[tbl < min_cluster_size]
            if (length(small)) {
                memb_chr <- as.character(memberships)
                memberships[memb_chr %in% small] <- NA
            }
        }
        tbl <- table(memberships, useNA = "no")
        stable_clusters <- names(tbl)[tbl >= min_cluster_size]
    }
    list(membership = memberships, stable_clusters = stable_clusters)
}

# ---------------------------------------------------------------------------
# Hotspot-style hierarchical clustering helpers
# ---------------------------------------------------------------------------
#' Merge small hotspot clusters into nearest neighbours.
#'
#' @param membership Integer vector of cluster assignments.
#' @param similarity_dense Dense similarity matrix.
#' @param min_size Minimum allowed cluster size.
#'
#' @return Updated membership vector.
#' @noRd
.merge_hotspot_clusters <- function(membership, similarity_dense, min_size) {
    membership <- as.integer(membership)
    names(membership) <- rownames(similarity_dense)
    repeat {
        freq <- table(membership, useNA = "no")
        small_clusters <- as.integer(names(freq)[freq < min_size])
        if (!length(small_clusters)) break
        small_clusters <- small_clusters[order(as.integer(freq[as.character(small_clusters)]))]
        changed <- FALSE
        for (cid in small_clusters) {
            current_idx <- which(membership == cid)
            other_ids <- setdiff(as.integer(names(freq)), cid)
            if (!length(other_ids)) next
            best_target <- NA_integer_
            best_value <- -Inf
            for (target in other_ids) {
                target_idx <- which(membership == target)
                if (!length(target_idx)) next
                avg_sim <- mean(similarity_dense[current_idx, target_idx, drop = FALSE], na.rm = TRUE)
                if (!is.finite(avg_sim)) avg_sim <- 0
                if (avg_sim > best_value) {
                    best_value <- avg_sim
                    best_target <- target
                }
            }
            if (!is.na(best_target)) {
                membership[current_idx] <- best_target
                changed <- TRUE
            }
        }
        if (!changed) break
    }
    unique_ids <- sort(na.omit(unique(membership)))
    if (length(unique_ids)) {
        id_map <- setNames(seq_along(unique_ids), unique_ids)
        idx <- !is.na(membership)
        membership[idx] <- id_map[as.character(membership[idx])]
    }
    membership
}

#' Build a binary co-membership consensus matrix.
#'
#' @param membership Integer vector of cluster labels (NA allowed).
#' @param gene_names Gene identifiers.
#'
#' @return Symmetric sparseMatrix with 1 for within-cluster pairs.
#' @noRd
.build_membership_consensus <- function(membership, gene_names) {
    n <- length(membership)
    if (n == 0) {
        return(Matrix::Matrix(0, 0, 0, sparse = TRUE))
    }
    membership <- as.integer(membership)
    groups <- split(seq_len(n), membership)
    ii <- jj <- integer(0)
    for (idx in groups) {
        idx <- idx[!is.na(idx)]
        s <- length(idx)
        if (s > 1) {
            cmb <- utils::combn(idx, 2)
            ii <- c(ii, cmb[1, ])
            jj <- c(jj, cmb[2, ])
        }
    }
    consensus <- Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, n), symmetric = TRUE)
    Matrix::diag(consensus) <- 0
    consensus <- Matrix::drop0(consensus)
    if (!is.null(gene_names) && length(gene_names) == n) {
        dimnames(consensus) <- list(gene_names, gene_names)
    }
    consensus
}

#' Perform hotspot-like hierarchical clustering on a similarity matrix.
#'
#' @param similarity_matrix Sparse similarity matrix (already thresholded).
#' @param cluster_count Desired cluster count.
#' @param min_module_size Minimum module size after merge.
#' @param apply_log1p Whether to apply log1p before normalisation.
#'
#' @return List with `membership` and `consensus` entries.
#' @noRd
.run_hotspot_clustering <- function(similarity_matrix,
                                    cluster_count,
                                    min_module_size,
                                    apply_log1p = TRUE) {
    if (cluster_count <= 0) cluster_count <- 1L
    if (!inherits(similarity_matrix, "Matrix") && !is.matrix(similarity_matrix)) {
        similarity_matrix <- as.matrix(similarity_matrix)
    }
    matrix_work <- similarity_matrix
    if (apply_log1p && inherits(matrix_work, "sparseMatrix")) {
        matrix_work <- methods::as(matrix_work, "dgCMatrix")
        if (length(matrix_work@x)) matrix_work@x <- log1p(matrix_work@x)
    } else if (apply_log1p) {
        matrix_work <- log1p(matrix_work)
    }
    dense_matrix <- if (inherits(matrix_work, "Matrix")) as.matrix(matrix_work) else matrix_work
    if (!is.double(dense_matrix)) dense_matrix <- matrix(as.numeric(dense_matrix), nrow = nrow(dense_matrix), dimnames = dimnames(dense_matrix))
    max_val <- suppressWarnings(max(dense_matrix, na.rm = TRUE))
    if (!is.finite(max_val) || max_val <= 0) max_val <- 1
    dense_matrix <- dense_matrix / max_val
    diag(dense_matrix) <- 1
    gene_names <- rownames(dense_matrix)
    cluster_count <- max(1L, min(as.integer(cluster_count), nrow(dense_matrix)))
    dist_matrix <- stats::as.dist(1 - dense_matrix)
    hc <- try(if (requireNamespace("fastcluster", quietly = TRUE)) fastcluster::hclust(dist_matrix, method = "ward.D2") else stats::hclust(dist_matrix, method = "ward.D2"), silent = TRUE)
    if (inherits(hc, "try-error")) hc <- stats::hclust(dist_matrix, method = "average")
    membership <- stats::cutree(hc, k = cluster_count)
    if (is.null(names(membership))) names(membership) <- gene_names
    if (!is.null(min_module_size) && is.finite(min_module_size) && min_module_size > 1) {
        membership <- .merge_hotspot_clusters(membership, dense_matrix, as.integer(min_module_size))
    }
    consensus <- .build_membership_consensus(membership, gene_names)
    list(membership = membership, consensus = consensus)
}

# ---------------------------------------------------------------------------
# NetworKit integration helpers
# ---------------------------------------------------------------------------
#' Internal state for NetworKit runtime configuration.
#' @noRd
.networkit_runtime_state <- new.env(parent = emptyenv())

#' Configure reticulate to use a specific Python environment for NetworKit.
#'
#' @param conda_env Optional conda environment name.
#' @param conda_bin Optional explicit conda binary path.
#' @param python_path Optional absolute python interpreter path.
#' @param required Whether failure should raise.
#'
#' @return Logical invisibly indicating whether configuration was attempted.
#' @noRd
.configure_networkit_python <- function(conda_env = NULL,
                                        conda_bin = NULL,
                                        python_path = NULL,
                                        required = FALSE) {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(FALSE))

    if (!is.null(python_path) && nzchar(python_path)) {
        .use_python_interpreter(python_path, required = required)
        return(invisible(TRUE))
    }

    if (!is.null(conda_env) && nzchar(conda_env)) {
        .use_conda_environment(conda_env, conda_bin = conda_bin, required = required)
        return(invisible(TRUE))
    }

    # fallback: try CONDA_PREFIX
    conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
    if (nzchar(conda_prefix)) {
        default_python <- file.path(conda_prefix, "bin", "python")
        try(.use_python_interpreter(default_python, required = FALSE), silent = TRUE)
        return(invisible(TRUE))
    }
    invisible(FALSE)
}

#' Request a specific Python interpreter for NetworKit via reticulate.
#'
#' @param python_path Absolute path to python executable.
#' @param required Whether to error if binding fails.
#'
#' @return Invisible TRUE when successful.
#' @noRd
.use_python_interpreter <- function(python_path, required = FALSE) {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(FALSE))
    reticulate::use_python(python_path, required = required)
    assign("configured", TRUE, envir = .networkit_runtime_state)
    invisible(TRUE)
}

#' Request a conda environment for NetworKit via reticulate.
#'
#' @param env_name Conda environment name.
#' @param conda_bin Optional explicit conda binary path.
#' @param required Whether to error if binding fails.
#'
#' @return Invisible TRUE when successful.
#' @noRd
.use_conda_environment <- function(env_name, conda_bin = NULL, required = FALSE) {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(FALSE))
    if (is.null(conda_bin)) {
        reticulate::use_condaenv(env_name, required = required)
    } else {
        reticulate::use_condaenv(env_name, conda = conda_bin, required = required)
    }
    assign("configured", TRUE, envir = .networkit_runtime_state)
    invisible(TRUE)
}

#' Request a virtualenv for NetworKit via reticulate.
#'
#' @param env_name Virtualenv name.
#' @param required Whether to error if binding fails.
#'
#' @return Invisible TRUE when successful.
#' @noRd
.use_virtual_environment <- function(env_name, required = FALSE) {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(FALSE))
    reticulate::use_virtualenv(env_name, required = required)
    assign("configured", TRUE, envir = .networkit_runtime_state)
    invisible(TRUE)
}

#' Autoconfigure NetworKit runtime from options/env vars.
#'
#' @return Logical invisibly indicating whether configuration was attempted.
#' @noRd
.autoconfigure_networkit_runtime <- function() {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(invisible(FALSE))

    python_opt <- getOption("clusterGenes.python", default = getOption("geneSCOPE.python", default = NULL))
    if (is.character(python_opt) && nzchar(python_opt)) {
        try(.use_python_interpreter(python_opt, required = FALSE), silent = TRUE)
    } else {
        conda_opt <- getOption("clusterGenes.condaenv", default = getOption("geneSCOPE.condaenv", default = NULL))
        if (is.character(conda_opt) && nzchar(conda_opt)) {
            try(.use_conda_environment(conda_opt, required = FALSE), silent = TRUE)
        }
        venv_opt <- getOption("clusterGenes.virtualenv", default = getOption("geneSCOPE.virtualenv", default = NULL))
        if (is.character(venv_opt) && nzchar(venv_opt)) {
            try(.use_virtual_environment(venv_opt, required = FALSE), silent = TRUE)
        }
    }

    env_python <- Sys.getenv("NETWORKIT_PYTHON", unset = Sys.getenv("NK_PYTHON", unset = Sys.getenv("RETICULATE_PYTHON", unset = "")))
    if (nzchar(env_python)) {
        try(.use_python_interpreter(env_python, required = FALSE), silent = TRUE)
    }
    invisible(TRUE)
}

#' Check whether NetworKit is available via reticulate.
#'
#' @return TRUE when NetworKit Python module can be imported.
#' @noRd
.networkit_available <- function() {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
    if (!isTRUE(get0("configured", envir = .networkit_runtime_state, ifnotfound = FALSE))) {
        .autoconfigure_networkit_runtime()
        assign("configured", TRUE, envir = .networkit_runtime_state)
    }
    ok <- tryCatch({ reticulate::py_module_available("networkit") }, error = function(e) FALSE)
    isTRUE(ok)
}

#' Run NetworKit Parallel Leiden across multiple seeds.
#'
#' @param vertex_names Character vector of vertex names.
#' @param edge_table Data frame with columns `from`, `to`.
#' @param edge_weights Numeric vector of edge weights.
#' @param gamma Resolution parameter for Leiden.
#' @param iterations Number of parallel Leiden iterations.
#' @param randomize Logical; randomize the solver.
#' @param threads Thread count.
#' @param seeds Integer vector of seeds.
#'
#' @return Integer matrix of memberships (rows = vertices, cols = seeds).
#' @noRd
.networkit_parallel_leiden_runs <- function(vertex_names, edge_table, edge_weights,
                                            gamma = 1.0, iterations = 10L,
                                            randomize = TRUE,
                                            threads = 1L, seeds = 1L) {
    if (!.networkit_available()) {
        g_local <- igraph::graph_from_data_frame(cbind(edge_table, weight = edge_weights), directed = FALSE, vertices = vertex_names)
        membership_list <- lapply(seeds, function(seed_id) {
            res <- try(igraph::cluster_leiden(g_local, weights = igraph::E(g_local)$weight, resolution_parameter = gamma), silent = TRUE)
            if (inherits(res, "try-error")) res <- try(igraph::cluster_leiden(g_local, weights = igraph::E(g_local)$weight, resolution = gamma), silent = TRUE)
            if (inherits(res, "try-error")) res <- igraph::cluster_leiden(g_local, weights = igraph::E(g_local)$weight)
            memberships <- igraph::membership(res)
            if (is.null(names(memberships))) names(memberships) <- igraph::V(g_local)$name
            as.integer(memberships[vertex_names])
        })
        membership_matrix <- do.call(cbind, membership_list)
        rownames(membership_matrix) <- vertex_names
        return(membership_matrix)
    }

    nk <- reticulate::import("networkit", delay_load = TRUE)
    threads <- as.integer(threads)
    if (!is.na(threads) && threads > 0) try(nk$setNumberOfThreads(threads), silent = TRUE)

    vertex_count <- length(vertex_names)
    nk_graph <- nk$Graph(as.integer(vertex_count), weighted = TRUE, directed = FALSE)
    if (nrow(edge_table)) {
        vertex_ids <- seq_len(vertex_count) - 1L
        name_to_id <- structure(vertex_ids, names = vertex_names)
        u <- name_to_id[as.character(edge_table$from)]
        v <- name_to_id[as.character(edge_table$to)]
        w <- as.numeric(edge_weights)
        keep <- !is.na(u) & !is.na(v) & is.finite(w)
        if (any(keep)) {
            u <- as.integer(u[keep]); v <- as.integer(v[keep]); w <- as.numeric(w[keep])
            for (idx in seq_along(u)) nk_graph$addEdge(u[[idx]], v[[idx]], w[[idx]])
        }
    }

    seeds <- as.integer(seeds)
    membership_list <- vector("list", length(seeds))
    for (i in seq_along(seeds)) {
        if (!is.na(seeds[i])) try(nk$random$setSeed(seeds[i], TRUE), silent = TRUE)
        pl <- nk$community$ParallelLeiden(nk_graph,
            iterations = as.integer(iterations),
            randomize = as.logical(randomize),
            gamma = as.numeric(gamma))
        pl$run()
        partition <- pl$getPartition()
        memberships_zero_based <- as.integer(unlist(partition$getVector()))
        membership_list[[i]] <- memberships_zero_based + 1L
    }
    membership_matrix <- do.call(cbind, membership_list)
    rownames(membership_matrix) <- vertex_names
    colnames(membership_matrix) <- paste0("r", seq_along(seeds))
    membership_matrix
}

#' Run NetworKit PLM (multi-run) across seeds.
#'
#' @param vertex_names Character vector of vertex names.
#' @param edge_table Data frame with columns `from`, `to`.
#' @param edge_weights Numeric vector of edge weights.
#' @param gamma Modularity resolution parameter.
#' @param refine Whether to enable PLM refinement.
#' @param threads Thread count.
#' @param seeds Integer vector of seeds.
#'
#' @return Integer matrix of memberships.
#' @noRd
.networkit_plm_runs <- function(vertex_names, edge_table, edge_weights,
                                gamma = 1.0, refine = TRUE,
                                threads = 1L, seeds = 1L) {
    if (!.networkit_available()) {
        g_local <- igraph::graph_from_data_frame(cbind(edge_table, weight = edge_weights), directed = FALSE, vertices = vertex_names)
        membership_list <- lapply(seeds, function(seed_id) {
            res <- try(igraph::cluster_louvain(g_local, weights = igraph::E(g_local)$weight, resolution = gamma), silent = TRUE)
            if (inherits(res, "try-error")) res <- igraph::cluster_louvain(g_local, weights = igraph::E(g_local)$weight)
            memberships <- igraph::membership(res)
            if (is.null(names(memberships))) names(memberships) <- igraph::V(g_local)$name
            as.integer(memberships[vertex_names])
        })
        membership_matrix <- do.call(cbind, membership_list)
        rownames(membership_matrix) <- vertex_names
        return(membership_matrix)
    }

    nk <- reticulate::import("networkit", delay_load = TRUE)
    threads <- as.integer(threads)
    if (!is.na(threads) && threads > 0) try(nk$setNumberOfThreads(threads), silent = TRUE)

    vertex_count <- length(vertex_names)
    nk_graph <- nk$Graph(as.integer(vertex_count), weighted = TRUE, directed = FALSE)
    if (nrow(edge_table)) {
        vertex_ids <- seq_len(vertex_count) - 1L
        name_to_id <- structure(vertex_ids, names = vertex_names)
        u <- name_to_id[as.character(edge_table$from)]
        v <- name_to_id[as.character(edge_table$to)]
        w <- as.numeric(edge_weights)
        keep <- !is.na(u) & !is.na(v) & is.finite(w)
        if (any(keep)) {
            u <- as.integer(u[keep]); v <- as.integer(v[keep]); w <- as.numeric(w[keep])
            for (idx in seq_along(u)) nk_graph$addEdge(u[[idx]], v[[idx]], w[[idx]])
        }
    }

    seeds <- as.integer(seeds)
    membership_list <- vector("list", length(seeds))
    for (i in seq_along(seeds)) {
        if (!is.na(seeds[i])) try(nk$random$setSeed(seeds[i], TRUE), silent = TRUE)
        plm_algo <- nk$community$PLM(nk_graph, refine = as.logical(refine), gamma = as.numeric(gamma))
        plm_algo$run()
        partition <- plm_algo$getPartition()
        memberships_zero_based <- as.integer(unlist(partition$getVector()))
        membership_list[[i]] <- memberships_zero_based + 1L
    }
    membership_matrix <- do.call(cbind, membership_list)
    rownames(membership_matrix) <- vertex_names
    colnames(membership_matrix) <- paste0("r", seq_along(seeds))
    membership_matrix
}

# Legacy aliases retained for backward compatibility with older scripts.
.nk_available <- .networkit_available
.nk_parallel_leiden_mruns <- .networkit_parallel_leiden_runs
.nk_plm_mruns <- .networkit_plm_runs

# ---------------------------------------------------------------------------
# Step 1: Threshold Lee's L matrix and build base adjacency
# ---------------------------------------------------------------------------
#' Threshold similarity matrix and prepare adjacency structures.
#'
#' @param similarity_matrix Square matrix of similarity scores.
#' @param significance_matrix Optional matching significance matrix.
#' @param config Configuration list produced by `.build_pipeline_config()`.
#' @param verbose Logical; emit progress messages.
#'
#' @return A list containing filtered adjacency matrices and gene sets.
#' @noRd
.threshold_similarity_graph <- function(similarity_matrix,
                                  significance_matrix = NULL,
                                  config,
                                  verbose = TRUE) {
    if (isTRUE(verbose)) message("[cluster] Step1: thresholding similarity matrix…")

    genes_all <- rownames(similarity_matrix)
    if (is.null(genes_all)) stop("Input similarity matrix must have rownames/colnames for node IDs.")

    similarity_filtered <- similarity_matrix
    similarity_filtered[abs(similarity_filtered) < config$min_cutoff] <- 0
    if (config$use_significance && !is.null(significance_matrix)) {
        similarity_filtered <- tryCatch(
            .align_and_filter_FDR(similarity_filtered, similarity_matrix, significance_matrix, config$significance_max),
            error = function(e) {
                warning("[cluster] significance filter disabled: ", conditionMessage(e))
                similarity_filtered
            }
        )
    }
    similarity_filtered[similarity_filtered < 0] <- 0
    similarity_filtered <- pmax(similarity_filtered, t(similarity_filtered))
    diag(similarity_filtered) <- 0
    similarity_filtered <- .filter_matrix_by_quantile(similarity_filtered, config$pct_min, "q100")
    if (!inherits(similarity_filtered, "sparseMatrix")) {
        similarity_filtered <- Matrix::drop0(Matrix::Matrix(similarity_filtered, sparse = TRUE))
    }

    keep <- if (config$drop_isolated) Matrix::rowSums(similarity_filtered != 0) > 0 else rep(TRUE, nrow(similarity_filtered))
    kept_genes <- genes_all[keep]
    if (length(kept_genes) < 2) stop("Fewer than 2 genes after filtering.")

    similarity_kept <- similarity_filtered[kept_genes, kept_genes, drop = FALSE]
    if (isTRUE(verbose)) {
        message(sprintf("[cluster]   retained %d/%d nodes; |E|=%d",
                        length(kept_genes), length(genes_all), Matrix::nnzero(similarity_kept) / 2))
    }

    list(
        similarity = similarity_filtered,
        similarity_sub = similarity_kept,
        kept_genes = kept_genes,
        genes_all = genes_all,
        A = similarity_filtered,                 # backward compatibility
        A_sub = similarity_kept           # backward compatibility
    )
}

# ---------------------------------------------------------------------------
# Step 2: Stage-1 multi-run community detection + consensus
# ---------------------------------------------------------------------------
#' Run Stage-1 multi-restart clustering and build consensus.
#'
#' @param similarity_matrix Thresholded similarity matrix.
#' @param threshold Output list from `.threshold_similarity_graph()`.
#' @param config Configuration list.
#' @param verbose Logical; emit progress messages.
#' @param networkit_config Optional list describing NetworKit runtime configuration.
#'
#' @return A list containing Stage-1 memberships, consensus matrix, and graphs.
#' @noRd
.stage1_consensus_workflow <- function(similarity_matrix,
                                      threshold,
                                      config,
                                      verbose = TRUE,
                                      networkit_config = list()) {
    if (isTRUE(verbose)) message("[cluster] Step2: running Stage-1 consensus…")

    kept_genes <- threshold$kept_genes
    similarity_sub <- threshold$similarity_sub
    genes_all <- threshold$genes_all

    if (!is.null(networkit_config) && length(networkit_config)) {
        config_conda_env <- networkit_config[["conda_env"]]
        if (is.null(config_conda_env)) config_conda_env <- networkit_config[["nk_condaenv"]]
        config_conda_bin <- networkit_config[["conda_bin"]]
        if (is.null(config_conda_bin)) config_conda_bin <- networkit_config[["nk_conda_bin"]]
        config_python_path <- networkit_config[["python_path"]]
        if (is.null(config_python_path)) config_python_path <- networkit_config[["nk_python"]]
        try(.configure_networkit_python(
            conda_env = config_conda_env,
            conda_bin = config_conda_bin,
            python_path = config_python_path,
            required = FALSE
        ), silent = TRUE)
    }

    adjacency_triplet <- as.data.frame(summary(similarity_sub))
    if (!all(c("i", "j", "x") %in% colnames(adjacency_triplet))) {
        colnames(adjacency_triplet) <- c("i", "j", "x")[seq_len(ncol(adjacency_triplet))]
    }
    adjacency_triplet <- adjacency_triplet[adjacency_triplet$i < adjacency_triplet$j & adjacency_triplet$x > 0, , drop = FALSE]
    if (!nrow(adjacency_triplet)) stop("All edges filtered out at Stage-1.")

    edge_table <- data.frame(
        from = kept_genes[adjacency_triplet$i],
        to = kept_genes[adjacency_triplet$j],
        sim_edge = adjacency_triplet$x,
        stringsAsFactors = FALSE
    )
    edge_table$weight <- if (config$use_log1p_weight) log1p(edge_table$sim_edge) else edge_table$sim_edge
    similarity_graph <- igraph::graph_from_data_frame(edge_table[, c("from", "to", "weight")],
        directed = FALSE, vertices = kept_genes
    )

    algo <- match.arg(config$algo, c("leiden", "louvain", "hotspot-like"))
    objective <- match.arg(config$objective, c("CPM", "modularity"))
    mode_selected <- switch(config$mode,
        auto = if (length(kept_genes) > config$large_n_threshold) "aggressive" else "safe",
        safe = "safe",
        aggressive = "aggressive"
    )
    if (isTRUE(verbose)) {
        message(sprintf("[cluster]   Stage-1 mode=%s (nodes=%d)", mode_selected, length(kept_genes)))
    }
    if (identical(mode_selected, "aggressive") && identical(objective, "CPM")) {
        warning("[cluster] CPM objective not supported in aggressive mode; falling back to 'modularity'.")
        objective <- "modularity"
    }

    stage1_consensus_matrix <- NULL
    stage1_membership_labels <- NULL
    restart_memberships <- NULL
    stage1_backend <- "igraph"
    stage1_backend <- "igraph"

    if (identical(algo, "hotspot-like")) {
        if (isTRUE(verbose)) message("[cluster]   Stage-1 backend: hotspot-like")
        hotspot_k <- config$hotspot_k
        if (is.null(hotspot_k) || !is.finite(hotspot_k)) {
            hotspot_k <- max(1L, min(20L, length(kept_genes)))
        }
        hotspot_result <- .run_hotspot_clustering(similarity_sub, hotspot_k, config$hotspot_min_module_size, config$use_log1p_weight)
        stage1_membership_labels <- hotspot_result$membership
        stage1_consensus_matrix <- hotspot_result$consensus
        restart_memberships <- matrix(as.integer(stage1_membership_labels), ncol = 1,
            dimnames = list(kept_genes, "hotspot"))
        stage1_backend <- "hotspot"
    } else {
        n_threads <- tryCatch({
            if (is.null(config$n_threads)) {
                if (exists("safe_thread_count", mode = "function", inherits = TRUE)) safe_thread_count() else 1L
            } else {
                as.integer(config$n_threads)
            }
        }, error = function(e) 1L)

        edge_weights <- igraph::E(similarity_graph)$weight
        vertex_names <- igraph::V(similarity_graph)$name
        edge_pairs <- igraph::as_data_frame(similarity_graph, what = "edges")[, c("from", "to")]

        res_param <- if (identical(objective, "modularity")) as.numeric(config$gamma) else as.numeric(config$resolution)
        algo_per_run <- if (identical(mode_selected, "aggressive") && identical(algo, "leiden")) "louvain" else algo

        old_max_size <- getOption("future.globals.maxSize")
        min_max_size <- 2 * 1024^3
        new_max_size <- old_max_size
        if (is.null(new_max_size) || (!is.infinite(new_max_size) && new_max_size < min_max_size)) {
            new_max_size <- min_max_size
        }
        options(future.globals.maxSize = new_max_size)
        on.exit(options(future.globals.maxSize = old_max_size), add = TRUE)

        aggressive_requested <- identical(mode_selected, "aggressive")
        if (aggressive_requested && !.networkit_available()) {
            stop("Aggressive mode requires the NetworKit backend; please configure reticulate/networkit before running clusterGenes_new.", call. = FALSE)
        }
        use_networkit_leiden <- aggressive_requested && identical(algo_per_run, "leiden")
        use_networkit_plm <- aggressive_requested && identical(algo_per_run, "louvain")

        if (use_networkit_leiden || use_networkit_plm) {
            if (isTRUE(verbose)) {
                message("[cluster]   Stage-1 backend: NetworKit (",
                    if (use_networkit_leiden) "Leiden" else "PLM/louvain", ")")
            }
            stage1_backend <- "networkit"
            seeds <- seq_len(config$n_restart)
            chunk_size <- config$aggr_batch_size
            if (is.null(chunk_size) || chunk_size <= 0) {
                chunk_size <- ceiling(config$n_restart / max(1L, as.integer(config$aggr_future_workers)))
            }
            split_chunks <- function(x, k) {
                if (k <= 0) return(list(x))
                split(x, ceiling(seq_along(x) / k))
            }
            chunks <- split_chunks(seeds, chunk_size)
            threads_per_worker <- max(1L, floor(n_threads / max(1L, as.integer(config$aggr_future_workers))))

            worker_fun <- function(ss) {
                if (use_networkit_leiden) {
                    .networkit_parallel_leiden_runs(
                        vertex_names = vertex_names,
                        edge_table = edge_pairs,
                        edge_weights = edge_weights,
                        gamma = as.numeric(res_param),
                        iterations = as.integer(config$nk_leiden_iterations),
                        randomize = isTRUE(config$nk_leiden_randomize),
                        threads = threads_per_worker,
                        seeds = as.integer(ss)
                    )
                } else {
                    .networkit_plm_runs(
                        vertex_names = vertex_names,
                        edge_table = edge_pairs,
                        edge_weights = edge_weights,
                        gamma = as.numeric(res_param),
                        refine = TRUE,
                        threads = threads_per_worker,
                        seeds = as.integer(ss)
                    )
                }
            }
            membership_blocks <- if (config$aggr_future_workers > 1L) {
                future.apply::future_lapply(chunks, worker_fun, future.seed = TRUE)
            } else {
                lapply(chunks, worker_fun)
            }
            restart_memberships <- do.call(cbind, membership_blocks)
        } else {
            if (isTRUE(verbose)) {
                message("[cluster]   Stage-1 backend: igraph (", algo_per_run, ")")
            }
            # Aggressive mode no longer falls back to igraph; the original igraph fallback is retained below for reference.
            restart_memberships <- future.apply::future_sapply(
                seq_len(config$n_restart),
                function(ii) {
                    g_local <- igraph::graph_from_data_frame(
                        cbind(edge_pairs, weight = edge_weights),
                        directed = FALSE,
                        vertices = vertex_names
                    )
                    res <- try(.run_graph_algorithm(g_local, algo_per_run, res_param, objective), silent = TRUE)
                    if (inherits(res, "try-error")) {
                        rep.int(1L, length(vertex_names))
                    } else {
                        mm <- try(igraph::membership(res), silent = TRUE)
                        if (inherits(mm, "try-error") || is.null(mm)) rep.int(1L, length(vertex_names)) else mm
                    }
                },
                future.seed = TRUE
            )
        }
        if (is.null(dim(restart_memberships))) restart_memberships <- matrix(restart_memberships, nrow = length(kept_genes))
        if (!is.matrix(restart_memberships) && !inherits(restart_memberships, "Matrix")) {
            restart_memberships <- matrix(restart_memberships, nrow = length(kept_genes))
        }
        if (is.null(nrow(restart_memberships)) || nrow(restart_memberships) != length(kept_genes)) {
            len <- length(restart_memberships)
            if (len %% length(kept_genes) == 0) {
                restart_memberships <- matrix(restart_memberships, nrow = length(kept_genes))
            } else {
                stop("Stage-1 membership shape mismatch: ", len, " elements vs ", length(kept_genes))
            }
        }
        rownames(restart_memberships) <- kept_genes

        if (!config$use_consensus) {
            stage1_membership_labels <- restart_memberships[, 1]
            stage1_consensus_matrix <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                dims = c(length(kept_genes), length(kept_genes)),
                dimnames = list(kept_genes, kept_genes))
        } else {
            stage1_consensus_matrix <- try(consensus_sparse(restart_memberships,
                thr = config$consensus_thr,
                n_threads = n_threads
            ), silent = TRUE)
            if (!inherits(stage1_consensus_matrix, "try-error")) {
                if (isTRUE(verbose)) message("[cluster]   Stage-1 consensus backend: C++ (consensus_sparse)")
            } else {
                if (isTRUE(verbose)) message("[cluster]   consensus_sparse failed; using R fallback.")
                stage1_consensus_matrix <- .compute_consensus_sparse_R(restart_memberships, thr = config$consensus_thr)
            }
            try({ if (length(stage1_consensus_matrix@i)) diag(stage1_consensus_matrix) <- 0 }, silent = TRUE)
            if (Matrix::nnzero(stage1_consensus_matrix) == 0) {
                stage1_membership_labels <- restart_memberships[, 1]
            } else {
                consensus_graph_view <- igraph::graph_from_adjacency_matrix(stage1_consensus_matrix,
                    mode = "undirected", weighted = TRUE, diag = FALSE)
                consensus_backend <- if (identical(stage1_backend, "networkit")) "networkit" else "igraph"
                if (isTRUE(verbose)) {
                    backend_label <- if (identical(consensus_backend, "networkit")) {
                        paste0("NetworKit (", algo, ")")
                    } else {
                        paste0("igraph (", algo, ")")
                    }
                    message("[cluster]   Stage-1 consensus graph backend: ", backend_label)
                }
                nk_opts <- list(
                    iterations = config$nk_leiden_iterations,
                    randomize = config$nk_leiden_randomize,
                    refine = TRUE,
                    seeds = 1L
                )
                stage1_membership_labels <- .run_single_partition(
                    consensus_graph_view,
                    backend = consensus_backend,
                    algo_name = algo,
                    res_param = res_param,
                    objective = objective,
                    threads = n_threads,
                    nk_opts = nk_opts
                )
            }
        }
    }

    # Stage-1 consensus graph weights for reporting
    stage1_edge_triplet <- summary(stage1_consensus_matrix)
    if (!all(c("i", "j", "x") %in% colnames(stage1_edge_triplet))) {
        colnames(stage1_edge_triplet) <- c("i", "j", "x")[seq_len(ncol(stage1_edge_triplet))]
    }
    stage1_edge_triplet <- stage1_edge_triplet[stage1_edge_triplet$i < stage1_edge_triplet$j & stage1_edge_triplet$x > 0, , drop = FALSE]
    stage1_cons_graph <- NULL
    if (nrow(stage1_edge_triplet)) {
        consensus_edge_df <- data.frame(
            from = kept_genes[stage1_edge_triplet$i],
            to = kept_genes[stage1_edge_triplet$j],
            weight = stage1_edge_triplet$x,
            stringsAsFactors = FALSE
        )
        stage1_cons_graph <- igraph::graph_from_data_frame(
            consensus_edge_df,
            directed = FALSE,
            vertices = kept_genes
        )
    } else {
        stage1_cons_graph <- igraph::make_empty_graph(n = length(kept_genes), directed = FALSE)
        stage1_cons_graph <- igraph::set_vertex_attr(stage1_cons_graph, "name", value = kept_genes)
    }

    list(
        kept_genes = kept_genes,
        genes_all = genes_all,
        similarity_full = threshold$similarity,
        similarity_sub = similarity_sub,
        A = threshold$similarity,          # backward compatibility
        A_sub = similarity_sub,            # backward compatibility
        g_stage1 = similarity_graph,
        mode_selected = mode_selected,
        stage1_membership_matrix = restart_memberships,
        stage1_consensus = stage1_consensus_matrix,
        stage1_membership = stage1_membership_labels,
        stage1_consensus_graph = stage1_cons_graph,
        stage1_backend = stage1_backend,
        stage1_algo = algo,
        stage1_algo_per_run = if (exists("algo_per_run")) algo_per_run else algo
    )
}

# ---------------------------------------------------------------------------
# Step 3: Stage-2 corrections + optional subclustering
# ---------------------------------------------------------------------------
#' Apply Stage-2 corrections and optional subclustering.
#'
#' @param stage1 List returned by `.stage1_consensus_workflow()`.
#' @param similarity_matrix Full similarity matrix.
#' @param config Configuration list.
#' @param FDR Optional FDR matrix aligned with similarity.
#' @param mh_object Optional MH weighting object.
#' @param aux_stats Auxiliary statistics list.
#' @param pearson_matrix Optional Pearson matrix for CI95 filter.
#' @param verbose Logical; emit progress messages.
#'
#' @return List with refined memberships, consensus graph, and diagnostics.
#' @noRd
.stage2_refine_workflow <- function(stage1,
                                   similarity_matrix,
                                   config,
                                   FDR = NULL,
                                   mh_object = NULL,
                                   aux_stats = NULL,
                                   pearson_matrix = NULL,
                                   verbose = TRUE) {
    if (isTRUE(verbose)) message("[cluster] Step3: Stage-2 correction & clustering…")

    kept_genes <- stage1$kept_genes
    genes_all <- stage1$genes_all
    stage1_membership_labels <- stage1$stage1_membership
    stage1_consensus_matrix <- stage1$stage1_consensus
    stage1_similarity_graph <- stage1$g_stage1
    stage1_mode_selected <- stage1$mode_selected
    if (is.null(stage1_mode_selected)) stage1_mode_selected <- "safe"
    stage1_membership_matrix <- stage1$stage1_membership_matrix
    stage1_backend <- stage1$stage1_backend %||% "igraph"
    stage1_algo_effective <- stage1$stage1_algo_per_run %||% stage1$stage1_algo %||% config$algo

    stage2_algo_final <- tryCatch({
        if (is.null(config$stage2_algo)) config$algo else match.arg(config$stage2_algo, c("leiden", "louvain", "hotspot-like"))
    }, error = function(e) config$algo)
    if (isTRUE(verbose)) message("[cluster]   Stage-2 base algorithm: ", stage2_algo_final)
    stage2_backend <- if (!is.null(config$stage2_algo)) {
        if (identical(stage2_algo_final, "hotspot-like")) "hotspot" else "igraph"
    } else {
        if (identical(stage2_algo_final, "hotspot-like")) {
            "hotspot"
        } else if (identical(stage1_backend, "networkit")) {
            "networkit"
        } else {
            "igraph"
        }
    }
    if (isTRUE(verbose)) message("[cluster]   Stage-2 target backend: ", stage2_backend)

    if (!config$use_mh_weight && !config$CI95_filter) {
        if (isTRUE(verbose)) message("[cluster]   No Stage-2 corrections requested; returning Stage-1 result.")
        return(list(
            membership = stage1_membership_labels,
            genes_all = genes_all,
            kept_genes = kept_genes,
            consensus_graph = stage1$stage1_consensus_graph,
            membership_matrix = stage1$stage1_membership_matrix,
            stage1_consensus = stage1_consensus_matrix,
            metrics_before = NULL,
            metrics_after = NULL,
            genes_before = NULL,
            genes_after = NULL,
            W = igraph::as_adjacency_matrix(stage1_similarity_graph, attr = "weight", sparse = TRUE)
        ))
    }

    L_post <- similarity_matrix[kept_genes, kept_genes, drop = FALSE]
    L_post[abs(L_post) < config$min_cutoff] <- 0
    if (config$use_significance && !is.null(FDR)) {
        L_post <- tryCatch(
            .align_and_filter_FDR(L_post, similarity_matrix, FDR, config$significance_max),
            error = function(e) {
                if (isTRUE(verbose)) message("[cluster]   significance subset disabled: ", conditionMessage(e))
                L_post
            }
        )
    }
    L_post[L_post < 0] <- 0
    L_post <- pmax(L_post, t(L_post))
    diag(L_post) <- 0
    L_post <- .filter_matrix_by_quantile(L_post, config$pct_min, "q100")
    if (!inherits(L_post, "sparseMatrix")) {
        L_post <- Matrix::drop0(Matrix::Matrix(L_post, sparse = TRUE))
    }

    # Optional CI95 filter
    if (config$CI95_filter) {
        if (is.null(aux_stats)) stop("CI95_filter requires aux_stats input.")
        curve_obj <- aux_stats[[config$curve_layer]]
        if (is.null(curve_obj)) stop("CI95 curve missing in stats layer.")
        curve_mat <- if (inherits(curve_obj, "big.matrix")) {
            bigmemory::as.matrix(curve_obj)
        } else if (is.matrix(curve_obj)) curve_obj else as.matrix(as.data.frame(curve_obj))
        xp <- as.numeric(curve_mat[, "Pear"])
        lo <- as.numeric(curve_mat[, "lo95"])
        hi <- as.numeric(curve_mat[, "hi95"])
        if (anyNA(c(xp, lo, hi))) stop("LR curve contains NA")
        f_lo <- approxfun(xp, lo, rule = 2)
        f_hi <- approxfun(xp, hi, rule = 2)
        if (is.null(pearson_matrix)) stop("CI95_filter requires pearson_matrix input.")
        rMat <- as.matrix(pearson_matrix[kept_genes, kept_genes, drop = FALSE])
        mask <- if (config$CI_rule == "remove_within") (L_post <= f_hi(rMat)) else (L_post < f_lo(rMat) | L_post > f_hi(rMat))
        L_post[mask] <- 0
    }

    # Optional MH re-weighting
    if (config$use_mh_weight) {
        if (is.null(mh_object)) stop("MH weighting requested but mh_object is NULL.")
        if (inherits(mh_object, "igraph")) {
            ed_sim <- igraph::as_data_frame(mh_object, what = "edges")
            wcol <- if ("MH" %in% names(ed_sim)) "MH" else if ("CMH" %in% names(ed_sim)) "CMH" else if ("weight" %in% names(ed_sim)) "weight" else stop("MH graph lacks weight attribute.")
            key_sim <- paste(pmin(ed_sim$from, ed_sim$to), pmax(ed_sim$from, ed_sim$to), sep = "|")
            sim_map <- setNames(ed_sim[[wcol]], key_sim)
            idx <- which(upper.tri(L_post) & L_post > 0, arr.ind = TRUE)
            if (nrow(idx)) {
                g1_names <- kept_genes[idx[, 1]]
                g2_names <- kept_genes[idx[, 2]]
                k <- paste(pmin(g1_names, g2_names), pmax(g1_names, g2_names), sep = "|")
                s <- sim_map[k]
                s[is.na(s)] <- stats::median(sim_map, na.rm = TRUE)
                L_post[cbind(idx[, 1], idx[, 2])] <- L_post[cbind(idx[, 1], idx[, 2])] * s
                L_post[cbind(idx[, 2], idx[, 1])] <- L_post[cbind(idx[, 2], idx[, 1])] * s
            }
        } else if (inherits(mh_object, "Matrix") || is.matrix(mh_object)) {
            M <- mh_object
            if (inherits(M, "Matrix") && !inherits(M, "CsparseMatrix")) M <- methods::as(M, "CsparseMatrix")
            if (!(identical(rownames(M), rownames(similarity_matrix)) && identical(colnames(M), colnames(similarity_matrix)))) {
                if (!is.null(rownames(M)) && !is.null(colnames(M)) &&
                    all(rownames(similarity_matrix) %in% rownames(M)) &&
                    all(colnames(similarity_matrix) %in% colnames(M))) {
                    M <- M[rownames(similarity_matrix), colnames(similarity_matrix), drop = FALSE]
                } else {
                    stop("MH matrix cannot be aligned to similarity matrix dimensions.")
                }
            }
            ridx <- match(kept_genes, rownames(M))
            M_sub <- M[ridx, ridx, drop = FALSE]
            idx <- which(upper.tri(L_post) & L_post > 0, arr.ind = TRUE)
            if (nrow(idx)) {
                s <- M_sub[cbind(idx[, 1], idx[, 2])]
                vec_all <- if (inherits(M_sub, "sparseMatrix")) M_sub@x else as.numeric(M_sub)
                vec_all <- vec_all[is.finite(vec_all)]
                s[is.na(s)] <- stats::median(vec_all, na.rm = TRUE)
                L_post[cbind(idx[, 1], idx[, 2])] <- L_post[cbind(idx[, 1], idx[, 2])] * s
                L_post[cbind(idx[, 2], idx[, 1])] <- L_post[cbind(idx[, 2], idx[, 1])] * s
            }
        } else {
            stop("Unsupported mh_object class: ", paste(class(mh_object), collapse = ","))
        }
    }

    L_post <- pmax(L_post, t(L_post))
    diag(L_post) <- 0

    if (isTRUE(config$keep_stage1_backbone)) {
        pos_vals <- as.numeric(L_post[L_post > 0])
        if (length(pos_vals) == 0) {
            w_floor <- 1e-6
        } else {
            qv <- stats::quantile(pos_vals, probs = min(max(config$backbone_floor_q, 0), 0.25), na.rm = TRUE)
            w_floor <- max(1e-8, as.numeric(qv))
        }
        sim_pos_full <- pmax(similarity_matrix, 0)
        cl_ids_stage1 <- sort(na.omit(unique(stage1_membership_labels)))
        for (cid in cl_ids_stage1) {
            genes_c <- kept_genes[stage1_membership_labels[kept_genes] == cid]
            if (length(genes_c) < 2) next
            M <- L_post[genes_c, genes_c, drop = FALSE]
            if (!any(M > 0, na.rm = TRUE)) {
                M <- sim_pos_full[genes_c, genes_c, drop = FALSE]
            }
            idx_back <- which(upper.tri(M) & M > 0, arr.ind = TRUE)
            if (nrow(idx_back) == 0) next
            mst_edge_table <- data.frame(from = genes_c[idx_back[, 1]], to = genes_c[idx_back[, 2]], w = M[idx_back], stringsAsFactors = FALSE)
            gtmp <- igraph::graph_from_data_frame(mst_edge_table, directed = FALSE, vertices = genes_c)
            if (igraph::ecount(gtmp) == 0) next
            mst_c <- igraph::mst(gtmp, weights = 1 / (igraph::E(gtmp)$w + 1e-9))
            if (igraph::ecount(mst_c) == 0) next
            ep <- igraph::as_data_frame(mst_c, what = "edges")
            for (k in seq_len(nrow(ep))) {
                i <- ep$from[k]; j <- ep$to[k]
                L_post[i, j] <- max(L_post[i, j], w_floor)
                L_post[j, i] <- L_post[i, j]
            }
        }
    }

    ed1 <- summary(L_post)
    ed1 <- ed1[ed1$i < ed1$j & ed1$x > 0, , drop = FALSE]
    if (!nrow(ed1)) stop("No edges remain after Stage-2 corrections.")
    edges_corr <- data.frame(
        from = kept_genes[ed1$i],
        to = kept_genes[ed1$j],
        L_corr = ed1$x,
        stringsAsFactors = FALSE
    )
    use_log1p_weight <- isTRUE(config$use_log1p_weight)
    edges_corr$w_raw <- if (use_log1p_weight) log1p(edges_corr$L_corr) else edges_corr$L_corr
    if (isTRUE(config$post_smooth)) {
        q <- stats::quantile(edges_corr$w_raw, probs = config$post_smooth_quant, na.rm = TRUE)
        w_clipped <- pmin(pmax(edges_corr$w_raw, q[1]), q[2])
        rng <- q[2] - q[1]
        if (!is.finite(rng) || rng <= 0) rng <- max(1e-8, max(w_clipped) - min(w_clipped))
        w01 <- (w_clipped - min(w_clipped)) / (rng + 1e-12)
        if (!is.null(config$post_smooth_power) && is.finite(config$post_smooth_power) && config$post_smooth_power != 1) {
            w01 <- w01 ^ config$post_smooth_power
        }
        edges_corr$weight <- w01
    } else {
        edges_corr$weight <- edges_corr$w_raw
    }

    corrected_similarity_graph <- igraph::graph_from_data_frame(
        edges_corr[, c("from", "to", "weight")],
        directed = FALSE,
        vertices = kept_genes
    )
    W <- as.matrix(igraph::as_adjacency_matrix(corrected_similarity_graph, attr = "weight", sparse = TRUE))
    rownames(W) <- colnames(W) <- kept_genes

    use_consensus <- isTRUE(config$use_consensus)
    consensus_thr <- config$consensus_thr
    n_restart <- max(1L, as.integer(config$n_restart))
    n_threads <- if (is.null(config$n_threads)) 1L else max(1L, as.integer(config$n_threads))
    objective <- config$objective
    res_param_base <- if (identical(objective, "modularity")) as.numeric(config$gamma) else as.numeric(config$resolution)
    if (!is.finite(res_param_base)) res_param_base <- 1
    mode_setting <- config$mode
    large_n_threshold <- if (is.null(config$large_n_threshold)) length(kept_genes) else as.integer(config$large_n_threshold)
    aggr_future_workers <- if (is.null(config$aggr_future_workers)) 1L else max(1L, as.integer(config$aggr_future_workers))
    aggr_batch_size <- config$aggr_batch_size
    nk_leiden_iterations <- if (is.null(config$nk_leiden_iterations)) 10L else as.integer(config$nk_leiden_iterations)
    nk_leiden_randomize <- if (is.null(config$nk_leiden_randomize)) TRUE else isTRUE(config$nk_leiden_randomize)
    min_cluster_keep <- if (is.null(config$min_cluster_size)) 2L else max(1L, as.integer(config$min_cluster_size))
    sub_min_size <- max(1L, as.integer(config$sub_min_size))
    sub_min_child_size <- max(1L, as.integer(config$sub_min_child_size))
    sub_resolution_factor <- if (is.null(config$sub_resolution_factor)) 1.3 else config$sub_resolution_factor
    sub_within_cons_max <- config$sub_within_cons_max
    sub_conductance_min <- config$sub_conductance_min
    sub_improve_within_cons_min <- config$sub_improve_within_cons_min
    sub_max_groups <- max(1L, as.integer(config$sub_max_groups))
    enable_qc_filter <- isTRUE(config$enable_qc_filter)
    qc_gene_intra_cons_min <- config$qc_gene_intra_cons_min
    qc_gene_best_out_cons_min <- config$qc_gene_best_out_cons_min
    qc_gene_intra_weight_q <- config$qc_gene_intra_weight_q
    keep_cross_stable <- isTRUE(config$keep_cross_stable)
    min_cutoff <- config$min_cutoff

    ng <- length(kept_genes)
    memb_final <- rep(NA_integer_, ng)
    names(memb_final) <- kept_genes
    cons <- Matrix::Matrix(0, ng, ng, sparse = TRUE, dimnames = list(kept_genes, kept_genes))

    cl_ids_stage1 <- sort(na.omit(unique(stage1_membership_labels)))
    next_global_id <- 1L
    for (cid in cl_ids_stage1) {
        idx <- which(stage1_membership_labels == cid)
        genes_c <- kept_genes[idx]
        if (!length(genes_c)) next
        g_sub <- igraph::induced_subgraph(corrected_similarity_graph, vids = genes_c)
        if (length(genes_c) <= 1 || igraph::ecount(g_sub) == 0) {
            memb_final[genes_c] <- next_global_id
            cons[genes_c, genes_c] <- 0
            next_global_id <- next_global_id + 1L
            next
        }

        if (identical(stage2_algo_final, "hotspot-like")) {
            S_sub <- igraph::as_adjacency_matrix(g_sub, attr = "weight", sparse = TRUE)
            S_sub <- methods::as(S_sub, "dgCMatrix")
            hotspot_k <- config$hotspot_k
            if (is.null(hotspot_k) || !is.finite(hotspot_k)) hotspot_k <- length(genes_c)
            hotspot_res <- .run_hotspot_clustering(S_sub, as.integer(max(1L, min(hotspot_k, nrow(S_sub)))), config$hotspot_min_module_size, use_log1p_weight)
            memb_sub <- hotspot_res$membership
            if (is.null(names(memb_sub))) names(memb_sub) <- genes_c
            cons_sub <- hotspot_res$consensus
            if (is.null(cons_sub)) {
                cons_sub <- Matrix::Matrix(0, length(memb_sub), length(memb_sub), sparse = TRUE,
                    dimnames = list(names(memb_sub), names(memb_sub)))
            }
        } else {
            mode_sub <- if (identical(mode_setting, "auto")) {
                if (length(genes_c) > large_n_threshold) "aggressive" else stage1_mode_selected
            } else stage1_mode_selected
            algo_per_run_sub <- stage2_algo_final
            ew <- igraph::E(g_sub)$weight
            el <- igraph::as_data_frame(g_sub, what = "edges")[, c("from", "to")]
            res_param_sub <- res_param_base

            backend_sub <- if (identical(stage2_backend, "networkit")) "networkit" else "igraph"
            if (backend_sub == "networkit" && !stage2_algo_final %in% c("leiden", "louvain")) {
                backend_sub <- "igraph"
            }
            if (backend_sub == "networkit" && !.nk_available()) {
                stop("Stage-2 NetworKit backend requested but NetworKit is unavailable.", call. = FALSE)
            }

            if (backend_sub == "networkit") {
                seeds <- seq_len(n_restart)
                chunk_size <- aggr_batch_size
                if (is.null(chunk_size) || chunk_size <= 0) {
                    chunk_size <- ceiling(n_restart / max(1L, aggr_future_workers))
                }
                split_chunks <- function(x, k) {
                    if (k <= 0) return(list(x))
                    split(x, ceiling(seq_along(x) / k))
                }
                chunks <- split_chunks(seeds, chunk_size)
                threads_per_worker <- max(1L, floor(n_threads / max(1L, aggr_future_workers)))
                vertex_names_sub <- igraph::V(g_sub)$name
                worker_fun <- function(ss) {
                    if (identical(stage2_algo_final, "leiden")) {
                        .nk_parallel_leiden_mruns(
                            vertex_names = vertex_names_sub,
                            edge_table = el,
                            edge_weights = ew,
                            gamma = as.numeric(res_param_sub),
                            iterations = nk_leiden_iterations,
                            randomize = nk_leiden_randomize,
                            threads = threads_per_worker,
                            seeds = as.integer(ss)
                        )
                    } else {
                        .nk_plm_mruns(
                            vertex_names = vertex_names_sub,
                            edge_table = el,
                            edge_weights = ew,
                            gamma = as.numeric(res_param_sub),
                            refine = TRUE,
                            threads = threads_per_worker,
                            seeds = as.integer(ss)
                        )
                    }
                }
                memb_blocks <- if (aggr_future_workers > 1L) {
                    future.apply::future_lapply(chunks, worker_fun, future.seed = TRUE)
                } else {
                    lapply(chunks, worker_fun)
                }
                memb_mat_sub <- do.call(cbind, memb_blocks)
            } else {
                memb_mat_sub <- future.apply::future_sapply(
                    seq_len(n_restart),
                    function(ii) {
                        g_local <- igraph::graph_from_data_frame(
                            cbind(el, weight = ew),
                            directed = FALSE,
                            vertices = igraph::V(g_sub)$name
                        )
                        res <- try(.run_graph_algorithm(g_local, algo_per_run_sub, res_param_sub, objective), silent = TRUE)
                        if (inherits(res, "try-error")) {
                            rep.int(1L, length(genes_c))
                        } else {
                            mm <- try(igraph::membership(res), silent = TRUE)
                            if (inherits(mm, "try-error") || is.null(mm)) rep.int(1L, length(genes_c)) else mm
                        }
                    },
                    future.seed = TRUE
                )
            }
            if (is.null(dim(memb_mat_sub))) {
                memb_mat_sub <- matrix(memb_mat_sub, ncol = 1)
            }
            rownames(memb_mat_sub) <- igraph::V(g_sub)$name

            if (!use_consensus) {
                if (isTRUE(verbose)) message("[cluster]   Stage-2 consensus backend: single-run (no aggregation)")
                memb_sub <- memb_mat_sub[, 1]
                labs <- as.integer(memb_sub)
                n_loc <- length(labs)
                split_idx <- split(seq_len(n_loc), labs)
                if (!length(split_idx)) {
                    cons_sub <- Matrix::Matrix(0, n_loc, n_loc, sparse = TRUE, dimnames = list(names(memb_sub), names(memb_sub)))
                } else {
                    ii <- integer(0); jj <- integer(0)
                    for (idv in split_idx) {
                        s <- length(idv)
                        if (s > 1) {
                            comb <- utils::combn(idv, 2)
                            ii <- c(ii, comb[1, ])
                            jj <- c(jj, comb[2, ])
                        }
                    }
                    cons_sub <- if (length(ii)) {
                        Matrix::sparseMatrix(i = ii, j = jj, x = rep(1, length(ii)), dims = c(n_loc, n_loc), symmetric = TRUE)
                    } else {
                        Matrix::Matrix(0, n_loc, n_loc, sparse = TRUE)
                    }
                    if (!is.null(names(memb_sub)) && length(names(memb_sub)) == n_loc) {
                        dimnames(cons_sub) <- list(names(memb_sub), names(memb_sub))
                    }
                }
            } else {
                cons_sub <- try(consensus_sparse(memb_mat_sub,
                    thr = consensus_thr,
                    n_threads = n_threads
                ), silent = TRUE)
                if (!inherits(cons_sub, "try-error")) {
                    if (isTRUE(verbose)) message("[cluster]   Stage-2 consensus backend: C++ (consensus_sparse)")
                } else {
                    if (isTRUE(verbose)) {
                        message("[cluster]   Stage-2 consensus backend: R fallback (.compute_consensus_sparse_R)")
                    }
                    cons_sub <- .compute_consensus_sparse_R(memb_mat_sub, thr = consensus_thr)
                }
                if (length(cons_sub@i)) Matrix::diag(cons_sub) <- 0
                if (Matrix::nnzero(cons_sub) == 0) {
                    memb_sub <- memb_mat_sub[, 1]
                } else {
                    g_cons_sub <- igraph::graph_from_adjacency_matrix(cons_sub, mode = "undirected", weighted = TRUE, diag = FALSE)
                    backend_consensus_sub <- backend_sub
                    if (!backend_consensus_sub %in% c("networkit", "igraph")) backend_consensus_sub <- "igraph"
                    if (isTRUE(verbose)) {
                        backend_label <- if (identical(backend_consensus_sub, "networkit")) {
                            paste0("NetworKit (", stage2_algo_final, ")")
                        } else {
                            paste0("igraph (", stage2_algo_final, ")")
                        }
                        message("[cluster]   Stage-2 consensus graph backend: ", backend_label)
                    }
                    nk_opts_cons <- list(
                        iterations = config$nk_leiden_iterations,
                        randomize = config$nk_leiden_randomize,
                        refine = TRUE,
                        seeds = 1L
                    )
                    memb_sub <- .run_single_partition(
                        g_cons_sub,
                        backend = backend_consensus_sub,
                        algo_name = stage2_algo_final,
                        res_param = res_param_base,
                        objective = objective,
                        threads = n_threads,
                        nk_opts = nk_opts_cons
                    )
                }
            }
        }

        sub_ids <- sort(unique(memb_sub))
        if (length(sub_ids) <= 1) {
            memb_final[genes_c] <- next_global_id
            next_global_id <- next_global_id + 1L
        } else {
            for (sid in sub_ids) {
                genes_sid <- names(memb_sub)[memb_sub == sid]
                memb_final[genes_sid] <- next_global_id
                next_global_id <- next_global_id + 1L
            }
        }

        block_names <- intersect(rownames(cons_sub), genes_c)
        if (length(block_names)) {
            cons[block_names, block_names] <- as.matrix(cons_sub[block_names, block_names])
        } else {
            cons[genes_c, genes_c] <- 0
        }
    }

    metrics_before <- .summarise_cluster_metrics(memb_final, cons, W)
    genes_before <- .summarise_gene_membership(memb_final, cons, W, kept_genes)

    if (enable_qc_filter) {
        q_w <- stats::quantile(genes_before$w_in, probs = qc_gene_intra_weight_q, na.rm = TRUE)
        bad <- which(!is.na(genes_before$p_in) &
            genes_before$p_in < qc_gene_intra_cons_min &
            genes_before$p_best_out > qc_gene_best_out_cons_min &
            genes_before$w_in <= q_w)
        if (length(bad)) {
            memb_final[match(genes_before$gene[bad], kept_genes)] <- NA
        }
    }

    if (config$enable_subcluster && !identical(stage2_algo_final, "hotspot-like")) {
        next_id <- if (length(na.omit(memb_final))) max(na.omit(memb_final)) + 1L else 1L
        for (row_idx in seq_len(nrow(metrics_before))) {
            cid <- metrics_before$cluster[row_idx]
            size <- metrics_before$size[row_idx]
            wcons <- metrics_before$within_cons[row_idx]
            condc <- metrics_before$conductance[row_idx]
            if (is.na(size) || size < sub_min_size || is.na(wcons) || is.na(condc)) next
            if (!(wcons < sub_within_cons_max && condc > sub_conductance_min)) next

            genes_c <- kept_genes[memb_final == cid]
            if (!length(genes_c)) next
            g_sub <- igraph::induced_subgraph(corrected_similarity_graph, vids = genes_c)
            if (igraph::ecount(g_sub) == 0) next
            ew <- igraph::E(g_sub)$weight
            el <- igraph::as_data_frame(g_sub, what = "edges")[, c("from", "to")]
            memb_mat_sub <- future.apply::future_sapply(
                seq_len(min(n_restart, 50)),
                function(ii) {
                    g_local <- igraph::graph_from_data_frame(cbind(el, weight = ew), directed = FALSE, vertices = igraph::V(g_sub)$name)
                    res <- .run_graph_algorithm(
                        g_local,
                        stage2_algo_final,
                        if (identical(objective, "modularity")) res_param_base * sub_resolution_factor else res_param_base * sub_resolution_factor,
                        objective
                    )
                    igraph::membership(res)
                },
                future.seed = TRUE
            )
            if (is.null(dim(memb_mat_sub))) memb_mat_sub <- matrix(memb_mat_sub, ncol = 1)
            rownames(memb_mat_sub) <- igraph::V(g_sub)$name
            cons_sub <- .compute_consensus_sparse_R(memb_mat_sub, thr = consensus_thr)
            if (Matrix::nnzero(cons_sub) == 0) next
            g_cons_sub <- igraph::graph_from_adjacency_matrix(cons_sub, mode = "undirected", weighted = TRUE, diag = FALSE)
            comm_sub <- .run_graph_algorithm(
                g_cons_sub,
                stage2_algo_final,
                if (identical(objective, "modularity")) res_param_base * sub_resolution_factor else res_param_base * sub_resolution_factor,
                objective
            )
            memb_sub <- igraph::membership(comm_sub)

            sub_ids <- sort(unique(memb_sub))
            if (length(sub_ids) <= 1 || length(sub_ids) > sub_max_groups) next
            sizes <- vapply(sub_ids, function(sid) sum(memb_sub == sid), integer(1))
            if (any(sizes < sub_min_child_size)) next

            wcons_sub <- vapply(sub_ids, function(sid) {
                idx_local <- which(memb_sub == sid)
                if (length(idx_local) > 1) mean(cons_sub[idx_local, idx_local][upper.tri(cons_sub[idx_local, idx_local])], na.rm = TRUE) else 0
            }, numeric(1))
            wcons_old <- wcons
            wcons_new <- sum(wcons_sub * sizes) / sum(sizes)
            if (!(is.finite(wcons_new) && wcons_new >= (wcons_old + sub_improve_within_cons_min))) next

            for (sid in sub_ids) {
                genes_sid <- names(memb_sub)[memb_sub == sid]
                memb_final[genes_sid] <- next_id
                next_id <- next_id + 1L
            }
            block_names <- intersect(rownames(cons_sub), genes_c)
            if (length(block_names)) {
                cons[block_names, block_names] <- as.matrix(cons_sub[block_names, block_names])
            }
        }
    }

    di <- Matrix::rowSums((L_post > 0)[kept_genes, kept_genes, drop = FALSE] &
        outer(memb_final, memb_final, "=="), na.rm = TRUE)
    memb_final[di == 0] <- NA
    tbl <- table(memb_final, useNA = "no")
    small <- names(tbl)[tbl < min_cluster_keep]
    if (length(small)) {
        memb_final[memb_final %in% small] <- NA
    }
    tbl <- table(memb_final, useNA = "no")
    if (length(tbl)) {
        unique_ids <- sort(as.integer(names(tbl)))
        id_map <- setNames(seq_along(unique_ids), unique_ids)
        idx_keep <- which(!is.na(memb_final))
        if (length(idx_keep)) {
            memb_final[idx_keep] <- as.integer(id_map[as.character(memb_final[idx_keep])])
        }
        tbl <- table(memb_final, useNA = "no")
        stable_cl <- as.integer(names(tbl)[tbl >= min_cluster_keep])
    } else {
        stable_cl <- integer(0)
    }
    names(memb_final) <- kept_genes

    metrics_after <- .summarise_cluster_metrics(memb_final, cons, W)
    genes_after <- .summarise_gene_membership(memb_final, cons, W, kept_genes)

    cons_ts <- methods::as(cons, "TsparseMatrix")
    if (length(cons_ts@x)) {
        sm2 <- data.frame(
            i = cons_ts@i + 1L,
            j = cons_ts@j + 1L,
            x = cons_ts@x,
            stringsAsFactors = FALSE
        )
    } else {
        sm2 <- data.frame(i = integer(0), j = integer(0), x = numeric(0), stringsAsFactors = FALSE)
    }
    sm2 <- sm2[sm2$i < sm2$j & sm2$x > 0, , drop = FALSE]

    key_corr <- paste(pmin(edges_corr$from, edges_corr$to), pmax(edges_corr$from, edges_corr$to), sep = "|")
    w_map <- setNames(edges_corr$weight, key_corr)
    intra_keys <- paste(
        pmin(kept_genes[sm2$i], kept_genes[sm2$j]),
        pmax(kept_genes[sm2$i], kept_genes[sm2$j]),
        sep = "|"
    )
    w_vec <- as.numeric(w_map[intra_keys])
    keep_mask <- !is.na(w_vec) & w_vec > 0
    intra_df <- data.frame(
        from = kept_genes[sm2$i[keep_mask]],
        to = kept_genes[sm2$j[keep_mask]],
        weight = w_vec[keep_mask],
        sign = ifelse(similarity_matrix[cbind(
            match(kept_genes[sm2$i[keep_mask]], rownames(similarity_matrix)),
            match(kept_genes[sm2$j[keep_mask]], colnames(similarity_matrix))
        )] < 0, "neg", "pos"),
        stringsAsFactors = FALSE
    )

    cross_df <- NULL
    if (keep_cross_stable && length(stable_cl) > 1) {
        idx_stable <- memb_final %in% stable_cl
        genes_stable <- kept_genes[idx_stable]
        if (length(genes_stable) > 1) {
            L_sub_raw <- L_post[genes_stable, genes_stable, drop = FALSE]
            ij <- which(L_sub_raw >= min_cutoff & upper.tri(L_sub_raw), arr.ind = TRUE)
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
                    wL <- if (use_log1p_weight) log1p(Lval) else Lval
                    mask <- wL > 0
                    if (any(mask)) {
                        cross_df <- data.frame(
                            from = g1[mask],
                            to = g2[mask],
                            weight = wL[mask],
                            sign = ifelse(Lval[mask] < 0, "neg", "pos")
                        )
                    }
                }
            }
        }
    }

    edges_df_final <- if (!is.null(cross_df)) unique(rbind(intra_df, cross_df)) else intra_df
    if (!nrow(edges_df_final)) {
        cons_graph_final <- igraph::make_empty_graph(n = length(kept_genes), directed = FALSE)
        cons_graph_final <- igraph::set_vertex_attr(cons_graph_final, "name", value = kept_genes)
    } else {
        cons_graph_final <- igraph::graph_from_data_frame(edges_df_final,
            directed = FALSE, vertices = kept_genes
        )
    }

    list(
        membership = memb_final,
        kept_genes = kept_genes,
        genes_all = genes_all,
        consensus_graph = cons_graph_final,
        membership_matrix = stage1_membership_matrix,
        stage1_consensus = stage1_consensus_matrix,
        metrics_before = metrics_before,
        metrics_after = metrics_after,
        genes_before = genes_before,
        genes_after = genes_after,
        corrected_graph = corrected_similarity_graph,
        weights_matrix = W,
        stage2_backend = stage2_backend
    )
}

# ---------------------------------------------------------------------------
# Wrapper: scope_obj compatible clustering (public entry)
# ---------------------------------------------------------------------------
#' Cluster genes using the multi-stage consensus pipeline.
#'
#' @param scope_obj CosMx scope object containing grid statistics.
#' @param grid_name Character grid identifier.
#' @param stats_layer Statistics layer name (default `"LeeStats_Xz"`).
#' @param lee_stats_layer Optional legacy LeeStats layer argument.
#' @param similarity_slot Slot containing similarity matrix (default `"L"`).
#' @param significance_slot Slot containing significance matrix (default `"FDR"`).
#' @param min_cutoff Minimum similarity cutoff (aliases `L_min`).
#' @param use_significance Logical; apply significance filtering.
#' @param significance_max Maximum allowed FDR (aliases `FDR_max`).
#' @param pct_min Quantile filter applied to similarity weights.
#' @param drop_isolated Logical; drop isolated nodes post-thresholding.
#' @param algo Stage-1 clustering algorithm (`"leiden"`, `"louvain"`, `"hotspot-like"`).
#' @param stage2_algo Optional Stage-2 algorithm override.
#' @param resolution CPM resolution parameter (for `"CPM"` objective).
#' @param gamma Modularity resolution parameter (for `"modularity"` objective).
#' @param objective Clustering objective (`"CPM"` or `"modularity"`).
#' @param cluster_name Column name written into `scope_obj@meta.data`.
#' @param use_log1p_weight Logical; log1p transform edge weights.
#' @param use_consensus Logical; enable consensus aggregation.
#' @param graph_slot_name Stats-layer slot to store consensus graph.
#' @param n_restart Number of Stage-1 restarts.
#' @param consensus_thr Consensus threshold applied to co-occurrence matrix.
#' @param keep_cross_stable Logical; retain cross-cluster stable edges.
#' @param CI95_filter Logical; enable CI95-based edge filtering.
#' @param curve_layer Name of curve layer for CI95 filtering.
#' @param CI_rule Rule applied to CI95 mask (`"remove_within"` or `"remove_outside"`).
#' @param use_mh_weight Logical; enable Morisita-Horn weighting.
#' @param mh_slot Slot containing MH weights (aliases `cmh_slot`).
#' @param K Desired number of modules when `algo = "hotspot-like"`.
#' @param min_module_size Minimum module size for hotspot-like clustering.
#' @param post_smooth Logical; smooth corrected weights.
#' @param post_smooth_quant Quantiles used for smoothing clamp.
#' @param post_smooth_power Exponent applied to re-scaled weights.
#' @param enable_subcluster Logical; run conservative Stage-2 subclustering.
#' @param sub_min_size Minimum cluster size eligible for subclustering.
#' @param sub_min_child_size Minimum size for subclusters.
#' @param sub_resolution_factor Resolution multiplier for Stage-2 runs.
#' @param sub_within_cons_max Maximum within-consensus to trigger splitting.
#' @param sub_conductance_min Minimum conductance to trigger splitting.
#' @param sub_improve_within_cons_min Minimum improvement required to accept split.
#' @param sub_max_groups Maximum resulting subclusters per parent.
#' @param enable_qc_filter Logical; enable QC-based gene pruning.
#' @param qc_gene_intra_cons_min Minimum within-consensus for genes.
#' @param qc_gene_best_out_cons_min Maximum best-out consensus for genes.
#' @param qc_gene_intra_weight_q Weight quantile for QC heuristic.
#' @param min_cluster_size Minimum cluster size retained after Stage-2 cleanup.
#' @param keep_stage1_backbone Logical; preserve Stage-1 MST backbone.
#' @param backbone_floor_q Quantile floor applied when restoring backbone edges.
#' @param return_report Logical; return diagnostics alongside modified object.
#' @param verbose Logical; emit progress messages.
#' @param ncores Parallel threads for consensus routines.
#' @param mode Clustering mode (`"auto"`, `"safe"`, `"aggressive"`).
#' @param large_n_threshold Node threshold for switching to aggressive mode.
#' @param nk_condaenv/nk_conda_bin/nk_python Optional NetworKit configuration.
#' @param nk_leiden_iterations Iterations for NetworKit Leiden backend.
#' @param nk_leiden_randomize Logical; randomize NetworKit Leiden.
#' @param aggr_future_workers Future workers for batched restarts.
#' @param aggr_batch_size Batch size for restart aggregation.
#' @param profile_timing Logical; collect timing metrics.
#'
#' @return Modified `scope_obj`, or list with `scope_obj` and diagnostics when
#'   `return_report = TRUE`.
#' @export
clusterGenes <- function(
    scope_obj,
    grid_name = NULL,
    stats_layer = "LeeStats_Xz",
    lee_stats_layer = NULL,
    similarity_slot = "L",
    significance_slot = "FDR",
    min_cutoff = 0,
    L_min = NULL,
    use_significance = TRUE,
    use_FDR = NULL,
    significance_max = 0.05,
    FDR_max = NULL,
    pct_min = "q0",
    drop_isolated = TRUE,
    algo = c("leiden", "louvain", "hotspot-like"),
    stage2_algo = NULL,
    resolution = 1,
    gamma = 1,
    objective = c("CPM", "modularity"),
    cluster_name = NULL,
    use_log1p_weight = TRUE,
    use_consensus = TRUE,
    graph_slot_name = "g_consensus",
    n_restart = 100,
    consensus_thr = 0.95,
    K = NULL,
    min_module_size = 5,
    CI95_filter = FALSE,
    curve_layer = "LR_curve",
    CI_rule = c("remove_within", "remove_outside"),
    use_mh_weight = FALSE,
    use_cmh_weight = NULL,
    mh_slot = NULL,
    cmh_slot = "g_morisita",
    post_smooth = TRUE,
    post_smooth_quant = c(0.05, 0.95),
    post_smooth_power = 0.5,
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
    keep_cross_stable = TRUE,
    min_cluster_size = 2,
    keep_stage1_backbone = TRUE,
    backbone_floor_q = 0.02,
    return_report = FALSE,
    verbose = TRUE,
    ncores = NULL,
    mode = c("auto", "safe", "aggressive"),
    large_n_threshold = 1000L,
    nk_condaenv = NULL,
    nk_conda_bin = NULL,
    nk_python = NULL,
    nk_leiden_iterations = 10L,
    nk_leiden_randomize = TRUE,
    aggr_future_workers = 2L,
    aggr_batch_size = NULL,
    profile_timing = FALSE) {

    if (length(algo) == 1L && identical(tolower(algo), "hotspot")) {
        algo <- "hotspot-like"
    }
    algo <- match.arg(algo)
    objective <- match.arg(objective)
    CI_rule <- match.arg(CI_rule)
    mode <- match.arg(mode)

    if (!is.null(lee_stats_layer)) stats_layer <- lee_stats_layer
    if (!is.null(L_min)) min_cutoff <- L_min
    if (!is.null(use_FDR)) use_significance <- use_FDR
    if (!is.null(FDR_max)) significance_max <- FDR_max
    if (!is.null(use_cmh_weight)) use_mh_weight <- use_cmh_weight
    if (!is.null(cmh_slot)) mh_slot <- cmh_slot
    if (!is.null(stage2_algo) && length(stage2_algo) == 1L && identical(tolower(stage2_algo), "hotspot")) {
        stage2_algo <- "hotspot-like"
    }

    cfg <- .build_pipeline_config(
        min_cutoff = min_cutoff,
        use_significance = use_significance,
        significance_max = significance_max,
        pct_min = pct_min,
        drop_isolated = drop_isolated,
        algo = algo,
        stage2_algo = stage2_algo,
        resolution = resolution,
        gamma = gamma,
        objective = objective,
        use_log1p_weight = use_log1p_weight,
        use_consensus = use_consensus,
        consensus_thr = consensus_thr,
        n_restart = n_restart,
        n_threads = ncores,
        mode = mode,
        large_n_threshold = large_n_threshold,
        nk_leiden_iterations = nk_leiden_iterations,
        nk_leiden_randomize = nk_leiden_randomize,
        aggr_future_workers = aggr_future_workers,
        aggr_batch_size = aggr_batch_size,
        enable_subcluster = enable_subcluster,
        sub_min_size = sub_min_size,
        sub_min_child_size = sub_min_child_size,
        sub_resolution_factor = sub_resolution_factor,
        sub_within_cons_max = sub_within_cons_max,
        sub_conductance_min = sub_conductance_min,
        sub_improve_within_cons_min = sub_improve_within_cons_min,
        sub_max_groups = sub_max_groups,
        enable_qc_filter = enable_qc_filter,
        qc_gene_intra_cons_min = qc_gene_intra_cons_min,
        qc_gene_best_out_cons_min = qc_gene_best_out_cons_min,
        qc_gene_intra_weight_q = qc_gene_intra_weight_q,
        keep_cross_stable = keep_cross_stable,
        min_cluster_size = min_cluster_size,
        CI95_filter = CI95_filter,
        CI_rule = CI_rule,
        curve_layer = curve_layer,
        similarity_slot = similarity_slot,
        significance_slot = significance_slot,
        use_mh_weight = use_mh_weight,
        mh_slot = mh_slot,
        post_smooth = post_smooth,
        post_smooth_quant = post_smooth_quant,
        post_smooth_power = post_smooth_power,
        keep_stage1_backbone = keep_stage1_backbone,
        backbone_floor_q = backbone_floor_q,
        hotspot_k = if (is.null(K)) NULL else K,
        hotspot_min_module_size = min_module_size
    )

    timing <- list()
    .tic <- function() proc.time()
    .toc <- function(t0) as.numeric((proc.time() - t0)["elapsed"])
    t_all0 <- .tic()

    inputs <- .extract_scope_layers(scope_obj, grid_name,
        stats_layer = stats_layer,
        similarity_slot = cfg$similarity_slot,
        significance_slot = cfg$significance_slot,
        use_significance = cfg$use_significance,
        verbose = verbose)

    th <- .threshold_similarity_graph(inputs$similarity, inputs$significance, cfg, verbose = verbose)
    stage1 <- .stage1_consensus_workflow(inputs$similarity, th, cfg, verbose = verbose,
        networkit_config = list(conda_env = nk_condaenv, conda_bin = nk_conda_bin, python_path = nk_python))
    timing$stage1 <- .toc(t_all0)

    mh_object <- NULL
    if (cfg$use_mh_weight) {
        lstats <- inputs$stats_object
        mh_object <- lstats[[cfg$mh_slot]]
        if (is.null(mh_object)) stop("MH slot missing: ", cfg$mh_slot)
    }

    pearson_mat <- NULL
    if (cfg$CI95_filter) {
        pearson_mat <- try(.getPearsonMatrix(scope_obj, grid_name = inputs$grid_name, level = "cell"), silent = TRUE)
        if (inherits(pearson_mat, "try-error") || is.null(pearson_mat)) {
            pearson_mat <- try(.getPearsonMatrix(scope_obj, grid_name = inputs$grid_name, level = "grid"), silent = TRUE)
        }
        if (inherits(pearson_mat, "try-error") || is.null(pearson_mat)) {
            stop("CI95_filter requested but Pearson matrix not found at cell/grid level.")
        }
    }

    stage2 <- .stage2_refine_workflow(stage1,
        similarity_matrix = inputs$similarity,
        config = cfg,
        FDR = inputs$significance,
        mh_object = mh_object,
        aux_stats = inputs$stats_object,
        pearson_matrix = pearson_mat,
        verbose = verbose)
    timing$total <- .toc(t_all0)

    # Update scope_obj meta.data
    genes_all <- stage2$genes_all
    kept_genes <- stage2$kept_genes
    memb_final <- stage2$membership

    if (is.null(cluster_name) || !nzchar(cluster_name)) {
        cluster_name <- if (identical(algo, "hotspot-like")) sprintf("modHS%.2f", min_cutoff) else sprintf("modL%.2f", min_cutoff)
    }

    md <- scope_obj@meta.data
    if (is.null(md) || !is.data.frame(md)) md <- try(as.data.frame(md), silent = TRUE)
    if (inherits(md, "try-error") || is.null(md) || !is.data.frame(md) || is.null(rownames(md))) {
        md <- data.frame(row.names = genes_all)
    } else {
        miss_rows <- setdiff(genes_all, rownames(md))
        if (length(miss_rows)) md <- rbind(md, data.frame(row.names = miss_rows))
    }
    md[genes_all, cluster_name] <- NA_integer_
    md[kept_genes, cluster_name] <- as.integer(memb_final[kept_genes])
    scope_obj@meta.data <- md

    if (isTRUE(verbose)) {
        nz <- sum(!is.na(scope_obj@meta.data[genes_all, cluster_name]))
        message("[clusterGenes] Written clustering column '", cluster_name, "': ",
            nz, "/", length(genes_all), " genes (",
            sprintf("%.1f%%", 100 * nz / length(genes_all)), ")")
    }

    # Update stats layer with final consensus graph
    gname <- inputs$grid_name
    if (is.null(scope_obj@stats[[gname]][[stats_layer]])) {
        scope_obj@stats[[gname]][[stats_layer]] <- list()
    }
    scope_obj@stats[[gname]][[stats_layer]][[graph_slot_name]] <- stage2$consensus_graph

    report <- NULL
    if (isTRUE(return_report)) {
        report <- list()
        report$timing <- timing
        report$stage1 <- list(
            mode = stage1$mode_selected,
            kept_genes = length(stage1$kept_genes)
        )
        report$metrics_before <- stage2$metrics_before
        report$metrics_after <- stage2$metrics_after
        report$genes_before <- stage2$genes_before
        report$genes_after <- stage2$genes_after
    }

    if (isTRUE(return_report)) {
        return(list(scope_obj = scope_obj, report = report))
    }
    scope_obj
}
