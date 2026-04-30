#' Cluster genes on a consensus graph.
#' @description
#' Runs the stage1/stage2 gene clustering pipeline on a `scope_object`, using a
#' similarity/statistics layer (typically produced by `computeL()` / `computeMH()`)
#' and writing cluster assignments back into `scope_obj@meta.data`.
#' @param scope_obj A \code{scope_object} containing graph/stat layers.
#' @param grid_name Character; name of the grid layer to operate on. If \code{NULL}
#'   and only one grid layer exists, it is auto-selected. If multiple layers exist,
#'   the name must be specified.
#' @param stats_layer Statistics layer name (e.g. `"LeeStats_Xz"`).
#' @param lee_stats_layer Optional alias of `stats_layer` (preferred for clarity).
#' @param similarity_slot Name of similarity matrix slot within the stats layer (default `"L"`).
#' @param significance_slot Name of significance/FDR matrix slot within the stats layer (default `"FDR"`).
#' @param min_cutoff Minimum similarity cutoff (alias: `L_min`).
#' @param L_min Deprecated alias of `min_cutoff`.
#' @param use_significance Whether to apply significance/FDR gating.
#' @param use_FDR Deprecated alias of `use_significance`.
#' @param significance_max Maximum allowed significance/FDR (alias: `FDR_max`).
#' @param FDR_max Deprecated alias of `significance_max`.
#' @param pct_min Gene prevalence filter used in stage1 (e.g. `"q0"`).
#' @param drop_isolated Whether to drop isolated nodes before clustering.
#' @param algo Stage1 clustering algorithm (`leiden`, `louvain`, `hotspot-like`).
#' @param stage2_algo Optional override of stage2 algorithm (defaults to `algo`).
#' @param resolution Resolution parameter passed to Leiden/Louvain. Must be a
#'   single finite numeric value greater than 0. **Default is 1**, but note that
#'   high resolution values (>= 0.5) combined with strict consensus thresholds
#'   (consensus_thr >= 0.9) may produce empty or near-empty clustering results
#'   in some datasets. For typical spatial gene clustering workflows, values
#'   in the range 0.05-0.1 often yield more robust results. A warning is emitted
#'   when the combination of resolution and consensus_thr appears likely to
#'   produce sparse outputs.
#' @param gamma Gamma parameter used by some clustering backends.
#' @param objective Optimization objective (backend-specific; e.g. `"CPM"` or `"modularity"`).
#' @param cluster_name Column name to store final memberships under `scope_obj@meta.data`.
#' @param use_log1p_weight Whether to log1p-transform edge weights before clustering.
#' @param use_consensus Whether to build/use a consensus graph across restarts.
#' @param graph_slot_name Graph slot used to drive clustering (e.g., consensus graph).
#' @param n_restart Number of random restarts for consensus building.
#' @param consensus_thr Consensus threshold applied when building the final graph.
#' @param K Hotspot-like neighbourhood parameter (backend-specific; optional).
#' @param min_module_size Minimum module size for hotspot-like backend.
#' @param CI95_filter Deprecated flag; use `CI_rule` instead.
#' @param CI_rule Optional confidence-interval filtering rule (`remove_within`, `remove_outside`).
#' @param curve_layer Name of curve layer used by CI filtering (required when `CI_rule` is set).
#' @param use_mh_weight Whether to use MH weighting (auto-detected from `mh_slot` when NULL).
#' @param use_cmh_weight Deprecated alias of `use_mh_weight`.
#' @param mh_slot Name of MH weight slot.
#' @param cmh_slot Deprecated alias of `mh_slot`.
#' @param post_smooth Whether to apply post-smoothing of weights/scores.
#' @param post_smooth_quant Quantile range for post-smoothing.
#' @param post_smooth_power Power used in post-smoothing.
#' @param enable_subcluster Whether to enable the subclustering pass.
#' @param sub_min_size Minimum parent-cluster size eligible for subclustering.
#' @param sub_min_child_size Minimum subcluster size.
#' @param sub_resolution_factor Resolution multiplier for subclustering.
#' @param sub_within_cons_max Maximum within-consensus threshold for subclustering.
#' @param sub_conductance_min Minimum conductance threshold for subclustering.
#' @param sub_improve_within_cons_min Minimum improvement threshold for within-consensus to accept a subcluster.
#' @param sub_max_groups Maximum number of subclusters produced per parent.
#' @param enable_qc_filter Whether to apply QC-based gene/module filtering.
#' @param qc_gene_intra_cons_min QC: minimum within-module consensus threshold.
#' @param qc_gene_best_out_cons_min QC: minimum "best outside" consensus threshold.
#' @param qc_gene_intra_weight_q QC: quantile for intra-module weight thresholding.
#' @param keep_cross_stable Whether to retain cross-stable edges/modules.
#' @param min_cluster_size Minimum allowed final cluster size. Must be a single
#'   finite numeric value greater than or equal to 1.
#' @param keep_stage1_backbone Whether to keep stage1 backbone edges.
#' @param backbone_floor_q Quantile floor for backbone retention.
#' @param return_report When `TRUE`, returns both updated scope object and a report list.
#' @param verbose Logical; whether to emit compact public progress messages.
#'   Default is \code{TRUE}. Set to \code{FALSE} for silent operation, or use
#'   \code{options(geneSCOPE.verbose = FALSE)} for global control.
#' @param ncores Integer; number of threads to use for parallel computation.
#' @param mode Runtime mode (`auto`, `safe`, `safe_sequential`, `aggressive`, `fast`).
#'   `safe_sequential` is the supported large-graph single-worker path for
#'   avoiding oversized future exports in Visium/CosMx runs.
#' @param backend Graph backend contract. `auto` keeps the existing runtime
#'   heuristics, `igraph` forces the safe igraph path, and `networkit` forces
#'   the aggressive NetworKit path (erroring if NetworKit is unavailable).
#' @param q_guard_method Q-distribution guardrail policy. `auto` (default)
#'   checks the similarity spread before clustering and switches to the
#'   compatibility `k_perms` fallback contract when the distribution is too
#'   narrow, `k_perms` forces the same fallback decision when triggered, and
#'   `none` disables the guard.
#' @param q_span_threshold Minimum allowed q95-q05 similarity span before the
#'   q-distribution guard recommends the `k_perms` fallback contract.
#' @param large_n_threshold Threshold controlling large-N heuristics.
#' @param nk_condaenv Optional conda environment name for Networkit/Python backends.
#' @param nk_conda_bin Optional path to conda binary.
#' @param nk_python Optional explicit python path for Networkit backends.
#' @param nk_leiden_iterations Leiden iterations (Networkit backend).
#' @param nk_leiden_randomize Whether to randomize Leiden (Networkit backend).
#' @param aggr_future_workers Number of future workers used by aggressive modes.
#' @param aggr_batch_size Batch size used in aggressive modes.
#' @param future_globals_min_bytes Minimum globals size threshold for futures.
#' @param profile_timing Whether to collect timing measurements for the returned report.
#' @return
#' If `return_report = FALSE` (default), returns the updated `scope_object`.
#' If `return_report = TRUE`, returns `list(scope_obj = <scope_object>, report = <list>)`.
#' @examples
#' \dontrun{
#' scope_obj <- createSCOPE(data_dir = "/path/to/output_dir",　grid_length = 30)
#' scope_obj <- normalizeMoleculesInGrid(scope_obj, grid_name = "grid30")
#' scope_obj <- computeL(scope_obj, grid_name = "grid30")
#' scope_obj <- computeMH(scope_obj, grid_name = "grid30")
#' scope_obj <- clusterGenes(scope_obj, grid_name = "grid30", pct_min = "q50", n_start = 100, consensus_thr = 0.9, resolution = 0.5)
#' }
#' @seealso `computeL()`, `computeMH()`, `plotNetwork()`, `plotDendroNetwork()`
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
    CI95_filter = NULL,
    CI_rule = NULL,
    curve_layer = NULL,
    use_mh_weight = NULL,
    use_cmh_weight = NULL,
    mh_slot = NULL,
    cmh_slot = NULL,
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
    mode = c("auto", "safe", "safe_sequential", "aggressive", "fast"),
    backend = c("auto", "igraph", "networkit"),
    q_guard_method = c("auto", "k_perms", "none"),
    q_span_threshold = 0.1,
    large_n_threshold = 1000L,
    nk_condaenv = NULL,
    nk_conda_bin = NULL,
    nk_python = NULL,
    nk_leiden_iterations = 10L,
    nk_leiden_randomize = TRUE,
    aggr_future_workers = 2L,
    aggr_batch_size = NULL,
    future_globals_min_bytes = 2 * 1024^3,
    profile_timing = FALSE) {
    verbose <- .resolve_runtime_verbose(verbose)
    
    # Emit warning for high resolution + strict consensus combination
    if (resolution >= 0.5 && consensus_thr >= 0.9) {
        warning(
            sprintf(
                "High resolution (%.2f) combined with strict consensus_thr (%.2f) may produce sparse clustering results. Consider using resolution in range 0.05-0.1 for typical spatial gene clustering workflows.",
                resolution, consensus_thr
            ),
            call. = FALSE
        )
    }
    
    .with_log_session("clusterGenes", verbose = verbose, {
        .log_start("clusterGenes", verbose = verbose)
        
        requested_return_report <- isTRUE(return_report)
        backend <- match.arg(backend)
        q_guard_method <- match.arg(q_guard_method)

        .log_key_values("clusterGenes", list(
            grid = if (!is.null(grid_name)) as.character(grid_name)[1] else NA_character_,
            stats_layer = if (!is.null(lee_stats_layer) && nzchar(lee_stats_layer)) lee_stats_layer else stats_layer,
            algo = as.character(algo)[1],
            resolution = resolution,
            pct_min = pct_min,
            backend = backend,
            q_guard_method = q_guard_method
        ))

        result <- .cluster_genes(
            scope_obj = scope_obj,
            grid_name = grid_name,
            stats_layer = stats_layer,
            lee_stats_layer = lee_stats_layer,
            similarity_slot = similarity_slot,
            significance_slot = significance_slot,
            min_cutoff = min_cutoff,
            L_min = L_min,
            use_significance = use_significance,
            use_FDR = use_FDR,
            significance_max = significance_max,
            FDR_max = FDR_max,
            pct_min = pct_min,
            drop_isolated = drop_isolated,
            algo = algo,
            stage2_algo = stage2_algo,
            resolution = resolution,
            gamma = gamma,
            objective = objective,
            cluster_name = cluster_name,
            use_log1p_weight = use_log1p_weight,
            use_consensus = use_consensus,
            graph_slot_name = graph_slot_name,
            n_restart = n_restart,
            consensus_thr = consensus_thr,
            K = K,
            min_module_size = min_module_size,
            CI95_filter = CI95_filter,
            CI_rule = CI_rule,
            curve_layer = curve_layer,
            use_mh_weight = use_mh_weight,
            use_cmh_weight = use_cmh_weight,
            mh_slot = mh_slot,
            cmh_slot = cmh_slot,
            post_smooth = post_smooth,
            post_smooth_quant = post_smooth_quant,
            post_smooth_power = post_smooth_power,
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
            keep_stage1_backbone = keep_stage1_backbone,
            backbone_floor_q = backbone_floor_q,
            return_report = TRUE,
            verbose = verbose,
            ncores = ncores,
            mode = mode,
            backend = backend,
            q_guard_method = q_guard_method,
            q_span_threshold = q_span_threshold,
            large_n_threshold = large_n_threshold,
            nk_condaenv = nk_condaenv,
            nk_conda_bin = nk_conda_bin,
            nk_python = nk_python,
            nk_leiden_iterations = nk_leiden_iterations,
            nk_leiden_randomize = nk_leiden_randomize,
            aggr_future_workers = aggr_future_workers,
            aggr_batch_size = aggr_batch_size,
            future_globals_min_bytes = future_globals_min_bytes,
            profile_timing = profile_timing
        )

        scope_out <- result$scope_obj
        report <- result$report
        cluster_name_out <- if (!is.null(cluster_name) && nzchar(cluster_name)) {
            cluster_name
        } else if (!is.null(scope_out@meta.data) && ncol(scope_out@meta.data) > 0L) {
            tail(colnames(scope_out@meta.data), 1L)
        } else {
            NA_character_
        }
        assigned <- if (!is.null(scope_out@meta.data) &&
            !is.na(cluster_name_out) &&
            cluster_name_out %in% colnames(scope_out@meta.data)) {
            sum(!is.na(scope_out@meta.data[[cluster_name_out]]))
        } else {
            NA_integer_
        }
        total_genes <- if (!is.null(scope_out@meta.data)) nrow(scope_out@meta.data) else NA_integer_

        .log_key_values("clusterGenes", list(
            cluster_name = cluster_name_out,
            assigned = assigned,
            total_genes = total_genes,
            kept_genes = if (!is.null(report$stage1)) report$stage1$kept_genes else NA_integer_,
            stage1_mode = if (!is.null(report$stage1)) report$stage1$mode else NA_character_,
            q_fallback = if (!is.null(report$q_distribution)) report$q_distribution$fallback_method else NA_character_
        ))

        if (isTRUE(requested_return_report)) {
            return(result)
        }
        .log_done("clusterGenes", verbose = verbose)
        scope_out
    })
}

#' Enumerate dendrogram walk paths.
#' @description
#' Extracts walk paths from a dendrogram network object produced by geneSCOPE
#' clustering. This is commonly used for downstream plotting/inspection.
#' @param dnet_obj A dendrogram network list (typically produced by clustering/plotting helpers).
#' @param gene Primary gene identifier to trace.
#' @param gene2 Optional second gene identifier (two-ended walk mode).
#' @param all_shortest Deprecated placeholder; retained for API compatibility.
#' @param cutoff Maximum walk length to consider when enumerating simple paths.
#' @param max_paths Maximum number of paths to return (guards combinatorial explosion).
#' @param verbose Logical; whether to emit compact public progress messages.
#'   Default is \code{TRUE}. Set to \code{FALSE} for silent operation, or use
#'   \code{options(geneSCOPE.verbose = FALSE)} for global control.
#' @return A list containing `paths` (character vectors) and `paths_df` (summary data.frame).
#' @examples
#' \dontrun{
#' scope_obj <- clusterGenes(scope_obj, grid_name = "grid30", pct_min = "q50", n_start = 100, consensus_thr = 0.9, resolution = 0.5)
#' dnet_obj <- plotDendroNetwork(scope_obj, grid_name = "grid30", cluster_name = "gene_cluster", graph_slot_name = "g_consensus")
#' paths <- getDendroWalkPaths(dnet_obj, gene = "ACTA2", cutoff = 6)
#' head(paths$paths_df)
#' }
#' @seealso `plotDendroNetwork()`
#' @export
getDendroWalkPaths <- function(
    dnet_obj,
    gene,
    gene2 = NULL,
    all_shortest = FALSE,
    cutoff = 8,
    max_paths = 50000,
    verbose = getOption("geneSCOPE.verbose", TRUE)) {
    .get_dendro_walk_paths(
        dnet_obj = dnet_obj,
        gene = gene,
        gene2 = gene2,
        all_shortest = all_shortest,
        cutoff = cutoff,
        max_paths = max_paths,
        verbose = verbose
    )
}
