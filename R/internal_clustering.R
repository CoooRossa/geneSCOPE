#' ANSI Red Helper
#' @description
#' Wraps a message in ANSI red for terminal emphasis.
#' @param text Message string.
#' @return Scalar character string.
#' @keywords internal
.ansi_red <- function(text) {
    paste0("\033[31m", text, "\033[39m")
}

#' Extract Similarity Values For Cluster Diagnostics
#' @description
#' Collects finite upper-triangle similarity values without densifying sparse
#' matrices.
#' @param similarity_matrix Similarity matrix.
#' @return Numeric vector.
#' @keywords internal
.cluster_similarity_values <- function(similarity_matrix) {
    if (is.null(similarity_matrix)) return(numeric())
    if (inherits(similarity_matrix, "sparseMatrix")) {
        triplet <- summary(similarity_matrix)
        if (!all(c("i", "j", "x") %in% colnames(triplet))) {
            colnames(triplet) <- c("i", "j", "x")[seq_len(ncol(triplet))]
        }
        vals <- triplet$x[triplet$i < triplet$j]
    } else {
        vals <- similarity_matrix[upper.tri(similarity_matrix, diag = FALSE)]
    }
    vals <- as.numeric(vals)
    vals[is.finite(vals)]
}

#' Summarize Cluster Q Distribution
#' @description
#' Computes quantile-span diagnostics used to guard unstable clustering inputs.
#' @param similarity_matrix Similarity matrix.
#' @param q_span_threshold Minimum allowed q95-q05 spread.
#' @return Diagnostic list.
#' @keywords internal
.cluster_q_distribution_summary <- function(similarity_matrix,
                                           q_span_threshold = 0.1) {
    values <- .cluster_similarity_values(similarity_matrix)
    nonzero_values <- values[values != 0]
    if (!length(nonzero_values)) nonzero_values <- values
    if (!length(nonzero_values)) {
        return(list(
            ok = FALSE,
            reason = "no_finite_similarity_values",
            n_values = 0L,
            q05 = NA_real_,
            q10 = NA_real_,
            q90 = NA_real_,
            q95 = NA_real_,
            span90_10 = NA_real_,
            span95_05 = NA_real_,
            threshold = as.numeric(q_span_threshold)
        ))
    }
    qs <- stats::quantile(
        nonzero_values,
        probs = c(0.05, 0.10, 0.90, 0.95),
        na.rm = TRUE,
        names = FALSE,
        type = 8
    )
    span90_10 <- unname(qs[[3]] - qs[[2]])
    span95_05 <- unname(qs[[4]] - qs[[1]])
    list(
        ok = is.finite(span95_05) && span95_05 >= as.numeric(q_span_threshold),
        reason = if (is.finite(span95_05) && span95_05 >= as.numeric(q_span_threshold)) "ok" else "narrow_q_distribution",
        n_values = length(nonzero_values),
        q05 = unname(qs[[1]]),
        q10 = unname(qs[[2]]),
        q90 = unname(qs[[3]]),
        q95 = unname(qs[[4]]),
        span90_10 = span90_10,
        span95_05 = span95_05,
        threshold = as.numeric(q_span_threshold)
    )
}

#' Resolve Cluster Q Guard
#' @description
#' Decides whether a narrow similarity distribution should trigger the
#' compatibility fallback contract labelled as `k_perms`.
#' @param similarity_matrix Similarity matrix.
#' @param algo Requested stage-1 algorithm.
#' @param stage2_algo Optional requested stage-2 algorithm.
#' @param K Optional hotspot neighbourhood size.
#' @param use_consensus Whether consensus graphing is enabled.
#' @param n_restart Number of consensus restarts.
#' @param q_guard_method Guardrail policy (`auto`, `k_perms`, `none`).
#' @param q_span_threshold Minimum allowed q95-q05 spread.
#' @return List containing effective runtime knobs and diagnostics.
#' @keywords internal
.resolve_cluster_q_guard <- function(similarity_matrix,
                                     algo,
                                     stage2_algo = NULL,
                                     K = NULL,
                                     use_consensus = TRUE,
                                     n_restart = 100L,
                                     q_guard_method = c("auto", "k_perms", "none"),
                                     q_span_threshold = 0.1) {
    q_guard_method <- match.arg(q_guard_method)
    diag <- .cluster_q_distribution_summary(
        similarity_matrix,
        q_span_threshold = q_span_threshold
    )
    triggered <- !identical(q_guard_method, "none") && !isTRUE(diag$ok)
    effective_algo <- algo
    effective_stage2_algo <- if (is.null(stage2_algo)) NULL else stage2_algo
    effective_K <- K
    effective_use_consensus <- use_consensus
    effective_n_restart <- n_restart
    fallback_method <- "none"

    if (isTRUE(triggered)) {
        fallback_method <- "k_perms"
        effective_algo <- "hotspot-like"
        effective_stage2_algo <- "hotspot-like"
        if (is.null(effective_K) || !is.finite(effective_K) || effective_K < 2) {
            graph_size <- tryCatch(nrow(similarity_matrix), error = function(e) 10L)
            effective_K <- max(2L, min(30L, as.integer(floor(sqrt(max(4L, graph_size))))))
        } else {
            effective_K <- as.integer(effective_K)
        }
        effective_use_consensus <- FALSE
        effective_n_restart <- 1L
    }

    diag$requested_method <- q_guard_method
    diag$triggered <- triggered
    diag$fallback_method <- fallback_method
    diag$effective_algo <- effective_algo
    diag$effective_stage2_algo <- if (is.null(effective_stage2_algo)) NA_character_ else effective_stage2_algo
    diag$effective_K <- if (is.null(effective_K)) NA_integer_ else as.integer(effective_K)
    diag$effective_use_consensus <- isTRUE(effective_use_consensus)
    diag$effective_n_restart <- as.integer(effective_n_restart)

    list(
        diag = diag,
        algo = effective_algo,
        stage2_algo = effective_stage2_algo,
        K = effective_K,
        use_consensus = effective_use_consensus,
        n_restart = effective_n_restart
    )
}

#' Cluster genes on a consensus graph
#' @description
#' Internal helper for `.cluster_genes`.
#' Runs the stage1/stage2 clustering pipeline on a scope object.
#' @param scope_obj A `scope_object` containing graph/stat layers.
#' @param grid_name Optional grid layer name; defaults to the active layer.
#' @param stats_layer Lee statistic layer name (overridden by `lee_stats_layer`).
#' @param similarity_slot Slot name holding similarity matrix (default `L`).
#' @param significance_slot Slot name holding significance/FDR matrix.
#' @param min_cutoff Minimum similarity cutoff (alias `L_min`).
#' @param use_significance Whether to apply significance/FDR gating.
#' @param resolution Clustering resolution parameter passed to Leiden/Louvain.
#' @param algo Clustering algorithm choice (`leiden`, `louvain`, `hotspot-like`).
#' @param graph_slot_name Graph slot used to drive clustering (e.g., consensus graph).
#' @param return_report When `TRUE`, returns both updated scope object and a report list.
#' @param lee_stats_layer Layer name.
#' @param L_min Numeric threshold.
#' @param use_FDR Logical flag.
#' @param significance_max Numeric threshold.
#' @param FDR_max Numeric threshold.
#' @param pct_min Numeric threshold.
#' @param drop_isolated Parameter value.
#' @param stage2_algo Parameter value.
#' @param gamma Parameter value.
#' @param objective Parameter value.
#' @param cluster_name Parameter value.
#' @param use_log1p_weight Logical flag.
#' @param use_consensus Logical flag.
#' @param n_restart Parameter value.
#' @param consensus_thr Parameter value.
#' @param K Parameter value.
#' @param min_module_size Numeric threshold.
#' @param CI95_filter Parameter value.
#' @param CI_rule Parameter value.
#' @param curve_layer Layer name.
#' @param use_mh_weight Logical flag.
#' @param use_cmh_weight Logical flag.
#' @param mh_slot Slot name.
#' @param cmh_slot Slot name.
#' @param post_smooth Parameter value.
#' @param post_smooth_quant Parameter value.
#' @param post_smooth_power Parameter value.
#' @param enable_subcluster Logical flag.
#' @param sub_min_size Parameter value.
#' @param sub_min_child_size Parameter value.
#' @param sub_resolution_factor Parameter value.
#' @param sub_within_cons_max Numeric threshold.
#' @param sub_conductance_min Numeric threshold.
#' @param sub_improve_within_cons_min Numeric threshold.
#' @param sub_max_groups Parameter value.
#' @param enable_qc_filter Logical flag.
#' @param qc_gene_intra_cons_min Numeric threshold.
#' @param qc_gene_best_out_cons_min Numeric threshold.
#' @param qc_gene_intra_weight_q Parameter value.
#' @param keep_cross_stable Logical flag.
#' @param min_cluster_size Numeric threshold.
#' @param keep_stage1_backbone Logical flag.
#' @param backbone_floor_q Parameter value.
#' @param verbose Logical; whether to emit progress messages.
#' @param ncores Number of cores/threads to use.
#' @param mode Parameter value.
#' @param large_n_threshold Parameter value.
#' @param nk_condaenv Parameter value.
#' @param nk_conda_bin Parameter value.
#' @param nk_python Parameter value.
#' @param nk_leiden_iterations Parameter value.
#' @param nk_leiden_randomize Parameter value.
#' @param aggr_future_workers Parameter value.
#' @param aggr_batch_size Parameter value.
#' @param future_globals_min_bytes Parameter value.
#' @param profile_timing Parameter value.
#' @param backend Graph backend contract.
#' @param q_guard_method Q-distribution guardrail policy.
#' @param q_span_threshold Minimum q95-q05 span allowed before fallback.
#' @return Return value used internally.
#' @keywords internal
.cluster_genes <- function(
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
    large_n_threshold = 1000L,
    nk_condaenv = NULL,
    nk_conda_bin = NULL,
    nk_python = NULL,
    nk_leiden_iterations = 10L,
    nk_leiden_randomize = TRUE,
    aggr_future_workers = 2L,
    aggr_batch_size = NULL,
    future_globals_min_bytes = 2 * 1024^3,
    profile_timing = FALSE,
    backend = c("auto", "igraph", "networkit"),
    q_guard_method = c("auto", "k_perms", "none"),
    q_span_threshold = 0.1) {
    step01 <- .log_step("clusterGenes", "S01", "normalize inputs and build config", verbose)
    step01$enter()

    if (length(algo) == 1L && identical(tolower(algo), "hotspot")) {
        algo <- "hotspot-like"
    }
    algo <- match.arg(algo)
    objective <- match.arg(objective)
    valid_ci_rules <- c("remove_within", "remove_outside")
    ci_rule_value <- NULL
    if (!is.null(CI_rule) && length(CI_rule)) {
        ci_rule_value <- match.arg(CI_rule, valid_ci_rules)
    }
    mode <- match.arg(mode)
    backend <- match.arg(backend)
    q_guard_method <- match.arg(q_guard_method)

    requested_mode <- mode
    if (identical(mode, "safe_sequential")) {
        mode <- "safe"
        ncores <- 1L
        aggr_future_workers <- 1L
        if (requireNamespace("future", quietly = TRUE)) {
            old_future_plan <- tryCatch(future::plan(), error = function(e) NULL)
            future::plan(future::sequential)
            on.exit({
                if (!is.null(old_future_plan)) {
                    try(future::plan(old_future_plan), silent = TRUE)
                }
            }, add = TRUE)
        }
        .log_key(
            "clusterGenes",
            .ansi_red("safe_sequential mode selected; using single-worker safe clustering to avoid future export pressure on large graphs."),
            level = "WARN"
        )
    }
    if (identical(backend, "igraph")) {
        mode <- if (identical(requested_mode, "fast")) "fast" else "safe"
    } else if (identical(backend, "networkit")) {
        if (!.networkit_available()) {
            stop("clusterGenes backend='networkit' requested but NetworKit is unavailable.")
        }
        mode <- "aggressive"
    }

    # Use requested cores directly (clamped to available logical cores)
    ncores_requested <- ncores
    avail_cores <- max(1L, detectCores(logical = TRUE))
    ncores_use <- if (is.null(ncores) || is.na(ncores)) {
        avail_cores
    } else if (identical(requested_mode, "safe_sequential")) {
        1L
    } else {
        max(1L, min(as.integer(ncores), avail_cores))
    }
    ncores <- ncores_use

    if (!is.null(lee_stats_layer)) stats_layer <- lee_stats_layer
    if (!is.null(L_min)) min_cutoff <- L_min
    if (!is.null(use_FDR)) use_significance <- use_FDR
    if (!is.null(FDR_max)) significance_max <- FDR_max
    if (!is.null(use_cmh_weight)) use_mh_weight <- use_cmh_weight
    if (!is.null(cmh_slot)) mh_slot <- cmh_slot
    normalize_name <- function(val) {
        if (is.null(val)) return(NULL)
        val <- as.character(val)[1]
        if (is.na(val) || !nzchar(val)) return(NULL)
        val
    }
    mh_slot <- normalize_name(mh_slot)
    curve_layer <- normalize_name(curve_layer)
    if (is.null(use_mh_weight)) {
        use_mh_weight <- !is.null(mh_slot)
    } else {
        use_mh_weight <- isTRUE(use_mh_weight)
        if (use_mh_weight && is.null(mh_slot)) {
            stop("mh_slot must be provided when enabling MH weighting.")
        }
    }
    if (!use_mh_weight) {
        mh_slot <- NULL
    }
    if (is.null(ci_rule_value) && !is.null(CI95_filter)) {
        if (isTRUE(CI95_filter)) {
            warning("CI95_filter is deprecated; please specify CI_rule instead.", call. = FALSE)
            ci_rule_value <- "remove_within"
        }
    }
    CI95_filter_flag <- !is.null(ci_rule_value)
    if (CI95_filter_flag && is.null(curve_layer)) {
        stop("curve_layer must be provided when specifying CI_rule.")
    }
    if (!CI95_filter_flag) {
        curve_layer <- NULL
    }
    if (!is.null(stage2_algo) && length(stage2_algo) == 1L && identical(tolower(stage2_algo), "hotspot")) {
        stage2_algo <- "hotspot-like"
    }

    sys_mem_bytes <- .detect_system_memory_bytes()

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
        n_threads = ncores_use,
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
        CI95_filter = CI95_filter_flag,
        CI_rule = ci_rule_value,
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
        hotspot_min_module_size = min_module_size,
        future_globals_min_bytes = as.numeric(future_globals_min_bytes),
        system_memory_bytes = sys_mem_bytes,
        execution_mode = requested_mode
    )
    cfg$system_memory_bytes <- if (is.finite(sys_mem_bytes) && sys_mem_bytes > 0) sys_mem_bytes else cfg$system_memory_bytes
    if (.log_enabled(verbose)) {
        .log_backend(
            "clusterGenes",
            "S01",
            "graph_backend",
            backend,
            reason = paste0("requested_mode=", requested_mode, " effective_mode=", cfg$mode),
            verbose = verbose
        )
    }
    step01$done(sprintf("algo=%s objective=%s mode=%s", algo, objective, cfg$mode))

    step02 <- .log_step("clusterGenes", "S02", "configure threads and resources", verbose)
    step02$enter()
    if (.log_enabled(verbose)) {
        if (is.null(ncores_requested) || is.na(ncores_requested)) {
            reason <- paste0("requested=auto available=", avail_cores)
        } else if (as.integer(ncores_requested) != ncores_use) {
            reason <- paste0("requested=", as.integer(ncores_requested), " available=", avail_cores)
        } else {
            reason <- paste0("requested=", as.integer(ncores_requested))
        }
        .log_backend("clusterGenes", "S02", "openmp_threads", ncores_use, reason = reason, verbose = verbose)
    }
    sys_mem_gb <- if (is.finite(sys_mem_bytes) && sys_mem_bytes > 0) sys_mem_bytes / 1024^3 else .get_system_memory_gb()
    if (.log_enabled(verbose)) {
        .log_info("clusterGenes", "S02", sprintf("system_memory_gb=%.1f", sys_mem_gb), verbose)
    }
    step02$done(sprintf("ncores=%d", ncores_use))

    timing <- list()
    .tic <- function() proc.time()
    .toc <- function(t0) as.numeric((proc.time() - t0)["elapsed"])
    t_all0 <- .tic()

    step03 <- .log_step("clusterGenes", "S03", "extract similarity inputs", verbose)
    step03$enter()
    inputs <- .extract_scope_layers(scope_obj, grid_name,
        stats_layer = stats_layer,
        similarity_slot = cfg$similarity_slot,
        significance_slot = cfg$significance_slot,
        use_significance = cfg$use_significance,
        verbose = verbose)
    sim_dim <- nrow(inputs$similarity)
    if (.log_enabled(verbose)) {
        sig_state <- if (is.null(inputs$significance)) "none" else "present"
        .log_info("clusterGenes", "S03", sprintf("similarity_dim=%d significance=%s", sim_dim, sig_state), verbose)
    }
    q_guard <- .resolve_cluster_q_guard(
        similarity_matrix = inputs$similarity,
        algo = cfg$algo,
        stage2_algo = cfg$stage2_algo,
        K = cfg$hotspot_k,
        use_consensus = cfg$use_consensus,
        n_restart = cfg$n_restart,
        q_guard_method = q_guard_method,
        q_span_threshold = q_span_threshold
    )
    if (.log_enabled(verbose)) {
        .log_info(
            "clusterGenes",
            "S03",
            sprintf(
                "q05=%.4f q95=%.4f span95_05=%.4f threshold=%.4f",
                q_guard$diag$q05,
                q_guard$diag$q95,
                q_guard$diag$span95_05,
                q_guard$diag$threshold
            ),
            verbose
        )
    }
    if (isTRUE(q_guard$diag$triggered)) {
        cfg$algo <- q_guard$algo
        cfg$stage2_algo <- q_guard$stage2_algo
        cfg$hotspot_k <- q_guard$K
        cfg$use_consensus <- q_guard$use_consensus
        cfg$n_restart <- q_guard$n_restart
        .log_key(
            "clusterGenes",
            .ansi_red(
                paste0(
                    "q distribution span95_05=",
                    sprintf("%.4f", q_guard$diag$span95_05),
                    " < ",
                    sprintf("%.4f", q_guard$diag$threshold),
                    "; using k_perms fallback contract (mapped to hotspot-like runtime path). ",
                    "Interpret cluster/modules cautiously."
                )
            ),
            level = "WARN"
        )
    }
    step03$done()

    # Memory guard: estimate dense footprint of similarity/significance matrices per thread
    step04 <- .log_step("clusterGenes", "S04", "select fast/safe mode", verbose)
    step04$enter()
    base_gb <- (sim_dim^2 * 8 * 2) / (1024^3) # two matrices (similarity + significance)
    est_total_gb <- base_gb * ncores_use
    if (est_total_gb > sys_mem_gb && !identical(requested_mode, "safe_sequential")) {
        stop(
            "[.cluster_genes] Estimated memory requirement (", round(est_total_gb, 1),
            " GB) exceeds system capacity (", round(sys_mem_gb, 1),
            " GB). Reduce ncores, gene count, or run clusterGenes(mode='safe_sequential')."
        )
    }

    fast_limit <- cfg$system_memory_bytes
    if (!is.finite(fast_limit) || fast_limit <= 0) fast_limit <- cfg$future_globals_min_bytes
    fast_est <- if (identical(cfg$mode, "fast") || identical(cfg$mode, "auto")) .estimate_fast_mode_memory(inputs$similarity, cfg) else NA_real_
    margin <- getOption("genescope.fast_mode_mem_margin", 0.9)
    if (!is.finite(margin) || margin <= 0 || margin > 1) margin <- 0.9
    finite_limit <- is.finite(fast_limit) && fast_limit > 0
    fast_allowed <- (identical(cfg$mode, "fast") || identical(cfg$mode, "auto")) &&
        (!finite_limit || !is.finite(fast_est) || fast_est <= fast_limit * margin)
    cfg$prefer_fast <- FALSE
    if (fast_allowed && (identical(cfg$mode, "fast") || identical(cfg$mode, "auto"))) {
        cfg$prefer_fast <- TRUE
        cfg$mode <- "fast"
        if (is.finite(fast_est) && fast_est > 0) {
            fg_multiplier <- getOption("genescope.fast_mode_fg_multiplier", 1.2)
            if (!is.finite(fg_multiplier) || fg_multiplier <= 0) fg_multiplier <- 1.2
            cfg$future_globals_min_bytes <- max(cfg$future_globals_min_bytes, fast_est * fg_multiplier)
        }
    } else if (identical(cfg$mode, "fast") && !fast_allowed) {
        limit_gb <- if (is.finite(fast_limit)) fast_limit / 1024^3 else NA_real_
        est_gb <- if (is.finite(fast_est)) fast_est / 1024^3 else NA_real_
        warning(sprintf("[cluster] fast mode requires %.2f GB but limit is %.2f GB; falling back to safe mode.", est_gb, limit_gb))
        cfg$mode <- "safe"
    }
    step04$done(sprintf("mode=%s prefer_fast=%s", cfg$mode, cfg$prefer_fast))

    step05 <- .log_step("clusterGenes", "S05", "stage1 consensus workflow", verbose)
    step05$enter()
    stage1 <- .stage1_consensus_workflow_v2(
        inputs$similarity,
        .threshold_similarity_graph(inputs$similarity, inputs$significance, cfg,
            verbose = verbose
        ),
        cfg,
        verbose = verbose,
        networkit_config = list(conda_env = nk_condaenv, conda_bin = nk_conda_bin, python_path = nk_python)
    )
    timing$stage1 <- .toc(t_all0)
    step05$done(sprintf("kept_genes=%d", length(stage1$kept_genes)))

    step06 <- .log_step("clusterGenes", "S06", "prepare MH/CI inputs", verbose)
    step06$enter()
    if (.log_enabled(verbose)) {
        .log_info("clusterGenes", "S06", paste0("use_mh_weight=", cfg$use_mh_weight, " CI95_filter=", cfg$CI95_filter), verbose)
    }
    mh_object <- NULL
    if (cfg$use_mh_weight) {
        lstats <- inputs$stats_object
        mh_object <- lstats[[cfg$mh_slot]]
        if (is.null(mh_object)) stop("MH slot missing: ", cfg$mh_slot)
        if (.log_enabled(verbose)) {
            mh_type <- if (inherits(mh_object, "igraph")) {
                "igraph"
            } else if (inherits(mh_object, "sparseMatrix")) {
                "sparse_matrix"
            } else if (is.matrix(mh_object)) {
                "dense_matrix"
            } else {
                class(mh_object)[1]
            }
            .log_backend("clusterGenes", "S06", "mh_weight",
                paste0("slot=", cfg$mh_slot, " type=", mh_type),
                verbose = verbose
            )
        }
    } else if (.log_enabled(verbose)) {
        .log_backend("clusterGenes", "S06", "mh_weight", "disabled",
            reason = "use_mh_weight=FALSE",
            verbose = verbose
        )
    }

    pearson_mat <- NULL
    if (cfg$CI95_filter) {
        pearson_mat <- try(.get_pearson_matrix(scope_obj, grid_name = inputs$grid_name, level = "cell"), silent = TRUE)
        if (inherits(pearson_mat, "try-error") || is.null(pearson_mat)) {
            pearson_mat <- try(.get_pearson_matrix(scope_obj, grid_name = inputs$grid_name, level = "grid"), silent = TRUE)
        }
        if (inherits(pearson_mat, "try-error") || is.null(pearson_mat)) {
            stop("CI95_filter requested but Pearson matrix not found at cell/grid level.")
        }
    }
    step06$done()

    step07 <- .log_step("clusterGenes", "S07", "stage2 refinement", verbose)
    step07$enter()
    stage2 <- .stage2_refine_workflow_v2(stage1,
        similarity_matrix = inputs$similarity,
        config = cfg,
        FDR = inputs$significance,
        mh_object = mh_object,
        aux_stats = inputs$stats_object,
        pearson_matrix = pearson_mat,
        verbose = verbose)
    timing$total <- .toc(t_all0)
    n_clusters <- length(unique(na.omit(stage2$membership)))
    step07$done(sprintf("clusters=%d kept_genes=%d", n_clusters, length(stage2$kept_genes)))

    step08 <- .log_step("clusterGenes", "S08", "write clustering outputs", verbose)
    step08$enter()
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

    if (.log_enabled(verbose)) {
        assigned <- sum(!is.na(scope_obj@meta.data[genes_all, cluster_name]))
        pct <- if (length(genes_all)) sprintf("%.1f%%", 100 * assigned / length(genes_all)) else "0.0%"
        .log_info("clusterGenes", "S08", paste0("cluster_name=", cluster_name, " assigned=", assigned, "/", length(genes_all), " (", pct, ")"), verbose)
    }

    # Update stats layer with final consensus graph
    gname <- inputs$grid_name
    if (is.null(scope_obj@stats[[gname]]) || !is.list(scope_obj@stats[[gname]])) {
        scope_obj@stats[[gname]] <- list()
    }
    if (is.null(scope_obj@stats[[gname]][[stats_layer]])) {
        scope_obj@stats[[gname]][[stats_layer]] <- list()
    } else if (is.matrix(scope_obj@stats[[gname]][[stats_layer]]) ||
        inherits(scope_obj@stats[[gname]][[stats_layer]], "Matrix")) {
        stats_mat <- scope_obj@stats[[gname]][[stats_layer]]
        scope_obj@stats[[gname]][[stats_layer]] <- list()
        scope_obj@stats[[gname]][[stats_layer]][[cfg$similarity_slot]] <- stats_mat
    }
    scope_obj@stats[[gname]][[stats_layer]][[graph_slot_name]] <- stage2$consensus_graph
    step08$done()

    step09 <- .log_step("clusterGenes", "S09", "finalize report and return", verbose)
    step09$enter()
    report <- NULL
    if (isTRUE(return_report)) {
        report <- list()
        report$timing <- timing
        report$stage1 <- list(
            mode = stage1$mode_selected,
            kept_genes = length(stage1$kept_genes),
            algo = cfg$algo,
            stage2_algo = if (is.null(cfg$stage2_algo)) cfg$algo else cfg$stage2_algo,
            backend_requested = backend,
            execution_mode_requested = requested_mode,
            execution_mode_effective = cfg$mode,
            ncores_effective = ncores_use
        )
        report$q_distribution <- q_guard$diag
        report$metrics_before <- stage2$metrics_before
        report$metrics_after <- stage2$metrics_after
        report$genes_before <- stage2$genes_before
        report$genes_after <- stage2$genes_after
    }

    result <- if (isTRUE(return_report)) {
        list(scope_obj = scope_obj, report = report)
    } else {
        scope_obj
    }
    step09$done(sprintf("return_report=%s", isTRUE(return_report)))
    result
}

#' Extract walk paths from a dendrogram network object
#' @description
#' Internal helper for `.get_dendro_walk_paths`.
#' @param dnet_obj A dendrogram network list (output from clustering).
#' @param gene Primary gene identifier to trace.
#' @param gene2 Optional second gene for two-ended walks.
#' @param cutoff Maximum walk length to consider.
#' @param max_paths Maximum number of paths to return.
#' @param verbose Emit progress messages when TRUE.
#' @param all_shortest Parameter value.
#' @return List of walk paths for downstream plotting/inspection.
#' @keywords internal
.get_dendro_walk_paths <- function(dnet_obj,
                               gene,
                               gene2 = NULL,
                               all_shortest = FALSE, # deprecated; placeholder only
                               cutoff = 8,
                               max_paths = 50000,
                               verbose = getOption("geneSCOPE.verbose", TRUE)) {
    ## 1. Locate the igraph object --------------------------------------------------
    g <- NULL
    for (nm in c("graph", "g", "igraph")) {
        if (!is.null(dnet_obj[[nm]])) {
            g <- dnet_obj[[nm]]
            break
        }
    }
    if (is.null(g) || !inherits(g, "igraph")) stop("No igraph graph found (fields graph/g/igraph)")

    vn <- igraph::V(g)$name
    if (is.null(vn)) stop("Graph vertices are missing names (V(g)$name)")

    gene <- as.character(gene)[1]
    if (!(gene %in% vn)) stop("gene not found in graph: ", gene)
    if (!is.null(gene2)) {
        gene2 <- as.character(gene2)[1]
        if (!(gene2 %in% vn)) stop("gene2 not found in graph: ", gene2)
        if (gene2 == gene) stop("gene and gene2 must be different")
    }

    ## 2. Normalize random-walk paths ----------------------------------------------
    has_walks <- FALSE
    walks_std <- NULL
    if (!is.null(dnet_obj$walks) && is.list(dnet_obj$walks)) {
        rw <- dnet_obj$walks
        rw <- rw[vapply(rw, function(x) length(x) >= 2, logical(1))]
        if (length(rw)) {
            walks_std <- lapply(rw, function(w) {
                if (is.factor(w)) w <- as.character(w)
                if (is.numeric(w)) {
                    if (all(w >= 1 & w <= length(vn))) vn[w] else as.character(w)
                } else {
                    as.character(w)
                }
            })
            has_walks <- TRUE
        }
    }

    ## 3. Helper: enumerate simple paths and filter --------------------------------
    enumerate_paths_filter <- function(filter_fun) {
        # filter_fun: function(char_vector_path) -> TRUE/FALSE
        nV <- length(vn)
        kept <- list()
        k <- 0L
        # source<target to avoid duplicates from reversed direction
        for (s in seq_len(nV - 1L)) {
            for (t in (s + 1L):nV) {
                ps <- igraph::all_simple_paths(g, from = s, to = t, mode = "all", cutoff = cutoff)
                if (length(ps)) {
                    for (pp in ps) {
                        path_chr <- vn[pp]
                        if (filter_fun(path_chr)) {
                            k <- k + 1L
                            kept[[k]] <- path_chr
                            if (k >= max_paths) {
                                if (verbose) message("[geneSCOPE] Reached max_paths=", max_paths, "; stopping enumeration")
                                return(kept)
                            }
                        }
                    }
                }
            }
        }
        kept
    }

    ## 4A. Single-gene mode: all paths containing `gene` ----------------------------
    if (is.null(gene2)) {
        if (has_walks) {
            sel <- walks_std[vapply(walks_std, function(p) gene %in% p, logical(1))]
            if (verbose) message("[geneSCOPE] Selected ", length(sel), " paths containing ", gene, " from random walks")
        } else {
            if (verbose) message("[geneSCOPE] No random walks found; enumerating graph paths (cutoff=", cutoff, ", max_paths=", max_paths, ")")
            sel <- enumerate_paths_filter(function(p) gene %in% p)
            if (verbose) message("[geneSCOPE] Enumerated and filtered ", length(sel), " paths")
        }
        # de-duplicate
        if (length(sel)) {
            key <- vapply(sel, paste, collapse = "->", FUN.VALUE = character(1))
            sel <- sel[!duplicated(key)]
        }
        paths_df <- if (length(sel)) {
            data.frame(
                path_id = seq_along(sel),
                length = vapply(sel, length, integer(1)),
                path_str = vapply(sel, paste, collapse = "->", FUN.VALUE = character(1)),
                stringsAsFactors = FALSE
            )
        } else {
            data.frame(path_id = integer(0), length = integer(0), path_str = character(0))
        }
        return(list(
            mode = "single_gene",
            genes = gene,
            paths = sel,
            paths_df = paths_df,
            n_paths = nrow(paths_df),
            source_has_walks = has_walks,
            cutoff_used = if (has_walks) NA_integer_ else cutoff
        ))
    }

    ## 4B. Dual-gene mode: all paths containing both `gene` and `gene2` -------------
    if (has_walks) {
        sel <- walks_std[vapply(walks_std, function(p) all(c(gene, gene2) %in% p), logical(1))]
        if (verbose) message("[geneSCOPE] Selected ", length(sel), " paths containing both ", gene, " and ", gene2, " from random walks")
    } else {
        if (verbose) message("[geneSCOPE] No random walks found; enumerating and filtering paths containing both genes (cutoff=", cutoff, ", max_paths=", max_paths, ")")
        sel <- enumerate_paths_filter(function(p) (gene %in% p) && (gene2 %in% p))
        if (verbose) message("[geneSCOPE] Enumerated and filtered ", length(sel), " paths")
    }
    # de-duplicate
    if (length(sel)) {
        key <- vapply(sel, paste, collapse = "->", FUN.VALUE = character(1))
        sel <- sel[!duplicated(key)]
    }
    paths_df <- if (length(sel)) {
        data.frame(
            path_id = seq_along(sel),
            length = vapply(sel, length, integer(1)),
            path_str = vapply(sel, paste, collapse = "->", FUN.VALUE = character(1)),
            stringsAsFactors = FALSE
        )
    } else {
        data.frame(path_id = integer(0), length = integer(0), path_str = character(0))
    }

    list(
        mode = "dual_gene_all_paths",
        genes = c(gene, gene2),
        paths = sel,
        paths_df = paths_df,
        n_paths = nrow(paths_df),
        source_has_walks = has_walks,
        cutoff_used = if (has_walks) NA_integer_ else cutoff
    )
}
