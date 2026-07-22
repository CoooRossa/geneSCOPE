#' Resolve weight_style attribute for Lee's L weights matrix
#' @description
#'   Internal helper that resolves the weight_style attribute from a spatial
#'   weights matrix W. If the attribute is missing or unsupported, infers the
#'   style from the matrix values (binary vs row-standardised).
#' @param W A spatial weights matrix (sparse or dense).
#' @param caller Character. Name of the calling function for warning messages.
#' @param verbose Logical. Whether to emit progress messages.
#' @return A list with components: weight_style ("B" or "W"), source ("attribute" or "inferred"),
#'   and inferred (logical indicating if style was inferred).
#' @importFrom stats sd pnorm p.adjust coef lm
#' @importFrom igraph graph_from_adjacency_matrix simplify degree cluster_louvain cluster_leiden components modularity
#' @keywords internal
.resolve_leesl_weight_style <- function(W,
                                        caller = "computeL",
                                        verbose = TRUE) {
    raw_style <- attr(W, "weight_style")
    style <- .transcripts_attestation_scalar(raw_style)
    style_upper <- toupper(style)

    if (identical(style_upper, "B") || identical(style_upper, "W")) {
        return(list(
            weight_style = style_upper,
            source = "attribute",
            inferred = FALSE
        ))
    }

    w_nonzero <- if (inherits(W, "Matrix")) W@x else W[W != 0]
    inferred <- if (length(w_nonzero) && all(w_nonzero %in% c(0, 1, 0L, 1L))) "B" else "W"
    if (is.na(style)) {
        warning(
            caller, ": weight_style attribute missing; inferred style='", inferred,
            "' from the W matrix values.",
            call. = FALSE
        )
    } else {
        warning(
            caller, ": unsupported weight_style='", style,
            "'; inferred style='", inferred, "' from the W matrix values.",
            call. = FALSE
        )
    }
    .log_info(caller, "S04", paste0(
        "weight_style resolved by inference: ", inferred
    ), verbose)
    list(
        weight_style = inferred,
        source = "inferred",
        inferred = TRUE
    )
}

#' Validate an integer-valued Lee workflow argument without truncation.
#' @keywords internal
.lee_validate_integer_scalar <- function(x, arg, lower = 0L) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        x < lower || x != floor(x) || x > .Machine$integer.max) {
        qualifier <- if (lower == 0L) "a non-negative integer" else paste0("an integer >= ", lower)
        stop(arg, " must be ", qualifier, ".", call. = FALSE)
    }
    as.integer(x)
}

#' Approximate q-values from low-resolution permutation p-values.
#'
#' Uses beta posterior-smoothed permutation p-values before Storey-style
#' monotone q-value estimation. When `perms` is supplied, p-values are treated
#' as Monte-Carlo values on the `(k + 1) / (B + 1)` grid and converted to
#' `(k + 1) / (B + 2)`, the posterior mean under a uniform beta-binomial model.
#' This avoids zero/one boundary artifacts and reduces staircase behavior for
#' low permutation budgets.
#'
#' @keywords internal
approximate_q_values <- function(p_values,
                                 perms = NULL,
                                 method = c("beta"),
                                 lambdas = seq(0.5, 0.95, by = 0.05)) {
    method <- match.arg(method)
    dims <- dim(p_values)
    dimn <- dimnames(p_values)
    pairwise <- length(dims) == 2L && dims[1L] == dims[2L]
    pair_mask <- if (pairwise) upper.tri(p_values, diag = FALSE) else NULL
    p <- if (pairwise) as.numeric(p_values[pair_mask]) else as.numeric(p_values)
    restore_shape <- function(values) {
        if (!pairwise) {
            return(if (is.null(dims)) values else matrix(values, nrow = dims[1L], dimnames = dimn))
        }
        ans <- matrix(NA_real_, dims[1L], dims[2L], dimnames = dimn)
        ans[pair_mask] <- values
        ans[lower.tri(ans)] <- t(ans)[lower.tri(ans)]
        diag(ans) <- NA_real_
        ans
    }
    out <- rep(NA_real_, length(p))
    finite <- is.finite(p)
    if (!any(finite)) {
        return(list(
            q_values = restore_shape(out),
            adjusted_p_values = restore_shape(out),
            pi0_hat = NA_real_,
            method = method,
            warning = "no finite p-values"
        ))
    }

    warn <- character()
    if (any(p[finite] < 0 | p[finite] > 1)) {
        warn <- c(warn, "p-values outside [0,1] were clamped")
    }
    p[finite] <- pmin(pmax(p[finite], 0), 1)

    p_adj <- p
    B <- if (is.null(perms)) NULL else .lee_validate_integer_scalar(perms, "perms", lower = 0L)
    if (!is.null(B) && B > 0L) {
        k_est <- round(p[finite] * (B + 1) - 1)
        k_est <- pmin(pmax(k_est, 0L), B)
        p_adj[finite] <- (k_est + 1) / (B + 2)
    } else {
        tiny <- .Machine$double.xmin
        p_adj[finite] <- pmin(pmax(p_adj[finite], tiny), 1 - tiny)
        warn <- c(warn, "perms missing or zero; beta grid smoothing skipped")
    }

    p_work <- p_adj[finite]
    lambdas <- lambdas[is.finite(lambdas) & lambdas > 0 & lambdas < 1]
    if (!length(lambdas)) lambdas <- seq(0.5, 0.95, by = 0.05)
    pi0_vals <- vapply(lambdas, function(lam) {
        mean(p_work > lam) / (1 - lam)
    }, numeric(1))
    pi0_hat <- suppressWarnings(min(1, min(pi0_vals, na.rm = TRUE)))
    if (!is.finite(pi0_hat) || pi0_hat <= 0) {
        pi0_hat <- 1
        warn <- c(warn, "pi0 estimate was unstable; using pi0=1")
    }

    m <- length(p_work)
    ord <- order(p_work, na.last = NA)
    ranked <- seq_along(ord)
    q_ord <- pi0_hat * m * p_work[ord] / pmax(ranked, 1)
    if (length(q_ord) > 1L) {
        for (i in seq(from = length(q_ord) - 1L, to = 1L, by = -1L)) {
            if (q_ord[i] > q_ord[i + 1L]) q_ord[i] <- q_ord[i + 1L]
        }
    }
    q_work <- numeric(m)
    q_work[ord] <- pmin(pmax(q_ord, 0), 1)
    out[finite] <- q_work

    list(
        q_values = restore_shape(out),
        adjusted_p_values = restore_shape(p_adj),
        pi0_hat = pi0_hat,
        method = method,
        warning = paste(unique(warn), collapse = "; ")
    )
}

#' Adjust a pairwise p-value matrix over unique unordered gene pairs.
#' @keywords internal
.lee_pairwise_p_adjust <- function(p_mat, method = "BH") {
    p_mat <- as.matrix(p_mat)
    out <- matrix(NA_real_, nrow(p_mat), ncol(p_mat), dimnames = dimnames(p_mat))
    if (nrow(p_mat) == ncol(p_mat)) {
        ut <- upper.tri(p_mat, diag = FALSE)
        vals <- p_mat[ut]
        good <- is.finite(vals)
        adj <- rep(NA_real_, length(vals))
        if (any(good)) adj[good] <- p.adjust(vals[good], method = method)
        out[ut] <- adj
        out[lower.tri(out)] <- t(out)[lower.tri(out)]
        diag(out) <- NA_real_
    } else {
        vals <- as.numeric(p_mat)
        good <- is.finite(vals)
        adj <- rep(NA_real_, length(vals))
        if (any(good)) adj[good] <- p.adjust(vals[good], method = method)
        out[] <- adj
    }
    out
}

#' Exploratory Storey q-values over the unique pairwise test family.
#' @keywords internal
.lee_pairwise_storey <- function(p_mat,
                                  lambdas = seq(0.5, 0.95, by = 0.05)) {
    p_mat <- as.matrix(p_mat)
    pairwise <- nrow(p_mat) == ncol(p_mat)
    mask <- if (pairwise) upper.tri(p_mat, diag = FALSE) else matrix(TRUE, nrow(p_mat), ncol(p_mat))
    p_all <- as.numeric(p_mat[mask])
    good <- is.finite(p_all)
    p_vec <- p_all[good]
    q_all <- rep(NA_real_, length(p_all))
    pi0_hat <- NA_real_

    if (length(p_vec)) {
        p_vec <- pmin(pmax(p_vec, 0), 1)
        pi0_vals <- vapply(lambdas, function(lam) {
            mean(p_vec > lam) / (1 - lam)
        }, numeric(1))
        pi0_hat <- min(1, min(pi0_vals, na.rm = TRUE))
        if (!is.finite(pi0_hat) || pi0_hat <= 0) pi0_hat <- 1

        ord <- order(p_vec)
        q_ord <- pi0_hat * length(p_vec) * p_vec[ord] / seq_along(ord)
        if (length(q_ord) > 1L) {
            for (i in seq.int(length(q_ord) - 1L, 1L, by = -1L)) {
                q_ord[i] <- min(q_ord[i], q_ord[i + 1L])
            }
        }
        q_vec <- numeric(length(p_vec))
        q_vec[ord] <- pmin(pmax(q_ord, 0), 1)
        q_all[good] <- q_vec
    }

    out <- matrix(NA_real_, nrow(p_mat), ncol(p_mat), dimnames = dimnames(p_mat))
    out[mask] <- q_all
    if (pairwise) {
        out[lower.tri(out)] <- t(out)[lower.tri(out)]
        diag(out) <- NA_real_
    }
    list(q_values = out, pi0_hat = pi0_hat, n_tests = sum(good))
}

#' Convert exact two-sided permutation p-values to signed normal scores.
#' @keywords internal
.lee_empirical_signed_z <- function(L, P) {
    L <- as.matrix(L)
    P <- as.matrix(P)
    if (!identical(dim(L), dim(P))) stop("L and P dimensions must match.", call. = FALSE)
    out <- sign(L) * stats::qnorm(pmax(P / 2, .Machine$double.xmin), lower.tail = FALSE)
    out[!is.finite(P)] <- NA_real_
    out[P >= 1] <- 0
    dimnames(out) <- dimnames(L)
    out
}

#' Validate that stored Lee statistics still match their normalized data and W.
#' @keywords internal
.lee_validate_current_provenance <- function(scope_obj,
                                             grid_name,
                                             lee_stats_layer,
                                             context,
                                             block_side = NULL) {
    rerun_hint <- if (grepl("clusterGenes", context, fixed = TRUE)) {
        "Rerun computeL() before clusterGenes()."
    } else if (grepl("computeLvsRCurve", context, fixed = TRUE)) {
        "Rerun computeL() before computeLvsRCurve()."
    } else {
        "Rerun computeL()."
    }
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    resolved_grid <- if (!is.null(grid_name)) {
        as.character(grid_name)[1L]
    } else {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)][1L]
    }
    if (is.na(resolved_grid) || !nzchar(resolved_grid)) {
        stop(context, ": unable to resolve the selected grid layer.", call. = FALSE)
    }

    stats_obj <- NULL
    if (!is.null(scope_obj@stats[[resolved_grid]])) {
        stats_obj <- scope_obj@stats[[resolved_grid]][[lee_stats_layer]]
    }
    if (is.null(stats_obj)) stats_obj <- g_layer[[lee_stats_layer]]
    meta <- if (is.list(stats_obj)) stats_obj$meta else NULL
    if (!is.list(meta) || !.lee_formula_id_is_supported(meta$formula_id)) {
        stop(
            context, ": LeeStats formula provenance is missing or incompatible. ",
            rerun_hint, " The canonical S2/W'W implementation is required.",
            call. = FALSE
        )
    }
    norm_layer <- meta$norm_layer
    if (!is.character(norm_layer) || length(norm_layer) != 1L ||
        is.na(norm_layer) || !nzchar(norm_layer)) {
        stop(context, ": LeeStats metadata has no valid norm_layer. ", rerun_hint, call. = FALSE)
    }
    if (is.null(g_layer[[norm_layer]])) {
        stop(
            context, ": the normalized layer '", norm_layer,
            "' used for Lee's L is no longer present; fallback to Xz is prohibited. ",
            rerun_hint,
            call. = FALSE
        )
    }
    X_current <- .validate_lee_norm_layer(g_layer[[norm_layer]], norm_layer, context = context)
    W_current <- g_layer$W
    grid_info <- g_layer$grid_info
    if (is.null(W_current)) {
        stop(context, ": spatial weights W are missing; rerun computeWeights() and computeL().", call. = FALSE)
    }
    if (is.null(attr(W_current, "weight_style", exact = TRUE))) {
        stored_style <- g_layer$weights_meta$weight_style
        if (!is.null(stored_style) && length(stored_style) == 1L &&
            !is.na(stored_style) && toupper(as.character(stored_style)) %in% c("B", "W")) {
            attr(W_current, "weight_style") <- toupper(as.character(stored_style))
        }
    }
    target_ids <- as.character(grid_info$grid_id)
    if (!length(target_ids) || anyNA(target_ids) || any(!nzchar(target_ids)) || anyDuplicated(target_ids)) {
        stop(context, ": grid_info$grid_id must be complete, non-empty, and unique.", call. = FALSE)
    }
    align_ids <- function(source_ids, label) {
        if (is.null(source_ids) || anyNA(source_ids) || anyDuplicated(source_ids)) {
            stop(context, ": ", label, " grid IDs must be complete and unique.", call. = FALSE)
        }
        idx <- match(target_ids, as.character(source_ids))
        if (anyNA(idx) || length(source_ids) != length(target_ids)) {
            stop(context, ": ", label, " grid IDs do not match grid_info$grid_id.", call. = FALSE)
        }
        idx
    }
    X_current <- X_current[align_ids(rownames(X_current), norm_layer), , drop = FALSE]
    current_weight_style <- .lee_weight_style(W_current)
    W_current <- W_current[
        align_ids(rownames(W_current), "W rows"),
        align_ids(colnames(W_current), "W columns"),
        drop = FALSE
    ]
    attr(W_current, "weight_style") <- current_weight_style
    rownames(X_current) <- target_ids
    dimnames(W_current) <- list(target_ids, target_ids)

    fingerprint_block_side <- block_side %||% meta$block_side
    fingerprint_block_side <- .lee_assert_positive_integer(
        fingerprint_block_side,
        "block_side"
    )
    current_fingerprint <- .lee_input_fingerprint(
        X_current,
        W_current,
        grid_info,
        norm_layer = norm_layer,
        block_side = fingerprint_block_side
    )
    if (!is.list(meta$input_fingerprint) ||
        !identical(meta$data_fingerprint, meta$input_fingerprint$data) ||
        !identical(meta$permutation_fingerprint, meta$input_fingerprint$permutation)) {
        stop(context, ": LeeStats fingerprint metadata are internally inconsistent. ", rerun_hint, call. = FALSE)
    }
    if (!.lee_same_observed_data(meta$input_fingerprint, current_fingerprint)) {
        stop(
            context, ": current X/W/grid data do not match the inputs that produced ",
            "stored Lee's L. ", rerun_hint,
            call. = FALSE
        )
    }
    if (!identical(meta$weight_style, current_weight_style)) {
        stop(context, ": current weight_style does not match LeeStats metadata; rerun computeL().", call. = FALSE)
    }
    current_S2 <- .lee_s2_value(W_current)
    if (!is.numeric(meta$S2) || length(meta$S2) != 1L || is.na(meta$S2) ||
        !is.finite(meta$S2) || meta$S2 <= 0 ||
        abs(current_S2 - meta$S2) > 1e-12 * max(1, abs(meta$S2))) {
        stop(context, ": current Lee S2 does not match LeeStats metadata; rerun computeL().", call. = FALSE)
    }

    list(
        grid_name = resolved_grid,
        grid_layer = g_layer,
        stats = stats_obj,
        meta = meta,
        norm_layer = norm_layer,
        X = X_current,
        W = W_current,
        grid_info = grid_info,
        input_fingerprint = current_fingerprint,
        S2 = current_S2,
        weight_style = current_weight_style
    )
}

#' @keywords internal
.compute_l <- function(scope_obj,
                     grid_name = NULL,
                     genes = NULL,
                     within = TRUE,
                     ncores = 1,
                     block_side = 8,
                     perms = 1000,
                     approximate_q = FALSE,
                     approximate_q_method = c("beta"),
                     approximate_q_threshold = 1000L,
                     block_size = 64,
                     L_min = 0,
                     norm_layer = "Xz",
                     lee_stats_layer_name = NULL,
                     legacy_formula = FALSE,
                     mem_limit_GB = 2,
                     chunk_size = 32L,
                     use_bigmemory = FALSE,
                     backing_path = tempdir(),
                     cache_inputs = TRUE,
                     verbose = TRUE,
                     ncore = NULL,
                     .user_set_bigmemory = NULL) {
    parent <- "computeL"
    user_set_bigmemory <- if (is.null(.user_set_bigmemory)) {
        !missing(use_bigmemory)
    } else {
        if (!is.logical(.user_set_bigmemory) || length(.user_set_bigmemory) != 1L || is.na(.user_set_bigmemory)) {
            stop(".user_set_bigmemory must be NULL or one non-missing logical value.", call. = FALSE)
        }
        isTRUE(.user_set_bigmemory)
    }
    approximate_q_method <- match.arg(approximate_q_method)
    approximate_q_threshold <- .lee_validate_integer_scalar(
        approximate_q_threshold,
        "approximate_q_threshold",
        lower = 1L
    )
    perms <- .lee_validate_integer_scalar(perms, "perms", lower = 0L)
    if (!is.logical(legacy_formula) || length(legacy_formula) != 1L || is.na(legacy_formula)) {
        stop("legacy_formula must be one non-missing logical value.", call. = FALSE)
    }
    if (legacy_formula) {
        stop(
            "legacy_formula=TRUE is no longer supported: computeL implements only the canonical Lee L denominator n/S2.",
            call. = FALSE
        )
    }
    within <- .lee_assert_flag(within, "within")
    if (perms > 0L && !isTRUE(within)) {
        stop(
            "Permutation inference for within=FALSE is not supported because it produces a rectangular gene-pair matrix; ",
            "use perms=0 for observed-L output or within=TRUE for inference.",
            call. = FALSE
        )
    }
    if (!is.logical(use_bigmemory) || length(use_bigmemory) != 1L || is.na(use_bigmemory)) {
        stop("use_bigmemory must be one non-missing logical value.", call. = FALSE)
    }
    if (use_bigmemory) {
        stop(
            "use_bigmemory=TRUE is disabled: the current route is RAM-backed and ",
            "does not provide a safe file-backed Lee L/FDR contract.",
            call. = FALSE
        )
    }

    ## --- 0. Thread-safe preprocessing: automatic thread management and error recovery ---

    # Handle deprecated parameter
    if (!is.null(ncore)) {
        warning("'ncore' is deprecated, please use 'ncores' instead. Will use its value this time.",
            call. = FALSE, immediate. = TRUE
        )
        ncores <- ncore
    }
    ncores <- .lee_assert_positive_integer(ncores, "ncores")
    block_side <- .lee_assert_positive_integer(block_side, "block_side")
    block_size <- .lee_assert_positive_integer(block_size, "block_size")
    mem_limit_GB <- .lee_assert_positive_finite(mem_limit_GB, "mem_limit_GB")
    chunk_size <- .lee_assert_positive_integer(chunk_size, "chunk_size")
    cache_inputs <- .lee_assert_flag(cache_inputs, "cache_inputs")
    verbose <- .lee_assert_flag(verbose, "verbose")
    norm_layer <- .lee_assert_layer_name(norm_layer, "norm_layer")
    if (!is.null(grid_name)) grid_name <- .lee_assert_layer_name(grid_name, "grid_name")
    if (!is.null(lee_stats_layer_name)) {
        lee_stats_layer_name <- .lee_assert_layer_name(lee_stats_layer_name, "lee_stats_layer_name")
    }
    if (!is.numeric(L_min) || length(L_min) != 1L || is.na(L_min) || !is.finite(L_min)) {
        stop("L_min must be a single finite number.", call. = FALSE)
    }
    if (!is.null(genes) && (!is.character(genes) || !length(genes) || anyNA(genes) ||
        any(!nzchar(genes)) || anyDuplicated(genes))) {
        stop("genes must be NULL or a non-empty vector of unique gene names.", call. = FALSE)
    }
    if (!is.character(backing_path) || length(backing_path) != 1L || is.na(backing_path) || !nzchar(backing_path)) {
        stop("backing_path must be a single non-empty character string.", call. = FALSE)
    }

    user_set_chunk <- !missing(chunk_size)

    # Get system information and aggressive thread count (use all visible cores up to request)
    os_type <- .detect_os()
    avail_cores <- suppressWarnings(detectCores(logical = TRUE))
    if (!is.numeric(avail_cores) || length(avail_cores) != 1L || is.na(avail_cores) ||
        !is.finite(avail_cores) || avail_cores < 1L) {
        avail_cores <- 1L
    }
    avail_cores <- as.integer(avail_cores)
    ncores <- min(ncores, avail_cores)
    grid_label <- if (is.null(grid_name)) "auto" else as.character(grid_name)[1]
    step01 <- .log_step(parent, "S01", "resolve inputs and weights", verbose)
    step01$enter(paste0("grid_name=", grid_label, " ncores=", ncores))

    # Use thread configuration for mixed OpenMP/BLAS operations
    thread_config <- .configure_threads_for("mixed", ncores, restore_after = TRUE)
    on.exit({
        restore_fn <- attr(thread_config, "restore_function")
        if (!is.null(restore_fn)) restore_fn()
    })

    # Use OpenMP threads for C++ operations
    ncores_cpp <- thread_config$openmp_threads
    .log_info(parent, "S01", paste0(
        "ncores_use=", ncores, "/", avail_cores,
        " openmp_threads=", ncores_cpp,
        " os=", os_type
    ), verbose)

    ## --- 1. Check if spatial weights exist, auto-compute if missing ---
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    grid_name <- if (is.null(grid_name)) {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    } else {
        grid_name
    }
    grid_info_n <- tryCatch(nrow(g_layer$grid_info), error = function(e) NA_integer_)
    if (!is.finite(grid_info_n) || grid_info_n < 2L) {
        stop("n_cells must be >= 2", call. = FALSE)
    }

    # Check if spatial weight matrix exists
    if (is.null(g_layer$W) || sum(g_layer$W) == 0) {
        warning(
            "No spatial weight matrix found; computing weights with style='B' (binary). ",
            "Binary weights affect Z-score variance computation. ",
            "Consider calling computeWeights() explicitly with your preferred style ",
            "before computeL() to ensure correct statistical properties.",
            call. = FALSE, immediate. = TRUE
        )
        .log_info(parent, "S01", "W missing; running computeWeights(style=B)", verbose)
        .log_backend(parent, "S01", "weights", "computeWeights(style=B)", verbose = verbose)
        scope_obj <- .compute_weights(scope_obj,
            grid_name = grid_name,
            style = "B",
            store_mat = TRUE,
            store_listw = FALSE,
            verbose = verbose,
            backend = "cpp"
        )
        # Re-extract layer after computing weights
        g_layer <- .select_grid_layer(scope_obj, grid_name)
    }

    step01$done(paste0("grid_name=", grid_name))

    ## --- 2. Memory guard and mode selection ---
    step02 <- .log_step(parent, "S02", "memory guard and mode selection", verbose)
    step02$enter(paste0(
        "mem_limit_GB=", mem_limit_GB,
        " use_bigmemory=", use_bigmemory,
        " chunk_size=", chunk_size
    ))
    # Cross-platform memory calculation
    all_genes <- if (!is.null(g_layer[[norm_layer]])) {
        ncol(g_layer[[norm_layer]])
    } else {
        stop("Normalized layer '", norm_layer, "' not found")
    }

    n_genes_use <- if (!is.null(genes) && isTRUE(within)) length(genes) else all_genes
    matrix_size_gb <- (as.double(n_genes_use)^2 * 8) / (1024^3)
    input_size_gb <- (as.double(nrow(g_layer[[norm_layer]])) * n_genes_use * 8) / (1024^3)

    # The observed kernel shares X/WX across workers. The permutation kernel
    # keeps one g x g accumulator per OpenMP worker to avoid data races.
    sys_mem_gb <- .get_system_memory_gb()
    observed_est_gb <- matrix_size_gb + 2 * input_size_gb
    permutation_est_gb <- if (perms > 0L) {
        (2 + ncores_cpp) * matrix_size_gb +
            (1 + 2 * ncores_cpp) * input_size_gb
    } else {
        0
    }
    est_total_gb <- max(observed_est_gb, permutation_est_gb)
    .log_info(parent, "S02", paste0(
        "matrix_size_gb=", round(matrix_size_gb, 1),
        " est_total_gb=", round(est_total_gb, 1),
        " sys_mem_gb=", round(sys_mem_gb, 1)
    ), verbose)
    if (est_total_gb > sys_mem_gb) {
        stop(
            "[geneSCOPE::.compute_l] Estimated memory requirement (",
            round(est_total_gb, 1), " GB) exceeds system capacity (",
            round(sys_mem_gb, 1), " GB). Reduce ncores or gene set size."
        )
    }

    mem_reason <- "within_limit"
    if (matrix_size_gb > mem_limit_GB) {
        stop(
            "Lee L output exceeds mem_limit_GB, but automatic bigmemory mode is disabled ",
            "because the current implementation is RAM-backed. Reduce the gene set or ",
            "raise mem_limit_GB only after confirming available RAM.",
            call. = FALSE
        )
    }

    mem_mode <- "inmemory"
    .log_backend(parent, "S02", "mem_mode", mem_mode, reason = mem_reason, verbose = verbose)

    # Cross-platform temporary directory
    if (use_bigmemory) {
        if (backing_path == tempdir()) {
            # Use platform-specific temporary directory
            backing_path <- switch(os_type,
                windows = Sys.getenv("TEMP", tempdir()),
                linux = Sys.getenv("TMPDIR", "/tmp"),
                macos = Sys.getenv("TMPDIR", "/tmp")
            )
        }

        if (!dir.exists(backing_path)) {
            dir.create(backing_path, recursive = TRUE)
        }
    }

    step02$done(paste0("mem_mode=", mem_mode))

    ## --- 3. Lee's L computation with error recovery ---
    step03 <- .log_step(parent, "S03", "compute Lee's L", verbose)
    step03$enter(paste0("mem_mode=", mem_mode, " ncores=", ncores))
    .log_backend(parent, "S03", "lee_L", paste0("C++ lee_L threads=", ncores_cpp), verbose = verbose)

    # Multi-retry mechanism with reduced thread count
    current_cores <- ncores
    min_cores <- 1
    success <- FALSE
    attempt <- 1
    last_error <- NULL

    while (!success && current_cores >= min_cores) {
        if (verbose && attempt > 1) {
            .log_info(parent, "S03", paste0("retry #", attempt, " with ", current_cores, " cores"), verbose)
            .log_backend(parent, "S03", "lee_L",
                paste0("C++ lee_L threads=", current_cores),
                reason = "retry",
                verbose = verbose
            )
        }

        result <- tryCatch(
            {
                # Execute computation
                t_start <- Sys.time()
                res <- .compute_lee_l(scope_obj,
                    grid_name = grid_name,
                    norm_layer = norm_layer,
                    genes = genes,
                    within = within,
                    ncores = current_cores, # Use current core count
                    mem_limit_GB = mem_limit_GB,
                    chunk_size = chunk_size,
                    use_bigmemory = use_bigmemory,
                    backing_path = backing_path,
                    block_side = block_side,
                    cache_inputs = cache_inputs,
                    input_cache = if (cache_inputs && !is.null(scope_obj@stats[[grid_name]])) {
                        scope_obj@stats[[grid_name]][["_lee_input_cache"]]
                    } else {
                        NULL
                    }
                )
                t_end <- Sys.time()

                if (verbose) {
                    time_msg <- if (attempt == 1) "completed" else "retry successful"
                    .log_info(parent, "S03", paste0(
                        "Lee's L ", time_msg, " (", format(t_end - t_start), ")"
                    ), verbose)
                }

                list(success = TRUE, object = res)
            },
            error = function(e) {
                if (verbose && attempt > 1) {
                    .log_info(parent, "S03", paste0("attempt failed: ", conditionMessage(e)), verbose)
                }
                list(success = FALSE, error = e)
            }
        )

        # Update status
        success <- result$success

        if (success) {
            res <- result$object
        } else {
            last_error <- conditionMessage(result$error)
            if (current_cores <= min_cores) {
                if (verbose) {
                    .log_info(parent, "S03", paste0("final attempt failed at ", current_cores, " core(s): ", last_error), verbose)
                }
                break
            }
            attempt <- attempt + 1
            # Reduce thread count
            current_cores <- max(floor(current_cores / 2), min_cores)
            if (verbose) {
                .log_info(parent, "S03", paste0("reducing cores to ", current_cores, " and retrying"), verbose)
                .log_backend(parent, "S03", "lee_L",
                    paste0("C++ lee_L threads=", current_cores),
                    reason = "retry",
                    verbose = verbose
                )
            }
            # Give system some time to recover
            Sys.sleep(3)
            # Force garbage collection
            gc(verbose = FALSE)
        }
    }

    if (!success) {
        if (!is.null(last_error) && nzchar(last_error)) {
            stop(last_error, call. = FALSE)
        }
        stop("Unable to compute Lee's L statistics even with minimal thread count", call. = FALSE)
    }

    step03$done(paste0("threads_final=", current_cores))

    L <- res$Lmat
    X_full <- res$X_full
    X_used <- res$X_used
    W <- res$W
    grid_inf <- res$grid_info
    gname <- res$grid_name
    n <- nrow(X_full)

    ## --- 4. Inference placeholder ---
    ## Moran's-I moments are not the moments of Lee's L.  Z is therefore
    ## derived from the exact two-sided permutation P matrix after Step 5.
    step04 <- .log_step(parent, "S04", "prepare empirical inference", verbose)
    step04$enter(paste0("within=", within, " mem_mode=", mem_mode))
    weight_style_info <- .resolve_leesl_weight_style(W, caller = parent, verbose = verbose)
    w_style <- weight_style_info$weight_style
    Z_mat <- matrix(NA_real_, nrow(L), ncol(L), dimnames = dimnames(L))
    step04$done("analytic Moran moments disabled")

    block_id <- res$block_id
    input_cache <- res$input_cache

    ## --- 5. Monte Carlo p-values with BLAS control and error recovery ---
    step05 <- .log_step(parent, "S05", "permutation p-values", verbose)
    step05$enter(paste0("perms=", perms, " block_size=", block_size))
    P <- if (perms > 0) {
        # Temporarily disable BLAS threads for permutation tests
        if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
            RhpcBLASctl::blas_set_num_threads(1)
        }

        if (inherits(L, "big.matrix")) {
            if (verbose) .log_info(parent, "S05", "converting to regular matrix for permutation tests", verbose)
            L_reg <- as.matrix(L)
        } else {
            L_reg <- L
        }

        # Permutation test with error recovery
        if (verbose) .log_info(parent, "S05", "running permutation tests", verbose)
        perm_success <- FALSE
        perm_cores <- current_cores
        perm_attempt <- 1
        perm_last_error <- NULL
        .log_backend(parent, "S05", "permutation",
            paste0("C++ lee_l_perm_block threads=", perm_cores),
            verbose = verbose
        )

        while (!perm_success && perm_cores >= 1) {
            if (verbose && perm_attempt > 1) {
                .log_info(parent, "S05", paste0("retry #", perm_attempt, " with ", perm_cores, " cores"), verbose)
                .log_backend(parent, "S05", "permutation",
                    paste0("C++ lee_l_perm_block threads=", perm_cores),
                    reason = "retry",
                    verbose = verbose
                )
            }

            # Sanity check to avoid C++ OOB when a gene subset is used
            if (ncol(X_used) != nrow(L_reg) || ncol(X_used) != ncol(L_reg)) {
                stop(sprintf("Permutation input mismatch: ncol(X_used)=%d but L_ref is %dx%d",
                    ncol(X_used), nrow(L_reg), ncol(L_reg)))
            }

            perm_result <- tryCatch(
                {
                    t_start <- Sys.time()
                    # Use X_used to match the dimensions of L_reg when a gene subset is used
                    p_result <- .lee_l_perm_block(X_used, W, L_reg,
                        block_id   = block_id,
                        perms      = perms,
                        block_size = block_size,
                        n_threads  = perm_cores
                    )
                    t_end <- Sys.time()
                    if (verbose) {
                        time_msg <- if (perm_attempt == 1) "completed" else "retry successful"
                        .log_info(parent, "S05", paste0(
                            "permutation test ", time_msg, " (", format(t_end - t_start), ")"
                        ), verbose)
                    }
                    list(success = TRUE, object = p_result)
                },
                error = function(e) {
                    if (verbose && perm_attempt > 1) {
                        .log_info(parent, "S05", paste0("test failed: ", conditionMessage(e)), verbose)
                    }
                    list(success = FALSE, error = e)
                }
            )

            if (perm_result$success) {
                perm_success <- TRUE
                P <- perm_result$object
            } else {
                perm_last_error <- conditionMessage(perm_result$error)
                if (perm_cores <= 1L) {
                    if (verbose) {
                        .log_info(
                            parent,
                            "S05",
                            paste0("final permutation attempt failed at ", perm_cores, " core(s): ", perm_last_error),
                            verbose
                        )
                    }
                    break
                }
                perm_attempt <- perm_attempt + 1
                perm_cores <- max(floor(perm_cores / 2), 1)
                if (verbose) {
                    .log_info(parent, "S05", paste0("reducing cores to ", perm_cores, " and retrying"), verbose)
                    .log_backend(parent, "S05", "permutation",
                        paste0("C++ lee_l_perm_block threads=", perm_cores),
                        reason = "retry",
                        verbose = verbose
                    )
                }
                Sys.sleep(2)
                gc(verbose = FALSE)
            }
        }

        if (!perm_success) {
            stop(
                "Lee's L permutation test failed",
                if (!is.null(perm_last_error) && nzchar(perm_last_error)) paste0(": ", perm_last_error) else ".",
                call. = FALSE
            )
        }
        P
    } else {
        .log_backend(parent, "S05", "permutation", "skipped", reason = "perms<=0", verbose = verbose)
        NULL
    }
    if (!is.null(P)) {
        dimnames(P) <- dimnames(L_reg)
        zero_var <- colSums(X_used^2) <= .Machine$double.eps
        if (any(zero_var)) {
            P[zero_var, ] <- NA_real_
            P[, zero_var] <- NA_real_
        }
        Z_mat <- .lee_empirical_signed_z(L_reg, P)
    } else {
        Z_mat <- NULL
    }
    step05$done(if (!is.null(P)) "empirical P and signed Z ready" else "no inference")

    use_approximate_q <- isTRUE(approximate_q) && !is.null(P) &&
        is.finite(perms) && perms > 0L && perms < approximate_q_threshold
    if (isTRUE(approximate_q) && is.null(P)) {
        warning("approximate_q=TRUE requested but perms=0; no inferential P/FDR/q matrices are computed.", call. = FALSE)
    }

    ## --- 6. FDR: Three Monte Carlo smoothing strategies as in .get_top_l_vs_r ---
    step06 <- .log_step(parent, "S06", "FDR and q-values", verbose)
    step06$enter(paste0("mem_mode=", mem_mode, " perms=", perms, " approximate_q=", use_approximate_q))
    # Explanation:
    #   If P exists: P = (k+1)/(perms+1); infer k = round(P*(perms+1)-1)
    #   p_beta    = (k+1)/(perms+2)        (Jeffreys, default more robust/slightly conservative)
    #   p_mid     = (k+0.5)/(perms+1)      (mid-p, less conservative)
    #   p_uniform = (k+U)/(perms+1)        (random jitter, expected equal to discrete grid, breaks staircase)
    #   Main output FDR = BH(p_beta)
    if (is.null(P)) {
        FDR_out_disc <- NULL
        FDR_out_beta <- NULL
        FDR_out_mid <- NULL
        FDR_out_uniform <- NULL
    } else if (inherits(Z_mat, "big.matrix")) {
        if (verbose) .log_info(parent, "S06", "computing FDR corrections", verbose)
        n_genes <- ncol(Z_mat)
        n_rows <- nrow(Z_mat)
        idx_chunks <- seq(1, n_genes, by = chunk_size)

        # Allocate big.matrices for final FDR output
        FDR_disc <- bigmemory::big.matrix(
            nrow = n_rows, ncol = n_genes,
            type = "double", init = NA_real_
        )
        FDR_beta <- bigmemory::big.matrix(
            nrow = n_rows, ncol = n_genes,
            type = "double", init = NA_real_
        )
        FDR_mid <- bigmemory::big.matrix(
            nrow = n_rows, ncol = n_genes,
            type = "double", init = NA_real_
        )
        FDR_uniform <- bigmemory::big.matrix(
            nrow = n_rows, ncol = n_genes,
            type = "double", init = NA_real_
        )

        # PASS 1: Collect all p-values (global BH requires all p-values simultaneously)
        all_p_disc <- numeric(0)
        all_p_beta <- numeric(0)
        all_p_mid <- numeric(0)
        all_p_uniform <- numeric(0)
        chunk_info <- list()

        for (start in idx_chunks) {
            end <- min(start + chunk_size - 1L, n_genes)
            chunk_n <- (end - start + 1L) * n_rows
            chunk_info[[length(chunk_info) + 1L]] <- list(start = start, end = end, n = chunk_n)

            if (!is.null(P)) {
                p_disc <- P[, start:end]
                k_est <- round(p_disc * (perms + 1) - 1)
                k_est[k_est < 0] <- 0
                k_est[k_est > perms] <- perms
                p_beta <- (k_est + 1) / (perms + 2)
                p_mid <- (k_est + 0.5) / (perms + 1)
                p_uniform <- (k_est + matrix(runif(length(k_est)), nrow = nrow(k_est))) / (perms + 1)
            } else {
                Z_chunk <- Z_mat[, start:end]
                p_disc <- 2 * pnorm(-abs(Z_chunk))
                p_beta <- p_mid <- p_uniform <- p_disc
            }

            all_p_disc <- c(all_p_disc, as.numeric(p_disc))
            all_p_beta <- c(all_p_beta, as.numeric(p_beta))
            all_p_mid <- c(all_p_mid, as.numeric(p_mid))
            all_p_uniform <- c(all_p_uniform, as.numeric(p_uniform))
        }

        # Apply global BH correction to all p-values at once (correct implementation)
        all_q_disc <- p.adjust(all_p_disc, method = "BH")
        all_q_beta <- p.adjust(all_p_beta, method = "BH")
        all_q_mid <- p.adjust(all_p_mid, method = "BH")
        all_q_uniform <- p.adjust(all_p_uniform, method = "BH")

        # PASS 2: Redistribute q-values to big.matrix chunks
        idx_offset <- 1L
        for (ci in chunk_info) {
            start <- ci$start
            end <- ci$end
            n_chunk <- ci$n

            q_disc_chunk <- all_q_disc[idx_offset:(idx_offset + n_chunk - 1L)]
            q_beta_chunk <- all_q_beta[idx_offset:(idx_offset + n_chunk - 1L)]
            q_mid_chunk <- all_q_mid[idx_offset:(idx_offset + n_chunk - 1L)]
            q_uniform_chunk <- all_q_uniform[idx_offset:(idx_offset + n_chunk - 1L)]

            FDR_disc[, start:end] <- matrix(q_disc_chunk, nrow = n_rows)
            FDR_beta[, start:end] <- matrix(q_beta_chunk, nrow = n_rows)
            FDR_mid[, start:end] <- matrix(q_mid_chunk, nrow = n_rows)
            FDR_uniform[, start:end] <- matrix(q_uniform_chunk, nrow = n_rows)

            idx_offset <- idx_offset + n_chunk
        }

        FDR_out_disc <- FDR_disc
        FDR_out_beta <- FDR_beta
        FDR_out_mid <- FDR_mid
        FDR_out_uniform <- FDR_uniform
    } else {
        # In-memory mode: compute in one go (exact BH)
        p_disc <- P
        k_est <- round(p_disc * (perms + 1) - 1)
        k_est[k_est < 0] <- 0
        k_est[k_est > perms] <- perms
        p_beta <- (k_est + 1) / (perms + 2)
        p_mid <- (k_est + 0.5) / (perms + 1)
        p_uniform <- (k_est + runif(length(k_est))) / (perms + 1)
        # Adjust the intended family: unique unordered, off-diagonal gene pairs.
        FDR_out_disc <- .lee_pairwise_p_adjust(p_disc, "BH")
        FDR_out_beta <- .lee_pairwise_p_adjust(p_beta, "BH")
        FDR_out_mid <- .lee_pairwise_p_adjust(p_mid, "BH")
        FDR_out_uniform <- .lee_pairwise_p_adjust(p_uniform, "BH")
    }

    ## --- 5.x Extra: Storey q-value main FDR (in-memory mode only, big.matrix fallback) ---
    FDR_storey <- NULL
    FDR_approx <- NULL
    approximate_q_report <- NULL
    pi0_hat <- NA_real_
    FDR_main_method <- "not computed"
    if (!inherits(Z_mat, "big.matrix") && !is.null(P)) {
        storey <- .lee_pairwise_storey(p_beta)
        FDR_storey <- storey$q_values
        pi0_hat <- storey$pi0_hat
    }
    if (use_approximate_q) {
        approximate_q_report <- approximate_q_values(
            p_values = P,
            perms = perms,
            method = approximate_q_method
        )
        FDR_approx <- approximate_q_report$q_values
        if (!identical(dim(FDR_approx), dim(P))) {
            FDR_approx <- matrix(FDR_approx, nrow = nrow(P), dimnames = dimnames(P))
        }
        FDR_main_method <- paste0("approximate Storey q-value (", approximate_q_method, " posterior-smoothed permutation p)")
        pi0_hat <- approximate_q_report$pi0_hat
        warning(
            sprintf(
                "computeL: approximate_q=TRUE and perms=%d < %d; using %s q-value approximation as a screening aid only, not confirmatory inference. Increase permutations for final inference.",
                as.integer(perms), as.integer(approximate_q_threshold), approximate_q_method
            ),
            call. = FALSE
        )
        if (!is.null(approximate_q_report$warning) && nzchar(approximate_q_report$warning)) {
            warning(
                "computeL: approximate_q diagnostic: ", approximate_q_report$warning,
                ". Treat this as a boundary condition and increase permutations where feasible.",
                call. = FALSE
            )
        }
        .log_key(
            parent,
            sprintf(
                "WARN: \033[31mapproximate_q active: perms=%d < threshold=%d, method=%s. Screening only; not confirmatory.\033[39m",
                as.integer(perms), as.integer(approximate_q_threshold), approximate_q_method
            )
        )
    }
    # The confirmatory output uses exact Monte-Carlo P and BH over unique pairs.
    # The explicitly requested low-budget approximation remains opt-in only.
    FDR_main <- if (!is.null(FDR_approx)) FDR_approx else FDR_out_disc
    if (is.null(P)) {
        FDR_main_method <- "unavailable (perms=0)"
    } else if (is.null(FDR_approx)) {
        FDR_main_method <- "BH(exact empirical P; unique unordered off-diagonal pairs)"
    }
    fdr_main_finite <- if (is.null(FDR_main)) numeric(0) else if (nrow(FDR_main) == ncol(FDR_main)) {
        FDR_main[upper.tri(FDR_main, diag = FALSE)]
    } else as.numeric(FDR_main)
    fdr_main_finite <- fdr_main_finite[is.finite(fdr_main_finite)]
    n_sig_005 <- if (is.null(P)) NA_integer_ else sum(fdr_main_finite < 0.05)
    p_test_values <- if (is.null(P)) numeric(0) else if (nrow(P) == ncol(P)) {
        P[upper.tri(P, diag = FALSE)]
    } else as.numeric(P)
    n_tests <- sum(is.finite(p_test_values))
    min_p_possible <- if (!is.null(P)) 1 / (perms + 1) else NA_real_

    step06$done(paste0("FDR_main_method=", FDR_main_method))

    ## --- 7. gradients and QC metrics ---
    step07 <- .log_step(parent, "S07", "gradients and QC metrics", verbose)
    step07$enter(paste0("within=", within, " L_min=", L_min))
    .log_backend(parent, "S07", "qc_graph", "igraph louvain (fallback=leiden)", verbose = verbose)

    ## --- 7. βx / βy ---
    centres <- with(
        grid_inf[match(res$cells, grid_inf$grid_id), ],
        data.frame(
            x = (xmin + xmax) / 2,
            y = (ymin + ymax) / 2
        )
    )
    betas <- t(apply(X_full, 2, function(v) {
        c(beta_x = coef(lm(v ~ centres$x + centres$y))[2:3])
    }))
    if (!is.null(genes)) betas <- betas[genes, , drop = FALSE]

    ## --- 7. QC computation ---
    qc <- NULL
    if (within || is.null(genes)) {
        # Convert to regular matrix for QC if needed
        L_qc <- if (inherits(L, "big.matrix")) {
            .log_info(parent, "S07", "converting subset of big.matrix for QC computation", verbose)
            # Sample subset for QC to avoid memory issues
            n_sample <- min(1000, nrow(L))
            sample_idx <- sample(nrow(L), n_sample)
            as.matrix(L[sample_idx, sample_idx])
        } else {
            L
        }

        A_bin <- abs(L_qc) >= L_min
        A_bin <- A_bin | t(A_bin)
        diag(A_bin) <- FALSE

        A_num <- (A_bin | t(A_bin)) * 1
        g_tmp <- igraph::simplify(
            igraph::graph_from_adjacency_matrix(A_num,
                mode = "max",
                diag = FALSE
            ),
            remove.multiple = TRUE,
            remove.loops = TRUE
        )

        deg <- igraph::degree(g_tmp)
        memb <- tryCatch(
            igraph::cluster_louvain(g_tmp)$membership,
            error = function(e) {
                igraph::cluster_leiden(g_tmp)$membership
            }
        )

        qc <- list(
            edge_density = 2 * sum(A_bin[upper.tri(A_bin)]) /
                (ncol(A_bin) * (ncol(A_bin) - 1)),
            components = igraph::components(g_tmp)$no,
            modularity_Q = igraph::modularity(g_tmp, membership = memb),
            mean_degree = mean(deg),
            sd_degree = sd(deg),
            hub_ratio = mean(deg > 2 * median(deg)),
            sig_edge_frac = if (is.null(P)) {
                NA
            } else {
                if (inherits(P, "big.matrix")) {
                    # Sample for significance fraction
                    p_sample <- as.matrix(P[sample_idx, sample_idx])
                    mean(p_sample[upper.tri(p_sample)] < 0.05, na.rm = TRUE)
                } else {
                    mean(P[upper.tri(P)] < 0.05, na.rm = TRUE)
                }
            }
        )
    }

    step07$done()

    ## --- 8. Write to stats (main FDR replacement; retain old matrices; add FDR_storey) ---
    step08 <- .log_step(parent, "S08", "store outputs", verbose)
    layer_name <- if (is.null(lee_stats_layer_name)) paste0("LeeStats_", norm_layer) else lee_stats_layer_name
    step08$enter(paste0("layer_name=", layer_name))
    if (is.null(scope_obj@stats)) scope_obj@stats <- list()
    if (is.null(scope_obj@stats[[gname]])) scope_obj@stats[[gname]] <- list()

    scope_obj@stats[[gname]][[layer_name]] <- list(
        L = L,
        Z = Z_mat,
        P = P,
        grad = betas,
        L_min = L_min,
        qc = qc,
        FDR = FDR_main, # Exact empirical-P BH, unless approximate_q was explicitly requested
        FDR_storey = FDR_storey, # Exploratory unique-pair diagnostic
        FDR_approx = FDR_approx,
        FDR_disc = FDR_out_disc, # Exact empirical-P BH over unique unordered pairs
        FDR_beta = FDR_out_beta,
        FDR_mid = FDR_out_mid,
        FDR_uniform = FDR_out_uniform,
        weight_style = w_style,
        meta = list(
            formula_id = .lee_formula_id(),
            norm_layer = norm_layer,
            input_fingerprint = res$input_fingerprint,
            data_fingerprint = res$input_fingerprint$data,
            permutation_fingerprint = res$input_fingerprint$permutation,
            weight_style = w_style,
            weight_style_source = weight_style_info$source,
            weight_style_inferred = weight_style_info$inferred,
            S2 = res$S2,
            perms = perms,
            block_side = block_side,
            block_size = block_size,
            ncores = current_cores,
            ncores_requested = ncores,
            permutation_ncores = if (perms > 0L) perm_cores else NA_integer_,
            mem_mode = if (inherits(L, "big.matrix")) "bigmemory" else "inmemory",
            p_source = if (!is.null(P)) "permutation" else "unavailable",
            z_source = if (!is.null(P)) "signed inverse-normal transform of exact two-sided permutation P" else "unavailable",
            p_resolution = if (!is.null(P)) paste0("~", format(1 / (perms + 1), scientific = TRUE)) else "unavailable",
            pi0_hat = pi0_hat,
            n_tests = n_tests,
            n_sig_FDR_lt_0.05 = n_sig_005,
            min_p_possible = min_p_possible,
            FDR_main_method = FDR_main_method,
            approximate_q = isTRUE(use_approximate_q),
            approximate_q_status = if (isTRUE(use_approximate_q)) "active_screening_only" else if (isTRUE(approximate_q)) "requested_but_not_active" else "disabled",
            approximate_q_method = if (isTRUE(use_approximate_q)) approximate_q_method else NA_character_,
            approximate_q_threshold = approximate_q_threshold,
            approximate_q_pi0_hat = if (!is.null(approximate_q_report)) approximate_q_report$pi0_hat else NA_real_,
            approximate_q_warning = if (!is.null(approximate_q_report)) approximate_q_report$warning else NA_character_,
            approximate_q_fallback = if (isTRUE(approximate_q) && !isTRUE(use_approximate_q)) {
                if (is.null(P)) "permutation_p_matrix_unavailable" else "perms_not_below_threshold"
            } else NA_character_,
            approximate_q_note = if (isTRUE(use_approximate_q)) "Low-permutation q-values were estimated from beta posterior-smoothed permutation p-values; use for screening, not final confirmatory inference." else NA_character_,
            pval_modes = c("exact=(k+1)/(B+1)", "beta=(k+1)/(B+2)", "mid=(k+0.5)/(B+1)", "uniform=(k+U)/(B+1)"),
            note = "FDR main is BH on exact empirical P over unique unordered off-diagonal pairs unless approximate_q is explicitly enabled; beta/mid/uniform/Storey matrices are exploratory diagnostics."
        )
    )

    if (cache_inputs) {
        if (is.null(scope_obj@stats[[gname]])) scope_obj@stats[[gname]] <- list()
        scope_obj@stats[[gname]][["_lee_input_cache"]] <- input_cache
    }

    step08$done(paste0("layer_name=", layer_name))

    invisible(scope_obj)
}

#' Bootstrap Lee's L vs Pearson r relationship
#' @description
#' Internal helper for `.compute_l_vs_r_curve`.
#' Computes a LOESS-based relationship between Lee's L and Pearson correlation
#' with bootstrap confidence intervals.
#' @param scope_obj A `scope_object` containing Lee's L and correlation matrices.
#' @param grid_name Grid layer name to operate on.
#' @param level Correlation level (`cell` or `grid`).
#' @param lee_stats_layer Lee statistics layer name.
#' @param span LOESS span parameter.
#' @param B Number of bootstrap iterations.
#' @param deg Degree for the LOESS fit.
#' @param ncores Number of threads to use.
#' @param length_out Grid size for the fitted curve.
#' @param downsample Optional downsampling rate for points.
#' @param n_strata Number of strata used when drawing bootstrap samples.
#' @param jitter_eps Optional jitter applied to correlation values.
#' @param ci_method Confidence interval method to use.
#' @param ci_adjust Optional analytic adjustment for CI width.
#' @param verbose Emit progress messages when TRUE.
#' @param k_max Numeric threshold.
#' @param min_rel_width Numeric threshold.
#' @param widen_span Parameter value.
#' @param curve_name Parameter value.
#' @return The modified `scope_object` with stored curve data.
#' @keywords internal
.compute_l_vs_r_curve <- function(scope_obj,
                        grid_name,
                        level = c("grid", "cell"),
                        lee_stats_layer = "LeeStats_Xz",
                        span = 0.45,
                        B = 1000,
                        deg = 1,
                        ncores = max(1, detectCores() - 1),
                        length_out = 1000,
                        downsample = 1,
                        n_strata = 50,
                        k_max = Inf,
                        jitter_eps = 0,
                        ci_method = c("percentile", "basic", "bc"),
                        ci_adjust = c("none", "analytic"),
                        min_rel_width = 0,
                        widen_span = 0.1,
                        curve_name = "LR_curve2",
                        backend = "cpp",
                        verbose = TRUE) {
  ci_method <- match.arg(ci_method)
  ci_adjust <- match.arg(ci_adjust)
  level <- match.arg(level)
  backend <- if (missing(backend)) {
    "cpp"
  } else {
    .normalize_core_backend(backend, arg = "backend", allow_auto = TRUE)
  }
  backend_policy <- .core_backend_policy(
    backend,
    cpp_supported = FALSE,
    python_supported = FALSE,
    r_supported = TRUE
  )
  if (identical(backend, "python")) {
    stop(
      "[geneSCOPE::computeLvsRCurve] backend='", backend,
      "' is not implemented for this R-only curve-fitting layer. ",
      "Use the default C++->R fallback policy or backend='r'.",
      call. = FALSE
    )
  } else if (backend %in% c("auto", "cpp") && verbose) {
    .log_info("computeLvsRCurve", "S01", "Curve fitting is R-only; using R fallback under C++->R policy.", verbose)
  }
  ncores <- max(1L, min(as.integer(ncores), detectCores(logical = TRUE)))
  if (min_rel_width < 0) stop("min_rel_width cannot be negative")
  parent <- "computeLvsRCurve"
  if (B < 20 && verbose) {
    .log_info(parent, "S01", "Warning: B < 20 may be unstable", verbose)
  }

  step01 <- .log_step(parent, "S01", "extract L and correlation", verbose)
  step01$enter(paste0("grid_name=", grid_name, " level=", level))
  .log_info(parent, "S01", paste0(
    "bootstrap_iterations=", B,
    " ci_method=", ci_method,
    " ncores=", ncores
  ), verbose)

  # 1. Extract matrices
  if (verbose) .log_info(parent, "S01", "extracting Lee's L and Pearson correlation matrices", verbose)
  provenance <- .lee_validate_current_provenance(
    scope_obj,
    grid_name = grid_name,
    lee_stats_layer = lee_stats_layer,
    context = "[geneSCOPE::computeLvsRCurve]"
  )
  g_layer <- provenance$grid_layer
  grid_name <- provenance$grid_name
  Lmat <- .get_lee_matrix(scope_obj, grid_name, lee_layer = lee_stats_layer)
  rmat <- .get_pearson_matrix(scope_obj, grid_name, level = ifelse(level == "grid", "grid", "cell"))

  common <- intersect(rownames(Lmat), rownames(rmat))
  if (length(common) < 2) stop("Insufficient common genes")

  if (verbose) {
    .log_info(parent, "S01", paste0("common_genes=", length(common)), verbose)
    .log_info(parent, "S01", paste0("matrix_dim=", nrow(Lmat), "x", ncol(Lmat)), verbose)
  }

  step01$done(paste0("common_genes=", length(common)))

  step02 <- .log_step(parent, "S02", "prepare data vectors", verbose)
  step02$enter(paste0("downsample=", downsample, " jitter_eps=", jitter_eps))

  # Memory guard: estimate per-thread footprint (two dense matrices) and stop if over system RAM
  sys_mem_gb <- .get_system_memory_gb()
  per_thread_gb <- (length(common)^2 * 8 * 2) / (1024^3)
  est_total_gb <- per_thread_gb * ncores
  if (verbose) {
    .log_info(parent, "S02", paste0(
      "est_total_gb=", round(est_total_gb, 1),
      " sys_mem_gb=", round(sys_mem_gb, 1)
    ), verbose)
  }
  if (est_total_gb > sys_mem_gb) {
    stop(
      "[geneSCOPE::.compute_l_vs_r_curve] Estimated memory requirement (",
      round(est_total_gb, 1), " GB) exceeds system capacity (",
      round(sys_mem_gb, 1), " GB). Reduce ncores, downsample, or gene set size."
    )
  }

  Lmat <- Lmat[common, common, drop = FALSE]
  rmat <- rmat[common, common, drop = FALSE]
  ut <- upper.tri(Lmat, diag = FALSE)
  Lv <- Lmat[ut]
  rv <- rmat[ut]

  # 2. Clean / Downsample
  if (verbose) .log_info(parent, "S02", "cleaning and preprocessing data points", verbose)
  ok <- is.finite(Lv) & is.finite(rv)
  Lv <- Lv[ok]
  rv <- rv[ok]
  if (!length(Lv)) stop("No valid points")

  if (verbose) .log_info(parent, "S02", paste0("initial_points=", length(Lv)), verbose)

  if (is.numeric(downsample) && downsample < 1) {
    keep <- sample.int(length(Lv), max(1L, floor(downsample * length(Lv))))
    Lv <- Lv[keep]
    rv <- rv[keep]
    if (verbose) .log_info(parent, "S02", paste0(
      "downsampled_points=", length(Lv), " ratio=", downsample
    ), verbose)
  } else if (is.numeric(downsample) && downsample >= 1 && length(Lv) > downsample) {
    keep <- sample.int(length(Lv), downsample)
    Lv <- Lv[keep]
    rv <- rv[keep]
    if (verbose) .log_info(parent, "S02", paste0(
      "downsampled_points=", length(Lv), " target=", downsample
    ), verbose)
  }
  if (jitter_eps > 0) {
    rv <- jitter(rv, factor = jitter_eps)
    if (verbose) .log_info(parent, "S02", paste0("jitter_eps=", jitter_eps), verbose)
  }

  step02$done(paste0("points=", length(Lv)))

  # 3. Stratify (by r quantiles) -- fallback if failed
  step03 <- .log_step(parent, "S03", "stratify and build grid", verbose)
  step03$enter(paste0("n_strata=", n_strata, " length_out=", length_out))
  if (verbose) .log_info(parent, "S03", "setting up stratified sampling", verbose)
  uniq_r <- sort(unique(rv))
  if (length(uniq_r) < 3) stop("r values too discrete")
  n_strata_eff <- min(n_strata, length(uniq_r) - 1)

  if (verbose) {
    .log_info(parent, "S03", paste0("unique_r=", length(uniq_r)), verbose)
    .log_info(parent, "S03", paste0("effective_strata=", n_strata_eff), verbose)
  }

  probs <- seq(0, 1, length.out = n_strata_eff + 1)
  brks <- unique(quantile(rv, probs, na.rm = TRUE))
  if (length(brks) < 2) {
    # Fallback: uniform slices
    brks <- seq(min(rv), max(rv), length.out = n_strata_eff + 1)
    brks <- unique(brks)
    if (verbose) .log_info(parent, "S03", "using fallback uniform stratification", verbose)
  }
  if (length(brks) < 2) stop("Cannot establish strata")
  strat <- cut(rv, breaks = brks, include.lowest = TRUE, labels = FALSE)
  ok2 <- !is.na(strat)
  rv <- rv[ok2]
  Lv <- Lv[ok2]
  strat <- strat[ok2]

  # 4. xgrid
  if (verbose) .log_info(parent, "S03", "preparing analysis grid and fitting LOESS model", verbose)
  xr <- range(rv)
  if (diff(xr) <= 0) stop("r has no span")
  xgrid <- seq(xr[1], xr[2], length.out = length_out)

  if (verbose) {
    .log_info(parent, "S03", paste0(
      "analysis_range=[", round(xr[1], 3), ", ", round(xr[2], 3), "]"
    ), verbose)
    .log_info(parent, "S03", paste0("grid_points=", length_out), verbose)
    .log_info(parent, "S03", paste0("bootstrap_cores=", ncores), verbose)
  }

  step03$done(paste0("xgrid_points=", length_out))

  # Directly call unified interface (no longer check old version)
  step04 <- .log_step(parent, "S04", "bootstrap LOESS curve", verbose)
  step04$enter(paste0("B=", B, " span=", span, " deg=", deg))
  .log_backend(parent, "S04", "loess_bootstrap",
    paste0("C++ loess_residual_bootstrap threads=", ncores),
    verbose = verbose
  )
  if (verbose) .log_info(parent, "S04", "running LOESS residual bootstrap analysis", verbose)
  keep_boot <- TRUE
  adjust_mode <- if (ci_method == "percentile" && ci_adjust == "analytic") 1L else 0L
  res <- .loess_residual_bootstrap(
    x = rv, y = Lv, strat = as.integer(strat),
    grid = xgrid,
    B = as.integer(B),
    span = span,
    deg = as.integer(deg),
    n_threads = as.integer(max(1, ncores)),
    k_max = if (is.finite(k_max)) as.integer(k_max) else -1L,
    keep_boot = keep_boot,
    adjust_mode = adjust_mode,
    ci_type = 0L,
    level = 0.95
  )

  fit <- res$fit
  lo <- res$lo
  hi <- res$hi
  if (ci_method == "basic") {
    lo <- res$lo_basic
    hi <- res$hi_basic
  } else if (ci_method == "bc") {
    lo <- res$lo_bc
    hi <- res$hi_bc
  }

  step04$done(paste0("B=", res$B))

  step05 <- .log_step(parent, "S05", "adjust confidence intervals", verbose)
  step05$enter(paste0("ci_method=", ci_method, " ci_adjust=", ci_adjust))
  if (verbose) {
    .log_info(parent, "S05", "bootstrap completed, processing confidence intervals", verbose)
    .log_info(parent, "S05", paste0("ci_method_applied=", ci_method), verbose)
  }

  # 6. floor widen
  if (verbose && min_rel_width > 0) {
    .log_info(parent, "S05", "applying minimum relative width constraints", verbose)
  }
  if (min_rel_width > 0) {
    rng_fit <- diff(range(fit, finite = TRUE))
    if (!is.finite(rng_fit) || rng_fit <= 0) rng_fit <- 1
    width <- hi - lo
    target_w <- min_rel_width * rng_fit
    width <- pmax(width, target_w)
    # Smooth width (LOESS)
    if (is.finite(widen_span) && widen_span > 0 && length(width) > 10) {
      sm <- tryCatch(loess(width ~ xgrid, span = widen_span), error = function(e) NULL)
      if (!is.null(sm)) {
        wp <- tryCatch(predict(sm, xgrid), error = function(e) NULL)
        if (!is.null(wp) && all(is.finite(wp))) {
          width <- pmax(wp, target_w)
        }
      }
    }
    center <- (lo + hi) / 2
    lo <- center - width / 2
    hi <- center + width / 2
  }

  # ---- New: Local residual scale adaptive lower bound (no new parameters) ----------------------------
  # Purpose: Avoid "points more dispersed but CI narrower"; force width >= 2*1.96*local MAD
  if (verbose) .log_info(parent, "S05", "applying local residual scale adaptation", verbose)
  {
    if (length(Lv) > 20) {
      # Fit values interpolated to original r
      fit_at_rv <- tryCatch(approx(xgrid, fit, xout = rv, rule = 2, ties = "ordered")$y,
        error = function(e) rep(mean(fit), length(rv))
      )
      res_raw <- Lv - fit_at_rv
      ord_r <- order(rv)
      rv_sorted <- rv[ord_r]
      res_sorted <- res_raw[ord_r]
      n_pts <- length(rv_sorted)
      # Fixed K (internal constant, no new parameters)
      K <- min(200L, n_pts)
      if (K > 5) {
        # Nearest neighbor index function (bidirectional expansion)
        get_knn_idx <- function(x0) {
          pos <- findInterval(x0, rv_sorted)
          l <- pos
          r <- pos + 1
          out <- integer(0)
          while (length(out) < K && (l >= 1 || r <= n_pts)) {
            dl <- if (l >= 1) abs(rv_sorted[l] - x0) else Inf
            dr <- if (r <= n_pts) abs(rv_sorted[r] - x0) else Inf
            if (dl <= dr) {
              if (l >= 1) {
                out <- c(out, l)
                l <- l - 1
              }
            } else {
              if (r <= n_pts) {
                out <- c(out, r)
                r <- r + 1
              }
            }
          }
          out
        }
        # Calculate local MAD
        loc_mad <- vapply(xgrid, function(xx) {
          idx <- get_knn_idx(xx)
          if (!length(idx)) {
            return(NA_real_)
          }
          rr <- res_sorted[idx]
          med <- median(rr)
          1.4826 * median(abs(rr - med))
        }, numeric(1))
        # Global benchmark (prevent all NA)
        med_mad <- median(loc_mad[is.finite(loc_mad) & loc_mad > 0], na.rm = TRUE)
        if (is.finite(med_mad) && med_mad > 0) {
          width <- hi - lo
          # Required minimum width
          w_need <- 2 * 1.96 * loc_mad
          # If local MAD is NA or 0, no restriction
          bad <- !is.finite(w_need) | w_need <= 0
          if (any(!bad)) {
            width_new <- width
            width_new[!bad] <- pmax(width[!bad], w_need[!bad])
            if (!identical(width_new, width)) {
              center <- (lo + hi) / 2
              lo <- center - width_new / 2
              hi <- center + width_new / 2
              # Record diagnostics (add to meta)
              local_mad_diag <- list(
                K = K,
                mad_global_med = med_mad,
                frac_expanded = mean(width_new > width)
              )
            }
          }
        }
      }
    }
  }

  # ---- New: CI edge smoothing (post-processing only, does not change algorithm) ---------------------------
  if (verbose) .log_info(parent, "S05", "applying confidence interval edge smoothing", verbose)
  {
    if (length(lo) >= 15 && all(is.finite(xgrid))) {
      span_s <- 0.06 + 10 / length(lo) # Adaptive small span, smaller with longer length
      lo_s <- tryCatch(
        predict(loess(lo ~ xgrid,
          span = span_s, degree = 1,
          surface = "direct", family = "gaussian"
        )),
        error = function(e) lo
      )
      hi_s <- tryCatch(
        predict(loess(hi ~ xgrid,
          span = span_s, degree = 1,
          surface = "direct", family = "gaussian"
        )),
        error = function(e) hi
      )
      # Use smoothed width as main, keep original center to avoid overall drift
      center_old <- (lo + hi) / 2
      width_new <- pmax(hi_s - lo_s, 0)
      lo_new <- center_old - width_new / 2
      hi_new <- center_old + width_new / 2
      # Prevent numerical error causing reversal
      swap_idx <- which(hi_new < lo_new)
      if (length(swap_idx)) {
        tmp <- lo_new[swap_idx]
        lo_new[swap_idx] <- hi_new[swap_idx]
        hi_new[swap_idx] <- tmp
      }
      # If NA produced, fallback
      if (all(is.finite(lo_new)) && all(is.finite(hi_new))) {
        edge_smooth_info <- list(
          edge_smooth = TRUE,
          edge_smooth_span = span_s,
          edge_smooth_expanded = mean(width_new > (hi - lo))
        )
        lo <- lo_new
        hi <- hi_new
      } else {
        edge_smooth_info <- list(edge_smooth = FALSE)
      }
    } else {
      edge_smooth_info <- list(edge_smooth = FALSE)
    }
  }

  # ---- CI edge smoothing completed, enter length check ----
  step05$done()

  step06 <- .log_step(parent, "S06", "store curve", verbose)
  step06$enter(paste0("curve_name=", curve_name))
  if (verbose) .log_info(parent, "S06", "finalizing curve data and storing results", verbose)
  ## ---- 6. Write back to @stats -------------------------------------------------------
  n_grid <- length(xgrid)
  if (!all(
    length(fit) == n_grid,
    length(lo) == n_grid,
    length(hi) == n_grid
  )) {
    stop(sprintf(
      "Internal length mismatch: fit=%d lo=%d hi=%d xgrid=%d",
      length(fit), length(lo), length(hi), n_grid
    ))
  }

  if (is.null(scope_obj@stats)) scope_obj@stats <- list()
  if (is.null(scope_obj@stats[[grid_name]])) {
    scope_obj@stats[[grid_name]] <- list()
  }
  if (is.null(scope_obj@stats[[grid_name]][[lee_stats_layer]])) {
    scope_obj@stats[[grid_name]][[lee_stats_layer]] <- list()
  }

  ## ---- meta update ----
  meta_obj <- list(
    B = res$B,
    span = span,
    deg = deg,
    backend_requested = backend,
    backend_selected = "r",
    backend_policy = backend_policy,
    ci_method = ci_method,
    ci_adjust = ci_adjust,
    min_rel_width = min_rel_width,
    n_strata = n_strata,
    k_max = k_max,
    edf = res$edf,
    sigma2_raw = res$sigma2_raw,
    sigma2_edf = res$sigma2_edf,
    resid_global_mad = res$resid_mad,
    adjust_mode = adjust_mode,
    note = "Generated by .compute_l_vs_r_curve with residual-MAD floor + edge smoothing",
    local_mad_diag = if (exists("local_mad_diag")) local_mad_diag else NULL,
    edge_smooth = if (exists("edge_smooth_info")) edge_smooth_info else NULL
  )
  meta_col <- rep(list(meta_obj), length(fit))

  scope_obj@stats[[grid_name]][[lee_stats_layer]][[curve_name]] <-
    data.frame(Pear = xgrid, fit = fit, lo95 = lo, hi95 = hi, meta = I(meta_col))

  if (verbose) {
    .log_info(parent, "S06", "analysis completed successfully", verbose)
    .log_info(parent, "S06", paste0("curve_name=", curve_name), verbose)
    .log_info(parent, "S06", paste0("curve_points=", length(xgrid)), verbose)
    .log_info(parent, "S06", paste0("bootstrap_iterations=", res$B), verbose)
  }

  step06$done(paste0("curve_name=", curve_name))

  invisible(scope_obj)
}

#' Rank gene pairs by Lee's L vs Pearson r
#' @description
#' Internal helper for `.get_top_l_vs_r`.
#' Retrieves top gene pairs with large deviations between Lee's L and Pearson
#' correlation, with optional confidence interval filtering and permutation
#' testing.
#' @param scope_obj A `scope_object` containing Lee's L and correlation matrices.
#' @param grid_name Grid layer name.
#' @param pear_level Correlation level (`cell` or `grid`).
#' @param lee_stats_layer Lee statistics layer name.
#' @param curve_layer Optional precomputed curve layer from `.compute_l_vs_r_curve`.
#' @param direction Direction of Delta selection (`both`, `pos`, `neg`).
#' @param top_n Number of top pairs to return.
#' @param ncores Number of threads to use for permutation testing.
#' @param perms Number of permutations for p-value estimation.
#' @param verbose Emit progress messages when TRUE.
#' @param expr_layer Layer name.
#' @param pear_range Parameter value.
#' @param L_range Parameter value.
#' @param do_perm Parameter value.
#' @param block_side Parameter value.
#' @param use_blocks Logical flag.
#' @param clamp_mode Parameter value.
#' @param p_adj_mode Parameter value.
#' @param mem_limit_GB Parameter value.
#' @param pval_mode Parameter value.
#' @param CI_rule Parameter value.
#' @return A data.frame of top pairs with statistics.
#' @keywords internal
.get_top_l_vs_r <- function(scope_obj,
                       grid_name,
                       pear_level = c("cell", "grid"),
                       lee_stats_layer = "LeeStats_Xz",
                       expr_layer = NULL,
                       pear_range = c(-1, 1),
                       L_range = c(-1, 1),
                       top_n = 10,
                       direction = c("largest", "smallest", "both"),
                       do_perm = TRUE,
                       perms = 1000,
                       block_side = 8,
                       use_blocks = TRUE,
                       ncores = 1,
                       clamp_mode = c("none", "ref_only", "both"),
                       p_adj_mode = c("BH", "BY", "BH_universe", "BY_universe", "bonferroni"),
                       mem_limit_GB = 2,
                       pval_mode = c("exact", "beta", "mid", "uniform"),
                       curve_layer = NULL,
                       support_pct_scale = c("auto", "0-1", "0-100"),
                       backend = "cpp",
                       CI_rule = c("remove_within", "remove_outside", "none"),
                       verbose = TRUE) {
  pear_level <- match.arg(pear_level)
  direction <- match.arg(direction)
  p_adj_mode <- match.arg(p_adj_mode)
  clamp_mode <- match.arg(clamp_mode)
  pval_mode <- match.arg(pval_mode)
  grid_name <- .lee_assert_layer_name(grid_name, "grid_name")
  lee_stats_layer <- .lee_assert_layer_name(lee_stats_layer, "lee_stats_layer")
  do_perm <- .lee_assert_flag(do_perm, "do_perm")
  use_blocks <- .lee_assert_flag(use_blocks, "use_blocks")
  verbose <- .lee_assert_flag(verbose, "verbose")
  ncores <- .lee_assert_positive_integer(ncores, "ncores")
  block_side <- .lee_assert_positive_integer(block_side, "block_side")
  mem_limit_GB <- .lee_assert_positive_finite(mem_limit_GB, "mem_limit_GB")
  top_n <- .lee_assert_positive_integer(top_n, "top_n")
  perms <- .lee_validate_integer_scalar(perms, "perms", lower = 0L)
  if (perms == 0L) do_perm <- FALSE
  if (!is.null(expr_layer)) expr_layer <- .lee_assert_layer_name(expr_layer, "expr_layer")
  if (!is.null(curve_layer)) curve_layer <- .lee_assert_layer_name(curve_layer, "curve_layer")
  validate_range <- function(x, name) {
    if (!is.numeric(x) || length(x) != 2L || anyNA(x) ||
        any(!is.finite(x)) || x[1L] > x[2L]) {
      stop(name, " must contain two finite numbers in non-decreasing order.", call. = FALSE)
    }
    as.numeric(x)
  }
  pear_range <- validate_range(pear_range, "pear_range")
  L_range <- validate_range(L_range, "L_range")
  support_pct_scale <- match.arg(support_pct_scale)
  backend <- if (missing(backend)) {
    "cpp"
  } else {
    .normalize_core_backend(backend, arg = "backend", allow_auto = TRUE)
  }
  if (isTRUE(do_perm) && identical(backend, "r")) {
    stop(
      "getTopLvsR backend='r' does not implement the fixed-Pearson Delta permutation null; ",
      "use backend='cpp' or set do_perm=FALSE.",
      call. = FALSE
    )
  }
  if (isTRUE(do_perm) && .native_permutation_disabled()) {
    stop(
      "getTopLvsR Delta permutation requires the native fixed-Pearson kernel, which is disabled. ",
      "Explicitly opt in with options(geneSCOPE.disable_native_all=FALSE, ",
      "geneSCOPE.disable_native_permutation=FALSE), or set do_perm=FALSE.",
      call. = FALSE
    )
  }
  backend_policy <- .core_backend_policy(
    backend,
    cpp_supported = FALSE,
    python_supported = FALSE,
    r_supported = TRUE
  )
  if (identical(backend, "python")) {
    stop(
      "[geneSCOPE::getTopLvsR] backend='", backend,
      "' is not implemented for this R-only ranking layer. ",
      "Use the default C++->R fallback policy or backend='r'.",
      call. = FALSE
    )
  } else if (backend %in% c("auto", "cpp") && verbose) {
    .log_info("getTopLvsR", "S01", "Top-pair ranking is R-only; using R fallback under C++->R policy.", verbose)
  }
  CI_rule <- match.arg(CI_rule)

  parent <- "getTopLvsR"
  step01 <- .log_step(parent, "S01", "load matrices and compute Delta", verbose)
  step01$enter(paste0("grid_name=", grid_name, " pear_level=", pear_level))
  if (verbose) {
    .log_info(parent, "S01", "starting top Delta analysis", verbose)
    .log_info(parent, "S01", paste0("direction=", direction, " top_n=", top_n), verbose)
  }

  logi <- suppressWarnings(detectCores(TRUE))
  if (!is.numeric(logi) || length(logi) != 1L || is.na(logi) ||
      !is.finite(logi) || logi < 1L) {
    logi <- 1L
  }
  logi <- as.integer(logi)
  safe_cores <- min(ncores, logi)
  if (ncores > safe_cores) {
    .log_info(parent, "S01", paste0(
      "adjusting ncores: requested=", ncores,
      " capped at available=", safe_cores
    ), verbose)
    ncores <- safe_cores
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
  }
  Sys.setenv(OMP_NUM_THREADS = ncores)

  provenance <- .lee_validate_current_provenance(
    scope_obj,
    grid_name = grid_name,
    lee_stats_layer = lee_stats_layer,
    context = "[geneSCOPE::getTopLvsR]",
    block_side = block_side
  )
  grid_name <- provenance$grid_name
  if (!is.null(expr_layer) && !identical(expr_layer, provenance$norm_layer)) {
    stop(
      "expr_layer ('", expr_layer, "') does not match the norm_layer ('",
      provenance$norm_layer, "') that produced observed Lee's L.",
      call. = FALSE
    )
  }
  expr_layer <- provenance$norm_layer

  L_mat <- .get_lee_matrix(scope_obj, grid_name, lee_layer = lee_stats_layer)
  r_mat <- .get_pearson_matrix(scope_obj, grid_name, level = pear_level)
  common <- intersect(rownames(L_mat), rownames(r_mat))
  if (length(common) < 2) stop("Insufficient common genes")

  if (verbose) {
    .log_info(parent, "S01", paste0("common_genes=", length(common)), verbose)
    .log_info(parent, "S01", paste0("matrix_dim=", nrow(L_mat), "x", ncol(L_mat)), verbose)
  }

  L_mat <- L_mat[common, common]
  r_mat <- r_mat[common, common]
  diag(L_mat) <- NA
  diag(r_mat) <- NA
  ut <- upper.tri(L_mat)
  LeesL_vec <- L_mat[ut]
  Pear_vec <- r_mat[ut]

  if (verbose) .log_info(parent, "S01", "computing Delta values and applying filters", verbose)
  internal_clamp_mode <- clamp_mode
  if (clamp_mode == "both") {
    .log_info(parent, "S01",
      "clamp_mode='both' currently equivalent to ref_only (reference truncation only)",
      verbose
    )
    internal_clamp_mode <- "ref_only"
  }
  Pear_for_delta <- if (internal_clamp_mode == "ref_only") pmax(Pear_vec, 0) else Pear_vec

  warning(
    "Delta (Lee's L - Pearson r) is an exploratory metric. ",
    "Cross-scale comparison of spatial vs. non-spatial correlation ",
    "lacks formal statistical theory for inference. ",
    "Use Delta for heuristic prioritization / hypothesis generation only, ",
    "not for confirmatory inference.",
    call. = FALSE, immediate. = TRUE
  )

  df <- data.frame(
    gene1 = rep(common, each = length(common))[ut],
    gene2 = rep(common, length(common))[ut],
    LeesL = LeesL_vec,
    Pear  = Pear_vec,
    Delta = LeesL_vec - Pear_for_delta,
    delta_method = "L_minus_r",
    delta_caveat = "exploratory: heuristic prioritization across spatial and non-spatial scales; permutation-assisted ranking, not confirmatory inference"
  )

  # Expression coverage percentages
  {
    expr_support_map <- setNames(rep(0, length(common)), common)
    g_layer_try <- tryCatch(.select_grid_layer(scope_obj, grid_name), error = function(e) NULL)
    coverage_done <- FALSE

    if (!is.null(g_layer_try) && !is.null(g_layer_try$counts)) {
      ct <- g_layer_try$counts
      if (is.data.frame(ct) && all(c("gene", "grid_id") %in% colnames(ct))) {
        total_cells <- if (!is.null(g_layer_try$grid_info)) nrow(g_layer_try$grid_info) else length(unique(ct$grid_id))
        if (total_cells <= 0) total_cells <- NA_real_
        if (inherits(ct, "data.table")) {
          gene_cells <- ct[gene %in% common, .(cells = uniqueN(grid_id)), by = gene]
        } else if (requireNamespace("data.table", quietly = TRUE)) {
          dct <- as.data.table(ct)
          gene_cells <- dct[gene %in% common, .(cells = uniqueN(grid_id)), by = gene]
        } else {
          keep <- ct$gene %in% common
          if (any(keep)) {
            spl <- split(ct$grid_id[keep], ct$gene[keep])
            gene_cells <- data.frame(
              gene = names(spl),
              cells = vapply(spl, function(v) length(unique(v)), integer(1)),
              row.names = NULL
            )
          } else {
            gene_cells <- data.frame(gene = character(0), cells = integer(0))
          }
        }
        if (nrow(gene_cells)) {
          expr_support_map[gene_cells$gene] <- gene_cells$cells / total_cells
        }
        coverage_done <- TRUE
      }
    }

    if (!coverage_done) {
      pick_matrix <- function(layer) {
        if (is.null(layer)) {
          return(NULL)
        }
        for (nm in c("counts", "raw_counts", "expr", "X", "data", "logCPM", "Xz")) {
          m <- layer[[nm]]
          if (!is.null(m) && (is.matrix(m) || inherits(m, "dgCMatrix"))) {
            return(m)
          }
        }
        NULL
      }
      expr_mat <- pick_matrix(g_layer_try)
      if (is.null(expr_mat) && pear_level == "cell") {
        cell_env <- tryCatch(scope_obj@cell, error = function(e) NULL)
        if (is.null(cell_env)) cell_env <- tryCatch(scope_obj@cells, error = function(e) NULL)
        expr_mat <- pick_matrix(cell_env)
      }
      if (!is.null(expr_mat)) {
        rn <- rownames(expr_mat)
        cn <- colnames(expr_mat)
        inter_col <- intersect(cn, common)
        inter_row <- intersect(rn, common)
        if (length(inter_col) == 0 && length(inter_row) == 0) {
          # keep zeros
        } else {
          gene_in_col <- (length(inter_col) >= length(inter_row))
          if (!gene_in_col) {
            expr_mat <- if (inherits(expr_mat, "dgCMatrix")) t(expr_mat) else t(expr_mat)
            cn <- colnames(expr_mat)
            inter_col <- intersect(cn, common)
          }
          if (length(inter_col)) {
            nz_counts <- if (inherits(expr_mat, "dgCMatrix")) {
              colSums(expr_mat[, inter_col, drop = FALSE] > 0)
            } else {
              colSums(expr_mat[, inter_col, drop = FALSE] > 0)
            }
            expr_support_map[inter_col] <- nz_counts / nrow(expr_mat)
          }
        }
      }
    }

    n_pair <- nrow(df)
    gene1_support <- expr_support_map[df$gene1]
    gene2_support <- expr_support_map[df$gene2]
    gene1_pct <- gene1_support * 100
    gene2_pct <- gene2_support * 100
    df$gene1_expr_support <- round(gene1_support, 6)
    df$gene2_expr_support <- round(gene2_support, 6)
    df$gene1_expr_pct <- round(gene1_pct, 3)
    df$gene2_expr_pct <- round(gene2_pct, 3)
    attr(df, "support_pct_requested_scale") <- support_pct_scale
    attr(df, "support_pct_input_scale") <- "0-1"
    attr(df, "support_pct_scale_used") <- "0-1"
    attr(df, "support_pct_internal_scale") <- "0-1"
    attr(df, "support_pct_legacy_scale") <- "0-100"
  }

  step01$done(paste0("pairs_total=", nrow(df)))

  # Curve filtering
  step02 <- .log_step(parent, "S02", "apply curve/CI filters", verbose)
  step02$enter(paste0("curve_layer=", if (is.null(curve_layer)) "none" else curve_layer, " CI_rule=", CI_rule))
  if (!is.null(curve_layer) || CI_rule != "none") {
    if (is.null(curve_layer)) stop("curve_layer must be provided when CI_rule != 'none'")
    curve_obj <- scope_obj@stats[[grid_name]][[lee_stats_layer]][[curve_layer]]
    if (is.null(curve_obj)) stop("Curve layer '", curve_layer, "' not found in stats")
    required_cols <- c("Pear", "lo95", "hi95")
    if (!all(required_cols %in% colnames(curve_obj))) {
      stop("Curve layer must contain columns: ", paste(required_cols, collapse = ", "))
    }
    lo_fun <- approxfun(curve_obj$Pear, curve_obj$lo95, rule = 2)
    hi_fun <- approxfun(curve_obj$Pear, curve_obj$hi95, rule = 2)
    df$curve_lo <- lo_fun(df$Pear)
    df$curve_hi <- hi_fun(df$Pear)
    df$outside_ci <- is.na(df$curve_lo) | is.na(df$curve_hi) | df$LeesL < df$curve_lo | df$LeesL > df$curve_hi
    if (CI_rule == "remove_within") {
      df <- df[df$outside_ci, , drop = FALSE]
    } else if (CI_rule == "remove_outside") {
      df <- df[!df$outside_ci, , drop = FALSE]
    }
    if (!nrow(df)) stop("No gene pairs satisfy the requested CI_rule")
  }

  if (verbose && !is.null(df$outside_ci) && CI_rule != "none") {
    if (CI_rule == "remove_within") {
      .log_info(parent, "S02", paste0("pairs_outside_ci95=", sum(df$outside_ci)), verbose)
    } else if (CI_rule == "remove_outside") {
      .log_info(parent, "S02", paste0("pairs_inside_ci95=", sum(!df$outside_ci)), verbose)
    }
  }

  step02$done(paste0("pairs_after_ci=", nrow(df)))

  step03 <- .log_step(parent, "S03", "filter ranges and select pairs", verbose)
  step03$enter(paste0(
    "pear_range=[", pear_range[1], ",", pear_range[2], "]",
    " L_range=[", L_range[1], ",", L_range[2], "]"
  ))
  if (verbose) .log_info(parent, "S03", "applying threshold filters", verbose)
  df <- df[df$Pear >= pear_range[1] & df$Pear <= pear_range[2] &
    df$LeesL >= L_range[1] & df$LeesL <= L_range[2], ]
  if (!nrow(df)) stop("Thresholds remove all pairs")
  total_universe <- nrow(df)

  if (verbose) {
    .log_info(parent, "S03", paste0("pairs_after_filtering=", total_universe), verbose)
    .log_info(parent, "S03", paste0(
      "pear_range=[", pear_range[1], ", ", pear_range[2], "]",
      " L_range=[", L_range[1], ", ", L_range[2], "]"
    ), verbose)
  }

  if (verbose) .log_info(parent, "S03", "selecting top gene pairs by Delta values", verbose)
  sel <- switch(direction,
    largest = slice_max(df, Delta, n = top_n),
    smallest = slice_min(df, Delta, n = top_n),
    both = bind_rows(
      slice_max(df, Delta, n = top_n),
      slice_min(df, Delta, n = top_n)
    )
  )
  sel <- distinct(sel, gene1, gene2, .keep_all = TRUE)
  rownames(sel) <- NULL

  if (verbose) {
    .log_info(parent, "S03", paste0("selected_pairs=", nrow(sel)), verbose)
    if (nrow(sel) > 0) {
      delta_range <- range(sel$Delta)
      .log_info(parent, "S03", paste0(
        "delta_range=[", round(delta_range[1], 4), ", ", round(delta_range[2], 4), "]"
      ), verbose)
    }
  }

  step03$done(paste0("selected_pairs=", nrow(sel)))

  if (!nrow(sel) || !do_perm) {
    step04 <- .log_step(parent, "S04", "permutation testing", verbose)
    step04$enter(paste0("perms=", perms, " do_perm=", do_perm))
    .log_backend(parent, "S04", "permutation", "skipped",
      reason = if (!nrow(sel)) "no_pairs" else "do_perm=FALSE",
      verbose = verbose
    )
    step04$done("skipped")

    step05 <- .log_step(parent, "S05", "return results", verbose)
    step05$enter("no permutation results")
    sel$FDR <- NA_real_
    out <- transmute(
      sel,
      gene1,
      gene2,
      L = LeesL,
      r = Pear,
      delta = Delta,
      delta_method = delta_method,
      delta_caveat = delta_caveat,
      support1 = gene1_expr_support,
      support2 = gene2_expr_support,
      pct1 = gene1_expr_pct,
      pct2 = gene2_expr_pct,
      fdr = FDR
    )
    attr(out, "backend") <- list(requested = backend, selected = "r", policy = backend_policy)
    attr(out, "support_pct_scale") <- list(
      requested = attr(df, "support_pct_requested_scale"),
      input = attr(df, "support_pct_input_scale"),
      used = attr(df, "support_pct_scale_used"),
      internal = attr(df, "support_pct_internal_scale"),
      legacy = attr(df, "support_pct_legacy_scale")
    )
    step05$done(paste0("pairs_returned=", nrow(out)))
    return(out)
  }

  step04 <- .log_step(parent, "S04", "permutation testing", verbose)
  step04$enter(paste0("perms=", perms, " pval_mode=", pval_mode))
  if (verbose) {
    .log_info(parent, "S04", "preparing permutation analysis", verbose)
    .log_info(parent, "S04", paste0("permutations=", perms), verbose)
    .log_info(parent, "S04", paste0("pval_mode=", pval_mode), verbose)
    .log_info(parent, "S04", paste0("p_adj_mode=", p_adj_mode), verbose)
  }
  Xz <- provenance$X
  W <- provenance$W
  grid_info <- provenance$grid_info
  genes_top <- unique(c(sel$gene1, sel$gene2))
  gene_map <- match(genes_top, colnames(Xz))
  if (any(is.na(gene_map))) stop("Selected genes not found in Xz")
  gene_pairs <- cbind(
    match(sel$gene1, genes_top) - 1L,
    match(sel$gene2, genes_top) - 1L
  )
  delta_ref <- sel$Delta
  # This is exactly the reference used to construct observed Delta, including
  # ref_only/both clamping.  A common row permutation keeps it fixed.
  pearson_ref <- sel$LeesL - delta_ref
  if (any(!is.finite(delta_ref)) || any(!is.finite(pearson_ref))) {
    stop("Selected Delta/Pearson reference values must be finite before permutation.", call. = FALSE)
  }

  backend_label <- if (use_blocks) "C++ delta_l_fixed_r_perm_csr_block" else "C++ delta_l_fixed_r_perm_csr"
  .log_backend(parent, "S04", "permutation_backend", paste0(
    backend_label,
    " threads=", ncores,
    " mem_limit_GB=", mem_limit_GB
  ), verbose = verbose)

  block_id <- if (use_blocks) {
    bx <- (grid_info$gx - 1L) %/% block_side
    by <- (grid_info$gy - 1L) %/% block_side
    max_by <- max(by)
    bx * (max_by + 1L) + by + 1L
  } else {
    seq_len(nrow(Xz))
  }
  split_rows <- split(seq_along(block_id), block_id)

  W_rows <- lapply(seq_len(nrow(W)), function(i) {
    nz <- which(W[i, ] != 0)
    if (length(nz)) list(indices = nz - 1L, values = as.numeric(W[i, nz])) else list(indices = integer(0), values = numeric(0))
  })
  W_row_lengths <- vapply(W_rows, function(x) length(x$indices), integer(1))
  W_row_ptr <- c(0L, cumsum(W_row_lengths))
  W_indices <- unlist(lapply(W_rows, `[[`, "indices"), use.names = FALSE)
  W_values <- unlist(lapply(W_rows, `[[`, "values"), use.names = FALSE)
  if (!length(W_indices)) stop("Weight matrix has no non-zero entries")
  Xz_sub <- Xz[, gene_map, drop = FALSE]
  n_cells <- nrow(Xz_sub)

  target_batch <- min(100L, perms)
  max_idx_bytes <- mem_limit_GB * 1024^3 * 0.30
  est_bytes <- function(bs) n_cells * bs * 4
  while (target_batch > 1L && est_bytes(target_batch) > max_idx_bytes) {
    target_batch <- max(1L, floor(target_batch / 2))
  }
  if (verbose) {
    .log_info(parent, "S04", paste0(
      "planned_batch_size=", target_batch,
      " (est_idx_mat_mb=", sprintf("%.2f", est_bytes(target_batch) / 1024^2),
      " limit_mb=", sprintf("%.2f", max_idx_bytes / 1024^2), ")"
    ), verbose)
    .log_info(parent, "S04", "starting permutation testing loop", verbose)
  }
  remaining <- perms
  exceed_count <- rep(0L, nrow(sel))
  perm_threads <- ncores
  attempt <- 1
  while (remaining > 0) {
    bsz <- min(target_batch, remaining)
    success <- FALSE
    while (!success && perm_threads >= 1) {
      if (verbose) .log_info(parent, "S04", sprintf(
        "permutation attempt #%d threads=%d batch=%d remain=%d clamp_mode=%s",
        attempt, perm_threads, bsz, remaining, clamp_mode
      ), verbose)
      attempt <- attempt + 1
      idx_mat <- matrix(NA_integer_, nrow = n_cells, ncol = bsz)
      for (p in seq_len(bsz)) {
        if (use_blocks) {
          idx <- seq_len(n_cells)
          for (rows in split_rows) {
            idx[rows] <- rows[sample.int(length(rows), length(rows), replace = FALSE)]
          }
          if (!identical(sort(idx), seq_len(n_cells)) || any(block_id[idx] != block_id)) {
            stop("Internal block permutation invariant failed", call. = FALSE)
          }
        } else {
          idx <- sample.int(n_cells, n_cells, replace = FALSE)
        }
        idx_mat[, p] <- idx - 1L
      }
      res <- tryCatch(
        {
          if (use_blocks) {
            .delta_l_fixed_r_perm_csr_block(
              Xz_sub, W_indices, W_values, W_row_ptr, idx_mat,
              as.integer(block_id) - 1L, gene_pairs, delta_ref, pearson_ref,
              perm_threads
            )
          } else {
            .delta_l_fixed_r_perm_csr(
              Xz_sub, W_indices, W_values, W_row_ptr, idx_mat,
              gene_pairs, delta_ref, pearson_ref,
              perm_threads
            )
          }
        },
        error = function(e) e
      )
      if (inherits(res, "error")) {
        if (verbose) .log_info(parent, "S04", paste0("batch failed: ", conditionMessage(res)), verbose)
        if (perm_threads > 1) {
          perm_threads <- max(1, floor(perm_threads / 2))
          .log_backend(parent, "S04", "permutation_backend", paste0(
            backend_label,
            " threads=", perm_threads,
            " mem_limit_GB=", mem_limit_GB
          ), reason = "retry", verbose = verbose)
          next
        } else {
          stop("Permutation failed at single-thread: ", conditionMessage(res))
        }
      } else {
        exceed_count <- exceed_count + res
        success <- TRUE
      }
    }
    remaining <- remaining - bsz
  }

  if (verbose) .log_info(parent, "S04", "computing p-values and applying multiple testing correction", verbose)
  N <- perms
  k <- exceed_count
  p_values <- switch(pval_mode,
    exact   = (k + 1) / (N + 1),
    beta    = (k + 1) / (N + 2),
    mid     = (k + 0.5) / (N + 1),
    uniform = (k + runif(length(k))) / (N + 1)
  )
  p_values[p_values > 1] <- 1
  mc_se <- sqrt(p_values * (1 - p_values) / N)
  p_ci_lo <- qbeta(0.025, k + 1, N - k + 1)
  p_ci_hi <- qbeta(0.975, k + 1, N - k + 1)

  FDR <- switch(p_adj_mode,
    BH          = p.adjust(p_values, "BH"),
    BY          = p.adjust(p_values, "BY"),
    BH_universe = p.adjust(p_values, "BH", n = total_universe),
    BY_universe = p.adjust(p_values, "BY", n = total_universe),
    bonferroni  = p.adjust(p_values, "bonferroni", n = total_universe)
  )

  if (verbose) {
    sig_count <- sum(p_values < 0.05, na.rm = TRUE)
    fdr_sig_count <- sum(FDR < 0.05, na.rm = TRUE)
    .log_info(parent, "S04", paste0(
      "significant_pairs_p_lt_0.05=", sig_count, "/", length(p_values)
    ), verbose)
    .log_info(parent, "S04", paste0(
      "significant_pairs_fdr_lt_0.05=", fdr_sig_count, "/", length(FDR)
    ), verbose)
  }

  sel$raw_p <- p_values
  sel$mc_se <- mc_se
  sel$p_ci_lo <- p_ci_lo
  sel$p_ci_hi <- p_ci_hi
  sel$FDR <- FDR
  sel$stat_type <- switch(clamp_mode,
    none     = "Delta",
    ref_only = "Delta_refClamp",
    both     = "Delta_clamp_Ronly"
  )
  sel$pval_mode <- pval_mode
  if (p_adj_mode == "bonferroni") {
    sel <- sel[sel$FDR < 0.05, , drop = FALSE]
    if (!nrow(sel)) {
      .log_info(parent, "S04", "no Bonferroni-significant pairs (FDR < 0.05)", verbose)
    }
  }

  step04$done(paste0("pairs_with_fdr=", nrow(sel)))

  step05 <- .log_step(parent, "S05", "return results", verbose)
  step05$enter(paste0("pairs_final=", nrow(sel)))
  out <- transmute(
    sel,
    gene1,
    gene2,
    L = LeesL,
    r = Pear,
    delta = Delta,
    delta_method = delta_method,
    delta_caveat = delta_caveat,
    support1 = gene1_expr_support,
    support2 = gene2_expr_support,
    pct1 = gene1_expr_pct,
    pct2 = gene2_expr_pct,
    fdr = FDR
  )
  attr(out, "backend") <- list(requested = backend, selected = "r", policy = backend_policy)
  attr(out, "support_pct_scale") <- list(
    requested = attr(df, "support_pct_requested_scale"),
    input = attr(df, "support_pct_input_scale"),
    used = attr(df, "support_pct_scale_used"),
    internal = attr(df, "support_pct_internal_scale"),
    legacy = attr(df, "support_pct_legacy_scale")
  )
  step05$done(paste0("pairs_returned=", nrow(out)))
  out
}

#' Compute Lee's L for Visium (internal implementation)
#' @description
#' Performs Visium-specific prescreening and setup (spatial weights, optional Idelta
#' prescreening, and gene filtering) before dispatching to `computeL()`.
#' This is the internal implementation used by the computeL_visium wrapper.
#' @param scope_obj A `scope_object` containing Visium data.
#' @param grid_name Grid layer name to operate on.
#' @param norm_layer Normalised expression layer name.
#' @param use_idelta Use Idelta prescreening when TRUE.
#' @param S_target Target number of genes after prescreening.
#' @param min_detect Minimum detection rate when filtering genes.
#' @param winsor_high Winsorization quantile for expression scaling.
#' @param ncores Number of threads to use (defaults to safe thread count).
#' @param verbose Emit progress messages when TRUE.
#' @param backend Lee's L backend policy. Uses the same C++ -> R executable
#'   contract as `computeL()`; `python` errors explicitly because no equivalent
#'   Python Lee's L backend is shipped.
#' @param ... Additional arguments forwarded to `.compute_l()`, such as
#'   `mem_limit_GB`, `use_bigmemory`, `chunk_size`, `perms`, or
#'   `lee_stats_layer_name`.
#' @return The modified `scope_object`.
#' @seealso `computeL()`, `computeWeights()`, `computeIDelta()`
#' @keywords internal
.computeL_visium_internal <- function(
    scope_obj,
    grid_name = NULL,
    norm_layer = "Xz",
    use_idelta = TRUE,
    S_target = 8000L,
    min_detect = 0.10,
    winsor_high = 0.95,
    ncores = NULL,
    backend = "cpp",
    verbose = getOption("geneSCOPE.verbose", TRUE),
    ...) {

    parent <- "computeL_visium"
    backend <- if (missing(backend)) {
        "cpp"
    } else {
        .normalize_core_backend(backend, arg = "backend", allow_auto = TRUE)
    }
    if (identical(backend, "python")) {
        stop(.lee_l_python_backend_unavailable(), call. = FALSE)
    }

    ## ---- 0. Thread count and grid layer selection ----
    if (is.null(ncores)) {
        ncores <- .get_safe_thread_count(default = 8L)
    }
    grid_label <- if (is.null(grid_name)) "auto" else as.character(grid_name)[1]
    step01 <- .log_step(parent, "S01", "select grid layer and validate", verbose)
    step01$enter(paste0("grid_name=", grid_label, " norm_layer=", norm_layer))

    g_layer <- .select_grid_layer(scope_obj, grid_name)
    grid_name <- if (is.null(grid_name)) {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    } else grid_name

    if (is.null(g_layer[[norm_layer]])) {
        stop("computeL_visium: normalized layer '", norm_layer, "' not found.")
    }
    X <- g_layer[[norm_layer]]  # n x g
    if (!is.matrix(X)) X <- as.matrix(X)
    n <- nrow(X); G <- ncol(X)
    if (n < 2L || G < 2L) stop("Insufficient data size to compute Lee's L.")
    .log_info(parent, "S01", paste0("cells=", n, " genes=", G, " ncores=", ncores), verbose)
    step01$done(paste0("grid_name=", grid_name, " genes=", G))

    ## ---- 1. Check/compute spatial weights W ----
    step02 <- .log_step(parent, "S02", "ensure spatial weights", verbose)
    step02$enter(paste0("grid_name=", grid_name))
    if (is.null(g_layer$W) || sum(g_layer$W) == 0) {
        .log_info(parent, "S02", "W missing; running computeWeights(style=B)", verbose)
        .log_backend(parent, "S02", "weights", "computeWeights(style=B)", verbose = verbose)
        scope_obj <- .compute_weights(
            scope_obj,
            grid_name = grid_name,
            style = "B",
            store_mat = TRUE,
            verbose = verbose
        )
        g_layer <- .select_grid_layer(scope_obj, grid_name)
    }
    step02$done(paste0("W=", if (is.null(g_layer$W)) "missing" else "ready"))
    W <- g_layer$W  # dgCMatrix (ensures W exists for .compute_l later)

    ## ---- 2. Idelta and gene prescreening ----
    step03 <- .log_step(parent, "S03", "optional I-delta prescreen", verbose)
    step03$enter(paste0("use_idelta=", use_idelta))
    if (use_idelta) {
        meta_col <- paste0(grid_name, "_iDelta")
        if (is.null(scope_obj@meta.data) || !(meta_col %in% colnames(scope_obj@meta.data))) {
            .log_info(parent, "S03", "computing I-delta", verbose)
            .log_backend(parent, "S03", "idelta", "C++ idelta_sparse_cpp", verbose = verbose)
            scope_obj <- .compute_idelta(
                scope_obj,
                grid_name = grid_name,
                level = "grid",
                ncores = min(8L, ncores),
                verbose = FALSE
            )
        } else {
            .log_info(parent, "S03", "I-delta already present; skipping recompute", verbose)
        }
    } else {
        .log_info(parent, "S03", "I-delta prescreen disabled", verbose)
    }
    step03$done()
    genes_all <- colnames(X)
    if (is.null(genes_all)) genes_all <- rownames(scope_obj@meta.data)
    if (is.null(genes_all)) genes_all <- as.character(seq_len(G))

    ## 2.1 Detection rate (based on counts)
    detect_rate <- NULL
    total_count <- NULL
    if (!is.null(g_layer$counts)) {
        dt <- g_layer$counts
        dt <- dt[dt$count > 0, c("gene", "grid_id", "count")]
        # Fast path when data.table is available
        if (requireNamespace("data.table", quietly = TRUE)) {
            dtt <- as.data.table(dt)
            n_spots <- length(g_layer$grid_info$grid_id)
            dr <- dtt[, .N, by = gene]
            detect_rate <- setNames(as.numeric(dr$N) / n_spots, dr$gene)
            total_count <- dtt[, sum(count), by = gene]
            total_count <- setNames(as.numeric(total_count$V1), total_count$gene)
        } else {
            n_spots <- length(unique(dt$grid_id))
            detect_rate <- tapply(dt$grid_id, dt$gene, function(x) length(unique(x)) / n_spots)
            total_count <- tapply(dt$count, dt$gene, sum)
        }
    }

    ## 2.2 Determine candidate gene set
    keep <- rep(TRUE, length(genes_all))
    names(keep) <- genes_all
    if (!is.null(detect_rate)) {
        low <- detect_rate < min_detect
        keep[names(low)] <- !low
    }
    if (!is.null(total_count)) {
        total_count <- total_count[names(keep)]
        total_count[is.na(total_count)] <- 0
    }
    genes_keep <- names(keep)[keep]

    if (!is.null(total_count)) {
        if (length(genes_keep) > S_target) {
            ord <- order(total_count[genes_keep], decreasing = TRUE)
            genes_keep <- genes_keep[ord][seq_len(S_target)]
        }
    }

    ## ---- 3. Winsorize & compute L ----
    if (!is.null(g_layer[[norm_layer]])) {
        X_use <- g_layer[[norm_layer]]
        if (!is.matrix(X_use)) X_use <- as.matrix(X_use)
        if (is.numeric(winsor_high) && winsor_high > 0 && winsor_high < 1) {
            qv <- apply(X_use, 2, stats::quantile, probs = winsor_high, na.rm = TRUE)
            X_use <- sweep(X_use, 2, qv, pmin)
            g_layer[[norm_layer]] <- X_use
            scope_obj@grid[[grid_name]][[norm_layer]] <- X_use
        }
    }

    scope_obj <- .with_lee_l_backend(
        backend,
        .compute_l(
            scope_obj,
            grid_name = grid_name,
            norm_layer = norm_layer,
            genes = genes_keep,
            within = TRUE,
            ncores = ncores,
            verbose = verbose,
            ...
        )
    )
    scope_obj
}
