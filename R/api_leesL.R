#' Compute Lee's L statistic.
#' @description
#' Computes Lee's L for gene pairs using a grid-level expression layer, stores
#' results under a statistics layer in `scope_obj@stats`, and returns the updated
#' object.
#' @param scope_obj A \code{scope_object} with at least one populated \code{@grid} slot.
#' @param grid_name Character; name of the grid layer to operate on. If \code{NULL}
#'   and only one grid layer exists, it is auto-selected. If multiple layers exist,
#'   the name must be specified.
#' @param genes Optional subset of genes to include.
#' @param within If `TRUE`, restricts analysis to the selected gene set on both axes.
#' @param ncores Number of cores for parallel processing.
#' @param block_side Number of grid cells per side for block partitioning.
#' @param perms Number of permutations for Monte-Carlo p-values.
#' @param approximate_q Use beta posterior-smoothed permutation p-values plus
#'   Storey q-value estimation when `perms` is below `approximate_q_threshold`.
#'   This is intended for large gene panels with low permutation budgets and is
#'   not enabled by default. It is a screening aid only and must not be treated
#'   as a replacement for adequately powered permutation-based confirmatory
#'   inference.
#' @param approximate_q_method Approximation method. Currently `"beta"` uses
#'   `(k + 1) / (B + 2)` posterior-smoothed permutation p-values before Storey
#'   q-value estimation.
#' @param approximate_q_threshold Only activate `approximate_q` when
#'   `0 < perms < approximate_q_threshold`.
#' @param block_size Number of permutations processed per batch.
#' @param L_min Similarity threshold used when building QC similarity graphs.
#' @param norm_layer Name of the normalised expression layer (default `"Xz"`).
#' @param lee_stats_layer_name Output statistics layer name (auto-generated when NULL).
#' @param legacy_formula Use legacy denominator for compatibility.
#' @param mem_limit_GB RAM threshold that triggers streaming mode.
#' @param chunk_size Number of columns processed per chunk in streaming mode.
#' @param use_bigmemory Whether to use file-backed matrices for large computations.
#'   Defaults to `FALSE`; if a matrix exceeds `mem_limit_GB` and the user did
#'   not explicitly set this argument, geneSCOPE may still enable chunked
#'   bigmemory mode as a safety fallback.
#' @param backing_path Directory for temporary files (default `tempdir()`).
#' @param cache_inputs Whether to .cache preprocessed inputs for reuse across calls.
#' @param verbose Logical; whether to emit compact public progress messages.
#'   Default is \code{TRUE}. Set to \code{FALSE} for silent operation, or use
#'   \code{options(geneSCOPE.verbose = FALSE)} for global control.
#' @param backend Deprecated compatibility argument. The runtime now uses the
#'   package policy `C++ first, R fallback`; `python` still errors explicitly
#'   because no Python Lee's L backend is shipped. Legacy alias `native` maps to
#'   `cpp`.
#'   Inputs must provide at least two observations (`n_cells >= 2`); non-finite
#'   values error explicitly, and extremely large finite values (`> 1e300`)
#'   trigger a precision warning.
#' @param ncore Deprecated alias of `ncores`.
#' @return The modified `scope_object`.
#' @examples
#' \dontrun{
#' scope_obj <- createSCOPE(data_dir = "/path/to/output_dir", grid_length = 30)
#' scope_obj <- normalizeMoleculesInGrid(scope_obj, grid_name = "grid30")
#' scope_obj <- computeWeights(scope_obj, grid_name = "grid30")
#' scope_obj <- computeL(scope_obj, grid_name = "grid30", ncores = 16)
#' }
#' @seealso `computeWeights()`, `computeMH()`, `computeLvsRCurve()`, `getTopLvsR()`
#' @export
computeL <- function(
    scope_obj,
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
    backend = "cpp",
    ncore = NULL) {
    approximate_q_method <- match.arg(approximate_q_method)
    backend <- if (missing(backend)) {
        "cpp"
    } else {
        .normalize_core_backend(backend, arg = "backend", allow_auto = TRUE)
    }
    verbose <- .resolve_runtime_verbose(verbose)

    .with_log_session("computeL", verbose = verbose, {
        .log_start("computeL", verbose = verbose)
        
        resolved_grid_name <- if (!is.null(grid_name) && nzchar(grid_name)) {
            as.character(grid_name)[1]
        } else if (length(scope_obj@grid) == 1L) {
            names(scope_obj@grid)[1]
        } else {
            NULL
        }
        layer_name <- if (is.null(lee_stats_layer_name) || !nzchar(lee_stats_layer_name)) {
            paste0("LeeStats_", norm_layer)
        } else {
            lee_stats_layer_name
        }

        .log_key_values("computeL", list(
            grid = resolved_grid_name,
            layer = layer_name,
            perms = perms,
            approximate_q = approximate_q,
            approximate_q_method = approximate_q_method,
            approximate_q_threshold = approximate_q_threshold,
            backend = backend,
            use_bigmemory = use_bigmemory
        ))

        result <- .with_lee_l_backend(
            backend,
            .compute_l(
                scope_obj = scope_obj,
                grid_name = grid_name,
                genes = genes,
                within = within,
                ncores = ncores,
                block_side = block_side,
                perms = perms,
                approximate_q = approximate_q,
                approximate_q_method = approximate_q_method,
                approximate_q_threshold = approximate_q_threshold,
                block_size = block_size,
                L_min = L_min,
                norm_layer = norm_layer,
                lee_stats_layer_name = lee_stats_layer_name,
                legacy_formula = legacy_formula,
                mem_limit_GB = mem_limit_GB,
                chunk_size = chunk_size,
                use_bigmemory = use_bigmemory,
                backing_path = backing_path,
                cache_inputs = cache_inputs,
                verbose = verbose,
                ncore = ncore
            )
        )

        grid_name_after <- resolved_grid_name
        if (is.null(grid_name_after) && length(result@grid) == 1L) {
            grid_name_after <- names(result@grid)[1]
        }
        lee_layer <- if (!is.null(grid_name_after) &&
            !is.null(result@stats[[grid_name_after]]) &&
            !is.null(result@stats[[grid_name_after]][[layer_name]])) {
            result@stats[[grid_name_after]][[layer_name]]
        } else {
            NULL
        }
        fdr_mat <- if (is.list(lee_layer)) lee_layer$FDR else NULL

        .log_key_values("computeL", list(
            grid = grid_name_after,
            layer = layer_name,
            dims = if (is.list(lee_layer) && !is.null(lee_layer$L)) paste(dim(lee_layer$L), collapse = "x") else NA_character_,
            weight_style = if (is.list(lee_layer)) .log_first_scalar(lee_layer$weight_style, lee_layer$meta$weight_style) else NA_character_,
            sig_fdr_005 = if (!is.null(fdr_mat)) sum(fdr_mat < 0.05, na.rm = TRUE) else NA_integer_
        ))

        .log_done("computeL", verbose = verbose)
        result
    })
}

#' Compute Lee's L using the exact reference R backend.
#' @description
#' Convenience wrapper that preserves the current Lee's L statistics while
#' forcing the exact R reference backend regardless of the active option or
#' `GENESCOPE_LEE_L_BACKEND` environment variable.
#' @inheritParams computeL
#' @param ... Additional arguments passed to `safe_computeL()`.
#' @param return_report Whether to attach a `genescope_computeL_report`
#'   attribute to the returned object.
#' @return A `scope_object` computed with the reference R Lee's L backend.
#' @seealso `computeL()`, `safe_computeL()`
#' @export
reference_computeL <- function(scope_obj,
                               ...,
                               return_report = FALSE) {
    safe_computeL(
        scope_obj = scope_obj,
        ...,
        backend = "r",
        return_report = return_report,
        on_failure = "error"
    )
}

#' Safely compute Lee's L statistic.
#' @description
#' Executes `computeL()` under an explicit Lee's L backend policy and attaches
#' a structured report. Use `backend = "r"` to avoid native Lee's L helpers
#' when auditing suspicious inputs or running on a branch with unresolved
#' native instability.
#' @inheritParams computeL
#' @param ... Additional arguments passed to `computeL()`.
#' @param backend Backend policy for Lee's L core computation:
#'   `"auto"` tries C++ first, skips the unavailable Python tier for Lee's L,
#'   and retries the R reference backend when a regular R error is raised by
#'   the native path. `"cpp"` requires the native helper, `"r"` avoids the
#'   native helper, and `"python"` fails explicitly because no Python Lee's L
#'   implementation is shipped. Legacy alias `"native"` maps to `"cpp"`.
#' @param return_report Whether to attach a `genescope_computeL_report`
#'   attribute to the returned object.
#' @param on_failure Whether to return the input object with an attached report
#'   (`"return_input"`) or rethrow the error (`"error"`).
#' @return A `scope_object`. When `return_report = TRUE`, the result carries a
#'   `genescope_computeL_report` attribute describing the backend used and any
#'   fallback/error state.
#' @seealso `computeL()`
#' @export
safe_computeL <- function(scope_obj,
                          ...,
                          backend = "cpp",
                          return_report = TRUE,
                          on_failure = c("return_input", "error")) {
    backend <- if (missing(backend)) {
        "cpp"
    } else {
        .normalize_core_backend(backend, arg = "backend", allow_auto = TRUE)
    }
    on_failure <- match.arg(on_failure)
    backend_policy <- .core_backend_policy(
        backend,
        python_supported = FALSE,
        cpp_supported = TRUE,
        r_supported = TRUE
    )
    report <- list(
        ok = FALSE,
        backend = backend,
        requested_backend = backend,
        backend_policy = backend_policy,
        selected_backend = NULL,
        fallback_tier = NULL,
        unavailable_tiers = character(0),
        unavailable_tier_notes = character(0),
        fallback_used = FALSE,
        error = NULL,
        fallback_error = NULL,
        timestamp = as.character(Sys.time())
    )

    attach_report <- function(obj) {
        if (isTRUE(return_report)) {
            attr(obj, "genescope_computeL_report") <- report
        }
        obj
    }

    run_compute <- function() {
        withCallingHandlers(
            computeL(scope_obj = scope_obj, ..., backend = backend),
            warning = function(w) {
                msg <- conditionMessage(w)
                if (grepl("Python Lee's L backend is not implemented", msg, fixed = TRUE)) {
                    report$unavailable_tiers <<- unique(c(report$unavailable_tiers, "python"))
                    report$unavailable_tier_notes <<- unique(c(report$unavailable_tier_notes, msg))
                }
                if (grepl("falling back to the reference R backend", msg, fixed = TRUE)) {
                    report$fallback_used <<- TRUE
                    report$fallback_tier <<- "r"
                    report$selected_backend <<- "r"
                    if (is.null(report$fallback_error)) {
                        report$fallback_error <<- msg
                    }
                }
            }
        )
    }

    result <- tryCatch(
        {
            obj <- run_compute()
            report$ok <- TRUE
            if (is.null(report$selected_backend)) {
                report$selected_backend <- if (identical(backend, "auto")) "cpp" else backend
            }
            attach_report(obj)
        },
        error = function(e) {
            report$error <<- conditionMessage(e)

            if (backend %in% c("auto", "cpp")) {
                fallback <- tryCatch(
                    computeL(scope_obj = scope_obj, ..., backend = "r"),
                    error = function(e2) e2
                )
                if (!inherits(fallback, "error")) {
                    report$ok <<- TRUE
                    report$fallback_used <<- TRUE
                    report$fallback_tier <<- "r"
                    report$selected_backend <<- "r"
                    report$unavailable_tiers <<- unique(c(report$unavailable_tiers, "python"))
                    report$unavailable_tier_notes <<- unique(c(
                        report$unavailable_tier_notes,
                        .lee_l_python_backend_unavailable()
                    ))
                    if (is.null(report$fallback_error)) {
                        report$fallback_error <<- report$error
                    }
                    return(attach_report(fallback))
                }
                report$fallback_used <<- TRUE
                report$fallback_tier <<- "r"
                report$fallback_error <<- conditionMessage(fallback)
            }

            if (identical(on_failure, "return_input")) {
                return(attach_report(scope_obj))
            }
            stop(e)
        }
    )

    result
}

#' Compute L vs R curve.
#' @description
#' Computes a smoothed relationship between Lee's L and Pearson correlation using
#' bootstrapping and stores the fitted curve and confidence intervals into the
#' `scope_object`.
#' @param scope_obj A \code{scope_object} containing Lee's L and correlation matrices.
#' @param grid_name Character; name of the grid layer to operate on.
#' @param level Correlation level used for Pearson r (`grid` or `cell`).
#' @param lee_stats_layer Lee statistics layer name.
#' @param span LOESS span parameter.
#' @param B Number of bootstrap iterations.
#' @param deg Degree for the LOESS fit.
#' @param ncores Number of threads to use.
#' @param length_out Grid size for the fitted curve.
#' @param downsample Downsampling control (ratio < 1, or target count >= 1).
#' @param n_strata Number of strata used when drawing bootstrap samples.
#' @param k_max Maximum number of points used per stratum (guards runtime).
#' @param jitter_eps Optional jitter factor applied to correlation values.
#' @param ci_method Confidence interval method (`percentile`, `basic`, `bc`).
#' @param ci_adjust Optional analytic adjustment for CI width.
#' @param min_rel_width Minimum relative CI width.
#' @param widen_span Additional span widening used when enforcing `min_rel_width`.
#' @param curve_name Name under which the curve is stored.
#' @param backend Backend policy for the Lee's L vs Pearson curve layer. This
#'   layer is currently implemented in R only. `auto` records the package-wide
#'   C++ → R policy and selects R; forcing `python` errors explicitly because no
#'   Python core compute tier is shipped.
#' @param verbose Logical; whether to emit compact public progress messages.
#'   Default is \code{TRUE}. Set to \code{FALSE} for silent operation, or use
#'   \code{options(geneSCOPE.verbose = FALSE)} for global control.
#' @return The modified `scope_object` with stored curve data.
#' @examples
#' \dontrun{
#' scope_obj <- computeLvsRCurve(scope_obj, grid_name = "grid30", level = "cell", downsample = 0.1, ncores = 16)
#' }
#' @seealso `computeL()`, `computeCorrelation()`, `getTopLvsR()`
#' @export
computeLvsRCurve <- function(
    scope_obj,
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
    backend <- .normalize_core_backend(backend, arg = "backend", allow_auto = TRUE)
    .compute_l_vs_r_curve(
        scope_obj = scope_obj,
        grid_name = grid_name,
        level = level,
        lee_stats_layer = lee_stats_layer,
        span = span,
        B = B,
        deg = deg,
        ncores = ncores,
        length_out = length_out,
        downsample = downsample,
        n_strata = n_strata,
        k_max = k_max,
        jitter_eps = jitter_eps,
        ci_method = ci_method,
        ci_adjust = ci_adjust,
        min_rel_width = min_rel_width,
        widen_span = widen_span,
        curve_name = curve_name,
        backend = backend,
        verbose = verbose
    )
}

#' Extract top L vs R pairs.
#' @description
#' Ranks gene pairs by the deviation between Lee's L and Pearson correlation and
#' returns a table of top pairs. The returned Delta-related columns are
#' explicitly marked as exploratory because they compare spatial and
#' non-spatial correlation scales and should be used for heuristic
#' prioritization rather than confirmatory inference.
#' @param scope_obj A `scope_object` containing Lee's L and correlation matrices.
#' @param grid_name Grid layer name.
#' @param pear_level Correlation level (`cell` or `grid`).
#' @param lee_stats_layer Lee statistics layer name.
#' @param expr_layer Optional expression layer name for gene prevalence filtering.
#' @param pear_range Range of Pearson r values to include.
#' @param L_range Range of Lee's L values to include.
#' @param top_n Number of top pairs to return.
#' @param direction Which tail to select (`largest`, `smallest`, `both`).
#' @param do_perm Whether to run permutation testing.
#' @param perms Number of permutations for p-value estimation.
#' @param block_side Block side length used when building permutation blocks.
#' @param use_blocks Whether to use spatial blocks for permutations.
#' @param ncores Number of threads to use.
#' @param clamp_mode Whether to clamp Delta using reference-only or both ends.
#' @param p_adj_mode Multiple-testing adjustment mode.
#' @param mem_limit_GB Memory threshold controlling chunking behavior.
#' @param pval_mode P-value mode used when transforming permutation counts.
#' @param curve_layer Optional curve layer from `computeLvsRCurve()`.
#' @param support_pct_scale Output support-scale contract for the returned gene
#'   coverage diagnostics. Canonical support columns are recorded on `[0,1]`
#'   while legacy percent columns remain available on `0-100` for backward
#'   compatibility.
#' @param backend Backend policy for the ranking layer. This layer is currently
#'   implemented in R only. `auto` records the package-wide C++ → R policy and
#'   selects R; forcing `python` errors explicitly because no Python core compute
#'   tier is shipped.
#' @param CI_rule Confidence-interval filtering rule (`remove_within`, `remove_outside`, `none`).
#' @param verbose Logical; whether to emit compact public progress messages.
#'   Default is \code{TRUE}. Set to \code{FALSE} for silent operation, or use
#'   \code{options(geneSCOPE.verbose = FALSE)} for global control.
#' @return A data.frame with ranked gene pairs and associated statistics,
#' including exploratory `delta`, `delta_method`, and `delta_caveat` columns
#' for heuristic prioritization rather than confirmatory inference.
#' @examples
#' \dontrun{
#' scope_obj <- computeL(scope_obj, grid_name = "grid30", ncores = 16)
#' scope_obj <- computeCorrelation(scope_obj, level = "cell", layer = "logCPM", ncores = 16)
#' top_pairs <- getTopLvsR(scope_obj, grid_name = "grid30", top_n = 200, ncores = 16, L_range = c(0.1, 1))
#' head(top_pairs)
#' }
#' @seealso `computeL()`, `computeCorrelation()`, `computeLvsRCurve()`
#' @export
getTopLvsR <- function(
    scope_obj,
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
    pval_mode = c("beta", "mid", "uniform"),
    curve_layer = NULL,
    support_pct_scale = c("auto", "0-1", "0-100"),
    backend = "cpp",
    CI_rule = c("remove_within", "remove_outside", "none"),
    verbose = TRUE) {
    support_pct_scale <- match.arg(support_pct_scale)
    backend <- .normalize_core_backend(backend, arg = "backend", allow_auto = TRUE)
    .get_top_l_vs_r(
        scope_obj = scope_obj,
        grid_name = grid_name,
        pear_level = pear_level,
        lee_stats_layer = lee_stats_layer,
        expr_layer = expr_layer,
        pear_range = pear_range,
        L_range = L_range,
        top_n = top_n,
        direction = direction,
        do_perm = do_perm,
        perms = perms,
        block_side = block_side,
        use_blocks = use_blocks,
        ncores = ncores,
        clamp_mode = clamp_mode,
        p_adj_mode = p_adj_mode,
        mem_limit_GB = mem_limit_GB,
        pval_mode = pval_mode,
        curve_layer = curve_layer,
        support_pct_scale = support_pct_scale,
        backend = backend,
        CI_rule = CI_rule,
        verbose = verbose
    )
}
