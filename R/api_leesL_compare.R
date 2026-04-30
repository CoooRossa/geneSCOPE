#' Compare Lee's L across samples
#'
#' @description
#' Compare Lee's L across samples without redefining the Lee's L statistic.
#' This function preserves the original Lee's L (and delta_Lr = L - Pearson r)
#' and performs cross-sample comparison through:
#' (1) a comparability gate based on support/overlap and estimation stability,
#' (2) precision-weighted standardized differences, and
#' (3) sample-level L-network background summaries.
#'
#' It does not modify existing `computeL()`, `lee_L()`, `getTopLvsR()`, or
#' `computeLvsRCurve()` behavior, and it is not a global rank/min-max normalization.
#'
#' @param x Either a list of `scope_object` samples, a list of pair-level
#'   Lee's L tables, a single `scope_object`, or a single pair-level table.
#'   Pair-level tables should include `gene_i`, `gene_j`, `L`, and optionally
#'   `Pearson_r`, `delta_Lr`, uncertainty, and support fields.
#' @param sample_ids Optional character vector of sample identifiers. When `x`
#'   is a list, this is used to label samples if names are not present.
#' @param grid_name Character; name of the grid layer to operate on. If \code{NULL}
#'   and only one grid layer exists, it is auto-selected. If multiple layers exist,
#'   the name must be specified.
#' @param lee_stats_layer Character; name of the Lee statistics layer (e.g., \code{"LeeStats_Xz"}).
#' @param pear_level Pearson correlation level (`cell` or `grid`).
#' @param pairs Optional data.frame/matrix with gene pair columns used to
#'   restrict comparisons. If omitted, all available pairs are used.
#' @param sample_info Optional data.frame with `sample_id` and optional
#'   `group`, `batch`, `donor`, and L-network QC metrics.
#' @param mode Comparison mode: `pairwise`, `group`, or `both`.
#' @param group_col Column name in `sample_info` (or pair tables) used for
#'   group comparisons.
#' @param min_support_n Minimum per-gene support count to pass the gate.
#' @param min_support_pct Minimum per-gene support percentage to pass the gate.
#' @param support_pct_range Declared support-percentage scale. geneSCOPE now
#'   keeps a canonical internal gate scale on `[0,1]`, while preserving legacy
#'   `0-100` percent columns in the returned table for backward compatibility.
#'   `auto` detects percent-like inputs and normalizes them to the canonical
#'   support scale before gating.
#' @param support_pct_scale Deprecated alias for `support_pct_range`.
#' @param backend Deprecated compatibility argument. This comparison layer is
#'   currently implemented in R only and therefore uses the R fallback under
#'   the package-wide C++->R policy. `python` errors explicitly.
#' @param min_overlap_n Minimum pair overlap count to pass the gate.
#' @param min_overlap_ratio Minimum pair overlap ratio to pass the gate.
#' @param max_rel_se Maximum relative SE (SE / |L|) allowed by the gate.
#' @param max_ci_width Maximum CI width allowed by the gate (if CI available).
#' @param require_uncertainty When TRUE, pairs without uncertainty estimates
#'   are marked non-comparable.
#' @param min_samples_per_group Minimum samples per group for group comparison.
#' @param p_adj_mode P-value adjustment for standardized differences.
#' @param verbose Logical; whether to emit compact public progress messages.
#'   Default is \code{TRUE}. Set to \code{FALSE} for silent operation, or use
#'   \code{options(geneSCOPE.verbose = FALSE)} for global control.
#' @return A list with:
#'   - `pair_table`: standardized per-sample pair table with gate results
#'   - `pairwise`: sample-to-sample comparison results (if requested)
#'   - `group`: group comparison results (if requested)
#'   - `background`: sample-level L-network background summary
#'   - `diagnostics`: parameters and filtering diagnostics
#' @export
compareLeesL <- function(
    x,
    sample_ids = NULL,
    grid_name = NULL,
    lee_stats_layer = NULL,
    pear_level = c("cell", "grid"),
    pairs = NULL,
    sample_info = NULL,
    mode = c("pairwise", "group", "both"),
    group_col = "group",
    min_support_n = 5,
    min_support_pct = 0,
    support_pct_range = c("auto", "0-100", "0-1"),
    support_pct_scale = NULL,
    backend = "cpp",
    min_overlap_n = 3,
    min_overlap_ratio = 0.1,
    max_rel_se = Inf,
    max_ci_width = Inf,
    require_uncertainty = TRUE,
    min_samples_per_group = 1,
    p_adj_mode = c("BH", "none"),
    verbose = getOption("geneSCOPE.verbose", TRUE)
) {
    mode <- match.arg(mode)
    p_adj_mode <- match.arg(p_adj_mode)
    pear_level <- match.arg(pear_level)
    verbose <- .resolve_runtime_verbose(verbose)

    .with_log_session("compareLeesL", verbose = verbose, {
        .log_start("compareLeesL", verbose = verbose)
        
        if (!is.null(support_pct_scale)) {
            support_pct_scale <- as.character(support_pct_scale)[1]
            if (!support_pct_scale %in% c("auto", "0-100", "0-1")) {
                stop("support_pct_scale must be one of 'auto', '0-100', or '0-1'")
            }
            if (!missing(support_pct_range)) {
                support_pct_range_arg <- as.character(support_pct_range)[1]
                if (!identical(support_pct_range_arg, support_pct_scale)) {
                    stop("support_pct_scale and support_pct_range disagree; use only support_pct_range or provide matching values")
                }
            } else {
                support_pct_range <- support_pct_scale
            }
            warning("support_pct_scale is deprecated; use support_pct_range instead.", call. = FALSE)
        }
        support_pct_range <- match.arg(support_pct_range)
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
                "[geneSCOPE::compareLeesL] backend='", backend,
                "' is not implemented for this R-only comparison statistic. ",
                "Use the default C++->R fallback policy or backend='r'.",
                call. = FALSE
            )
        } else if (backend %in% c("auto", "cpp") && .log_enabled(verbose)) {
            .log_info("compareLeesL", "S01", "Comparison statistic is R-only; using R fallback under C++->R policy.", verbose)
        }

        n_inputs <- if (.is_scope_object(x) || is.data.frame(x) || is.matrix(x)) 1L else length(x)
        .log_key_values("compareLeesL", list(
            mode = mode,
            backend = backend,
            inputs = n_inputs,
            support_pct_range = support_pct_range,
            min_support_pct = min_support_pct
        ))

        pairs_df <- .leesl_compare_standardize_pairs(pairs)

        # Normalize input to list
        input_list <- list()
        if (.is_scope_object(x)) {
            input_list <- list(x)
        } else if (is.data.frame(x) || is.matrix(x)) {
            input_list <- list(as.data.frame(x, stringsAsFactors = FALSE))
        } else if (is.list(x)) {
            input_list <- x
        } else {
            stop("x must be a scope_object, data.frame, or list")
        }

        # Resolve sample IDs
        if (is.null(sample_ids) && length(input_list) > 1) {
            if (!is.null(names(input_list)) && all(nzchar(names(input_list)))) {
                sample_ids <- names(input_list)
            }
        }

        pair_tables <- list()
        sample_meta_list <- list()
        sample_qc_list <- list()

        add_sample_entry <- function(df, sid) {
            pair_tables[[length(pair_tables) + 1]] <<- df
            sample_meta_list[[length(sample_meta_list) + 1]] <<- data.frame(sample_id = sid, stringsAsFactors = FALSE)
            sample_qc_list[[length(sample_qc_list) + 1]] <<- data.frame(sample_id = sid, stringsAsFactors = FALSE)
        }

    for (i in seq_along(input_list)) {
        entry <- input_list[[i]]
        sid <- if (!is.null(sample_ids) && length(sample_ids) >= i) {
            as.character(sample_ids[i])
        } else if (!is.null(names(input_list)) && nzchar(names(input_list)[i])) {
            as.character(names(input_list)[i])
        } else {
            paste0("sample", i)
        }

        if (.is_scope_object(entry)) {
            res <- .leesl_compare_extract_from_scope(
                scope_obj = entry,
                sample_id = sid,
                grid_name = grid_name,
                lee_stats_layer = lee_stats_layer,
                pear_level = pear_level,
                pairs = pairs_df,
                verbose = verbose
            )
            pair_tables[[length(pair_tables) + 1]] <- res$pairs
            sample_meta_list[[length(sample_meta_list) + 1]] <- res$sample_meta
            sample_qc_list[[length(sample_qc_list) + 1]] <- res$sample_qc
        } else if (is.data.frame(entry) || is.matrix(entry)) {
            df <- as.data.frame(entry, stringsAsFactors = FALSE)
            if ("sample_id" %in% names(df)) {
                split_df <- split(df, df$sample_id)
                for (sid2 in names(split_df)) {
                    tbl <- .leesl_compare_standardize_pair_table(split_df[[sid2]], sample_id = sid2, sample_info = sample_info)
                    pair_tables[[length(pair_tables) + 1]] <- tbl
                    sample_meta_list[[length(sample_meta_list) + 1]] <- data.frame(sample_id = sid2, stringsAsFactors = FALSE)
                    sample_qc_list[[length(sample_qc_list) + 1]] <- data.frame(sample_id = sid2, stringsAsFactors = FALSE)
                }
            } else {
                tbl <- .leesl_compare_standardize_pair_table(df, sample_id = sid, sample_info = sample_info)
                add_sample_entry(tbl, sid)
            }
        } else {
            stop("Unsupported element in x at index ", i)
        }
    }

        if (!length(pair_tables)) stop("No input samples detected")

    # Build reference pair list
    if (is.null(pairs_df) || !nrow(pairs_df)) {
        ref_pairs <- do.call(rbind, lapply(pair_tables, function(df) {
            unique(df[, c("gene_i", "gene_j", "pair_id"), drop = FALSE])
        }))
        if (!nrow(ref_pairs)) stop("No gene pairs found in inputs")
        ref_pairs <- unique(ref_pairs)
    } else {
        ref_pairs <- pairs_df
    }

    # Expand to full pair set per sample
    expanded <- list()
    for (tbl in pair_tables) {
        sid <- unique(tbl$sample_id)
        sid <- sid[!is.na(sid)][1]
        expanded[[length(expanded) + 1]] <- .leesl_compare_expand_pairs(tbl, ref_pairs, sid)
    }
    pair_table <- do.call(rbind, expanded)

    # Attach sample-level info
        pair_table <- .leesl_compare_attach_sample_info(pair_table, sample_info)

    # Merge sample metadata and infer uncertainty
    sample_meta_all <- if (length(sample_meta_list)) do.call(rbind, sample_meta_list) else data.frame()
    sample_qc_all <- if (length(sample_qc_list)) do.call(rbind, sample_qc_list) else data.frame()

    pair_table <- .leesl_compare_infer_uncertainty(pair_table, sample_meta = sample_meta_all)
    pair_table <- .leesl_compare_normalize_support_pct(
        pair_table,
        support_pct_range = support_pct_range,
        caller = "compareLeesL()"
    )
    min_support_prop <- .normalize_support_pct_threshold(
        min_support_pct,
        support_pct_scale = support_pct_range,
        detected_scale = attr(pair_table, "support_pct_input_scale"),
        caller = "compareLeesL()"
    )

    # Apply comparability gate
    pair_table <- .leesl_compare_apply_gate(
        pair_table,
        min_support_n = min_support_n,
        min_support_prop = min_support_prop,
        min_overlap_n = min_overlap_n,
        min_overlap_ratio = min_overlap_ratio,
        max_rel_se = max_rel_se,
        max_ci_width = max_ci_width,
        require_uncertainty = require_uncertainty
    )

    # Build background summary
        background <- .leesl_compare_build_background(sample_qc_all, sample_info = sample_info)

    # Pairwise comparisons
    pairwise <- data.frame()
    if (mode %in% c("pairwise", "both")) {
        pairwise <- .leesl_compare_pairwise(
            pair_table,
            sample_info = sample_info,
            background = background,
            p_adj_mode = p_adj_mode
        )
    }

    # Group comparisons
    group_res <- data.frame()
    if (mode %in% c("group", "both") && group_col %in% names(pair_table)) {
        group_res <- .leesl_compare_group(
            pair_table,
            group_col = group_col,
            min_samples_per_group = min_samples_per_group,
            p_adj_mode = p_adj_mode
        )
    }

    # Diagnostics
    gate_summary <- aggregate(
        cbind(gate_L = pair_table$gate_L, gate_delta = pair_table$gate_delta),
        by = list(sample_id = pair_table$sample_id),
        FUN = function(v) sum(v, na.rm = TRUE)
    )
    total_pairs <- aggregate(pair_id ~ sample_id, data = pair_table, FUN = length)
    names(total_pairs)[2] <- "n_pairs"
    gate_summary <- merge(gate_summary, total_pairs, by = "sample_id", all.x = TRUE, sort = FALSE)

    diagnostics <- list(
        params = list(
            mode = mode,
            group_col = group_col,
            min_support_n = min_support_n,
            min_support_pct = min_support_pct,
            min_support_prop = min_support_prop,
            support_pct_range = support_pct_range,
            support_pct_scale = support_pct_range,
            min_overlap_n = min_overlap_n,
            min_overlap_ratio = min_overlap_ratio,
            max_rel_se = max_rel_se,
            max_ci_width = max_ci_width,
            require_uncertainty = require_uncertainty,
            min_samples_per_group = min_samples_per_group,
            p_adj_mode = p_adj_mode
        ),
        backend = list(
            requested = backend,
            selected = "r",
            policy = backend_policy
        ),
        uncertainty_sources = list(
            L_se_source = table(pair_table$L_se_source, useNA = "ifany"),
            delta_se_source = table(pair_table$delta_se_source, useNA = "ifany")
        ),
        gate_summary = gate_summary,
        gate_reason_L = table(pair_table$gate_reason_L, useNA = "ifany"),
        gate_reason_delta = table(pair_table$gate_reason_delta, useNA = "ifany"),
        background_enabled = nrow(background) > 0,
        support_pct_scale = list(
            requested = attr(pair_table, "support_pct_requested_scale"),
            input = attr(pair_table, "support_pct_input_scale"),
            used = attr(pair_table, "support_pct_scale_used"),
            normalized = attr(pair_table, "support_pct_normalized_scale"),
            internal = attr(pair_table, "support_pct_internal_scale"),
            legacy = attr(pair_table, "support_pct_legacy_scale")
        ),
        support_pct_auto_normalized = attr(pair_table, "support_pct_auto_normalized"),
        support_pct_invalid_values = attr(pair_table, "support_pct_invalid_values")
    )

    result <- list(
        pair_table = pair_table,
        pairwise = pairwise,
        group = group_res,
        background = background,
        diagnostics = diagnostics
    )

    .log_key_values("compareLeesL", list(
        samples = length(unique(pair_table$sample_id)),
        pairs = nrow(pair_table),
        gate_L_pass = sum(pair_table$gate_L, na.rm = TRUE),
        gate_delta_pass = sum(pair_table$gate_delta, na.rm = TRUE),
        pairwise_rows = if (is.data.frame(pairwise)) nrow(pairwise) else 0L,
        group_rows = if (is.data.frame(group_res)) nrow(group_res) else 0L,
        support_scale_used = diagnostics$support_pct_scale$used
    ))

        .log_done("compareLeesL", verbose = verbose)
        result
    })
}
