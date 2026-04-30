#' Score module quality using STRINGdb (module coherence)
#'
#' Computes module-level coherence against STRINGdb by delegating to
#' `evaluateModuleQuality_STRING()` and (optionally) writing results back to
#' `scope_obj@stats[[grid_name]][[stats_layer]]$module_quality`.
#'
#' This function supersedes the previous Hotspot-style module-quality scoring:
#' module-quality evaluation in geneSCOPE is now STRINGdb-only.
#'
#' @param x A `scope_object`, an `igraph`, or a list with `graph` and optional `membership`.
#' @param membership Optional membership vector named by gene (module ID per gene).
#' @param membership_cols Optional character vector of `x@meta.data` column names to score.
#'   When provided, `membership` is ignored and the function scores each column separately.
#' @param grid_name Grid layer name for `scope_object` input (default `"grid30"`).
#' @param stats_layer Stats layer name for `scope_object` input (default `"LeeStats_Xz"`).
#' @param graph_slot Graph slot name for `scope_object` input (default `"g_consensus"`).
#' @param edge_weight Edge attribute used for weighting (default `"weight"`).
#' @param species NCBI taxonomy ID for STRING (default 9606).
#' @param string_score_threshold STRING combined score threshold for positives (default 700).
#' @param n_perm Number of membership label shuffles for a within/between delta null (default 200).
#' @param seed Random seed for permutations.
#' @param writeback When TRUE and `x` is a `scope_object`, store results under
#'   `scope_obj@stats[[grid_name]][[stats_layer]]$module_quality$STRINGdb`.
#' @param verbose Emit progress messages when TRUE.
#' @param ... Passed to `evaluateModuleQuality_STRING()` (reserved).
#'
#' @return A list with `summary`, `per_module`, and `meta`. If `writeback=TRUE`
#'   and `x` is a `scope_object`, returns the same list plus `scope_obj`.
#' @export
ScoreModuleQuality <- function(
    x,
    membership = NULL,
    membership_cols = NULL,
    grid_name = "grid30",
    stats_layer = "LeeStats_Xz",
    graph_slot = "g_consensus",
    edge_weight = "weight",
    species = 9606,
    string_score_threshold = 700,
    n_perm = 200,
    seed = 1,
    writeback = FALSE,
    verbose = getOption("geneSCOPE.verbose", TRUE),
    ...
) {
    parent <- "ScoreModuleQuality"

    eval_one <- function(membership_vec, membership_col = NA_character_) {
        if (!is.null(membership_vec) && length(membership_vec)) {
            is_unassigned <- suppressWarnings(as.character(membership_vec)) == "-1"
            membership_vec[is_unassigned] <- NA
        }
        out <- evaluateModuleQuality_STRING(
            x = x,
            membership = membership_vec,
            grid_name = grid_name,
            stats_layer = stats_layer,
            graph_slot = graph_slot,
            edge_weight = edge_weight,
            species = species,
            string_score_threshold = string_score_threshold,
            n_perm = n_perm,
            perm_mode = "label_shuffle",
            seed = seed,
            verbose = verbose,
            ...
        )
        if (is.data.frame(out$per_module) && nrow(out$per_module)) {
            out$per_module$membership_col <- membership_col
        }
        if (is.data.frame(out$summary) && nrow(out$summary)) {
            out$summary$membership_col <- membership_col
        }
        out
    }

    outputs <- list()
    if (!is.null(membership_cols)) {
        if (!inherits(x, "scope_object")) stop("membership_cols requires x to be a scope_object.")
        cols <- as.character(membership_cols)
        cols <- cols[!is.na(cols) & nzchar(cols)]
        if (!length(cols)) stop("membership_cols was provided but empty.")
        missing <- setdiff(cols, colnames(x@meta.data))
        if (length(missing)) stop("membership_cols not found in x@meta.data: ", paste(missing, collapse = ", "))
        for (col in cols) {
            memb <- x@meta.data[[col]]
            names(memb) <- rownames(x@meta.data)
            is_unassigned <- suppressWarnings(as.character(memb)) == "-1"
            memb[is_unassigned] <- NA
            outputs[[col]] <- eval_one(memb, membership_col = col)
        }
    } else {
        outputs[["__single__"]] <- eval_one(membership, membership_col = NA_character_)
    }

    per_module_rows <- lapply(outputs, `[[`, "per_module")
    per_module <- if (length(per_module_rows)) {
        do.call(rbind, Filter(function(x) is.data.frame(x) && nrow(x), per_module_rows))
    } else {
        data.frame()
    }

    summary_rows <- lapply(outputs, `[[`, "summary")
    summary <- if (length(summary_rows)) {
        do.call(rbind, Filter(function(x) is.data.frame(x) && nrow(x), summary_rows))
    } else {
        data.frame()
    }

    meta <- list(
        schema = "STRINGdb_module_quality_v1",
        created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
        grid_name = grid_name,
        stats_layer = stats_layer,
        graph_slot = graph_slot,
        edge_weight = edge_weight,
        species = species,
        string_score_threshold = string_score_threshold,
        n_perm = as.integer(n_perm),
        seed = as.integer(seed),
        membership_cols = if (!is.null(membership_cols)) as.character(membership_cols) else NA_character_
    )

    payload <- list(summary = summary, per_module = per_module, meta = meta)

    if (isTRUE(writeback) && inherits(x, "scope_object")) {
        scope_obj <- x
        if (is.null(scope_obj@stats)) scope_obj@stats <- list()
        if (is.null(scope_obj@stats[[grid_name]])) scope_obj@stats[[grid_name]] <- list()
        if (is.null(scope_obj@stats[[grid_name]][[stats_layer]])) scope_obj@stats[[grid_name]][[stats_layer]] <- list()
        mq <- scope_obj@stats[[grid_name]][[stats_layer]]$module_quality
        if (is.null(mq) || !is.list(mq)) mq <- list()
        if (is.null(mq$STRINGdb) || !is.list(mq$STRINGdb)) mq$STRINGdb <- list()
        run_id <- paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"))
        mq$STRINGdb[[run_id]] <- payload
        scope_obj@stats[[grid_name]][[stats_layer]]$module_quality <- mq
        payload$scope_obj <- scope_obj
        payload$run_id <- run_id
    }

    payload
}

#' Extract stored STRINGdb module-quality scores.
#'
#' Reads `module_quality$STRINGdb` payloads from `scope_obj@stats` and returns a
#' flattened table.
#'
#' @param x A `scope_object`.
#' @param grid_name Optional grid filter.
#' @param stats_layer Optional stats layer filter.
#' @return A data.frame with one row per module per run (or empty if not found).
#' @export
GetModuleQualityScores <- function(x, grid_name = NULL, stats_layer = NULL) {
    if (!inherits(x, "scope_object")) stop("GetModuleQualityScores requires a scope_object.")
    stats_root <- x@stats
    rows <- list()
    for (gname in names(stats_root)) {
        if (!is.null(grid_name) && !identical(gname, grid_name)) next
        layer_root <- stats_root[[gname]]
        for (lname in names(layer_root)) {
            if (!is.null(stats_layer) && !identical(lname, stats_layer)) next
            entry <- layer_root[[lname]]$module_quality
            if (is.null(entry) || is.null(entry$STRINGdb) || !length(entry$STRINGdb)) next
            for (run_id in names(entry$STRINGdb)) {
                payload <- entry$STRINGdb[[run_id]]
                if (is.null(payload$per_module) || !is.data.frame(payload$per_module)) next
                df <- payload$per_module
                df$grid_name <- gname
                df$stats_layer <- lname
                df$run_id <- run_id
                rows[[length(rows) + 1]] <- df
            }
        }
    }
    if (!length(rows)) return(data.frame())
    do.call(rbind, rows)
}

#' Summarize stored STRINGdb module-quality scores.
#'
#' @param x A `scope_object`.
#' @param grid_name Optional grid filter.
#' @param stats_layer Optional stats layer filter.
#' @return A list with `per_module`, `per_run`, and `summary_text`.
#' @export
SummarizeModuleQualityScores <- function(x, grid_name = NULL, stats_layer = NULL) {
    per_module <- GetModuleQualityScores(x = x, grid_name = grid_name, stats_layer = stats_layer)
    per_run <- data.frame()
    if (nrow(per_module)) {
        groups <- interaction(per_module$grid_name, per_module$stats_layer, per_module$run_id, drop = TRUE)
        rows <- lapply(split(per_module, groups), function(df) {
            data.frame(
                grid_name = df$grid_name[1],
                stats_layer = df$stats_layer[1],
                run_id = df$run_id[1],
                n_modules = nrow(df),
                mean_within_string = mean(df$within_string_mean, na.rm = TRUE),
                mean_between_string = mean(df$between_string_mean, na.rm = TRUE),
                mean_delta_within_between = mean(df$delta_within_between, na.rm = TRUE),
                stringsAsFactors = FALSE
            )
        })
        per_run <- do.call(rbind, rows)
    }
    summary_text <- paste0("Summarized STRINGdb module quality for ", nrow(per_run), " run(s).")
    list(per_module = per_module, per_run = per_run, summary_text = summary_text)
}
