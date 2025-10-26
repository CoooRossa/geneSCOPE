#' computeL_visium: Visium-friendly wrapper around computeL with light prescreen
#'
#' Goal: For G≈40k genes and N≈500–2000 spots, first perform a fast gene
#' prescreen using detection rate and Morisita's Iδ to select S≈several thousand
#' informative genes, then delegate to the classic computeL() so that
#' permutations, statistics, and the output structure are exactly identical to
#' the legacy algorithm.
#'
#' @param scope_obj A scope_object.
#' @param grid_name Character. If NULL and only one layer exists, it is chosen.
#' @param norm_layer Character. Normalized layer name (default "Xz"). Rows=spots,
#'        cols=genes.
#' @param use_idelta Logical. If Iδ is missing, call computeIDelta() (default TRUE).
#' @param S_target Integer. Target number of genes after prescreening (default 8000).
#' @param min_detect Numeric. Baseline detection-rate threshold for the first
#'        prescreen step (default 0.10).
#' @param winsor_high Numeric. Upper winsorization quantile for Iδ (default 0.95).
#' @param ncores Integer. Thread hint for prescreen only (computeL has its own
#'        threading); default uses getSafeThreadCount().
#' @param ... Additional arguments forwarded verbatim to computeL(), e.g.
#'        perms, block_size, L_min, lee_stats_layer_name, ncores, etc.
#' @param verbose Logical. Print progress (default TRUE).
#' @return The modified scope_obj (invisible).
#' @export
computeL_visium <- function(
    scope_obj,
    grid_name = NULL,
    norm_layer = "Xz",
    use_idelta = TRUE,
    S_target = 8000L,
    min_detect = 0.10,
    winsor_high = 0.95,
    ncores = NULL,
    verbose = getOption("geneSCOPE.verbose", TRUE),
    ...) {

    ## ---- 0. Thread count and grid layer selection ----
    if (is.null(ncores)) {
        ncores <- getSafeThreadCount(default = 8L)
    }
    g_layer <- .selectGridLayer(scope_obj, grid_name)
    grid_name <- if (is.null(grid_name)) {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    } else grid_name

    if (is.null(g_layer[[norm_layer]])) {
        stop("computeL_visium: normalized layer '", norm_layer, "' not found.")
    }
    X <- g_layer[[norm_layer]]  # n × g
    if (!is.matrix(X)) X <- as.matrix(X)
    n <- nrow(X); G <- ncol(X)
    if (n < 2L || G < 2L) stop("Insufficient data size to compute Lee's L.")

    ## ---- 1. Check/compute spatial weights W ----
    if (is.null(g_layer$W) || sum(g_layer$W) == 0) {
        if (verbose) message("[computeL_visium] W not found; auto computeWeights(style='B')…")
        scope_obj <- computeWeights(scope_obj, grid_name = grid_name, style = "B", store_mat = TRUE)
        g_layer <- .selectGridLayer(scope_obj, grid_name)
    }
    W <- g_layer$W  # dgCMatrix (ensures W exists for computeL later)

    ## ---- 2. Iδ and gene prescreening ----
    if (use_idelta) {
        meta_col <- paste0(grid_name, "_iDelta")
        if (is.null(scope_obj@meta.data) || !(meta_col %in% colnames(scope_obj@meta.data))) {
            if (verbose) message("[computeL_visium] Computing Iδ…")
            scope_obj <- computeIDelta(scope_obj, grid_name = grid_name, level = "grid", ncore = min(8L, ncores), verbose = verbose)
        }
    }
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
            dtt <- data.table::as.data.table(dt)
            n_spots <- length(g_layer$grid_info$grid_id)
            dr <- dtt[, .N, by = gene]
            detect_rate <- dr$N / n_spots; names(detect_rate) <- dr$gene
            tc <- dtt[, sum(count), by = gene]
            total_count <- tc$V1; names(total_count) <- tc$gene
        } else {
            n_spots <- length(g_layer$grid_info$grid_id)
            det_list <- split(dt$grid_id, dt$gene)
            detect_rate <- vapply(det_list, function(v) length(unique(v)) / n_spots, numeric(1))
            cnt_list <- split(dt$count, dt$gene)
            total_count <- vapply(cnt_list, sum, numeric(1))
        }
    } else {
        # Fallback: infer from X (may be biased)
        detect_rate <- colMeans(X > 0)
        total_count <- colSums(pmax(X, 0))
        names(detect_rate) <- names(total_count) <- colnames(X)
    }

    # Iδ fetch and winsorization
    idelta_vec <- rep(NA_real_, G); names(idelta_vec) <- genes_all
    meta_col <- paste0(grid_name, "_iDelta")
    if (!is.null(scope_obj@meta.data) && (meta_col %in% colnames(scope_obj@meta.data))) {
        idelta_vec[rownames(scope_obj@meta.data)] <- scope_obj@meta.data[, meta_col]
    }
    idelta_vec[is.infinite(idelta_vec)] <- NA_real_
    idelta_vec_w <- idelta_vec
    if (all(is.na(idelta_vec_w))) {
        if (verbose) message("[computeL_visium] Iδ not found; skip Iδ-based ranking and use detection rate only.")
        idelta_vec_w <- rep(0, G); names(idelta_vec_w) <- genes_all
    } else {
        qh <- stats::quantile(idelta_vec_w, probs = winsor_high, na.rm = TRUE)
        idelta_vec_w <- pmin(idelta_vec_w, qh)
    }

    # 2.2 Prescreen: rank by Iδ within detection-rate filter; keep top S_target
    keep0 <- intersect(names(detect_rate)[detect_rate >= min_detect], genes_all)
    ord <- order(idelta_vec_w[keep0], decreasing = TRUE, na.last = NA)
    keep <- keep0[ord]
    if (length(keep) > S_target) keep <- keep[seq_len(S_target)]
    S <- length(keep)
    if (S < 100L) stop("Too few genes after prescreen (<100); lower min_detect or increase S_target.")
    if (verbose) message("[computeL_visium] Prescreened genes S=", S, " (target ", S_target, ")")

    ## ---- 3. Delegate to computeL with the filtered gene set ----
    if (verbose) message("[computeL_visium] Delegating to computeL() with prescreened genes …")
    scope_obj <- computeL(
        scope_obj,
        grid_name = grid_name,
        norm_layer = norm_layer,
        genes = keep,
        verbose = verbose,
        ...
    )

    invisible(scope_obj)
}
