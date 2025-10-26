##' computeL_visium: Fast (≤1 hour) Visium-oriented Lee's L sparse network
##'
##' Goal: For G≈40k genes and N≈500–2000 spots, produce a high-confidence
##' sparse network (Top-K edges per gene) and multiple-testing results within
##' about an hour, without forming the full G×G matrix.
##'
##' Pipeline:
##' 1) Gene prescreening (Iδ + detection/abundance) to keep S≈several thousand;
##' 2) Sparse–dense multiply WX = W %*% X subset (C++ accelerated);
##' 3) Candidate generation via sign random projection LSH (C++ signatures,
##'    R-side bucketing);
##' 4) Exact refinement only on candidates to compute L = XᵀWX (C++
##'    per-row Top-K on candidate CSR);
##' 5) Analytical Z, p, BH/FDR; write results into stats layer LeeStats_visium.
##'
##' @param scope_obj A scope_object.
##' @param grid_name Character. If NULL and only one layer exists, it is chosen.
##' @param norm_layer Character. Normalized layer name (default "Xz"). Rows=spots,
##'        cols=genes.
##' @param use_idelta Logical. If Iδ is missing, call computeIDelta() (default TRUE).
##' @param S_target Integer. Target number of genes after prescreening (default 8000).
##' @param min_detect Numeric. Baseline detection-rate threshold for the first
##'        prescreen step (default 0.10).
##' @param winsor_high Numeric. Upper winsorization quantile for Iδ (default 0.95).
##' @param lsh_bits Integer. Bits per hash table b (default 12; avg bucket size ~ S/2^b).
##' @param lsh_tables Integer. Number of hash tables L (default 6).
##' @param K_candidate Integer. Max candidates per gene after merging LSH buckets
##'        (default 600).
##' @param K_keep Integer. Final Top-K edges per gene to keep (default 100).
##' @param ncores Integer. Number of threads (default: safe auto).
##' @param chunk_size Integer. Column-block size for WX etc. (default 1024L).
##' @param precision Character. "float32" or "float64" (default "float32").
##' @param mode Character. "network": output sparse network only;
##'        "full_if_small": if S≤3000, switch to exact full matrix (upper triangle).
##' @param calibrate Logical. Enable global permutation tail calibration (placeholder,
##'        default FALSE).
##' @param verbose Logical. Print progress (default TRUE).
##' @return The modified scope_obj (invisible).
##' @export
computeL_visium <- function(
    scope_obj,
    grid_name = NULL,
    norm_layer = "Xz",
    use_idelta = TRUE,
    S_target = 8000L,
    min_detect = 0.10,
    winsor_high = 0.95,
    lsh_bits = 12L,
    lsh_tables = 6L,
    K_candidate = 600L,
    K_keep = 100L,
    ncores = NULL,
    chunk_size = 1024L,
    precision = c("float32", "float64"),
    mode = c("network", "full_if_small"),
    calibrate = FALSE,
    verbose = getOption("geneSCOPE.verbose", TRUE)) {

    precision <- match.arg(precision)
    mode <- match.arg(mode)

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
    W <- g_layer$W  # dgCMatrix
    S0 <- sum(W)

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

    Xs <- X[, keep, drop = FALSE]        # n × S
    colnames(Xs) <- keep

    ## ---- 3. Sparse–dense multiply WX = W %*% Xs (C++ accelerated) ----
    if (verbose) message("[computeL_visium] Computing WX (sparse × dense)…")
    WX <- spmm_dgc_dense(Xs, W, n_threads = ncores, precision = precision)

    ## ---- 4. LSH candidate generation (sign RP; C++ signatures; R bucketing) ----
    if (verbose) message("[computeL_visium] Building LSH candidates (bits=", lsh_bits, ", tables=", lsh_tables, ")…")
    # Compute signatures on Xs (more stable); could switch to WX if desired
    sig <- rp_sign_bits(Xs, bits = lsh_bits, n_tables = lsh_tables, seed = 20251026L, n_threads = ncores)
    cand <- .build_lsh_candidates_from_sign(sig, K_cap = as.integer(K_candidate))

    ## If candidates are too few and S small, fall back to dense neighbors
    if (mode == "full_if_small" && S <= 3000L) {
        if (verbose) message("[computeL_visium] S≤3000: fallback to exact upper-triangular block Top-K (skip LSH)…")
        cand <- .build_dense_candidates(S)
    }

    ## ---- 5. Exact L on candidate CSR and per-row Top-K (C++ accelerated) ----
    if (verbose) message("[computeL_visium] Refining candidate edges (L) and keeping per-row Top-", K_keep, "…")
    topk <- leeL_topk_candidates(Xs, WX, cand$row_ptr, cand$indices, K_keep = as.integer(K_keep), n_threads = ncores)

    ## ---- 6. Z, p, FDR and edge table ----
    EZ <- -1 / (n - 1)
    Var <- (n^2 * (n - 2)) / ((n - 1)^2 * (n - 3) * S0)
    Zvals <- (topk$values - EZ) / sqrt(Var)
    pvals <- 2 * stats::pnorm(-abs(Zvals))
    qvals <- stats::p.adjust(pvals, method = "BH")

    # Expand to edge table (deduplicate symmetric pairs; keep max)
    ii <- rep.int(seq_len(S), diff(topk$row_ptr))
    jj <- topk$indices
    g_i <- keep[ii]; g_j <- keep[jj]
    # Undirected key
    key <- ifelse(g_i < g_j, paste0(g_i, "|", g_j), paste0(g_j, "|", g_i))
    o <- order(key, -topk$values)
    key <- key[o]
    df <- data.frame(
        from = g_i[o], to = g_j[o], L = topk$values[o],
        Z = Zvals[o], p = pvals[o], q = qvals[o],
        stringsAsFactors = FALSE
    )
    # Keep the strongest entry per pair
    keep_idx <- !duplicated(key)
    edges <- df[keep_idx, ]

    ## ---- 7. Write to stats layer ----
    layer_name <- paste0("LeeStats_visium_", norm_layer)
    if (is.null(scope_obj@stats)) scope_obj@stats <- list()
    if (is.null(scope_obj@stats[[grid_name]])) scope_obj@stats[[grid_name]] <- list()
    scope_obj@stats[[grid_name]][[layer_name]] <- list(
        edges = edges,
        meta = list(
            n = n, G = G, S = S,
            W_nonzero = length(W@x), S0 = S0,
            S_target = S_target,
            min_detect = min_detect,
            winsor_high = winsor_high,
            lsh_bits = lsh_bits, lsh_tables = lsh_tables,
            K_candidate = K_candidate, K_keep = K_keep,
            precision = precision,
            mode = mode,
            calibrate = calibrate,
            time = Sys.time()
        )
    )

    if (verbose) {
        message("[computeL_visium] Done: edges=", nrow(edges), 
                ", edges with FDR<0.05=", sum(edges$q < 0.05, na.rm = TRUE))
    }
    invisible(scope_obj)
}

## ---- Internal helper: build candidate CSR from signature matrix ----
.build_lsh_candidates_from_sign <- function(sig, K_cap = 600L) {
    # sig: S × L character/integer; each column is a 64-bit bucket ID per table
    S <- nrow(sig); L <- ncol(sig)
    # Collect candidate sets per row
    cand_list <- vector("list", S)
    for (t in seq_len(L)) {
        keys <- sig[, t]
        ord <- order(keys)
        k_sorted <- keys[ord]
        # Linear scan equal-key segments
        start <- 1L
        while (start <= S) {
            end <- start
            while (end < S && k_sorted[end + 1L] == k_sorted[start]) end <- end + 1L
            if ((end - start + 1L) >= 2L) {
                idx <- ord[start:end]
                for (u in idx) {
                    cand_list[[u]] <- c(cand_list[[u]], setdiff(idx, u))
                }
            }
            start <- end + 1L
        }
    }
    # Deduplicate and cap
    for (i in seq_len(S)) {
        if (length(cand_list[[i]]) == 0L) next
        cand_list[[i]] <- unique(cand_list[[i]])
        if (length(cand_list[[i]]) > K_cap) {
            cand_list[[i]] <- cand_list[[i]][seq_len(K_cap)]
        }
    }
    # Convert to CSR
    nnz <- sum(vapply(cand_list, length, integer(1)))
    row_ptr <- integer(S + 1L)
    indices <- integer(nnz)
    pos <- 1L
    for (i in seq_len(S)) {
        row_ptr[i] <- pos
        ci <- cand_list[[i]]
        if (length(ci)) {
            k <- length(ci)
            indices[pos:(pos + k - 1L)] <- ci
            pos <- pos + k
        }
    }
    row_ptr[S + 1L] <- pos
    list(row_ptr = row_ptr, indices = indices)
}

## ---- Internal helper: dense candidates (upper triangle) for small S ----
.build_dense_candidates <- function(S) {
    # i<j full upper triangle CSR
    lens <- (S - 1L):0L
    row_ptr <- c(1L, 1L + cumsum(lens))
    nnz <- row_ptr[length(row_ptr)] - 1L
    indices <- integer(nnz)
    pos <- 1L
    for (i in seq_len(S)) {
        if (i < S) {
            jidx <- (i + 1L):S
            k <- length(jidx)
            indices[pos:(pos + k - 1L)] <- jidx
            pos <- pos + k
        }
    }
    list(row_ptr = row_ptr, indices = indices)
}
