#' Resolve gene indices for Lee's L computations.
#' @description
#' Internal helper for `.compute_lee_l()`. Filters missing/empty gene names and
#' drops unknown genes with a warning.
#' @param genes Character vector of requested genes (or `NULL` for all genes).
#' @param all_genes Character vector of available gene names (colnames of Xz).
#' @return Integer vector of indices into `all_genes`.
#' @keywords internal
.resolve_leeL_gene_indices <- function(genes, all_genes) {
    stopifnot(is.character(all_genes), length(all_genes) > 0L)

    if (is.null(genes)) {
        return(seq_along(all_genes))
    }

    genes <- as.character(genes)
    genes <- genes[!is.na(genes)]
    genes <- trimws(genes)
    genes <- genes[nzchar(genes)]

    if (!length(genes)) {
        stop("`genes` was provided but is empty after filtering.", call. = FALSE)
    }

    genes <- unique(genes)
    idx <- match(genes, all_genes)
    missing <- genes[is.na(idx)]

    if (length(missing)) {
        warning(
            paste0(
                "Dropping ", length(missing),
                " requested genes not found in norm_layer columns. Examples: ",
                paste(utils::head(missing, 10L), collapse = ", "),
                if (length(missing) > 10L) " ..." else ""
            ),
            call. = FALSE
        )
    }

    idx <- idx[!is.na(idx)]
    if (!length(idx)) {
        stop("None of the requested genes were found in norm_layer column names.", call. = FALSE)
    }

    idx
}

.resolve_lee_l_backend <- function(preferred = getOption("geneSCOPE.lee_l.backend", "cpp"),
                                    respect_env = TRUE) {
    if (isTRUE(respect_env)) {
        env_backend <- Sys.getenv("GENESCOPE_LEE_L_BACKEND", "")
        if (nzchar(env_backend)) {
            preferred <- env_backend
        }
    }
    .normalize_core_backend(preferred, arg = "backend", allow_auto = TRUE)
}

.with_lee_l_backend <- function(backend, expr) {
    backend <- .resolve_lee_l_backend(backend, respect_env = FALSE)
    if (identical(backend, "python")) {
        stop(.lee_l_python_backend_unavailable(), call. = FALSE)
    }

    old_backend <- getOption("geneSCOPE.lee_l.backend")
    old_env <- Sys.getenv("GENESCOPE_LEE_L_BACKEND", unset = NA_character_)

    options(geneSCOPE.lee_l.backend = backend)
    Sys.setenv(GENESCOPE_LEE_L_BACKEND = backend)
    on.exit({
        options(geneSCOPE.lee_l.backend = old_backend)
        if (is.na(old_env)) {
            Sys.unsetenv("GENESCOPE_LEE_L_BACKEND")
        } else {
            Sys.setenv(GENESCOPE_LEE_L_BACKEND = old_env)
        }
    }, add = TRUE)

    force(expr)
}

.lee_l_python_backend_unavailable <- function() {
    .core_backend_python_unavailable_message("computeL")
}

.lee_l_auto_fallback_warning <- function(error,
                                         chunked = FALSE) {
    msg <- conditionMessage(error)
    if (grepl("disabled on Darwin", msg, fixed = TRUE)) {
        return(invisible(NULL))
    }
    prefix <- if (isTRUE(chunked)) {
        "Native Lee's L chunk backend failed"
    } else {
        "Native Lee's L backend failed"
    }
    warning(
        prefix,
        "; Python Lee's L backend is not implemented, so falling back to the reference R backend: ",
        msg,
        call. = FALSE,
        immediate. = TRUE
    )
}

.validate_lee_l_native_inputs <- function(Xz, W, caller = ".compute_lee_l") {
    Xz <- as.matrix(Xz)
    W <- methods::as(W, "dgCMatrix")

    if (!nrow(Xz) || !ncol(Xz)) {
        stop(caller, ": Xz must contain at least one row and one column.", call. = FALSE)
    }
    if (!is.numeric(Xz)) {
        stop(caller, ": Xz must be numeric.", call. = FALSE)
    }
    .validate_core_numeric_input(Xz, caller = caller, observation_axis = "rows")
    if (nrow(W) != nrow(Xz) || ncol(W) != nrow(Xz)) {
        stop(
            caller, ": W must be square with nrow(W) == nrow(Xz); got ",
            nrow(W), "x", ncol(W), " versus ", nrow(Xz), " rows in Xz.",
            call. = FALSE
        )
    }
    if (any(!is.finite(Xz))) {
        stop(caller, ": Xz contains non-finite values.", call. = FALSE)
    }
    if (length(W@x) && any(!is.finite(W@x))) {
        stop(caller, ": W contains non-finite values.", call. = FALSE)
    }

    list(Xz = Xz, W = W)
}

#' Unified compute block for Lee's L chunked computation
#' @description
#' Internal helper for chunked Lee's L computation. Abstracts the compute_block
#' logic that was previously defined separately in within=TRUE and within=FALSE branches.
#' @param Xz Numeric matrix (genes in columns).
#' @param W Sparse weight matrix.
#' @param chunk_pos Integer vector of column indices (1-based).
#' @param backend Character. Backend to use ("r" or "cpp").
#' @param ncores Integer. Number of threads for native backend.
#' @return Numeric matrix of Lee's L values for the specified columns.
#' @keywords internal
.compute_block_unified <- function(Xz, W, chunk_pos, backend = "cpp", ncores = 1L) {
    if (identical(backend, "r")) {
        .lee_l_cols_r(Xz, W, chunk_pos)
    } else {
        tryCatch(
            .lee_l_cols(Xz, W, cols0 = as.integer(chunk_pos - 1L), n_threads = ncores),
            error = function(e) {
                .lee_l_auto_fallback_warning(e, chunked = TRUE)
                .lee_l_cols_r(Xz, W, chunk_pos)
            }
        )
    }
}

.lee_l_cache_r <- function(Xz, W) {
    validated <- .validate_lee_l_native_inputs(Xz, W, caller = ".lee_l_cache_r")
    Xz <- validated$Xz
    W <- validated$W

    S0 <- sum(W)
    if (!is.finite(S0)) {
        stop(".lee_l_cache_r: W sum is non-finite.", call. = FALSE)
    }
    if (S0 == 0) {
        out <- matrix(0, ncol(Xz), ncol(Xz))
        dimnames(out) <- NULL
        return(out)
    }

    dz2 <- colSums(Xz^2)
    WZ <- as.matrix(W %*% Xz)
    num <- crossprod(Xz, WZ)
    den <- sqrt(outer(dz2, dz2))
    out <- (nrow(Xz) / S0) * (num / den)
    out[!is.finite(out)] <- 0
    dimnames(out) <- NULL
    out
}

.lee_l_cols_r <- function(Xz, W, col_indices) {
    validated <- .validate_lee_l_native_inputs(Xz, W, caller = ".lee_l_cols_r")
    Xz <- validated$Xz
    W <- validated$W

    col_indices <- as.integer(col_indices)
    if (!length(col_indices)) {
        stop(".lee_l_cols_r: col_indices must contain at least one value.", call. = FALSE)
    }
    if (anyNA(col_indices) || any(col_indices < 1L) || any(col_indices > ncol(Xz))) {
        stop(
            ".lee_l_cols_r: col_indices must be within [1, ", ncol(Xz), "].",
            call. = FALSE
        )
    }

    S0 <- sum(W)
    if (!is.finite(S0)) {
        stop(".lee_l_cols_r: W sum is non-finite.", call. = FALSE)
    }
    if (S0 == 0) {
        out <- matrix(0, ncol(Xz), length(col_indices))
        dimnames(out) <- NULL
        return(out)
    }

    dz2_all <- colSums(Xz^2)
    Xz_sub <- Xz[, col_indices, drop = FALSE]
    WZ_sub <- as.matrix(W %*% Xz_sub)
    num <- crossprod(Xz, WZ_sub)
    den <- sqrt(outer(dz2_all, dz2_all[col_indices]))
    out <- (nrow(Xz) / S0) * (num / den)
    out[!is.finite(out)] <- 0
    dimnames(out) <- NULL
    out
}

#' Compute Lee L
#' @description
#' Internal helper for `.compute_lee_l`.
#' @param scope_obj A `scope_object`.
#' @param grid_name Grid layer name (or `NULL` to use the active layer).
#' @param norm_layer Layer name.
#' @param genes Parameter value.
#' @param within Parameter value.
#' @param ncores Number of cores/threads to use.
#' @param mem_limit_GB Parameter value.
#' @param chunk_size Parameter value.
#' @param use_bigmemory Logical flag.
#' @param backing_path Filesystem path.
#' @param block_side Parameter value.
#' @param cache_inputs Parameter value.
#' @param input_cache Parameter value.
#' @return Return value used internally.
#' @keywords internal
.compute_lee_l <- function(scope_obj,
                         grid_name = NULL,
                         norm_layer = "Xz",
                         genes = NULL,
                         within = TRUE,
                         ncores = 1,
                         mem_limit_GB = 16,
                         chunk_size = 256L,
                         use_bigmemory = FALSE,
                         backing_path = tempdir(),
                         block_side = 8,
                         cache_inputs = FALSE,
                         input_cache = NULL) {
    if (use_bigmemory && !requireNamespace("bigmemory", quietly = TRUE)) {
        use_bigmemory <- FALSE
    }

    ## ---- 0. Get grid layer (new helper) ---------------------------------
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    # Fill back layer name string, for later writing back to stats
    if (is.null(grid_name)) {
        grid_name <- names(scope_obj@grid)[
            vapply(scope_obj@grid, identical, logical(1), g_layer)
        ]
    }

    ## ---- 1. Extract expression matrix and weight matrix (with optional reuse) ----
    grid_info <- g_layer$grid_info
    cache_valid <- cache_inputs &&
        !is.null(input_cache) &&
        identical(input_cache$norm_layer, norm_layer) &&
        identical(input_cache$grid_id, grid_info$grid_id) &&
        identical(input_cache$block_side, block_side)

    if (cache_valid) {
        Xz_full <- input_cache$Xz_full
        W <- input_cache$W
        block_id <- input_cache$block_id
    } else {
        Xz_full <- as.matrix(g_layer[[norm_layer]])
        ord <- match(grid_info$grid_id, rownames(Xz_full))
        Xz_full <- Xz_full[ord, , drop = FALSE]
        w_style_attr <- attr(g_layer$W, "weight_style", exact = TRUE)
        W <- g_layer$W[ord, ord]
        if (!is.null(w_style_attr) && is.null(attr(W, "weight_style", exact = TRUE))) {
            attr(W, "weight_style") <- w_style_attr
        }
        block_id <- .assign_block_id(grid_info, block_side = block_side)
    }

    block_id <- .compact_block_id(block_id)

    all_genes <- colnames(Xz_full)
    idx_keep <- .resolve_leeL_gene_indices(genes, all_genes)
    validated <- .validate_lee_l_native_inputs(Xz_full, W, caller = ".compute_lee_l")
    Xz_full <- validated$Xz
    W <- validated$W
    backend <- .resolve_lee_l_backend()
    if (identical(backend, "python")) {
        stop(.lee_l_python_backend_unavailable(), call. = FALSE)
    }

    n_g <- length(idx_keep)
    n_cells <- nrow(Xz_full)
    n_g_full <- ncol(Xz_full)
    output_ncol <- if (within) n_g else n_g_full
    bytes_output <- (n_g * output_ncol * 8) / 1024^3
    bytes_native_full <- ((n_g_full * n_g_full * 8) + (n_cells * n_g_full * 8)) / 1024^3

    if (!use_bigmemory && bytes_native_full > mem_limit_GB) {
        stop(
            ".compute_lee_l: estimated Lee's L working set (~",
            sprintf("%.2f", bytes_native_full),
            " GB) exceeds mem_limit_GB=", sprintf("%.2f", mem_limit_GB),
            " with use_bigmemory=FALSE. Enable chunked mode or reduce the gene set.",
            call. = FALSE
        )
    }

    need_stream <- use_bigmemory && bytes_native_full > mem_limit_GB

    ## ======================================================
    ## =============  A. One‑shot computation (fits in RAM) ===============
    ## ======================================================
    if (!need_stream) {
        RhpcBLASctl::blas_set_num_threads(1)
        native_call <- function() .lee_l_cache(Xz_full, W, n_threads = ncores)
        if (identical(backend, "r")) {
            L_full <- .lee_l_cache_r(Xz_full, W)
        } else {
            L_full <- tryCatch(
                native_call(),
                error = function(e) {
                    .lee_l_auto_fallback_warning(e)
                    .lee_l_cache_r(Xz_full, W)
                }
            )
        }

        if (within) {
            Lmat <- L_full[idx_keep, idx_keep, drop = FALSE]
            Xuse <- Xz_full[, idx_keep, drop = FALSE]
        } else {
            Lmat <- L_full[idx_keep, , drop = FALSE]
            Xuse <- Xz_full[, idx_keep, drop = FALSE]
        }

        ## ======================================================
        ## =============  B. Chunked computation with file mapping ================
        ## ======================================================
    } else {
        ## Use in-memory shared big.matrix to avoid unsupported 'shared' arg on filebacked in some versions
        ## Note: this stores the chunked result in RAM. Ensure sufficient memory is available.
        L_bm <- bigmemory::big.matrix(
            nrow = n_g, ncol = output_ncol, type = "double",
            init = NA_real_,
            dimnames = NULL,
            shared = TRUE
        )

        RhpcBLASctl::blas_set_num_threads(1)

        ## ---- Process column blocks and write ----
        if (within) {
            Xz_work <- Xz_full[, idx_keep, drop = FALSE]

            for (start in seq(1L, n_g, by = chunk_size)) {
                chunk_pos <- start:min(n_g, start + chunk_size - 1L)
                L_block <- .compute_block_unified(Xz_work, W, chunk_pos, backend = backend, ncores = ncores)

                ## Write column block + symmetric columns
                L_bm[, chunk_pos] <- L_block
                L_bm[chunk_pos, ] <- t(L_block)

                rm(L_block)
                gc(verbose = FALSE)
            }
        } else {
            for (start in seq(1L, n_g, by = chunk_size)) {
                chunk_pos <- start:min(n_g, start + chunk_size - 1L)
                full_idx_chunk <- idx_keep[chunk_pos]
                L_block <- .compute_block_unified(Xz_full, W, full_idx_chunk, backend = backend, ncores = ncores)
                L_bm[chunk_pos, ] <- t(L_block)

                rm(L_block)
                gc(verbose = FALSE)
            }
        }

        ## ---- Add gene names (only two vectors) ----
        dimnames(L_bm) <- list(
            all_genes[idx_keep],
            if (within) all_genes[idx_keep] else all_genes
        )

        Lmat <- L_bm
        Xuse <- Xz_full[, idx_keep, drop = FALSE]
    }

	    dimnames(Lmat) <- if (within) {
	        list(row = colnames(Xuse), col = colnames(Xuse))
	    } else {
	        list(row = colnames(Xuse), col = all_genes)
	    }

    list(
        Lmat      = Lmat,
        X_used    = Xuse,
        X_full    = Xz_full,
        cells     = rownames(Xz_full),
        W         = W,
        grid_info = g_layer$grid_info,
        grid_name = grid_name, # ← Additional return, convenient for .compute_l
        block_id  = block_id,
        input_cache = if (cache_inputs) {
            list(
                norm_layer = norm_layer,
                grid_id = grid_info$grid_id,
                block_side = block_side,
                Xz_full = Xz_full,
                W = W,
                block_id = block_id
            )
        } else {
            NULL
        }
    )
}

#' Lee L Perm Block
#' @description
#' Internal helper for `.lee_l_perm_block`.
#' @param Xz Parameter value.
#' @param W Parameter value.
#' @param L_ref Parameter value.
#' @param block_id Parameter value.
#' @param perms Parameter value.
#' @param block_size Parameter value.
#' @param n_threads Number of threads to use.
#' @return Return value used internally.
#' @keywords internal
.lee_l_perm_block <- function(Xz, W, L_ref,
                             block_id,
                             perms = 999,
                             block_size = 64,
                             n_threads = 1) {
    stopifnot(
        is.matrix(Xz), inherits(W, "dgCMatrix"),
        length(block_id) == nrow(Xz)
    )

    RhpcBLASctl::blas_set_num_threads(1)
    ngen <- ncol(Xz)
    geCnt <- matrix(0, ngen, ngen)
    done <- 0L

    ## ---- Pre-group row indices by block ----
    split_rows <- split(seq_along(block_id), block_id)

    while (done < perms) {
        bsz <- min(block_size, perms - done)

        idx_mat <- replicate(bsz,
            {
                # Preserve block positions and shuffle rows within each block.
                # This implements the documented constrained block-wise null.
                unlist(lapply(split_rows, function(x) x[sample.int(length(x))]), use.names = FALSE)
            },
            simplify = "matrix"
        ) - 1L # 0-based for C++

        storage.mode(idx_mat) <- "integer"
        # Use correct export function name
        geCnt <- geCnt + .lee_perm_block(Xz, W, idx_mat, as.integer(block_id) - 1L, L_ref, n_threads)
        done <- done + bsz
    }
    (geCnt + 1) / (perms + 1)
}

#' Lee L Full
#' @description
#' Internal helper for `.lee_l_full`.
#' @param Xz Parameter value.
#' @param W Parameter value.
#' @param n_threads Number of threads to use.
#' @param .cache Parameter value.
#' @return Return value used internally.
#' @keywords internal
.lee_l_full <- function(Xz, W, n_threads = 1L, .cache = TRUE) {
    if (.cache) .lee_l_cache(Xz, W, n_threads) else .lee_l(Xz, W, n_threads)
}

#' Lee L Subset
#' @description
#' Internal helper for `.lee_l_subset`.
#' @param Xz Parameter value.
#' @param W Parameter value.
#' @param cols0 Parameter value.
#' @param n_threads Number of threads to use.
#' @return Return value used internally.
#' @keywords internal
.lee_l_subset <- function(Xz, W, cols0, n_threads = 1L) {
    .lee_l_cols(Xz, W, cols0, n_threads)
}
