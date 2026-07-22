#' Identifier for the canonical Lee's L implementation.
#' @keywords internal
.lee_formula_id <- function() "Lee_S2_WtW_v1"

#' Recognize formula identifiers emitted by corrected geneSCOPE releases.
#' @keywords internal
.lee_formula_id_is_supported <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) &&
        x %in% c(.lee_formula_id(), "Lee2009_S2_v1")
}

#' Validate a positive integer Lee workflow parameter.
#' @keywords internal
.lee_assert_positive_integer <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        x < 1 || x > .Machine$integer.max || x != floor(x)) {
        stop(name, " must be a single positive integer.", call. = FALSE)
    }
    as.integer(x)
}

#' Validate a non-negative integer Lee workflow parameter.
#' @keywords internal
.lee_assert_nonnegative_integer <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
        x < 0 || x > .Machine$integer.max || x != floor(x)) {
        stop(name, " must be a single non-negative integer.", call. = FALSE)
    }
    as.integer(x)
}

#' Validate a positive finite Lee workflow parameter.
#' @keywords internal
.lee_assert_positive_finite <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
        stop(
            name, " must be a single positive finite number.",
            if (identical(name, "mem_limit_GB")) {
                " An invalid limit is fail-closed because the requested output exceeds mem_limit_GB."
            } else {
                ""
            },
            call. = FALSE
        )
    }
    as.numeric(x)
}

#' Validate a scalar logical Lee workflow parameter.
#' @keywords internal
.lee_assert_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop(name, " must be a single non-missing logical value.", call. = FALSE)
    }
    x
}

#' Validate a scalar layer name.
#' @keywords internal
.lee_assert_layer_name <- function(x, name) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
        stop(name, " must be a single non-empty character string.", call. = FALSE)
    }
    x
}

#' Resolve the weight style recorded on a Lee weight matrix.
#' @keywords internal
.lee_weight_style <- function(W) {
    style <- attr(W, "weight_style", exact = TRUE)
    if (is.null(style)) {
        provenance <- attr(W, "weight_provenance", exact = TRUE)
        if (is.list(provenance)) style <- provenance$weight_style
    }
    style <- if (is.null(style) || length(style) != 1L || is.na(style)) {
        ""
    } else {
        toupper(as.character(style))
    }
    if (style %in% c("B", "W")) return(style)
    nonzero <- if (inherits(W, "sparseMatrix")) W@x else as.matrix(W)[as.matrix(W) != 0]
    if (length(nonzero) && all(nonzero == 1)) "B" else "W"
}

#' Validate a normalized expression layer used by Lee's L.
#' @description
#' Lee's formula assumes centred variables. Unit variance is not required
#' because the denominator removes positive per-column scale factors.
#' @keywords internal
.validate_lee_norm_layer <- function(X, norm_layer, context = "Lee's L") {
    norm_layer <- .lee_assert_layer_name(norm_layer, "norm_layer")
    X <- tryCatch(as.matrix(X), error = function(e) {
        stop(
            context, " expression layer '", norm_layer,
            "' cannot be converted to a matrix: ", conditionMessage(e),
            call. = FALSE
        )
    })
    if (!is.numeric(X) || length(dim(X)) != 2L || any(dim(X) < 1L)) {
        stop(
            context, " expression layer '", norm_layer,
            "' must be a non-empty numeric matrix.", call. = FALSE
        )
    }
    if (any(!is.finite(X))) {
        stop(
            context, " expression layer '", norm_layer,
            "' contains non-finite values.", call. = FALSE
        )
    }
    centre <- colMeans(X)
    rms <- sqrt(colMeans(X * X))
    relative_offset <- abs(centre) / pmax(rms, .Machine$double.eps)
    bad <- which(!is.finite(relative_offset) | relative_offset > 1e-7)
    if (length(bad)) {
        labels <- colnames(X)
        labels <- if (is.null(labels)) as.character(bad) else labels[bad]
        stop(
            context, " expression layer '", norm_layer,
            "' must be column-centred (abs(mean)/RMS <= 1e-7); offending columns: ",
            paste(utils::head(labels, 5L), collapse = ", "),
            if (length(bad) > 5L) " ..." else "",
            ". Unit variance is not required.",
            call. = FALSE
        )
    }
    X
}

#' Compute Lee's canonical spatial normalizer from a weight matrix.
#' @keywords internal
.lee_s2_value <- function(W) {
    if (is.null(dim(W)) || length(dim(W)) != 2L || nrow(W) != ncol(W) || nrow(W) < 1L) {
        stop("W must be a non-empty square matrix before computing Lee S2.", call. = FALSE)
    }
    row_totals <- as.numeric(W %*% rep.int(1, ncol(W)))
    S2 <- sum(row_totals * row_totals)
    if (length(S2) != 1L || !is.finite(S2) || S2 <= 0) {
        stop("W has non-finite or non-positive Lee S2; Lee's L is undefined.", call. = FALSE)
    }
    as.numeric(S2)
}

#' Fingerprint the inputs that determine observed Lee's L and its permutations.
#' @keywords internal
.lee_hash_object <- function(x) {
    # Use one dependency-free algorithm in every environment. Mixing an
    # optional xxhash path with this MD5 path would make identical scope
    # objects fail provenance checks after moving between installations.
    path <- tempfile("genescope-lee-fingerprint-", fileext = ".rdsbin")
    con <- file(path, open = "wb")
    con_open <- TRUE
    on.exit({
        if (con_open) try(close(con), silent = TRUE)
        unlink(path)
    }, add = TRUE)
    serialize(x, con, version = 3L)
    close(con)
    con_open <- FALSE
    unname(tools::md5sum(path)[[1L]])
}

#' Fingerprint the inputs that determine observed Lee's L and its permutations.
#' @keywords internal
.lee_input_fingerprint <- function(Xz, W, grid_info, norm_layer, block_side) {
    norm_layer <- .lee_assert_layer_name(norm_layer, "norm_layer")
    block_side <- .lee_assert_positive_integer(block_side, "block_side")
    hash <- .lee_hash_object
    w_payload <- if (inherits(W, "sparseMatrix")) {
        Wc <- if (inherits(W, "dgCMatrix")) {
            W
        } else {
            methods::as(methods::as(W, "generalMatrix"), "CsparseMatrix")
        }
        list(
            class = "dgCMatrix", dim = dim(Wc), dimnames = dimnames(Wc),
            p = Wc@p, i = Wc@i, x = Wc@x
        )
    } else {
        list(class = class(W), dim = dim(W), dimnames = dimnames(W), values = as.matrix(W))
    }
    grid_frame <- as.data.frame(grid_info, stringsAsFactors = FALSE)
    grid_cols <- intersect(
        c("grid_id", "gx", "gy", "xmin", "xmax", "ymin", "ymax"),
        names(grid_frame)
    )
    grid_payload <- grid_frame[, grid_cols, drop = FALSE]
    data_components <- list(
        X = hash(list(dim = dim(Xz), dimnames = dimnames(Xz), values = Xz)),
        W = hash(w_payload),
        grid = hash(grid_payload),
        weight_style = hash(.lee_weight_style(W))
    )
    permutation_components <- list(block_side = block_side)
    list(
        schema = "lee_input_fingerprint_v2",
        norm_layer = norm_layer,
        data = c(
            list(schema = "lee_observed_data_v1", hash = hash(data_components)),
            data_components
        ),
        permutation = c(
            list(schema = "lee_permutation_config_v1", hash = hash(permutation_components)),
            permutation_components
        )
    )
}

#' Compare only the observed-data part of two Lee fingerprints.
#' @keywords internal
.lee_same_observed_data <- function(observed, current) {
    is.list(observed) && is.list(current) &&
        identical(observed$schema, "lee_input_fingerprint_v2") &&
        identical(current$schema, "lee_input_fingerprint_v2") &&
        identical(observed$norm_layer, current$norm_layer) &&
        identical(observed$data, current$data)
}

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
                if (isTRUE(getOption("geneSCOPE.lee_l.strict_cpp", FALSE))) {
                    stop("Strict C++ Lee's L backend failed: ", conditionMessage(e), call. = FALSE)
                }
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

    S2 <- .lee_s2_value(W)

    dz2 <- colSums(Xz^2)
    if (any(!is.finite(dz2))) {
        stop(".lee_l_cache_r: gene norms are non-finite.", call. = FALSE)
    }
    zero_norm <- dz2 <= 0
    WZ <- as.matrix(W %*% Xz)
    if (any(!is.finite(WZ))) {
        stop(".lee_l_cache_r: W %*% Xz produced non-finite values.", call. = FALSE)
    }
    num <- crossprod(WZ)
    den <- sqrt(outer(dz2, dz2))
    out <- (nrow(Xz) / S2) * (num / den)
    zero_pairs <- outer(zero_norm, zero_norm, `|`)
    if (any(!is.finite(out[!zero_pairs]))) {
        stop(".lee_l_cache_r: Lee's L produced unexpected non-finite values.", call. = FALSE)
    }
    out[zero_pairs] <- 0
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

    S2 <- .lee_s2_value(W)

    dz2_all <- colSums(Xz^2)
    if (any(!is.finite(dz2_all))) {
        stop(".lee_l_cols_r: gene norms are non-finite.", call. = FALSE)
    }
    zero_norm_all <- dz2_all <= 0
    Xz_sub <- Xz[, col_indices, drop = FALSE]
    WZ_all <- as.matrix(W %*% Xz)
    WZ_sub <- as.matrix(W %*% Xz_sub)
    if (any(!is.finite(WZ_all)) || any(!is.finite(WZ_sub))) {
        stop(".lee_l_cols_r: W %*% Xz produced non-finite values.", call. = FALSE)
    }
    num <- crossprod(WZ_all, WZ_sub)
    den <- sqrt(outer(dz2_all, dz2_all[col_indices]))
    out <- (nrow(Xz) / S2) * (num / den)
    zero_pairs <- outer(zero_norm_all, zero_norm_all[col_indices], `|`)
    if (any(!is.finite(out[!zero_pairs]))) {
        stop(".lee_l_cols_r: Lee's L produced unexpected non-finite values.", call. = FALSE)
    }
    out[zero_pairs] <- 0
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
    within <- .lee_assert_flag(within, "within")
    ncores <- .lee_assert_positive_integer(ncores, "ncores")
    block_side <- .lee_assert_positive_integer(block_side, "block_side")
    cache_inputs <- .lee_assert_flag(cache_inputs, "cache_inputs")
    use_bigmemory <- .lee_assert_flag(use_bigmemory, "use_bigmemory")
    norm_layer <- .lee_assert_layer_name(norm_layer, "norm_layer")
    if (isTRUE(use_bigmemory)) {
        stop(
            "use_bigmemory=TRUE is disabled: the native Lee L chunk route allocates a shared RAM-backed big.matrix, not a safe file-backed result.",
            call. = FALSE
        )
    }
    mem_limit_GB <- .lee_assert_positive_finite(mem_limit_GB, "mem_limit_GB")

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
    if (is.null(g_layer[[norm_layer]])) {
        stop("Grid layer is missing norm_layer='", norm_layer, "'.", call. = FALSE)
    }
    Xz_source <- .validate_lee_norm_layer(
        g_layer[[norm_layer]], norm_layer,
        context = "[geneSCOPE::.compute_lee_l]"
    )
    W_source <- g_layer$W
    if (is.null(W_source)) stop("Grid layer is missing W; run computeWeights() first.", call. = FALSE)
    source_weight_style <- attr(W_source, "weight_style", exact = TRUE) %||%
        g_layer$weights_meta$weight_style %||% "unknown"
    if (is.null(attr(W_source, "weight_style", exact = TRUE))) {
        attr(W_source, "weight_style") <- as.character(source_weight_style)[1L]
    }
    source_fingerprint <- .lee_input_fingerprint(
        Xz = Xz_source,
        W = W_source,
        grid_info = grid_info,
        norm_layer = norm_layer,
        block_side = block_side
    )
    alignment_schema <- "grid_id_independent_fingerprint_v4"
    cache_valid <- cache_inputs &&
        !is.null(input_cache) &&
        identical(input_cache$alignment_schema, alignment_schema) &&
        identical(input_cache$norm_layer, norm_layer) &&
        identical(input_cache$grid_id, grid_info$grid_id) &&
        identical(input_cache$block_side, block_side) &&
        identical(input_cache$source_fingerprint, source_fingerprint)

    target_ids <- as.character(grid_info$grid_id)
    if (!length(target_ids) || anyNA(target_ids) || any(!nzchar(target_ids)) || anyDuplicated(target_ids)) {
        stop("grid_info$grid_id must be present, non-empty, and unique before Lee's L alignment.", call. = FALSE)
    }
    if (cache_valid) {
        Xz_full <- input_cache$Xz_full
        W <- input_cache$W
        block_id <- input_cache$block_id
        cached_fingerprint <- tryCatch(
            .lee_input_fingerprint(
                Xz_full, W, grid_info,
                norm_layer = norm_layer,
                block_side = block_side
            ),
            error = function(e) NULL
        )
        expected_block_id <- tryCatch(
            .compact_block_id(.assign_block_id(grid_info, block_side = block_side)),
            error = function(e) NULL
        )
        cached_block_id <- tryCatch(
            as.integer(.compact_block_id(block_id)),
            error = function(e) NULL
        )
        cache_valid <-
            is.matrix(Xz_full) &&
            identical(rownames(Xz_full), target_ids) &&
            identical(rownames(W), target_ids) &&
            identical(colnames(W), target_ids) &&
            length(block_id) == length(target_ids) &&
            identical(cached_block_id, as.integer(expected_block_id)) &&
            identical(cached_fingerprint, input_cache$input_fingerprint)
    }
    if (!cache_valid) {
        Xz_full <- as.matrix(Xz_source)

        .alignment_index <- function(source_ids, label) {
            if (is.null(source_ids)) {
                stop(label, " has no grid_id dimnames; rerun the producing step.", call. = FALSE)
            }
            source_ids <- as.character(source_ids)
            if (anyNA(source_ids) || anyDuplicated(source_ids)) {
                stop(label, " grid_id dimnames contain missing or duplicated values.", call. = FALSE)
            }
            idx <- match(target_ids, source_ids)
            if (anyNA(idx) || length(source_ids) != length(target_ids)) {
                stop(label, " grid_id set does not exactly match grid_info$grid_id.", call. = FALSE)
            }
            idx
        }

        x_ord <- .alignment_index(rownames(Xz_full), norm_layer)
        w_row_ord <- .alignment_index(rownames(W_source), "W rows")
        w_col_ord <- .alignment_index(colnames(W_source), "W columns")

        Xz_full <- Xz_full[x_ord, , drop = FALSE]
        w_style_attr <- attr(W_source, "weight_style", exact = TRUE)
        W <- W_source[w_row_ord, w_col_ord, drop = FALSE]
        if (!is.null(w_style_attr) && is.null(attr(W, "weight_style", exact = TRUE))) {
            attr(W, "weight_style") <- w_style_attr
        }
        rownames(Xz_full) <- target_ids
        dimnames(W) <- list(target_ids, target_ids)

        if (!identical(rownames(Xz_full), target_ids) ||
            !identical(rownames(W), target_ids) ||
            !identical(colnames(W), target_ids)) {
            stop("Failed to align Xz and W independently to grid_info$grid_id.", call. = FALSE)
        }
        block_id <- .assign_block_id(grid_info, block_side = block_side)
    }

    block_id <- .compact_block_id(block_id)

    input_fingerprint <- .lee_input_fingerprint(
        Xz = Xz_full,
        W = W,
        grid_info = grid_info,
        norm_layer = norm_layer,
        block_side = block_side
    )
    weight_style <- .lee_weight_style(W)
    S2 <- .lee_s2_value(W)

    all_genes <- colnames(Xz_full)
    if (is.null(all_genes) || anyNA(all_genes) || any(!nzchar(all_genes)) || anyDuplicated(all_genes)) {
        stop("The normalized Lee input must have unique, non-empty gene column names.", call. = FALSE)
    }
    idx_keep <- .resolve_leeL_gene_indices(genes, all_genes)
    validated <- .validate_lee_l_native_inputs(Xz_full, W, caller = ".compute_lee_l")
    Xz_full <- validated$Xz
    W <- validated$W
    backend <- .resolve_lee_l_backend()
    if (identical(backend, "python")) {
        stop(.lee_l_python_backend_unavailable(), call. = FALSE)
    }

    Xz_compute <- if (isTRUE(within)) Xz_full[, idx_keep, drop = FALSE] else Xz_full
    n_cells <- nrow(Xz_compute)
    n_g_compute <- ncol(Xz_compute)
    bytes_native_full <- ((n_g_compute * n_g_compute * 8) + (n_cells * n_g_compute * 8)) / 1024^3

    if (bytes_native_full > mem_limit_GB) {
        stop(
            ".compute_lee_l: estimated Lee's L working set (~",
            sprintf("%.2f", bytes_native_full),
            " GB) exceeds mem_limit_GB=", sprintf("%.2f", mem_limit_GB),
            " with use_bigmemory=FALSE. Reduce the gene set or raise mem_limit_GB only after confirming available RAM.",
            call. = FALSE
        )
    }

    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        RhpcBLASctl::blas_set_num_threads(1)
    }
    native_call <- function() .lee_l_cache(Xz_compute, W, n_threads = ncores)
    if (identical(backend, "r")) {
        L_full <- .lee_l_cache_r(Xz_compute, W)
    } else {
        L_full <- tryCatch(
            native_call(),
            error = function(e) {
                if (isTRUE(getOption("geneSCOPE.lee_l.strict_cpp", FALSE))) {
                    stop("Strict C++ Lee's L backend failed: ", conditionMessage(e), call. = FALSE)
                }
                .lee_l_auto_fallback_warning(e)
                .lee_l_cache_r(Xz_compute, W)
            }
        )
    }

    if (within) {
        Lmat <- L_full
        Xuse <- Xz_compute
    } else {
        Lmat <- L_full[idx_keep, , drop = FALSE]
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
        input_fingerprint = input_fingerprint,
        weight_style = weight_style,
        S2 = S2,
        input_cache = if (cache_inputs) {
            list(
                alignment_schema = alignment_schema,
                norm_layer = norm_layer,
                grid_id = grid_info$grid_id,
                block_side = block_side,
                source_fingerprint = source_fingerprint,
                input_fingerprint = input_fingerprint,
                weight_style = weight_style,
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
    perms <- .lee_validate_integer_scalar(perms, "perms", lower = 0L)

    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        RhpcBLASctl::blas_set_num_threads(1)
    }
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
                # Sampling each block and concatenating the results is not
                # equivalent when block positions are non-contiguous.
                idx <- seq_len(nrow(Xz))
                for (rows in split_rows) {
                    idx[rows] <- rows[sample.int(length(rows), length(rows), replace = FALSE)]
                }
                idx
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
