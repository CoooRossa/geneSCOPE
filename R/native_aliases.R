# Internal aliases for RcppExports functions with dot-prefixed names.

# Round 4: Darwin native safety - unified check for all native backends
.native_all_disabled <- function() {
    isTRUE(getOption("geneSCOPE.disable_native_all", FALSE))
}

.native_all_explicitly_enabled <- function() {
    identical(getOption("geneSCOPE.disable_native_all", NULL), FALSE)
}

.native_backend_disabled <- function(option) {
    if (.native_all_explicitly_enabled()) {
        return(FALSE)
    }
    isTRUE(getOption(option, FALSE)) || .native_all_disabled()
}

.native_permutation_disabled <- function() {
    .native_backend_disabled("geneSCOPE.disable_native_permutation")
}

.lee_l_native_disabled <- function() {
    .native_backend_disabled("geneSCOPE.disable_native_lee_l_backend")
}

.validate_native_thread_count <- function(n_threads, caller) {
    if (!is.numeric(n_threads) || length(n_threads) != 1L || is.na(n_threads) ||
        !is.finite(n_threads) || n_threads < 1 || n_threads > .Machine$integer.max ||
        n_threads != floor(n_threads)) {
        stop(caller, ": n_threads must be at least 1 and integer-valued.", call. = FALSE)
    }
    as.integer(n_threads)
}

.validate_lee_l_native_thread_count <- function(n_threads, caller) {
    .validate_native_thread_count(n_threads, caller)
}

.native_lee_L_call <- function(Xz, W, n_threads = 1L) {
    .Call(`_geneSCOPE_lee_L`, Xz, W, n_threads)
}

.native_lee_L_cache_call <- function(Xz, W, n_threads = 1L) {
    .Call(`_geneSCOPE_lee_L_cache`, Xz, W, n_threads)
}

.native_lee_L_cols_call <- function(Xz, W, cols0, n_threads = 1L) {
    .Call(`_geneSCOPE_lee_L_cols`, Xz, W, cols0, n_threads)
}

.native_lee_perm_call <- function(Xz, W, idx_mat, L_ref, n_threads = 1L) {
    .Call(`_geneSCOPE_lee_perm`, Xz, W, idx_mat, L_ref, n_threads)
}

.native_lee_perm_block_call <- function(Xz, W, idx_mat, block_ids, L_ref, n_threads = 1L) {
    .Call(`_geneSCOPE_lee_perm_block`, Xz, W, idx_mat, block_ids, L_ref, n_threads)
}

.validate_lee_permutation_indices <- function(idx_mat, n, caller) {
    idx_mat <- as.matrix(idx_mat)
    if (!is.numeric(idx_mat) || !ncol(idx_mat) || nrow(idx_mat) != n ||
        anyNA(idx_mat) || any(!is.finite(idx_mat)) || any(idx_mat != floor(idx_mat))) {
        stop(caller, ": idx_mat must be an n x B integer matrix with B >= 1.", call. = FALSE)
    }
    if (any(idx_mat < 0) || any(idx_mat >= n)) {
        stop(caller, ": idx_mat contains an out-of-range index.", call. = FALSE)
    }
    expected <- seq.int(0L, n - 1L)
    for (b in seq_len(ncol(idx_mat))) {
        if (!identical(sort(as.integer(idx_mat[, b])), expected)) {
            stop(caller, ": every idx_mat column must be a bijection of 0:(n-1).", call. = FALSE)
        }
    }
    storage.mode(idx_mat) <- "integer"
    idx_mat
}

.validate_lee_reference <- function(L_ref, g, caller) {
    L_ref <- as.matrix(L_ref)
    if (!identical(dim(L_ref), c(g, g))) {
        stop(caller, ": L_ref must be a square matrix with dimensions matching ncol(Xz).", call. = FALSE)
    }
    if (any(!is.finite(L_ref))) {
        stop(caller, ": L_ref contains a non-finite value.", call. = FALSE)
    }
    L_ref
}

lee_L <- function(Xz, W, n_threads = 1L) {
    n_threads <- .validate_lee_l_native_thread_count(n_threads, "lee_L")
    # Round 4: Darwin native safety - check before any native call
    if (.lee_l_native_disabled()) {
        return(.lee_l_cache_r(Xz, W))
    }
    
    validated <- .validate_lee_l_native_inputs(Xz, W, caller = "lee_L")
    .native_lee_L_call(validated$Xz, validated$W, n_threads)
}

lee_L_cache <- function(Xz, W, n_threads = 1L) {
    n_threads <- .validate_lee_l_native_thread_count(n_threads, "lee_L_cache")
    # Round 4: Darwin native safety - check before any native call
    if (.lee_l_native_disabled()) {
        return(.lee_l_cache_r(Xz, W))
    }
    
    validated <- .validate_lee_l_native_inputs(Xz, W, caller = "lee_L_cache")
    .native_lee_L_cache_call(validated$Xz, validated$W, n_threads)
}

lee_L_cols <- function(Xz, W, cols0, n_threads = 1L) {
    n_threads <- .validate_lee_l_native_thread_count(n_threads, "lee_L_cols")
    validated <- .validate_lee_l_native_inputs(Xz, W, caller = "lee_L_cols")
    if (!is.numeric(cols0) || !length(cols0) || anyNA(cols0) ||
        any(!is.finite(cols0)) || any(cols0 != floor(cols0))) {
        stop("cols0 must contain at least one integer value.", call. = FALSE)
    }
    cols0 <- as.integer(cols0)
    if (any(cols0 < 0L | cols0 >= ncol(validated$Xz))) {
        stop("cols0 contains an out-of-range index.", call. = FALSE)
    }
    if (.lee_l_native_disabled()) {
        return(.lee_l_cols_r(validated$Xz, validated$W, cols0 + 1L))
    }
    .native_lee_L_cols_call(validated$Xz, validated$W, cols0, n_threads)
}

lee_perm <- function(Xz, W, idx_mat, L_ref, n_threads = 1L) {
    n_threads <- .validate_lee_l_native_thread_count(n_threads, "lee_perm")
    validated <- .validate_lee_l_native_inputs(Xz, W, caller = "lee_perm")
    idx_mat <- .validate_lee_permutation_indices(idx_mat, nrow(validated$Xz), "lee_perm")
    L_ref <- .validate_lee_reference(L_ref, ncol(validated$Xz), "lee_perm")
    if (.native_permutation_disabled()) {
        return(.lee_perm_r(validated$Xz, validated$W, idx_mat, L_ref))
    }
    .native_lee_perm_call(
        validated$Xz,
        validated$W,
        idx_mat,
        L_ref,
        n_threads
    )
}

lee_perm_block <- function(Xz, W, idx_mat, block_ids, L_ref, n_threads = 1L) {
    n_threads <- .validate_lee_l_native_thread_count(n_threads, "lee_perm_block")
    validated <- .validate_lee_l_native_inputs(Xz, W, caller = "lee_perm_block")
    idx_mat <- .validate_lee_permutation_indices(idx_mat, nrow(validated$Xz), "lee_perm_block")
    if (!is.numeric(block_ids) || length(block_ids) != nrow(validated$Xz) ||
        anyNA(block_ids) || any(!is.finite(block_ids)) || any(block_ids != floor(block_ids))) {
        stop("block_ids must be an integer vector with length nrow(Xz).", call. = FALSE)
    }
    block_ids <- as.integer(block_ids)
    for (b in seq_len(ncol(idx_mat))) {
        if (any(block_ids[idx_mat[, b] + 1L] != block_ids)) {
            stop("lee_perm_block: idx_mat moves an observation across blocks.", call. = FALSE)
        }
    }
    L_ref <- .validate_lee_reference(L_ref, ncol(validated$Xz), "lee_perm_block")
    if (.native_permutation_disabled()) {
        return(.lee_perm_r(validated$Xz, validated$W, idx_mat, L_ref))
    }
    .native_lee_perm_block_call(
        validated$Xz,
        validated$W,
        idx_mat,
        block_ids,
        L_ref,
        n_threads
    )
}

.lee_l <- function(...) {
    lee_L(...)
}
.lee_l_cache <- function(...) {
    lee_L_cache(...)
}
.lee_l_cols <- function(...) {
    lee_L_cols(...)
}
.lee_perm <- function(...) {
    lee_perm(...)
}
.lee_perm_block <- function(...) {
    lee_perm_block(...)
}

# Exact R fallback for the common-row Lee permutation count.  This is used by
# the macOS safety policy so observed L and permutation inference cannot select
# incompatible backends.
.lee_perm_r <- function(Xz, W, idx_mat, L_ref) {
    counts <- matrix(0, nrow(L_ref), ncol(L_ref), dimnames = dimnames(L_ref))
    for (b in seq_len(ncol(idx_mat))) {
        L_perm <- .lee_l_cache_r(Xz[idx_mat[, b] + 1L, , drop = FALSE], W)
        counts <- counts + (abs(L_perm) >= abs(L_ref))
    }
    counts
}

.idelta_sparse_cpp <- function(...) idelta_sparse_cpp(...)

.pearson_native_disabled <- function() {
    .native_backend_disabled("geneSCOPE.disable_native_correlation_backend")
}

.native_pearson_block_cpp_call <- function(X, bs = 2000L, n_threads = 1L) {
    .Call(`_geneSCOPE_pearson_block_cpp`, X, bs, n_threads)
}

.native_pearson_cor_call <- function(X, bs = 2000L, n_threads = 1L) {
    .Call(`_geneSCOPE_pearson_cor`, X, bs, n_threads)
}

pearson_block_cpp <- function(X, bs = 2000L, n_threads = 1L) {
    X <- as.matrix(X)
    .validate_core_numeric_input(X, caller = "pearson_block_cpp", observation_axis = "rows")
    bs <- suppressWarnings(as.integer(bs)[1L])
    if (is.na(bs) || bs < 1L) {
        stop("pearson_block_cpp: bs must be >= 1.", call. = FALSE)
    }
    if (.pearson_native_disabled()) {
        return(.compute_correlation_dense_pearson_r(X))
    }
    .native_pearson_block_cpp_call(X, bs, .validate_native_thread_count(n_threads, "pearson_block_cpp"))
}

pearson_cor <- function(X, bs = 2000L, n_threads = 1L) {
    X <- as.matrix(X)
    .validate_core_numeric_input(X, caller = "pearson_cor", observation_axis = "rows")
    bs <- suppressWarnings(as.integer(bs)[1L])
    if (is.na(bs) || bs < 1L) {
        stop("pearson_cor: bs must be >= 1.", call. = FALSE)
    }
    if (.pearson_native_disabled()) {
        return(.compute_correlation_dense_pearson_r(X))
    }
    .native_pearson_cor_call(X, bs, .validate_native_thread_count(n_threads, "pearson_cor"))
}

.pearson_block_cpp <- function(...) {
    if (.pearson_native_disabled()) {
        stop("pearson_block_cpp native backend disabled by option; using R Pearson backend.", call. = FALSE)
    }
    pearson_block_cpp(...)
}
.pearson_cor <- function(...) {
    if (.pearson_native_disabled()) {
        stop("pearson_cor native backend disabled by option; using R Pearson backend.", call. = FALSE)
    }
    pearson_cor(...)
}

.morisita_horn_sparse <- function(...) morisita_horn_sparse(...)
.loess_residual_bootstrap <- function(...) loess_residual_bootstrap(...)

.delta_lr_perm_tiny <- function(...) delta_lr_perm_tiny(...)
.delta_lr_perm_block_tiny <- function(...) delta_lr_perm_block_tiny(...)
.delta_lr_perm <- function(...) delta_lr_perm(...)
.delta_lr_perm_block <- function(...) delta_lr_perm_block(...)
.delta_lr_perm_csr <- function(...) delta_lr_perm_csr(...)
.delta_lr_perm_csr_block <- function(...) delta_lr_perm_csr_block(...)
.delta_l_fixed_r_perm <- function(...) delta_l_fixed_r_perm(...)
.delta_l_fixed_r_perm_block <- function(...) delta_l_fixed_r_perm_block(...)
.delta_l_fixed_r_perm_csr <- function(...) delta_l_fixed_r_perm_csr(...)
.delta_l_fixed_r_perm_csr_block <- function(...) delta_l_fixed_r_perm_csr_block(...)
