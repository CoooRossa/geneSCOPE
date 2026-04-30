# Internal aliases for RcppExports functions with dot-prefixed names.

# Round 4: Darwin native safety - unified check for all native backends
.native_all_disabled <- function() {
    isTRUE(getOption("geneSCOPE.disable_native_all", FALSE))
}

.native_permutation_disabled <- function() {
    isTRUE(getOption("geneSCOPE.disable_native_permutation", FALSE)) || .native_all_disabled()
}

.lee_l_native_disabled <- function() {
    isTRUE(getOption("geneSCOPE.disable_native_lee_l_backend", FALSE)) || .native_all_disabled()
}

.validate_native_thread_count <- function(n_threads, caller) {
    n_threads <- suppressWarnings(as.integer(n_threads)[1L])
    if (is.na(n_threads) || n_threads < 1L) {
        stop(caller, ": n_threads must be >= 1.", call. = FALSE)
    }
    n_threads
}

.validate_lee_l_native_thread_count <- function(n_threads, caller) {
    n_threads <- .validate_native_thread_count(n_threads, caller)
    sysname <- tolower(Sys.info()[["sysname"]] %||% "")
    if (
        identical(sysname, "darwin") &&
        n_threads > 1L &&
        !isTRUE(getOption("geneSCOPE.allow_darwin_native_lee_l_threads", FALSE))
    ) {
        return(1L)
    }
    n_threads
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

lee_L <- function(Xz, W, n_threads = 1L) {
    # Round 4: Darwin native safety - check before any native call
    if (.lee_l_native_disabled()) {
        return(.lee_l_cache_r(Xz, W))
    }
    
    validated <- .validate_lee_l_native_inputs(Xz, W, caller = "lee_L")
    .native_lee_L_call(validated$Xz, validated$W, .validate_lee_l_native_thread_count(n_threads, "lee_L"))
}

lee_L_cache <- function(Xz, W, n_threads = 1L) {
    # Round 4: Darwin native safety - check before any native call
    if (.lee_l_native_disabled()) {
        return(.lee_l_cache_r(Xz, W))
    }
    
    validated <- .validate_lee_l_native_inputs(Xz, W, caller = "lee_L_cache")
    .native_lee_L_cache_call(validated$Xz, validated$W, .validate_lee_l_native_thread_count(n_threads, "lee_L_cache"))
}

lee_L_cols <- function(Xz, W, cols0, n_threads = 1L) {
    validated <- .validate_lee_l_native_inputs(Xz, W, caller = "lee_L_cols")
    cols0 <- suppressWarnings(as.integer(cols0))
    if (!length(cols0) || anyNA(cols0)) {
        stop("cols0 must be a non-empty integer vector.", call. = FALSE)
    }
    if (any(cols0 < 0L | cols0 >= ncol(validated$Xz))) {
        stop("cols0 must be within [0, ncol(Xz)-1].", call. = FALSE)
    }
    if (.lee_l_native_disabled()) {
        return(.lee_l_cols_r(validated$Xz, validated$W, cols0 + 1L))
    }
    .native_lee_L_cols_call(validated$Xz, validated$W, cols0, .validate_lee_l_native_thread_count(n_threads, "lee_L_cols"))
}

lee_perm <- function(Xz, W, idx_mat, L_ref, n_threads = 1L) {
    # Round 4: Darwin native safety - check before any native call
    if (.native_permutation_disabled()) {
        stop("lee_perm native backend disabled by option (Darwin safety); permutation path unavailable. ",
             "Set options(geneSCOPE.disable_native_all=FALSE, geneSCOPE.disable_native_permutation=FALSE) to enable.",
             call. = FALSE)
    }
    
    validated <- .validate_lee_l_native_inputs(Xz, W, caller = "lee_perm")
    idx_mat <- as.matrix(idx_mat)
    if (!is.integer(idx_mat)) storage.mode(idx_mat) <- "integer"
    if (nrow(idx_mat) != nrow(validated$Xz)) {
        stop("idx_mat must have nrow(Xz) rows.", call. = FALSE)
    }
    L_ref <- as.matrix(L_ref)
    if (!identical(dim(L_ref), c(ncol(validated$Xz), ncol(validated$Xz)))) {
        stop("L_ref must be a square matrix with dimensions matching ncol(Xz).", call. = FALSE)
    }
    .native_lee_perm_call(
        validated$Xz,
        validated$W,
        idx_mat,
        L_ref,
        .validate_lee_l_native_thread_count(n_threads, "lee_perm")
    )
}

lee_perm_block <- function(Xz, W, idx_mat, block_ids, L_ref, n_threads = 1L) {
    # Round 4: Darwin native safety - check before any native call
    if (.native_permutation_disabled()) {
        stop("lee_perm_block native backend disabled by option (Darwin safety); permutation path unavailable. ",
             "Set options(geneSCOPE.disable_native_all=FALSE, geneSCOPE.disable_native_permutation=FALSE) to enable.",
             call. = FALSE)
    }
    
    validated <- .validate_lee_l_native_inputs(Xz, W, caller = "lee_perm_block")
    idx_mat <- as.matrix(idx_mat)
    if (!is.integer(idx_mat)) storage.mode(idx_mat) <- "integer"
    if (nrow(idx_mat) != nrow(validated$Xz)) {
        stop("idx_mat must have nrow(Xz) rows.", call. = FALSE)
    }
    block_ids <- suppressWarnings(as.integer(block_ids))
    if (length(block_ids) != nrow(validated$Xz) || anyNA(block_ids)) {
        stop("block_ids must be an integer vector with length nrow(Xz).", call. = FALSE)
    }
    L_ref <- as.matrix(L_ref)
    if (!identical(dim(L_ref), c(ncol(validated$Xz), ncol(validated$Xz)))) {
        stop("L_ref must be a square matrix with dimensions matching ncol(Xz).", call. = FALSE)
    }
    .native_lee_perm_block_call(
        validated$Xz,
        validated$W,
        idx_mat,
        block_ids,
        L_ref,
        .validate_lee_l_native_thread_count(n_threads, "lee_perm_block")
    )
}

.lee_l <- function(...) {
    if (.lee_l_native_disabled()) {
        stop("lee_L native backend disabled by option; using R reference backend.", call. = FALSE)
    }
    lee_L(...)
}
.lee_l_cache <- function(...) {
    if (.lee_l_native_disabled()) {
        stop("lee_L_cache native backend disabled by option; using R reference backend.", call. = FALSE)
    }
    lee_L_cache(...)
}
.lee_l_cols <- function(...) {
    if (.lee_l_native_disabled()) {
        stop("lee_L_cols native backend disabled by option; using R reference backend.", call. = FALSE)
    }
    lee_L_cols(...)
}
.lee_perm <- function(...) {
    if (.lee_l_native_disabled()) {
        stop("lee_perm native backend disabled by option; permutation path unavailable.", call. = FALSE)
    }
    lee_perm(...)
}
.lee_perm_block <- function(...) {
    if (.lee_l_native_disabled()) {
        stop("lee_perm_block native backend disabled by option; permutation path unavailable.", call. = FALSE)
    }
    lee_perm_block(...)
}

.idelta_sparse_cpp <- function(...) idelta_sparse_cpp(...)

.pearson_native_disabled <- function() {
    isTRUE(getOption("geneSCOPE.disable_native_correlation_backend", FALSE))
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
