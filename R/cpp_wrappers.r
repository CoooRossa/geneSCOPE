#' @title Compute Grid Neighbors
#' @description
#'   Creates a spatial neighborhood list for a rectangular grid using either
#'   Queen (8-neighbor) or Rook (4-neighbor) connectivity.
#' @param nrow Integer. Number of rows in the grid
#' @param ncol Integer. Number of columns in the grid
#' @param queen Logical. TRUE for Queen moves, FALSE for Rook
#' @return A list of integer vectors representing neighbors
#' @export
grid_nb <- function(nrow, ncol, queen = TRUE) {
    .Call(`_geneSCOPE_grid_nb`, nrow, ncol, queen, PACKAGE = "geneSCOPE")
}

#' @title Compute Grid Neighbors with Parallelism
#' @description
#'   Creates a spatial neighborhood list for a rectangular grid using parallel processing.
#' @param nrow Integer. Number of rows in the grid
#' @param ncol Integer. Number of columns in the grid
#' @param queen Logical. TRUE for Queen moves, FALSE for Rook
#' @return A list of integer vectors representing neighbors
#' @export
grid_nb_omp <- function(nrow, ncol, queen = TRUE) {
    .Call(`_geneSCOPE_grid_nb_omp`, nrow, ncol, queen, PACKAGE = "geneSCOPE")
}

#' @title Convert Neighbor List to Sparse Matrix
#' @description
#'   Transforms a neighborhood list into a binary sparse matrix (dgCMatrix).
#' @param nb Neighborhood list object
#' @return A dgCMatrix sparse matrix
#' @export
nb2mat <- function(nb) {
    .Call(`_geneSCOPE_nb2mat`, nb, PACKAGE = "geneSCOPE")
}

#' @title Compute Lee's L Statistic
#' @description
#'   Computes Lee's L spatial correlation statistic for gene pairs.
#' @param Xz Numeric matrix of z-scored expression (rows=cells, cols=genes)
#' @param W Sparse weight matrix (dgCMatrix)
#' @param n_threads Integer. Number of OpenMP threads
#' @return Dense matrix of Lee's L values
#' @export
lee_L <- function(Xz, W, n_threads = 1) {
    .Call(`_geneSCOPE_lee_L`, Xz, W, n_threads, PACKAGE = "geneSCOPE")
}

#' @title Compute Lee's L for Column Subset
#' @description
#'   Computes Lee's L between specified columns and all genes.
#' @param Xz Expression matrix
#' @param W Weight matrix
#' @param cols0 Zero-based column indices
#' @param n_threads Number of threads
#' @return Matrix of Lee's L values
#' @export
lee_L_cols <- function(Xz, W, cols0, n_threads = 1) {
    .Call(`_geneSCOPE_lee_L_cols`, Xz, W, cols0, n_threads, PACKAGE = "geneSCOPE")
}

#' @title Compute Lee's L with Caching
#' @description
#'   Same as lee_L but caches spatially lagged matrix for efficiency.
#' @param Xz n × g numeric matrix of z-scored gene expression (rows = cells)
#' @param W n × n sparse weight matrix in dgCMatrix format
#' @param n_threads Integer. Number of OpenMP threads (default 1)
#' @return Dense matrix of Lee's L values
#' @export
lee_L_cache <- function(Xz, W, n_threads = 1) {
    .Call(`_geneSCOPE_lee_L_cache`, Xz, W, n_threads, PACKAGE = "geneSCOPE")
}

#' @title Monte-Carlo Permutation for Lee's L
#' @description
#'   Performs permutation tests for Lee's L statistics.
#' @param Xz Expression matrix
#' @param W Weight matrix
#' @param idx_mat Permutation index matrix
#' @param L_ref Reference Lee's L matrix
#' @param n_threads Number of threads
#' @return Matrix of exceedance counts
#' @export
lee_perm <- function(Xz, W, idx_mat, L_ref, n_threads = 1) {
    .Call(`_geneSCOPE_lee_perm`, Xz, W, idx_mat, L_ref, n_threads, PACKAGE = "geneSCOPE")
}

#' @title Block-wise Permutation for Lee's L
#' @description
#'   Performs block-constrained permutation tests.
#' @param Xz Expression matrix
#' @param W Weight matrix
#' @param idx_mat Permutation index matrix
#' @param block_ids Block identifier vector
#' @param L_ref Reference Lee's L matrix
#' @param n_threads Number of threads
#' @return Matrix of exceedance counts
#' @export
lee_perm_block <- function(Xz, W, idx_mat, block_ids, L_ref, n_threads = 1) {
    .Call(`_geneSCOPE_lee_perm_block`, Xz, W, idx_mat, block_ids, L_ref, n_threads, PACKAGE = "geneSCOPE")
}

#' @title Compute Morisita's Index Delta
#' @description
#'   Computes Morisita's index of dispersion for sparse count matrix.
#' @param G Sparse count matrix (genes x grids)
#' @param n_threads Number of threads
#' @return Vector of Morisita's delta values
#' @export
idelta <- function(G, n_threads = 1) {
    .Call(`_geneSCOPE_idelta`, G, n_threads, PACKAGE = "geneSCOPE")
}

#' @title Compute Pearson Correlation
#' @description
#'   Computes block-wise Pearson correlation matrix in parallel.
#' @param X Numeric matrix (already centered)
#' @param bs Block size for processing
#' @param n_threads Number of threads
#' @return Correlation matrix
#' @export
pearson_cor <- function(X, bs = 2000, n_threads = 1) {
    .Call(`_geneSCOPE_pearson_cor`, X, bs, n_threads, PACKAGE = "geneSCOPE")
}

#' @title LOESS Bootstrap
#' @description
#'   Performs bootstrap LOESS regression with confidence intervals.
#' @param x Predictor vector
#' @param y Response vector
#' @param strata Stratification vector
#' @param xgrid Grid points for prediction
#' @param B Number of bootstrap samples
#' @param span LOESS span parameter
#' @param deg Polynomial degree
#' @param n_threads Number of threads
#' @param k_max Maximum number of iterations
#' @param keep_boot Logical. Keep bootstrap samples
#' @param adjust_mode Adjustment mode
#' @param ci_type Confidence interval type
#' @param level Confidence level
#' @return List with fit, confidence intervals
#' @export
loess_residual_bootstrap <- function(x, y, strat, grid,
                                     B = 1000L,
                                     span = 0.45,
                                     deg = 1L,
                                     n_threads = 1L,
                                     k_max = -1L,
                                     keep_boot = TRUE,
                                     adjust_mode = 0L,
                                     ci_type = 0L,
                                     level = 0.95) {
    .Call(`_geneSCOPE_loess_residual_bootstrap`,
        x, y, strat, grid,
        as.integer(B),
        span,
        as.integer(deg),
        as.integer(n_threads),
        as.integer(k_max),
        as.logical(keep_boot),
        as.integer(adjust_mode),
        as.integer(ci_type),
        level,
        PACKAGE = "geneSCOPE"
    )
}

#' @title LOESS Bootstrap (Compatibility Wrapper)
#' @description
#'   Compatibility wrapper for LOESS bootstrap with old parameter names.
#' @param x Predictor vector
#' @param y Response vector
#' @param strata Stratification vector
#' @param xgrid Grid points for prediction
#' @param B Number of bootstrap samples
#' @param span LOESS span parameter
#' @param degree Polynomial degree
#' @param nthreads Number of threads
#' @param ... Additional arguments
#' @return List with fit, confidence intervals
#' @export
loess_residual_bootstrap_compat <- function(x, y, strat, grid,
                                            B = 1000L, span = 0.45,
                                            degree = 1L, nthreads = 1L,
                                            ...) {
    loess_residual_bootstrap(x, y, strat, grid,
        B = B, span = span,
        deg = degree,
        n_threads = nthreads, 
    )
}

# Unified thin wrappers (for future logging or checks)

lr_bootstrap <- function(x, y, strat, grid,
                         B = 1000L, span = 0.45, deg = 1L,
                         n_threads = 1L, k_max = -1L,
                         keep_boot = TRUE, adjust_mode = 0L,
                         ci_type = 0L, level = 0.95) {
    stopifnot(
        is.numeric(x), is.numeric(y), length(x) == length(y),
        length(strat) == length(x), is.numeric(grid)
    )
    loess_residual_bootstrap(x, y, strat, grid,
        B = B, span = span, deg = deg,
        n_threads = n_threads, k_max = k_max,
        keep_boot = keep_boot, adjust_mode = adjust_mode,
        ci_type = ci_type, level = level
    )
}

leeL_full <- function(Xz, W, n_threads = 1L, cache = TRUE) {
    if (cache) lee_L_cache(Xz, W, n_threads) else lee_L(Xz, W, n_threads)
}

leeL_subset <- function(Xz, W, cols0, n_threads = 1L) {
    lee_L_cols(Xz, W, cols0, n_threads)
}

safe_thread_count <- function() {
    # Base heuristic first
    core_guess <- tryCatch(parallel::detectCores(), error = function(e) 1L)
    n_req <- max(1L, min(8L, as.integer(core_guess) - 1L))

    if (exists("getSafeThreadCount", mode = "function", inherits = TRUE)) {
        f <- get("getSafeThreadCount", mode = "function")
        an <- tryCatch(names(formals(f)), error = function(e) character(0))
        # Support both old and new signatures
        if ("max_requested" %in% an) {
            return(tryCatch(f(max_requested = n_req), error = function(e) n_req))
        } else if ("requested_cores" %in% an) {
            return(tryCatch(f(requested_cores = n_req), error = function(e) n_req))
        } else {
            return(tryCatch(f(n_req), error = function(e) n_req))
        }
    }
    n_req
}

#' @title Consensus co-occurrence (sparse COO via C++)
#' @description
#'   Works in two modes:
#'   - Package mode: calls registered symbol `_geneSCOPE_consensus_coo_cpp`.
#'   - Sourced mode: if you ran `Rcpp::sourceCpp('src/9.ConsensusAccel.cpp')`,
#'     calls the Rcpp-exposed function `consensus_coo_cpp()` directly.
#' @param memb Integer matrix (genes × runs), each column a run's cluster labels
#' @param thr  Numeric in [0,1]. Threshold on fraction across runs (default 0)
#' @param n_threads Integer. Threads (ignored if OpenMP not available)
#' @return list(i,j,x,n) with 1-based COO
#' @export
consensus_coo <- function(memb, thr = 0, n_threads = NULL) {
    if (is.null(n_threads)) n_threads <- safe_thread_count()
    # Try package-registered symbol
    call_sym_ok <- FALSE
    res <- try({
        .Call(`_geneSCOPE_consensus_coo_cpp`, memb, thr, as.integer(n_threads), PACKAGE = "geneSCOPE")
    }, silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    # Try direct Rcpp function if sourceCpp() was used
    if (exists("consensus_coo_cpp", mode = "function", inherits = TRUE)) {
        return(consensus_coo_cpp(memb, thr, as.integer(n_threads)))
    }
    stop("consensus_coo: C++ backend not found. Install/load geneSCOPE or run Rcpp::sourceCpp('genescope/src/9.ConsensusAccel.cpp').")
}

#' @title Consensus sparse matrix (dgCMatrix)
#' @description
#'   Convenience wrapper building a symmetric sparse matrix with edge weights
#'   equal to fraction of co-occurrence across runs for pairs meeting `thr`.
#' @inheritParams consensus_coo
#' @return A symmetric `dgCMatrix` with zero diagonal
#' @export
consensus_sparse <- function(memb, thr = 0, n_threads = NULL) {
    if (is.null(n_threads)) n_threads <- safe_thread_count()
    coo <- consensus_coo(memb, thr = thr, n_threads = n_threads)
    # Build symmetric sparse matrix
    M <- if (length(coo$i) == 0L) {
        Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(coo$n, coo$n), symmetric = TRUE)
    } else {
        Matrix::sparseMatrix(i = coo$i, j = coo$j, x = coo$x, dims = c(coo$n, coo$n), symmetric = TRUE)
    }
    # Propagate dimnames from membership matrix when available to enable name-based indexing
    rn <- try(if (is.matrix(memb) || inherits(memb, "Matrix")) rownames(memb) else rownames(as.matrix(memb)), silent = TRUE)
    if (!inherits(rn, "try-error") && !is.null(rn) && length(rn) == nrow(M)) {
        # Guard against rare cases where dimnames<- is applied to non-arrays
        try({
            if (is.matrix(M) || inherits(M, "Matrix")) dimnames(M) <- list(rn, rn)
        }, silent = TRUE)
    }
    M
}

delta_perm_pairs <- function(Xz, W, idx_mat, gene_pairs, delta_ref,
                             block_ids = NULL, csr = NULL,
                             n_threads = 1L,
                             chunk_size = 1000L,
                             clamp_nonneg_r = FALSE,
                             tiny = FALSE) {
    if (!is.null(csr)) {
        if (is.null(block_ids)) {
            delta_lr_perm_csr(Xz, csr$indices, csr$values, csr$row_ptr,
                idx_mat, gene_pairs, delta_ref,
                n_threads = n_threads, clamp_nonneg_r = clamp_nonneg_r
            )
        } else {
            delta_lr_perm_csr_block(Xz, csr$indices, csr$values, csr$row_ptr,
                idx_mat, block_ids, gene_pairs, delta_ref,
                n_threads = n_threads, clamp_nonneg_r = clamp_nonneg_r
            )
        }
    } else if (tiny) {
        if (is.null(block_ids)) {
            delta_lr_perm_tiny(Xz, W, idx_mat, gene_pairs, delta_ref, n_threads)
        } else {
            delta_lr_perm_block_tiny(Xz, W, idx_mat, block_ids, gene_pairs, delta_ref, n_threads)
        }
    } else {
        if (is.null(block_ids)) {
            delta_lr_perm(Xz, W, idx_mat, gene_pairs, delta_ref,
                n_threads = n_threads, chunk_size = chunk_size
            )
        } else {
            delta_lr_perm_block(Xz, W, idx_mat, block_ids, gene_pairs, delta_ref,
                n_threads = n_threads, chunk_size = chunk_size
            )
        }
    }
}
# Thin R wrappers for Visium-specific C++ accelerators

#' Sparse × dense multiply: WX = W %*% X
#' @param X n×S numeric matrix
#' @param W n×n dgCMatrix
#' @param n_threads Integer, OpenMP threads
#' @param precision "float32" or "float64"
#' @keywords internal
spmm_dgc_dense <- function(X, W, n_threads = 1L, precision = c("float32", "float64")) {
    precision <- match.arg(precision)
    if (precision == "float32") {
        .Call(`_geneSCOPE_spmm_dgc_dense_f32`, X, W, as.integer(n_threads), PACKAGE = "geneSCOPE")
    } else {
        .Call(`_geneSCOPE_spmm_dgc_dense_f64`, X, W, as.integer(n_threads), PACKAGE = "geneSCOPE")
    }
}

#' Sign random projection LSH: per-gene bucket IDs (hex) per table
#' @param X n×S numeric matrix
#' @param bits Bits per table b (<= 60)
#' @param n_tables Number of tables L
#' @param seed Random seed
#' @param n_threads Threads
#' @return Character matrix S×L, each entry like "0xABCD..." (64-bit bucket)
#' @keywords internal
rp_sign_bits <- function(X, bits = 12L, n_tables = 6L, seed = 1L, n_threads = 1L) {
    .Call(`_geneSCOPE_rp_sign_bits`, X, as.integer(bits), as.integer(n_tables), as.integer(seed), as.integer(n_threads), PACKAGE = "geneSCOPE")
}

#' Compute L on candidate CSR and keep per-row Top-K
#' @param X n×S matrix
#' @param WX n×S matrix (W%*%X)
#' @param row_ptr Integer vector length S+1 (1-based CSR start/end)
#' @param indices Integer vector (1-based column indices)
#' @param K_keep Keep Top-K per row
#' @param n_threads Threads
#' @return list(row_ptr, indices, values)
#' @keywords internal
leeL_topk_candidates <- function(X, WX, row_ptr, indices, K_keep = 100L, n_threads = 1L) {
    .Call(`_geneSCOPE_leeL_topk_candidates`, X, WX, as.integer(row_ptr), as.integer(indices), as.integer(K_keep), as.integer(n_threads), PACKAGE = "geneSCOPE")
}
