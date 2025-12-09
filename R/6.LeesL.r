#' @title Add Lee's L statistics to a scope_object
#' @description
#'   High-level wrapper that computes Lee's L, analytical Z-scores, empirical
#'   p-values via block permutations, FDR, spatial gradients, quality-control
#'   metrics, and stores everything under a new layer in \code{@grid}.
#' @param scope_obj A \code{scope_object} with at least one populated \code{@grid} slot.
#' @param grid_name Character. Name of the grid sub-layer to process. If
#'   \code{NULL} and only one sub-layer exists, it is selected automatically.
#' @param genes Optional character vector of genes to include; if \code{NULL}
#'   all genes are used.
#' @param within Logical. If \code{TRUE} restrict analysis to the selected
#'   gene set on both axes (default). Otherwise compute gene × all.
#' @param ncores Integer. Number of cores for parallel processing (default 1).
#' @param block_side Integer. Number of grid cells per side for block partitioning (default 8).
#' @param perms Integer. Number of permutations for Monte-Carlo p-values (default 1000).
#' @param block_size Integer. Number of permutations processed per batch (default 64).
#' @param L_min Numeric threshold used when building the QC similarity graph.
#' @param norm_layer Character. Name of the normalised expression layer (default "Xz").
#' @param lee_stats_layer_name Character. Name for the output statistics layer.
#' @param legacy_formula Logical. Use legacy denominator for compatibility.
#' @param mem_limit_GB Numeric. RAM threshold that triggers streaming mode (default 8).
#' @param chunk_size Integer. Number of columns processed per chunk (default 128).
#' @param use_bigmemory Logical. Use file-backed matrices for large computations (default TRUE).
#' @param backing_path Character. Directory for temporary files (default tempdir()).
#' @param cache_inputs Logical. Cache preprocessed X/Z/W and block IDs for reuse across calls (default TRUE).
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#' @param ncore Deprecated. Use \code{ncores} instead.
#' @return The modified \code{scope_object}.
#' @importFrom RhpcBLASctl blas_set_num_threads
#' @importFrom stats sd pnorm p.adjust coef lm
#' @importFrom igraph graph_from_adjacency_matrix simplify degree cluster_louvain cluster_leiden components modularity
#' @export
computeL <- function(scope_obj,
                     grid_name = NULL,
                     genes = NULL,
                     within = TRUE,
                     ncores = 1,
                     block_side = 8,
                     perms = 1000,
                     block_size = 64,
                     L_min = 0,
                     norm_layer = "Xz",
                     lee_stats_layer_name = NULL,
                     legacy_formula = FALSE,
                     mem_limit_GB = 2,
                     chunk_size = 32L,
                     use_bigmemory = TRUE,
                     backing_path = tempdir(),
                     cache_inputs = TRUE,
                     verbose = TRUE,
                     ncore = NULL) {
    ## --- 0. Thread-safe preprocessing: automatic thread management and error recovery ---

    # Handle deprecated parameter
    if (!is.null(ncore)) {
        warning("'ncore' is deprecated, please use 'ncores' instead. Will use its value this time.",
            call. = FALSE, immediate. = TRUE
        )
        ncores <- ncore
    }

    user_set_bigmemory <- !missing(use_bigmemory)
    user_set_chunk <- !missing(chunk_size)

    # Get system information and aggressive thread count (use all visible cores up to request)
    os_type <- detectOS()
    avail_cores <- max(1L, parallel::detectCores(logical = TRUE))
    ncores <- max(1L, min(ncores, avail_cores))
    if (verbose) {
        message("[geneSCOPE::computeL] Core configuration: using ", ncores, "/", avail_cores, " logical cores")
    }

    # Use thread configuration for mixed OpenMP/BLAS operations
    thread_config <- configureThreadsFor("mixed", ncores, restore_after = TRUE)
    on.exit({
        restore_fn <- attr(thread_config, "restore_function")
        if (!is.null(restore_fn)) restore_fn()
    })

    # Use OpenMP threads for C++ operations
    ncores_cpp <- thread_config$openmp_threads

    ## --- 1. Check if spatial weights exist, auto-compute if missing ---
    g_layer <- .selectGridLayer(scope_obj, grid_name)
    grid_name <- if (is.null(grid_name)) {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    } else {
        grid_name
    }

    # Check if spatial weight matrix exists
    if (is.null(g_layer$W) || sum(g_layer$W) == 0) {
        message("[geneSCOPE::computeL] Computing spatial weights matrix")
        scope_obj <- computeWeights(scope_obj,
            grid_name = grid_name,
            style = "B",
            store_mat = TRUE,
            store_listw = FALSE
        )
        # Re-extract layer after computing weights
        g_layer <- .selectGridLayer(scope_obj, grid_name)
    }

    # Cross-platform memory calculation
    all_genes <- if (!is.null(g_layer[[norm_layer]])) {
        ncol(g_layer[[norm_layer]])
    } else {
        stop("Normalized layer '", norm_layer, "' not found")
    }

    n_genes_use <- if (!is.null(genes)) length(genes) else all_genes
    matrix_size_gb <- (n_genes_use^2 * 8) / (1024^3)

    # Respect user bigmemory/chunk overrides while still hinting when streaming is advisable
    if (matrix_size_gb > mem_limit_GB) {
        if (use_bigmemory) {
            if (verbose) {
                message("[geneSCOPE::computeL] Large matrix detected (", round(matrix_size_gb, 1), " GB); staying in bigmemory/streaming mode.",
                    if (!user_set_chunk) " You may tune chunk_size to trade IO vs RAM." else "")
            }
        } else if (user_set_bigmemory) {
            if (verbose) message("[geneSCOPE::computeL] Warning: requested use_bigmemory=FALSE with large matrix (", round(matrix_size_gb, 1), " GB); proceeding in-memory as requested.")
        } else if (os_type == "windows" && !requireNamespace("bigmemory", quietly = TRUE)) {
            if (verbose) message("[geneSCOPE::computeL] !!! Warning: bigmemory not available on Windows; using regular matrices !!!")
            use_bigmemory <- FALSE
        } else {
            use_bigmemory <- TRUE
            if (verbose) message("[geneSCOPE::computeL] Large matrix detected (", round(matrix_size_gb, 1), " GB); enabling bigmemory/streaming.")
        }
    }

    # Cross-platform temporary directory
    if (use_bigmemory) {
        if (backing_path == tempdir()) {
            # Use platform-specific temporary directory
            backing_path <- switch(os_type,
                windows = Sys.getenv("TEMP", tempdir()),
                linux = Sys.getenv("TMPDIR", "/tmp"),
                macos = Sys.getenv("TMPDIR", "/tmp")
            )
        }

        if (!dir.exists(backing_path)) {
            dir.create(backing_path, recursive = TRUE)
        }
    }

    ## --- 2. Lee's L computation with error recovery ---
    if (verbose) message("[geneSCOPE::computeL] Computing Lee's L statistics")

    # Multi-retry mechanism with reduced thread count
    current_cores <- ncores
    min_cores <- 1
    success <- FALSE
    attempt <- 1

    while (!success && current_cores >= min_cores) {
        if (verbose && attempt > 1) message("[geneSCOPE::computeL]  Retry #", attempt, " with ", current_cores, " cores")

        result <- tryCatch(
            {
                # Execute computation
                t_start <- Sys.time()
                res <- .computeLeeL(scope_obj,
                    grid_name = grid_name,
                    norm_layer = norm_layer,
                    genes = genes,
                    within = within,
                    ncore = current_cores, # Use current core count
                    mem_limit_GB = mem_limit_GB,
                    chunk_size = chunk_size,
                    use_bigmemory = use_bigmemory,
                    backing_path = backing_path,
                    block_side = block_side,
                    cache_inputs = cache_inputs,
                    input_cache = if (cache_inputs && !is.null(scope_obj@stats[[grid_name]])) {
                        scope_obj@stats[[grid_name]][["_lee_input_cache"]]
                    } else {
                        NULL
                    }
                )
                t_end <- Sys.time()

                if (verbose) {
                    time_msg <- if (attempt == 1) "completed" else "retry successful"
                    message("[geneSCOPE::computeL]  Lee's L ", time_msg, " (", format(t_end - t_start), ")")
                }

                list(success = TRUE, object = res)
            },
            error = function(e) {
                if (verbose && attempt > 1) message("[geneSCOPE::computeL]  Attempt failed: ", conditionMessage(e))
                list(success = FALSE, error = e)
            }
        )

        # Update status
        success <- result$success

        if (success) {
            res <- result$object
        } else {
            attempt <- attempt + 1
            # Reduce thread count
            current_cores <- max(floor(current_cores / 2), min_cores)
            if (verbose) message("[geneSCOPE::computeL]  Reducing cores to ", current_cores, " and retrying")
            # Give system some time to recover
            Sys.sleep(3)
            # Force garbage collection
            gc(verbose = FALSE)
        }
    }

    if (!success) {
        stop("Unable to compute Lee's L statistics even with minimal thread count")
    }

    L <- res$Lmat
    X_full <- res$X_full
    X_used <- res$X_used
    W <- res$W
    grid_inf <- res$grid_info
    gname <- res$grid_name
    n <- nrow(X_full)

    ## --- 3. Analytical Z computation with memory optimization ---
    if (within || is.null(genes)) {
        S0 <- sum(W)
        EZ <- -1 / (n - 1)
        Var <- (n^2 * (n - 2)) / ((n - 1)^2 * (n - 3) * S0)

        # Handle large matrices more carefully
        if (inherits(L, "big.matrix")) {
            if (verbose) message("[geneSCOPE::computeL]  Computing Z-scores in chunks")
            Z_mat <- L # Use the same backing store
            # Process in chunks to avoid memory issues
            n_genes <- ncol(L)
            chunk_genes <- seq(1, n_genes, by = chunk_size)
            for (start in chunk_genes) {
                end <- min(start + chunk_size - 1, n_genes)
                L_chunk <- L[, start:end]
                Z_chunk <- (L_chunk - EZ) / sqrt(Var)
                Z_mat[, start:end] <- Z_chunk
            }
        } else {
            Z_mat <- (L - EZ) / sqrt(Var)
        }
    } else {
        if (inherits(L, "big.matrix")) {
            stop("Asymmetric Z-score computation not yet supported for big.matrix")
        }
        Z_mat <- t(apply(L, 1, function(v) {
            sdv <- stats::sd(v)
            if (sdv == 0) rep(0, length(v)) else (v - mean(v)) / sdv
        }))
        dimnames(Z_mat) <- dimnames(L)
    }

    block_id <- res$block_id
    input_cache <- res$input_cache

    ## --- 4. Monte Carlo p-values with BLAS control and error recovery ---
    P <- if (perms > 0) {
        # Temporarily disable BLAS threads for permutation tests
        if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
            RhpcBLASctl::blas_set_num_threads(1)
        }

        if (inherits(L, "big.matrix")) {
            if (verbose) message("[geneSCOPE::computeL]  Converting to regular matrix for permutation tests")
            L_reg <- as.matrix(L)
        } else {
            L_reg <- L
        }

        # Permutation test with error recovery
        if (verbose) message("[geneSCOPE::computeL] Running permutation tests")
        perm_success <- FALSE
        perm_cores <- current_cores
        perm_attempt <- 1

        while (!perm_success && perm_cores >= 1) {
            if (verbose && perm_attempt > 1) message("[geneSCOPE::computeL]  Retry #", perm_attempt, " with ", perm_cores, " cores")

            # Sanity check to avoid C++ OOB when a gene subset is used
            if (ncol(X_used) != nrow(L_reg) || ncol(X_used) != ncol(L_reg)) {
                stop(sprintf("Permutation input mismatch: ncol(X_used)=%d but L_ref is %dx%d",
                    ncol(X_used), nrow(L_reg), ncol(L_reg)))
            }

            perm_result <- tryCatch(
                {
                    t_start <- Sys.time()
                    # Use X_used to match the dimensions of L_reg when a gene subset is used
                    p_result <- .leeL_perm_block(X_used, W, L_reg,
                        block_id   = block_id,
                        perms      = perms,
                        block_size = block_size,
                        n_threads  = perm_cores
                    )
                    t_end <- Sys.time()
                    if (verbose) {
                        time_msg <- if (perm_attempt == 1) "completed" else "retry successful"
                        message("[geneSCOPE::computeL]  Permutation test ", time_msg, " (", format(t_end - t_start), ")")
                    }
                    list(success = TRUE, object = p_result)
                },
                error = function(e) {
                    if (verbose && perm_attempt > 1) message("[geneSCOPE::computeL]  Test failed: ", conditionMessage(e))
                    list(success = FALSE, error = e)
                }
            )

            if (perm_result$success) {
                perm_success <- TRUE
                P <- perm_result$object
            } else {
                perm_attempt <- perm_attempt + 1
                perm_cores <- max(floor(perm_cores / 2), 1)
                if (verbose) message("[geneSCOPE::computeL]  Reducing cores to ", perm_cores, " and retrying")
                Sys.sleep(2)
                gc(verbose = FALSE)
            }
        }

        if (!perm_success) {
            if (verbose) message("[geneSCOPE::computeL] !!! Warning: Permutation test failed; setting P = NULL !!!")
            NULL
        } else {
            P
        }
    } else {
        NULL
    }

    ## --- 5. FDR: Three Monte Carlo smoothing strategies as in getTopLvsR ---
    # Explanation:
    #   If P exists: P = (k+1)/(perms+1); infer k = round(P*(perms+1)-1)
    #   p_beta    = (k+1)/(perms+2)        (Jeffreys, default more robust/slightly conservative)
    #   p_mid     = (k+0.5)/(perms+1)      (mid-p, less conservative)
    #   p_uniform = (k+U)/(perms+1)        (random jitter, expected equal to discrete grid, breaks staircase)
    #   Main output FDR = BH(p_beta)
    if (inherits(Z_mat, "big.matrix")) {
        if (verbose) message("[geneSCOPE::computeL] Computing FDR corrections")
        n_genes <- ncol(Z_mat)
        idx_chunks <- seq(1, n_genes, by = chunk_size)

        FDR_disc <- bigmemory::big.matrix(
            nrow = nrow(Z_mat), ncol = n_genes,
            type = "double", init = NA_real_
        )
        FDR_beta <- bigmemory::big.matrix(
            nrow = nrow(Z_mat), ncol = n_genes,
            type = "double", init = NA_real_
        )
        FDR_mid <- bigmemory::big.matrix(
            nrow = nrow(Z_mat), ncol = n_genes,
            type = "double", init = NA_real_
        )
        FDR_uniform <- bigmemory::big.matrix(
            nrow = nrow(Z_mat), ncol = n_genes,
            type = "double", init = NA_real_
        )

        for (start in idx_chunks) {
            end <- min(start + chunk_size - 1L, n_genes)
            if (!is.null(P)) {
                p_disc <- P[, start:end]
                k_est <- round(p_disc * (perms + 1) - 1)
                k_est[k_est < 0] <- 0
                k_est[k_est > perms] <- perms
                p_beta <- (k_est + 1) / (perms + 2)
                p_mid <- (k_est + 0.5) / (perms + 1)
                p_uniform <- (k_est + matrix(runif(length(k_est)), nrow = nrow(k_est))) / (perms + 1)
            } else {
                # No permutation → analytical p
                Z_chunk <- Z_mat[, start:end]
                p_disc <- 2 * pnorm(-abs(Z_chunk))
                p_beta <- p_mid <- p_uniform <- p_disc
            }
            FDR_disc[, start:end] <- matrix(p.adjust(c(p_disc), "BH"), nrow = nrow(p_disc))
            FDR_beta[, start:end] <- matrix(p.adjust(c(p_beta), "BH"), nrow = nrow(p_beta))
            FDR_mid[, start:end] <- matrix(p.adjust(c(p_mid), "BH"), nrow = nrow(p_mid))
            FDR_uniform[, start:end] <- matrix(p.adjust(c(p_uniform), "BH"), nrow = nrow(p_uniform))
        }
        FDR_out_disc <- FDR_disc
        FDR_out_beta <- FDR_beta
        FDR_out_mid <- FDR_mid
        FDR_out_uniform <- FDR_uniform
    } else {
        # In-memory mode: compute in one go (exact BH)
        p_disc <- if (!is.null(P)) P else 2 * pnorm(-abs(Z_mat))
        if (!is.null(P)) {
            k_est <- round(p_disc * (perms + 1) - 1)
            k_est[k_est < 0] <- 0
            k_est[k_est > perms] <- perms
            p_beta <- (k_est + 1) / (perms + 2)
            p_mid <- (k_est + 0.5) / (perms + 1)
            p_uniform <- (k_est + runif(length(k_est))) / (perms + 1)
        } else {
            p_beta <- p_mid <- p_uniform <- p_disc
        }
        # BH on flattened vectors, then reshape
        shp <- dim(p_disc)
        FDR_out_disc <- matrix(p.adjust(p_disc, "BH"), nrow = shp[1], dimnames = dimnames(p_disc))
        FDR_out_beta <- matrix(p.adjust(p_beta, "BH"), nrow = shp[1], dimnames = dimnames(p_disc))
        FDR_out_mid <- matrix(p.adjust(p_mid, "BH"), nrow = shp[1], dimnames = dimnames(p_disc))
        FDR_out_uniform <- matrix(p.adjust(p_uniform, "BH"), nrow = shp[1], dimnames = dimnames(p_disc))
    }

    ## --- 5.x Extra: Storey q-value main FDR (in-memory mode only, big.matrix fallback) ---
    FDR_storey <- NULL
    pi0_hat <- NA_real_
    FDR_main_method <- if (inherits(Z_mat, "big.matrix")) "BH(beta p) (chunked, no q-value)" else "Storey q-value (beta p)"
    if (!inherits(Z_mat, "big.matrix")) {
        # Flatten to vector
        p_vec <- as.numeric(if (!is.null(P)) {
            # p_beta already computed above
            p_beta
        } else {
            p_disc
        })
        m_tot <- length(p_vec)
        if (m_tot > 0) {
            # 1) Estimate pi0
            lambdas <- seq(0.5, 0.95, by = 0.05)
            pi0_vals <- sapply(lambdas, function(lam) {
                mean(p_vec > lam) / (1 - lam)
            })
            pi0_hat <- min(1, min(pi0_vals, na.rm = TRUE))
            if (!is.finite(pi0_hat) || pi0_hat <= 0) pi0_hat <- 1
            # 2) Compute raw q
            o <- order(p_vec, na.last = NA)
            ro <- integer(m_tot)
            ro[o] <- seq_along(o)
            q_raw <- pi0_hat * m_tot * p_vec / pmax(ro, 1)
            # 3) Monotonic decreasing adjustment (tail to head)
            q_ord <- q_raw[o]
            for (i in seq_along(q_ord)[(length(q_ord) - 1):1]) {
                if (q_ord[i] > q_ord[i + 1]) q_ord[i] <- q_ord[i + 1]
            }
            q_final <- numeric(m_tot)
            q_final[o] <- pmin(q_ord, 1)
            # 4) Restore matrix
            FDR_storey <- matrix(q_final, nrow = nrow(FDR_out_beta), dimnames = dimnames(FDR_out_beta))
        }
    }
    # Main FDR choice: use FDR_storey if available, otherwise fallback to FDR_out_beta
    FDR_main <- if (!is.null(FDR_storey)) FDR_storey else FDR_out_beta
    n_sig_005 <- sum(FDR_main < 0.05, na.rm = TRUE)
    min_p_possible <- if (!is.null(P)) 1 / (perms + 1) else NA_real_

    ## --- 6. βx / βy ---
    centres <- with(
        grid_inf[match(res$cells, grid_inf$grid_id), ],
        data.frame(
            x = (xmin + xmax) / 2,
            y = (ymin + ymax) / 2
        )
    )
    betas <- t(apply(X_full, 2, function(v) {
        c(beta_x = coef(lm(v ~ centres$x + centres$y))[2:3])
    }))
    if (!is.null(genes)) betas <- betas[genes, , drop = FALSE]

    ## --- 7. QC computation ---
    qc <- NULL
    if (within || is.null(genes)) {
        # Convert to regular matrix for QC if needed
        L_qc <- if (inherits(L, "big.matrix")) {
            message("Converting subset of big.matrix for QC computation...")
            # Sample subset for QC to avoid memory issues
            n_sample <- min(1000, nrow(L))
            sample_idx <- sample(nrow(L), n_sample)
            as.matrix(L[sample_idx, sample_idx])
        } else {
            L
        }

        A_bin <- abs(L_qc) >= L_min
        A_bin <- A_bin | t(A_bin)
        diag(A_bin) <- FALSE

        A_num <- (A_bin | t(A_bin)) * 1
        g_tmp <- igraph::simplify(
            igraph::graph_from_adjacency_matrix(A_num,
                mode = "undirected",
                diag = FALSE
            ),
            remove.multiple = TRUE,
            remove.loops = TRUE
        )

        deg <- igraph::degree(g_tmp)
        memb <- tryCatch(
            igraph::cluster_louvain(g_tmp)$membership,
            error = function(e) {
                igraph::cluster_leiden(g_tmp)$membership
            }
        )

        qc <- list(
            edge_density = 2 * sum(A_bin[upper.tri(A_bin)]) /
                (ncol(A_bin) * (ncol(A_bin) - 1)),
            components = igraph::components(g_tmp)$no,
            modularity_Q = igraph::modularity(g_tmp, membership = memb),
            mean_degree = mean(deg),
            sd_degree = stats::sd(deg),
            hub_ratio = mean(deg > 2 * median(deg)),
            sig_edge_frac = if (is.null(P)) {
                NA
            } else {
                if (inherits(P, "big.matrix")) {
                    # Sample for significance fraction
                    p_sample <- as.matrix(P[sample_idx, sample_idx])
                    mean(p_sample[upper.tri(p_sample)] < 0.05, na.rm = TRUE)
                } else {
                    mean(P[upper.tri(P)] < 0.05, na.rm = TRUE)
                }
            }
        )
    }

    ## --- 8. Write to stats (main FDR replacement; retain old matrices; add FDR_storey) ---
    layer_name <- if (is.null(lee_stats_layer_name)) paste0("LeeStats_", norm_layer) else lee_stats_layer_name
    if (is.null(scope_obj@stats)) scope_obj@stats <- list()
    if (is.null(scope_obj@stats[[gname]])) scope_obj@stats[[gname]] <- list()

    scope_obj@stats[[gname]][[layer_name]] <- list(
        L = L,
        Z = Z_mat,
        P = P,
        grad = betas,
        L_min = L_min,
        qc = qc,
        FDR = FDR_main, # Main recommendation (q-value or beta-BH fallback)
        FDR_storey = FDR_storey, # Available in in-memory mode only
        FDR_disc = FDR_out_disc, # Original discrete BH
        FDR_beta = FDR_out_beta,
        FDR_mid = FDR_out_mid,
        FDR_uniform = FDR_out_uniform,
        meta = list(
            perms = perms,
            block_size = block_size,
            ncores = ncores,
            mem_mode = if (inherits(L, "big.matrix")) "bigmemory" else "inmemory",
            p_source = if (!is.null(P)) "permutation" else "analytic",
            p_resolution = if (!is.null(P)) paste0("~", format(1 / (perms + 1), scientific = TRUE)) else "analytic",
            pi0_hat = pi0_hat,
            n_tests = length(FDR_main),
            n_sig_FDR_lt_0.05 = n_sig_005,
            min_p_possible = min_p_possible,
            FDR_main_method = FDR_main_method,
            pval_modes = c("disc=(k+1)/(B+1)", "beta=(k+1)/(B+2)", "mid=(k+0.5)/(B+1)", "uniform=(k+U)/(B+1)"),
            note = "FDR: main=Storey q (beta p) if possible; fallback=BH(beta). All matrices retained for diagnostics."
        )
    )

    if (cache_inputs) {
        if (is.null(scope_obj@stats[[gname]])) scope_obj@stats[[gname]] <- list()
        scope_obj@stats[[gname]][["_lee_input_cache"]] <- input_cache
    }

    invisible(scope_obj)
}


#' @noRd
.computeLeeL <- function(scope_obj,
                         grid_name = NULL,
                         norm_layer = "Xz",
                         genes = NULL,
                         within = TRUE,
                         ncore = 1,
                         mem_limit_GB = 16,
                         chunk_size = 256L,
                         use_bigmemory = TRUE,
                         backing_path = tempdir(),
                         block_side = 8,
                         cache_inputs = FALSE,
                         input_cache = NULL) {
    if (use_bigmemory && !requireNamespace("bigmemory", quietly = TRUE)) {
        stop("Package 'bigmemory' is required.")
    }

    ## ---- 0. Get grid layer (new helper) ---------------------------------
    g_layer <- .selectGridLayer(scope_obj, grid_name)
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
        W <- g_layer$W[ord, ord]
        block_id <- .assign_block_id(grid_info, block_side = block_side)
    }

    block_id <- .compact_block_id(block_id)


    all_genes <- colnames(Xz_full)
    idx_keep <- if (is.null(genes)) {
        seq_along(all_genes)
    } else {
        m <- match(genes, all_genes, nomatch = 0L)
        if (any(m == 0L)) stop("Some genes were not found in column names")
        m
    }

    n_g <- length(idx_keep)
    bytes_L <- n_g^2 * 8 / 1024^3 # in GB
    need_stream <- use_bigmemory && bytes_L > mem_limit_GB

    ## ======================================================
    ## =============  A. One‑shot computation (fits in RAM) ===============
    ## ======================================================
    if (!need_stream) {
        RhpcBLASctl::blas_set_num_threads(1)
        Sys.setenv(OMP_NUM_THREADS = ncore)
        # Use correct export function name
        L_full <- lee_L_cache(Xz_full, W, n_threads = ncore)

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
        message(sprintf(
            "[geneSCOPE::computeL] Estimated matrix %.1f GB > limit %.1f GB; switching to chunked, streamed write",
            bytes_L, mem_limit_GB
        ))

        ## Use in-memory shared big.matrix to avoid unsupported 'shared' arg on filebacked in some versions
        ## Note: this stores the chunked result in RAM. Ensure sufficient memory is available.
        L_bm <- bigmemory::big.matrix(
            nrow = n_g, ncol = n_g, type = "double",
            init = NA_real_,
            dimnames = NULL,
            shared = TRUE
        )

        RhpcBLASctl::blas_set_num_threads(1)
        Sys.setenv(OMP_NUM_THREADS = ncore)

        ## ---- Process column blocks and write ----
        for (start in seq(1L, n_g, by = chunk_size)) {
            idx_chunk <- idx_keep[start:min(n_g, start + chunk_size - 1L)]
            # Use correct export function name
            L_block <- lee_L_cols(
                Xz_full, W,
                cols0 = as.integer(idx_chunk - 1L),
                n_threads = ncore
            )

            ## Write column block + symmetric columns
            L_bm[, idx_chunk] <- L_block
            L_bm[idx_chunk, ] <- t(L_block)

            rm(L_block)
            gc(verbose = FALSE)
        }

        ## ---- Add gene names (only two vectors) ----
        dimnames(L_bm) <- list(
            all_genes[idx_keep],
            all_genes[idx_keep]
        )

        Lmat <- L_bm
        Xuse <- Xz_full[, idx_keep, drop = FALSE]
    }

    dimnames(Lmat) <- list(
        row = colnames(Xuse),
        col = colnames(Xuse)
    )

    list(
        Lmat      = Lmat,
        X_used    = Xuse,
        X_full    = Xz_full,
        cells     = rownames(Xz_full),
        W         = W,
        grid_info = g_layer$grid_info,
        grid_name = grid_name, # ← Additional return, convenient for computeL
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

#' @noRd
.leeL_perm_block <- function(Xz, W, L_ref,
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
    blk_keys <- names(split_rows)
    n_blk <- length(split_rows)

    while (done < perms) {
        bsz <- min(block_size, perms - done)

        idx_mat <- replicate(bsz,
            {
                new_order <- sample(n_blk) # shuffle block order
                unlist(split_rows[new_order], use.names = FALSE)
            },
            simplify = "matrix"
        ) - 1L # 0-based for C++

        storage.mode(idx_mat) <- "integer"
        # Use correct export function name
        geCnt <- geCnt + lee_perm_block(Xz, W, idx_mat, as.integer(block_id) - 1L, L_ref, n_threads)
        done <- done + bsz
    }
    (geCnt + 1) / (perms + 1)
}
