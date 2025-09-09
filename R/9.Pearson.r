#' @title Compute gene correlation matrix for grid or cell data
#' @description
#'   High-level wrapper that extracts expression data from the specified level
#'   and layer, applies column centering, and computes the Pearson correlation
#'   matrix using the C++ backend \code{pearson_block_cpp()}. The result is
#'   stored in the coordObj and the modified object is returned.
#' @param coordObj A \code{CoordObj} with populated \code{@grid} or \code{@cells}.
#' @param level Character. \code{"grid"} or \code{"cell"}.
#' @param grid_name Character (for grid level). Grid sub-layer name.
#' @param layer Character. Layer name to extract.
#' @param method Character. Correlation method (currently only \code{"pearson"}).
#' @param blocksize Integer. Block size for the C++ correlation computation.
#' @param ncores Integer. Number of OpenMP threads.
#' @param memory_limit_gb Numeric. Maximum memory limit for dense operations (default 32GB).
#' @param force_sparse Logical. Force sparse computation even for small matrices (default FALSE).
#' @param chunk_size Integer. Chunk size for large matrix processing (default 1000).
#' @param gene_subset Character vector or NULL. Specific genes to include in correlation analysis.
#' @param max_genes Integer. Maximum number of genes to include (triggers automatic selection if exceeded).
#' @param gene_selection_method Character. Method for automatic gene selection: "variance", "mean", "dropout", "random".
#' @param save_subset_info Logical. Whether to save information about gene selection.
#' @param force_compute Logical. Force computation even for very large matrices (default FALSE).
#' @param use_bigmemory Logical. Use file-backed matrices for large computations (default TRUE).
#' @param backing_path Character. Directory for temporary files (default tempdir()).
#' @param correation_slot_name Character. Name of the slot to store correlation matrix (default based on level and layer).
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#' @return The modified \code{CoordObj} (invisibly).
#' @importFrom Matrix t
#' @export
geneCorrelation <- function(coordObj,
                            level = c("cell", "grid"),
                            layer = "logCPM",
                            method = c("pearson", "spearman", "kendall"),
                            blocksize = 2000,
                            ncores = 16,
                            chunk_size = 1000,
                            memory_limit_gb = 16,
                            store_layer = "pearson_cor",
                            compute_fdr = TRUE,
                            fdr_method = "BH",
                            verbose = TRUE) {
  level <- match.arg(level)
  method <- match.arg(method)

  ## ---- Simple thread count control --------------------------------
  if (verbose) message("[geneSCOPE] Configuring thread settings...")
  # Get system information
  max_cores <- parallel::detectCores()
  if (is.na(max_cores) || max_cores <= 0) {
    max_cores <- 1
    if (verbose) message("[geneSCOPE] !!! Warning: Could not detect cores, using single core !!!")
  }

  # Simple safety limit
  ncores_safe <- min(ncores, max_cores - 1, 16) # Maximum 16 cores
  if (ncores_safe <= 0) {
    ncores_safe <- 1
  }

  if (ncores_safe < ncores && verbose) {
    message("[geneSCOPE] Reduced ncores from ", ncores, " to ", ncores_safe, " for safety")
  }

  ## ---- Set BLAS to single thread to avoid conflicts ----------------------
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
  }

  # Set environment variables to ensure BLAS single threading
  old_blas_env <- Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA)
  old_mkl_env <- Sys.getenv("MKL_NUM_THREADS", unset = NA)

  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1")

  # Restore environment variables on exit
  on.exit({
    if (is.na(old_blas_env)) {
      Sys.unsetenv("OPENBLAS_NUM_THREADS")
    } else {
      Sys.setenv(OPENBLAS_NUM_THREADS = old_blas_env)
    }
    if (is.na(old_mkl_env)) {
      Sys.unsetenv("MKL_NUM_THREADS")
    } else {
      Sys.setenv(MKL_NUM_THREADS = old_mkl_env)
    }
  })

  ## ---- Get expression matrix ------------------------------------
  if (verbose) message("[geneSCOPE] Loading expression matrix from level: ", level)
  if (level == "cell") {
    if (is.null(coordObj@cells) || is.null(coordObj@cells[[layer]])) {
      stop("Cell layer '", layer, "' not found in @cells")
    }
    expr_mat <- coordObj@cells[[layer]]
  } else {
    stop("Grid level correlation not yet implemented")
  }

  if (!inherits(expr_mat, "dgCMatrix")) {
    message("[geneSCOPE] Converting to sparse matrix...")
    expr_mat <- as(expr_mat, "dgCMatrix")
  }

  n_genes <- nrow(expr_mat)
  n_cells <- ncol(expr_mat)

  if (verbose) message("[geneSCOPE] Data dimensions: ", n_genes, " genes × ", n_cells, " cells")

  ## ---- Memory check and chunking strategy -------------------------------
  # Estimate correlation matrix size (n_genes × n_genes × 8 bytes)
  cor_matrix_gb <- (n_genes^2 * 8) / (1024^3)

  if (verbose) message("[geneSCOPE] Estimated correlation matrix size: ", round(cor_matrix_gb, 2), " GB")

  if (cor_matrix_gb > memory_limit_gb) {
    stop(
      "Correlation matrix too large (", round(cor_matrix_gb, 1), " GB) exceeds limit (",
      memory_limit_gb, " GB). Please reduce gene number or increase memory limit."
    )
  }

  ## ---- Data preprocessing --------------------------------------
  if (verbose) message("[geneSCOPE] Preprocessing expression matrix...")

  # Transpose matrix to cells × genes (pearson_cor expects this)
  expr_dense <- as.matrix(Matrix::t(expr_mat))

  # Center data (pearson_cor expects centered data)
  if (verbose) message("[geneSCOPE] Centering data...")
  expr_centered <- scale(expr_dense, center = TRUE, scale = FALSE)

  # Remove NaN/Inf
  expr_centered[!is.finite(expr_centered)] <- 0

  if (verbose) message("[geneSCOPE] Data preprocessing completed")

  ## ---- Call C++ correlation function --------------------------------
  if (verbose) message("[geneSCOPE] Starting correlation matrix computation (using ", ncores_safe, " threads)...")

  start_time <- Sys.time()

  # Use safe thread count to call C++ function
  cor_matrix <- tryCatch(
    {
      pearson_cor(
        X = expr_centered,
        bs = blocksize,
        n_threads = ncores_safe
      )
    },
    error = function(e) {
      # If error, try single thread
      message("[geneSCOPE] !!! Warning: Multithreaded computation failed, trying single thread: ", e$message, " !!!")
      pearson_cor(
        X = expr_centered,
        bs = blocksize,
        n_threads = 1
      )
    }
  )

  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))

  if (verbose) message("[geneSCOPE] Correlation matrix computation completed, time elapsed: ", round(elapsed, 2), " minutes")

  ## ---- Compute FDR (if needed) ----------------------------------
  fdr_matrix <- NULL
  if (compute_fdr) {
    if (verbose) message("[geneSCOPE] Computing FDR...")

    # Avoid large conversion from sparse to dense matrix
    # Directly extract upper triangle values from correlation matrix, avoid as.vector()
    upper_tri <- upper.tri(cor_matrix, diag = FALSE)

    # Use matrix indexing instead of converting entire matrix
    upper_indices <- which(upper_tri, arr.ind = TRUE)
    cor_values <- cor_matrix[upper_indices]

    # Remove NaN and Inf values
    valid_mask <- is.finite(cor_values) & !is.na(cor_values)
    cor_values_clean <- cor_values[valid_mask]
    valid_indices <- upper_indices[valid_mask, , drop = FALSE]

    if (length(cor_values_clean) == 0) {
      message("[geneSCOPE] !!! Warning: All correlation coefficients are NaN or Inf, cannot compute FDR !!!")
      fdr_matrix <- matrix(1, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
    } else {
      # Compute two-tailed p-values (only for valid values)
      n_samples <- n_cells
      df <- n_samples - 2

      if (df > 0) {
        # Compute t-statistic (avoid division by zero)
        cor_values_safe <- pmax(pmin(cor_values_clean, 0.9999), -0.9999)
        t_stat <- cor_values_safe * sqrt(df / (1 - cor_values_safe^2))

        # Compute two-tailed p-values
        p_values <- 2 * pt(-abs(t_stat), df = df)

        # Apply FDR correction
        fdr_values <- p.adjust(p_values, method = fdr_method)

        # Reconstruct symmetric matrix (initialize to 1)
        fdr_matrix <- matrix(1, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))

        # Fill only valid positions
        fdr_matrix[valid_indices] <- fdr_values
        fdr_matrix[valid_indices[, c(2, 1)]] <- fdr_values # Symmetric fill
        diag(fdr_matrix) <- 0 # Set diagonal to 0

        # Set row and column names
        rownames(fdr_matrix) <- rownames(cor_matrix)
        colnames(fdr_matrix) <- colnames(cor_matrix)

        if (verbose) message("[geneSCOPE] FDR computation completed, using method: ", fdr_method)
        if (verbose) message("[geneSCOPE] Valid correlation coefficients: ", length(cor_values_clean), " / ", length(cor_values))
      } else {
        message("[geneSCOPE] !!! Warning: Insufficient sample size, cannot compute FDR !!!")
        fdr_matrix <- matrix(1, nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
      }
    }
  }

  ## ---- Set gene names ------------------------------------
  gene_names <- rownames(expr_mat)
  rownames(cor_matrix) <- gene_names
  colnames(cor_matrix) <- gene_names

  ## ---- Store results ----------------------------------------
  coordObj@cells[[store_layer]] <- cor_matrix

  if (!is.null(fdr_matrix)) {
    fdr_layer_name <- paste0(store_layer, "_FDR")
    coordObj@cells[[fdr_layer_name]] <- fdr_matrix
    if (verbose) message("[geneSCOPE] FDR matrix stored in @cells$", fdr_layer_name)
  }

  if (verbose) message("[geneSCOPE] Correlation matrix stored in @cells$", store_layer)

  invisible(coordObj)
}

#' @title Use bigmemory for chunked correlation computation
#' @description For extremely large matrices, use file-backed storage to avoid memory overflow
#' @param mat Expression matrix
#' @param method Correlation method
#' @param chunk_size Chunk size
#' @param ncores Number of threads
#' @param force_compute Whether to force computation
#' @param memory_limit_gb Memory limit
#' @param backing_path File storage path
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#' @return File-backed big.matrix object or error message
.computeCorrelationBigMatrix <- function(mat, method, chunk_size, ncores,
                                         force_compute, memory_limit_gb, backing_path,
                                         verbose = TRUE) {
  n_genes <- if (inherits(mat, "dgCMatrix")) nrow(mat) else ncol(mat)
  result_size_gb <- (as.double(n_genes) * as.double(n_genes) * 8) / (1024^3)

  if (!force_compute && result_size_gb > 1000) {
    # Provide intelligent suggestions
    suggested_genes <- floor(sqrt(1000 * 1024^3 / 8)) # Approx. 310k genes for 1000GB

    stop(
      "Result matrix would be ", round(result_size_gb, 1),
      "GB which is too large. To force computation, set force_compute = TRUE.\n",
      "Alternatively, consider:\n",
      "1. Subsampling genes (e.g., select top ", suggested_genes, " variable genes)\n",
      "2. Using gene_subset parameter to specify genes of interest\n",
      "3. Computing correlations for gene subsets separately\n",
      "4. Using approximate correlation methods"
    )
  }

  if (force_compute && result_size_gb > 5000) { # 5TB limit
    stop(
      "Even with bigmemory, ", round(result_size_gb, 1),
      "GB is extremely large. Consider reducing the gene set."
    )
  }

  # Check bigmemory availability
  if (!requireNamespace("bigmemory", quietly = TRUE)) {
    stop("bigmemory package is required for large matrix computation but not available.")
  }

  # Cross-platform temporary directory
  os_type <- detectOS()
  if (backing_path == tempdir()) {
    backing_path <- switch(os_type,
      windows = Sys.getenv("TEMP", tempdir()),
      linux = Sys.getenv("TMPDIR", "/tmp"),
      macos = Sys.getenv("TMPDIR", "/tmp")
    )
  }

  if (!dir.exists(backing_path)) {
    dir.create(backing_path, recursive = TRUE)
  }

  # Create unique temporary file name to avoid conflicts
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  session_id <- sample.int(10000, 1)
  random_id <- sample(letters, 8, replace = TRUE)
  file_prefix <- paste0(
    "correlation_", timestamp, "_", session_id, "_",
    paste(random_id, collapse = ""), "_", Sys.getpid()
  )

  bm_file <- file.path(backing_path, paste0(file_prefix, ".bin"))
  bm_desc <- file.path(backing_path, paste0(file_prefix, ".desc"))

  # Ensure files do not exist, if they do add more randomness
  attempt <- 1
  while ((file.exists(bm_file) || file.exists(bm_desc)) && attempt <= 10) {
    extra_random <- sample.int(99999, 1)
    file_prefix_new <- paste0(file_prefix, "_", extra_random)
    bm_file <- file.path(backing_path, paste0(file_prefix_new, ".bin"))
    bm_desc <- file.path(backing_path, paste0(file_prefix_new, ".desc"))
    attempt <- attempt + 1
  }

  if (file.exists(bm_file) || file.exists(bm_desc)) {
    stop("Unable to create unique temporary files after 10 attempts")
  }

  # Clean up any existing files with the same name
  if (file.exists(bm_file)) {
    message("Removing existing backing file: ", basename(bm_file))
    unlink(bm_file)
  }
  if (file.exists(bm_desc)) {
    message("Removing existing descriptor file: ", basename(bm_desc))
    unlink(bm_desc)
  }

  # Set up automatic cleanup function
  cleanup_files <- function() {
    if (file.exists(bm_file)) {
      tryCatch(unlink(bm_file), error = function(e) NULL)
    }
    if (file.exists(bm_desc)) {
      tryCatch(unlink(bm_desc), error = function(e) NULL)
    }
  }

  # Automatically clean up at the end of the R session (as a backup)
  reg.finalizer(environment(), function(e) cleanup_files(), onexit = TRUE)

  if (verbose) message("[geneSCOPE] Creating file-backed correlation matrix (", round(result_size_gb, 1), "GB)...")
  if (verbose) message("[geneSCOPE] Files: ", basename(bm_file), " & ", basename(bm_desc))
  if (verbose) message("[geneSCOPE] This may take several hours for very large matrices...")

  # Use tryCatch to ensure files are cleaned up on error
  cor_bm <- tryCatch(
    {
      bigmemory::filebacked.big.matrix(
        nrow = n_genes, ncol = n_genes, type = "double",
        backingfile = basename(bm_file),
        descriptorfile = basename(bm_desc),
        backingpath = backing_path,
        init = 0,
        dimnames = list(
          if (inherits(mat, "dgCMatrix")) rownames(mat) else colnames(mat),
          if (inherits(mat, "dgCMatrix")) rownames(mat) else colnames(mat)
        )
      )
    },
    error = function(e) {
      cleanup_files()
      stop("Failed to create bigmemory matrix: ", e$message)
    }
  )

  # Set diagonal to 1
  if (verbose) message("[geneSCOPE] Initializing diagonal elements...")
  for (i in seq_len(n_genes)) {
    cor_bm[i, i] <- 1
  }

  if (verbose) message("[geneSCOPE] Computing correlations in chunks...")

  # Use smaller chunk size to avoid memory issues
  file_chunk_size <- min(chunk_size, 50)

  # Preprocessing: convert to dense matrix and center
  if (inherits(mat, "dgCMatrix")) {
    if (verbose) message("[geneSCOPE] Converting sparse to dense matrix in chunks...")
    mat_dense <- .sparseToDenseChunked(mat, file_chunk_size * 4)
  } else {
    mat_dense <- as.matrix(mat)
  }

  # Column centering
  if (verbose) message("[geneSCOPE] Column centering...")
  mat_centered <- scale(mat_dense, center = TRUE, scale = FALSE)
  mat_centered[is.na(mat_centered)] <- 0

  # Chunked correlation computation and direct write to file
  total_blocks <- ceiling(n_genes / file_chunk_size)
  current_block <- 0

  if (verbose) message("[geneSCOPE] Processing ", total_blocks^2, " correlation blocks...")

  # Add progress tracking
  progress_interval <- max(1, floor(total_blocks^2 / 20)) # 5% interval reporting

  for (i_start in seq(1, n_genes, by = file_chunk_size)) {
    i_end <- min(i_start + file_chunk_size - 1, n_genes)

    for (j_start in seq(i_start, n_genes, by = file_chunk_size)) {
      j_end <- min(j_start + file_chunk_size - 1, n_genes)

      current_block <- current_block + 1
      if (current_block %% progress_interval == 0) {
        percent_done <- round(100 * current_block / total_blocks^2, 1)
        if (verbose) message("[geneSCOPE]   Progress: ", percent_done, "% (", current_block, "/", total_blocks^2, " blocks)")
      }

      if (j_start < i_start) next # Only compute upper triangle

      # Extract blocks
      block_i <- mat_centered[, i_start:i_end, drop = FALSE]
      block_j <- if (j_start == i_start) block_i else mat_centered[, j_start:j_end, drop = FALSE]

      # Compute block correlation (using C++ acceleration)
      tryCatch(
        {
          if (j_start == i_start) {
            # Diagonal block - compute full correlation matrix
            block_cor <- pearson_block_cpp(block_i,
              bs = min(1000, ncol(block_i)),
              n_threads = ncores
            )
          } else {
            # Off-diagonal block - compute inter-block correlation
            combined_block <- cbind(block_i, block_j)
            full_cor <- pearson_block_cpp(combined_block,
              bs = min(1000, ncol(combined_block)),
              n_threads = ncores
            )

            # Extract relevant submatrix
            ni <- ncol(block_i)
            nj <- ncol(block_j)
            block_cor <- full_cor[1:ni, (ni + 1):(ni + nj)]
          }

          # Directly write to big.matrix
          cor_bm[i_start:i_end, j_start:j_end] <- block_cor
          if (j_start != i_start) {
            cor_bm[j_start:j_end, i_start:i_end] <- t(block_cor)
          }
        },
        error = function(e) {
          warning(
            "Error computing block [", i_start, ":", i_end, ", ",
            j_start, ":", j_end, "]: ", e$message
          )
          # Continue processing other blocks
        }
      )

      # Periodically force garbage collection
      if (current_block %% 50 == 0) {
        gc(verbose = FALSE)
      }
    }
  }

  if (verbose) message("[geneSCOPE] Correlation computation completed. Result stored in file-backed matrix.")
  if (verbose) message("[geneSCOPE] File locations:")
  if (verbose) message("[geneSCOPE]   Data: ", bm_file)
  if (verbose) message("[geneSCOPE]   Descriptor: ", bm_desc)
  if (verbose) message("[geneSCOPE] Note: These files will be automatically cleaned up when R session ends.")

  # Add cleanup information to return object
  attr(cor_bm, "cleanup_function") <- cleanup_files
  attr(cor_bm, "temp_files") <- c(bm_file, bm_desc)

  return(cor_bm)
}

#' @title Compute correlation between two sparse matrix blocks
#' @param block_i First block
#' @param block_j Second block
#' @param means_i Column means of first block
#' @param means_j Column means of second block
#' @param sds_i Column standard deviations of first block
#' @param sds_j Column standard deviations of second block
#' @return Correlation matrix block
.computeBlockCorrelation <- function(block_i, block_j, means_i, means_j, sds_i, sds_j) {
  n_obs <- nrow(block_i)
  n_i <- ncol(block_i)
  n_j <- ncol(block_j)

  cor_block <- matrix(0, n_i, n_j)

  for (i in seq_len(n_i)) {
    for (j in seq_len(n_j)) {
      if (sds_i[i] == 0 || sds_j[j] == 0) {
        cor_block[i, j] <- 0
        next
      }

      # Safely extract column vectors and convert to numeric vectors
      x_col <- block_i[, i, drop = FALSE]
      y_col <- block_j[, j, drop = FALSE]

      # Convert to numeric vectors
      x_vec <- as.numeric(x_col)
      y_vec <- as.numeric(y_col)

      # Compute covariance - use standard formula
      if (all(x_vec == 0) && all(y_vec == 0)) {
        # Both columns are zero
        cov_xy <- 0
      } else {
        # Use standard covariance formula
        mean_x <- means_i[i]
        mean_y <- means_j[j]
        cov_xy <- sum((x_vec - mean_x) * (y_vec - mean_y)) / (n_obs - 1)
      }

      # Compute correlation coefficient
      cor_block[i, j] <- cov_xy / (sds_i[i] * sds_j[j])
    }
  }

  return(cor_block)
}

#' @title Chunked conversion of sparse matrix to dense matrix
#' @param sparse_mat Sparse matrix
#' @param chunk_size Chunk size
#' @return Dense matrix
.sparseToDenseChunked <- function(sparse_mat, chunk_size) {
  n_cols <- ncol(sparse_mat)
  dense_mat <- matrix(0, nrow(sparse_mat), n_cols)

  for (start in seq(1, n_cols, by = chunk_size)) {
    end <- min(start + chunk_size - 1, n_cols)
    dense_mat[, start:end] <- as.matrix(sparse_mat[, start:end])
  }

  return(dense_mat)
}

#' @title Intelligent gene subset selection
#' @description Select gene subset intelligently based on dataset size and available memory
#' @param mat Expression matrix
#' @param max_genes Maximum number of genes
#' @param method Selection method: "variance", "mean", "dropout", "random"
#' @return Selected gene names vector
.selectTopGenes <- function(mat, max_genes = 50000, method = "variance") {
  if (ncol(mat) <= max_genes) {
    return(colnames(mat))
  }

  message("[geneSCOPE] Selecting top ", max_genes, " genes by ", method, "...")

  if (method == "variance") {
    # Compute gene variance
    gene_vars <- Matrix::colMeans(mat^2) - Matrix::colMeans(mat)^2
    top_idx <- order(gene_vars, decreasing = TRUE)[1:max_genes]
  } else if (method == "mean") {
    # Select by mean expression
    gene_means <- Matrix::colMeans(mat)
    top_idx <- order(gene_means, decreasing = TRUE)[1:max_genes]
  } else if (method == "dropout") {
    # Select by dropout rate (proportion of cells expressing)
    gene_dropout <- Matrix::colMeans(mat > 0)
    top_idx <- order(gene_dropout, decreasing = TRUE)[1:max_genes]
  } else if (method == "random") {
    # Random selection
    top_idx <- sample(ncol(mat), max_genes)
  } else {
    stop("Unknown selection method: ", method)
  }

  return(colnames(mat)[top_idx])
}
