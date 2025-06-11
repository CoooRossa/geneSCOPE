#' @title Compute Lee’s L from Pre-constructed Grid Data (No Global Shifting)
#' @description
#'   This function computes Lee’s L statistics directly using either a column-centered
#'   or Z-score–normalized matrix without applying any “+1” or row-wise shifting. This
#'   preserves symmetry around zero in the input data so that Lee’s L can take both
#'   positive and negative values. Given a spatial weights matrix W, for each gene the
#'   spatial lag vector \(\tilde{x}_i = \sum_j W_{ij} x_j\) is column-centered and then
#'   Pearson correlation is computed and multiplied by the corresponding spatial autocorrelation
#'   coefficient, following the definition in Lee (2001).
#'
#' @param coordObj      A coordObj processed through createCoordObj + countMPG + norMPG
#'                      (or equivalent layers) + computeSpatialWeights.
#' @param grid_name     Character or NULL. The name of the grid sublayer to use (e.g., "grid_lenGrid50").
#'                      If NULL and there is only one sublayer under coordObj\$grid, that layer is used;
#'                      otherwise this must be provided explicitly.
#' @param norm_layer    Character or NULL. The name of the normalization layer to use:
#'                      - If "Xz", reads coordObj\$grid[[grid_name]][["Xz"]], which is already Z-score normalized
#'                        (column means ≈ 0).
#'                      - If NULL, uses the raw counts layer and applies column-centered normalization.
#' @param genes         Character vector or NULL. If non-NULL, compute Lee’s L only for these genes:
#'                      - If within = TRUE, returns a |genes|×|genes| submatrix.
#'                      - If within = FALSE, returns a |genes|×|all genes| matrix.
#' @param within        Logical (default TRUE):
#'                      - TRUE: compute Lee’s L within the provided genes list (|genes|×|genes|).
#'                      - FALSE: compute Lee’s L between the provided genes list and all genes (|genes|×|all genes|).
#' @param ncore         Integer. Number of parallel cores to use when use_cpp = FALSE.
#' @param use_cpp       Logical (default TRUE). Whether to call the C++ interface leeL_cpp for acceleration.
#'
#' @return A list containing:
#'   \item{Lmat}{Lee’s L submatrix (if within = TRUE, a |genes|×|genes| matrix; otherwise |genes|×|all genes|).}
#'   \item{X_used}{The subsetted matrix used for Lee’s L (rows = grid IDs; columns = selected genes).}
#'   \item{X_full}{The full column-centered or Z-score matrix (rows = grid IDs; columns = all genes).}
#'   \item{grid_info}{coordObj\$grid[[grid_name]]\$grid_info (contains each cell’s spatial boundary info).}
#'   \item{cells}{All grid IDs (in the same order as rows of X_full).}
#'   \item{W}{The Queen adjacency weights matrix coordObj\$grid[[grid_name]]\$W (|cells|×|cells|, entries 0/1).}
#'
#' @importFrom Matrix sparseMatrix t
#' @importFrom spdep cell2nb nb2listw listw2mat
#' @importFrom RhpcBLASctl blas_set_num_threads
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom foreach foreach %dopar%
#' @export
.computeLeeL <- function(coordObj,
                                  grid_name  = NULL,
                                  norm_layer = "Xz",
                                  genes      = NULL,
                                  within     = TRUE,
                                  ncore      = 1,
                                  use_cpp    = TRUE) {
  ## ========== 0. Pre-checks ==========
  stopifnot(!is.null(coordObj@grid))
  # If grid_name is NULL but only one sublayer exists, auto-select it
  if (is.null(grid_name)) {
    sub_layers <- names(coordObj@grid)
    if (length(sub_layers) == 1L) {
      grid_layer_name <- sub_layers[[1]]
    } else {
      stop("Multiple sublayers exist under coordObj@grid; specify grid_name explicitly.")
    }
  } else {
    grid_layer_name <- as.character(grid_name)
  }
  stopifnot(grid_layer_name %in% names(coordObj@grid))
  
  # Extract the chosen grid sublayer
  grid_layer     <- coordObj@grid[[grid_layer_name]]
  stopifnot("grid_info" %in% names(grid_layer), "W" %in% names(grid_layer))
  grid_info_full <- grid_layer$grid_info    # Raw grid information (includes boundaries, gx, gy, etc.)
  W_full         <- grid_layer$W            # Full Queen adjacency sparse matrix W
  
  ## ========== 1. Construct X_full: Column-centered or Z-score Matrix ==========
  if (!is.null(norm_layer)) {
    # ---- 1.1 Read user-provided Z-score matrix (assumed column means ≈ 0, variance ≈ 1) ----
    stopifnot(norm_layer %in% names(grid_layer))
    X_in <- grid_layer[[norm_layer]]
    stopifnot(is.matrix(X_in) || inherits(X_in, "Matrix"))
    # Convert sparse or dense to a regular matrix
    X_dense_in    <- as.matrix(X_in)
    # Replace any NA with 0 to ensure completeness
    X_dense_in[is.na(X_dense_in)] <- 0
    # Extract all genes and all grid IDs (ensuring order matches grid_info)
    genes_all  <- colnames(X_dense_in)
    cells_all  <- grid_info_full$grid_id
    # Build a zero matrix of dimension (cells_all × genes_all) and copy values from X_dense_in
    X_full <- matrix(
      0,
      nrow = length(cells_all),
      ncol = length(genes_all),
      dimnames = list(cells_all, genes_all)
    )
    common_cells <- intersect(rownames(X_dense_in), cells_all)
    common_genes <- intersect(colnames(X_dense_in), genes_all)
    if (length(common_cells) > 0L && length(common_genes) > 0L) {
      X_full[common_cells, common_genes] <- X_dense_in[common_cells, common_genes]
    }
    # At this point, each column’s mean ≈ 0 (since it’s copied from a Z-score matrix);
    # values may be positive or negative, satisfying Lee’s L requirement for Pearson correlation.
  } else {
    # ---- 1.2 Use raw counts layer and apply column-centered normalization ----
    stopifnot("counts" %in% names(grid_layer))
    counts_dt <- grid_layer$counts
    stopifnot(all(c("grid_id", "gene", "count") %in% colnames(counts_dt)))
    # Extract all genes and all grid IDs
    genes_all <- sort(unique(counts_dt$gene))
    cells_all <- grid_info_full$grid_id
    # Build a sparse matrix with rows = genes, columns = grid IDs
    P_mat <- sparseMatrix(
      i        = match(counts_dt$gene,    genes_all),
      j        = match(counts_dt$grid_id, cells_all),
      x        = counts_dt$count,
      dims     = c(length(genes_all), length(cells_all)),
      dimnames = list(genes_all, cells_all)
    )
    # Transpose to get the “grid × gene” raw counts matrix
    X_counts <- t(P_mat)  # Dimensions: |cells_all| × |genes_all|
    X_counts <- as.matrix(X_counts)
    # Column-center each gene: X_centered[, f] = X_counts[, f] - mean(X_counts[, f])
    col_means <- colMeans(X_counts, na.rm = TRUE)
    # X_full    <- sweep(X_counts, 2, col_means, FUN = "-")
    # For columns that are all zero or NaN, replace with zero
    X_full[is.na(X_full)] <- 0
    # Now each column’s mean = 0, with positive/negative structure intact.
  }
  
  ## ========== 2. Subset Columns According to 'genes' Parameter ==========
  all_genes <- colnames(X_full)
  cells     <- rownames(X_full)
  if (!is.null(genes)) {
    keep <- intersect(genes, all_genes)
    if (length(keep) == 0L) {
      stop("Provided genes do not overlap with the complete gene set; cannot compute.")
    }
    if (within) {
      # Keep only the specified genes
      X_used <- X_full[, keep, drop = FALSE]
    } else {
      # Keep all genes
      X_used <- X_full
    }
  } else {
    # If genes = NULL, treat as within = TRUE but keep all genes
    keep   <- all_genes
    X_used <- X_full
    within <- TRUE
  }
  
  ## ========== 3. Subset Spatial Weights W and grid_info ==========
  # Ensure W_full rows/columns match the order of 'cells'
  W         <- W_full[cells, cells, drop = FALSE]
  grid_info <- grid_info_full[match(cells, grid_info_full$grid_id), ]
  
  ## ========== 4. Compute Full Lee’s L Matrix ==========
  ncell    <- nrow(X_full)
  all_f    <- length(all_genes)
  # Compute prefactor: prefac = n / sum( (∑_j w_{ij})^2 )
  # Here W_full is binary (0/1) and not row-normalized; compute row sums and their squared sum
  wij_sum_vec <- rowSums(W)                   # For each i, ∑_j w_{ij}
  S0           <- sum(wij_sum_vec^2)          # ∑_i (∑_j w_{ij})^2
  prefac       <- ncell / S0                  # Prefactor as in Lee (2001)
  
  # Preallocate L_full (dimensions: all_genes × all_genes)
  L_full <- matrix(NA_real_, all_f, all_f, dimnames = list(all_genes, all_genes))
  
  if (use_cpp) {
    # ---- 4.1 Call C++ interface for parallel computation ----
    L_full <- FG2CLI:::leeL_cpp(as.matrix(X_full), W)
    dimnames(L_full) <- list(all_genes, all_genes)
  } else {
    # ---- 4.2 R-side iterative computation per gene ----
    # Denominator vector: ∑_i X_full[i, f]^2 for each gene f
    den_vec <- colSums(X_full^2)
    # If multithreading via BLAS is active, set to 1 to avoid conflicts
    if (ncore > 1L) {
      ## Enable OpenMP in C++ while forcing BLAS to single‑thread
      Sys.setenv(OMP_NUM_THREADS = ncore)
      RhpcBLASctl::blas_set_num_threads(1)
    } else {
      Sys.unsetenv("OMP_NUM_THREADS")
    }
    if (ncore > 1L) {
      cl <- parallel::makeCluster(ncore)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, c("X_full", "W", "den_vec", "prefac"), envir = environment())
      parallel::clusterEvalQ(cl, {
        Sys.setenv(OMP_NUM_THREADS = ncore)
        library(Matrix)
        RhpcBLASctl::blas_set_num_threads(1)
      })
      # Parallel loop: each iteration computes one column f of L_full
      L_full <- foreach::foreach(
        f = seq_len(all_f), 
        .combine = "cbind", 
        .packages = "Matrix"
      ) %dopar% {
        zf  <- X_full[, f]         # Values of gene f across all grid cells
        Wzf <- W %*% zf            # Spatial lag \tilde{z}_f = W z_f
        # Numerator: ∑_i z_f[i] * Wz_f[i] (since column-centered, mean = 0)
        num <- as.numeric(crossprod(zf, Wzf))
        # Denominator: sqrt( ∑_i z_f[i]^2 * den_vec[f] )
        den        <- sqrt(sum(zf^2) * den_vec[f])
        den[den == 0] <- NA        # If denominator is 0 (constant column), L = NA
        prefac * (num / den)
      }
      colnames(L_full) <- all_genes
    } else {
      # Single-threaded mode
      for (f in seq_len(all_f)) {
        zf  <- X_full[, f]
        Wzf <- W %*% zf
        num <- as.numeric(crossprod(zf, Wzf))
        den        <- sqrt(sum(zf^2) * den_vec[f])
        den[den == 0] <- NA
        L_full[, f] <- prefac * (num / den)
      }
      colnames(L_full) <- all_genes
    }
  }
  
  ## ========== 5. Extract Submatrix Based on 'within' ==========
  if (within) {
    Lmat <- L_full[keep, keep, drop = FALSE]
  } else {
    Lmat <- L_full[keep, , drop = FALSE]
  }
  
  ## ========== 6. Return Results ==========
  list(
    Lmat      = Lmat,       # Lee’s L submatrix
    X_used    = X_used,     # Subsetted input matrix
    X_full    = X_full,     # Full column-centered or Z-score matrix
    grid_info = grid_info,  # Spatial grid information
    cells     = cells,      # All grid IDs
    W         = W           # Subsetted spatial weights matrix
  )
}
#' @title Compute and Save Lee’s L Statistics with QC from Grid Data
#' @description
#'   This function computes Lee’s L statistics for a given grid layer in a coordObj,
#'   using either a column-centered or Z-score–normalized matrix.
#'   It saves the results as a new layer in coordObj@grid[[grid_name]].
#' @param coordObj      A coordObj processed through createCoordObj + norMPG
#'                      (or equivalent layers) + computeSpatialWeights.
#' @param grid_name     Character or NULL. The name of the grid sublayer to use (e.g., "grid_lenGrid50").
#'                      If NULL and there is only one sublayer under coordObj\$grid, that layer is used;
#'                      otherwise this must be provided explicitly.
#' @param genes         Character vector or NULL. If non-NULL, compute Lee’s L only for these genes:
#'                       - If within = TRUE, returns a |genes|×|genes| submatrix.
#'                       - If within = FALSE, returns a |genes|×|all genes| matrix.
#' @param within        Logical (default TRUE):
#'                       - TRUE: compute Lee’s L within the provided genes list (|genes|×|genes|).
#'                       - FALSE: compute Lee’s L between the provided genes list and all genes (|genes|×|all genes|).
#' @param ncore         Integer. Number of parallel cores to use when use_cpp = FALSE.
#' @param perms         Integer. Number of Monte-Carlo permutations to perform for significance testing.
#' @param L_min         Numeric. Minimum Lee’s L value to consider for QC metrics (default 0.15).
#' @param use_cpp       Logical (default TRUE). Whether to call the C++ interface for acceleration.
#' @param norm_layer    Character or NULL. The name of the normalization layer to use:
#'                       - If "Xz", reads coordObj\$grid[[grid_name]][["Xz"]], which is already Z-score normalized
#'                         (column means ≈ 0).
#'                       - If NULL, uses the raw counts layer and applies column-centered normalization.
#' @param lee_stats_layer_name Character or NULL. If provided, the results will be saved under this name;
#                              otherwise defaults to "LeeStats_Xz" or "LeeStats_Xcounts" based on norm_layer.
#' @return The modified `coordObj`, with a new entry:
#'          - `coordObj@grid[[grid_name]][[lee_stats_layer_name]]`: a list containing:
#'          - `L`: Lee’s L submatrix (if within = TRUE, a |genes|×|genes| matrix; otherwise |genes|×|all genes|).
#'          - `Z`: Analytical Z-score matrix (same dimensions as L).
#'          - `P`: Permutation p-values matrix (same dimensions as L, NA if perms = 0).
#'          - `FDR`: BH-adjusted p-values matrix (same dimensions as L).
#'          - `betas`: Spatial gradient coefficients (|genes|×2 matrix, columns = beta_x, beta_y).
#'          - `qc`: QC metrics list (if within = TRUE or genes is NULL, contains edge_density, components, modularity_Q,
#'             mean_degree, sd_degree, hub_ratio, sig_edge_frac; otherwise NULL).
#' @importFrom Matrix Matrix
#' @importFrom spdep cell2nb nb2listw listw2mat
#' @importFrom RhpcBLASctl blas_set_num_threads
#' @importFrom foreach foreach %dopar%
#' @importFrom igraph graph_from_adjacency_matrix components modularity cluster_louvain
#' @importFrom stats p.adjust
#' @importFrom data.table as.data.table setDT 
#' @export
addLeeStats <- function(coordObj,
                                 grid_name            = NULL,
                                 genes                = NULL,
                                 within               = TRUE,
                                 ncore                = 1,
                                 perms                = 0,
                                 L_min                = 0.15,
                                 use_cpp              = TRUE,
                                 norm_layer           = "Xz",
                                 lee_stats_layer_name = NULL) {  # ---- 新增参数 ----
  # ---------- 1. Pre-checks ----------
  stopifnot(!is.null(coordObj@grid))
  if (is.null(grid_name)) {
    sub_layers <- names(coordObj@grid)
    if (length(sub_layers) == 1L) {
      grid_layer_name <- sub_layers[[1]]
    } else {
      stop("Multiple sublayers exist under coordObj@grid; please specify grid_name explicitly.")
    }
  } else {
    grid_layer_name <- as.character(grid_name)
  }
  stopifnot(grid_layer_name %in% names(coordObj@grid))

  # ---------- Set OpenMP/BLAS threading for Lee's L calculation ----------
  if (ncore > 1L) {
    Sys.setenv(OMP_NUM_THREADS = ncore)
    RhpcBLASctl::blas_set_num_threads(1)
  } else {
    Sys.unsetenv("OMP_NUM_THREADS")
  }

  # ---------- 2. Extract X_full and W ----------
  res <- .computeLeeL(coordObj,
                               grid_name  = grid_layer_name,
                               norm_layer = norm_layer,
                               genes      = genes,
                               within     = within,
                               ncore      = ncore,
                               use_cpp    = use_cpp)
  L       <- res$Lmat
  X_used  <- res$X_used
  X_full  <- res$X_full
  W       <- res$W
  n       <- nrow(X_full)

  grid_info <- coordObj@grid[[grid_layer_name]]$grid_info

  # ---------- 3. Analytical Z ----------
  if (within || is.null(genes)) {
    nbins <- max(grid_info$gx)
    nb    <- spdep::cell2nb(nbins, nbins, type = "queen")
    S0    <- sum(spdep::listw2mat(spdep::nb2listw(nb, style = "B")))
    E_L   <- -1 / (n - 1)
    VarL  <- (n^2 * (n - 2)) / ((n - 1)^2 * (n - 3) * S0)
    Z_mat <- (L - E_L) / sqrt(VarL)
  } else {
    Z_mat <- t(apply(L, 1, function(v) {
      mu  <- mean(v, na.rm = TRUE)
      sdv <- stats::sd(v, na.rm = TRUE)
      if (sdv == 0) rep(0, length(v)) else (v - mu) / sdv
    }))
    rownames(Z_mat) <- rownames(L)
    colnames(Z_mat) <- colnames(L)
  }

  # ---------- 4. Monte-Carlo permutations 〈STREAMING VERSION〉 ----------
  P <- NULL
  if (perms > 0L) {
    idx_mat <- replicate(perms, sample.int(nrow(X_full)), simplify = "matrix")
    geCnt   <- FG2CLI:::leeL_perm_ge_cpp(X_full, W, idx_mat - 1L, L)  # Armadillo 用 0-based
    P       <- (geCnt + 1) / (perms + 1)
  }
  # ---------- 5. FDR adjustment ----------
  # If P is NULL, compute p-values from Z-scores
  # If P is provided, it is assumed to be the raw p-values.
  if (!is.null(P)) {
    p_raw <- P
  } else {
    p_raw <- 2 * stats::pnorm(-abs(Z_mat))
  }
  p_vec        <- as.vector(p_raw)
  p_vec[is.na(p_vec)] <- 1
  q_vec        <- stats::p.adjust(p_vec, method = "BH")
  FDR_mat      <- matrix(q_vec, nrow = nrow(p_raw), ncol = ncol(p_raw),
                         dimnames = dimnames(p_raw))

  # ---------- 6. Compute spatial gradient coefficients ----------
  centres <- data.frame(
    x = (grid_info$xmin + grid_info$xmax) / 2,
    y = (grid_info$ymin + grid_info$ymax) / 2
  )[match(res$cells, grid_info$grid_id), ]
  if (use_cpp) {
    betas_full <- FG2CLI:::grad_betas_cpp(as.matrix(X_full), centres$x, centres$y)
  } else {
    betas_full <- t(apply(X_full, 2, function(v) {
      coefs <- stats::lm(v ~ centres$x + centres$y)$coefficients
      c(beta_x = coefs[2], beta_y = coefs[3])
    }))
  }
  betas <- if (!is.null(genes)) betas_full[intersect(genes, rownames(betas_full)), , drop = FALSE] else betas_full

  # ---------- 7. Quality control metrics ----------
  qc <- NULL
  if (within || is.null(genes)) {
    A_bin <- abs(L) >= L_min
    g_tmp <- igraph::graph_from_adjacency_matrix(A_bin, mode = "undirected", diag = FALSE)
    qc <- list(
      edge_density  = 2 * sum(A_bin[upper.tri(A_bin)]) /
                      (ncol(A_bin) * (ncol(A_bin) - 1)),
      components    = igraph::components(g_tmp)$no,
      modularity_Q  = igraph::modularity(g_tmp,
                       membership = igraph::cluster_louvain(g_tmp)$membership),
      mean_degree   = mean(igraph::degree(g_tmp)),
      sd_degree     = stats::sd(igraph::degree(g_tmp)),
      hub_ratio     = mean(igraph::degree(g_tmp) >
                           2 * median(igraph::degree(g_tmp))),
      sig_edge_frac = if (is.null(P)) NA
                      else mean(P[upper.tri(P)] < 0.05, na.rm = TRUE)
    )
  }

  # ---------- 8. Choose layer name ----------
  chosen_layer_name <- if (!is.null(lee_stats_layer_name)) {
    lee_stats_layer_name
  } else if (!is.null(norm_layer)) {
    paste0("LeeStats_", norm_layer)
  } else {
    "LeeStats_Xcounts"
  }

  dimnames(FDR_mat) <- dimnames(L)

  # ---------- 9. Save results to coordObj ----------
  out_list <- list(
    L      = L,
    Z      = Z_mat,
    P      = P,
    grad   = betas,
    L_min  = L_min,
    qc     = qc,
    FDR    = FDR_mat 
  )
  coordObj@grid[[grid_layer_name]][[chosen_layer_name]] <- out_list
  return(coordObj)
}
