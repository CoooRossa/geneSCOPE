#' @title Compute gene × gene correlation and save to CoordObj
#' Compute gene × gene correlation and (optionally) save to CoordObj
#'
#' @inheritParams geneCorrelation
#' @param cor_layer_name Character; name of the new sub-layer to save.
#'                       Default = paste0(method, "_cor")
#' @param write_back     Logical; if TRUE, write matrix to
#'                       coordObj@grid[[grid_name]][[cor_layer_name]]
#'
#' @return If write_back = TRUE, the modified CoordObj; otherwise the
#'         correlation matrix (dense or dgCMatrix)
#' @importFrom Matrix bandSparse sparseMatrix
#' @importFrom stats cor
#' @importFrom stats setNames
#' 
#' @export
geneCorrelation <- function(coordObj,
                                  grid_name       = "grid_lenGrid16",
                                  layer           = "Xz",
                                  genes           = NULL,
                                  method          = "pearson",
                                  blocksize       = 2000,
                                  as_sparse       = TRUE,
                                  thr             = 0.3,
                                  use_parallel    = FALSE,
                                  ncores          = 15,
                                  cor_layer_name  = NULL,
                                  write_back      = TRUE) {

  stopifnot(method %in% c("pearson", "spearman", "kendall"))

  g <- coordObj@grid[[grid_name]]

  ## ---------- fetch matrix ----------
  mat <- switch(layer,
    "Xz"     = g$Xz,
    "counts" = {
      requireNamespace("Matrix"); requireNamespace("data.table")
      dt <- data.table::as.data.table(g$counts)[count > 0]
      dt[, rn := match(grid_id, g$grid_info$grid_id)]
      i <- dt$rn
      gene_levels <- unique(dt$gene)
      j <- match(dt$gene, gene_levels)
      Matrix::sparseMatrix(
        i = i, j = j, x = dt$count,
        dims = c(nrow(g$grid_info), length(gene_levels)),
        dimnames = list(g$grid_info$grid_id, gene_levels)
      )
    },
    stop("Unsupported layer: ", layer)
  )

  ## ---------- gene filter ----------
  if (!is.null(genes)) {
    keep <- intersect(genes, colnames(mat))
    if (length(keep) == 0)
      stop("None of the requested genes found.")
    if (length(keep) < length(genes))
      warning("Dropped ", length(setdiff(genes, keep)), " genes not present.")
    mat <- mat[, keep, drop = FALSE]
  }

  nr <- nrow(mat); nc <- ncol(mat)
  message("Grids (obs): ", nr, " | Genes (vars): ", nc)

  ## ---------- helper ---------------
  corBlock <- function(j1, j2 = NULL) {
    m1 <- mat[, j1, drop = FALSE]
    m2 <- if (is.null(j2)) m1 else mat[, j2, drop = FALSE]
    stats::cor(m1, m2, method = method)
  }

  ## ---------- main loop ------------
  if (is.na(blocksize) || nc <= blocksize) {
    C <- corBlock(seq_len(nc))
  } else {
    idx <- split(seq_len(nc), ceiling(seq_len(nc) / blocksize))
    n   <- length(idx)

    C <- Matrix::bandSparse(nc, nc, k = 0,
                            diagonals = list(rep(1, nc)),
                            dimnames = list(colnames(mat), colnames(mat)))

    run_pair <- function(pair) {
      a <- pair[1]; b <- pair[2]
      cab <- corBlock(idx[[a]], idx[[b]])
      list(a = a, b = b, cab = cab)
    }

    pairs <- do.call(rbind,
                     lapply(seq_along(idx),
                            function(a) cbind(a, b = a:length(idx))))
    pairlist <- split(t(pairs), seq_len(nrow(pairs)))

    if (use_parallel && ncores > 1) {
      res <- parallel::mclapply(pairlist, run_pair,
                                mc.cores = ncores, mc.preschedule = FALSE)
    } else {
      res <- lapply(pairlist, run_pair)
    }

    for (r in res) {
      ja <- idx[[r$a]]; jb <- idx[[r$b]]
      if (r$a == r$b) {
        C[ja, jb] <- r$cab
      } else {
        C[ja, jb] <- r$cab
        C[jb, ja] <- t(r$cab)
      }
    }
  }

  ## ---------- sparsify -------------
  if (as_sparse) {
    C[abs(C) < thr] <- 0
    C <- as(C, "dgCMatrix")
  }

  ## ---------- write back -----------
  if (write_back) {
    if (is.null(cor_layer_name))
      cor_layer_name <- paste0(method, "_cor")
    coordObj@grid[[grid_name]][[cor_layer_name]] <- C
    return(coordObj)
  } else {
    return(C)
  }
}


#' @title Compute Δ = Lee's L − Pearson *r* with permutation p‑values **and QC**
#'
#' @description
#' For a given grid layer:
#' 1. Load Lee's L (`lee_stats_layer`) and Pearson correlation (`cor_layer_name`).
#' 2. Compute the difference Δ = **L − r**.
#' 3. Perform `perms` row‑randomisation permutations via an OpenMP C++ backend
#'    (`leeL_minusR_perm_ge_cpp()`), return two‑tailed empirical p‑values and
#'    adjust FDR.
#' 4. Classify gene pairs:
#'    * *L > 0 & Δ > 0 & FDR < 0.05* → "Spatial‑driven"
#'    * *L < 0 & Δ < 0*              → "False‑pair"
#' 5. **QC metrics** are included in the output (see *Return*).
#' 6. Optionally write the result back to
#'    `coordObj@grid[[grid_name]][["LeeMinusR_stats"]]`.
#'
#' @inheritParams computeDeltaLee
#' @param grid_name        Character; name of the grid layer to use.
#' @param lee_stats_layer  Character; name of the Lee's L sub‑layer to use.
#' @param cor_layer_name   Character; name of the Pearson correlation sub‑layer to use.
#' @param perms            Integer; number of permutations to perform.
#' @param fdr_method       Character; FDR adjustment method (default "BY").
#' @param seed             Integer; random seed for reproducibility.
#' @param write_back       Logical; if TRUE, write the results back to
#'                        `coordObj@grid[[grid_name]][["LeeMinusR_stats"]]`.
#' 
#' 
#' @importFrom Matrix as.matrix 
#' @importFrom stats p.adjust
#' @export
computeDeltaLee <- function(coordObj,
                                 grid_name        = "grid_lenGrid50",
                                 lee_stats_layer  = "LeeStats_Xz",
                                 cor_layer_name   = "pearson_Xz",
                                 perms            = 1000,
                                 fdr_method       = "BY",
                                 seed             = 1) {
  ## ---------- A. fetch matrices -----------------------------------------
  stopifnot(grid_name %in% names(coordObj@grid))
  L_mat <- coordObj@grid[[grid_name]][[lee_stats_layer]]$L
  R_mat <- coordObj@grid[[grid_name]][[cor_layer_name]]
  stopifnot(!is.null(L_mat), !is.null(R_mat))

  if (inherits(L_mat, "Matrix")) L_mat <- as.matrix(L_mat)
  if (inherits(R_mat, "Matrix")) R_mat <- as.matrix(R_mat)

  common <- intersect(colnames(L_mat), colnames(R_mat))
  L_mat  <- L_mat[common, common, drop = FALSE]
  R_mat  <- R_mat[common, common, drop = FALSE]
  Delta  <- L_mat - R_mat                                   # observed Δ

  ## ---------- B. permutation test --------------------------------------
  if (perms > 0L) {
    Xz <- coordObj@grid[[grid_name]]$Xz[, common, drop = FALSE]
    W  <- coordObj@grid[[grid_name]]$W
    ncell <- nrow(Xz)

    set.seed(seed)
    idx_mat <- replicate(perms, sample.int(ncell) - 1L, simplify = "matrix")

    # C++ backend returns exceedance counts (|Δ_perm| ≥ |Δ_obs|)
    geCnt <- FG2CLI:::leeL_minusR_perm_ge_cpp(Xz, W, idx_mat, Delta)

    P_mat <- (geCnt + 1) / (perms + 1)                    # empirical p‑values
  } else {
    warning("perms = 0 → no significance test, returning raw Δ only")
    P_mat <- matrix(1, nrow = nrow(Delta), ncol = ncol(Delta),
                    dimnames = dimnames(Delta))
    geCnt <- NULL
  }

  ## ---------- C. FDR ----------------------------------------------------
  fdr_vec <- p.adjust(as.vector(P_mat), method = fdr_method)
  FDR_mat <- matrix(fdr_vec, nrow = nrow(P_mat), dimnames = dimnames(P_mat))

  ## ---------- D. pair classification -----------------------------------
  class_mat <- matrix("Other", nrow = nrow(Delta), ncol = ncol(Delta),
                      dimnames = dimnames(Delta))
  class_mat[(L_mat > 0) & (Delta > 0) & (FDR_mat < 0.05)] <- "Spatial-driven"
  class_mat[(L_mat < 0) & (Delta < 0)]                    <- "False-pair"

  ## ---------- E. QC diagnostics ----------------------------------------
  qc <- list(
    summary_stats = list(
      Delta = summary(as.vector(Delta)),
      P     = summary(as.vector(P_mat)),
      FDR   = summary(as.vector(FDR_mat))
    ),
    class_counts = table(factor(class_mat,
                                levels = c("Spatial-driven",
                                           "False-pair",
                                           "Other"))),
    is_symmetric = isTRUE(all.equal(Delta, t(Delta))),
    diag_zero    = all(diag(Delta) == 0)
  )
  if (!is.null(geCnt)) qc$tail_cnt <- geCnt

  ## ---------- F. assemble output ---------------------------------------
  out <- list(
    Delta = Delta,
    P     = P_mat,
    FDR   = FDR_mat,
    class = class_mat,
    meta  = list(perms = perms, fdr_method = fdr_method),
    qc    = qc
  )
  
    coordObj@grid[[grid_name]][["LeeMinusR_stats"]] <- out
    return(coordObj)
  
}
