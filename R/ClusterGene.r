#' @title Two-step gene clustering: Lee’s L macro-modules → Δ-L refined sub-modules
#'
#' @description
#' **Step 1 Macro modules** — keep only positive Lee’s L edges (`L ≥ L_min`) and apply Louvain/Leiden to obtain macro-module IDs *macro*.
#'
#' **Step 2 Δ-L refinement** (within each macro module)
#'   * **Conflict edges** = `Δ ≥ delta_split_min & FDR_Δ ≤ delta_FDR_max` **or** `L < 0`.
#'   * **Candidate genes** = genes incident to at least one conflict edge.
#'   * Non-candidate genes are labeled `macro.0`.
#'   * **Glue edges** = `L > 0 & Δ < 0`; connected components of the glue graph define sub-modules `macro.1`, `macro.2`, …
#'     If no glue edges exist, each candidate gene forms its own sub-module.
#'
#' **Output** (stored in `coordObj@meta.data`)
#'   * `{cluster_name}`       — final module ID (macro module unless a sub-module contains >2 genes, in which case `macro.sub` overwrites)
#'   * `{cluster_name}_sub`   — full Δ-L sub-module ID (`macro.sub`)
#'
#' @inheritParams graphClusterGenes
#' @param lee_stats_layer Sub-layer that stores Lee’s L / FDR (default "LeeStats_Xz")
#' @param delta_layer     Sub-layer that stores Δ / FDR (default "LeeMinusR_stats")
#' @param L_min           Minimum positive Lee’s L to retain an edge
#' @param algo            Community detection algorithm: "louvain" or "leiden"
#' @param resolution      Resolution parameter for Louvain/Leiden
#' @param delta_split_min Minimum Δ to trigger a conflict edge
#' @param delta_FDR_max   FDR threshold for Δ
#' @param split_min_size  Attempt sub-splitting only if a macro module contains at least this many genes
#' @param cluster_name    Prefix of output column; generated automatically if NULL
#' @param use_delta_L     Whether to perform the Δ-L second step (default TRUE)
#'
#' @importFrom igraph cluster_louvain cluster_leiden graph_from_adjacency_matrix
#' @importFrom igraph components membership
#' @importFrom igraph E V
#' @importFrom stats setNames
#' 
#' @return The modified `coordObj`
#' @export

graphClusterGenes <- function(
  coordObj,
  grid_name       = NULL,
  lee_stats_layer = "LeeStats_Xz",
  delta_layer     = "LeeMinusR_stats",
  L_min           = 0,
  algo            = "louvain",
  resolution      = 1,
  pct_min         = "q0",
  pct_max         = "q100",
  use_delta_L     = FALSE,
  delta_split_min = 0,
  delta_FDR_max   = 0.05,
  split_min_size  = 20,
  cluster_name    = NULL,
  drop_isolated   = TRUE
) {
  ## ---------- 0. Select grid layer ----------
  algo <- match.arg(algo)
  if (is.null(grid_name)) {
    subs <- names(coordObj@grid)
    if (length(subs) == 1L) grid_name <- subs else
      stop("Multiple grid sub-layers detected; please specify grid_name.")
  }
  stopifnot(grid_name %in% names(coordObj@grid))
  lay <- coordObj@grid[[grid_name]]

  ## ---------- 1. Read Lee’s L ----------
  if (!lee_stats_layer %in% names(lay))
    stop("Layer '", lee_stats_layer, "' not found under grid '", grid_name, "'.")
  Lmat <- as.matrix(lay[[lee_stats_layer]]$L)
  if (is.null(dimnames(Lmat))) stop("Lee’s L matrix must have dimnames.")
  genes <- rownames(Lmat)

  ## ---------- 2. Macro-module clustering ----------
  A <- Lmat; A[A < L_min] <- 0; diag(A) <- 0; A[A < 0] <- 0
  W <- 1 - A
  if (exists(".filter_matrix_by_quantile", mode = "function"))
    W <- .filter_matrix_by_quantile(W, pct_min, pct_max)
  W[W == 1] <- 0

  keep <- if (drop_isolated) rowSums(W != 0) > 0 else rep(TRUE, length(genes))
  kept_genes <- genes[keep]
  if (length(kept_genes) < 2L)
    stop("Fewer than two genes remain after filtering; relax thresholds or set drop_isolated = FALSE.")

  g_macro <- igraph::graph_from_adjacency_matrix(W[kept_genes, kept_genes], "undirected", weighted = TRUE)
  comm_macro <- if (algo == "leiden")
    igraph::cluster_leiden(g_macro, weights = igraph::E(g_macro)$weight, resolution_parameter = resolution)
  else igraph::cluster_louvain(g_macro, weights = igraph::E(g_macro)$weight)
  macro <- setNames(rep(NA_integer_, length(genes)), genes)
  macro[names(igraph::membership(comm_macro))] <- igraph::membership(comm_macro)

  ## ---------- 3. Δ-L sub-modules ----------
  submod <- setNames(paste0(macro, ".0"), genes)  # default .0

  if (use_delta_L) {
    if (!delta_layer %in% names(lay))
      stop("Layer '", delta_layer, "' not found under grid '", grid_name, "'.")
    Dmat <- as.matrix(lay[[delta_layer]]$Delta)
    Fmat <- as.matrix(lay[[delta_layer]]$FDR)
    dimnames(Dmat) <- dimnames(Fmat) <- list(genes, genes)

    for (m in na.omit(unique(macro[kept_genes]))) {
      ## genes in this macro
      g_mod <- names(macro)[macro == m & !is.na(macro)]
      g_mod <- intersect(g_mod, rownames(Dmat))
      if (length(g_mod) < split_min_size) next

      ## conflict edges
      conflict <- (Dmat[g_mod, g_mod] >= delta_split_min & Fmat[g_mod, g_mod] <= delta_FDR_max) |
                  (Lmat[g_mod, g_mod] < 0)
      g_conf <- unique(c(rownames(conflict)[rowSums(conflict) > 0],
                         colnames(conflict)[colSums(conflict) > 0]))
      if (length(g_conf) == 0) next

      ## retain genes that still have positive L to others
      posL <- Lmat[g_conf, g_conf] > 0; diag(posL) <- FALSE
      g_pos <- rownames(posL)[rowSums(posL) > 0]
      if (length(g_pos) == 0) {
        submod[g_conf] <- paste0(m, ".", seq_along(g_conf))
        next
      }

      ## glue graph: positive L & Δ < delta_split_min
      glue_mask <- (posL[g_pos, g_pos]) & (Dmat[g_pos, g_pos] < delta_split_min)
      diag(glue_mask) <- FALSE

      if (!any(glue_mask)) {
        submod[g_pos] <- paste0(m, ".", seq_along(g_pos))
      } else {
        g_glue <- igraph::graph_from_adjacency_matrix(glue_mask, "undirected")
        cc <- igraph::components(g_glue)$membership
        submod[g_pos] <- paste0(m, ".", cc)
      }
    }
  }

  ## ---------- 4. Compose final cluster column ----------
  sub_sizes <- table(submod)
  large_sub <- names(sub_sizes)[sub_sizes > 2]
  cluster_final <- as.character(macro); names(cluster_final) <- genes
  cluster_final[submod %in% large_sub] <- submod[submod %in% large_sub]
  cluster_final[grepl("NA", cluster_final)] <- NA  # keep NAs explicit

  ## ---------- 5. Write to meta.data ----------
  if (is.null(coordObj@meta.data) || nrow(coordObj@meta.data) == 0)
    coordObj@meta.data <- data.frame(row.names = genes)
  miss <- setdiff(genes, rownames(coordObj@meta.data))
  if (length(miss)) coordObj@meta.data[miss, ] <- NA

  if (is.null(cluster_name))
    cluster_name <- paste0("modL", gsub("\\.", "", format(L_min, trim = TRUE)))

  coordObj@meta.data[genes, cluster_name]            <- cluster_final
  coordObj@meta.data[genes, paste0(cluster_name, "_sub")] <- submod

  invisible(coordObj)
}