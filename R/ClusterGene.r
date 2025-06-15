#' @title Gene macro‑module clustering on Lee's *L* with optional Inter/Intra filter
#'
#' @description
#' Clusters genes into **macro‑modules** using an undirected graph built from
#' positive Lee's \emph{L} edges.  Three sequential filters are applied to the
#' \eqn{L} matrix:
#'
#' \enumerate{
#'   \item **Inter/Intra class mask** — If \code{cluster_by_delta_class} is set
#'         to "Inter" or "Intra", all edges whose corresponding entry in the
#'         \code{class} matrix (taken from \code{delta_layer}) is *not*
#'         “Inter‑cellular” or “Intra‑cellular” respectively are zeroed out.
#'   \item **Hard threshold** — Edges with \eqn{|L| < L_min} are removed; if
#'         \code{use_FDR = TRUE}, edges with \code{FDR > FDR_max} are also
#'         removed.
#'   \item **Quantile filter** — Optionally drop the weakest positive edges by
#'         calling \code{.filter_matrix_by_quantile()}.
#' }
#'
#' The remaining positive edges define a weighted undirected graph that is
#' partitioned by Louvain (default) or Leiden community detection.  Resulting
#' cluster IDs are stored in \code{coordObj@meta.data[[cluster_name]]}.  No
#' further sub‑clustering is performed.
#'
#' @param coordObj       A **CoordObj** returned by FG²CLI helpers.
#' @param grid_name      Name of the grid sub‑layer.  If \code{NULL} and only a
#'                       single grid layer exists, that layer is used
#'                       automatically.
#' @param lee_stats_layer List entry storing Lee's \emph{L} results
#'                       (default "LeeStats_Xz").  Must contain matrices
#'                       \code{$L} and (optionally) \code{$FDR}.
#' @param L_min          Minimum absolute Lee's \emph{L} for an edge to be kept
#'                       (default 0).
#' @param use_FDR        Apply FDR filter (\code{FDR_max}) on \eqn{L} edges?
#'                       (default \code{TRUE}).
#' @param FDR_max        Maximum FDR allowed when \code{use_FDR = TRUE}
#'                       (default 0.05).
#' @param pct_min        Quantile filter to drop weakest positive edges before
#'                       graph construction; pass "q0" to disable.
#' @param drop_isolated  Drop genes without any retained positive edge
#'                       (default \code{TRUE}).
#' @param algo           Community algorithm: "louvain" (default) or "leiden".
#' @param resolution     Leiden resolution parameter (ignored for Louvain).
#' @param delta_layer Name of the list entry that stores the \code{class}
#'                       matrix (default "LeeMinusR_bycell_stats").  Must
#'                       contain a character matrix \code{$class} with the same
#'                       dimnames as \code{$L}.
#' @param cluster_by_delta_class  One of "none" (default), "Inter", or "Intra".
#'                       When set to "Inter" or "Intra", only edges labelled as
#'                       “Inter‑cellular” or “Intra‑cellular” in the class
#'                       matrix are retained before further filtering.
#' @param cluster_name   Column name for the macro‑module IDs; a sensible
#'                       default is generated from \code{L_min} and the chosen
#'                       Inter/Intra filter.
#'
#' @return The modified \code{CoordObj} (invisibly).
#' @export
clusterGenes <- function(
  coordObj,
  grid_name              = NULL,
  lee_stats_layer        = "LeeStats_Xz",
  L_min                  = 0,
  use_FDR               = TRUE,
  FDR_max                = 0.05,
  drop_isolated          = TRUE,
  algo                   = c("louvain", "leiden"),
  resolution             = 1,
  pct_min                = "q0",
  filter_by_delta_class  = TRUE,
  delta_layer            = "LeeMinusR_bycell_stats",
  cluster_by_delta_class = c("Other", "Inter", "Intra"),
  cluster_name           = NULL
) {
  ## ---------- 0. Argument checks ----------
  algo  <- match.arg(algo)
  cluster_by_delta_class <- match.arg(cluster_by_delta_class)

  ## ---------- 1. Select grid layer ----------
  if (is.null(grid_name)) {
    subs <- names(coordObj@grid)
    if (length(subs) == 1L) grid_layer_name <- subs else
      stop("Multiple grid sub‑layers detected; please specify grid_name.")
  } else {
    grid_layer_name <- as.character(grid_name)
  }
  stopifnot(grid_layer_name %in% names(coordObj@grid))
  layer <- coordObj@grid[[grid_layer_name]]

  ## ---------- 2. Read Lee's L (+FDR) ----------
  if (!lee_stats_layer %in% names(layer))
    stop("Layer '", lee_stats_layer, "' not found under grid '", grid_layer_name, "'.")
  LeeStats <- layer[[lee_stats_layer]]
  Lmat <- LeeStats$L
  stopifnot(is.matrix(Lmat))
  genes_all <- rownames(Lmat)

  if (use_FDR) {
    FDRmat <- LeeStats$FDR
    stopifnot(is.matrix(FDRmat), all(dim(FDRmat) == dim(Lmat)))
  } else {
    FDRmat <- NULL
  }

  ## ---------- 3. Optional Inter/Intra mask ----------
  if (filter_by_delta_class) {
    if (!delta_layer %in% names(layer))
      stop("Layer '", delta_layer, "' not found under grid '", grid_layer_name, "'.")
    classMat <- layer[[delta_layer]]$class
    stopifnot(is.matrix(classMat), all(dim(classMat) == dim(Lmat)))

    target <- if (cluster_by_delta_class == "Inter") "Inter-cellular" else "Intra-cellular"

    keep_mask <- classMat == target
    keep_mask[is.na(keep_mask)] <- FALSE

    # Ensure symmetry by intersecting with its transpose
    keep_mask <- keep_mask & t(keep_mask)

    Lmat[!keep_mask] <- 0
    if (use_FDR) {
      FDRmat[!keep_mask] <- 1
    } # force drop later
  }

  ## ---------- 4. Build adjacency ----------
  A <- Lmat
  A[abs(A) < L_min] <- 0
  if (use_FDR) A[FDRmat > FDR_max] <- 0
  diag(A) <- 0
  A[A < 0] <- 0
  W <- .filter_matrix_by_quantile(A, pct_min, "q100")

  keep <- if (drop_isolated) rowSums(W != 0) > 0 else rep(TRUE, length(genes_all))
  kept_genes <- genes_all[keep]
  if (length(kept_genes) < 2L)
    stop("Fewer than two genes remain after filtering; relax thresholds or set drop_isolated = FALSE.")

  g_macro <- igraph::graph_from_adjacency_matrix(W[kept_genes, kept_genes],
                                                 mode = "undirected",
                                                 weighted = TRUE,
                                                 diag = FALSE)
  comm <- if (algo == "leiden") {
    igraph::cluster_leiden(g_macro,
      weights = igraph::E(g_macro)$weight,
      resolution_parameter = resolution
    )
  } else
    igraph::cluster_louvain(g_macro,
      weights = igraph::E(g_macro)$weight,
      resolution = resolution
    )

  macro <- setNames(rep(NA_integer_, length(genes_all)), genes_all)
  macro[names(igraph::membership(comm))] <- igraph::membership(comm)

  ## ---------- 5. Write back to meta.data ----------
  if (is.null(cluster_name)) {
    suffix <- switch(cluster_by_delta_class, none = "", Inter = "_Inter", Intra = "_Intra")
    cluster_name <- paste0("modL", gsub("\\.", "", format(L_min, trim = TRUE)), suffix)
  }

  if (is.null(coordObj@meta.data) || nrow(coordObj@meta.data) == 0)
    coordObj@meta.data <- data.frame(row.names = genes_all)
  miss <- setdiff(genes_all, rownames(coordObj@meta.data))
  if (length(miss)) coordObj@meta.data[miss, ] <- NA

  coordObj@meta.data[genes_all, cluster_name] <- macro

  invisible(coordObj)
}



#' @noRd
.parse_q <- function(qstr) {
  if (!is.character(qstr) || length(qstr) != 1 || !grepl("^[qQ][0-9]+\\.?[0-9]*$", qstr)) {
    stop("parse_q(): qstr must be a single character string in the format 'q0', 'q50', 'q100', etc.")
  }
  pct_num <- as.numeric(sub("^[qQ]", "", qstr)) / 100
  if (is.na(pct_num) || pct_num < 0 || pct_num > 1) {
    stop("parse_q(): Parsed percentage must be between 0 and 1.")
  }
  return(pct_num)
}

#' @noRd
.filter_matrix_by_quantile <- function(mat, pct_min = "q0", pct_max = "q100") {
  ## Chekk input matrix -----------------------------------------------
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("Input 'mat' must be a numeric matrix.")
  }
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("Input matrix must have rownames and colnames.")
  }

  ## Argument checks for quantiles -----------------------------------
  pmin <- .parse_q(pct_min)
  pmax <- .parse_q(pct_max)
  if (pmin > pmax)
    stop("'pct_min' cannot be greater than 'pct_max'.")

  vec <- as.vector(mat)
  lower_cut <- as.numeric(stats::quantile(vec, pmin, type = 7, na.rm = TRUE))
  upper_cut <- as.numeric(stats::quantile(vec, pmax, type = 7, na.rm = TRUE))

  ## Filter indices based on quantiles -----------------------------
  sel_idx <- which(vec >= lower_cut & vec <= upper_cut)

  ## Create new matrix with filtered values ----------------------
  new_mat <- matrix(0, nrow(mat), ncol(mat),
                    dimnames = dimnames(mat)) 
  new_mat[sel_idx] <- mat[sel_idx]

  new_mat
}


