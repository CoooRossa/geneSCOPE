#' @title Plot Gene Co-distribution Network with Cluster Coloring
#'
#' @description
#' Generates a publication-quality network plot in which edges represent
#' significant Lee’s _L_ spatial correlations between genes, and nodes are
#' coloured by a user-supplied cluster assignment (e.g. Louvain clusters).
#' Positive and negative correlations can be visualised with different
#' linetypes; hub genes (high-degree nodes) are outlined for emphasis, and
#' an optional QC caption summarises network statistics is printed below the
#' figure.
#'
#' @param coordObj       A **`CoordObj`** containing spatial‐grid statistics
#'                       produced by downstream FG²CLI helpers.
#' @param grid_name      Character. Name of the grid sub-layer (e.g.
#'                       `"grid_lenGrid50"`). If `NULL`, the single available
#'                       grid layer is used.
#' @param lee_stats_layer Character. Slot name that stores Lee’s _L_ results
#'                       inside the chosen grid sub-layer. Default `"LeeStats"`.
#' @param gene_subset    Optional character vector of gene names to keep.
#' @param cluster_vec    Either a named vector mapping genes to clusters, or a
#'                       single string giving the column name inside
#'                       `coordObj@meta.data` that contains cluster labels.
#' @param cluster_palette Character vector of colours (named or unnamed) used
#'                       to paint clusters. Unnamed palettes are re-cycled
#'                       with `colorRampPalette()`.
#' @param L_min          Minimum |L| retained when `weight_abs = TRUE`.
#' @param L_min_neg      Separate threshold for negative edges when
#'                       `weight_abs = FALSE`.  Default uses `L_min`.
#' @param p_cut          Optional *P*-value cutoff (used if a `P` matrix exists
#'                       inside the LeeStats layer).
#' @param drop_isolated  Logical. Drop nodes with zero degree after thresholding.
#' @param weight_abs     Logical. If `TRUE`, edge sign is ignored when filtering
#'                       (`|L| >= L_min`).  If `FALSE`, separate thresholds are
#'                       applied for positive/negative edges.
#' @param vertex_size    Base node size (passed to `ggplot2`).
#' @param base_edge_mult Multiplier controlling the thickest edge width.
#' @param label_cex      Node-label text size.
#' @param layout_niter   Number of iterations for the Fruchterman-Reingold
#'                       layout.
#' @param seed           Random seed for reproducible layouts.
#' @param hub_factor     Nodes with degree > `hub_factor × median(deg)` are
#'                       treated as hubs and drawn with a border.
#' @param max.overlaps   Passed to `geom_node_text(repel = TRUE)` to limit
#'                       label collisions.
#' @param L_max          Absolute |L| values larger than this are discarded.
#' @param hub_border_col Colour of hub borders.
#' @param hub_border_size Stroke width of hub borders.
#' @param show_sign      Logical. If `TRUE`, positive/negative edges are drawn
#'                       with different linetypes.
#' @param neg_linetype   Linetype used for negative edges when `show_sign = TRUE`.
#' @param neg_legend_lab Legend label for negative edges.
#' @param pos_legend_lab Legend label for positive edges.
#' @param show_qc_caption Logical. Appends a one-line QC caption if the LeeStats
#'                       layer contains a `$qc` list.
#' @param title          Optional plot title.
#' @import igraph
#' @import ggraph
#' @import ggplot2
#' @import dplyr
#'
#' @return A **`ggplot`** object (class `ggraph_tbl_graph`).
#'
#' @examples
#' \dontrun{
#' p <- plotNetworkClu(
#'   coordObj     = P5.coord,
#'   grid_name    = "grid_lenGrid50",
#'   cluster_vec  = P5.coord@meta.data$louvain_Xz_top0.1_grid50,
#'   gene_subset  = head(colnames(P5.coord@grid$grid_lenGrid50$Counts), 300)
#' )
#' print(p)
#' }
plotNetworkGene <- function(coordObj,
                           grid_name        = NULL,
                           lee_stats_layer  = "LeeStats",
                           gene_subset      = NULL,
                           cluster_vec      = NULL,
                           cluster_palette   = c(
                             "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00",
                             "#FFFF33", "#A65628", "#984EA3",
                             "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                             "#A6D854", "#FFD92F", "#E5C494"
                           ),
                           L_min            = 0.15,
                           L_min_neg        = NULL,
                           p_cut            = NULL,
                           drop_isolated    = TRUE,
                           weight_abs       = TRUE,
                           vertex_size      = 8,
                           base_edge_mult   = 12,
                           label_cex        = 3,
                           layout_niter     = 1000,
                           seed             = 1,
                           hub_factor       = 3,
                           max.overlaps     = 20,
                           L_max            = 1,
                           hub_border_col   = "#4B4B4B",
                           hub_border_size  = 0.8,
                           show_sign        = FALSE,
                           neg_linetype     = "dashed",
                           neg_legend_lab   = "Negative",
                           pos_legend_lab   = "Positive",
                           show_qc_caption  = TRUE,
                           title = NULL) {

  ## ——— Dependencies ———————————————————————————————————————————————— ##
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(dplyr)

  ## ——— 0. Determine the grid sub-layer ——————————————— ##
  if (is.null(coordObj@grid) || length(coordObj@grid) == 0L)
    stop("coordObj@grid is empty; cannot draw network.")

  if (is.null(grid_name)) {
    sub_layers <- names(coordObj@grid)
    if (length(sub_layers) == 1L) {
      grid_layer_name <- sub_layers[[1]]
    } else {
      stop("coordObj@grid contains multiple sublayers; please specify `grid_name`.")
    }
  } else {
    grid_layer_name <- as.character(grid_name)
  }
  if (!(grid_layer_name %in% names(coordObj@grid)))
    stop("Specified `grid_name` '", grid_layer_name, "' does not exist.")

  ## ——— 0.5. Normalise L_min_neg ——————————————— ##
  L_min_neg <- if (is.null(L_min_neg)) L_min else abs(L_min_neg)

  ## ——— 1. Check LeeStats layer ——————————————— ##
  if (is.null(coordObj@grid[[grid_layer_name]][[lee_stats_layer]]))
    stop("LeeStats layer ‘", lee_stats_layer, "’ not found. Did you run addLeeStats_fromGrid()?")

  res <- coordObj@grid[[grid_layer_name]][[lee_stats_layer]]

  ## ——— 2. Gene subset ————————————————————————————— ##
  keep_genes <- colnames(res$L)
  if (!is.null(gene_subset))
    keep_genes <- intersect(keep_genes, gene_subset)

  if (length(keep_genes) < 2L)
    stop("After subsetting, fewer than two genes remain; nothing to plot.")

  ## ——— 3. Threshold the L-matrix ————————————————————— ##
  A <- res$L[keep_genes, keep_genes, drop = FALSE]

  ## 3.1 P-value cut
  if (!is.null(p_cut) && !is.null(res$P)) {
    idx_keep <- match(keep_genes, rownames(res$L))
    Pmat <- res$P[idx_keep, idx_keep, drop = FALSE]
    A[Pmat >= p_cut | is.na(Pmat)] <- 0
  }

  ## 3.2 Remove extreme |L|
  A[abs(A) > L_max] <- 0

  ## 3.3 Apply L_min / L_min_neg
  if (weight_abs) {
    A[abs(A) < L_min] <- 0
  } else {
    A[A > 0  &  A <  L_min]      <- 0       # weak positives
    A[A < 0  &  abs(A) < L_min_neg] <- 0    # weak negatives
  }

  ## 3.4 Symmetrise
  A2 <- (A + t(A)) / 2
  diag(A2) <- 0

  ## ——— 4. Drop isolated nodes ——————————————— ##
  if (drop_isolated) {
    keep_nodes <- rowSums(abs(A2) > 0) > 0
    A2 <- A2[keep_nodes, keep_nodes, drop = FALSE]
  }
  if (nrow(A2) < 2L || all(A2 == 0))
    stop("No edges remain after filtering; relax thresholds.")

  ## ——— 5. Build igraph; mark edge sign ——————————— ##
  g <- graph_from_adjacency_matrix(abs(A2), mode = "undirected",
                                   weighted = TRUE, diag = FALSE)
  e_idx          <- get.edgelist(g, names = FALSE)
  E(g)$sign      <- ifelse(A2[e_idx] < 0, "neg", "pos")

  ## ——— 6. Flags for positive/negative neighbours ——— ##
  Vnames  <- V(g)$name
  has_pos <- has_neg <- setNames(rep(FALSE, length(Vnames)), Vnames)
  for (v in Vnames) {
    sigs <- E(g)$sign[incident(g, v)]
    has_pos[v] <- any(sigs == "pos")
    has_neg[v] <- any(sigs == "neg")
  }

  ## ——— 7. Cluster assignment ———————————————— ##
  clu              <- setNames(rep(NA_character_, length(Vnames)), Vnames)
  is_factor_input  <- FALSE
  factor_levels    <- NULL

  get_meta_cluster <- function(col) {
    if (is.null(coordObj@meta.data) || !(col %in% colnames(coordObj@meta.data)))
      stop("Cluster column ‘", col, "’ not found in coordObj@meta.data.")
    tmp <- coordObj@meta.data[[col]]
    names(tmp) <- rownames(coordObj@meta.data)
    tmp[Vnames]
  }

  if (!is.null(cluster_vec)) {
    if (length(cluster_vec) > 1L) {          # named vector
      if (is.factor(cluster_vec)) {
        is_factor_input <- TRUE; factor_levels <- levels(cluster_vec)
      }
      tmp <- as.character(cluster_vec[Vnames])
      clu[!is.na(tmp) & has_pos] <- tmp[!is.na(tmp) & has_pos]
    } else {                                 # column name
      meta_clu <- get_meta_cluster(cluster_vec)
      if (is.factor(coordObj@meta.data[[cluster_vec]])) {
        is_factor_input <- TRUE; factor_levels <- levels(meta_clu)
      }
      clu[!is.na(meta_clu) & has_pos] <- as.character(meta_clu[!is.na(meta_clu) & has_pos])
    }
  }

  ## ——— 8. Keep assigned nodes only ————————————— ##
  keep_nodes <- names(clu)[!is.na(clu)]
  if (length(keep_nodes) < 2L)
    stop("Fewer than two nodes are assigned to a cluster; cannot plot.")

  g      <- induced_subgraph(g, keep_nodes)
  Vnames <- V(g)$name
  clu    <- clu[Vnames]

  ## ——— 9–11. Degree, palette, edge colours ————————— ##
  deg_vec <- degree(g)

  uniq_clu <- if (is_factor_input && !is.null(factor_levels))
                intersect(factor_levels, unique(clu))
              else sort(unique(clu))
  n_clu <- length(uniq_clu)

  if (is.null(cluster_palette)) {
    palette_vals <- setNames(rainbow(n_clu), uniq_clu)
  } else if (is.null(names(cluster_palette))) {
    palette_vals <- setNames(colorRampPalette(cluster_palette)(n_clu), uniq_clu)
  } else {
    palette_vals <- cluster_palette
    missing <- setdiff(uniq_clu, names(palette_vals))
    if (length(missing) > 0L) {
      extra <- colorRampPalette(cluster_palette)(length(missing))
      palette_vals <- c(palette_vals, setNames(extra, missing))
    }
    palette_vals <- palette_vals[uniq_clu]
  }

  basecol <- setNames(rep("gray80", length(Vnames)), Vnames)
  basecol[!is.na(clu)] <- palette_vals[clu[!is.na(clu)]]

  ## edge colours (weight-based shading)
  e_idx       <- get.edgelist(g, names = FALSE)
  w_norm      <- E(g)$weight / max(E(g)$weight)
  cluster_ord <- setNames(seq_along(uniq_clu), uniq_clu)
  edge_cols   <- character(ecount(g))

  for (i in seq_len(ecount(g))) {
    v1 <- Vnames[e_idx[i, 1]]; v2 <- Vnames[e_idx[i, 2]]
    if (E(g)$sign[i] == "neg") {
      edge_cols[i] <- "gray40"
    } else {
      # choose reference colour
      ref_col <- {
        if (!is.na(clu[v1]) && !is.na(clu[v2]) && clu[v1] == clu[v2]) {
          basecol[v1]
        } else {
          d1 <- deg_vec[v1]; d2 <- deg_vec[v2]
          if (d1 > d2) basecol[v1]
          else if (d2 > d1) basecol[v2]
          else {
            if (is.na(clu[v1]) && !is.na(clu[v2]))       basecol[v2]
            else if (is.na(clu[v2]) && !is.na(clu[v1]))  basecol[v1]
            else if (!is.na(clu[v1]) && !is.na(clu[v2])) {
              if (cluster_ord[clu[v1]] <= cluster_ord[clu[v2]])
                basecol[v1] else basecol[v2]
            } else "gray80"
          }
        }
      }
      rgb_ref <- col2rgb(ref_col) / 255
      tval    <- w_norm[i]
      edge_cols[i] <- rgb(1 - tval + tval * rgb_ref[1],
                          1 - tval + tval * rgb_ref[2],
                          1 - tval + tval * rgb_ref[3])
    }
  }
  E(g)$edge_col <- edge_cols
  E(g)$linetype <- if (show_sign) E(g)$sign else "solid"

  ## ——— 12. Layout —————————————————————————— ##
  set.seed(seed)
  lay <- create_layout(g, layout = "fr", niter = layout_niter)
  lay$basecol <- basecol
  lay$deg     <- deg_vec
  lay$hub     <- deg_vec > hub_factor * median(deg_vec)

  ## ——— 13. QC caption —————————————————————— ##
  qc_txt <- if (show_qc_caption && !is.null(res$qc))
              with(res$qc, sprintf(
                "density = %.3f | comp = %d | Q = %.2f | mean-deg = %.1f ± %.1f | hubs = %.1f%% | sig.edge = %.1f%%",
                edge_density, components, modularity_Q,
                mean_degree, sd_degree, hub_ratio * 100, sig_edge_frac * 100))
            else NULL

  ## ——— 14. Draw ——————————————————————————— ##
  if (show_sign) {
    p <- ggraph(lay) +
      geom_edge_link(
        aes(width    = weight,
            colour   = edge_col,
            linetype = linetype),
        lineend = "round",
        show.legend = c(width = TRUE, linetype = TRUE, colour = FALSE)
      ) +
      scale_edge_width(name = "|L|", range = c(0.3, base_edge_mult)) +
      scale_edge_colour_identity() +
      ## ——— FIX: ensure correct mapping & order ——— ##
      scale_edge_linetype_manual(
        name   = "Direction",
        values = c(pos = "solid", neg = neg_linetype),
        breaks = c("pos", "neg"),
        labels = c(pos = pos_legend_lab, neg = neg_legend_lab)
      )
  } else {
    p <- ggraph(lay) +
      geom_edge_link(
        aes(width = weight, colour = edge_col),
        lineend = "round",
        show.legend = c(width = TRUE, colour = FALSE)
      ) +
      scale_edge_width(name = "|L|", range = c(0.3, base_edge_mult)) +
      scale_edge_colour_identity()
  }

  p <- p +
    geom_node_point(
      data = ~ dplyr::filter(.x, !hub),
      aes(size = deg, fill = basecol),
      shape = 21, stroke = 0,
      show.legend = c(fill = TRUE, size = TRUE)
    ) +
    geom_node_point(
      data = ~ dplyr::filter(.x, hub),
      aes(size = deg, fill = basecol),
      shape = 21, colour = hub_border_col,
      stroke = hub_border_size,
      show.legend = FALSE
    ) +
    geom_node_text(
      aes(label = name), size = label_cex, repel = TRUE,
      vjust = 1.4, max.overlaps = max.overlaps
    ) +
    scale_size_continuous(
      name  = "Node size",
      range = c(vertex_size / 2, vertex_size * 1.5)
    ) +
    scale_fill_identity(
      name  = "Cluster",
      guide = guide_legend(
        override.aes = list(shape = 21, size = vertex_size * 0.5),
        nrow = 1, order = 1
      ),
      breaks = unname(palette_vals),
      labels = names(palette_vals)
    ) +
    labs(title = title) +
    guides(
      fill = guide_legend(
        title = "Cluster",
        override.aes = list(shape = 21, size = vertex_size * 0.5),
        nrow = 1, order = 1
      ),
      size = guide_legend(
        title = "Node size",
        override.aes = list(shape = 21, fill = "grey80"),
        nrow = 1, order = 2
      ),
      edge_width = guide_legend(
        title = "|L|", direction = "horizontal", nrow = 1, order = 3
      ),
      linetype = if (show_sign)
        guide_legend(title = "Direction", nrow = 1, order = 4) else FALSE
    ) +
    coord_fixed() +
    theme_void(base_size = 15) +
    theme(
      legend.text       = element_text(size = 15),
      legend.title      = element_text(size = 20),
      plot.title        = element_text(size = 16, hjust = 0.5),
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.box.just   = "center",
      legend.spacing.y  = unit(0.2, "cm"),
      legend.box.margin = margin(t = 5, b = 5),
      plot.caption      = element_text(hjust = 0, size = 9,
                                       margin = margin(t = 6))
    ) +
    labs(caption = qc_txt)

  return(p)
}

#' @title Plot Gene Co-distribution Network with Cluster Coloring
#' @description
#' This function is a wrapper around `plotNetworkGene()` that uses a cluster
#' vector to colour nodes. It is intended for use with cluster assignments
#' such as Louvain clusters.
#' @param coordObj       A **`CoordObj`** containing spatial‐grid statistics
#'                      produced by downstream FG²CLI helpers.
#' @param grid_name      Character. Name of the grid sub-layer (e.g.
#'                     `"grid_lenGrid50"`). If `NULL`, the single available
#' #                      grid layer is used.
#' @param lee_stats_layer Character. Slot name that stores Lee’s _L_ results
#'                      inside the chosen grid sub-layer. Default `"LeeStats"`.
#' @param gene_subset    Optional character vector of gene names to keep.
#' @param cluster_vec    Either a named vector mapping genes to clusters, or a
#'                     single string giving the column name inside
#'                      `
#'  
#' 
#' @param cluster_palette Character vector of colours (named or unnamed) used
#'                     to paint clusters. Unnamed palettes are re-cycled
#' #                     with `colorRampPalette()`.
#' @param L_min          Minimum |L| retained when `weight_abs = TRUE`.
#' @param L_min_neg      Separate threshold for negative edges when
#'                     `weight_abs = FALSE`.  Default uses `L_min`.
#' @param p_cut          Optional *P*-value cutoff (used if a `P` matrix exists
#' #                     inside the LeeStats layer).
#' @param drop_isolated  Logical. Drop nodes with zero degree after thresholding.
#' @param weight_abs     Logical. If `TRUE`, edge sign is ignored when filtering
#' #                     (`|L| >= L_min`).  If `FALSE`, separate thresholds are
#' 
#' @import igraph
#' @import ggraph
#' @import ggplot2
#' @import dplyr
#' @param vertex_size    Base node size (passed to `ggplot2`).
#' @param base_edge_mult Multiplier controlling the thickest edge width.
#' @param edge_width_max Maximum edge width (default 5).
#' @param label_cex      Node-label text size.
#' @param layout_niter   Number of iterations for the Fruchterman-Reingold
#'                       layout.
#' @param seed           Random seed for reproducible layouts.
#' @param max.overlaps   Passed to `geom_node_text(repel = TRUE)` to limit
#' #                     label collisions.
#' @param L_max          Absolute |L| values larger than this are discarded.
#' @param show_sign      Logical. If `TRUE`, positive/negative edges are drawn
#' #                     with different linetypes.
#' @param pos_edge_col   Colour for positive edges.
#' @param neg_edge_col   Colour for negative edges.
#' @param neg_linetype   Linetype used for negative edges when `show_sign = TRUE`.
#' @param neg_legend_lab Legend label for negative edges.
#' @param pos_legend_lab Legend label for positive edges.
#' @param show_qc_caption Logical. Appends a one-line QC caption if the LeeStats
#' #                     layer contains a `$qc` list.
#' @param title          Optional plot title.
#' @return A **`ggplot`** object (class `ggraph_tbl_graph`).
#' @export
plotNetworkClu <- function(coordObj,
                               grid_name        = NULL,
                               lee_stats_layer  = "LeeStats",
                               gene_subset      = NULL,
                               cluster_vec      = NULL,
                               cluster_palette  = c(
                                 "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00",
                                 "#FFFF33", "#A65628", "#984EA3", 
                                 "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                                 "#A6D854", "#FFD92F", "#E5C494"
                               ),
                               L_min            = 0.15,
                               L_min_neg        = NULL,
                               p_cut            = NULL,
                               drop_isolated    = TRUE,
                               weight_abs       = TRUE,
                               vertex_size      = 8,
                               base_edge_mult   = 12,   # kept for backward‑compat (ignored)
                               edge_width_max   = 5,
                               label_cex        = 3,
                               layout_niter     = 1000,
                               seed             = 1,
                               max.overlaps     = 20,
                               L_max            = 1,
                               show_sign        = FALSE,
                               pos_edge_col     = "#E41A1C",
                               neg_edge_col     = "#377EB8",
                               neg_linetype     = "dashed",
                               neg_legend_lab   = "Negative",
                               pos_legend_lab   = "Positive",
                               show_qc_caption  = TRUE,
                               title = NULL) {
  library(igraph)
  library(ggraph)
  library(dplyr)

  ## 0. Resolve grid sub‑layer ------------------------------------------------
  if (is.null(coordObj@grid) || length(coordObj@grid) == 0L)
    stop("coordObj@grid is empty; cannot draw network.")
  if (is.null(grid_name)) {
    sub_layers <- names(coordObj@grid)
    if (length(sub_layers) == 1L) {
      grid_layer_name <- sub_layers[[1]]
    } else {
      stop("coordObj@grid contains multiple sub‑layers; please set `grid_name`.")
    }
  } else {
    grid_layer_name <- as.character(grid_name)
  }
  if (!(grid_layer_name %in% names(coordObj@grid)))
    stop("`grid_name` ", grid_layer_name, " not found in coordObj@grid.")

  ## 0.5 Default L_min_neg ----------------------------------------------------
  if (is.null(L_min_neg)) L_min_neg <- L_min

  ## 1. Retrieve Lee’s L ------------------------------------------------------
  res <- coordObj@grid[[grid_layer_name]][[lee_stats_layer]]
  if (is.null(res))
    stop("LeeStats layer ‘", lee_stats_layer, "’ not found in grid sub‑layer.")
  Lmat <- res$L
  gene_names <- colnames(Lmat)
  dimnames(res$P) <- dimnames(res$L)  # ensure P matrix matches Lmat

  ## 2. Optional gene subset --------------------------------------------------
  if (!is.null(gene_subset)) {
    keep_genes <- intersect(gene_names, gene_subset)
    if (length(keep_genes) < 2L)
      stop("After subsetting, fewer than 2 genes remain; cannot draw network.")
    Lmat <- Lmat[keep_genes, keep_genes, drop = FALSE]
    gene_names <- colnames(Lmat)
  }

  ## 3. Cluster assignment ----------------------------------------------------
  gene_clu <- rep(NA_character_, length(gene_names))
  names(gene_clu) <- gene_names

  if (!is.null(cluster_vec)) {
    if (length(cluster_vec) == 1L && is.character(cluster_vec)) {
      colname <- cluster_vec
      if (is.null(coordObj@meta.data) || !(colname %in% colnames(coordObj@meta.data)))
        stop("Column ‘", colname, "’ not found in coordObj@meta.data.")
      tmp <- coordObj@meta.data[[colname]]
      names(tmp) <- rownames(coordObj@meta.data)
      gene_clu[gene_names] <- as.character(tmp[gene_names])
    } else {  # named vector
      gene_clu[gene_names] <- as.character(cluster_vec[gene_names])
    }
  }
  keep <- !is.na(gene_clu)
  if (sum(keep) < 2L)
    stop("Fewer than 2 genes have valid cluster labels.")
  gene_clu <- gene_clu[keep]
  Lmat     <- Lmat[keep, keep, drop = FALSE]


  ## 4. Optional p‑value filter ----------------------------------------------

  if (!is.null(p_cut) && !is.null(res$P)) {
    Pmat <- res$P[rownames(Lmat), colnames(Lmat), drop = FALSE]
    Lmat[Pmat >= p_cut | is.na(Pmat)] <- NA
  }

  ## 5. L_max filter ----------------------------------------------------------
  Lmat[abs(Lmat) > L_max] <- NA

  ## 6. L_min / L_min_neg threshold ------------------------------------------
  if (weight_abs) {
    Lmat[abs(Lmat) < L_min] <- NA
  } else {
    Lmat[Lmat > 0 & Lmat <  L_min]          <- NA
    Lmat[Lmat < 0 & abs(Lmat) < L_min_neg]  <- NA
  }

  ## 7. Build cluster‑pair mean matrix ----------------------------------------
  df <- as.data.frame(as.table(Lmat), stringsAsFactors = FALSE)
  colnames(df) <- c("g1", "g2", "L")
  df$clu1 <- gene_clu[df$g1]
  df$clu2 <- gene_clu[df$g2]
  df <- df[df$g1 != df$g2 & !is.na(df$L), ]

  meanL_df <- df %>%
    group_by(clu1, clu2) %>%
    summarise(meanL = mean(L, na.rm = TRUE), .groups = "drop")

  clusters <- sort(unique(c(meanL_df$clu1, meanL_df$clu2)))
  Cmat <- matrix(0, length(clusters), length(clusters),
                 dimnames = list(clusters, clusters))
  for (i in seq_len(nrow(meanL_df)))
    Cmat[meanL_df$clu1[i], meanL_df$clu2[i]] <- meanL_df$meanL[i]
  Cmat <- (Cmat + t(Cmat)) / 2   # symmetrise
  Cmat[is.na(Cmat)] <- 0

  if (all(Cmat == 0))
    stop("All cluster‑pair mean L values are zero after filtering.")

  ## 8. Optionally drop isolated clusters ------------------------------------
  if (drop_isolated) {
    keep_idx <- rowSums(abs(Cmat) > 0) > 0
    Cmat <- Cmat[keep_idx, keep_idx, drop = FALSE]
    clusters <- rownames(Cmat)
  }
  if (nrow(Cmat) < 2L) stop("Fewer than 2 clusters with connections remain.")

  ## 9. igraph object ---------------------------------------------------------
  W <- abs(Cmat)
  g <- graph_from_adjacency_matrix(W, mode = "undirected",
                                   weighted = TRUE, diag = FALSE)
  e_idx <- get.edgelist(g, names = FALSE)
  E(g)$sign <- ifelse(Cmat[e_idx] < 0, "neg", "pos")

  ## 10. Layout prep ----------------------------------------------------------
  set.seed(seed)
  lay <- create_layout(g, layout = "fr", niter = layout_niter)

  # palette
  if (is.null(names(cluster_palette))) {
    clu_cols <- setNames(colorRampPalette(cluster_palette)(length(clusters)),
                         clusters)
  } else {
    clu_cols <- cluster_palette[clusters]
  }
  lay$fill <- clu_cols[lay$name]

  # node size = gene counts (numeric)
  gene_counts <- table(gene_clu)
  lay$size <- as.numeric(gene_counts[lay$name])

  ## 11. Plot -----------------------------------------------------------------
  if (show_sign) {
    edge_layer <- geom_edge_link(
      aes(width = weight,
          linetype = sign,
          colour   = sign),
      lineend     = "round",
      show.legend = c(width = TRUE, linetype = TRUE, colour = TRUE)
    )
    linetype_scale <- scale_edge_linetype_manual(
      name   = "Direction",
      values = c(pos = "solid", neg = neg_linetype),
      breaks = c("pos", "neg"),
      labels = c(pos_legend_lab, neg_legend_lab)
    )
    colour_scale <- scale_edge_colour_manual(
      name   = "Direction",
      values = c(pos = pos_edge_col, neg = neg_edge_col),
      breaks = c("pos", "neg"),
      labels = c(pos_legend_lab, neg_legend_lab)
    )
  } else {
    edge_layer <- geom_edge_link(
      aes(width = weight),
      colour      = "grey50",
      lineend     = "round",
      show.legend = c(width = TRUE)
    )
    linetype_scale <- NULL
  }

  p <- ggraph(lay) +
    edge_layer +
    labs(title = title) +
    scale_edge_width(name = "|mean L|",
                     range = c(0.2, edge_width_max)) +
    linetype_scale +
    colour_scale +
    geom_node_point(aes(size = size, fill = fill),
                    shape = 21, colour = "black", stroke = 0.3,
                    show.legend = c(fill = TRUE, size = TRUE)) +
    geom_node_text(aes(label = name),
                   repel = TRUE, vjust = 1.4,
                   size = label_cex, max.overlaps = max.overlaps) +
    scale_size_continuous(name = "Genes in cluster",
                          range = c(vertex_size / 2, vertex_size * 1.5)) +
    scale_fill_identity(name = "Cluster") +
    theme_void(base_size = 15) +
    theme(
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.box.just   = "center",
      legend.text       = element_text(size = 15),
      legend.title      = element_text(size = 20),
      plot.title        = element_text(size = 16, hjust = 0.5),
    ) +
    coord_fixed()

  return(p)
}


