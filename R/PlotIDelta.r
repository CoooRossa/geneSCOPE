#' @title Plot Iδ by Cluster
#' @param coordObj        A `coordObj` containing either
#'                        `coordObj@grid[[grid_name]]$iDeltaStats$delta_raw` or `delta_nor`,
#'                        or a column named `<grid_name>_iDelta_raw` / `<grid_name>_iDelta_nor` in
#'                        `coordObj@meta.data`. Additionally, `coordObj@meta.data` must have a column
#'                        indicating each gene’s cluster labels (`cluster_col`).
#' @param grid_name       Character, name of the grid sublayer to use (e.g., `"grid_lenGrid50"`).
#'                        If NULL and `coordObj@grid` has exactly one sublayer, that sublayer is used;
#'                        otherwise, this parameter is required.
#' @param cluster_col     Character, column name in `coordObj@meta.data` that indicates gene cluster labels.
#' @param suffix          Character, either `"raw"` or `"nor"` (default `"raw"`). Determines which Iδ values
#'                        to use (`delta_raw` or `delta_nor`).
#' @param top_n           Integer (optional). If provided, only the top `top_n` genes by Iδ value in each
#'                        cluster are displayed (default NULL, no filtering).
#' @param cluster_palette Named character vector (optional), mapping cluster names to colors. If not provided,
#'                        colors are generated automatically and extended if needed.
#' @param point_size      Numeric, size of the points (default 3).
#' @param line_size       Numeric, thickness of the connecting lines (default 0.5).
#' @param label_size      Numeric, font size for gene labels (default 3).
#'
#' @return A `ggplot2` object with one horizontal row of panels—each panel corresponding to one cluster—where
#'         panel widths scale proportionally to the number of genes in that cluster. Points and connecting lines
#'         show each gene’s Iδ value, and each gene is labeled by name. Gene labels at the far left/right are not clipped.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_text facet_grid scale_color_manual scale_x_discrete labs theme_minimal theme element_text element_blank element_rect unit coord_cartesian
#' @importFrom dplyr filter group_by arrange mutate ungroup tally
#' @export
plotIDeltaByCluster <- function(coordObj,
                                grid_name = NULL,
                                cluster_col,
                                suffix = c("raw", "nor"),
                                top_n = NULL,
                                cluster_palette   = c(
                                  "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00",
                                  "#FFFF33", "#A65628", "#984EA3",
                                  "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
                                  "#A6D854", "#FFD92F", "#E5C494"
                                ),
                                point_size = 3,
                                line_size = 0.5,
                                label_size = 3,
                                subCluster = NULL) {
  library(ggplot2)
  library(dplyr)
  library(scales)  # for colorRampPalette

  # —— 0. Validate parameters & determine suffix —— 
  suffix <- match.arg(suffix) 
  delta_key <- paste0("delta_", suffix) 
  meta_key  <- paste0(grid_name, "_iDelta_", suffix) 

  # —— 1. Determine which grid sublayer to use —— 
  if (is.null(grid_name)) {
    sub_layers <- names(coordObj@grid)
    if (length(sub_layers) == 1L) {
      grid_layer_name <- sub_layers[[1]]
    } else {
      stop("coordObj@grid contains multiple sublayers; please specify `grid_name` explicitly.")
    }
  } else {
    grid_layer_name <- as.character(grid_name)
  }
  if (!(grid_layer_name %in% names(coordObj@grid))) {
    stop("Specified `grid_name` '", grid_layer_name, "' does not exist in coordObj@grid.")
  }

  # —— 2. Read Iδ values from iDeltaStats or fallback to meta.data —— 
  iDeltaStats <- coordObj@grid[[grid_layer_name]]$iDeltaStats
  if (!is.null(iDeltaStats) && !is.null(iDeltaStats[[delta_key]])) {
    # Prefer to read directly from iDeltaStats if available
    delta_vals <- iDeltaStats[[delta_key]]
    genes <- names(delta_vals)
  } else {
    # Otherwise, read from coordObj@meta.data
    if (is.null(coordObj@meta.data) || !(meta_key %in% colnames(coordObj@meta.data))) {
      stop("Cannot find column '", delta_key, "' in coordObj@grid[['", grid_layer_name,
           "']]$iDeltaStats, nor column '", meta_key, "' in coordObj@meta.data.")
    }
    delta_vals <- coordObj@meta.data[[meta_key]]
    genes <- rownames(coordObj@meta.data)
  }

  # Verify that the specified cluster_col exists in meta.data
  if (is.null(coordObj@meta.data) || !(cluster_col %in% colnames(coordObj@meta.data))) {
    stop("coordObj@meta.data is missing column: ", cluster_col)
  }

  # —— 3. Construct a data frame with gene, delta, and cluster information —— 
  meta <- coordObj@meta.data
  if (!all(genes %in% rownames(meta))) {
    stop("coordObj@meta.data must have rownames matching all genes (names of `delta_vals`).")
  }
  df <- data.frame(
    gene    = genes,
    delta   = as.numeric(delta_vals),
    cluster = meta[genes, cluster_col],
    stringsAsFactors = FALSE
  )
  
  # ---- Ensure cluster is treated as DISCRETE ----
  # If the supplied cluster column is numeric (or integer), convert it to factor
  # so that ggplot treats it as a discrete aesthetic for scale_color_manual().
  if (is.numeric(df$cluster) || is.integer(df$cluster)) {
    df$cluster <- factor(df$cluster)
  } else {
    df$cluster <- as.factor(df$cluster)
  }

    if(!is.null(subCluster)){
        df[df$cluster %in% subCluster, ] -> df
    }



  # —— 4. Exclude genes with missing cluster labels —— 
  df <- df %>% filter(!is.na(cluster))

  # —— 5. If top_n is specified, keep only the top_n genes by delta in each cluster —— 
  if (!is.null(top_n) && top_n > 0) {
    df <- df %>%
      group_by(cluster) %>%
      arrange(desc(delta)) %>%
      slice_head(n = top_n) %>%
      ungroup()
  }

  # —— 6. Remove clusters that contain only a single gene —— 
  count_per_clu <- df %>%
    group_by(cluster) %>%
    tally(name = "n_genes") %>%
    ungroup()
  valid_clu <- count_per_clu$cluster[count_per_clu$n_genes > 1]
  df <- df %>% filter(cluster %in% valid_clu)

  if (nrow(df) == 0) {
    stop("No cluster has more than one gene; nothing to plot.")
  }

  # —— 7. Within each remaining cluster, sort genes by descending delta and create an ordered factor —— 
  df <- df %>%
    group_by(cluster) %>%
    arrange(desc(delta)) %>%
    mutate(gene_ord = factor(gene, levels = unique(gene))) %>%
    ungroup()

  # —— 8. Determine the ordering of clusters (uniq_clu) —— 
  #     Priority: if original cluster_col is a factor, use its levels;
  #               else if all cluster names can be coerced to numeric, sort numerically;
  #               otherwise, sort alphabetically.

  # Extract the original cluster vector from meta.data
  orig_cluster_vec <- coordObj@meta.data[[cluster_col]]
  is_factor_input <- is.factor(orig_cluster_vec)
  factor_levels   <- if (is_factor_input) levels(orig_cluster_vec) else NULL

  observed_clu <- unique(as.character(df$cluster))  # current clusters present

  if (is_factor_input && !is.null(factor_levels)) {
    # Keep only those factor levels that actually appear in df, in the same order
    uniq_clu <- intersect(factor_levels, observed_clu)
  } else {
    # Attempt to coerce cluster names to numeric
    suppress_ex <- function(x) suppressWarnings(as.numeric(x))
    numeric_vals <- suppress_ex(observed_clu)
    if (!any(is.na(numeric_vals))) {
      # All cluster names are numeric → sort by numeric value
      sorted_idx <- order(numeric_vals)
      uniq_clu <- observed_clu[sorted_idx]
    } else {
      # Not a factor and not purely numeric → alphabetical sort
      uniq_clu <- sort(observed_clu)
    }
  }
  n_clu <- length(uniq_clu)

  # Build a color palette in which the names of palette_vals follow uniq_clu
  if (is.null(cluster_palette)) {
    palette_vals <- setNames(rainbow(n_clu), uniq_clu)
  } else {
    if (is.null(names(cluster_palette))) {
      palette_vals <- setNames(colorRampPalette(cluster_palette)(n_clu), uniq_clu)
    } else {
      provided <- names(cluster_palette)
      if (!all(uniq_clu %in% provided)) {
        missing_clu <- setdiff(uniq_clu, provided)
        extra_cols  <- colorRampPalette(cluster_palette)(length(missing_clu))
        palette_vals <- c(cluster_palette, setNames(extra_cols, missing_clu))
        palette_vals <- palette_vals[uniq_clu]
      } else {
        palette_vals <- cluster_palette[uniq_clu]
      }
    }
  }
  df$col <- palette_vals[as.character(df$cluster)]

  # —— 9. Build the ggplot object —— 
  #     Use facet_grid(. ~ cluster, scales="free_x", space="free_x") so that each cluster panel
  #     scales horizontally to the number of genes. Use coord_cartesian(clip="off") to avoid
  #     clipping labels at the panel edges.
  y_label   <- "Iδ"
  plot_title <- paste0("Iδ by Cluster (", grid_layer_name, ")")

  p <- ggplot(df, aes(x = gene_ord, y = delta, group = cluster)) +
    # Draw points, using cluster color and including in legend
    geom_point(aes(color = cluster), size = point_size, show.legend = TRUE) +
    # Draw lines connecting points within each cluster
    geom_line(aes(color = cluster), size = line_size, show.legend = TRUE) +
    # Draw gene labels (not included in legend)
    geom_text(aes(label = gene, color = cluster),
              vjust = -1, size = label_size,
              check_overlap = TRUE,
              show.legend = FALSE) +
    facet_grid(. ~ cluster, scales = "free_x", space = "free_x", switch = "x") +
    scale_color_manual(values = palette_vals, name = "Cluster") +
    scale_x_discrete(expand = c(0, 0)) +
    labs(
      x = "Cluster",
      y = y_label,
      title = plot_title
    ) +
    theme_minimal(base_size = 8) +
    theme(
      strip.placement    = "outside",               # place facet labels outside the panels
      strip.text.x       = element_text(face = "bold"),
      axis.text.x        = element_blank(),          # hide x-axis text
      panel.title       = element_text(size = 10, face = "bold"),
      axis.ticks.x       = element_blank(),          # hide x-axis ticks
      panel.spacing.x    = unit(0.5, "lines"),       # small horizontal spacing between panels
      panel.grid.major.x = element_blank(),          # remove vertical grid lines
      panel.grid.minor.x = element_blank(),
      plot.margin        = unit(c(0.5, 1, 0.5, 1), "lines"),
      legend.position = "bottom"  # give extra margin for labels
    ) +
    coord_cartesian(clip = "off")  # ensure labels near panel borders are not clipped

  return(p)
}

