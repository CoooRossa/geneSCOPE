#' @title Plot Distribution of Lee’s L Values
#' @description
#'   Extract the Lee’s L matrix from a specified grid sublayer and LeeStats layer within `coordObj`.
#'   Bin the off-diagonal elements of the Lee’s L matrix into intervals and plot a histogram to
#'   visualize the distribution of Lee’s L values. Optionally, take the absolute value of Lee’s L
#'   before plotting. The plot adheres to publication-quality aesthetic standards.
#'
#' @param coordObj         An object in which Lee’s L has already been computed and stored under
#'                         `coordObj@grid[[grid_name]][[lee_stats_layer]]$L`.
#' @param grid_name        Character (optional). The name of the grid sublayer to use (e.g., `"grid_lenGrid50"`).
#'                         If NULL and there is only one sublayer under `coordObj@grid`, that layer is used.
#'                         Otherwise, this parameter is required.
#' @param lee_stats_layer  Character (optional). Name of the LeeStats layer to use (e.g., `"LeeStats_Xz"`).
#'                         If NULL, attempts to automatically detect a single layer whose name starts
#'                         with `"LeeStats_"`. If multiple matching layers are found, an error is thrown
#'                         and the user must specify `lee_stats_layer`.
#' @param bins            Integer or numeric vector. Passed to `ggplot2::geom_histogram()` as bin settings.
#'                         If a single integer, that number of equal-width bins is used. If a numeric
#'                         vector of breakpoints, those defines the histogram bins.
#' @param xlim            Numeric vector of length 2 (optional). X-axis display limits (c(min, max)).
#'                         Defaults to NULL (automatic scaling).
#' @param title           Character (optional). The plot title. Defaults to `"Lee’s L Distribution"` or
#'                         `"Absolute Lee’s L Distribution"` if `use_abs = TRUE`.
#' @param use_abs         Logical. Whether to take the absolute value of Lee’s L before plotting.
#'                         Defaults to FALSE (plot original Lee’s L).
#'
#' @return A `ggplot2` object displaying the histogram of Lee’s L values (or their absolute values).
#'
#' @importFrom ggplot2 ggplot aes geom_histogram labs theme_bw theme element_rect element_line element_text coord_cartesian
#' @importFrom grid unit
#' @export
plotLeeLDistribution <- function(coordObj,
                                 grid_name       = NULL,
                                 lee_stats_layer = NULL,
                                 bins            = 30,
                                 xlim            = NULL,
                                 title           = NULL,
                                 use_abs         = FALSE) {
  # —— 0. Determine grid_layer_name —— 
  if (is.null(coordObj@grid) || length(coordObj@grid) == 0) {
    stop("coordObj@grid is empty; cannot extract Lee’s L information.")
  }
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
  if (!(grid_layer_name %in% names(coordObj@grid))) {
    stop("Specified grid_name does not exist in coordObj@grid: ", grid_layer_name)
  }
  grid_layer <- coordObj@grid[[grid_layer_name]]

  # —— 1. Determine the LeeStats layer to use —— 
  if (is.null(lee_stats_layer)) {
    candidate_layers <- grep("^LeeStats_", names(grid_layer), value = TRUE)
    if (length(candidate_layers) == 1L) {
      lee_stats_layer <- candidate_layers
    } else if (length(candidate_layers) == 0L) {
      stop("No layers starting with 'LeeStats_' found under coordObj@grid[['", grid_layer_name, "']]. Please run Lee’s L first.")
    } else {
      stop("Multiple LeeStats layers found: ", paste(candidate_layers, collapse = ", "),
           ". Please specify lee_stats_layer explicitly.")
    }
  }
  if (!(lee_stats_layer %in% names(grid_layer))) {
    stop("Specified lee_stats_layer does not exist: ", lee_stats_layer)
  }

  # —— 2. Check that LeeStats and L matrix exist —— 
  LeeStats <- grid_layer[[lee_stats_layer]]
  if (is.null(LeeStats) || is.null(LeeStats$L)) {
    stop("coordObj@grid[['", grid_layer_name, "']][['", lee_stats_layer, "']]$L is empty; run Lee’s L calculation first.")
  }
  Lmat <- LeeStats$L
  if (!is.matrix(Lmat)) {
    stop("LeeStats$L must be a matrix.")
  }

  # —— 3. Extract off-diagonal elements and apply absolute value if requested —— 
  raw_vals <- as.numeric(Lmat[upper.tri(Lmat)])
  if (use_abs) {
    Lvals <- abs(raw_vals)
  } else {
    Lvals <- raw_vals
  }
  df <- data.frame(LeeL = Lvals)

  # —— 4. Set default title if not provided —— 
  if (is.null(title)) {
    if (use_abs) {
      title <- "Absolute Lee’s L Distribution"
    } else {
      title <- "Lee’s L Distribution"
    }
  }

  # —— 5. Build ggplot with publication-quality aesthetics —— 
  library(ggplot2)
  library(grid)

  p <- ggplot(df, aes(x = LeeL)) +
    {
      if (length(bins) == 1L && bins %% 1L == 0L) {
        geom_histogram(
          bins  = bins,
          fill  = "white",
          color = "black",
          size  = 0.3
        )
      } else if (is.numeric(bins) && length(bins) > 1L) {
        geom_histogram(
          breaks = bins,
          fill   = "white",
          color  = "black",
          size   = 0.3
        )
      } else {
        stop("Argument `bins` must be a single integer or a numeric vector of breakpoints.")
      }
    } +
    labs(
      title = title,
      x     = if (use_abs) expression("|Lee’s L|") else expression("Lee’s L"),
      y     = "Frequency"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.background   = element_rect(fill = "#c0c0c0", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "#c0c0c0", size = 0.2),
      panel.border       = element_rect(colour = "black", fill = NA, size = 0.5),
      axis.line          = element_line(colour = "black", size = 0.3),
      axis.ticks         = element_line(colour = "black", size = 0.3),
      axis.text          = element_text(size = 10, colour = "black"),
      axis.title         = element_text(size = 10, colour = "black", face = "plain"),
      plot.margin        = unit(c(1, 1, 1, 1), "cm"),
      legend.background  = element_rect(fill = "#c0c0c0", colour = NA),
      legend.key         = element_rect(fill = "#c0c0c0", colour = "black", size = 0.3),
      legend.title       = element_text(size = 11, colour = "black", face = "plain"),
      legend.text        = element_text(size = 10, colour = "black")
    ) +
    coord_cartesian(expand = FALSE)

  # —— 6. Apply x-axis limits if provided —— 
  if (!is.null(xlim)) {
    if (length(xlim) != 2L || !is.numeric(xlim)) {
      stop("`xlim` must be a numeric vector of length 2 (c(min, max)).")
    }
    p <- p + coord_cartesian(xlim = xlim)
  }

  return(p)
}
#' @title Plot Lee’s L vs. Gene Total Count Fold Change (All Gene Pairs)
#' @description
#'   For a given `coordObj`, extract the Lee’s L matrix from
#'   `coordObj@grid[[grid_name]][[lee_stats_layer]]$L`.
#'   Summarize each gene’s total count across all grid cells from
#'   `coordObj@grid[[grid_name]]$counts`. For every gene pair (i, j) with i < j,
#'   compute Lee’s L and the fold change of their total counts
#'   (fold change = max(total_i, total_j) / min(total_i, total_j), ≥ 1).
#'   Finally, plot a scatterplot with Lee’s L on the x-axis and fold change on
#'   the y-axis, using a 25% gray background and publication‐quality aesthetics.
#'
#' @param coordObj        An object containing
#'                        `coordObj@grid[[grid_name]]$counts` and
#'                        `coordObj@grid[[grid_name]][[lee_stats_layer]]$L`.
#' @param grid_name       Character (optional). The name of the grid sublayer to use.
#'                        If NULL and `coordObj@grid` has exactly one sublayer,
#'                        that sublayer is used automatically. Otherwise, it must be provided.
#' @param lee_stats_layer Character (optional). The name of the LeeStats layer (e.g., `"LeeStats_Xz"`).
#'                        If NULL, the function searches for exactly one layer under
#'                        `coordObj@grid[[grid_name]]` whose name starts with `"LeeStats_"`.
#'                        If multiple matches are found, an error is thrown and the user must
#'                        supply `lee_stats_layer` explicitly.
#' @param title           Character (optional). Plot title. Defaults to
#'                        `"Lee’s L vs. Fold Change (All Gene Pairs)"` if NULL.
#'
#' @return A `ggplot2` object showing a scatterplot of Lee’s L vs. fold change
#'         for all gene pairs (i < j).
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_rect element_line element_text coord_cartesian
#' @importFrom dplyr group_by summarise filter left_join rename rowwise ungroup
#' @importFrom grid unit
#' @export
plotLeeLScatter <- function(coordObj,
                                    grid_name       = NULL,
                                    lee_stats_layer = NULL,
                                    title           = NULL) {
  # —— 0. Determine grid_layer_name —— 
  if (is.null(coordObj@grid) || length(coordObj@grid) == 0L) {
    stop("coordObj@grid is empty; cannot extract data.")
  }
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
  if (!(grid_layer_name %in% names(coordObj@grid))) {
    stop("Specified grid_name does not exist: ", grid_layer_name)
  }
  grid_layer <- coordObj@grid[[grid_layer_name]]

  # —— 1. Determine lee_stats_layer —— 
  if (is.null(lee_stats_layer)) {
    candidate_layers <- grep("^LeeStats_", names(grid_layer), value = TRUE)
    if (length(candidate_layers) == 1L) {
      lee_stats_layer <- candidate_layers
    } else if (length(candidate_layers) == 0L) {
      stop("No layers starting with 'LeeStats_' found under coordObj@grid[['",
           grid_layer_name, "']]. Please run Lee’s L calculation first.")
    } else {
      stop("Multiple LeeStats layers found: ", paste(candidate_layers, collapse = ", "),
           ". Please specify lee_stats_layer explicitly.")
    }
  }
  if (!(lee_stats_layer %in% names(grid_layer))) {
    stop("Specified lee_stats_layer does not exist: ", lee_stats_layer)
  }

  # —— 2. Extract Lee’s L matrix —— 
  LeeStats <- grid_layer[[lee_stats_layer]]
  if (is.null(LeeStats) || is.null(LeeStats$L)) {
    stop("coordObj@grid[['", grid_layer_name, "']][['",
         lee_stats_layer, "']]$L is missing; please check.")
  }
  Lmat <- LeeStats$L
  if (!is.matrix(Lmat)) {
    stop("LeeStats$L must be a matrix.")
  }

  # —— 3. Summarize each gene’s total count across the entire grid —— 
  counts_df <- grid_layer$counts
  if (!all(c("gene", "count") %in% colnames(counts_df))) {
    stop("coordObj@grid[['", grid_layer_name, "']]$counts must contain columns 'gene' and 'count'.")
  }
  gene_totals <- counts_df |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      total_count = sum(count, na.rm = TRUE),
      .groups     = "drop"
    )

  # —— 4. Build all unique gene pairs (i < j) from the intersection of genes in Lmat and totals —— 
  genes_L <- rownames(Lmat)
  genes_T <- gene_totals$gene
  genes_common <- intersect(genes_L, genes_T)
  if (length(genes_common) < 2L) {
    stop("Fewer than 2 overlapping genes; cannot plot scatter.")
  }
  genes_sorted <- sort(genes_common)

  comb_df <- tidyr::crossing(geneA = genes_sorted, geneB = genes_sorted) |>
    dplyr::filter(geneA < geneB)

  # —— 5. Merge in total counts and Lee’s L, then compute fold change —— 
  comb_df <- comb_df |>
    dplyr::left_join(gene_totals, by = c("geneA" = "gene")) |>
    dplyr::rename(totalA = total_count) |>
    dplyr::left_join(gene_totals, by = c("geneB" = "gene")) |>
    dplyr::rename(totalB = total_count) |>
    dplyr::filter(!is.na(totalA), !is.na(totalB)) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      LeeL = Lmat[geneA, geneB],
      FoldRatio = ifelse(
        totalA >= totalB,
        totalA / totalB,
        totalB / totalA
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(LeeL), is.finite(FoldRatio))

  # —— 6. Set default title if not provided —— 
  if (is.null(title)) {
    title <- "Lee’s L vs. Fold Change (All Gene Pairs)"
  }

  # —— 7. Plot scatterplot —— 
  library(ggplot2)
  library(grid)

  p_all <- ggplot2::ggplot(comb_df, aes(x = LeeL, y = FoldRatio)) +
    ggplot2::geom_point(
      shape  = 21,
      size   = 1.2,
      fill   = "white",
      color  = "black",
      stroke = 0.2
    ) +
            scale_x_continuous(expand = expansion(mult = 0.05)) +   
         scale_y_continuous(expand = expansion(mult = 0.05)) +   
    ggplot2::labs(
      title = title,
      x     = expression("Lee’s L"),
      y     = "Fold Change (≥1)"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.background   = element_rect(fill = "#c0c0c0", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "#c0c0c0", size = 0.2),
      panel.border       = element_rect(colour = "black", fill = NA, size = 0.5),
      axis.line          = element_line(colour = "black", size = 0.3),
      axis.ticks         = element_line(colour = "black", size = 0.3),
      axis.text          = element_text(size = 8, colour = "black"),
      axis.title         = element_text(size = 10, colour = "black", face = "plain"),
      plot.margin        = unit(c(1, 1, 1, 1), "cm")
    ) +
    ggplot2::coord_cartesian(expand = T)

  return(p_all)
}


#' @title Plot –log10(FDR) vs. Gene Total Count Fold-Ratio (All Gene Pairs)
#' @description
#'       For a given `coordObj`, extract the Lee’s FDR matrix from
#'       `coordObj@grid[[grid_name]][[lee_stats_layer]]$FDR`.
#'       Summarize each gene’s total count across all grid cells from
#'       `coordObj@grid[[grid_name]]$counts`.
#'       For every gene pair (i, j) with i < j,
#'       compute –log10(FDR) and the fold change of their total counts
#'       (fold change = max(total_i, total_j) / min(total_i, total_j), ≥ 1).
#'       Finally, plot a scatterplot with –log10(FDR) on the x-axis and fold change on
#'       the y-axis, using a 25% gray background and publication‐quality aesthetics.
#' 
#'
#' @param coordObj        
#' @param grid_name       
#' @param lee_stats_layer 
#' @param title           
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_rect
#'                     element_line element_text coord_cartesian
#' @importFrom dplyr group_by summarise left_join rename rowwise ungroup filter
#' @importFrom grid unit
#' @export
plotLogFDR <- function(coordObj,
                                       grid_name       = NULL,
                                       lee_stats_layer = NULL,
                                       title           = NULL) {
  ## ---------- 0. 选择 grid 子层 ----------
  if (is.null(coordObj@grid) || length(coordObj@grid) == 0L)
    stop("coordObj@grid is empty.")
  if (is.null(grid_name)) {
    sub_layers <- names(coordObj@grid)
    if (length(sub_layers) == 1L) {
      grid_layer_name <- sub_layers[[1]]
    } else {
      stop("Multiple sublayers exist; please supply grid_name.")
    }
  } else {
    grid_layer_name <- as.character(grid_name)
  }
  if (!(grid_layer_name %in% names(coordObj@grid)))
    stop("grid_name not found: ", grid_layer_name)
  grid_layer <- coordObj@grid[[grid_layer_name]]

  # ---------- 1. ----------
  if (is.null(lee_stats_layer)) {
    cand <- grep("^LeeStats_", names(grid_layer), value = TRUE)
    if (length(cand) == 1L)       lee_stats_layer <- cand
    else if (length(cand) == 0L)  stop("No LeeStats layer found.")
    else                          stop("Multiple LeeStats layers: ",
                                       paste(cand, collapse = ", "),
                                       "; please specify lee_stats_layer.")
  }
  if (!(lee_stats_layer %in% names(grid_layer)))
    stop("lee_stats_layer not found: ", lee_stats_layer)

  LeeStats <- grid_layer[[lee_stats_layer]]
  if (is.null(LeeStats$FDR))
    stop("FDR matrix missing in the specified LeeStats layer. ",
         "Please ensure addLeeStats_fromGrid(…, keep_FDR = TRUE) 已执行。")
  FDRmat <- LeeStats$FDR
  if (!is.matrix(FDRmat))
    stop("LeeStats$FDR must be a matrix.")

  # ---------- 2. ----------
  counts_df <- grid_layer$counts
  if (!all(c("gene", "count") %in% colnames(counts_df)))
    stop("grid_layer$counts 需含 'gene' 与 'count' 列。")
  gene_totals <- counts_df |>
    dplyr::group_by(gene) |>
    dplyr::summarise(total_count = sum(count, na.rm = TRUE),
                     .groups = "drop")

  # ---------- 3. ----------
  genes_F <- rownames(FDRmat)
  genes_T <- gene_totals$gene
  genes_common <- intersect(genes_F, genes_T)
  if (length(genes_common) < 2L)
    stop("Less than two overlapping genes.")

  genes_sorted <- sort(genes_common)
  comb_df <- tidyr::crossing(geneA = genes_sorted, geneB = genes_sorted) |>
    dplyr::filter(geneA < geneB)

  # ---------- 4. ----------
  comb_df <- comb_df |>
    dplyr::left_join(gene_totals, by = c("geneA" = "gene")) |>
    dplyr::rename(totalA = total_count) |>
    dplyr::left_join(gene_totals, by = c("geneB" = "gene")) |>
    dplyr::rename(totalB = total_count) |>
    dplyr::filter(!is.na(totalA), !is.na(totalB)) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      FDR       = FDRmat[geneA, geneB],
      logFDR    = ifelse(is.na(FDR) | FDR <= 0,
                         NA_real_,
                         -log10(pmax(FDR, 1e-300))),
      FoldRatio = ifelse(totalA >= totalB,
                         totalA / totalB,
                         totalB / totalA)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(logFDR), is.finite(FoldRatio))

  # ---------- 5. ----------
  if (is.null(title))
    title <- "–log10(FDR) vs. Fold Change (All Gene Pairs)"

  # ---------- 6. ----------
  library(ggplot2)
  library(grid)

  p_all <- ggplot2::ggplot(comb_df, aes(x = logFDR, y = FoldRatio)) +
      ggplot2::geom_point(
          shape  = 21,
          size   = 1.2,
          fill   = "white",
          color  = "black",
          stroke = 0.2
      ) +
        scale_x_continuous(expand = expansion(mult = 0.05)) +   # 横坐标留 2%
         scale_y_continuous(expand = expansion(mult = 0.05)) +   # 纵坐标留 2%
    ggplot2::labs(
      title = title,
      x     = expression("\u2212log"[10]*"(FDR)"),
      y     = "Fold Change (≥1)"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      panel.background   = element_rect(fill = "#c0c0c0", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "#c0c0c0", size = 0.2),
      panel.border       = element_rect(colour = "black", fill = NA, size = 0.5),
      axis.line          = element_line(colour = "black", size = 0.3),
      axis.ticks         = element_line(colour = "black", size = 0.3),
      axis.text          = element_text(size = 8, colour = "black"),
      axis.title         = element_text(size = 10, colour = "black", face = "plain"),
      plot.margin        = unit(c(1, 1, 1, 1), "cm")
    ) +
    ggplot2::coord_cartesian(expand = T)

  return(p_all)
}
