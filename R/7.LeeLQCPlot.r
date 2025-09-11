#' @title Plot Distribution of Lee's L Values
#' @description
#'   Extract the Lee's L matrix from a specified grid sublayer and LeeStats layer within `scope_obj`.
#'   Bin the off-diagonal elements of the Lee's L matrix into intervals and plot a histogram to
#'   visualize the distribution of Lee's L values. Optionally, take the absolute value of Lee's L
#'   before plotting. The plot adheres to publication-quality aesthetic standards.
#'
#' @param scope_obj         An object in which Lee's L has already been computed and stored under
#'                         `scope_obj@grid[[grid_name]][[lee_stats_layer]]$L`.
#' @param grid_name        Character (optional). The name of the grid sublayer to use (e.g., `"grid50"`).
#'                         If NULL and there is only one sublayer under `scope_obj@grid`, that layer is used.
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
#' @param title           Character (optional). The plot title. Defaults to `"Lee's L Distribution"` or
#'                         `"Absolute Lee's L Distribution"` if `use_abs = TRUE`.
#' @param use_abs         Logical. Whether to take the absolute value of Lee's L before plotting.
#'                         Defaults to FALSE (plot original Lee's L).
#'
#' @return A `ggplot2` object displaying the histogram of Lee's L values (or their absolute values).
#'
#' @importFrom ggplot2 ggplot aes geom_histogram labs theme_bw theme element_rect element_line element_text coord_cartesian expansion
#' @importFrom grid unit
#' @export
plotLDistribution <- function(scope_obj,
                              grid_name = NULL,
                              lee_stats_layer = NULL,
                              bins = 30,
                              xlim = NULL,
                              title = NULL,
                              use_abs = FALSE) {
  ## ---- 1. Get Lee's L (helper will auto-match layer names) ---------------------------
  Lmat <- .getLeeMatrix(scope_obj,
    grid_name = grid_name,
    lee_layer = lee_stats_layer
  )

  vals <- as.numeric(Lmat[upper.tri(Lmat)])
  if (use_abs) vals <- abs(vals)
  df <- data.frame(LeeL = vals)

  if (is.null(title)) {
    title <- if (use_abs) "Absolute Lee's L Distribution" else "Lee's L Distribution"
  }

  library(ggplot2)
  library(grid)

  p <- ggplot(df, aes(x = LeeL)) +
    {
      if (length(bins) == 1L) {
        geom_histogram(
          bins = bins,
          fill = "white",
          colour = "black",
          size = .3
        )
      } else {
        geom_histogram(
          breaks = bins,
          fill = "white",
          colour = "black",
          size = .3
        )
      }
    } +
    labs(
      title = title,
      x = if (use_abs) expression("|Lee's L|") else expression("Lee's L"),
      y = "Frequency"
    ) +
    theme_bw(base_size = 8) +
    theme(
      panel.title        = element_text(size = 10, colour = "black"),
      panel.background   = element_rect(fill = "#c0c0c0", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "#c0c0c0", size = .2),
      panel.border       = element_rect(colour = "black", fill = NA, size = .5),
      axis.line          = element_line(colour = "black", size = .3),
      axis.ticks         = element_line(colour = "black", size = .3),
      axis.text          = element_text(size = 8, colour = "black"),
      axis.title         = element_text(size = 9, colour = "black"),
      plot.margin        = unit(c(1, 1, 1, 1), "cm")
    ) +
    coord_cartesian(expand = FALSE)

  if (!is.null(xlim)) {
    stopifnot(is.numeric(xlim) && length(xlim) == 2L)
    p <- p + coord_cartesian(xlim = xlim, expand = FALSE)
  }
  p
}

#' @title Plot Lee's L vs. Gene Total Count Fold Change (All Gene Pairs)
#' @description
#'   For a given `scope_obj`, extract the Lee's L matrix from
#'   `scope_obj@grid[[grid_name]][[lee_stats_layer]]$L`.
#'   Summarize each gene's total count across all grid cells from
#'   `scope_obj@grid[[grid_name]]$counts`. For every gene pair (i, j) with i < j,
#'   compute Lee's L and the fold change of their total counts
#'   (fold change = max(total_i, total_j) / min(total_i, total_j), ≥ 1).
#'   Finally, plot a scatterplot with Lee's L on the x-axis and fold change on
#'   the y-axis, using a 25% gray background and publication-quality aesthetics.
#'
#' @param scope_obj        An object containing
#'                        `scope_obj@grid[[grid_name]]$counts` and
#'                        `scope_obj@grid[[grid_name]][[lee_stats_layer]]$L`.
#' @param grid_name       Character (optional). The name of the grid sublayer to use.
#'                        If NULL and `scope_obj@grid` has exactly one sublayer,
#'                        that sublayer is used automatically. Otherwise, it must be provided.
#' @param lee_stats_layer Character (optional). The name of the LeeStats layer (e.g., `"LeeStats_Xz"`).
#'                        If NULL, the function searches for exactly one layer under
#'                        `scope_obj@grid[[grid_name]]` whose name starts with `"LeeStats_"`.
#'                        If multiple matches are found, an error is thrown and the user must
#'                        supply `lee_stats_layer` explicitly.
#' @param title           Character (optional). Plot title. Defaults to
#'                        `"Lee's L vs. Fold Change (All Gene Pairs)"` if NULL.
#'
#' @return A `ggplot2` object showing a scatterplot of Lee's L vs. fold change
#'         for all gene pairs (i < j).
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_bw theme element_rect element_line element_text coord_cartesian scale_x_continuous scale_y_continuous expansion
#' @importFrom dplyr group_by summarise filter left_join rename rowwise ungroup
#' @importFrom tidyr crossing
#' @importFrom grid unit
#' @export
plotLScatter <- function(scope_obj,
                         grid_name = NULL,
                         lee_stats_layer = NULL,
                         title = NULL) {
  ## ---- 1. Grid layer & counts ----------------------------------------------
  g_layer <- .selectGridLayer(scope_obj, grid_name)
  if (is.null(grid_name)) {
    grid_name <- names(scope_obj@grid)[
      vapply(scope_obj@grid, identical, logical(1), g_layer)
    ]
  }

  counts_df <- g_layer$counts
  if (!all(c("gene", "count") %in% names(counts_df))) {
    stop("counts must contain columns 'gene' and 'count'.")
  }

  ## ---- 2. Lee's L -------------------------------------------------------
  Lmat <- .getLeeMatrix(scope_obj,
    grid_name = grid_name,
    lee_layer = lee_stats_layer
  )

  ## ---- 3. Total & combine gene pairs ---------------------------------------------
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)

  gene_totals <- counts_df |>
    group_by(gene) |>
    summarise(total = sum(count), .groups = "drop")

  genes_common <- intersect(rownames(Lmat), gene_totals$gene)
  if (length(genes_common) < 2L) {
    stop("Need at least two overlapping genes.")
  }

  comb_df <- tidyr::crossing(
    geneA = genes_common,
    geneB = genes_common
  ) |>
    filter(geneA < geneB) |>
    left_join(gene_totals, by = c("geneA" = "gene")) |>
    rename(totalA = total) |>
    left_join(gene_totals, by = c("geneB" = "gene")) |>
    rename(totalB = total) |>
    rowwise() |>
    mutate(
      LeeL = Lmat[geneA, geneB],
      FoldRatio = max(totalA, totalB) /
        pmax(1, pmin(totalA, totalB))
    ) |>
    ungroup() |>
    filter(is.finite(LeeL), is.finite(FoldRatio))

  if (is.null(title)) {
    title <- "Lee's L vs. Fold Change (All Gene Pairs)"
  }

  p <- ggplot(comb_df, aes(x = LeeL, y = FoldRatio)) +
    geom_point(
      shape = 21, size = 1.2,
      fill = "white", colour = "black", stroke = .2
    ) +
    scale_x_continuous(expand = expansion(mult = .05)) +
    scale_y_continuous(expand = expansion(mult = .05)) +
    labs(
      title = title,
      x = expression("Lee's L"),
      y = "Fold Change (≥1)"
    ) +
    theme_bw(base_size = 8) +
    theme(
      panel.title        = element_text(size = 10, colour = "black"),
      panel.background   = element_rect(fill = "#c0c0c0", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(colour = "#c0c0c0", size = .2),
      panel.border       = element_rect(colour = "black", fill = NA, size = .5),
      axis.line          = element_line(colour = "black", size = .3),
      axis.ticks         = element_line(colour = "black", size = .3),
      axis.text          = element_text(size = 8, colour = "black"),
      axis.title         = element_text(size = 9, colour = "black"),
      plot.margin        = unit(c(1, 1, 1, 1), "cm")
    ) +
    coord_cartesian(expand = TRUE)

  p
}
