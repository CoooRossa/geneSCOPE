#' @title Plot Spatial Segmentation from coordObj with Optional Group Coloring
#' @description
#'   Draw tissue spatial contours based on the segmentation data stored in a `coordObj`.
#'   You can optionally pass a `data.frame` (with columns **cell.id** and **group**) to
#'   colour each cell polygon; when omitted, all polygons are rendered in grey.
#'
#' @param coordObj   A `coordObj` object that already contains
#'                   `coordObj@coord$segmentation` with columns `x`, `y`, and `cell`.
#' @param colors     Named character vector. Colors for each group; if `NULL`, all cells are grey.
#'                   If `colors` is provided, it should have names corresponding to the
#'                   `group` column in `coordObj@coord$segmentation`. If no `group` column exists,
#'                   all cells are rendered in grey.
#' @param grid_gap   Numeric. Grid line spacing (default `100`).
#' @param bar_len    Numeric. Length of the scale bar (default `400`, in coordinate units).
#' @param bar_offset Numeric. Vertical offset of the scale bar relative to the y‑range
#'                   (default `0.01`).
#' @param arrow_pt   Numeric. Arrow‑head size of the scale bar (pt; default `4`).
#' @param title      Character. Plot title; default is "Spatial Segmentation" when `NULL`.
#' @return           A `ggplot2` object containing only the plot.
#' @import ggplot2
#' @import ggforce
#' @import scales
#' @export
plotCoordObjSpatial <- function(coordObj,
                                colors      = NULL,
                                grid_gap    = 100,
                                bar_len     = 400,
                                bar_offset  = 0.01,
                                arrow_pt    = 4,
                                title       = NULL,
                                alpha       = 0.5) {
  # ---- 1. Check coordObj@coord$segmentation ----
  if (is.null(coordObj@coord) || is.null(coordObj@coord$segmentation)) {
    stop("coordObj@coord$segmentation does not exist. Please ensure segmentation data is present in coordObj.")
  }
  seg_df <- coordObj@coord$segmentation
  if (!all(c("x", "y", "cell") %in% colnames(seg_df))) {
    stop("coordObj@coord$segmentation must contain columns: x, y, cell.")
  }
  seg_df$cell <- as.character(seg_df$cell)
  seg_df$group <- factor("all")
  colors <- setNames(colors, "all")

  # ---- 3. Set title ----
  if (is.null(title)) title <- "Spatial Segmentation"

  # ---- 4. Draw base polygon map ----
  p <- ggplot() +
    ggforce::geom_shape(
      data      = seg_df,
      aes(x = x, y = y, group = cell, fill = group, colour = group),
      linewidth = 0.15,
      alpha     = alpha,
      radius    = 0
    ) +
    scale_fill_manual(values = colors, na.value = "grey") +
    scale_colour_manual(values = colors, guide = "none", na.value = "grey") +
    coord_fixed() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid      = element_blank(),
      panel.border    = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      legend.position = "bottom",
      legend.key.size = unit(0.3, "cm"),
      legend.title    = element_blank(),
      legend.text     = element_text(size = 9),
      plot.title      = element_text(size = 10, hjust = 0.5),
      plot.margin     = margin(1, 1, 1, 1, "cm")
    ) +
    labs(title = title)

  # ---- 5. Limit plot range and add border ----
  x_rng <- range(seg_df$x, na.rm = TRUE)
  y_rng <- range(seg_df$y, na.rm = TRUE)

  dx <- diff(x_rng)
  dy <- diff(y_rng)
  if (dx < dy) {
    pad <- (dy - dx) / 2
    x_rng <- x_rng + c(-pad, pad)
  } else {
    pad <- (dx - dy) / 2
    y_rng <- y_rng + c(-pad, pad)
  }
  x_mid  <- mean(x_rng)
  y_mid  <- mean(y_rng)
  x_half <- diff(x_rng)/2 * 0.9
  y_half <- diff(y_rng)/2 * 0.9
  x_lim  <- c(x_mid - x_half, x_mid + x_half)
  y_lim  <- c(y_mid - y_half, y_mid + y_half)

  p <- p +
    scale_x_continuous(limits = x_rng, expand = c(0, 0)) +
    scale_y_continuous(limits = y_rng, expand = c(0, 0)) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))

#   # ---- 6. Add grid lines; keep only left / bottom axes ----
  grid_x <- seq(x_rng[1], x_rng[2], by = grid_gap)
  grid_y <- seq(y_rng[1], y_rng[2], by = grid_gap)
  p <- p +
    theme(
      axis.line.x.bottom = element_line(colour = "black"),
      axis.line.y.left   = element_line(colour = "black"),
      axis.line.x.top    = element_blank(),
      axis.line.y.right  = element_blank(),
      axis.title         = element_blank(),
      axis.text          = element_blank(),
      axis.ticks.length  = unit(2, "pt")
    )

  # ---- 7. Add scale bar ----
  p <- p +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(20, 40, 40, 40, "pt")) +
    scale_x_continuous(
      limits = x_rng, expand = c(0, 0),
      breaks  = grid_x, oob = scales::oob_keep
    ) +
    scale_y_continuous(
      limits = y_rng, expand = c(0, 0),
      breaks  = grid_y, oob = scales::oob_keep
    )

  x0    <- x_rng[1]
  y_off <- bar_offset * diff(x_rng)
  p <- p +
    annotate(
      "segment",
      x    = x0, xend = x0 + bar_len,
      y    = y_rng[1] - y_off,
      yend = y_rng[1] - y_off,
      arrow= arrow(length = unit(arrow_pt, "pt"), ends = "both"),
      linewidth = 0.4
    ) +
    annotate(
      "text",
      x = x0 + bar_len / 2,
      y = y_rng[1] - y_off * 1.6,
      label = paste0(bar_len, " \u00B5m"),
      vjust = 1, size = 5.2
    )
  p <- p + theme(legend.position = "none")
  return(p)
}
