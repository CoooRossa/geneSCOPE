#' @title Plot Density Heatmap with Segmentation Overlay (Full Grid)
#' @description
#'   Based on `coordObj@grid[[grid_name]]`'s `grid_info` and `densityDF`:
#'   1. Draw a complete continuous density heatmap so that **all** grids (including those with no
#'      expression values) are rendered. Grids that are absent in `densityDF` or whose expression is NA
#'      are treated as zero and filled with the "zero‑expression" colour, ensuring the heatmap covers
#'      the entire frame with no gaps.
#'   2. Overlay the segmentation polygons (alpha adjustable).
#'   3. Re‑use the visual style of `plotCoordObjSpatial()`: black outer border, bottom/left axes only,
#'      equal aspect ratio, and whitespace padding around the panel (`plot.margin`).
#'
#' @param coordObj      A `coordObj` object that must contain:
#'                        `coordObj@grid[[grid_name]]$grid_info` (columns: `grid_id`, `xmin`, `xmax`, `ymin`, `ymax`)
#'                        `coordObj@grid[[grid_name]]$densityDF`  (rownames are `grid_id`)
#'                        `coordObj@coord$segmentation`           (`x`, `y`, `cell`)
#' @param grid_name     Character, e.g. "grid_lenGrid50".
#' @param density1_name Character. Column name of the first density layer in `densityDF`.
#' @param density2_name Character (optional). Column name of the second density layer in `densityDF`.
#' @param palette1      Colour for high values in the first density. Default `"#fc3d5d"`.
#' @param palette2      Colour for high values in the second density. Default `"#4753f8"`.
#' @param alpha1        Alpha for the first layer (0–1). Default `0.5`.
#' @param alpha2        Alpha for the second layer (0–1). Default `0.5`.
#' @param alpha_seg     Alpha of the segmentation overlay. Default `0.2`.
#' @param grid_gap      Spacing of grid lines. Default `100`.
#' @param bar_len       Length of the scale bar. Default `400`.
#' @param bar_offset    Vertical offset of the scale bar as a fraction of the y‑range. Default `0.01`.
#' @param arrow_pt      Arrow‑head size of the scale bar. Default `4`.
#' @return              A `ggplot` object.
#' @import ggplot2 ggforce viridis ggnewscale
#' @export

#' @title Plot Density Heatmap with Segmentation Overlay (Full Grid)
#' @description
#'   Based on `coordObj@grid[[grid_name]]`'s `grid_info` and `densityDF`:
#'   1. Draw a complete continuous density heatmap so that **all** grids (including those with no
#'      expression values) are rendered. Grids that are absent in `densityDF` or whose expression is NA
#'      are treated as zero and filled with the "zero‑expression" colour, ensuring the heatmap covers
#'      the entire frame with no gaps.
#'   2. Overlay the segmentation polygons (alpha adjustable).
#'   3. Re‑use the visual style of `plotCoordObjSpatial()`: black outer border, bottom/left axes only,
#'      equal aspect ratio, and whitespace padding around the panel (`plot.margin`).
#'
#' @param coordObj      A `coordObj` object that must contain:
#'                        `coordObj@grid[[grid_name]]$grid_info` (columns: `grid_id`, `xmin`, `xmax`, `ymin`, `ymax`)
#'                        `coordObj@grid[[grid_name]]$densityDF`  (rownames are `grid_id`)
#'                        `coordObj@coord$segmentation`           (`x`, `y`, `cell`)
#' @param grid_name     Character, e.g. "grid_lenGrid50".
#' @param density1_name Character. Column name of the first density layer in `densityDF`.
#' @param density2_name Character (optional). Column name of the second density layer in `densityDF`.
#' @param palette1      Colour for high values in the first density. Default `"#fc3d5d"`.
#' @param palette2      Colour for high values in the second density. Default `"#4753f8"`.
#' @param alpha1        Alpha for the first layer (0–1). Default `0.5`.
#' @param alpha2        Alpha for the second layer (0–1). Default `0.5`.
#' @param seg_type      Character string: one of "cell", "nucleus", or "both".
#' @param colour_cell   Line colour for cell segmentation polygons. Default "black".
#' @param colour_nucleus Line colour for nucleus segmentation polygons. Default `"#3182bd"`.
#' @param alpha_seg     Alpha of the segmentation overlay. Default `0.2`.
#' @param grid_gap      Spacing of grid lines. Default `100`.
#' @param scale_text_size Font size for scale‑bar text. Default `2.4`.
#' @param bar_len       Length of the scale bar. Default `400`.
#' @param bar_offset    Vertical offset of the scale bar as a fraction of the y‑range. Default `0.01`.
#' @param arrow_pt      Arrow‑head size of the scale bar. Default `4`.
#' @param scale_legend_colour Colour of the scale‑bar and legend text. Default "black".
#' @param max.cutoff1   Fraction (0–1) of the first density's maximum used for clipping. Default `1`.
#' @param max.cutoff2   Fraction (0–1) of the second density's maximum used for clipping. Default `1`.
#' @param legend_digits Number of decimal places in legend labels. Default `1`.
#' @return              A `ggplot` object.
#' @import ggplot2 ggforce viridis ggnewscale
#' @export

plotDensity <- function(coordObj,
                               grid_name,
                               density1_name,
                               density2_name = NULL,
                               palette1       = "#fc3d5d",
                               palette2       = "#4753f8",
                               alpha1         = 0.5,
                               alpha2         = 0.5,
                               seg_type       = c("cell", "nucleus", "both"),
                               colour_cell    = "black",
                               colour_nucleus = "#3182bd",
                               alpha_seg      = 0.2, 
                               grid_gap    = 100, 
                               scale_text_size = 2.4,
                               bar_len     = 400, 
                               bar_offset  = 0.01, 
                               arrow_pt    = 4,
                               scale_legend_colour = "black",
                               max.cutoff1 = 1,
                               max.cutoff2 = 1,
                               legend_digits = 1) {
  # ---- 0. Data validation ----
  seg_type <- match.arg(seg_type)
  accuracy_val <- 1 / (10 ^ legend_digits)
  if (is.null(coordObj@grid) || !(grid_name %in% names(coordObj@grid))) {
    stop("coordObj@grid does not contain sub‑layer:", grid_name)
  }
  layer      <- coordObj@grid[[grid_name]]
  grid_info  <- layer$grid_info
  densityDF  <- layer$densityDF
  if (is.null(grid_info) || is.null(densityDF)) {
    stop("layer lacks grid_info or densityDF")
  }
  if (!all(c("grid_id", "xmin", "xmax", "ymin", "ymax") %in% colnames(grid_info))) {
    stop("grid_info is missing required columns")
  }
  if (!(density1_name %in% colnames(densityDF))) {
    stop("densityDF does not contain column:", density1_name)
  }
  if (!is.null(density2_name) && !(density2_name %in% colnames(densityDF))) {
    stop("densityDF does not contain column:", density2_name)
  }

  # ---- 1. Merge & fill zeros ----
  ids <- grid_info$grid_id
  d1  <- data.frame(grid_id = rownames(densityDF), d = densityDF[[density1_name]])
  heat1 <- merge(grid_info, d1, by = "grid_id", all.x = TRUE)
  heat1$d[is.na(heat1$d)] <- 0
  dmax1 <- max(heat1$d, na.rm = TRUE)
  cut1  <- dmax1 * max.cutoff1
  heat1$d <- pmin(heat1$d, cut1)

  if (!is.null(density2_name)) {
    d2   <- data.frame(grid_id = rownames(densityDF), d = densityDF[[density2_name]])
    heat2 <- merge(grid_info, d2, by = "grid_id", all.x = TRUE)
    heat2$d[is.na(heat2$d)] <- 0
    dmax2 <- max(heat2$d, na.rm = TRUE)
    cut2  <- dmax2 * max.cutoff2
    heat2$d <- pmin(heat2$d, cut2)
  } else heat2 <- NULL

  # ---- 2. Continuous heatmap (geom_tile) ----
  tile_df <- function(df) transform(df,
                                    x = (xmin + xmax)/2,
                                    y = (ymin + ymax)/2,
                                    w = xmax - xmin,
                                    h = ymax - ymin)
  p <- ggplot()
  p <- p + geom_tile(data = tile_df(heat1),
                     aes(x = x, y = y, width = w, height = h, fill = d),
                     colour = NA, alpha = alpha1)
  p <- p + scale_fill_gradient(
    name     = density1_name,
    low      = "transparent",
    high     = palette1,
    na.value = "transparent",
    limits   = c(0, cut1),
    oob      = scales::squish,
    labels   = scales::number_format(accuracy = accuracy_val),
    guide = guide_colourbar(
      title.position = "left",
      title.hjust    = 0.5,
      barwidth       = unit(3, "cm"),   # shortened legend
      barheight      = unit(0.2, "cm"),
      label.theme    = element_text(angle = 90, size = 8),
      order          = 1
    )
  )

  if (!is.null(heat2)) {
    library(ggnewscale)
    p <- p + ggnewscale::new_scale_fill() +
      geom_tile(data = tile_df(heat2),
                aes(x = x, y = y, width = w, height = h, fill = d),
                colour = NA, alpha = alpha2) +
      scale_fill_gradient(
        name     = density2_name,
        low      = "transparent",
        high     = palette2,
        na.value = "transparent",
        limits   = c(0, cut2),
        oob      = scales::squish,
        labels   = scales::number_format(accuracy = accuracy_val)
      )
    p <- p + guides(fill = guide_colourbar(title.position = "left",
                                           title.hjust    = 0.5,
                                           barwidth       = unit(3, "cm"),
                                           barheight      = unit(0.2, "cm"),
                                           label.theme    = element_text(angle = 90, size = 8),
                                           order = 2))
  }

  # ---- 3. Segmentation overlay -----------------------------------------
  seg_layers <- switch(seg_type,
                       cell     = "segmentation_cell",
                       nucleus  = "segmentation_nucleus",
                       both     = c("segmentation_cell", "segmentation_nucleus"))
  seg_layers <- seg_layers[seg_layers %in% names(coordObj@coord)]
  if (length(seg_layers) == 0) {
    warning("No segmentation vertices found for seg_type = ", seg_type)
  } else {
    seg_dt <- data.table::rbindlist(coordObj@coord[seg_layers], use.names = TRUE, fill = TRUE)
    stopifnot(all(c("x", "y", "cell") %in% names(seg_dt)))
    seg_dt$cell <- as.character(seg_dt$cell)

    if (seg_type == "both") {
      seg_dt[, type := ifelse(cell %in% coordObj@coord$segmentation_cell$cell, "cell", "nucleus")]
      p <- p + ggforce::geom_shape(data = seg_dt[type == "cell"],
                                   aes(x = x, y = y, group = cell),
                                   fill = NA, colour = colour_cell,
                                   linewidth = 0.05, alpha = alpha_seg) +
               ggforce::geom_shape(data = seg_dt[type == "nucleus"],
                                   aes(x = x, y = y, group = cell),
                                   fill = NA, colour = colour_nucleus,
                                   linewidth = 0.05, alpha = alpha_seg)
    } else {
      col_use <- if (seg_type == "cell") colour_cell else colour_nucleus
      p <- p + ggforce::geom_shape(data = seg_dt,
                                   aes(x = x, y = y, group = cell),
                                   fill = NA, colour = col_use,
                                   linewidth = 0.05, alpha = alpha_seg)
    }
  }

  # ---- 4. Theme & border ----
  x_rng <- range(grid_info$xmin, grid_info$xmax)
  y_rng <- range(grid_info$ymin, grid_info$ymax)

  dx <- diff(x_rng)
  dy <- diff(y_rng)
  if (dx < dy) {
    pad <- (dy - dx) / 2
    x_rng <- x_rng + c(-pad, pad)
  } else {
    pad <- (dx - dy) / 2
    y_rng <- y_rng + c(-pad, pad)
  }

  # Grid line breaks
  grid_x <- seq(x_rng[1], x_rng[2], by = grid_gap)
  grid_y <- seq(y_rng[1], y_rng[2], by = grid_gap)

  p <- p +
    scale_x_continuous(limits = x_rng, expand = c(0, 0)) +
    scale_y_continuous(limits = y_rng, expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid      = element_blank(),
      panel.border    = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      legend.position = "bottom", 
      legend.direction = "horizontal", 
      legend.box       = "horizontal", 
      legend.key.size = unit(0.5, "cm"),
      legend.text     = element_text(size = 9, angle = 90),
      legend.title    = element_text(size = 8, hjust = 0, vjust = 1),
      plot.title      = element_text(size = 10, hjust = 0.5),
      plot.margin     = margin(1.2, 1, 1.5, 1, "cm"),
      axis.line.x.bottom = element_line(colour = "black"),
      axis.line.y.left   = element_line(colour = "black"),
      axis.line.x.top    = element_blank(),
      axis.line.y.right  = element_blank(),
      axis.title = element_blank()
    ) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))

  # ---- 5. Scale bar -----------------------------------------------------------

  p <- p +
    coord_cartesian(clip = "off") +
    theme(plot.margin = margin(20, 40, 40, 40, "pt"))

  x0    <- x_rng[1] + 0.001 * diff(x_rng)          # 3% padding from left
  y_bar <- y_rng[1] + bar_offset * diff(y_rng)    # inside plot area
  p <- p +
    annotate("segment",
             x = x0, xend = x0 + bar_len,
             y = y_bar,
             yend = y_bar,
             arrow = arrow(length = unit(arrow_pt, "pt"), ends = "both"),
              colour = scale_legend_colour,
             linewidth = 0.4) +
    annotate("text",
             x = x0 + bar_len / 2,
             y = y_bar + 0.025 * diff(y_rng),
             label = paste0(bar_len, " \u00B5m"),
             colour = scale_legend_colour,
             vjust = 1, size = scale_text_size)

  return(p)
}
