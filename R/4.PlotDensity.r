#' 4.PlotDensity.r (2025-01-06 rev)
#' @title plotDensity - visualise grid-level density stored in @density
#'
#' @param scope_obj  A \code{scope_object}.
#' @param grid_name Name of the grid layer to plot (e.g. "grid50").
#' @param density1_name,density2_name Column names to plot.
#' @param palette1,palette2 Color palettes for density1 and density2.
#' @param alpha1,alpha2 Alpha transparency for density layers.
#' @param seg_type Segmentation overlay type: "cell", "nucleus", or "both".
#' @param colour_cell,colour_nucleus Colors for cell and nucleus segmentation.
#' @param alpha_seg Alpha transparency for segmentation overlay.
#' @param grid_gap Spacing for background grid lines.
#' @param scale_text_size Font size for scale bar text.
#' @param bar_len Length of scale bar in micrometers.
#' @param bar_offset Vertical offset of scale bar as fraction of y-range.
#' @param arrow_pt Arrow point size for scale bar.
#' @param scale_legend_colour Color of scale bar and text.
#' @param max.cutoff1,max.cutoff2 Maximum cutoff fractions for density values.
#' @param legend_digits Number of decimal places in legend.
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot geom_tile aes scale_fill_gradient guide_colorbar geom_vline geom_hline scale_x_continuous scale_y_continuous coord_fixed theme_minimal theme element_blank element_rect element_text element_line unit margin annotate arrow
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggforce geom_shape
#' @importFrom scales alpha number_format squish
#' @importFrom data.table rbindlist
#' @export
plotDensity <- function(scope_obj,
                        grid_name,
                        density1_name,
                        density2_name = NULL,
                        palette1 = "#fc3d5d",
                        palette2 = "#4753f8",
                        alpha1 = 0.5,
                        alpha2 = 0.5,
                        overlay_image = FALSE,
                        image_path = NULL,
                        image_alpha = 0.6,
                        image_choice = c("auto", "hires", "lowres"),
                        seg_type = c("cell", "nucleus", "both"),
                        colour_cell = "black",
                        colour_nucleus = "#3182bd",
                        alpha_seg = 0.2,
                        grid_gap = 100,
                        scale_text_size = 2.4,
                        bar_len = 400,
                        bar_offset = 0.01,
                        arrow_pt = 4,
                        scale_legend_colour = "black",
                        max.cutoff1 = 1,
                        max.cutoff2 = 1,
                        legend_digits = 1) {
    seg_type <- match.arg(seg_type)
    image_choice <- match.arg(image_choice)
    accuracy_val <- 1 / (10^legend_digits)

    ## ------------------------------------------------------------------ 1
    ## validate grid layer & retrieve geometry
    g_layer <- .selectGridLayer(scope_obj, grid_name)
    grid_layer_name <- names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    grid_info <- g_layer$grid_info
    if (is.null(grid_info) ||
        !all(c("grid_id", "xmin", "xmax", "ymin", "ymax") %in% names(grid_info))) {
        stop("grid_info is missing required columns for layer '", grid_layer_name, "'.")
    }

    ## ------------------------------------------------------------------ 2
    ## locate density data frame (new slot first, old slot as fallback)
    densityDF <- scope_obj@density[[grid_layer_name]]
    if (is.null(densityDF)) {
        densityDF <- g_layer$densityDF
    } # legacy
    if (is.null(densityDF)) {
        stop("No density table found for grid '", grid_layer_name, "'.")
    }

    ## make sure rownames are grid IDs
    if (is.null(rownames(densityDF)) && "grid_id" %in% names(densityDF)) {
        rownames(densityDF) <- densityDF$grid_id
    }

    for (dcol in c(density1_name, density2_name)) {
        if (!is.null(dcol) && !(dcol %in% colnames(densityDF))) {
            stop(
                "Density column '", dcol, "' not found in table for grid '",
                grid_layer_name, "'."
            )
        }
    }

    ## helper to assemble heatmap data.frame
    build_heat <- function(d_col, cutoff_frac) {
        df <- data.frame(
            grid_id = rownames(densityDF),
            d = densityDF[[d_col]],
            stringsAsFactors = FALSE
        )
        heat <- merge(grid_info, df, by = "grid_id", all.x = TRUE)
        heat$d[is.na(heat$d)] <- 0
        maxv <- max(heat$d, na.rm = TRUE)
        heat$cut <- maxv * cutoff_frac
        heat$d <- pmin(heat$d, heat$cut)
        heat
    }

    heat1 <- build_heat(density1_name, max.cutoff1)
    heat2 <- if (!is.null(density2_name)) build_heat(density2_name, max.cutoff2)

    ## ------------------------------------------------------------------ 3
    ## tile geometry (centre & size)
    tile_df <- function(df) {
        transform(df,
            x = (xmin + xmax) / 2,
            y = (ymin + ymax) / 2,
            w = xmax - xmin,
            h = ymax - ymin
        )
    }

    library(ggplot2)
    p <- ggplot()

    ## ------------------------------------------------------------------ 3.0
    ## optional background image overlay (auto-detect from grid layer)
    if (isTRUE(overlay_image)) {
        img_info <- g_layer$image_info
        mpp_fullres <- g_layer$microns_per_pixel
        hires_path <- if (!is.null(img_info)) img_info$hires_path else NULL
        lowres_path <- if (!is.null(img_info)) img_info$lowres_path else NULL
        if (!is.null(image_path)) {
            hires_path <- image_path
            lowres_path <- NULL
        }
        sel_path <- switch(image_choice,
            auto = if (!is.null(hires_path)) hires_path else lowres_path,
            hires = hires_path,
            lowres = lowres_path
        )
        if (!is.null(sel_path) && file.exists(sel_path) && is.finite(mpp_fullres)) {
            ext <- tolower(tools::file_ext(sel_path))
            img <- NULL
            if (ext %in% c("png") && requireNamespace("png", quietly = TRUE)) {
                img <- png::readPNG(sel_path)
            } else if (ext %in% c("jpg","jpeg") && requireNamespace("jpeg", quietly = TRUE)) {
                img <- jpeg::readJPEG(sel_path)
            }
            if (!is.null(img)) {
                # Normalize to RGBA array and apply alpha (annotation_raster has no alpha arg)
                if (length(dim(img)) == 2L) {
                    # grayscale -> RGB
                    img_rgb <- array(NA_real_, dim = c(dim(img)[1], dim(img)[2], 3L))
                    img_rgb[,,1] <- img; img_rgb[,,2] <- img; img_rgb[,,3] <- img
                    img <- img_rgb
                }
                if (length(dim(img)) == 3L) {
                    if (dim(img)[3] == 3L) {
                        img_rgba <- array(NA_real_, dim = c(dim(img)[1], dim(img)[2], 4L))
                        img_rgba[,,1:3] <- img
                        img_rgba[,,4]   <- max(0, min(1, image_alpha))
                        img <- img_rgba
                    } else if (dim(img)[3] == 4L) {
                        img[,,4] <- pmin(1, pmax(0, img[,,4] * image_alpha))
                    }
                }

                # image size (width x height) in um
                wpx <- dim(img)[2]; hpx <- dim(img)[1]
                # Determine effective microns-per-pixel depending on image type and scalefactors
                hires_sf <- if (!is.null(img_info)) img_info$tissue_hires_scalef else NA_real_
                lowres_sf <- if (!is.null(img_info)) img_info$tissue_lowres_scalef else NA_real_
                eff_mpp <- mpp_fullres
                # Detect by filename if necessary
                is_hires <- grepl("tissue_hires_image", basename(sel_path), ignore.case = TRUE)
                is_lowres <- grepl("tissue_lowres_image", basename(sel_path), ignore.case = TRUE)
                is_aligned <- grepl("aligned_tissue_image", basename(sel_path), ignore.case = TRUE)
                if (is_hires && is.finite(hires_sf)) eff_mpp <- mpp_fullres / hires_sf
                if (is_lowres && is.finite(lowres_sf)) eff_mpp <- mpp_fullres / lowres_sf
                if (is_aligned) eff_mpp <- mpp_fullres
                # Fallback: if cannot infer, try to match grid width to image width by simple ratio
                if (!is.finite(eff_mpp) || eff_mpp <= 0) eff_mpp <- mpp_fullres
                w_um <- as.numeric(wpx) * as.numeric(eff_mpp)
                h_um <- as.numeric(hpx) * as.numeric(eff_mpp)
                y_origin <- if (!is.null(img_info$y_origin)) img_info$y_origin else "top-left"
                # If grid was not flipped (top-left), flip image vertically by swapping ymin/ymax
                if (identical(y_origin, "top-left")) {
                    p <- p + annotation_raster(img,
                        xmin = 0, xmax = w_um,
                        ymin = h_um, ymax = 0,
                        interpolate = TRUE
                    )
                } else {
                    p <- p + annotation_raster(img,
                        xmin = 0, xmax = w_um,
                        ymin = 0, ymax = h_um,
                        interpolate = TRUE
                    )
                }
            }
        }
    }

    p <- p +
        geom_tile(
            data = tile_df(heat1),
            aes(x = x, y = y, width = w, height = h, fill = d),
            colour = NA, alpha = alpha1
        ) +
        scale_fill_gradient(
            name = density1_name,
            low = "transparent",
            high = palette1,
            limits = c(0, unique(heat1$cut)),
            oob = scales::squish,
            labels = scales::number_format(accuracy = accuracy_val),
            na.value = "transparent",
            guide = guide_colorbar(order = 1)
        )

    ## second density with new fill scale
    if (!is.null(heat2)) {
        library(ggnewscale)
        p <- p + ggnewscale::new_scale_fill() +
            geom_tile(
                data = tile_df(heat2),
                aes(x = x, y = y, width = w, height = h, fill = d),
                colour = NA, alpha = alpha2
            ) +
            scale_fill_gradient(
                name = density2_name,
                low = "transparent",
                high = palette2,
                limits = c(0, unique(heat2$cut)),
                oob = scales::squish,
                labels = scales::number_format(accuracy = accuracy_val),
                na.value = "transparent",
                guide = guide_colorbar(order = 2)
            )
    }

    ## ------------------------------------------------------------------ 4
    ## segmentation overlay
    seg_layers <- switch(seg_type,
        cell    = "segmentation_cell",
        nucleus = "segmentation_nucleus",
        both    = c("segmentation_cell", "segmentation_nucleus")
    )
    seg_layers <- seg_layers[seg_layers %in% names(scope_obj@coord)]

    if (length(seg_layers)) {
        seg_dt <- data.table::rbindlist(scope_obj@coord[seg_layers],
            use.names = TRUE, fill = TRUE
        )
        if (seg_type == "both") {
            is_cell <- seg_dt$cell %in% scope_obj@coord$segmentation_cell$cell
            seg_dt[, segClass := ifelse(is_cell, "cell", "nucleus")]
            p <- p +
                ggforce::geom_shape(
                    data = seg_dt[segClass == "cell"],
                    aes(x = x, y = y, group = cell),
                    colour = scales::alpha(colour_cell, alpha_seg),
                    fill = NA,
                    linewidth = 0.05
                ) +
                ggforce::geom_shape(
                    data = seg_dt[segClass == "nucleus"],
                    aes(x = x, y = y, group = cell),
                    colour = scales::alpha(colour_nucleus, alpha_seg),
                    fill = NA,
                    linewidth = 0.05
                )
        } else {
            ## Single type overlay (stroke transparency baked into colour)
            col_use <- if (seg_type == "cell") colour_cell else colour_nucleus
            p <- p +
                ggforce::geom_shape(
                    data = seg_dt,
                    aes(x = x, y = y, group = cell),
                    colour = scales::alpha(col_use, alpha_seg),
                    fill = NA,
                    linewidth = 0.05
                )
        }
    }

    ## ------------------------------------------------------------------ 5
    ## aesthetics: square, gridlines, scale‑bar
    x_rng <- range(grid_info$xmin, grid_info$xmax)
    y_rng <- range(grid_info$ymin, grid_info$ymax)

    # pad shorter axis
    dx <- diff(x_rng)
    dy <- diff(y_rng)
    if (dx < dy) {
        pad <- (dy - dx) / 2
        x_rng <- x_rng + c(-pad, pad)
    }
    if (dy < dx) {
        pad <- (dx - dy) / 2
        y_rng <- y_rng + c(-pad, pad)
    }

    # gridlines
    grid_x <- seq(x_rng[1], x_rng[2], by = grid_gap)
    grid_y <- seq(y_rng[1], y_rng[2], by = grid_gap)
    p <- p +
        geom_vline(xintercept = grid_x, linewidth = 0.05, colour = "grey80") +
        geom_hline(yintercept = grid_y, linewidth = 0.05, colour = "grey80")

    p <- p +
        scale_x_continuous(limits = x_rng, expand = c(0, 0)) +
        scale_y_continuous(limits = y_rng, expand = c(0, 0)) +
        coord_fixed() +
        theme_minimal(base_size = 9) +
        theme(
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "horizontal",
            legend.key.width = unit(0.4, "cm"), # key 的“宽度”
            legend.key.height = unit(0.15, "cm"), # key 的“高度”
            legend.text = element_text(size = 9, angle = 90),
            legend.title = element_text(size = 8, hjust = 0, vjust = 1),
            plot.title = element_text(size = 10, hjust = 0.5),
            plot.margin = margin(1.2, 1, 1.5, 1, "cm"),
            axis.line.x.bottom = element_line(colour = "black"),
            axis.line.y.left = element_line(colour = "black"),
            axis.line.x.top = element_blank(),
            axis.line.y.right = element_blank(),
            axis.title = element_blank()
        ) +
        theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))

    ## scale‑bar
    x0 <- x_rng[1] + 0.01 * diff(x_rng) # 3% padding from left
    y_bar <- y_rng[1] + bar_offset * diff(y_rng) # inside plot area
    p <- p +
        annotate("segment",
            x = x0, xend = x0 + bar_len,
            y = y_bar,
            yend = y_bar,
            arrow = arrow(length = unit(arrow_pt, "pt"), ends = "both"),
            colour = scale_legend_colour,
            linewidth = 0.4
        ) +
        annotate("text",
            x = x0 + bar_len / 2,
            y = y_bar + 0.025 * diff(y_rng),
            label = paste0(bar_len, " \u00B5m"),
            colour = scale_legend_colour,
            vjust = 1, size = scale_text_size
        )
    return(p)
}

#' @title Mirror Plot of Two Gene Densities Across Grids
#'
#' @description
#'   Creates a mirrored area plot in which the density values of \code{gene1}
#'   are shown above the baseline and those of \code{gene2} are negated and
#'   shown below, producing an intuitive comparison across grids.  The grids
#'   can be ordered by the sum, difference, or individual gene density, and
#'   densities may be min-max rescaled to 0-1 for direct visual comparison.
#'
#' @param scope_obj A \code{scope_object} containing a \code{@grid} slot that
#'                  includes \code{densityDF}.
#' @param layer_name Character. The grid sub-layer to query.
#' @param gene1,gene2 Character. Names of the two genes to compare.
#' @param sort_by  Ordering criterion: \code{"sum"}, \code{"diff"},
#'                 \code{"gene1"}, or \code{"gene2"}.
#' @param rescale  Logical. Rescale densities to 0-1 before plotting.
#'
#' @return Invisibly returns the generated \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_area scale_y_continuous expansion labs theme_bw theme element_rect element_blank element_line element_text
#' @importFrom dplyr filter summarise collect select pull mutate arrange case_when if_else across all_of
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom scales rescale
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' p <- plotDensityMirror(P5.coord, "25um", "GAPDH", "MALAT1",
#'     sort_by = "diff"
#' )
#' plot(p)
#' }
#' @export
plotDensityMirror <- function(scope_obj,
                              layer_name,
                              gene1,
                              gene2,
                              sort_by = c("sum", "diff", "gene1", "gene2"),
                              rescale = TRUE) {
    sort_by <- match.arg(sort_by)
    genes <- c(gene1, gene2)

    # --- 0. Extract layer data ------------------------------------------------
    if (layer_name %in% names(scope_obj@grid)) {
        # Grid-based layer
        g_layer <- .selectGridLayer(scope_obj, layer_name)
        density_df <- scope_obj@density[[layer_name]]
        if (is.null(density_df) && !is.null(g_layer$densityDF)) {
            density_df <- g_layer$densityDF
        }
        if (is.null(density_df)) {
            stop("Density table not found for grid '", layer_name, "'.")
        }
        if (is.null(rownames(density_df)) && "grid_id" %in% colnames(density_df)) {
            rownames(density_df) <- density_df$grid_id
        }
        df <- tibble::rownames_to_column(density_df, "grid_id")[, c("grid_id", genes)]
    } else if (layer_name %in% names(scope_obj@cells)) {
        # Cell-based layer (counts or logCPM)
        mat <- scope_obj@cells[[layer_name]]
        if (methods::is(mat, "dgCMatrix")) {
            mat <- as.matrix(mat)
        }
        if (!all(genes %in% rownames(mat))) {
            stop("Specified gene(s) not found in cell layer '", layer_name, "'.")
        }
        df <- data.frame(
            grid_id = colnames(mat),
            gene1 = mat[gene1, , drop = TRUE],
            gene2 = mat[gene2, , drop = TRUE],
            stringsAsFactors = FALSE
        )
        colnames(df)[2:3] <- genes
    } else {
        stop("Layer '", layer_name, "' not found in scope_obj.")
    }

    # --- 1. Filter & optional rescaling ---------------------------------------
    df <- df[df[[gene1]] != 0 | df[[gene2]] != 0, ]
    if (rescale) {
        df <- df |>
            dplyr::mutate(
                across(all_of(genes), scales::rescale)
            )
    }

    # --- 2. Compute sort key & ordering ---------------------------------------
    df$sort_key <- dplyr::case_when(
        sort_by == "sum" ~ df[[gene1]] + df[[gene2]],
        sort_by == "diff" ~ df[[gene1]] - df[[gene2]],
        sort_by == "gene1" ~ df[[gene1]],
        TRUE ~ df[[gene2]]
    )
    df <- df |>
        dplyr::arrange(sort_key) |>
        dplyr::mutate(
            grid_order = factor(grid_id, levels = grid_id)
        )

    plot_df <- df |>
        dplyr::select(grid_order, all_of(genes)) |>
        tidyr::pivot_longer(
            -grid_order,
            names_to = "gene",
            values_to = "value"
        ) |>
        dplyr::mutate(
            value = if_else(gene == gene2, -value, value),
            gene  = factor(gene, levels = c(gene1, gene2))
        )

    # --- 3. Plot --------------------------------------------------------------
    p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = grid_order, y = value, fill = gene, group = gene)
    ) +
        ggplot2::geom_area(alpha = .65, size = .15, position = "identity") +
        ggplot2::scale_y_continuous(labels = abs, expand = ggplot2::expansion(mult = .02)) +
        ggplot2::labs(
            title = sprintf(
                "%s: %s vs %s",
                layer_name, gene1, gene2
            ),
            x = ifelse(layer_name %in% names(scope_obj@grid), "Grid", "Cell"),
            y = ifelse(rescale, "Rescaled expression", "Expression")
        ) +
        ggplot2::theme_bw(base_size = 8) +
        ggplot2::theme(
            panel.background   = ggplot2::element_rect(fill = "#c0c0c0", colour = NA),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor   = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(colour = "#c0c0c0", size = 0.2),
            panel.border       = ggplot2::element_rect(colour = "black", fill = NA, size = 0.5),
            axis.line.x        = ggplot2::element_blank(),
            axis.ticks.x       = ggplot2::element_blank(),
            axis.text.x        = ggplot2::element_blank(),
            axis.text.y        = ggplot2::element_text(size = 8),
            plot.title         = ggplot2::element_text(hjust = .5, size = 10)
        )

    p
}


#' plotGridBoundary
#'
#' Draw grid cell boundaries for a given grid layer in a coord object and
#' guarantee that the output image keeps the exact spatial aspect ratio
#' of the original coordinates (no stretching or squeezing), following the
#' same padding logic used in `plotDensity()`.
#'
#' @param scope_obj   A coord-style S4 object that contains a `@grid` slot,
#'                   e.g. produced by the STUtility pipeline.
#' @param grid_name  Character. The name of the grid layer to plot, e.g.
#'                   "grid31". Must exist inside `scope_obj@grid`.
#' @param colour     Border colour for grid rectangles. Default "black".
#' @param linewidth  Border line width for rectangles. Default 0.2.
#' @param panel_bg   Background colour of the panel. Default "#C0C0C0".
#' @param base_size  Base font size for `theme_minimal()`. Default 10.
#'
#' @return A `ggplot` object.
#' @importFrom ggplot2 ggplot geom_rect aes coord_fixed theme_minimal theme element_text element_rect element_blank element_line labs unit
#' @export
#'
#' @examples
#' p <- plotGridBoundary(P5.coord, "grid31")
#' plot(p)
#' @export
plotGridBoundary <- function(scope_obj,
                             grid_name,
                             colour = "black",
                             linewidth = 0.2,
                             panel_bg = "#C0C0C0",
                             base_size = 10) {
    ## --- 0. Retrieve grid layer & grid_info ------------------------------
    g_layer <- .selectGridLayer(scope_obj, grid_name)
    grid_info <- g_layer$grid_info
    if (is.null(grid_info) ||
        !all(c("grid_id", "xmin", "xmax", "ymin", "ymax") %in% names(grid_info))) {
        stop("grid_info is missing required columns.")
    }

    ## --- 1. Square padding (same logic as plotDensity) -------------------
    x_rng <- range(grid_info$xmin, grid_info$xmax, na.rm = TRUE)
    y_rng <- range(grid_info$ymin, grid_info$ymax, na.rm = TRUE)
    dx <- diff(x_rng)
    dy <- diff(y_rng)
    if (dx < dy) {
        pad <- (dy - dx) / 2
        x_rng <- x_rng + c(-pad, pad)
    }
    if (dy < dx) {
        pad <- (dx - dy) / 2
        y_rng <- y_rng + c(-pad, pad)
    }

    ## --- 2. Build plot ----------------------------------------------------
    library(ggplot2)
    p <- ggplot(grid_info) +
        geom_rect(
            aes(
                xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax
            ),
            fill = NA, colour = colour, linewidth = linewidth
        ) +
        coord_fixed(
            ratio = diff(y_rng) / diff(x_rng),
            xlim = x_rng, ylim = y_rng,
            expand = FALSE, clip = "off"
        ) +
        theme_minimal(base_size = base_size) +
        theme(
            panel.title        = element_text(size = 10, colour = "black"),
            panel.background   = element_rect(fill = panel_bg, colour = NA),
            panel.grid.major.x = element_blank(),
            panel.grid.minor   = element_blank(),
            panel.grid.major.y = element_line(colour = "#c0c0c0", size = 0.2),
            panel.border       = element_rect(colour = "black", fill = NA, size = 0.5),
            axis.line          = element_line(colour = "black", size = 0.3),
            axis.ticks         = element_line(colour = "black", size = 0.3),
            axis.text          = element_text(size = 8, colour = "black"),
            axis.title         = element_text(size = 9, colour = "black"),
            plot.margin        = unit(c(1, 1, 1, 1), "cm")
        ) +
        labs(
            x = "X", y = "Y",
            title = paste0(
                "Grid Cells: Spatial Boundary Distribution\nGrid Size ",
                gsub(".*Grid", "", grid_name)
            )
        )

    p
}

#' @title Plot Gene Centroid Expression with Segmentation Overlay
#' @description
#'   Plots cell centroids coloured by the normalised expression of two specified genes.
#'   Expression values are extracted from `scope_obj@cells$counts`, min-max scaled
#'   to 0-1 separately for each gene, and mapped to point transparency (`alpha`) so
#'   that higher expression appears more opaque. Centroids of cells expressing
#'   `gene1` are drawn with `palette1`; centroids of cells expressing `gene2` with
#'   `palette2`. Cells that express neither gene are not shown. Optionally overlays
#'   cell/nucleus segmentation polygons.
#'
#' @param scope_obj   A `scope_object` containing `coord$centroids`, segmentation
#'                   vertices, and `cells$counts`.
#' @param gene1_name Character. Name of the first gene.
#' @param gene2_name Character. Name of the second gene.
#' @param palette1   Colour used for cells expressing `gene1`. Default "#fc3d5d".
#' @param palette2   Colour used for cells expressing `gene2`. Default "#4753f8".
#' @param size1      Point size for `gene1` cells. Default 0.3.
#' @param size2      Point size for `gene2` cells. Default 0.3.
#' @param alpha1     Maximum alpha for `gene1` cells. Default 0.7.
#' @param alpha2     Maximum alpha for `gene2` cells. Default 0.7.
#' @param seg_type   Segmentation overlay type: "none", "cell", "nucleus",
#'                   or "both". Default "none".
#' @param colour_cell    Line colour for cell segmentation. Default "black".
#' @param colour_nucleus Line colour for nucleus segmentation. Default "#3182bd".
#' @param alpha_seg  Alpha for segmentation polygons. Default 0.2.
#' @param grid_gap      Spacing of background grid lines (µm). Default 100.
#' @param scale_text_size Font size for scale-bar text. Default 2.4.
#' @param bar_len       Length of the scale bar (µm). Default 400.
#' @param bar_offset    Vertical offset of the scale bar as a fraction of the y-range. Default 0.01.
#' @param arrow_pt      Arrow-head size of the scale bar. Default 4.
#' @param scale_legend_colour Colour of the scale-bar text and line. Default "black".
#' @param max.cutoff1 Fraction (0-1) of gene1's maximum expression used for clipping. Default 1.
#' @param max.cutoff2 Fraction (0-1) of gene2's maximum expression used for clipping. Default 1.
#' @return           A `ggplot` object.
#' @importFrom ggplot2 ggplot geom_point aes scale_alpha_continuous scale_x_continuous scale_y_continuous coord_fixed theme_minimal theme element_blank element_rect element_line element_text unit margin annotate arrow
#' @importFrom ggforce geom_shape
#' @importFrom ggnewscale new_scale
#' @importFrom data.table rbindlist
#' @importFrom scales alpha
#' @export
plotDensityCentroids <- function(scope_obj,
                                 gene1_name,
                                 gene2_name,
                                 palette1 = "#fc3d5d",
                                 palette2 = "#4753f8",
                                 size1 = 0.3,
                                 size2 = 0.3,
                                 alpha1 = 1,
                                 alpha2 = 0.5,
                                 seg_type = c("none", "cell", "nucleus", "both"),
                                 colour_cell = "black",
                                 colour_nucleus = "#3182bd",
                                 alpha_seg = 0.2,
                                 grid_gap = 100,
                                 scale_text_size = 2.4,
                                 bar_len = 400,
                                 bar_offset = 0.01,
                                 arrow_pt = 4,
                                 scale_legend_colour = "black",
                                 max.cutoff1 = 1,
                                 max.cutoff2 = 1) {
    seg_type <- match.arg(seg_type)

    ## ---- 0. counts / centroid ----------
    if (is.null(scope_obj@cells$counts)) {
        stop("scope_obj@cells$counts is missing.")
    }
    counts <- scope_obj@cells$counts
    if (!(gene1_name %in% rownames(counts))) stop("Gene1 not found in counts.")
    if (!(gene2_name %in% rownames(counts))) stop("Gene2 not found in counts.")

    ctd <- scope_obj@coord$centroids
    stopifnot(all(c("cell", "x", "y") %in% names(ctd)))

    expr1 <- counts[gene1_name, ]
    expr2 <- counts[gene2_name, ]

    if (max.cutoff1 < 1) expr1 <- pmin(expr1, max(expr1) * max.cutoff1)
    if (max.cutoff2 < 1) expr2 <- pmin(expr2, max(expr2) * max.cutoff2)

    sc1 <- if (max(expr1) == 0) expr1 else expr1 / max(expr1)
    sc2 <- if (max(expr2) == 0) expr2 else expr2 / max(expr2)

    df1 <- data.frame(
        cell = colnames(counts),
        x    = ctd$x[match(colnames(counts), ctd$cell)],
        y    = ctd$y[match(colnames(counts), ctd$cell)],
        e    = as.numeric(sc1)
    )[sc1 > 0, ]

    df2 <- data.frame(
        cell = colnames(counts),
        x    = ctd$x[match(colnames(counts), ctd$cell)],
        y    = ctd$y[match(colnames(counts), ctd$cell)],
        e    = as.numeric(sc2)
    )[sc2 > 0, ]

    ## ---- 1. points --------------------------------------------------------
    library(ggplot2)
    library(ggnewscale)
    library(ggforce)
    library(data.table)

    p <- ggplot() +
        geom_point(
            data = df1,
            aes(x = x, y = y, alpha = e),
            colour = palette1, size = size1, shape = 16
        ) +
        scale_alpha_continuous(range = c(0, alpha1), guide = "none") +
        ggnewscale::new_scale("alpha") +
        geom_point(
            data = df2,
            aes(x = x, y = y, alpha = e),
            colour = palette2, size = size2, shape = 16
        ) +
        scale_alpha_continuous(range = c(0, alpha2), guide = "none")

    ## ---- 2. segmentation overlay -----------------------------------------
    if (seg_type != "none") {
        seg_layers <- switch(seg_type,
            cell    = "segmentation_cell",
            nucleus = "segmentation_nucleus",
            both    = c("segmentation_cell", "segmentation_nucleus")
        )
        seg_layers <- seg_layers[seg_layers %in% names(scope_obj@coord)]

        if (length(seg_layers)) {
            seg_dt <- data.table::rbindlist(scope_obj@coord[seg_layers],
                use.names = TRUE, fill = TRUE
            )
            seg_dt$cell <- as.character(seg_dt$cell)

            if (seg_type == "both") {
                seg_dt[, type := ifelse(cell %in% scope_obj@coord$segmentation_cell$cell,
                    "cell", "nucleus"
                )]
                p <- p +
                    ggforce::geom_shape(
                        data = seg_dt[type == "cell"],
                        aes(x = x, y = y, group = cell),
                        fill = NA,
                        colour = scales::alpha(colour_cell, alpha_seg),
                        linewidth = .05
                    ) +
                    ggforce::geom_shape(
                        data = seg_dt[type == "nucleus"],
                        aes(x = x, y = y, group = cell),
                        fill = NA,
                        colour = scales::alpha(colour_nucleus, alpha_seg),
                        linewidth = .05
                    )
            } else {
                col_use <- if (seg_type == "cell") colour_cell else colour_nucleus
                p <- p +
                    ggforce::geom_shape(
                        data = seg_dt,
                        aes(x = x, y = y, group = cell),
                        fill = NA,
                        colour = scales::alpha(col_use, alpha_seg),
                        linewidth = .05
                    )
            }
        }
    }

    ## ---- 3. theme / grid / scalebar --------------------------------------
    x_rng <- range(ctd$x, na.rm = TRUE)
    y_rng <- range(ctd$y, na.rm = TRUE)
    dx <- diff(x_rng)
    dy <- diff(y_rng)
    if (dx < dy) {
        pad <- (dy - dx) / 2
        x_rng <- x_rng + c(-pad, pad)
    }
    if (dy < dx) {
        pad <- (dx - dy) / 2
        y_rng <- y_rng + c(-pad, pad)
    }

    grid_x <- seq(x_rng[1], x_rng[2], by = grid_gap)
    grid_y <- seq(y_rng[1], y_rng[2], by = grid_gap)

    p <- p +
        scale_x_continuous(limits = x_rng, expand = c(0, 0), breaks = grid_x) +
        scale_y_continuous(limits = y_rng, expand = c(0, 0), breaks = grid_y) +
        coord_fixed(
            ratio = diff(y_rng) / diff(x_rng),
            xlim = x_rng, ylim = y_rng,
            expand = FALSE, clip = "off"
        ) +
        theme_minimal(base_size = 10) +
        theme(
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = .8),
            axis.line.x.bottom = element_line(colour = "black"),
            axis.line.y.left = element_line(colour = "black"),
            axis.line.x.top = element_blank(),
            axis.line.y.right = element_blank(),
            axis.title = element_blank(),
            plot.margin = margin(20, 40, 40, 40, "pt")
        )

    # scale‑bar
    x0 <- x_rng[1] + 0.001 * diff(x_rng)
    y_bar <- y_rng[1] + bar_offset * diff(y_rng)
    p <- p +
        annotate("segment",
            x = x0, xend = x0 + bar_len,
            y = y_bar, yend = y_bar,
            arrow = arrow(length = unit(arrow_pt, "pt"), ends = "both"),
            colour = scale_legend_colour, linewidth = .4
        ) +
        annotate("text",
            x = x0 + bar_len / 2,
            y = y_bar + 0.025 * diff(y_rng),
            label = paste0(bar_len, " \u00B5m"),
            colour = scale_legend_colour,
            vjust = 1, size = scale_text_size
        )

    p
}
