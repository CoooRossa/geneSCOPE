#' 3.GeneDensity.r (2025-01-06 rev)
#' @title computeDensity - save density tables under @density
#'
#' @description
#'   For a chosen grid sub-layer, sum the counts of a user-defined gene set
#'   (specified directly or via a cluster in \code{coordObj@meta.data}),
#'   optionally normalise the counts, divide by the grid-cell area, and write
#'   the result into \code{coordObj@density[[grid_name]]}.
#'   Each call can add **multiple named density columns** to the same grid-layer
#'   table by using different \code{density_name}s.
#'
#' @param coordObj        A \code{CoordObj} produced by \code{createCoordObj()}.
#' @param grid_name       Character (optional). Target grid sub-layer; if
#'                        \code{NULL} and only one grid layer exists, that one
#'                        is used automatically.
#' @param density_name    Character. Column name to create / overwrite inside
#'                        \code{coordObj@density[[grid_name]]}.
#' @param genes           Character vector (optional). Gene symbols to include.
#'                        If supplied, this takes priority over
#'                        \code{cluster_col}/\code{cluster_num}.
#' @param cluster_col     Character (optional). Column in
#'                        \code{coordObj@meta.data} that stores cluster labels.
#' @param cluster_num     Scalar (optional). The cluster value whose member
#'                        genes are selected when \code{cluster_col} is given.
#' @param layer_name      Character. Which count layer inside the grid to read
#'                        (default \code{"counts"}).
#' @param normalize_method One of \code{"none"}, \code{"per_grid"},
#'                        \code{"global_gene"}.
#'                        * \strong{none}: raw counts.
#'                        * \strong{per_grid}: divide each gene count by the
#'                          total counts in that grid-cell.
#'                        * \strong{global_gene}: divide each gene count by its
#'                          global total across all grids.
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#'
#' @return The modified \code{CoordObj} (invisibly).
#' @importFrom data.table as.data.table copy fifelse
#' @export
computeDensity <- function(coordObj,
                           grid_name = NULL,
                           density_name = "density",
                           genes = NULL,
                           cluster_col = NULL,
                           cluster_num = NULL,
                           layer_name = "counts",
                           normalize_method = c("none", "per_grid", "global_gene"),
                           verbose = TRUE) {
    normalize_method <- match.arg(normalize_method)

    if (verbose) message("[geneSCOPE] Computing density for layer: ", density_name)

    ## --------------------------------------------------------------------- 1
    ## pick grid layer & sanity check
    if (verbose) message("[geneSCOPE] Selecting and validating grid layer...")
    g_layer <- .selectGridLayer(coordObj, grid_name)
    grid_layer_name <- names(coordObj@grid)[vapply(coordObj@grid, identical, logical(1), g_layer)]

    .checkGridContent(coordObj, grid_layer_name)

    grid_info <- g_layer$grid_info
    counts_dt <- g_layer[[layer_name]]
    if (is.null(counts_dt)) {
        stop("Layer '", layer_name, "' not found in grid '", grid_layer_name, "'.")
    }

    if (!all(c("grid_id", "gene", "count") %in% colnames(counts_dt))) {
        stop("Layer '", layer_name, "' must contain columns grid_id, gene, count.")
    }

    counts_dt <- as.data.table(counts_dt)

    ## --------------------------------------------------------------------- 2
    ## decide gene set
    if (verbose) message("[geneSCOPE] Determining gene subset...")
    sel_genes <- .getGeneSubset(coordObj,
        genes        = genes,
        cluster_col  = cluster_col,
        cluster_num  = cluster_num
    )

    if (verbose) message("[geneSCOPE] Selected ", length(sel_genes), " genes for density computation")

    ## --------------------------------------------------------------------- 3
    ## optional normalisation
    if (verbose && normalize_method != "none") {
        message("[geneSCOPE] Applying normalization method: ", normalize_method)
    }

    counts_proc <- copy(counts_dt)

    if (normalize_method == "per_grid") {
        tot_grid <- counts_proc[, .(total = sum(count)), by = grid_id]
        counts_proc <- counts_proc[tot_grid, on = "grid_id"]
        counts_proc[, count := fifelse(total == 0, 0, count / total)][, total := NULL]
    } else if (normalize_method == "global_gene") {
        tot_gene <- counts_proc[, .(g_tot = sum(count)), by = gene]
        counts_proc <- counts_proc[tot_gene, on = "gene"]
        counts_proc[, count := fifelse(g_tot == 0, 0, count / g_tot)][, g_tot := NULL]
    }
    ## else "none": keep raw counts

    ## --------------------------------------------------------------------- 4
    ## sum across selected genes → per‑grid counts
    if (verbose) message("[geneSCOPE] Aggregating counts across selected genes...")
    sel_counts <- counts_proc[gene %in% sel_genes,
        .(count = sum(count)),
        by = grid_id
    ]

    ## --------------------------------------------------------------------- 5
    ## attach cell area & compute density
    if (verbose) message("[geneSCOPE] Computing density values...")
    grid_dt <- as.data.table(grid_info)[
        ,
        .(grid_id, xmin, xmax, ymin, ymax,
            area = (xmax - xmin) * (ymax - ymin)
        )
    ]

    dens_dt <- merge(grid_dt, sel_counts, by = "grid_id", all.x = TRUE)
    dens_dt[is.na(count), count := 0]
    dens_dt[, density := count / area]

    ## --------------------------------------------------------------------- 6
    ## create / update @density entry for this grid layer
    if (verbose) message("[geneSCOPE] Storing density results...")
    all_ids <- grid_dt$grid_id
    dens_df <- if (!is.null(coordObj@density[[grid_layer_name]])) {
        coordObj@density[[grid_layer_name]]
    } else {
        data.frame(row.names = all_ids)
    }

    # ensure all grids represented
    if (nrow(dens_df) == 0) {
        dens_df <- data.frame(row.names = all_ids)
    }

    dens_vec <- dens_dt$density
    names(dens_vec) <- dens_dt$grid_id
    dens_df[[density_name]] <- dens_vec[rownames(dens_df)]
    dens_df[[density_name]][is.na(dens_df[[density_name]])] <- 0

    coordObj@density[[grid_layer_name]] <- dens_df

    if (verbose) message("[geneSCOPE] Density computation completed for column: ", density_name)

    invisible(coordObj)
}
