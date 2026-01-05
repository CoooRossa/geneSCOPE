#' @title Compute Morisita's Iδ per gene and store in a scope_object
#'
#' @description
#'   Builds a sparse gene × grid/cell count matrix, computes Morisita's
#'   index Iδ for every gene using the parallel C++ kernel
#'   \code{idelta_sparse_cpp()}, and writes the raw scores back to
#'   \code{scope_obj}.  The result is also copied into
#'   \code{scope_obj@meta.data} under a column named
#'   \code{<grid_name>_iDelta_raw} or \code{cell_iDelta_raw}.
#'
#' @param scope_obj  A \code{scope_object} that already contains a populated
#'                  \code{@grid} slot with \code{counts} and
#'                  \code{grid_info}, or \code{@cells} slot with count matrix.
#' @param grid_name Character. Name of the grid sub-layer to use when
#'                  \code{level = "grid"}. If \code{NULL} and only one layer
#'                  exists it is chosen automatically.
#' @param level     Character. Either \code{"grid"} (default) to compute Iδ
#'                  across spatial grids, or \code{"cell"} to compute across
#'                  single cells.
#' @param ncores    Integer. Number of OpenMP threads passed to the
#'                  C++ routine (default 4).
#' @param verbose   Logical. Whether to print progress messages (default TRUE).
#'
#' @return The modified \code{scope_object} is returned invisibly.
#'
#' @details
#'   When \code{level = "grid"}, only grid–gene pairs with positive molecule
#'   counts are considered. When \code{level = "cell"}, the function uses the
#'   count matrix from \code{scope_obj@cells$counts}. In both cases, memory
#'   usage scales with the number of non-zero entries. Genes with fewer than
#'   two total molecules receive \code{NaN}.
#'
#' @examples
#' \dontrun{
#' # Grid-level Iδ
#' coord <- computeIDelta(coord, grid_name = "25um", level = "grid", ncores = 8)
#'
#' # Single-cell Iδ
#' coord <- computeIDelta(coord, level = "cell", ncores = 8)
#' }
#' @importFrom Matrix sparseMatrix
#' @export
computeIDelta <- function(scope_obj,
                          grid_name = NULL,
                          level = c("grid", "cell"),
                          ncores = 4,
                          verbose = getOption("geneSCOPE.verbose", TRUE)) {
    level <- match.arg(level)
    ncores <- max(1L, min(as.integer(ncores), parallel::detectCores(logical = TRUE)))

    if (level == "grid") {
        if (verbose) message("[geneSCOPE::computeIDelta] Selecting grid layer and building sparse gene-by-grid matrix")
        ## ---- 1. Select grid layer ---------------------------------------------------
        g_layer <- .selectGridLayer(scope_obj, grid_name)
        # Parse actual grid_name (selectGridLayer might auto-select)
        grid_name <- names(scope_obj@grid)[
            vapply(scope_obj@grid, identical, logical(1), g_layer)
        ]

        ## ---- 2. Build gene × grid sparse matrix ------------------------------------
        dt <- g_layer$counts[g_layer$counts$count > 0, ] # Keep only positive counts
        genes <- sort(unique(dt$gene))
        grids <- g_layer$grid_info$grid_id

        Gsp <- Matrix::sparseMatrix(
            i = match(dt$gene, genes),
            j = match(dt$grid_id, grids),
            x = dt$count,
            dims = c(length(genes), length(grids)),
            dimnames = list(genes, grids)
        )

        storage_key <- grid_name
        col_suffix <- paste0(grid_name, "_iDelta")
    } else { # level == "cell"
        if (verbose) message("[geneSCOPE::computeIDelta] Using existing gene-by-cell sparse count matrix")
        ## ---- 1. Check for cell count matrix ----------------------------------------
        if (is.null(scope_obj@cells) || is.null(scope_obj@cells$counts)) {
            stop("No cell count matrix found. Please run addSingleCells() first.")
        }

        ## ---- 2. Use existing gene × cell sparse matrix -----------------------------
        Gsp <- scope_obj@cells$counts # Already genes × cells
        if (!inherits(Gsp, "dgCMatrix")) {
            stop("Cell count matrix must be a dgCMatrix.")
        }

        genes <- rownames(Gsp)

        storage_key <- "cell"
        col_suffix <- "cell_iDelta"
    }

    ## ---- 3. Call C++ kernel to compute Iδ -----------------------------------------
    # Memory guard: approximate sparse footprint (values + indices) per thread
    nnz <- length(Gsp@x)
    approx_gb <- (nnz * 16) / (1024^3) # rough estimate of numeric + index bytes
    sys_mem_gb <- .getSystemMemoryGB()
    est_total_gb <- approx_gb * ncores
    if (est_total_gb > sys_mem_gb) {
        stop(
            "[geneSCOPE::computeIDelta] Estimated memory requirement (",
            round(est_total_gb, 1), " GB) exceeds system capacity (",
            round(sys_mem_gb, 1), " GB). Reduce ncores or gene/cell counts."
        )
    }

    if (verbose) message("[geneSCOPE::computeIDelta] Computing Iδ for ", length(genes), " genes using ", ncores, " thread(s)")
    delta_raw <- idelta_sparse_cpp(Gsp, n_threads = ncores)
    names(delta_raw) <- genes

    ## ---- 4. Ensure delta_raw is a numeric vector, not matrix ----------------------
    delta_raw <- as.numeric(delta_raw)
    names(delta_raw) <- genes

    ## ---- 5. Write to scope_obj -----------------------------------------------------
    if (verbose) message("[geneSCOPE::computeIDelta] Writing Iδ results into scope_obj")
    # 5a. @stats new location
    if (is.null(scope_obj@stats)) scope_obj@stats <- list()
    if (is.null(scope_obj@stats[[storage_key]])) {
        scope_obj@stats[[storage_key]] <- list()
    }
    scope_obj@stats[[storage_key]]$iDeltaStats <- list(
        genes      = genes,
        delta_raw  = delta_raw, # Ensure this is a numeric vector
        level      = level
    )

    # 5b. @meta.data maintain old behavior (column name based on level)
    if (is.null(scope_obj@meta.data)) {
        scope_obj@meta.data <- data.frame(row.names = genes)
    }

    # Ensure all genes have rows in meta.data
    missing_genes <- setdiff(genes, rownames(scope_obj@meta.data))
    if (length(missing_genes) > 0) {
        new_rows <- data.frame(row.names = missing_genes)
        scope_obj@meta.data <- rbind(scope_obj@meta.data, new_rows)
    }

    # Assign scalar values (not matrix) to each gene
    scope_obj@meta.data[genes, col_suffix] <- delta_raw[genes]

    if (verbose) message("[geneSCOPE::computeIDelta] Iδ computation completed")
    invisible(scope_obj)
}






#' @title Plot Morisita Iδ peRr cluster
#'
#' @description
#'   Creates a faceted scatter-line plot showing Morisita's Iδ values for the
#'   top genes in each cluster of a chosen grid layer.  Clusters can be
#'   down-sampled, re-ordered to balance facet rows, and coloured by a custom
#'   palette.
#'
#' @param scope_obj        A \code{scope_object} containing \code{@grid} and
#'                        \code{@meta.data}.
#' @param grid_name       Character. Grid sub-layer name. If \code{NULL} and
#'                        only one layer exists, that layer is used.
#' @param cluster_col     Column name in \code{meta.data} that assigns each
#'                        gene to a cluster.
#' @param suffix          Either \code{"raw"} or \code{"nor"}; chooses between
#'                        raw or normalised Iδ.
#' @param top_n           Optional integer. Keep the top \code{n} genes
#'                        (highest Iδ) per cluster.
#' @param nrow            Integer. Number of facet rows to aim for.
#' @param min_genes       Integer. Minimum number of genes a cluster must keep
#'                        after filtering (default 1).
#' @param cluster_palette Character vector of colours to use for clusters. If
#'                        unnamed the order is taken as-is; otherwise names
#'                        are matched to cluster labels.
#' @param point_size,line_size,label_size Numeric. Aesthetics for the plot.
#' @param subCluster      Optional character vector. Restrict to these cluster
#'                        labels.
#'
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 ggplot aes geom_point geom_line facet_wrap labs theme_minimal theme element_blank element_rect margin coord_cartesian unit
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr group_by slice_max ungroup filter arrange desc mutate
#' @export
plotIDelta <- function(
    scope_obj,
    grid_name,
    cluster_col,
    top_n = NULL,
    min_genes = 1,
    nrow = 1,
    point_size = 3,
    line_size = 0.5,
    label_size = 2, # default reduced to 2
    subCluster = NULL) {
    meta_col <- paste0(grid_name, "_iDelta")
    if (!meta_col %in% colnames(scope_obj@meta.data)) {
        stop("Cannot find Iδ values: meta.data column '", meta_col, "' is missing.")
    }

    meta <- scope_obj@meta.data
    genes <- rownames(meta)
    df <- data.frame(
        gene = genes,
        delta = as.numeric(meta[genes, meta_col]),
        cluster = meta[genes, cluster_col],
        stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$cluster), ]

    if (!is.null(subCluster)) {
        df <- subset(df, cluster %in% subCluster)
    }
    if (!is.null(top_n) && top_n > 0) {
        df <- df %>%
            group_by(cluster) %>%
            slice_max(order_by = delta, n = top_n) %>%
            ungroup()
    }
    df <- df %>%
        group_by(cluster) %>%
        filter(n() >= min_genes) %>%
        ungroup()
    if (nrow(df) == 0) {
        stop("No cluster meets min_genes ≥ ", min_genes, ".")
    }
    df <- df %>%
        group_by(cluster) %>%
        arrange(desc(delta)) %>%
        mutate(gene = factor(gene, levels = gene)) %>%
        ungroup()

    p <- ggplot(df, aes(x = gene, y = delta, group = cluster, color = factor(cluster))) +
        geom_point(size = point_size) +
        geom_line(size = line_size) +
        # Use geom_text_repel to avoid overlap; if don't want ggrepel dependency, change back to geom_text()
        geom_text_repel(
            aes(label = gene),
            size = label_size,
            box.padding = unit(0.15, "lines"),
            point.padding = unit(0.15, "lines"),
            segment.size = 0.3,
            segment.color = "grey50",
            force = 0.5,
            max.overlaps = Inf
        ) +
        facet_wrap(~cluster, scales = "free_x", nrow = nrow, strip.position = "bottom") +
        labs(
            x     = NULL, # remove x-axis label as well
            y     = expression(I[delta]),
            title = paste0("Iδ by Cluster (", grid_name, ")")
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_blank(), # Remove x-axis text
            panel.grid = element_blank(), # Remove grid lines
            panel.border = element_rect(colour = "black", fill = NA), # Panel border
            axis.ticks.x = element_blank(), # Remove x-axis ticks (optional)
            plot.margin = margin(5, 20, 5, 5), # Increase right margin to prevent label clipping
            strip.background = element_rect(fill = "white", colour = NA)
        ) +
        coord_cartesian(clip = "off") # Allow labels to extend beyond plot area

    return(p)
}
