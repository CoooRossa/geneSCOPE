#' @title Add Morisita–Horn similarities to a gene–gene network
#'
#' @description
#'   Constructs or retrieves a gene-level network based on Lee's L
#'   (\code{L_min} threshold), computes Morisita–Horn similarity along the
#'   selected edges using the sparse C++ kernel
#'   \code{morisita_horn_sparse_cpp()}, and writes the updated network back
#'   into the chosen grid layer of a \code{CoordObj}.
#'
#' @param coordObj   A \code{CoordObj} with populated \code{@grid} slot.
#' @param grid_name  Character. Grid sub-layer to use. If \code{NULL} and only
#'                   one layer exists it is chosen automatically.
#' @param lee_layer  Name of the layer that contains Lee's statistics
#'                   (default \code{"LeeStats_Xz"}).
#' @param graph_slot Name of the element that stores, or will store, the
#'                   consensus graph built from Lee's L (default
#'                   \code{"g_consensus"}). If \code{NULL}, a temporary graph is
#'                   constructed on the fly.
#' @param graph_slot_mh Name of the element in which to save the graph with
#'                   Morisita–Horn edge weights (default \code{"g_morisita"}).
#' @param L_min      Numeric. Minimum Lee's L value required for an edge.
#' @param area_norm  Logical. Divide counts by grid area before similarity
#'                   calculation (default \code{TRUE}).
#' @param ncore      Integer. Number of OpenMP threads for the C++ kernel.
#' @param use_chao   Logical. Apply Chao bias correction when estimating
#'                   lambda terms (default \code{TRUE}).
#'
#' @return The modified \code{CoordObj} (invisibly).  The graph stored in
#'   \code{graph_slot_mh} gains a numeric edge attribute \code{CMH}.
#'
#' @seealso \code{\link{morisita_horn_sparse_cpp}}
#' @examples
#' \dontrun{
#' coord <- morisitaHorn_on_network2(coord,
#'     grid_name = "25um",
#'     L_min     = 0.2,
#'     ncore     = 8
#' )
#' igraph::edge_attr(coord@grid$`25um`$LeeStats_Xz$g_morisita, "CMH")[1:5]
#' }
#' @importFrom Matrix sparseMatrix Diagonal
#' @importFrom igraph as_edgelist edge_attr graph_from_data_frame
#' @export
morisitaHornOnNetwork <- function(coordObj,
                                  grid_name = NULL,
                                  lee_layer = "LeeStats_Xz",
                                  graph_slot = "g_consensus",
                                  graph_slot_mh = "g_morisita",
                                  L_min = 0,
                                  area_norm = TRUE,
                                  ncore = 8,
                                  use_chao = TRUE) {
    ## —— 1. Select and check grid layer —— ##
    g_layer <- .selectGridLayer(coordObj, grid_name)
    gname <- names(coordObj@grid)[
        vapply(coordObj@grid, identical, logical(1), g_layer)
    ]
    .checkGridContent(coordObj, gname)

    ## —— 2. Get Lee's L —— ##
    Lmat <- .getLeeMatrix(coordObj, grid_name = gname, lee_layer = lee_layer)

    ## —— 3. Try to read existing graph (only when graph_slot is not NULL) —— ##
    gnet <- NULL
    if (!is.null(graph_slot)) {
        # ① First check @stats
        if (!is.null(coordObj@stats) &&
            !is.null(coordObj@stats[[gname]]) &&
            !is.null(coordObj@stats[[gname]][[lee_layer]]) &&
            !is.null(coordObj@stats[[gname]][[lee_layer]][[graph_slot]])) {
            gnet <- coordObj@stats[[gname]][[lee_layer]][[graph_slot]]
        }
        # ② Then check legacy @grid
        else if (!is.null(g_layer[[lee_layer]]) &&
            !is.null(g_layer[[lee_layer]][[graph_slot]])) {
            gnet <- g_layer[[lee_layer]][[graph_slot]]
        }
    }

    ## —— 4. If still no network, build based on L_min —— ##
    if (is.null(gnet)) {
        keep <- (Lmat >= L_min)
        diag(keep) <- FALSE
        idx <- which(keep, arr.ind = TRUE)
        if (!nrow(idx)) {
            stop("`L_min` too high, no edges in Lee's L matrix meet criteria.")
        }
        edges_df <- data.frame(
            from = rownames(Lmat)[idx[, 1]],
            to = colnames(Lmat)[idx[, 2]],
            weight = Lmat[idx],
            stringsAsFactors = FALSE
        )
        gnet <- igraph::graph_from_data_frame(edges_df, directed = FALSE)
    }

    ## —— 5. gene × grid sparse count matrix (same as legacy) —— ##
    counts_dt <- g_layer$counts[count > 0]
    genes <- sort(unique(counts_dt$gene))
    grids <- g_layer$grid_info$grid_id
    Gsp <- Matrix::sparseMatrix(
        i = match(counts_dt$gene, genes),
        j = match(counts_dt$grid_id, grids),
        x = counts_dt$count,
        dims = c(length(genes), length(grids)),
        dimnames = list(genes, grids)
    )

    ## —— 5a. Optional area normalization —— ##
    if (area_norm) {
        gi <- g_layer$grid_info
        if (all(c("width", "height") %in% names(gi))) {
            gi <- gi[match(colnames(Gsp), gi$grid_id), ]
            area_vec <- gi$width * gi$height
            if (length(unique(area_vec)) == 1L) {
                Gsp@x <- Gsp@x / unique(area_vec)
            } else {
                Gsp <- Gsp %*% Matrix::Diagonal(x = 1 / area_vec)
            }
        } else {
            if (verbose) message("[geneSCOPE::morisitaHornOnNetwork] !!! Warning: grid_info lacks width/height; skipping area_norm !!!")
        }
    }

    ## —— 6. Morisita–Horn (C++) —— ##
    ed <- t(igraph::as_edgelist(gnet, names = FALSE) - 1L)
    storage.mode(ed) <- "integer"
    mh_vec <- morisita_horn_sparse(
        G = Gsp, edges = ed,
        use_chao = use_chao, nthreads = ncore
    )
    igraph::edge_attr(gnet, "CMH") <- mh_vec

    ## —— 7. Write to @stats —— ##
    if (is.null(coordObj@stats)) coordObj@stats <- list()
    if (is.null(coordObj@stats[[gname]])) coordObj@stats[[gname]] <- list()
    if (is.null(coordObj@stats[[gname]][[lee_layer]])) coordObj@stats[[gname]][[lee_layer]] <- list()
    coordObj@stats[[gname]][[lee_layer]][[graph_slot_mh]] <- gnet

    invisible(coordObj)
}
