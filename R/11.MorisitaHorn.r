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
#' @param verbose    Logical. Whether to print progress messages (default TRUE).
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
                                  use_chao = TRUE,
                                  verbose = TRUE) {
    if (verbose) {
        message("[geneSCOPE::morisitaHornOnNetwork] Starting Morisita-Horn similarity computation on gene network")
        message("[geneSCOPE::morisitaHornOnNetwork]   Lee layer: ", lee_layer)
        message("[geneSCOPE::morisitaHornOnNetwork]   L_min threshold: ", L_min)
        message("[geneSCOPE::morisitaHornOnNetwork]   Area normalization: ", area_norm)
        message("[geneSCOPE::morisitaHornOnNetwork]   Chao correction: ", use_chao)
        message("[geneSCOPE::morisitaHornOnNetwork]   OpenMP threads: ", ncore)
    }

    ## —— 1. Select and check grid layer —— ##
    if (verbose) message("[geneSCOPE::morisitaHornOnNetwork] Selecting and validating grid layer")
    g_layer <- .selectGridLayer(coordObj, grid_name)
    gname <- names(coordObj@grid)[
        vapply(coordObj@grid, identical, logical(1), g_layer)
    ]
    .checkGridContent(coordObj, gname)

    if (verbose) message("[geneSCOPE::morisitaHornOnNetwork]   Selected grid layer: ", gname)

    ## —— 2. Get Lee's L —— ##
    if (verbose) message("[geneSCOPE::morisitaHornOnNetwork] Extracting Lee's L matrix")
    Lmat <- .getLeeMatrix(coordObj, grid_name = gname, lee_layer = lee_layer)

    if (verbose) {
        message("[geneSCOPE::morisitaHornOnNetwork]   Lee's L matrix dimensions: ", nrow(Lmat), "×", ncol(Lmat))
        lee_stats <- summary(as.vector(Lmat[upper.tri(Lmat)]))
        message("[geneSCOPE::morisitaHornOnNetwork]   Lee's L range: [", round(lee_stats[1], 4), ", ", round(lee_stats[6], 4), "]")
    }

    ## —— 3. Try to read existing graph (only when graph_slot is not NULL) —— ##
    if (verbose) message("[geneSCOPE::morisitaHornOnNetwork] Searching for existing gene network")
    gnet <- NULL
    if (!is.null(graph_slot)) {
        if (verbose) message("[geneSCOPE::morisitaHornOnNetwork]   Looking for graph in slot: ", graph_slot)
        # ① First check @stats
        if (!is.null(coordObj@stats) &&
            !is.null(coordObj@stats[[gname]]) &&
            !is.null(coordObj@stats[[gname]][[lee_layer]]) &&
            !is.null(coordObj@stats[[gname]][[lee_layer]][[graph_slot]])) {
            gnet <- coordObj@stats[[gname]][[lee_layer]][[graph_slot]]
            if (verbose) message("[geneSCOPE::morisitaHornOnNetwork]   Found existing graph in @stats")
        }
        # ② Then check legacy @grid
        else if (!is.null(g_layer[[lee_layer]]) &&
            !is.null(g_layer[[lee_layer]][[graph_slot]])) {
            gnet <- g_layer[[lee_layer]][[graph_slot]]
            if (verbose) message("[geneSCOPE::morisitaHornOnNetwork]   Found existing graph in @grid (legacy)")
        }
    } else {
        if (verbose) message("[geneSCOPE::morisitaHornOnNetwork]   graph_slot is NULL, will build temporary graph")
    }

    ## —— 4. If still no network, build based on L_min —— ##
    if (is.null(gnet)) {
        if (verbose) message("[geneSCOPE::morisitaHornOnNetwork] Building new network from Lee's L matrix")
        keep <- (Lmat >= L_min)
        diag(keep) <- FALSE
        idx <- which(keep, arr.ind = TRUE)
        if (!nrow(idx)) {
            stop("`L_min` too high, no edges in Lee's L matrix meet criteria.")
        }

        if (verbose) {
            message("[geneSCOPE::morisitaHornOnNetwork]   Edges meeting L_min threshold: ", nrow(idx))
            message("[geneSCOPE::morisitaHornOnNetwork]   Edge weight range: [", round(min(Lmat[idx]), 4), ", ", round(max(Lmat[idx]), 4), "]")
        }

        edges_df <- data.frame(
            from = rownames(Lmat)[idx[, 1]],
            to = colnames(Lmat)[idx[, 2]],
            weight = Lmat[idx],
            stringsAsFactors = FALSE
        )
        gnet <- igraph::graph_from_data_frame(edges_df, directed = FALSE)
    }

    if (verbose) {
        n_vertices <- igraph::vcount(gnet)
        n_edges <- igraph::ecount(gnet)
        message("[geneSCOPE::morisitaHornOnNetwork]   Final network: ", n_vertices, " vertices, ", n_edges, " edges")
    }

    ## —— 5. gene × grid sparse count matrix (same as legacy) —— ##
    if (verbose) message("[geneSCOPE::morisitaHornOnNetwork] Constructing gene × grid sparse count matrix")
    counts_dt <- g_layer$counts[count > 0]
    genes <- sort(unique(counts_dt$gene))
    grids <- g_layer$grid_info$grid_id

    if (verbose) {
        message("[geneSCOPE::morisitaHornOnNetwork]   Unique genes: ", length(genes))
        message("[geneSCOPE::morisitaHornOnNetwork]   Grid cells: ", length(grids))
        message("[geneSCOPE::morisitaHornOnNetwork]   Non-zero counts: ", nrow(counts_dt))
    }

    Gsp <- Matrix::sparseMatrix(
        i = match(counts_dt$gene, genes),
        j = match(counts_dt$grid_id, grids),
        x = counts_dt$count,
        dims = c(length(genes), length(grids)),
        dimnames = list(genes, grids)
    )

    if (verbose) {
        sparsity <- (1 - length(Gsp@x) / (nrow(Gsp) * ncol(Gsp))) * 100
        message("[geneSCOPE::morisitaHornOnNetwork]   Sparse matrix: ", nrow(Gsp), "×", ncol(Gsp), " (", round(sparsity, 2), "% sparse)")
    }

    ## —— 5a. Optional area normalization —— ##
    if (area_norm) {
        if (verbose) message("[geneSCOPE::morisitaHornOnNetwork] Applying area normalization to count matrix")
        gi <- g_layer$grid_info
        if (all(c("width", "height") %in% names(gi))) {
            gi <- gi[match(colnames(Gsp), gi$grid_id), ]
            area_vec <- gi$width * gi$height
            if (length(unique(area_vec)) == 1L) {
                if (verbose) message("[geneSCOPE::morisitaHornOnNetwork]   Uniform grid area: ", unique(area_vec))
                Gsp@x <- Gsp@x / unique(area_vec)
            } else {
                if (verbose) {
                    area_stats <- summary(area_vec)
                    message("[geneSCOPE::morisitaHornOnNetwork]   Variable grid areas: [", round(area_stats[1], 3), ", ", round(area_stats[6], 3), "]")
                }
                Gsp <- Gsp %*% Matrix::Diagonal(x = 1 / area_vec)
            }
        } else {
            if (verbose) message("[geneSCOPE::morisitaHornOnNetwork] Warning: grid_info lacks width/height; skipping area_norm")
        }
    } else {
        if (verbose) message("[geneSCOPE::morisitaHornOnNetwork] Skipping area normalization (area_norm = FALSE)")
    }

    ## —— 6. Morisita–Horn (C++) —— ##
    if (verbose) message("[geneSCOPE::morisitaHornOnNetwork] Computing Morisita-Horn similarities using C++ kernel")
    ed <- t(igraph::as_edgelist(gnet, names = FALSE) - 1L)
    storage.mode(ed) <- "integer"

    if (verbose) {
        message("[geneSCOPE::morisitaHornOnNetwork]   Processing ", ncol(ed), " edges")
        message("[geneSCOPE::morisitaHornOnNetwork]   C++ threads: ", ncore)
    }

    mh_vec <- morisita_horn_sparse(
        G = Gsp, edges = ed,
        use_chao = use_chao, nthreads = ncore
    )

    if (verbose) {
        mh_stats <- summary(mh_vec)
        message("[geneSCOPE::morisitaHornOnNetwork]   Morisita-Horn similarity range: [", round(mh_stats[1], 4), ", ", round(mh_stats[6], 4), "]")
        message("[geneSCOPE::morisitaHornOnNetwork]   Mean similarity: ", round(mh_stats[4], 4))
    }

    igraph::edge_attr(gnet, "CMH") <- mh_vec

    ## —— 7. Write to @stats —— ##
    if (verbose) message("[geneSCOPE::morisitaHornOnNetwork] Storing network with Morisita-Horn similarities")
    if (is.null(coordObj@stats)) coordObj@stats <- list()
    if (is.null(coordObj@stats[[gname]])) coordObj@stats[[gname]] <- list()
    if (is.null(coordObj@stats[[gname]][[lee_layer]])) coordObj@stats[[gname]][[lee_layer]] <- list()
    coordObj@stats[[gname]][[lee_layer]][[graph_slot_mh]] <- gnet

    if (verbose) {
        message("[geneSCOPE::morisitaHornOnNetwork] Analysis completed successfully")
        message("[geneSCOPE::morisitaHornOnNetwork]   Network stored in slot: ", graph_slot_mh)
        message("[geneSCOPE::morisitaHornOnNetwork]   Edge attribute 'CMH' contains Morisita-Horn similarities")
    }

    invisible(coordObj)
}
