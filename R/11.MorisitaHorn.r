#' @title Add Morisita–Horn similarities to a gene–gene network
#'
#' @description
#'   Constructs or retrieves a gene-level network based on Lee's L
#'   (\code{L_min} threshold), computes Morisita–Horn similarity along the
#'   selected edges using the sparse C++ kernel
#'   \code{morisita_horn_sparse_cpp()}, and writes the updated network back
#'   into the chosen grid layer of a \code{scope_object}.
#'
#' @param scope_obj   A \code{scope_object} with populated \code{@grid} slot.
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
#' @param ncores     Integer. Number of OpenMP threads for the C++ kernel.
#' @param use_chao   Logical. Apply Chao bias correction when estimating
#'                   lambda terms (default \code{TRUE}).
#' @param verbose    Logical. Whether to print progress messages (default TRUE).
#'
#' @return The modified \code{scope_object} (invisibly).  The graph stored in
#'   \code{graph_slot_mh} gains a numeric edge attribute \code{CMH}.
#'
#' @seealso \code{\link{morisita_horn_sparse_cpp}}
#' @examples
#' \dontrun{
#' coord <- morisitaHorn_on_network2(coord,
#'     grid_name = "25um",
#'     L_min     = 0.2,
#'     ncores    = 8
#' )
#' igraph::edge_attr(coord@grid$`25um`$LeeStats_Xz$g_morisita, "CMH")[1:5]
#' }
#' @importFrom Matrix sparseMatrix Diagonal
#' @importFrom igraph as_edgelist edge_attr graph_from_data_frame
#' @export
computeMH <- function(scope_obj,
                      grid_name = NULL,
                      lee_layer = "LeeStats_Xz",
                      graph_slot = "g_consensus",
                      graph_slot_mh = "g_morisita",
                      matrix_slot_mh = NULL,
                      out = c("matrix", "igraph", "both"),
                      L_min = 0,
                      area_norm = TRUE,
                      ncores = 8,
                      use_chao = TRUE,
                      verbose = TRUE) {
    out <- match.arg(out)
    if (verbose) {
        message("[geneSCOPE::computeMH] Starting Morisita-Horn similarity computation on gene network")
        message("[geneSCOPE::computeMH]   Lee layer: ", lee_layer)
        message("[geneSCOPE::computeMH]   L_min threshold: ", L_min)
        message("[geneSCOPE::computeMH]   Area normalization: ", area_norm)
        message("[geneSCOPE::computeMH]   Chao correction: ", use_chao)
        message("[geneSCOPE::computeMH]   OpenMP threads: ", ncores)
        message("[geneSCOPE::computeMH]   Output: ", out)
    }

    ## —— 1. Select and check grid layer —— ##
    if (verbose) message("[geneSCOPE::computeMH] Selecting and validating grid layer")
    g_layer <- .selectGridLayer(scope_obj, grid_name)
    gname <- names(scope_obj@grid)[
        vapply(scope_obj@grid, identical, logical(1), g_layer)
    ]
    .checkGridContent(scope_obj, gname)

    if (verbose) message("[geneSCOPE::computeMH]   Selected grid layer: ", gname)

    ## —— 2. Get Lee's L —— ##
    if (verbose) message("[geneSCOPE::computeMH] Extracting Lee's L matrix")
    Lmat <- .getLeeMatrix(scope_obj, grid_name = gname, lee_layer = lee_layer)

    if (verbose) {
        message("[geneSCOPE::computeMH]   Lee's L matrix dimensions: ", nrow(Lmat), "×", ncol(Lmat))
        lee_stats <- summary(as.vector(Lmat[upper.tri(Lmat)]))
        message("[geneSCOPE::computeMH]   Lee's L range: [", round(lee_stats[1], 4), ", ", round(lee_stats[6], 4), "]")
    }

    ## —— 3. Try to read existing graph (only when graph_slot is not NULL) —— ##
    if (verbose) message("[geneSCOPE::computeMH] Searching for existing gene network")
    gnet <- NULL
    if (!is.null(graph_slot)) {
        if (verbose) message("[geneSCOPE::computeMH]   Looking for graph in slot: ", graph_slot)
        # ① First check @stats
        if (!is.null(scope_obj@stats) &&
            !is.null(scope_obj@stats[[gname]]) &&
            !is.null(scope_obj@stats[[gname]][[lee_layer]]) &&
            !is.null(scope_obj@stats[[gname]][[lee_layer]][[graph_slot]])) {
            gnet <- scope_obj@stats[[gname]][[lee_layer]][[graph_slot]]
            if (verbose) message("[geneSCOPE::computeMH]   Found existing graph in @stats")
        }
        # ② Then check legacy @grid
        else if (!is.null(g_layer[[lee_layer]]) &&
            !is.null(g_layer[[lee_layer]][[graph_slot]])) {
            gnet <- g_layer[[lee_layer]][[graph_slot]]
            if (verbose) message("[geneSCOPE::computeMH]   Found existing graph in @grid (legacy)")
        }
    } else {
        if (verbose) message("[geneSCOPE::computeMH]   graph_slot is NULL, will build temporary graph")
    }

    ## —— 4. If still no network, build based on L_min —— ##
    if (is.null(gnet)) {
        if (verbose) message("[geneSCOPE::computeMH] Building new network from Lee's L matrix")
        keep <- (Lmat >= L_min)
        diag(keep) <- FALSE
        idx <- which(keep, arr.ind = TRUE)
        if (!nrow(idx)) {
            stop("`L_min` too high, no edges in Lee's L matrix meet criteria.")
        }

        if (verbose) {
            message("[geneSCOPE::computeMH]   Edges meeting L_min threshold: ", nrow(idx))
            message("[geneSCOPE::computeMH]   Edge weight range: [", round(min(Lmat[idx]), 4), ", ", round(max(Lmat[idx]), 4), "]")
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
        message("[geneSCOPE::computeMH]   Final network: ", n_vertices, " vertices, ", n_edges, " edges")
    }

    ## —— 5. gene × grid sparse count matrix (same as legacy) —— ##
    if (verbose) message("[geneSCOPE::computeMH] Constructing gene × grid sparse count matrix")
    counts_dt <- g_layer$counts[count > 0]
    genes <- sort(unique(counts_dt$gene))
    grids <- g_layer$grid_info$grid_id

    if (verbose) {
        message("[geneSCOPE::computeMH]   Unique genes: ", length(genes))
        message("[geneSCOPE::computeMH]   Grid cells: ", length(grids))
        message("[geneSCOPE::computeMH]   Non-zero counts: ", nrow(counts_dt))
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
        message("[geneSCOPE::computeMH]   Sparse matrix: ", nrow(Gsp), "×", ncol(Gsp), " (", round(sparsity, 2), "% sparse)")
    }

    ## —— 5a. Optional area normalization —— ##
    if (area_norm) {
        if (verbose) message("[geneSCOPE::computeMH] Applying area normalization to count matrix")
        gi <- g_layer$grid_info
        if (all(c("width", "height") %in% names(gi))) {
            gi <- gi[match(colnames(Gsp), gi$grid_id), ]
            area_vec <- gi$width * gi$height
            if (length(unique(area_vec)) == 1L) {
                if (verbose) message("[geneSCOPE::computeMH]   Uniform grid area: ", unique(area_vec))
                Gsp@x <- Gsp@x / unique(area_vec)
            } else {
                if (verbose) {
                    area_stats <- summary(area_vec)
                    message("[geneSCOPE::computeMH]   Variable grid areas: [", round(area_stats[1], 3), ", ", round(area_stats[6], 3), "]")
                }
                Gsp <- Gsp %*% Matrix::Diagonal(x = 1 / area_vec)
            }
        } else {
            if (verbose) message("[geneSCOPE::computeMH] Warning: grid_info lacks width/height; skipping area_norm")
        }
    } else {
        if (verbose) message("[geneSCOPE::computeMH] Skipping area normalization (area_norm = FALSE)")
    }

    ## —— 6. Morisita–Horn (C++) —— ##
    if (verbose) message("[geneSCOPE::computeMH] Computing Morisita-Horn similarities using C++ kernel")
    ed <- t(igraph::as_edgelist(gnet, names = FALSE) - 1L)
    storage.mode(ed) <- "integer"

    if (verbose) {
        message("[geneSCOPE::computeMH]   Processing ", ncol(ed), " edges")
        message("[geneSCOPE::computeMH]   C++ threads: ", ncores)
    }

    mh_vec <- morisita_horn_sparse(
        G = Gsp, edges = ed,
        use_chao = use_chao, nthreads = ncores
    )

    if (verbose) {
        # Check if mh_vec contains valid numeric values
        if (is.numeric(mh_vec) && length(mh_vec) > 0) {
            finite_vals <- mh_vec[is.finite(mh_vec)]
            if (length(finite_vals) > 0) {
                mh_min <- min(finite_vals, na.rm = TRUE)
                mh_max <- max(finite_vals, na.rm = TRUE)
                mh_mean <- mean(finite_vals, na.rm = TRUE)
                message("[geneSCOPE::computeMH]   Morisita-Horn similarity range: [", round(mh_min, 4), ", ", round(mh_max, 4), "]")
                message("[geneSCOPE::computeMH]   Mean similarity: ", round(mh_mean, 4))
                if (length(finite_vals) < length(mh_vec)) {
                    n_invalid <- length(mh_vec) - length(finite_vals)
                    message("[geneSCOPE::computeMH]   Warning: ", n_invalid, " non-finite values found")
                }
            } else {
                message("[geneSCOPE::computeMH]   Warning: All similarity values are non-finite")
            }
        } else {
            message("[geneSCOPE::computeMH]   Warning: Similarity computation returned non-numeric results")
        }
    }

    if (out %in% c("igraph", "both")) {
        igraph::edge_attr(gnet, "CMH") <- mh_vec
    }

    ## —— 7. Build sparse CMH matrix aligned to Lee's L gene set —— ##
    gene_all <- rownames(Lmat)
    if (is.null(gene_all)) {
        # Fallback to vertex names if Lmat lacks dimnames
        gene_all <- igraph::V(gnet)$name
    }
    ed_names <- igraph::as_edgelist(gnet, names = TRUE)
    ii <- match(ed_names[,1], gene_all)
    jj <- match(ed_names[,2], gene_all)
    ok <- !is.na(ii) & !is.na(jj)
    if (!all(ok)) {
        if (verbose) message("[geneSCOPE::computeMH] Note: ", sum(!ok), " edges not in L gene set; skipped in matrix build.")
    }
    ii <- as.integer(ii[ok]); jj <- as.integer(jj[ok]); xv <- as.numeric(mh_vec[ok])
    MH <- Matrix::sparseMatrix(i = c(ii, jj), j = c(jj, ii), x = c(xv, xv),
                               dims = c(length(gene_all), length(gene_all)),
                               dimnames = list(gene_all, gene_all))
    MH <- Matrix::drop0(methods::as(MH, "dgCMatrix"))
    if (verbose) {
        message("[geneSCOPE::computeMH]   CMH matrix built: ", nrow(MH), "×", ncol(MH), "; nnz=", Matrix::nnzero(MH))
    }

    ## —— 8. Write to @stats —— ##
    if (verbose) message("[geneSCOPE::computeMH] Storing CMH outputs to stats layer")
    if (is.null(scope_obj@stats)) scope_obj@stats <- list()
    if (is.null(scope_obj@stats[[gname]])) scope_obj@stats[[gname]] <- list()
    if (is.null(scope_obj@stats[[gname]][[lee_layer]])) scope_obj@stats[[gname]][[lee_layer]] <- list()

    # Determine matrix slot name if not provided
    if (is.null(matrix_slot_mh) || !nzchar(matrix_slot_mh)) {
        if (startsWith(graph_slot_mh, "g_")) {
            matrix_slot_mh <- sub("^g_", "m_", graph_slot_mh)
        } else {
            matrix_slot_mh <- paste0("m_", graph_slot_mh)
        }
    }

    if (out %in% c("igraph", "both")) {
        scope_obj@stats[[gname]][[lee_layer]][[graph_slot_mh]] <- gnet
    }
    if (out %in% c("matrix", "both")) {
        scope_obj@stats[[gname]][[lee_layer]][[matrix_slot_mh]] <- MH
    }

    if (verbose) {
        message("[geneSCOPE::computeMH] Analysis completed successfully")
        if (out %in% c("igraph","both")) {
            message("[geneSCOPE::computeMH]   Graph stored in slot: ", graph_slot_mh, " (edge attr 'CMH')")
        }
        if (out %in% c("matrix","both")) {
            message("[geneSCOPE::computeMH]   Matrix stored in slot: ", matrix_slot_mh)
        }
    }

    invisible(scope_obj)
}
