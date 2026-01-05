#' @title Extract paths from dendrogram/network results
#' @description
#' Process a \code{plotDendroNetwork} output object to extract informative paths:
#' \itemize{
#'   \item Only \code{gene} provided: return all paths that contain the gene (prefer paths seen in random walks; if no walks available, enumerate simple paths in the graph up to \code{cutoff} and filter those containing the gene). The gene is not required to be an endpoint.
#'   \item Both \code{gene} and \code{gene2} provided: return all paths that contain both genes (prefer random-walk paths; otherwise enumerate simple paths and filter). The two genes are not forced to be endpoints.
#' }
#' @details
#' Full-graph enumeration considers all (source,target) with source<target using \code{igraph::all_simple_paths} (constrained by \code{cutoff} and \code{max_paths}). For large graphs this can be expensive; reduce \code{cutoff} and/or \code{max_paths} for safety.
#' @param dnet_obj A plotDendroNetwork output object
#' @param gene A gene symbol
#' @param gene2 Optional second gene symbol
#' @param all_shortest Deprecated; kept for compatibility but ignored
#' @param cutoff Maximum number of nodes per simple path considered
#' @param max_paths Maximum number of paths to return (stop once reached)
#' @param verbose Logical; print progress messages (default: getOption("geneSCOPE.verbose", TRUE))
#' @return A list with selected paths and summary data
#' @importFrom igraph V all_simple_paths
#' @export
getDendroWalkPaths <- function(dnet_obj,
                               gene,
                               gene2 = NULL,
                               all_shortest = FALSE, # deprecated; placeholder only
                               cutoff = 8,
                               max_paths = 50000,
                               verbose = getOption("geneSCOPE.verbose", TRUE)) {
    ## 1. Locate the igraph object --------------------------------------------------
    g <- NULL
    for (nm in c("graph", "g", "igraph")) {
        if (!is.null(dnet_obj[[nm]])) {
            g <- dnet_obj[[nm]]
            break
        }
    }
    if (is.null(g) || !inherits(g, "igraph")) stop("No igraph graph found (fields graph/g/igraph)")

    vn <- igraph::V(g)$name
    if (is.null(vn)) stop("Graph vertices are missing names (V(g)$name)")

    gene <- as.character(gene)[1]
    if (!(gene %in% vn)) stop("gene not found in graph: ", gene)
    pair_mode <- !is.null(gene2)
    if (pair_mode) {
        gene2 <- as.character(gene2)[1]
        if (!(gene2 %in% vn)) stop("gene2 not found in graph: ", gene2)
        if (gene2 == gene) stop("gene and gene2 must be different")
    }

    ## 2. Normalize random-walk paths ----------------------------------------------
    has_walks <- FALSE
    walks_std <- NULL
    if (!is.null(dnet_obj$walks) && is.list(dnet_obj$walks)) {
        rw <- dnet_obj$walks
        rw <- rw[vapply(rw, function(x) length(x) >= 2, logical(1))]
        if (length(rw)) {
            walks_std <- lapply(rw, function(w) {
                if (is.factor(w)) w <- as.character(w)
                if (is.numeric(w)) {
                    if (all(w >= 1 & w <= length(vn))) vn[w] else as.character(w)
                } else {
                    as.character(w)
                }
            })
            has_walks <- TRUE
        }
    }

    ## 3. Helper: enumerate simple paths and filter --------------------------------
    enumerate_paths_filter <- function(filter_fun) {
        # filter_fun: function(char_vector_path) -> TRUE/FALSE
        nV <- length(vn)
        kept <- list()
        k <- 0L
        # source<target to avoid duplicates from reversed direction
        for (s in seq_len(nV - 1L)) {
            for (t in (s + 1L):nV) {
                ps <- igraph::all_simple_paths(g, from = s, to = t, mode = "all", cutoff = cutoff)
                if (length(ps)) {
                    for (pp in ps) {
                        path_chr <- vn[pp]
                        if (filter_fun(path_chr)) {
                            k <- k + 1L
                            kept[[k]] <- path_chr
                            if (k >= max_paths) {
                                if (verbose) message("[geneSCOPE] Reached max_paths=", max_paths, "; stopping enumeration")
                                return(kept)
                            }
                        }
                    }
                }
            }
        }
        kept
    }

    ## 4A. Single-gene mode: all paths containing `gene` ----------------------------
    if (!pair_mode) {
        if (has_walks) {
            sel <- walks_std[vapply(walks_std, function(p) gene %in% p, logical(1))]
            if (verbose) message("[geneSCOPE] Selected ", length(sel), " paths containing ", gene, " from random walks")
        } else {
            if (verbose) message("[geneSCOPE] No random walks found; enumerating graph paths (cutoff=", cutoff, ", max_paths=", max_paths, ")")
            sel <- enumerate_paths_filter(function(p) gene %in% p)
            if (verbose) message("[geneSCOPE] Enumerated and filtered ", length(sel), " paths")
        }
        # de-duplicate
        if (length(sel)) {
            key <- vapply(sel, paste, collapse = "->", FUN.VALUE = character(1))
            sel <- sel[!duplicated(key)]
        }
        paths_df <- if (length(sel)) {
            data.frame(
                path_id = seq_along(sel),
                length = vapply(sel, length, integer(1)),
                path_str = vapply(sel, paste, collapse = "->", FUN.VALUE = character(1)),
                stringsAsFactors = FALSE
            )
        } else {
            data.frame(path_id = integer(0), length = integer(0), path_str = character(0))
        }
        return(list(
            mode = "single_gene",
            genes = gene,
            paths = sel,
            paths_df = paths_df,
            n_paths = nrow(paths_df),
            source_has_walks = has_walks,
            cutoff_used = if (has_walks) NA_integer_ else cutoff
        ))
    }

    ## 4B. Dual-gene mode: all paths containing both `gene` and `gene2` -------------
    if (has_walks) {
        sel <- walks_std[vapply(walks_std, function(p) all(c(gene, gene2) %in% p), logical(1))]
        if (verbose) message("[geneSCOPE] Selected ", length(sel), " paths containing both ", gene, " and ", gene2, " from random walks")
    } else {
        if (verbose) message("[geneSCOPE] No random walks found; enumerating and filtering paths containing both genes (cutoff=", cutoff, ", max_paths=", max_paths, ")")
        sel <- enumerate_paths_filter(function(p) (gene %in% p) && (gene2 %in% p))
        if (verbose) message("[geneSCOPE] Enumerated and filtered ", length(sel), " paths")
    }
    # de-duplicate
    if (length(sel)) {
        key <- vapply(sel, paste, collapse = "->", FUN.VALUE = character(1))
        sel <- sel[!duplicated(key)]
    }
    paths_df <- if (length(sel)) {
        data.frame(
            path_id = seq_along(sel),
            length = vapply(sel, length, integer(1)),
            path_str = vapply(sel, paste, collapse = "->", FUN.VALUE = character(1)),
            stringsAsFactors = FALSE
        )
    } else {
        data.frame(path_id = integer(0), length = integer(0), path_str = character(0))
    }

    list(
        mode = "dual_gene_all_paths",
        genes = c(gene, gene2),
        paths = sel,
        paths_df = paths_df,
        n_paths = nrow(paths_df),
        source_has_walks = has_walks,
        cutoff_used = if (has_walks) NA_integer_ else cutoff
    )
}
