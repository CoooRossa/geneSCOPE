#' @title Visualise the Lee's-L gene-gene network
#'
#' @description
#'   Builds a gene network from Lee's L statistics or from a pre-computed
#'   consensus graph, applies multiple optional filters (L thresholds, FDR,
#'   LR-curve confidence bands, quantile cut-off), and renders the graph with
#'   attractive aesthetics.  Edge widths encode |L|, edge colours fade toward
#'   the stronger node's cluster colour, and dashed lines mark negative
#'   correlations when \code{show_sign = TRUE}.  Nodes are coloured by supplied
#'   cluster labels or a meta-data column; hubs receive a thicker outline.
#'
#' @param scope_obj        A \code{scope_object} containing \code{@grid} and Lee's
#'                        statistics.
#' @param grid_name       Character. Grid sub-layer to plot. Defaults to the
#'                        first layer.
#' @param lee_stats_layer Name of the Lee's stats layer (default \code{"LeeStats_Xz"}).
#' @param gene_subset     Optional character vector of genes to keep.
#' @param L_min,L_min_neg Minimum positive / negative |L| thresholds.
#' @param p_cut,FDR_max   p-value or FDR cut-offs when the corresponding
#'                        matrices are available.
#' @param pct_min         Quantile string such as \code{"q80"} passed to the
#'                        internal quantile filter.
#' @param CI95_filter     Logical. Remove edges according to the 95 % LR-curve
#'                        band using \code{CI_rule}.
#' @param curve_layer     Name of the LR-curve layer for CI filtering.
#' @param CI_rule         \code{"remove_within"} or \code{"remove_outside"}.
#' @param drop_isolated   Logical. Discard isolated nodes.
#' @param weight_abs      Use absolute L for edge weight (default \code{TRUE}).
#' @param use_consensus_graph Logical. Use pre-computed consensus graph stored
#'                        in \code{graph_slot_name}.
#' @param graph_slot_name Name under which the consensus graph is stored.
#' @param cluster_vec     Either a named vector of cluster assignments or the
#'                        name of a column in \code{meta.data}.
#' @param cluster_palette Character vector of colours (named or unnamed).
#' @param vertex_size,base_edge_mult,label_cex Graph aesthetics.
#' @param layout_niter    Iterations for the Fruchterman-Reingold layout.
#' @param seed            Random seed for layout reproducibility.
#' @param hub_factor      Degree multiplier defining "hub" nodes.
#' @param length_scale    Global edge-length multiplier.
#' @param show_sign       Draw negative edges with distinct linetype.
#' @param neg_linetype    Linetype for negative edges.
#' @param title           Optional plot title.
#'
#' @return A \code{ggplot} object.
#' @importFrom ggraph ggraph create_layout geom_edge_link geom_node_point geom_node_text scale_edge_width scale_edge_colour_identity scale_edge_linetype_manual
#' @importFrom ggplot2 aes labs guides guide_legend coord_fixed theme_void theme element_text unit margin
#' @importFrom igraph V E induced_subgraph graph_from_adjacency_matrix get.edgelist subgraph.edges degree
#' @importFrom grDevices col2rgb rgb rainbow colorRampPalette
#' @importFrom dplyr filter
#' @importFrom stats approxfun quantile hclust as.dendrogram
#' @importFrom Matrix drop0
#' @export
plotNetwork <- function(
    scope_obj,
    ## ---------- Data layers ----------
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    gene_subset = NULL, # Now compatible with direct gene vector or (cluster_col, cluster_num)
    ## ---------- Filter thresholds ----------
    L_min = 0,
    L_min_neg = NULL,
    p_cut = NULL,
    use_FDR = TRUE,
    FDR_max = 0.05,
    pct_min = "q0",
    CI95_filter = FALSE,
    curve_layer = "LR_curve",
    CI_rule = c("remove_within", "remove_outside"),
    drop_isolated = TRUE,
    weight_abs = TRUE,
    ## ---------- Consensus network ----------
    use_consensus_graph = TRUE,
    graph_slot_name = "g_consensus",
    ## ---------- New FDR source ----------
    fdr_source = c("FDR", "FDR_storey", "FDR_beta", "FDR_mid", "FDR_uniform", "FDR_disc"),
    ## ---------- Plot parameters ----------
    cluster_vec = NULL,
    cluster_palette = c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33",
        "#A65628", "#984EA3", "#66C2A5", "#FC8D62", "#8DA0CB",
        "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"
    ),
    vertex_size = 8,
    base_edge_mult = 12,
    label_cex = 3,
    layout_niter = 1000,
    seed = 1,
    hub_factor = 3,
    length_scale = 1,
    max.overlaps = 20,
    L_max = 1,
    hub_border_col = "#4B4B4B",
    hub_border_size = 0.8,
    show_sign = FALSE,
    neg_linetype = "dashed",
    neg_legend_lab = "Negative",
    pos_legend_lab = "Positive",
    show_qc_caption = TRUE,
    title = NULL) {
    CI_rule <- match.arg(CI_rule)
    fdr_source <- match.arg(fdr_source)

    ## ===== 0. Lock and validate grid layer =====
    g_layer <- .selectGridLayer(scope_obj, grid_name) # If grid_name=NULL auto select unique layer
    grid_name <- if (is.null(grid_name)) {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    } else {
        as.character(grid_name)
    }

    .checkGridContent(scope_obj, grid_name) # Error directly if missing required fields

    ## ===== 1. Read LeeStats object and its matrices =====
    ##  First locate LeeStats layer (can be in @stats[[grid]] or @grid[[grid]])
    leeStat <- if (!is.null(scope_obj@stats[[grid_name]]) &&
        !is.null(scope_obj@stats[[grid_name]][[lee_stats_layer]])) {
        scope_obj@stats[[grid_name]][[lee_stats_layer]]
    } else if (!is.null(g_layer[[lee_stats_layer]])) {
        g_layer[[lee_stats_layer]]
    } else {
        stop(
            "Cannot find layer '", lee_stats_layer,
            "' in grid '", grid_name, "'."
        )
    }

    Lmat <- .getLeeMatrix(scope_obj,
        grid_name = grid_name,
        lee_layer = lee_stats_layer
    ) # Only extract aligned L

    ## ===== 2. Gene subset =====
    all_genes <- rownames(Lmat)
    if (is.null(gene_subset)) {
        keep_genes <- all_genes
    } else if (is.character(gene_subset)) {
        ## (A) Directly provide gene vector
        keep_genes <- intersect(
            all_genes,
            .getGeneSubset(scope_obj, genes = gene_subset)
        )
    } else if (is.list(gene_subset) &&
        all(c("cluster_col", "cluster_num") %in% names(gene_subset))) {
        ## (B) Use cluster column + number syntax: gene_subset = list(cluster_col = "..", cluster_num = ..)
        keep_genes <- intersect(
            all_genes,
            .getGeneSubset(scope_obj,
                cluster_col = gene_subset$cluster_col,
                cluster_num = gene_subset$cluster_num
            )
        )
    } else {
        stop("`gene_subset` must be a character vector or list(cluster_col, cluster_num)")
    }
    if (length(keep_genes) < 2) {
        stop("Less than two genes remain after sub‑setting.")
    }

    ## ===== 3. Retrieve consensus graph or rebuild =====
    use_consensus <- isTRUE(use_consensus_graph) &&
        !is.null(leeStat[[graph_slot_name]])

    if (use_consensus && isTRUE(use_FDR)) {
        # Consensus graphs are typically pre-filtered; skip redundant FDR filtering
        message("[geneSCOPE::plotNetwork] use_consensus_graph=TRUE ignoring FDR re-filtering (assuming graph already filtered).")
    }

    if (use_consensus) {
        g_raw <- leeStat[[graph_slot_name]]
        g <- igraph::induced_subgraph(
            g_raw,
            intersect(igraph::V(g_raw)$name, keep_genes)
        )
        if (igraph::ecount(g) == 0) {
            stop("Consensus graph has no edges under the chosen gene subset.")
        }
        ## Update edge weights / signs when required
    } else {
        ## ---------- Reconstruct graph using thresholds ----------
        idx <- match(keep_genes, rownames(Lmat))
        A <- Lmat[idx, idx, drop = FALSE]

        ## (i) Early filter by expression proportion (pct_min)
        A <- .filter_matrix_by_quantile(A, pct_min, "q100") # Fixed function name

        ## (ii) Cap the maximum absolute value
        A[abs(A) > L_max] <- 0

        ## (iii) Apply positive/negative thresholds
        L_min_neg <- if (is.null(L_min_neg)) L_min else abs(L_min_neg)
        if (weight_abs) {
            A[abs(A) < L_min] <- 0
        } else {
            A[A > 0 & A < L_min] <- 0
            A[A < 0 & abs(A) < L_min_neg] <- 0
        }

        ## (iv) 95% confidence interval filter
        if (CI95_filter) {
            curve <- leeStat[[curve_layer]]
            if (is.null(curve)) {
                stop("Cannot find `curve_layer = ", curve_layer, "` in LeeStats.")
            }
            f_lo <- approxfun(curve$Pear, curve$lo95, rule = 2)
            f_hi <- approxfun(curve$Pear, curve$hi95, rule = 2)

            rMat <- .getPearsonMatrix(scope_obj, level = "cell") # use single-cell Pearson correlations
            rMat <- rMat[keep_genes, keep_genes, drop = FALSE]

            mask <- if (CI_rule == "remove_within") {
                (A >= f_lo(rMat)) & (A <= f_hi(rMat))
            } else {
                (A < f_lo(rMat)) | (A > f_hi(rMat))
            }
            A[mask] <- 0
        }

        ## (v) p‑value & FDR
        if (!is.null(p_cut) && !is.null(leeStat$P)) {
            Pmat <- leeStat$P[idx, idx, drop = FALSE]
            A[Pmat >= p_cut | is.na(Pmat)] <- 0
        }
        if (isTRUE(use_FDR)) {
            ## -------- FDR matrix selection logic --------
            pref_order <- c(
                fdr_source,
                setdiff(
                    c("FDR", "FDR_storey", "FDR_beta", "FDR_mid", "FDR_uniform", "FDR_disc"),
                    fdr_source
                )
            )
            FDR_sel <- NULL
            FDR_used_name <- NULL
            for (nm in pref_order) {
                cand <- leeStat[[nm]]
                if (!is.null(cand)) {
                    FDR_sel <- cand
                    FDR_used_name <- nm
                    break
                }
            }
            if (is.null(FDR_sel)) {
                stop(
                    "use_FDR=TRUE but no usable FDR matrix found (expected: ",
                    paste(pref_order, collapse = ", "), ")."
                )
            }
            if (inherits(FDR_sel, "big.matrix")) {
                message(
                    "[geneSCOPE::plotNetwork] Using FDR source '", FDR_used_name,
                    "' (big.matrix) → converting to regular matrix for subset filtering"
                )
                FDR_sel <- bigmemory::as.matrix(FDR_sel)
            } else if (!is.matrix(FDR_sel)) {
                FDR_sel <- as.matrix(FDR_sel)
            }
            if (!identical(dim(FDR_sel), dim(Lmat))) {
                stop("FDR matrix dimensions do not match L: ", FDR_used_name)
            }
            FDRmat <- FDR_sel[idx, idx, drop = FALSE]
            A[FDRmat > FDR_max | is.na(FDRmat)] <- 0
            message(sprintf("[geneSCOPE::plotNetwork] Using FDR source '%s' (FDR_max=%.3g)", FDR_used_name, FDR_max))
        }

        ## (vi) Symmetrise and zero the diagonal
        A <- (A + t(A)) / 2
        diag(A) <- 0

        ## (vii) Drop isolated vertices
        if (drop_isolated) {
            keep <- which(rowSums(abs(A) > 0) > 0 | colSums(abs(A) > 0) > 0)
            A <- A[keep, keep, drop = FALSE]
        }
        if (nrow(A) < 2 || all(A == 0)) {
            stop("No edges remain after filtering thresholds.")
        }

        g <- igraph::graph_from_adjacency_matrix(abs(A),
            mode = "undirected",
            weighted = TRUE, diag = FALSE
        )
        e_idx <- igraph::as_edgelist(g, names = FALSE) # replacement for deprecated get.edgelist
        igraph::E(g)$sign <- ifelse(A[e_idx] < 0, "neg", "pos")
    }

    ## ===== 4. Weight adjustment & global scaling =====
    igraph::E(g)$weight <- abs(igraph::E(g)$weight)
    bad <- which(!is.finite(igraph::E(g)$weight) | igraph::E(g)$weight <= 0)
    if (length(bad)) {
        min_pos <- min(igraph::E(g)$weight[igraph::E(g)$weight > 0], na.rm = TRUE)
        igraph::E(g)$weight[bad] <- ifelse(is.finite(min_pos), min_pos, 1)
    }
    if (!is.null(length_scale) && length_scale != 1) {
        igraph::E(g)$weight <- igraph::E(g)$weight * length_scale
    }

    if (!is.null(L_min) && L_min > 0) {
        keep_e <- which(igraph::E(g)$weight >= L_min)
        if (length(keep_e) == 0) stop("`L_min` too strict; no edges to plot.")
        g <- igraph::subgraph.edges(g, keep_e, delete.vertices = TRUE)
    }

    ## ===== 5. Node color / cluster labels =====
    Vnames <- igraph::V(g)$name
    deg_vec <- igraph::degree(g)

    ## ① Derive clu vector from cluster_vec or meta.data
    clu <- setNames(rep(NA_character_, length(Vnames)), Vnames)
    is_factor_input <- FALSE
    factor_levels <- NULL

    ## helper – extract column from meta.data
    get_meta_cluster <- function(col) {
        if (is.null(scope_obj@meta.data) || !(col %in% colnames(scope_obj@meta.data))) {
            stop("Column ‘", col, "’ not found in scope_obj@meta.data.")
        }
        tmp <- scope_obj@meta.data[[col]]
        names(tmp) <- rownames(scope_obj@meta.data)
        tmp[Vnames]
    }

    if (!is.null(cluster_vec)) {
        if (length(cluster_vec) > 1) { # manual vector
            if (is.factor(cluster_vec)) {
                is_factor_input <- TRUE
                factor_levels <- levels(cluster_vec)
            }
            tmp <- as.character(cluster_vec[Vnames])
            clu[!is.na(tmp)] <- tmp[!is.na(tmp)]
        } else { # column name
            meta_clu <- get_meta_cluster(cluster_vec)
            if (is.factor(meta_clu)) {
                is_factor_input <- TRUE
                factor_levels <- levels(meta_clu)
            }
            clu[!is.na(meta_clu)] <- as.character(meta_clu[!is.na(meta_clu)])
        }
    }

    keep_nodes <- names(clu)[!is.na(clu)]
    if (length(keep_nodes) < 2) {
        stop("Less than two nodes have cluster labels.")
    }
    g <- igraph::induced_subgraph(g, keep_nodes)
    Vnames <- igraph::V(g)$name
    clu <- clu[Vnames]
    deg_vec <- deg_vec[Vnames]




    ## Step 2: generate palette
    uniq_clu <- if (is_factor_input && !is.null(factor_levels)) {
        intersect(factor_levels, unique(clu))
    } else {
        sort(unique(clu))
    }
    n_clu <- length(uniq_clu)
    palette_vals <- if (is.null(cluster_palette)) {
        setNames(rainbow(n_clu), uniq_clu)
    } else if (is.null(names(cluster_palette))) {
        setNames(colorRampPalette(cluster_palette)(n_clu), uniq_clu)
    } else {
        pal <- cluster_palette
        miss <- setdiff(uniq_clu, names(pal))
        if (length(miss)) {
            pal <- c(
                pal,
                setNames(colorRampPalette(cluster_palette)(length(miss)), miss)
            )
        }
        pal[uniq_clu]
    }
    basecol <- setNames(palette_vals[clu], Vnames)

    ## Step 3: edge colours (sign / intensity gradient)
    e_idx <- igraph::as_edgelist(g, names = FALSE)
    w_norm <- igraph::E(g)$weight / max(igraph::E(g)$weight)
    cluster_ord <- setNames(seq_along(uniq_clu), uniq_clu)
    edge_cols <- character(igraph::ecount(g))
    for (i in seq_along(edge_cols)) {
        v1 <- Vnames[e_idx[i, 1]]
        v2 <- Vnames[e_idx[i, 2]]
        if (show_sign && igraph::E(g)$sign[i] == "neg") {
            edge_cols[i] <- "gray40"
        } else { # use gradient for positive correlations
            ref_col <- if (!is.na(clu[v1]) && !is.na(clu[v2]) && clu[v1] == clu[v2]) {
                basecol[v1]
            } else {
                # decide primary colour by degree or cluster order
                if (deg_vec[v1] > deg_vec[v2]) {
                    basecol[v1]
                } else if (deg_vec[v2] > deg_vec[v1]) {
                    basecol[v2]
                } else if (!is.na(clu[v1]) && is.na(clu[v2])) {
                    basecol[v1]
                } else if (!is.na(clu[v2]) && is.na(clu[v1])) {
                    basecol[v2]
                } else if (!is.na(clu[v1]) && !is.na(clu[v2])) {
                    if (cluster_ord[clu[v1]] <= cluster_ord[clu[v2]]) {
                        basecol[v1]
                    } else {
                        basecol[v2]
                    }
                } else {
                    "gray80"
                }
            }
            rgb_ref <- grDevices::col2rgb(ref_col) / 255
            tval <- w_norm[i]
            edge_cols[i] <- grDevices::rgb(
                1 - tval + tval * rgb_ref[1],
                1 - tval + tval * rgb_ref[2],
                1 - tval + tval * rgb_ref[3]
            )
        }
    }
    igraph::E(g)$edge_col <- edge_cols
    igraph::E(g)$linetype <- if (show_sign) igraph::E(g)$sign else "solid"

    ## ===== 6. ggraph rendering =====
    library(ggraph)
    library(ggplot2)
    library(dplyr)
    set.seed(seed)
    lay <- ggraph::create_layout(g, layout = "fr", niter = layout_niter)
    lay$basecol <- basecol[lay$name]
    lay$deg <- deg_vec[lay$name]
    lay$hub <- lay$deg > hub_factor * median(lay$deg)

    qc_txt <- if (show_qc_caption && !is.null(leeStat$qc)) {
        with(leeStat$qc, sprintf(
            "density = %.3f | comp = %d | Q = %.2f | mean‑deg = %.1f ± %.1f | hubs = %.1f%% | sig.edge = %.1f%%",
            edge_density, components, modularity_Q,
            mean_degree, sd_degree, hub_ratio * 100, sig_edge_frac * 100
        ))
    } else {
        NULL
    }

    ## ---- edge width / colour / linetype layers ----
    p <- ggraph(lay) +
        geom_edge_link(
            aes(
                width = weight, colour = edge_col,
                linetype = linetype
            ),
            lineend = "round",
            show.legend = c(
                width = TRUE, colour = FALSE,
                linetype = show_sign
            )
        ) +
        scale_edge_width(name = "|L|", range = c(0.3, base_edge_mult)) +
        scale_edge_colour_identity() + {
            if (show_sign) {
                scale_edge_linetype_manual(
                    name   = "Direction",
                    values = c(pos = "solid", neg = neg_linetype),
                    breaks = c("pos", "neg"),
                    labels = c(pos = pos_legend_lab, neg = neg_legend_lab)
                )
            }
        }

    ## ---- node points / text ----
    p <- p +
        geom_node_point(
            data = ~ dplyr::filter(.x, !hub),
            aes(size = deg, fill = basecol), shape = 21, stroke = 0,
            show.legend = c(fill = TRUE, size = TRUE)
        ) +
        geom_node_point(
            data = ~ dplyr::filter(.x, hub),
            aes(size = deg, fill = basecol),
            shape = 21, colour = hub_border_col, stroke = hub_border_size,
            show.legend = FALSE
        ) +
        geom_node_text(aes(label = name),
            size = label_cex,
            repel = TRUE, vjust = 1.4,
            max.overlaps = max.overlaps
        ) +
        scale_size_continuous(
            name = "Node size",
            range = c(vertex_size / 2, vertex_size * 1.5),
            guide = guide_legend(override.aes = list(shape = 19))
        ) +
        scale_fill_identity(
            name = "Cluster",
            guide = guide_legend(
                override.aes = list(shape = 21, size = vertex_size * 0.5),
                nrow = 2, order = 1
            ),
            breaks = unname(palette_vals),
            labels = names(palette_vals)
        ) +
        coord_fixed() + theme_void(base_size = 12) +
        theme(
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 15),
            plot.title = element_text(size = 15, hjust = 0.5),
            legend.position = "bottom",
            legend.box = "vertical",
            legend.box.just = "center",
            legend.spacing.y = unit(0.2, "cm"),
            legend.box.margin = margin(t = 5, b = 5),
            plot.caption = element_text(
                hjust = 0, size = 9,
                margin = margin(t = 6)
            )
        ) +
        labs(title = title, caption = qc_txt)

    ## Optional legend order adjustment
    p <- p +
        guides(
            fill = guide_legend(order = 1),
            size = guide_legend(order = 2),
            edge_width = guide_legend(
                order = 3, direction = "horizontal",
                nrow = 1
            ),
            linetype = if (show_sign) guide_legend(order = 4, nrow = 1) else "none"
        )

    return(p)
}



#' @title Plot Dendrogram-style Network Layout
#'
#' @description
#'   Creates a tree-like visualization of a gene network using minimum spanning
#'   tree algorithms within and between clusters, optionally weighted by
#'   personalized PageRank scores derived from Morisita's Iδ values.
#'
#' @inheritParams plotNetwork
#' @param IDelta_col_name Character. Column name in \code{meta.data} containing
#'                        Morisita's Iδ values for PageRank weighting. If \code{NULL},
#'                        no PageRank weighting is applied.
#' @param damping         Numeric. Damping factor for PageRank algorithm (default 0.85).
#' @param weight_low_cut  Numeric. Minimum edge weight threshold after PageRank
#'                        weighting (default 0).
#' @param k_top           Integer. Maximum number of high-weight non-tree edges
#'                        to retain between clusters (default 1).
#' @param tree_mode       Character. Tree layout style: \code{"rooted"},
#'                        \code{"radial"}, or \code{"forest"} (default "rooted").
#'
#' @return Invisibly returns a list with components:
#'   \describe{
#'     \item{\code{graph}}{The processed \code{igraph} object.}
#'     \item{\code{pagerank}}{PageRank scores if \code{IDelta_col_name} was provided.}
#'     \item{\code{cross_edges}}{Data frame of inter-cluster edges in the final graph.}
#'   }
#'
#' @details
#'   The algorithm first constructs minimum spanning trees within each cluster,
#'   then connects clusters using a second MST on the cluster-level graph.
#'   Optionally, high-weight non-tree edges between clusters are preserved
#'   based on \code{k_top}. If \code{IDelta_col_name} is provided, edge weights
#'   are re-scaled by the mean PageRank of incident vertices.
#'
#' @importFrom igraph mst components subgraph.edges delete_edges delete_vertices degree E V as_data_frame ends edge_attr page_rank induced_subgraph ecount vcount
#' @importFrom stats quantile median
#' @export
plotDendroNetwork <- function(
    scope_obj,
    ## ---------- Data layers ----------
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    gene_subset = NULL,
    ## ---------- Filtering thresholds ----------
    L_min = 0,
    drop_isolated = TRUE,
    ## ---------- Network source ----------
    use_consensus_graph = TRUE,
    graph_slot_name = "g_consensus",
    ## ---------- Cluster labels ----------
    cluster_vec = NULL,
    ## ---------- Δ-PageRank options ----------
    IDelta_col_name = NULL, # New option; NULL disables Δ-PageRank
    damping = 0.85,
    weight_low_cut = 0,
    ## ---------- Tree construction ----------
    ## ---------- Visual details (re-use defaults) ----------
    cluster_palette = c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00",
        "#FFFF33", "#A65628", "#984EA3", "#66C2A5",
        "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
        "#FFD92F", "#E5C494"
    ),
    vertex_size = 8,
    base_edge_mult = 12,
    label_cex = 3,
    seed = 1,
    hub_factor = 3,
    length_scale = 1,
    max.overlaps = 20,
    hub_border_col = "#4B4B4B",
    hub_border_size = 0.8,
    show_sign = FALSE,
    neg_linetype = "dashed",
    neg_legend_lab = "Negative",
    pos_legend_lab = "Positive",
    show_qc_caption = TRUE,
    title = NULL,
    k_top = 1,
    tree_mode = c("rooted", "radial", "forest")) {
    tree_layout <- TRUE # keep tree layout
    ## ========= 0. Read consensus graph ========
    g_layer <- .selectGridLayer(scope_obj, grid_name)
    grid_name <- if (is.null(grid_name)) {
        names(scope_obj@grid)[vapply(scope_obj@grid, identical, logical(1), g_layer)]
    } else {
        grid_name
    }
    leeStat <- if (!is.null(scope_obj@stats[[grid_name]]) &&
        !is.null(scope_obj@stats[[grid_name]][[lee_stats_layer]])) {
        scope_obj@stats[[grid_name]][[lee_stats_layer]]
    } else {
        g_layer[[lee_stats_layer]]
    }

    g_raw <- leeStat[[graph_slot_name]]
    if (is.null(g_raw)) {
        stop(
            "plotDendroNetwork: consensus graph '", graph_slot_name,
            "' not found in leeStat[[\"", graph_slot_name, "\"]]; ensure clusterGenes has populated this layer."
        )
    }
    stopifnot(inherits(g_raw, "igraph"))

    ## ========= 1. Subset genes =========
    keep_genes <- rownames(scope_obj@meta.data)
    if (!is.null(gene_subset)) {
        keep_genes <- intersect(
            keep_genes,
            .getGeneSubset(scope_obj, genes = gene_subset)
        )
    }
    g <- igraph::induced_subgraph(g_raw, intersect(V(g_raw)$name, keep_genes))
    if (drop_isolated) g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))
    if (igraph::vcount(g) < 2) stop("Subgraph contains fewer than two vertices.")

    ## ========= 2. Δ-PageRank reweighting =========
    if (!is.null(IDelta_col_name)) {
        delta <- scope_obj@meta.data[V(g)$name, IDelta_col_name, drop = TRUE]
        delta[is.na(delta)] <- median(delta, na.rm = TRUE)
        q <- stats::quantile(delta, c(.1, .9))
        delta <- pmax(pmin(delta, q[2]), q[1]) # align with dendroRW behaviour
        pers <- exp(delta - max(delta))
        pers <- pers / sum(pers)
        pr <- igraph::page_rank(g,
            personalized = pers,
            damping = damping,
            weights = igraph::E(g)$weight,
            directed = FALSE
        )$vector
        et <- igraph::as_data_frame(g, "edges")
        w_rw <- igraph::E(g)$weight * ((pr[et$from] + pr[et$to]) / 2)
        w_rw[w_rw <= weight_low_cut] <- 0
        igraph::E(g)$weight <- w_rw
    }

    ## ========= 3. Global rescaling / L_min =========
    if (!is.null(length_scale) && length_scale != 1) {
        igraph::E(g)$weight <- igraph::E(g)$weight * length_scale
    }
    if (!is.null(L_min) && L_min > 0) {
        keep_e <- which(igraph::E(g)$weight >= L_min)
        g <- igraph::subgraph.edges(g, keep_e, delete.vertices = TRUE)
    }

    ## ========= 4. Retrieve cluster labels =========
    Vnames <- V(g)$name
    clu <- rep(NA_character_, length(Vnames))
    names(clu) <- Vnames
    if (!is.null(cluster_vec)) {
        cv <- if (length(cluster_vec) == 1) {
            scope_obj@meta.data[Vnames, cluster_vec, drop = TRUE]
        } else {
            cluster_vec[Vnames]
        }
        clu[!is.na(cv)] <- as.character(cv[!is.na(cv)])
    }
    g <- igraph::induced_subgraph(g, names(clu)[!is.na(clu)])
    Vnames <- V(g)$name
    clu <- clu[Vnames]

    ## ========= 5. Build cluster MST backbone =========
    ## —— 5.1 Intra-cluster MST ——
    all_edges <- igraph::as_data_frame(g, "edges")
    all_edges$key <- with(
        all_edges,
        ifelse(from < to, paste(from, to, sep = "|"),
            paste(to, from, sep = "|")
        )
    )
    keep_key <- character(0)

    for (cl in unique(clu)) {
        vsub <- Vnames[clu == cl]
        if (length(vsub) < 2) next
        g_sub <- igraph::induced_subgraph(g, vsub)
        if (ecount(g_sub) == 0) next
        mst_sub <- igraph::mst(g_sub, weights = 1 / (E(g_sub)$weight + 1e-9))
        ks <- igraph::as_data_frame(mst_sub, "edges")
        ks$key <- with(
            ks,
            ifelse(from < to, paste(from, to, sep = "|"),
                paste(to, from, sep = "|")
            )
        )
        keep_key <- c(keep_key, ks$key)
    }

    ## —— 5.2 Inter-cluster MST ——
    if (length(unique(clu)) > 1) {
        ed <- all_edges
        ed$cl1 <- clu[ed$from]
        ed$cl2 <- clu[ed$to]
        inter <- ed[ed$cl1 != ed$cl2, ]
        if (nrow(inter)) {
            inter$pair <- ifelse(inter$cl1 < inter$cl2,
                paste(inter$cl1, inter$cl2, sep = "|"),
                paste(inter$cl2, inter$cl1, sep = "|")
            )
            agg <- aggregate(weight ~ pair + cl1 + cl2, data = inter, max)
            g_clu <- igraph::graph_from_data_frame(
                agg[, c("cl1", "cl2", "weight")],
                directed = FALSE, vertices = unique(clu)
            )
            ## Obtain MST for each connected component
            cmp <- igraph::components(g_clu)$membership
            for (cc in unique(cmp)) {
                sub <- igraph::induced_subgraph(g_clu, which(cmp == cc))
                if (ecount(sub) == 0) next
                mstc <- igraph::mst(sub, weights = 1 / (E(sub)$weight + 1e-9))
                ks <- igraph::as_data_frame(mstc, "edges")
                ks$key <- with(
                    ks,
                    ifelse(from < to, paste(from, to, sep = "|"),
                        paste(to, from, sep = "|")
                    )
                )
                # retain the original edge with maximal weight
                for (k in ks$key) {
                    cand <- inter[inter$pair == k, ]
                    cand <- cand[order(cand$weight, decreasing = TRUE), ]
                    keep_key <- c(keep_key, cand$key[1])
                }
            }
        }
    }

    ## —— 5.3 Filter edges based on keep_key ——
    keep_eid <- which(all_edges$key %in% unique(keep_key))
    g <- igraph::subgraph.edges(g, keep_eid, delete.vertices = TRUE)
    g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))

    ## ============ Tree-layout helper steps ==============
    ## Retain original edge IDs for later mapping
    igraph::E(g)$eid <- seq_len(igraph::ecount(g))

    ## ---- 6.1  Intra-cluster MST ----
    keep_eid <- integer(0)
    for (cl in unique(clu)) {
        vsub <- Vnames[clu == cl]
        if (length(vsub) < 2) next
        g_sub <- igraph::induced_subgraph(g, vsub)
        if (igraph::ecount(g_sub) > 0) {
            mst_sub <- igraph::mst(g_sub, weights = igraph::E(g_sub)$weight)
            keep_eid <- c(keep_eid, igraph::edge_attr(mst_sub, "eid"))
        }
    }

    ## ---- 6.2  Inter-cluster MST on the cluster graph ----
    ## ---- 6.2  Add back high-weight off-tree edges ----
    ## Parameters:
    ##   k_top        : maximum number of off-tree edges to reproject (may be NULL)
    ##   w_extra_min  : minimum off-tree edge weight to reproject (may be NULL)
    if (length(unique(clu)) > 1) {
        ed_tab <- igraph::as_data_frame(g, what = "edges")
        ed_tab$eid <- igraph::E(g)$eid
        ed_tab$cl1 <- clu[ed_tab$from]
        ed_tab$cl2 <- clu[ed_tab$to]
        inter <- ed_tab[ed_tab$cl1 != ed_tab$cl2, ] # consider only inter-cluster edges

        if (nrow(inter)) {
            ## ---- 6.2a  Build the cluster-level graph and take its MST ----
            inter$pair <- ifelse(inter$cl1 < inter$cl2,
                paste(inter$cl1, inter$cl2, sep = "|"),
                paste(inter$cl2, inter$cl1, sep = "|")
            )
            agg <- aggregate(weight ~ pair + cl1 + cl2, data = inter, max)
            g_clu <- igraph::graph_from_data_frame(
                agg[, c("cl1", "cl2", "weight")],
                directed = FALSE, vertices = unique(clu)
            )

            ## For disconnected cases, take an MST per component
            mst_clu <- igraph::mst(
                g_clu,
                weights = 1 / (igraph::E(g_clu)$weight + 1e-9)
            )

            # Extract endpoint pairs from that MST‐forest:
            ep <- igraph::ends(mst_clu, igraph::E(mst_clu))
            keep_pairs <- ifelse(
                ep[, 1] < ep[, 2],
                paste(ep[, 1], ep[, 2], sep = "|"),
                paste(ep[, 2], ep[, 1], sep = "|")
            )
            inter_MST <- inter[inter$pair %in% keep_pairs, ]

            ## ---- 6.2b  Identify high-weight off-tree edges ----
            extra_edges <- inter[!(inter$pair %in% keep_pairs), ]
            if (nrow(extra_edges)) {
                if (!is.null(w_extra_min)) {
                    extra_edges <- extra_edges[extra_edges$weight >= w_extra_min, ]
                }
                if (!is.null(k_top) && nrow(extra_edges) > k_top) {
                    extra_edges <- extra_edges[order(extra_edges$weight, decreasing = TRUE)[seq_len(k_top)], ]
                }
            }
            ## ---- 6.2c  Collect edge IDs to retain ----
            keep_eid <- c(
                keep_eid,
                inter_MST$eid,
                if (nrow(extra_edges)) extra_edges$eid else integer(0)
            )
        }
    }
    keep_eid <- unique(keep_eid)
    g <- igraph::delete_edges(g, igraph::E(g)[!eid %in% keep_eid])
    ## Remove isolated vertices that remain
    g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))
    ## ============ Tree layout complete ==============
    ## —— 6.3  Derive inter-cluster edges for plotting ——
    # Convert the current graph into a data.frame of edges
    edf_final <- igraph::as_data_frame(g, what = "edges")
    # Keep only edges whose endpoints belong to different clusters
    cross_edges <- edf_final[clu[edf_final$from] != clu[edf_final$to], ]


    ## ===== 7. Recompute degree / colours =====
    Vnames <- igraph::V(g)$name
    deg_vec <- igraph::degree(g)


    is_factor_input <- FALSE
    factor_levels <- NULL

    ## helper – extract column from meta.data
    get_meta_cluster <- function(col) {
        if (is.null(scope_obj@meta.data) || !(col %in% colnames(scope_obj@meta.data))) {
            stop("Column ‘", col, "’ not found in scope_obj@meta.data.")
        }
        tmp <- scope_obj@meta.data[[col]]
        names(tmp) <- rownames(scope_obj@meta.data)
        tmp[Vnames]
    }

    if (!is.null(cluster_vec)) {
        if (length(cluster_vec) > 1) { # manual vector
            if (is.factor(cluster_vec)) {
                is_factor_input <- TRUE
                factor_levels <- levels(cluster_vec)
            }
            tmp <- as.character(cluster_vec[Vnames])
            clu[!is.na(tmp)] <- tmp[!is.na(tmp)]
        } else { # column name
            meta_clu <- get_meta_cluster(cluster_vec)
            if (is.factor(meta_clu)) {
                is_factor_input <- TRUE
                factor_levels <- levels(meta_clu)
            }
            clu[!is.na(meta_clu)] <- as.character(meta_clu[!is.na(meta_clu)])
        }
    }


    uniq_clu <- if (is_factor_input && !is.null(factor_levels)) {
        intersect(factor_levels, unique(clu))
    } else {
        sort(unique(clu))
    }
    n_clu <- length(uniq_clu)
    palette_vals <- if (is.null(cluster_palette)) {
        setNames(rainbow(n_clu), uniq_clu)
    } else if (is.null(names(cluster_palette))) {
        setNames(colorRampPalette(cluster_palette)(n_clu), uniq_clu)
    } else {
        pal <- cluster_palette
        miss <- setdiff(uniq_clu, names(pal))
        if (length(miss)) {
            pal <- c(
                pal,
                setNames(colorRampPalette(cluster_palette)(length(miss)), miss)
            )
        }
        pal[uniq_clu]
    }
    basecol <- setNames(palette_vals[clu], Vnames)

    ## ===== 8. Edge colours (legacy behaviour) =====
    e_idx <- igraph::as_edgelist(g, names = FALSE)
    w_norm <- igraph::E(g)$weight / max(igraph::E(g)$weight)
    cluster_ord <- setNames(seq_along(uniq_clu), uniq_clu)
    edge_cols <- character(igraph::ecount(g))
    for (i in seq_along(edge_cols)) {
        v1 <- Vnames[e_idx[i, 1]]
        v2 <- Vnames[e_idx[i, 2]]
        if (show_sign && igraph::E(g)$sign[i] == "neg") {
            edge_cols[i] <- "gray40"
        } else {
            ref_col <- if (!is.na(clu[v1]) && !is.na(clu[v2]) && clu[v1] == clu[v2]) {
                basecol[v1]
            } else {
                if (deg_vec[v1] > deg_vec[v2]) {
                    basecol[v1]
                } else if (deg_vec[v2] > deg_vec[v1]) {
                    basecol[v2]
                } else if (!is.na(clu[v1]) && is.na(clu[v2])) {
                    basecol[v1]
                } else if (!is.na(clu[v2]) && is.na(clu[v1])) {
                    basecol[v2]
                } else if (!is.na(clu[v1]) && !is.na(clu[v2])) {
                    if (cluster_ord[clu[v1]] <= cluster_ord[clu[v2]]) {
                        basecol[v1]
                    } else {
                        basecol[v2]
                    }
                } else {
                    "gray80"
                }
            }
            rgb_ref <- grDevices::col2rgb(ref_col) / 255
            tval <- w_norm[i]
            edge_cols[i] <- grDevices::rgb(
                1 - tval + tval * rgb_ref[1],
                1 - tval + tval * rgb_ref[2],
                1 - tval + tval * rgb_ref[3]
            )
        }
    }
    igraph::E(g)$edge_col <- edge_cols
    igraph::E(g)$linetype <- if (show_sign) igraph::E(g)$sign else "solid"

    ## ===== 9. ggraph rendering =====
    library(ggraph)
    library(ggplot2)
    library(dplyr)
    set.seed(seed)

    tree_mode <- match.arg(tree_mode) # recently added parameter; default is "rooted"

    if (tree_mode == "rooted") {
        root_v <- V(g)[which.max(deg_vec)] # same approach as before
        lay <- create_layout(g, layout = "tree", root = root_v)
    } else if (tree_mode == "radial") {
        lay <- create_layout(g, layout = "tree", circular = TRUE)
    } else if (tree_mode == "forest") {
        root_v <- V(g)[which.max(deg_vec)]
        lay <- create_layout(g, layout = "tree", root = root_v, circular = TRUE)
    }
    lay$basecol <- basecol[lay$name]
    lay$deg <- deg_vec[lay$name]
    lay$hub <- lay$deg > hub_factor * median(lay$deg)

    qc_txt <- if (show_qc_caption && !is.null(leeStat$qc)) {
        with(leeStat$qc, sprintf(
            "density = %.3f | comp = %d | Q = %.2f | mean‑deg = %.1f ± %.1f | hubs = %.1f%% | sig.edge = %.1f%%",
            edge_density, components, modularity_Q,
            mean_degree, sd_degree, hub_ratio * 100, sig_edge_frac * 100
        ))
    } else {
        NULL
    }

    p <- ggraph(lay) +
        geom_edge_link(aes(width = weight, colour = edge_col, linetype = linetype),
            lineend = "round",
            show.legend = c(
                width = TRUE, colour = FALSE,
                linetype = show_sign
            )
        ) +
        scale_edge_width(name = "|L|", range = c(0.3, base_edge_mult)) +
        scale_edge_colour_identity() +
        {
            if (show_sign) {
                scale_edge_linetype_manual(
                    name = "Direction",
                    values = c(pos = "solid", neg = neg_linetype),
                    breaks = c("pos", "neg"),
                    labels = c(pos = pos_legend_lab, neg = neg_legend_lab)
                )
            }
        } +
        geom_node_point(
            data = ~ dplyr::filter(.x, !hub),
            aes(size = deg, fill = basecol), shape = 21, stroke = 0,
            show.legend = c(fill = TRUE, size = TRUE)
        ) +
        geom_node_point(
            data = ~ dplyr::filter(.x, hub),
            aes(size = deg, fill = basecol),
            shape = 21, colour = hub_border_col, stroke = hub_border_size,
            show.legend = FALSE
        ) +
        geom_node_text(aes(label = name),
            size = label_cex,
            repel = TRUE,
            vjust = 1.4, max.overlaps = max.overlaps
        ) +
        scale_size_continuous(
            name = "Node size",
            range = c(vertex_size / 2, vertex_size * 1.5),
            guide = guide_legend(override.aes = list(shape = 19))
        ) +
        scale_fill_identity(
            name = "Cluster",
            guide = guide_legend(
                override.aes = list(shape = 21, size = vertex_size * 0.5),
                nrow = 2, order = 1
            ),
            breaks = unname(palette_vals),
            labels = names(palette_vals)
        ) +
        coord_fixed() + theme_void(base_size = 12) +
        theme(
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 15),
            plot.title = element_text(size = 15, hjust = 0.5),
            legend.position = "bottom",
            legend.box = "vertical",
            legend.box.just = "center",
            legend.spacing.y = unit(0.2, "cm"),
            legend.box.margin = margin(t = 5, b = 5),
            plot.caption = element_text(
                hjust = 0, size = 9,
                margin = margin(t = 6)
            )
        ) +
        labs(title = title, caption = qc_txt) +
        guides(
            fill        = guide_legend(order = 1),
            size        = guide_legend(order = 2),
            edge_width  = guide_legend(order = 3, direction = "horizontal", nrow = 1),
            linetype    = if (show_sign) guide_legend(order = 4, nrow = 1) else "none"
        )

    invisible(list(
        graph       = g,
        pagerank    = if (exists("pr")) pr else NULL,
        cross_edges = cross_edges,
        plot        = p
    ))
}

#' @title Compare cluster assignments produced by multiple CI³ parameter sets
#'
#' @description
#'   Generates a dot plot that juxtaposes gene–cluster memberships obtained
#'   under different CI³ thresholds or weighting options.  The function takes a
#'   set of membership columns already stored in \code{scope_obj@meta.data},
#'   reshapes them into long format, applies an adaptive discrete palette, and
#'   saves the figure to disk.
#'
#' @param scope_obj   A \code{scope_object} whose \code{@meta.data} contains the
#'                   membership columns listed in \code{method_cols}.
#' @param method_cols Character vector of column names, ordered from most
#'                   stringent to least stringent. **Each column must be a
#'                   factor or vector of cluster labels.**
#' @param method_labels Optional character vector of display labels for the
#'                   rows; defaults to \code{method_cols}.
#' @param point_size Size of the dot glyphs (default 3).
#' @param palette    Character vector of colours or a palette function.  When
#'                   the number of clusters exceeds 12, a \code{hue_pal()} ramp
#'                   is used automatically.
#' @param output_file Full path to the PNG file. When \code{NULL} the plot is
#'                   returned but not saved.
#' @param width,height,dpi Dimensions (in inches) and resolution for
#'                   \code{ggsave()}.
#'
#' @return A \code{ggplot} object invisibly.
#'
#' @examples
#' \dontrun{
#' p <- plotClusterComparison(
#'     scope_obj,
#'     method_cols = c(
#'         "All_res0.1_grid50",
#'         "q80_res0.1_grid50",
#'         "q80_res0.1_grid50_log1p",
#'         "q80_res0.1_grid50_log1p_pie0.8",
#'         "q80_res0.1_grid50_log1p_pie0.8_outCI95",
#'         "q80_res0.1_grid50_log1p_pie0.8_outCI95_cmh"
#'     ),
#'     method_labels = c(
#'         "All (Leiden)",
#'         "q80 (Leiden)",
#'         "q80 log1p (Leiden)",
#'         "q80 log1p pie0.8 (Leiden)",
#'         "q80 log1p pie0.8 outCI95 (Leiden)",
#'         "q80 log1p pie0.8 outCI95 CMH (Leiden)"
#'     ),
#'     output_file = "cluster_comparison.png"
#' )
#' }
#' @export
plotClusterComparison <- function(scope_obj,
                                  method_cols,
                                  method_labels = method_cols,
                                  point_size = 3,
                                  palette = RColorBrewer::brewer.pal) {
    stopifnot(length(method_cols) == length(method_labels))

    df <- scope_obj@meta.data |>
        tibble::rownames_to_column("gene") |>
        dplyr::select(gene, tidyselect::all_of(method_cols))

    # keep genes with at least one non‑NA assignment
    df <- df |>
        dplyr::filter(dplyr::if_any(-gene, ~ !is.na(.x)))

    # order genes by progressively less stringent columns
    df <- df |>
        dplyr::arrange(dplyr::across(dplyr::all_of(method_cols), ~ dplyr::desc(.x)))

    long_df <- df |>
        tidyr::pivot_longer(
            cols = tidyselect::all_of(method_cols),
            names_to = "method",
            values_to = "cluster"
        ) |>
        tidyr::drop_na(cluster) |>
        dplyr::mutate(
            method  = factor(method, levels = method_cols, labels = method_labels),
            gene    = factor(gene, levels = df$gene),
            cluster = factor(cluster)
        )

    n_col <- nlevels(long_df$cluster)
    pal <- if (is.function(palette)) {
        if (n_col <= 12) palette(n_col, "Set3") else scales::hue_pal()(n_col)
    } else {
        if (length(palette) < n_col) {
            grDevices::colorRampPalette(palette)(n_col)
        } else {
            palette[seq_len(n_col)]
        }
    }

    p <- ggplot2::ggplot(
        long_df,
        ggplot2::aes(x = gene, y = method, fill = cluster)
    ) +
        ggplot2::geom_point(shape = 21, size = point_size, colour = "black") +
        ggplot2::scale_fill_manual(values = pal) +
        ggplot2::labs(x = NULL, y = "Group", fill = "Cluster") +
        ggplot2::theme_minimal(base_size = 8) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 90, vjust = .5, hjust = 1,
                size = 5
            ),
            panel.title = element_text(size = 10, colour = "black"),
            panel.background = element_rect(fill = "#c0c0c0", colour = NA),
            panel.border = element_rect(colour = "black", fill = NA, size = .5),
            axis.line = element_line(colour = "black", size = .3),
            axis.ticks = element_line(colour = "black", size = .3),
            axis.text = element_text(size = 8, colour = "black"),
            axis.title = element_text(size = 9, colour = "black"),
            legend.position = "none"
        )

    p
}

.order_levels_numeric <- function(lev) {
    lev <- as.character(lev)
    num <- suppressWarnings(as.numeric(lev))
    if (all(!is.na(num))) return(lev[order(num)])
    rx <- regmatches(lev, regexpr("-?\\d+\\.?\\d*", lev))
    num2 <- suppressWarnings(as.numeric(rx))
    ord <- order(is.na(num2), num2, lev)
    lev[ord]
}

.make_pal_map <- function(levels_vec, cluster_palette, tip_palette) {
    levels_vec <- as.character(levels_vec)
    n <- length(levels_vec)
    if (is.null(tip_palette)) {
        if (is.null(cluster_palette)) {
            pal <- grDevices::rainbow(n)
            return(stats::setNames(pal, levels_vec))
        } else if (is.null(names(cluster_palette))) {
            pal <- grDevices::colorRampPalette(cluster_palette)(n)
            return(stats::setNames(pal, levels_vec))
        } else {
            pal <- cluster_palette
            miss <- setdiff(levels_vec, names(pal))
            if (length(miss)) {
                pal <- c(pal, stats::setNames(grDevices::colorRampPalette(cluster_palette)(length(miss)), miss))
            }
            return(pal[levels_vec])
        }
    } else if (is.function(tip_palette)) {
        pal <- tryCatch(tip_palette(n, "Set3"), error = function(...) NULL)
        if (is.null(pal)) pal <- grDevices::hcl.colors(n, palette = "Dynamic")
        return(stats::setNames(pal, levels_vec))
    } else {
        pal <- if (length(tip_palette) < n) grDevices::colorRampPalette(tip_palette)(n) else tip_palette[seq_len(n)]
        return(stats::setNames(pal, levels_vec))
    }
}

.format_label_display <- function(x) {
    gsub("_sub", ".", x, fixed = TRUE)
}

.draw_wrapped_legend <- function(
    x_start, y_start, labels, pch, bg, border,
    point_cex, text_cex, x_range, y_range,
    max_width, dot_gap_factor = 0.6, item_gap_factor = 0.5,
    row_spacing_factor = 1.2) {

    if (!length(labels)) return(invisible(NULL))

    width_user <- function(txt, cex) {
        strwidth(txt, cex = cex, units = "figure") * x_range
    }
    height_user <- function(txt, cex) {
        strheight(txt, cex = cex, units = "figure") * y_range
    }

    text_widths <- width_user(labels, text_cex)
    symbol_width <- width_user("M", point_cex) * 0.6
    dot_gap <- width_user("M", text_cex) * dot_gap_factor
    item_gap <- width_user(" ", text_cex) * item_gap_factor
    row_spacing <- height_user("M", text_cex) * row_spacing_factor

    if (!is.finite(max_width) || max_width <= 0) {
        max_width <- sum(symbol_width + dot_gap + text_widths) + item_gap * (length(labels) - 1)
    }
    min_required <- min(symbol_width + dot_gap + text_widths)
    if (max_width < min_required) {
        max_width <- min_required * 1.1
    }

    rows <- list()
    current <- integer(0)
    width_current <- 0
    for (i in seq_along(labels)) {
        item_width <- symbol_width + dot_gap + text_widths[i]
        item_total <- if (length(current)) item_width + item_gap else item_width
        if (length(current) && width_current + item_total > max_width) {
            rows[[length(rows) + 1]] <- current
            current <- i
            width_current <- item_width
        } else {
            current <- c(current, i)
            width_current <- width_current + item_total
        }
    }
    if (length(current)) rows[[length(rows) + 1]] <- current

    y_current <- y_start
    for (row_idx in seq_along(rows)) {
        idxs <- rows[[row_idx]]
        x_pos <- x_start
        for (j in idxs) {
            graphics::points(x_pos, y_current, pch = pch[j], bg = bg[j], col = border[j], cex = point_cex)
            x_pos <- x_pos + symbol_width
            graphics::text(x_pos + dot_gap, y_current, labels[j], adj = c(0, 0.5), cex = text_cex)
            x_pos <- x_pos + dot_gap + text_widths[j] + item_gap
        }
        y_current <- y_current - row_spacing
    }

    invisible(NULL)
}

#' @title Dendrogram of the plotDendroNetwork skeleton (plotRWDendrogram8-style)
#'
#' @description
#'   Builds the same MST/forest skeleton as in `plotDendroNetwork` (cluster-internal
#'   MST + cluster-level MST after optional PageRank reweighting) from the
#'   specified consensus graph, then renders a dendrogram in the identical
#'   visual style as `plotRWDendrogram8` with optional dual tip-dot rows,
#'   adjustable label offsets, and OLO leaf-order optimisation.
#'
#' @inheritParams plotRWDendrogram8
#' @param graph_slot Slot in the LeeStats layer that stores the base graph
#'                   (e.g. the consensus graph). Defaults to `"g_consensus"`.
#' @param title Custom title for the rendered dendrogram. Defaults to
#'   `"Graph-weighted Dendrogram"`.
#'
#' @return Invisibly returns a list containing `dend`, `hclust`, `dist`,
#'   `genes`, `cluster_map`, `tree_graph` (MST skeleton), `graph` (full graph
#'   on the selected genes), and `PageRank` if applied.
#'
#' @importFrom igraph V E as_data_frame induced_subgraph distances ecount vcount edge_attr_names page_rank ends
#' @importFrom stats hclust as.dendrogram quantile runmed
#' @importFrom graphics plot points text par
#' @export
plotDendro <- function(
    scope_obj,
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    graph_slot = "g_consensus",
    cluster_name,
    cluster_ids = NULL,
    IDelta_col_name = NULL,
    linkage = "average",
    plot_dend = TRUE,
    weight_low_cut = 0,
    damping = 0.85,
    cluster_palette = c(
        "#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33",
        "#A65628", "#984EA3", "#66C2A5", "#FC8D62", "#8DA0CB",
        "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"
    ),
    tip_dots = TRUE,
    tip_point_cex = 1.2,
    tip_palette = NULL,
    cluster_name2 = NULL,
    tip_palette2 = NULL,
    tip_row_offset2 = -0.4,
    tip_label_offset = -0.6,
    tip_label_cex = 0.6,
    tip_label_adj = 1,
    tip_label_srt = 90,
    tip_label_col = "black",
    leaf_order = c("OLO", "none"),
    length_mode = c("neg_log", "inverse", "inverse_sqrt"),
    weight_normalize = TRUE,
    weight_clip_quantile = 0.05,
    height_rescale = c("q95", "max", "none"),
    height_power = 1,
    height_scale = 1.2,
    distance_on = c("tree", "graph"),
    enforce_cluster_contiguity = TRUE,
    distance_smooth_power = 1,
    tip_shape2 = c("square", "circle", "diamond", "triangle"),
    tip_row1_label = NULL,
    tip_row2_label = NULL,
    tip_label_indent = 0.01,
    legend_inline = FALSE,
    legend_files = NULL,
    compose_outfile = NULL,
    compose_width = 2400,
    compose_height = 1600,
    compose_res = 200,
    legend_ncol1 = 1,
    legend_ncol2 = 1,
    title = "Graph-weighted Dendrogram") {

    leaf_order <- match.arg(leaf_order)
    length_mode <- match.arg(length_mode)
    height_rescale <- match.arg(height_rescale)
    distance_on <- match.arg(distance_on)
    tip_shape2 <- match.arg(tip_shape2)
    if (is.null(tip_row1_label)) tip_row1_label <- cluster_name
    if (is.null(tip_row2_label) && !is.null(cluster_name2)) tip_row2_label <- cluster_name2
    if (!exists("tip_label_indent", inherits = FALSE) || is.null(tip_label_indent)) tip_label_indent <- 0.01
    if (!exists("legend_inline", inherits = FALSE) || is.null(legend_inline)) legend_inline <- FALSE
    if (!is.null(distance_smooth_power) && distance_smooth_power <= 0) {
        stop("distance_smooth_power must be > 0 or NULL")
    }

    ## ---------- 0. Resolve grid layer & LeeStats ----------
    g_layer <- .selectGridLayer(scope_obj, grid_name)
    if (is.null(grid_name)) {
        grid_name <- names(scope_obj@grid)[
            vapply(scope_obj@grid, identical, logical(1), g_layer)
        ]
    }

    ## LeeStats object
    leeStat <- if (!is.null(scope_obj@stats[[grid_name]]) &&
        !is.null(scope_obj@stats[[grid_name]][[lee_stats_layer]])) {
        scope_obj@stats[[grid_name]][[lee_stats_layer]]
    } else if (!is.null(g_layer[[lee_stats_layer]])) {
        g_layer[[lee_stats_layer]]
    } else {
        stop("Cannot find layer '", lee_stats_layer, "' in grid '", grid_name, "'.")
    }

    ## ---------- 1. Extract base graph and subset to target genes ----------
    g_base <- leeStat[[graph_slot]]
    if (is.null(g_base) || !inherits(g_base, "igraph")) {
        stop("Graph slot '", graph_slot, "' not found or not an igraph in LeeStats layer.")
    }

    memb_all <- as.character(scope_obj@meta.data[[cluster_name]])
    if (is.null(cluster_ids)) cluster_ids <- sort(unique(na.omit(memb_all)))
    cluster_ids <- as.character(cluster_ids)
    genes_all <- rownames(scope_obj@meta.data)[memb_all %in% cluster_ids]
    if (length(genes_all) < 2) stop("Need at least two genes in the selected clusters.")

    g <- igraph::induced_subgraph(g_base, vids = intersect(igraph::V(g_base)$name, genes_all))
    if (igraph::vcount(g) < 2) stop("Too few vertices in induced subgraph.")

    ## Ensure edges have usable numeric weights (fallback to |L| if missing/invalid)
    Lmat <- .getLeeMatrix(scope_obj,
        grid_name = grid_name,
        lee_layer = lee_stats_layer
    )
    ed_l <- igraph::as_edgelist(g, names = TRUE)
    if (nrow(ed_l) > 0) {
        w_lmat <- mapply(function(a, b) abs(Lmat[a, b]), ed_l[, 1], ed_l[, 2], USE.NAMES = FALSE)
        if (!"weight" %in% igraph::edge_attr_names(g)) {
            igraph::E(g)$weight <- w_lmat
        } else {
            w_now <- igraph::E(g)$weight
            bad <- is.na(w_now) | !is.finite(w_now) | w_now <= 0
            if (any(bad)) w_now[bad] <- w_lmat[bad]
            igraph::E(g)$weight <- w_now
        }
        # Final NA guard
        igraph::E(g)$weight[is.na(igraph::E(g)$weight) | !is.finite(igraph::E(g)$weight)] <- 0
    }

    ## ---------- 2. Δ-PageRank reweighting (same as plotDendroNetwork) ----------
    pr <- NULL
    if (!is.null(IDelta_col_name)) {
        delta <- scope_obj@meta.data[igraph::V(g)$name, IDelta_col_name, drop = TRUE]
        delta[is.na(delta)] <- stats::median(delta, na.rm = TRUE)
        q <- stats::quantile(delta, c(.1, .9))
        delta <- pmax(pmin(delta, q[2]), q[1])
        pers <- exp(delta - max(delta)); pers <- pers / sum(pers)
        pr_tmp <- igraph::page_rank(
            g,
            personalized = pers,
            damping = damping,
            weights = igraph::E(g)$weight,
            directed = FALSE
        )$vector
        pr <- pr_tmp
        ed <- igraph::as_data_frame(g, "edges")
        w_new <- igraph::E(g)$weight * ((pr[ed$from] + pr[ed$to]) / 2)
        w_new[w_new <= weight_low_cut] <- 0
        igraph::E(g)$weight <- w_new
    }

    ## ---------- 3. Skeleton: intra-cluster MST + inter-cluster MST (same as plotDendroNetwork) ----------
    Vnames <- igraph::V(g)$name
    clu <- as.character(scope_obj@meta.data[Vnames, cluster_name, drop = TRUE])
    keep <- !is.na(clu)
    g <- igraph::induced_subgraph(g, Vnames[keep])
    Vnames <- igraph::V(g)$name
    clu <- clu[keep]

    # tag edge ids to track selection
    igraph::E(g)$eid <- seq_len(igraph::ecount(g))

    # 3.1 Intra-cluster MST: use max weight via length 1/(w + eps)
    keep_eid <- integer(0)
    for (cl in unique(clu)) {
        vsub <- Vnames[clu == cl]
        if (length(vsub) < 2) next
        g_sub <- igraph::induced_subgraph(g, vsub)
        if (igraph::ecount(g_sub) > 0) {
            mst_sub <- igraph::mst(g_sub, weights = 1 / (igraph::E(g_sub)$weight + 1e-9))
            keep_eid <- c(keep_eid, igraph::edge_attr(mst_sub, "eid"))
        }
    }

    # 3.2 Inter-cluster MST: aggregate by each pair's max-weight edge
    if (length(unique(clu)) > 1) {
        ed_tab <- igraph::as_data_frame(g, "edges"); ed_tab$eid <- igraph::E(g)$eid
        ed_tab$cl1 <- clu[ed_tab$from]; ed_tab$cl2 <- clu[ed_tab$to]
        inter <- ed_tab[ed_tab$cl1 != ed_tab$cl2, , drop = FALSE]
        if (nrow(inter) > 0) {
            inter$pair <- ifelse(inter$cl1 < inter$cl2,
                paste(inter$cl1, inter$cl2, sep = "|"), paste(inter$cl2, inter$cl1, sep = "|"))
            # drop rows with NA weights to avoid aggregate() zero-row error
            inter <- inter[!is.na(inter$weight) & is.finite(inter$weight), , drop = FALSE]
            if (nrow(inter) == 0) {
                agg <- inter
            } else {
                agg <- stats::aggregate(weight ~ pair + cl1 + cl2, data = inter, max)
            }
            g_clu <- igraph::graph_from_data_frame(
                agg[, c("cl1", "cl2", "weight")], directed = FALSE, vertices = unique(clu)
            )
            if (igraph::ecount(g_clu) > 0) {
                mstc <- igraph::mst(g_clu, weights = 1 / (igraph::E(g_clu)$weight + 1e-9))
                ep <- igraph::ends(mstc, igraph::E(mstc))
                keep_pairs <- ifelse(ep[, 1] < ep[, 2], paste(ep[, 1], ep[, 2], sep = "|"), paste(ep[, 2], ep[, 1], sep = "|"))
                inter_MST <- inter[inter$pair %in% keep_pairs, ]
                keep_eid <- unique(c(keep_eid, inter_MST$eid))
            }
        }
    }

    # 3.3 Keep skeleton edges only; drop isolates
    g_tree <- g
    if (length(keep_eid)) {
        g_tree <- igraph::delete_edges(g_tree, igraph::E(g_tree)[!eid %in% keep_eid])
    }
    g_tree <- igraph::delete_vertices(g_tree, which(igraph::degree(g_tree) == 0))
    if (igraph::vcount(g_tree) < 2) stop("MST-based tree has fewer than 2 vertices.")

    ## ---------- 4. Distance matrix (tree-based to merge within clusters first) ----------
    genes <- igraph::V(g_tree)$name
    # keep full-graph view for return, regardless of distance_on
    g_full <- igraph::induced_subgraph(g, vids = genes)
    if (distance_on == "graph") {
        ew <- igraph::E(g_full)$weight
        if (isTRUE(weight_normalize)) {
            ew <- (ew - min(ew, na.rm = TRUE)) / max(1e-12, diff(range(ew, na.rm = TRUE)))
        }
        if (!is.null(weight_clip_quantile) && weight_clip_quantile > 0) {
            ql <- stats::quantile(ew, probs = weight_clip_quantile, na.rm = TRUE)
            qu <- stats::quantile(ew, probs = 1 - weight_clip_quantile, na.rm = TRUE)
            ew <- pmin(pmax(ew, ql), qu)
        }
        ew <- ifelse(ew <= 0 | is.na(ew), 1e-6, ew)
        elen <- switch(length_mode,
            neg_log      = -log(pmax(ew, 1e-9)),
            inverse      = 1 / pmax(ew, 1e-9),
            inverse_sqrt = 1 / sqrt(pmax(ew, 1e-9))
        )
        elen <- (elen ^ height_power) * height_scale
        igraph::E(g_full)$length <- elen
        Dm <- igraph::distances(g_full, v = igraph::V(g_full), to = igraph::V(g_full), weights = igraph::E(g_full)$length)
        rownames(Dm) <- igraph::V(g_full)$name
        colnames(Dm) <- igraph::V(g_full)$name
        Dm <- Dm[genes, genes, drop = FALSE]
    } else {
        # compute lengths on the skeleton tree and use its geodesics
        ew_t <- igraph::E(g_tree)$weight
        if (isTRUE(weight_normalize)) {
            ew_t <- (ew_t - min(ew_t, na.rm = TRUE)) / max(1e-12, diff(range(ew_t, na.rm = TRUE)))
        }
        if (!is.null(weight_clip_quantile) && weight_clip_quantile > 0) {
            ql <- stats::quantile(ew_t, probs = weight_clip_quantile, na.rm = TRUE)
            qu <- stats::quantile(ew_t, probs = 1 - weight_clip_quantile, na.rm = TRUE)
            ew_t <- pmin(pmax(ew_t, ql), qu)
        }
        ew_t <- ifelse(ew_t <= 0 | is.na(ew_t), 1e-6, ew_t)
        elen_t <- switch(length_mode,
            neg_log      = -log(pmax(ew_t, 1e-9)),
            inverse      = 1 / pmax(ew_t, 1e-9),
            inverse_sqrt = 1 / sqrt(pmax(ew_t, 1e-9))
        )
        elen_t <- (elen_t ^ height_power) * height_scale
        igraph::E(g_tree)$length <- elen_t
        Dm <- igraph::distances(g_tree, v = genes, to = genes, weights = igraph::E(g_tree)$length)
        rownames(Dm) <- genes
        colnames(Dm) <- genes
    }

    # optional smoothing of the distance matrix
    if (!is.null(distance_smooth_power)) {
        k <- max(3L, as.integer(round(ncol(Dm) * (distance_smooth_power / 10))))
        if ((k %% 2L) == 0L) k <- k + 1L
        if (k > 3) {
            ut <- upper.tri(Dm, diag = FALSE)
            dv <- as.numeric(Dm[ut])
            ord <- order(dv); x <- dv[ord]
            xs <- stats::runmed(x, k = k, endrule = "keep")
            # guard NA from runmed at boundaries
            xs[!is.finite(xs)] <- 0
            dv[ord] <- pmax(xs, 0)
            Dm2 <- Dm; Dm2[ut] <- dv
            Dm2 <- Dm2 + t(Dm2); diag(Dm2) <- 0
            Dm <- Dm2
        }
    }

    # fill disconnected/invalid entries with a large finite distance
    if (any(!is.finite(Dm))) {
        maxf <- suppressWarnings(max(Dm[is.finite(Dm)], na.rm = TRUE))
        if (!is.finite(maxf) || is.na(maxf)) maxf <- 1
        Dm[!is.finite(Dm)] <- maxf * 1.2 + 1
    }
    diag(Dm) <- 0

        # Optionally force within-cluster merges to occur first
    if (isTRUE(enforce_cluster_contiguity)) {
        clv <- memb_all[genes]
        within_mask <- outer(clv, clv, "==")
        ut <- upper.tri(Dm)
        wvals <- Dm[ut & within_mask]
        bvals <- Dm[ut & !within_mask]
        w_max <- suppressWarnings(max(wvals[wvals > 0 & is.finite(wvals)], na.rm = TRUE))
        b_min <- suppressWarnings(min(bvals[bvals > 0 & is.finite(bvals)], na.rm = TRUE))
        if (is.finite(w_max) && is.finite(b_min) && b_min <= w_max) {
            off <- (w_max - b_min) + max(1e-6, 0.05 * w_max)
            Dm[!within_mask] <- Dm[!within_mask] + off
        }
    }

    # rescale heights if requested
    h <- if (height_rescale == "q95") {
        stats::quantile(Dm[upper.tri(Dm)], 0.95, na.rm = TRUE)
    } else if (height_rescale == "max") {
        max(Dm, na.rm = TRUE)
    } else {
        NA_real_
    }
    if (is.finite(h) && h > 0) Dm <- Dm / h

    hc <- stats::hclust(stats::as.dist(Dm), method = linkage)
    dend <- stats::as.dendrogram(hc)

    ## ---------- Optional OLO ----------
    if (leaf_order == "OLO") {
        ok_ser <- requireNamespace("seriation", quietly = TRUE)
        ok_den <- requireNamespace("dendextend", quietly = TRUE)
        if (ok_ser && ok_den) {
            ord <- seriation::get_order(seriation::seriate(stats::as.dist(Dm), method = "OLO"))
            desired <- rownames(Dm)[ord]
            dend <- dendextend::rotate(dend, order = desired)
        } else if (requireNamespace("gclus", quietly = TRUE)) {
            hc <- gclus::reorder.hclust(hc, dist = stats::as.dist(Dm))
            dend <- stats::as.dendrogram(hc)
        } else if (requireNamespace("dendsort", quietly = TRUE)) {
            dend <- dendsort::dendsort(dend)
        } else if (!ok_ser || !ok_den) {
            warning("Leaf-order optimization requested but required packages are missing; install one of: seriation + dendextend (preferred), gclus, or dendsort.")
        }
    }

    labs <- labels(dend)
    x <- seq_along(labs)

    tip1 <- list(show = isTRUE(tip_dots) && length(labs) > 0,
                 title = tip_row1_label)
    if (tip1$show) {
        cl_raw1 <- scope_obj@meta.data[labs, cluster_name, drop = TRUE]
        cl_chr1 <- as.character(cl_raw1)
        lev1 <- .order_levels_numeric(unique(cl_chr1))
        pal_map1 <- .make_pal_map(lev1, cluster_palette, tip_palette)
        cols1 <- unname(pal_map1[cl_chr1]); cols1[is.na(cols1)] <- "gray80"
        tip1$cols <- cols1
        tip1$legend_labels <- .format_label_display(lev1)
        tip1$legend_colors <- unname(pal_map1[lev1])
        tip1$pch <- rep(21, length(lev1))
        tip1$lev_order <- lev1
    } else {
        tip1$legend_labels <- character(0)
        tip1$legend_colors <- character(0)
        tip1$pch <- numeric(0)
    }
    if (is.null(tip1$title)) tip1$title <- ""

    tip2 <- list(show = !is.null(cluster_name2) && length(labs) > 0,
                 title = if (is.null(tip_row2_label)) "" else tip_row2_label)
    if (tip2$show) {
        cl_raw2 <- scope_obj@meta.data[labs, cluster_name2, drop = TRUE]
        cl_chr2 <- as.character(cl_raw2)
        lev2 <- .order_levels_numeric(unique(cl_chr2))
        pal_map2 <- .make_pal_map(lev2, cluster_palette, tip_palette2)
        cols2 <- unname(pal_map2[cl_chr2]); cols2[is.na(cols2)] <- "gray80"
        tip2$cols <- cols2
        tip2$legend_labels <- .format_label_display(lev2)
        tip2$legend_colors <- unname(pal_map2[lev2])
        tip2$pch <- rep(switch(tip_shape2, square = 22, circle = 21, diamond = 23, triangle = 24), length(lev2))
        tip2$lev_order <- lev2
    } else {
        tip2$legend_labels <- character(0)
        tip2$legend_colors <- character(0)
        tip2$pch <- numeric(0)
    }
    if (is.null(tip2$title)) tip2$title <- ""

    render_dend_panel <- function(draw_legends = FALSE, preserve_par = TRUE) {
        if (preserve_par) {
            op <- graphics::par(no.readonly = TRUE)
            on.exit(graphics::par(op), add = TRUE)
        } else {
            op <- graphics::par(c("mar", "mgp", "xpd"))
            on.exit(do.call(graphics::par, op), add = TRUE)
        }

        graphics::par(mar = c(3.9, 3.3, 3.1, 0.2), mgp = c(1.5, 0.25, 0), xpd = NA)

        max_h <- if (length(hc$height)) max(hc$height, na.rm = TRUE) else 1
        if (!is.finite(max_h) || is.na(max_h) || max_h <= 0) max_h <- 1
        y_gap <- 0.06 * max_h
        min_off <- min(0, tip_label_offset, if (tip2$show) tip_row_offset2 else 0)
        pad_off <- 0.1
        ylim_use <- c(y_gap * (min_off - pad_off), max_h)

        graphics::plot(dend,
            main = title,
            ylab = "Weighted distance",
            ylim = ylim_use,
            leaflab = "none"
        )

        labs_local <- labels(dend)
        x_local <- seq_along(labs_local)
        usr <- graphics::par("usr")
        x_range <- usr[2] - usr[1]
        y_range <- usr[4] - usr[3]

        if (tip1$show) {
            y1 <- rep(0, length(labs_local))
            graphics::points(x_local, y1, pch = 21, bg = tip1$cols, col = "black", cex = tip_point_cex)
            if (draw_legends && length(tip1$legend_labels)) {
                label_x1 <- usr[1] + tip_label_indent * x_range
                legend_y1 <- mean(y1) - 0.12 * y_gap
                legend_margin <- 0.02 * x_range
                max_width1 <- max(legend_margin, usr[2] - label_x1 - legend_margin)
                .draw_wrapped_legend(
                    x_start = label_x1,
                    y_start = legend_y1,
                    labels = tip1$legend_labels,
                    pch = tip1$pch,
                    bg = tip1$legend_colors,
                    border = rep("black", length(tip1$legend_labels)),
                    point_cex = tip_point_cex,
                    text_cex = 0.8,
                    x_range = x_range,
                    y_range = y_range,
                    max_width = max_width1,
                    dot_gap_factor = 0.6,
                    item_gap_factor = 0.4,
                    row_spacing_factor = 1.1
                )
            }
        }

        if (tip2$show) {
            y2 <- rep(tip_row_offset2 * y_gap, length(labs_local))
            pch2_sym <- if (length(tip2$pch)) tip2$pch[1] else 21
            graphics::points(x_local, y2, pch = pch2_sym, bg = tip2$cols, col = "black", cex = tip_point_cex)
            if (draw_legends && length(tip2$legend_labels)) {
                label_x2 <- usr[1] + tip_label_indent * x_range
                legend_y2 <- mean(y2) - 0.28 * y_gap
                legend_margin <- 0.02 * x_range
                max_width2 <- max(legend_margin, usr[2] - label_x2 - legend_margin)
                .draw_wrapped_legend(
                    x_start = label_x2,
                    y_start = legend_y2,
                    labels = tip2$legend_labels,
                    pch = rep(pch2_sym, length(tip2$legend_labels)),
                    bg = tip2$legend_colors,
                    border = rep("black", length(tip2$legend_labels)),
                    point_cex = tip_point_cex,
                    text_cex = 0.8,
                    x_range = x_range,
                    y_range = y_range,
                    max_width = max_width2,
                    dot_gap_factor = 0.6,
                    item_gap_factor = 0.4,
                    row_spacing_factor = 1.1
                )
            }
        }

        y_lab <- rep(tip_label_offset * y_gap, length(labs_local))
        graphics::text(x_local, y_lab, labels = labs_local,
            srt = tip_label_srt, xpd = NA,
            cex = tip_label_cex, adj = tip_label_adj,
            col = tip_label_col)
    }

    render_legend_panel <- function() {
        graphics::par(mar = c(0.5, 0.4, 0.5, 0.2))
        graphics::plot.new()
        usr <- graphics::par("usr")
        x_span <- usr[2] - usr[1]
        y_span <- usr[4] - usr[3]
        x_left <- usr[1] + 0.08 * x_span
        legend_width <- 0.28 * x_span
        legend_gap <- legend_width
        x_right <- x_left + 2 * legend_width + legend_gap
        if (x_right > usr[2] - 0.02 * x_span) {
            overflow <- x_right - (usr[2] - 0.02 * x_span)
            x_left <- x_left - overflow
            x_right <- x_right - overflow
        }
        y_top <- usr[4] - 0.15 * y_span
        if (tip1$show && length(tip1$legend_labels)) {
            graphics::legend(x_left, y_top, legend = tip1$legend_labels,
                title = tip1$title, pch = tip1$pch,
                pt.bg = tip1$legend_colors, col = "black",
                pt.cex = tip_point_cex, cex = 0.9, bty = "n",
                ncol = legend_ncol1, xjust = 0, yjust = 1,
                x.intersp = 0.6, y.intersp = 0.8)
        }
        if (tip2$show && length(tip2$legend_labels)) {
            pch2_sym <- if (length(tip2$pch)) tip2$pch[1] else 21
            graphics::legend(x_right, y_top, legend = tip2$legend_labels,
                title = tip2$title, pch = rep(pch2_sym, length(tip2$legend_labels)),
                pt.bg = tip2$legend_colors, col = "black",
                pt.cex = tip_point_cex, cex = 0.9, bty = "n",
                ncol = legend_ncol2, xjust = 1, yjust = 1,
                x.intersp = 0.6, y.intersp = 0.8)
        }
    }
    if (!is.null(compose_outfile)) {
        grDevices::png(filename = compose_outfile, width = compose_width, height = compose_height, res = compose_res)
        on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
        graphics::layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1), heights = 1)
        render_dend_panel(draw_legends = FALSE, preserve_par = FALSE)
        render_legend_panel()
        graphics::layout(1)
    }

    if (isTRUE(plot_dend)) {
        if (legend_inline) {
            render_dend_panel(draw_legends = TRUE)
        } else {
            op_layout <- graphics::par(no.readonly = TRUE)
            graphics::layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1), heights = 1)
            render_dend_panel(draw_legends = FALSE, preserve_par = FALSE)
            render_legend_panel()
            graphics::layout(1)
            graphics::par(op_layout)
        }
    }

    if (!is.null(legend_files)) {
        if (length(legend_files) >= 1 && !is.na(legend_files[1]) && tip1$show && length(tip1$legend_labels)) {
            grDevices::png(filename = legend_files[1], width = 1200, height = 800, res = 150)
            graphics::par(mar = c(0, 0, 0, 0))
            graphics::plot.new()
            graphics::legend(
                x = "center", y = "center", legend = tip1$legend_labels,
                title = tip1$title, pch = tip1$pch,
                pt.bg = tip1$legend_colors, col = "black",
                pt.cex = tip_point_cex, cex = 0.9, bty = "n",
                ncol = legend_ncol1, x.intersp = 0.6, y.intersp = 0.8
            )
            grDevices::dev.off()
        }
        if (length(legend_files) >= 2 && !is.na(legend_files[2]) && tip2$show && length(tip2$legend_labels)) {
            grDevices::png(filename = legend_files[2], width = 1200, height = 800, res = 150)
            graphics::par(mar = c(0, 0, 0, 0))
            graphics::plot.new()
            graphics::legend(
                x = "center", y = "center", legend = tip2$legend_labels,
                title = tip2$title, pch = rep(tip2$pch[1], length(tip2$legend_labels)),
                pt.bg = tip2$legend_colors, col = "black",
                pt.cex = tip_point_cex, cex = 0.9, bty = "n",
                ncol = legend_ncol2, x.intersp = 0.6, y.intersp = 0.8
            )
            grDevices::dev.off()
        }
    }

    invisible(list(
        dend        = dend,
        hclust      = hc,
        dist        = stats::as.dist(Dm),
        genes       = genes,
        cluster_map = setNames(memb_all[genes], genes),
        tree_graph  = g_tree,
        graph       = g_full,
        PageRank    = if (!is.null(pr)) pr[names(pr) %in% genes] else NULL,
        legend_row1 = list(
            labels = tip1$legend_labels,
            colors = tip1$legend_colors,
            pch = tip1$pch,
            title = tip1$title
        ),
        legend_row2 = list(
            labels = tip2$legend_labels,
            colors = tip2$legend_colors,
            pch = tip2$pch,
            title = tip2$title
        )
    ))
}
