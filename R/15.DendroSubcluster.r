#' @title Dendrogram-based Branch Subclustering Visualization
#' @description
#' Calls existing \code{plotDendroNetwork()} to build tree/skeleton network, then performs
#' "articulation point branching" subclustering refinement on specified clusters (or all clusters);
#' if no candidates and \code{fallback_community=TRUE} is set,
#' falls back to Louvain/Leiden community splitting. Outputs new network plot with subcluster coloring.
#' @inheritParams plotDendroNetwork
#' @param enable_subbranch Whether to enable branch subclustering refinement (default TRUE)
#' @param cluster_id Optional character vector: clusters to subcluster split only; NULL=all
#' @param include_root Whether articulation point method includes articulation points in each branch
#' @param max_subclusters Maximum number of subclusters to retain per cluster (greedy deduplication by size)
#' @param fallback_community Whether to fall back to community detection when articulation points yield no results
#' @param min_sub_size Minimum size threshold for subclusters in community fallback mode
#' @param community_method Priority order of community detection algorithms
#' @param subbranch_palette Subcluster color palette: 1st color = remaining/unrefined, others cycle for subclusters
#' @param downstream_min_size Downstream subtree size threshold (NULL = adaptive)
#' @param force_split Whether to attempt relaxed acceptance conditions when original community split fails
#' @param main_fraction_cap Maximum community proportion threshold (allows retention of smaller communities even when exceeded)
#' @param core_periph Allow core-periphery splitting
#' @param core_degree_quantile Core degree threshold quantile
#' @param core_min_fraction Minimum core fraction
#' @param degree_gini_threshold Degree Gini threshold that triggers core-periphery splitting
#' @return list(plot, graph, subclusters, subcluster_df, cluster_df, method_subcluster,
#'              base_network, branch_network, params)
#' @examples
#' \dontrun{
#' p <- plotDendroNetworkWithBranches(coordObj,
#'     grid_name = "grid_lenGrid30",
#'     cluster_vec = "modL0.15",
#'     enable_subbranch = TRUE
#' )
#' }
#' @export
plotDendroNetworkWithBranches <- function(
    coordObj,
    ## Base plotDendroNetwork parameters (maintain order for compatibility)
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    gene_subset = NULL,
    L_min = 0,
    drop_isolated = TRUE,
    use_consensus_graph = TRUE,
    graph_slot_name = "g_consensus",
    cluster_vec = NULL,
    IDelta_col_name = NULL,
    damping = 0.85,
    weight_low_cut = 0,
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
    tree_mode = c("rooted", "radial", "forest"),
    ## New subcluster-related parameters
    enable_subbranch = TRUE,
    cluster_id = NULL,
    include_root = TRUE,
    max_subclusters = 10,
    fallback_community = TRUE,
    min_sub_size = 3,
    community_method = c("louvain", "leiden"),
    subbranch_palette = c(
        "#999999", "#D55E00", "#0072B2", "#009E73",
        "#CC79A7", "#F0E442", "#56B4E9", "#E69F00"
    ),
    downstream_min_size = NULL,
    ## Additional control parameters (conservative defaults)
    force_split = TRUE,
    main_fraction_cap = 0.9,
    core_periph = TRUE,
    core_degree_quantile = 0.75,
    core_min_fraction = 0.05,
    degree_gini_threshold = 0.35,
    verbose = TRUE) {
    tree_mode <- match.arg(tree_mode)
    community_method <- match.arg(community_method, several.ok = TRUE)

    # 0. Call base network function
    if (verbose) message("[geneSCOPE] [Base] Constructing basic tree network...")
    base_args <- list(
        coordObj = coordObj,
        grid_name = grid_name,
        lee_stats_layer = lee_stats_layer,
        gene_subset = gene_subset,
        L_min = L_min,
        drop_isolated = drop_isolated,
        use_consensus_graph = use_consensus_graph,
        graph_slot_name = graph_slot_name,
        cluster_vec = cluster_vec,
        IDelta_col_name = IDelta_col_name,
        damping = damping,
        weight_low_cut = weight_low_cut,
        cluster_palette = cluster_palette,
        vertex_size = vertex_size,
        base_edge_mult = base_edge_mult,
        label_cex = label_cex,
        seed = seed,
        hub_factor = hub_factor,
        length_scale = length_scale,
        max.overlaps = max.overlaps,
        hub_border_col = hub_border_col,
        hub_border_size = hub_border_size,
        show_sign = show_sign,
        neg_linetype = neg_linetype,
        neg_legend_lab = neg_legend_lab,
        pos_legend_lab = pos_legend_lab,
        show_qc_caption = show_qc_caption,
        title = title,
        k_top = k_top,
        tree_mode = tree_mode
    )
    base_network <- do.call(plotDendroNetwork, base_args)
    g <- base_network$graph
    if (is.null(g) || !inherits(g, "igraph")) {
        stop("Base network construction failed or did not return igraph object")
    }

    vnames <- igraph::V(g)$name

    # 1. Parse cluster labels (align with base function logic)
    if (verbose) message("[geneSCOPE] [Subbranch] Preparing cluster labels...")
    Vnames <- vnames

    if (!is.null(cluster_vec)) {
        if (length(cluster_vec) == 1) {
            if (is.null(coordObj@meta.data) || !(cluster_vec %in% colnames(coordObj@meta.data))) {
                stop("meta.data does not contain column: ", cluster_vec)
            }
            clu_full <- coordObj@meta.data[[cluster_vec]]
            names(clu_full) <- rownames(coordObj@meta.data)
            clu <- as.character(clu_full[Vnames])
        } else {
            if (is.null(names(cluster_vec))) {
                stop("cluster_vec as vector must have names")
            }
            clu <- as.character(cluster_vec[Vnames])
        }
    } else {
        clu <- rep("C1", length(Vnames))
        if (verbose) message("[geneSCOPE] [Subbranch] No cluster_vec provided, using single cluster C1")
    }
    names(clu) <- Vnames
    cluster_df <- data.frame(gene = Vnames, cluster = clu, stringsAsFactors = FALSE)

    # 2. If subclustering is disabled, return directly
    if (!enable_subbranch) {
        if (verbose) message("[geneSCOPE] [Subbranch] enable_subbranch=FALSE, returning base network")
        return(list(
            plot = base_network$plot,
            graph = g,
            subclusters = list(),
            subcluster_df = data.frame(
                gene = Vnames, cluster = clu,
                subcluster = NA_character_, color = NA_character_
            ),
            cluster_df = cluster_df,
            method_subcluster = "disabled",
            base_network = base_network,
            params = list(enable_subbranch = FALSE)
        ))
    }

    # 3. Define single cluster detection function
    detect_one_cluster <- function(genes_target, cid) {
        tryCatch(
            {
                # Validate gene set
                genes_target <- unique(stats::na.omit(genes_target))
                if (any(!genes_target %in% igraph::V(g)$name)) {
                    genes_target <- intersect(genes_target, igraph::V(g)$name)
                }
                if (!length(genes_target)) {
                    return(list(method = "none", subclusters = list()))
                }

                # Create subgraph
                sg <- tryCatch(
                    igraph::induced_subgraph(g, vids = genes_target),
                    error = function(e) {
                        warning("[geneSCOPE] Subgraph construction failed for cluster ", cid, ": ", e$message)
                        return(NULL)
                    }
                )
                if (is.null(sg) || igraph::vcount(sg) == 0) {
                    return(list(method = "none", subclusters = list()))
                }

                if (verbose) {
                    message(
                        "[geneSCOPE] Cluster ", cid, " subgraph: nodes=",
                        igraph::vcount(sg), " edges=", igraph::ecount(sg)
                    )
                }

                method_used <- "none"
                subclusters <- list()
                vcount_sg <- igraph::vcount(sg)

                # Try articulation point method
                if (vcount_sg >= 3) {
                    tryCatch(
                        {
                            arts_vs <- igraph::articulation_points(sg)
                            if (length(arts_vs)) {
                                arts_ids <- as.integer(igraph::as_ids(arts_vs))
                                arts_ids <- arts_ids[is.finite(arts_ids) & arts_ids >= 1 & arts_ids <= vcount_sg]

                                cand <- list()
                                for (a_id in arts_ids) {
                                    root_name <- igraph::V(sg)$name[a_id]
                                    sg_minus <- tryCatch(igraph::delete_vertices(sg, a_id), error = function(e) NULL)
                                    if (is.null(sg_minus)) next

                                    comps <- tryCatch(igraph::components(sg_minus), error = function(e) NULL)
                                    if (is.null(comps)) next

                                    for (cid2 in seq_len(comps$no)) {
                                        idx_comp <- which(comps$membership == cid2)
                                        if (!length(idx_comp)) next

                                        sg_comp <- igraph::induced_subgraph(sg_minus, vids = idx_comp)
                                        if (!any(igraph::degree(sg_comp) >= 2)) next

                                        genes_comp <- igraph::V(sg_minus)$name[idx_comp]
                                        branch_genes <- if (include_root) {
                                            unique(c(root_name, genes_comp))
                                        } else {
                                            genes_comp
                                        }
                                        cand[[length(cand) + 1]] <- list(
                                            size = length(branch_genes),
                                            genes = branch_genes
                                        )
                                    }
                                }

                                if (length(cand)) {
                                    ord <- order(vapply(cand, `[[`, numeric(1), "size"), decreasing = TRUE)
                                    used <- character(0)
                                    kept <- list()
                                    for (i in ord) {
                                        gs <- cand[[i]]$genes
                                        if (!any(gs %in% used)) {
                                            kept[[length(kept) + 1]] <- cand[[i]]
                                            used <- c(used, gs)
                                            if (length(kept) >= max_subclusters) break
                                        }
                                    }
                                    if (length(kept)) {
                                        method_used <- "articulation"
                                        for (k in seq_along(kept)) {
                                            subclusters[[paste0(cid, "_sub", k)]] <- kept[[k]]$genes
                                        }
                                    }
                                }
                            }
                        },
                        error = function(e) {
                            if (verbose) {
                                message("[geneSCOPE] Articulation method failed for cluster ", cid, ": ", e$message)
                            }
                        }
                    )
                }

                # Fallback to community detection
                if (method_used == "none" && fallback_community) {
                    comm <- NULL
                    for (mtd in community_method) {
                        comm <- tryCatch(
                            switch(mtd,
                                louvain = igraph::cluster_louvain(sg),
                                leiden = igraph::cluster_leiden(sg),
                                NULL
                            ),
                            error = function(e) NULL
                        )
                        if (!is.null(comm) && length(unique(comm$membership)) > 1) break
                    }

                    if (!is.null(comm)) {
                        tab <- table(comm$membership)
                        if (length(tab) >= 2) {
                            main_c <- as.integer(names(tab)[which.max(tab)])
                            cand_ids <- as.integer(names(tab)[names(tab) != main_c & tab >= min_sub_size])
                            if (length(cand_ids)) {
                                method_used <- "community"
                                k <- 1
                                for (sid in cand_ids) {
                                    subclusters[[paste0(cid, "_sub", k)]] <- igraph::V(sg)$name[comm$membership == sid]
                                    k <- k + 1
                                    if (length(subclusters) >= max_subclusters) break
                                }
                            }
                        }
                    }
                }

                list(method = method_used, subclusters = subclusters)
            },
            error = function(e) {
                warning("[geneSCOPE] Error processing cluster ", cid, ": ", e$message)
                list(method = "none", subclusters = list())
            }
        )
    }

    # 4. Process clusters for subclustering
    target_clusters <- if (is.null(cluster_id)) {
        sort(unique(na.omit(clu)))
    } else {
        intersect(unique(clu), cluster_id)
    }

    if (!length(target_clusters)) {
        if (verbose) message("[geneSCOPE] [Subbranch] No available target clusters, skipping subdivision")
        enable_subbranch <- FALSE
    }

    sub_attr <- setNames(rep(NA_character_, length(Vnames)), Vnames)
    subclusters_all <- list()
    methods_seen <- character(0)

    if (enable_subbranch) {
        if (verbose) message("[geneSCOPE] [Subbranch] Processing ", length(target_clusters), " clusters")
        for (cid in target_clusters) {
            genes_target_raw <- names(clu)[clu == cid]
            genes_target <- intersect(unique(stats::na.omit(genes_target_raw)), Vnames)

            if (length(genes_target_raw) != length(genes_target) && verbose) {
                message(
                    "[geneSCOPE] Cluster ", cid, ": Filtered genes ",
                    length(genes_target_raw), " -> ", length(genes_target)
                )
            }

            if (length(genes_target) < 3) {
                if (verbose) message("[geneSCOPE] Cluster ", cid, ": Too few nodes, skipping")
                next
            }

            res_c <- detect_one_cluster(genes_target, cid)
            methods_seen <- c(methods_seen, res_c$method)

            if (length(res_c$subclusters)) {
                for (nm in names(res_c$subclusters)) {
                    gs <- intersect(res_c$subclusters[[nm]], Vnames)
                    if (!length(gs)) next
                    subclusters_all[[nm]] <- gs
                    sub_attr[gs] <- nm
                }
            }

            if (verbose) {
                message(
                    "[geneSCOPE] Cluster ", cid, ": method=", res_c$method,
                    " subclusters=", length(res_c$subclusters)
                )
            }
        }
    }

    method_final <- if (!length(subclusters_all)) {
        if (enable_subbranch) "none" else "disabled"
    } else {
        setdiff(unique(methods_seen), "none")[1]
    }

    igraph::V(g)$subcluster <- sub_attr

    # 5. Build ordered factor levels (parent cluster -> subclusters)
    numeric_sort_levels <- function(x) {
        ux <- unique(x)
        suppressWarnings({
            nx <- as.numeric(ux)
        })
        if (all(!is.na(nx))) as.character(sort(as.numeric(ux))) else sort(ux)
    }

    base_levels <- numeric_sort_levels(clu)

    if (method_final %in% c("none", "disabled")) {
        cluster_df$cluster <- factor(cluster_df$cluster, levels = base_levels)
    } else {
        sub_levels_raw <- sort(unique(na.omit(sub_attr)))
        parent_of <- sub("^([^_]+)_.*$", "\\1", sub_levels_raw)
        level_vec <- character(0)

        for (pv in base_levels) {
            level_vec <- c(level_vec, pv)
            if (pv %in% parent_of) {
                kids <- sub_levels_raw[parent_of == pv]
                level_vec <- c(level_vec, kids)
            }
        }

        orphan <- setdiff(sub_levels_raw, level_vec)
        if (length(orphan)) level_vec <- c(level_vec, orphan)

        cluster_df$cluster <- factor(cluster_df$cluster, levels = base_levels)
    }

    # 6. Build result tables and colors
    subcluster_df <- data.frame(
        gene = Vnames,
        cluster = cluster_df$cluster,
        subcluster = sub_attr,
        stringsAsFactors = FALSE
    )

    if (!(method_final %in% c("none", "disabled"))) {
        sub_levels_only <- level_vec[level_vec %in% subcluster_df$subcluster]
        subcluster_df$subcluster <- factor(subcluster_df$subcluster, levels = sub_levels_only)
    }

    if (method_final %in% c("none", "disabled")) {
        uniq_clu <- sort(unique(na.omit(clu)))
        col_map <- setNames(rep(cluster_palette, length.out = length(uniq_clu)), uniq_clu)
        node_cols <- col_map[clu]
    } else {
        sub_levels <- sort(unique(na.omit(sub_attr)))
        col_map <- setNames(rep(subbranch_palette[1], length(Vnames)), Vnames)
        if (length(sub_levels)) {
            cols_sub <- rep(subbranch_palette[-1], length.out = length(sub_levels))
            for (i in seq_along(sub_levels)) {
                gs <- names(sub_attr)[sub_attr == sub_levels[i]]
                col_map[gs] <- cols_sub[i]
            }
        }
        node_cols <- col_map[Vnames]
    }
    subcluster_df$color <- node_cols

    # 7. Generate final plot
    if (verbose) message("[geneSCOPE] [Subbranch] Generating final plot...")
    branch_network <- NULL

    if (method_final %in% c("none", "disabled")) {
        cluster_vec_base <- factor(clu, levels = base_levels)
        new_args <- base_args
        new_args$cluster_vec <- setNames(cluster_vec_base, names(clu))
        new_args$title <- if (is.null(title)) "Base Network" else title
        branch_network <- do.call(plotDendroNetwork, new_args)
        plt <- branch_network$plot
    } else {
        cluster_vec_sub <- setNames(ifelse(is.na(sub_attr), clu, sub_attr), names(clu))
        cluster_vec_sub <- factor(cluster_vec_sub, levels = level_vec)
        new_args <- base_args
        new_args$cluster_vec <- cluster_vec_sub
        new_args$title <- if (is.null(title)) {
            paste0("Sub-branch: ", method_final)
        } else {
            paste0(title, " | Sub-branch: ", method_final)
        }
        branch_network <- do.call(plotDendroNetwork, new_args)
        plt <- branch_network$plot
    }

    # 8. Return results
    list(
        plot = plt,
        graph = g,
        subclusters = subclusters_all,
        subcluster_df = subcluster_df,
        cluster_df = cluster_df,
        method_subcluster = method_final,
        base_network = base_network,
        branch_network = branch_network,
        params = list(
            enable_subbranch = enable_subbranch,
            cluster_id = cluster_id,
            include_root = include_root,
            max_subclusters = max_subclusters,
            fallback_community = fallback_community,
            min_sub_size = min_sub_size,
            downstream_min_size = downstream_min_size,
            community_method = community_method,
            subbranch_palette = subbranch_palette,
            force_split = force_split,
            main_fraction_cap = main_fraction_cap,
            core_periph = core_periph,
            core_degree_quantile = core_degree_quantile,
            core_min_fraction = core_min_fraction,
            degree_gini_threshold = degree_gini_threshold
        )
    )
}
