## ===== Helper utilities =====
.parse_q <- function(qstr) {
    if (!is.character(qstr) || length(qstr) != 1 ||
        !grepl("^[qQ][0-9]+\\.?[0-9]*$", qstr)) {
        stop("pct string must look like 'q90' / 'q99.5' …")
    }
    val <- as.numeric(sub("^[qQ]", "", qstr)) / 100
    if (is.na(val) || val < 0 || val > 1) {
        stop("pct out of range 0–100")
    }
    val
}

.filter_matrix_by_quantile <- function(mat, pct_min = "q0", pct_max = "q100") {
    stopifnot(is.matrix(mat))
    pmin <- .parse_q(pct_min)
    pmax <- .parse_q(pct_max)
    if (pmin > pmax) stop("pct_min > pct_max")
    vec <- as.vector(mat)
    thrL <- as.numeric(stats::quantile(vec, pmin, na.rm = TRUE))
    thrU <- as.numeric(stats::quantile(vec, pmax, na.rm = TRUE))
    mat[vec < thrL | vec > thrU] <- 0
    Matrix::drop0(mat)
}

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
#' @param coordObj        A \code{CoordObj} containing \code{@grid} and Lee's
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
plotNetworkGenes <- function(
    coordObj,
    ## ---------- 数据层 ----------
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    gene_subset = NULL, # 现在兼容直接给基因向量或 (cluster_col, cluster_num)
    ## ---------- 过滤阈值 ----------
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
    ## ---------- 共识网络 ----------
    use_consensus_graph = TRUE,
    graph_slot_name = "g_consensus",
    ## ---------- FDR 来源新增 ----------
    fdr_source = c("FDR", "FDR_storey", "FDR_beta", "FDR_mid", "FDR_uniform", "FDR_disc"),
    ## ---------- 绘图参数 ----------
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

    ## ===== 0. 锁定并校验网格层 =====
    g_layer <- .selectGridLayer(coordObj, grid_name) # 若 grid_name=NULL 自动唯一层
    grid_name <- if (is.null(grid_name)) {
        names(coordObj@grid)[vapply(coordObj@grid, identical, logical(1), g_layer)]
    } else {
        as.character(grid_name)
    }

    .checkGridContent(coordObj, grid_name) # 若缺必要字段将直接报错

    ## ===== 1. 读取 LeeStats 对象及其矩阵 =====
    ##  先定位 LeeStats 层（可在 @stats[[grid]] 或 @grid[[grid]]）
    leeStat <- if (!is.null(coordObj@stats[[grid_name]]) &&
        !is.null(coordObj@stats[[grid_name]][[lee_stats_layer]])) {
        coordObj@stats[[grid_name]][[lee_stats_layer]]
    } else if (!is.null(g_layer[[lee_stats_layer]])) {
        g_layer[[lee_stats_layer]]
    } else {
        stop(
            "Cannot find layer '", lee_stats_layer,
            "' in grid '", grid_name, "'."
        )
    }

    Lmat <- .getLeeMatrix(coordObj,
        grid_name = grid_name,
        lee_layer = lee_stats_layer
    ) # 仅提取对齐后的 L

    ## ===== 2. 基因子集 =====
    all_genes <- rownames(Lmat)
    if (is.null(gene_subset)) {
        keep_genes <- all_genes
    } else if (is.character(gene_subset)) {
        ## (A) 直接给定基因向量
        keep_genes <- intersect(
            all_genes,
            .getGeneSubset(coordObj, genes = gene_subset)
        )
    } else if (is.list(gene_subset) &&
        all(c("cluster_col", "cluster_num") %in% names(gene_subset))) {
        ## (B) 用聚类列 + 编号语法：gene_subset = list(cluster_col = "..", cluster_num = ..)
        keep_genes <- intersect(
            all_genes,
            .getGeneSubset(coordObj,
                cluster_col = gene_subset$cluster_col,
                cluster_num = gene_subset$cluster_num
            )
        )
    } else {
        stop("`gene_subset` 必须是字符向量或 list(cluster_col, cluster_num)")
    }
    if (length(keep_genes) < 2) {
        stop("Less than two genes remain after sub‑setting.")
    }

    ## ===== 3. 获取 (共识) 图或重新构图 =====
    use_consensus <- isTRUE(use_consensus_graph) &&
        !is.null(leeStat[[graph_slot_name]])

    if (use_consensus && isTRUE(use_FDR)) {
        # 共识图模式通常已预先过滤，这里不给二次 FDR 筛选
        message("[geneSCOPE] use_consensus_graph=TRUE ignoring FDR re-filtering (assuming graph already filtered).")
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
        ## 需要时更新权重 / 符号
    } else {
        ## ---------- 按阈值重新构图 ----------
        idx <- match(keep_genes, rownames(Lmat))
        A <- Lmat[idx, idx, drop = FALSE]

        ## (i) 先按表达占比 (pct_min) 做早期过滤
        A <- .filter_matrix_by_quantile(A, pct_min, "q100") # 修复函数名

        ## (ii) 限定最大绝对值
        A[abs(A) > L_max] <- 0

        ## (iii) 正/负阈值
        L_min_neg <- if (is.null(L_min_neg)) L_min else abs(L_min_neg)
        if (weight_abs) {
            A[abs(A) < L_min] <- 0
        } else {
            A[A > 0 & A < L_min] <- 0
            A[A < 0 & abs(A) < L_min_neg] <- 0
        }

        ## (iv) 95% 置信区间过滤
        if (CI95_filter) {
            curve <- leeStat[[curve_layer]]
            if (is.null(curve)) {
                stop("Cannot find `curve_layer = ", curve_layer, "` in LeeStats.")
            }
            f_lo <- approxfun(curve$Pear, curve$lo95, rule = 2)
            f_hi <- approxfun(curve$Pear, curve$hi95, rule = 2)

            rMat <- .getPearsonMatrix(coordObj, level = "cell") # 取单细胞 Pearson
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
            ## -------- 新增: FDR 矩阵选择逻辑 --------
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
                    "use_FDR=TRUE 但未找到任何可用 FDR 矩阵 (期待: ",
                    paste(pref_order, collapse = ", "), ")."
                )
            }
            if (inherits(FDR_sel, "big.matrix")) {
                message(
                    "[geneSCOPE] Using FDR source '", FDR_used_name,
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
            message(sprintf("[geneSCOPE] Using FDR source '%s' (FDR_max=%.3g)", FDR_used_name, FDR_max))
        }

        ## (vi) 对称化并去对角
        A <- (A + t(A)) / 2
        diag(A) <- 0

        ## (vii) 去掉孤立点
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
        e_idx <- igraph::as_edgelist(g, names = FALSE) # 替换 get.edgelist (弃用)
        igraph::E(g)$sign <- ifelse(A[e_idx] < 0, "neg", "pos")
    }

    ## ===== 4. 权重修正 & 全局缩放 =====
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

    ## ===== 5. 节点颜色 / 聚类标签 =====
    Vnames <- igraph::V(g)$name
    deg_vec <- igraph::degree(g)

    ## ① 由 cluster_vec 或 meta.data 派生 clu 向量
    clu <- setNames(rep(NA_character_, length(Vnames)), Vnames)
    is_factor_input <- FALSE
    factor_levels <- NULL

    ## helper – 从 meta.data 拿列
    get_meta_cluster <- function(col) {
        if (is.null(coordObj@meta.data) || !(col %in% colnames(coordObj@meta.data))) {
            stop("Column ‘", col, "’ not found in coordObj@meta.data.")
        }
        tmp <- coordObj@meta.data[[col]]
        names(tmp) <- rownames(coordObj@meta.data)
        tmp[Vnames]
    }

    if (!is.null(cluster_vec)) {
        if (length(cluster_vec) > 1) { # 手动向量
            if (is.factor(cluster_vec)) {
                is_factor_input <- TRUE
                factor_levels <- levels(cluster_vec)
            }
            tmp <- as.character(cluster_vec[Vnames])
            clu[!is.na(tmp)] <- tmp[!is.na(tmp)]
        } else { # 给列名
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




    ## ② 生成调色板
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

    ## ③ 边颜色（正负 / 强度渐变）
    e_idx <- igraph::as_edgelist(g, names = FALSE)
    w_norm <- igraph::E(g)$weight / max(igraph::E(g)$weight)
    cluster_ord <- setNames(seq_along(uniq_clu), uniq_clu)
    edge_cols <- character(igraph::ecount(g))
    for (i in seq_along(edge_cols)) {
        v1 <- Vnames[e_idx[i, 1]]
        v2 <- Vnames[e_idx[i, 2]]
        if (show_sign && igraph::E(g)$sign[i] == "neg") {
            edge_cols[i] <- "gray40"
        } else { # 正相关用渐变
            ref_col <- if (!is.na(clu[v1]) && !is.na(clu[v2]) && clu[v1] == clu[v2]) {
                basecol[v1]
            } else {
                # 以度数或 cluster 顺序来决定主色
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

    ## ===== 6. ggraph 绘图 =====
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

    ## ---- 边宽 / 颜色 / 线型图层 ----
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

    ## ---- 节点点位 / 文本 ----
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
            range = c(vertex_size / 2, vertex_size * 1.5)
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

    ## (可选) 调整图例顺序
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
#' @inheritParams plotNetworkGenes
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
    coordObj,
    ## ---------- 数据层 ----------
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    gene_subset = NULL,
    ## ---------- 过滤阈值 ----------
    L_min = 0,
    drop_isolated = TRUE,
    ## ---------- 网络来源 ----------
    use_consensus_graph = TRUE,
    graph_slot_name = "g_consensus",
    ## ---------- 聚类标签 ----------
    cluster_vec = NULL,
    ## ---------- Δ-PageRank 相关 ----------
    IDelta_col_name = NULL, # ★ 新增，NULL 表示不开启 Δ-PageRank
    damping = 0.85,
    weight_low_cut = 0,
    ## ---------- 树状化 ----------
    ## ---------- 可视化细节（其余保持不变） ----------
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
    tree_layout <- TRUE # 固定为树状布局
    ## ========= 0. 读共识图 ========
    g_layer <- .selectGridLayer(coordObj, grid_name)
    grid_name <- if (is.null(grid_name)) {
        names(coordObj@grid)[vapply(coordObj@grid, identical, logical(1), g_layer)]
    } else {
        grid_name
    }
    leeStat <- if (!is.null(coordObj@stats[[grid_name]]) &&
        !is.null(coordObj@stats[[grid_name]][[lee_stats_layer]])) {
        coordObj@stats[[grid_name]][[lee_stats_layer]]
    } else {
        g_layer[[lee_stats_layer]]
    }

    g_raw <- leeStat[[graph_slot_name]] # 修复变量名
    if (is.null(g_raw)) {
        stop(
            "plotDendroNetwork: 未找到共识图 graph_slot_name='", graph_slot_name,
            "' 于 leeStat[[\"", graph_slot_name, "\"]] 中。请确认 clusterGenes 已生成并写入。"
        )
    }
    stopifnot(inherits(g_raw, "igraph"))

    ## ========= 1. 子集基因 =========
    keep_genes <- rownames(coordObj@meta.data)
    if (!is.null(gene_subset)) {
        keep_genes <- intersect(
            keep_genes,
            .getGeneSubset(coordObj, genes = gene_subset)
        )
    }
    g <- igraph::induced_subgraph(g_raw, intersect(V(g_raw)$name, keep_genes))
    if (drop_isolated) g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))
    if (igraph::vcount(g) < 2) stop("子图节点不足 2 个。")

    ## ========= 2. Δ-PageRank 重新加权  =========
    if (!is.null(IDelta_col_name)) {
        delta <- coordObj@meta.data[V(g)$name, IDelta_col_name, drop = TRUE]
        delta[is.na(delta)] <- median(delta, na.rm = TRUE)
        q <- stats::quantile(delta, c(.1, .9))
        delta <- pmax(pmin(delta, q[2]), q[1]) # 与 dendroRW 同
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

    ## ========= 3. 全局缩放 / L_min =========
    if (!is.null(length_scale) && length_scale != 1) {
        igraph::E(g)$weight <- igraph::E(g)$weight * length_scale
    }
    if (!is.null(L_min) && L_min > 0) {
        keep_e <- which(igraph::E(g)$weight >= L_min)
        g <- igraph::subgraph.edges(g, keep_e, delete.vertices = TRUE)
    }

    ## ========= 4. 取 cluster 标签 =========
    Vnames <- V(g)$name
    clu <- rep(NA_character_, length(Vnames))
    names(clu) <- Vnames
    if (!is.null(cluster_vec)) {
        cv <- if (length(cluster_vec) == 1) {
            coordObj@meta.data[Vnames, cluster_vec, drop = TRUE]
        } else {
            cluster_vec[Vnames]
        }
        clu[!is.na(cv)] <- as.character(cv[!is.na(cv)])
    }
    g <- igraph::induced_subgraph(g, names(clu)[!is.na(clu)])
    Vnames <- V(g)$name
    clu <- clu[Vnames]

    ## ========= 5. ★ cluster-MST → 总树 ★ =========
    ## —— 5.1 各 cluster 内 MST ——
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

    ## —— 5.2 cluster 间 MST ——
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
            ## 连通分量逐个取 MST
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
                # 选回原图中对应 weight 最大的一条
                for (k in ks$key) {
                    cand <- inter[inter$pair == k, ]
                    cand <- cand[order(cand$weight, decreasing = TRUE), ]
                    keep_key <- c(keep_key, cand$key[1])
                }
            }
        }
    }

    ## —— 5.3 根据 keep_key 过滤边 ——
    keep_eid <- which(all_edges$key %in% unique(keep_key))
    g <- igraph::subgraph.edges(g, keep_eid, delete.vertices = TRUE)
    g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))

    ## ============ ★ 6. tree_layout 关键步骤 ★ ==============
    ## 为后续映射保留原始边 ID
    igraph::E(g)$eid <- seq_len(igraph::ecount(g))

    ## ---- 6.1  每个 cluster 的 MST ----
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

    ## ---- 6.2  cluster 之间再做 MST（在 cluster-graph 上） ----
    ## ---- 6.2  cluster 之间主干 + 高权重边重投影 ----
    ## 参数:
    ##   k_top        : 需要重投影的非树边数量上限 (可为 NULL)
    ##   w_extra_min  : 需要重投影的非树边权重下限 (可为 NULL)
    if (length(unique(clu)) > 1) {
        ed_tab <- igraph::as_data_frame(g, what = "edges")
        ed_tab$eid <- igraph::E(g)$eid
        ed_tab$cl1 <- clu[ed_tab$from]
        ed_tab$cl2 <- clu[ed_tab$to]
        inter <- ed_tab[ed_tab$cl1 != ed_tab$cl2, ] # 只看跨 cluster 边

        if (nrow(inter)) {
            ## ---- 6.2a  构建 cluster-level 图并取 MST ----
            inter$pair <- ifelse(inter$cl1 < inter$cl2,
                paste(inter$cl1, inter$cl2, sep = "|"),
                paste(inter$cl2, inter$cl1, sep = "|")
            )
            agg <- aggregate(weight ~ pair + cl1 + cl2, data = inter, max)
            g_clu <- igraph::graph_from_data_frame(
                agg[, c("cl1", "cl2", "weight")],
                directed = FALSE, vertices = unique(clu)
            )

            ## 可能不连通，逐分量取 MST
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

            ## ---- 6.2b  非树高权重边筛选 ----
            extra_edges <- inter[!(inter$pair %in% keep_pairs), ]
            if (nrow(extra_edges)) {
                if (!is.null(w_extra_min)) {
                    extra_edges <- extra_edges[extra_edges$weight >= w_extra_min, ]
                }
                if (!is.null(k_top) && nrow(extra_edges) > k_top) {
                    extra_edges <- extra_edges[order(extra_edges$weight, decreasing = TRUE)[seq_len(k_top)], ]
                }
            }
            ## ---- 6.2c  汇总需要保留的 edge ID ----
            keep_eid <- c(
                keep_eid,
                inter_MST$eid,
                if (nrow(extra_edges)) extra_edges$eid else integer(0)
            )
        }
    }
    keep_eid <- unique(keep_eid)
    g <- igraph::delete_edges(g, igraph::E(g)[!eid %in% keep_eid])
    ## 删除完全孤立点
    g <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))
    ## ============ tree_layout 结束 ==============
    ## —— 6.3  计算最终用于画图的跨簇边 ——
    # 把当前 g 的所有边转成 data.frame
    edf_final <- igraph::as_data_frame(g, what = "edges")
    # 挑出 endpoints 属于不同簇的那些
    cross_edges <- edf_final[clu[edf_final$from] != clu[edf_final$to], ]


    ## ===== 7. (重新) 计算度数 / 颜色等 =====
    Vnames <- igraph::V(g)$name
    deg_vec <- igraph::degree(g)


    is_factor_input <- FALSE
    factor_levels <- NULL

    ## helper – 从 meta.data 拿列
    get_meta_cluster <- function(col) {
        if (is.null(coordObj@meta.data) || !(col %in% colnames(coordObj@meta.data))) {
            stop("Column ‘", col, "’ not found in coordObj@meta.data.")
        }
        tmp <- coordObj@meta.data[[col]]
        names(tmp) <- rownames(coordObj@meta.data)
        tmp[Vnames]
    }

    if (!is.null(cluster_vec)) {
        if (length(cluster_vec) > 1) { # 手动向量
            if (is.factor(cluster_vec)) {
                is_factor_input <- TRUE
                factor_levels <- levels(cluster_vec)
            }
            tmp <- as.character(cluster_vec[Vnames])
            clu[!is.na(tmp)] <- tmp[!is.na(tmp)]
        } else { # 给列名
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

    ## ===== 8. 边颜色（与旧逻辑相同） =====
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

    ## ===== 9. ggraph 绘图 =====
    library(ggraph)
    library(ggplot2)
    library(dplyr)
    set.seed(seed)

    tree_mode <- match.arg(tree_mode) # 新增参数，默认 "rooted"

    if (tree_mode == "rooted") {
        root_v <- V(g)[which.max(deg_vec)] # 和之前一样
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
            range = c(vertex_size / 2, vertex_size * 1.5)
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
#' @title Random-walk-weighted multi-cluster dendrogram
#'
#' @description
#'   Constructs a dendrogram of selected gene clusters by propagating
#'   personalized PageRank scores over a Lee's-L consensus graph, rescaling
#'   edge weights, transforming similarities to distances, and applying
#'   hierarchical clustering.  The procedure highlights genes that are
#'   highly connected within and between user-specified clusters, producing
#'   an intuitive tree that can be plotted or further analysed.
#'
#' @param coordObj        A \code{CoordObj} containing \code{@grid} with a
#'                        consensus gene–gene graph and associated metadata.
#' @param grid_name       Character. Grid sub-layer to use; if \code{NULL}
#'                        the first layer is chosen automatically.
#' @param lee_stats_layer Name of the Lee's statistics layer
#'                        (default \code{"LeeStats_Xz"}).
#' @param graph_slot      Slot that stores the consensus graph
#'                        (default \code{"g_consensus"}).
#' @param cluster_name    Column in \code{meta.data} that assigns genes to
#'                        clusters.
#' @param cluster_ids     Character vector of cluster labels to include;
#'                        if \code{NULL} all labels in \code{cluster_name}
#'                        are used.
#' @param cmh_slot        Slot containing the Morisita–Horn graph
#'                        (default \code{"g_morisita"}).
#' @param use_mh_avg      Logical. If \code{TRUE} average Lee's L and
#'                        Morisita–Horn when both are present.
#' @param IDelta_col_name Column in \code{meta.data} with Δ statistics used
#'                        to build the personalised PageRank vector.
#' @param linkage         Linkage method passed to \code{hclust}
#'                        (e.g. \code{"average"}, \code{"complete"}).
#' @param dist_mode       Either \code{"one_minus"} to use
#'                        \eqn{1 - S / \max(S)} or \code{"inverse"} to use
#'                        \eqn{1/(S + 1e-6)} as the distance measure.
#' @param plot_dend       Logical. Draw the dendrogram if \code{TRUE}
#'                        (default).
#' @param weight_low_cut  Numeric. Edge weights at or below this value are
#'                        set to zero after PageRank scaling.
#' @param damping         Damping factor for personalised PageRank
#'                        (range 0–1; default 0.85).
#'
#' @return Invisibly returns a list with components:
#'   \describe{
#'     \item{\code{dend}}{A \code{dendrogram} object.}
#'     \item{\code{hclust}}{The underlying \code{hclust} object.}
#'     \item{\code{dist}}{The distance matrix used for clustering.}
#'     \item{\code{genes}}{Character vector of genes included.}
#'     \item{\code{cluster_map}}{Named vector mapping genes to clusters.}
#'     \item{\code{PageRank}}{PageRank scores for the included genes.}
#'   }
#'
#' @details
#'   Edge weights in the consensus graph are first re-weighted by the mean
#'   personalised PageRank of their incident vertices, emphasising genes that
#'   rank highly under the chosen \code{IDelta_col_name}.  Very weak edges
#'   (\code{<= weight_low_cut}) are removed.  The remaining similarity matrix
#'   is converted to a distance matrix according to \code{dist_mode} and fed
#'   to \code{hclust}.  The default \code{linkage = "average"} reproduces the
#'   legacy behaviour.  Setting \code{plot_dend = TRUE} calls the base
#'   graphics plotting method with a descriptive title.
#'
#' @examples
#' \dontrun{
#' dend_res <- buildMultiClusterDendrogramRW(
#'     coordObj        = P5.coord,
#'     grid_name       = "25um",
#'     cluster_name    = "modL0.15",
#'     cluster_ids     = c("C3", "C5", "C7"),
#'     IDelta_col_name = "25um_iDelta_raw",
#'     linkage         = "average",
#'     dist_mode       = "one_minus",
#'     plot_dend       = TRUE
#' )
#' }
#' @importFrom igraph edge_attr_names as_edgelist as_data_frame page_rank induced_subgraph ecount
#' @importFrom stats hclust as.dendrogram quantile
#' @importFrom graphics plot
#' @export
buildMultiClusterDendrogramRW <- function(
    coordObj,
    grid_name = NULL,
    lee_stats_layer = "LeeStats_Xz",
    graph_slot = "g_consensus",
    cluster_name,
    cluster_ids = NULL,
    cmh_slot = "g_morisita",
    use_mh_avg = TRUE,
    IDelta_col_name,
    linkage = "average",
    dist_mode = c("one_minus", "inverse"),
    plot_dend = TRUE,
    weight_low_cut = 0,
    damping = 0.85) {
    dist_mode <- match.arg(dist_mode)

    ## ---------- 0. 锁定网格层 & LeeStats ----------
    g_layer <- .selectGridLayer(coordObj, grid_name)
    if (is.null(grid_name)) {
        grid_name <- names(coordObj@grid)[
            vapply(coordObj@grid, identical, logical(1), g_layer)
        ]
    }

    ## LeeStats 对象（备用）
    leeStat <- if (!is.null(coordObj@stats[[grid_name]]) &&
        !is.null(coordObj@stats[[grid_name]][[lee_stats_layer]])) {
        coordObj@stats[[grid_name]][[lee_stats_layer]]
    } else if (!is.null(g_layer[[lee_stats_layer]])) {
        g_layer[[lee_stats_layer]]
    } else {
        stop(
            "Cannot find layer '", lee_stats_layer,
            "' in grid '", grid_name, "'."
        )
    }

    ## 统一用 .getLeeMatrix() 取 L，并对齐 Pearson 基因集（若存在）
    Lmat <- .getLeeMatrix(coordObj,
        grid_name = grid_name,
        lee_layer = lee_stats_layer
    )

    ## ---------- 1. 共识图 ----------
    g <- leeStat[[graph_slot]]
    stopifnot(inherits(g, "igraph"))

    ## 用 Lmat 重写 / 补全边权（Lee’s L 绝对值 → 正值）
    ed_l <- igraph::as_edgelist(g, names = TRUE)
    w_lmat <- mapply(function(a, b) abs(Lmat[a, b]),
        ed_l[, 1], ed_l[, 2],
        USE.NAMES = FALSE
    )

    ## 若 Lmat 没含此边（NA），保留原权或置 0
    if (!"weight" %in% igraph::edge_attr_names(g)) {
        igraph::E(g)$weight <- w_lmat
    } else {
        w_now <- igraph::E(g)$weight
        w_now[is.na(w_now) | w_now <= 0] <- w_lmat[is.na(w_now) | w_now <= 0]
        igraph::E(g)$weight <- w_now
    }

    ## ---------- 2. 解析聚类 ----------
    memb <- coordObj@meta.data[[cluster_name]]
    memb <- if (is.factor(memb)) as.character(memb) else memb
    if (is.null(cluster_ids)) cluster_ids <- sort(unique(na.omit(memb)))
    cluster_ids <- as.character(cluster_ids)

    gene_all <- rownames(coordObj@meta.data)[memb %in% cluster_ids]
    if (length(gene_all) < 2) {
        stop("Need at least two genes in the selected clusters.")
    }

    ## ---------- 3. 若仍缺权重则构造 ----------
    if (max(igraph::E(g)$weight, na.rm = TRUE) <= 1) {
        L_vals <- pmax(igraph::E(g)$weight, 0)
        w_L <- log1p(L_vals)

        cmh_vals <- if ("CMH" %in% igraph::edge_attr_names(g)) {
            igraph::E(g)$CMH
        } else if (!is.null(leeStat[[cmh_slot]])) {
            g_cmh <- leeStat[[cmh_slot]]
            ed_c <- igraph::as_data_frame(g_cmh, "edges")
            key_c <- paste(pmin(ed_c$from, ed_c$to),
                pmax(ed_c$from, ed_c$to),
                sep = "|"
            )
            cmh_map <- setNames(ed_c$CMH, key_c)
            key_now <- paste(pmin(ed_l[, 1], ed_l[, 2]),
                pmax(ed_l[, 1], ed_l[, 2]),
                sep = "|"
            )
            cmh_map[key_now]
        } else {
            1
        }

        if ("MH" %in% igraph::edge_attr_names(g) && use_mh_avg) {
            cmh_vals <- rowMeans(cbind(cmh_vals, igraph::E(g)$MH), na.rm = TRUE)
        }

        w_base <- w_L * cmh_vals
        w_base[is.na(w_base) | w_base <= weight_low_cut] <- 0
        igraph::E(g)$weight <- w_base
    }

    ## ---------- 4. Personalized PageRank ----------
    delta <- coordObj@meta.data[[IDelta_col_name]]
    delta[is.na(delta)] <- median(delta, na.rm = TRUE)
    q <- quantile(delta, c(.1, .9))
    delta <- pmax(pmin(delta, q[2]), q[1])
    pers <- exp(delta - max(delta))
    pers <- pers / sum(pers)
    names(pers) <- rownames(coordObj@meta.data)

    pers_g <- pers[V(g)$name]
    pers_g <- pers_g / sum(pers_g)

    pr <- igraph::page_rank(
        g,
        personalized = pers_g,
        damping      = damping,
        weights      = igraph::E(g)$weight,
        directed     = FALSE
    )$vector

    ## ---------- 5. 调整边权（乘 PR） ----------
    ed_tbl <- igraph::as_data_frame(g, "edges")
    w_new <- igraph::E(g)$weight * ((pr[ed_tbl$from] + pr[ed_tbl$to]) / 2)
    w_new[w_new <= weight_low_cut] <- 0
    igraph::E(g)$weight <- w_new

    ## ---------- 6. 相似度 → 距离 ----------
    subg <- igraph::induced_subgraph(g, vids = gene_all)
    n <- length(gene_all)
    S <- matrix(0, n, n, dimnames = list(gene_all, gene_all))

    if (igraph::ecount(subg)) {
        ed2 <- igraph::as_data_frame(subg, "edges")
        idx <- cbind(match(ed2$from, gene_all), match(ed2$to, gene_all))
        S[idx] <- ed2$weight
        S[cbind(idx[, 2], idx[, 1])] <- ed2$weight
    }

    max_w <- max(S)
    if (max_w == 0) stop("All edge weights are zero after random-walk weighting.")

    D <- if (dist_mode == "one_minus") {
        1 - S / max_w
    } else {
        M <- 1 / (S + 1e-6)
        diag(M) <- 0
        M
    }

    hc <- hclust(as.dist(D), method = linkage)
    dend <- as.dendrogram(hc)

    ## ---------- 7. 可视化 ----------
    if (plot_dend) {
        op <- par(no.readonly = TRUE)
        on.exit(par(op), add = TRUE)
        plot(dend,
            main = "Weighted Multi-Cluster Dendrogram",
            ylab = "Distance (1 − Weighted Similarity)"
        )
    }

    invisible(list(
        dend        = dend,
        hclust      = hc,
        dist        = as.dist(D),
        genes       = gene_all,
        cluster_map = memb[gene_all],
        PageRank    = pr[gene_all]
    ))
}


#' @title Compare cluster assignments produced by multiple CI³ parameter sets
#'
#' @description
#'   Generates a dot plot that juxtaposes gene–cluster memberships obtained
#'   under different CI³ thresholds or weighting options.  The function takes a
#'   set of membership columns already stored in \code{coordObj@meta.data},
#'   reshapes them into long format, applies an adaptive discrete palette, and
#'   saves the figure to disk.
#'
#' @param coordObj   A \code{CoordObj} whose \code{@meta.data} contains the
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
#'     coordObj,
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
plotClusterComparison <- function(coordObj,
                                  method_cols,
                                  method_labels = method_cols,
                                  point_size = 3,
                                  palette = RColorBrewer::brewer.pal) {
    stopifnot(length(method_cols) == length(method_labels))

    df <- coordObj@meta.data |>
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

    invisible(p)
}
#' @title 统一聚类结果数据类型
#' @description 将CoordObj中指定的聚类结果列统一转换为因子类型
#' @param coordObj CoordObj对象，包含聚类结果
#' @param grid_sizes 要处理的网格尺寸向量
#' @param verbose 是否打印详细信息
#' @return 修改后的CoordObj对象
unifyClusteringTypes <- function(coordObj, grid_sizes = c(10, 30, 55), verbose = TRUE) {
    # 获取所有列名
    all_cols <- colnames(coordObj@meta.data)

    # 计数统计
    converted_count <- 0
    already_factor <- 0

    # 对每个网格尺寸处理
    for (i in grid_sizes) {
        # 查找与当前网格尺寸相关的列
        grid_pattern <- paste0("_grid", i, "$|_grid", i, "_")
        matching_cols <- grep(grid_pattern, all_cols, value = TRUE)

        if (length(matching_cols) == 0) {
            if (verbose) message("[geneSCOPE] No clustering results found for grid size ", i)
            next
        }

        if (verbose) message("[geneSCOPE] Processing grid size ", i, ", found ", length(matching_cols), " clustering results")

        # 处理每列
        for (col in matching_cols) {
            # 跳过已经是因子类型的列
            if (is.factor(coordObj@meta.data[[col]])) {
                already_factor <- already_factor + 1
                next
            }

            # 转换为因子类型（保留NA值）
            values <- coordObj@meta.data[[col]]
            if (all(is.na(values))) {
                # 全部为NA的列，仍然创建一个空因子
                coordObj@meta.data[[col]] <- factor(values)
            } else {
                non_na <- values[!is.na(values)]
                # 确保有序的水平
                levels <- sort(unique(non_na))
                # 转换为因子
                coordObj@meta.data[[col]] <- factor(values, levels = levels)
            }
            converted_count <- converted_count + 1
        }
    }

    if (verbose) {
        message("[geneSCOPE] Converted ", converted_count, " columns to factor type")
        message("[geneSCOPE] Already had ", already_factor, " columns as factor type")
    }

    return(coordObj)
}

#' @title 获取有效的聚类列名
#' @description 根据网格尺寸和聚类参数模式生成有效的列名列表
#' @param coordObj CoordObj对象
#' @param grid_size 网格尺寸
#' @param include_cmh 是否包含cmh相关列
#' @param pct_filters 百分比过滤字符串向量
#' @return 存在于coordObj中的有效列名向量
getValidClusterColumns <- function(coordObj,
                                   grid_size = 10,
                                   include_cmh = FALSE,
                                   pct_filters = c("All", "q85.0", "q90.0", "q95.0", "q99.0", "q99.9")) {
    # 确保meta.data存在
    if (is.null(coordObj@meta.data) || ncol(coordObj@meta.data) == 0) {
        stop("CoordObj中的meta.data为空")
    }

    # 获取所有列名
    all_cols <- colnames(coordObj@meta.data)

    # 构建列名模式
    column_patterns <- c()

    # 循环构建基本列名
    for (pct in pct_filters) {
        # 基础列名
        base_name <- paste0(pct, "_res0.1_grid", grid_size)
        column_patterns <- c(column_patterns, base_name)

        # log1p列名
        log1p_name <- paste0(base_name, "_log1p")
        column_patterns <- c(column_patterns, log1p_name)

        # pie0.95列名
        pie_name <- paste0(base_name, "_log1p_pie0.95")
        column_patterns <- c(column_patterns, pie_name)

        # outCI95列名
        ci_name <- paste0(base_name, "_log1p_pie0.95_outCI95")
        column_patterns <- c(column_patterns, ci_name)

        # 可选：添加cmh列名
        if (include_cmh) {
            cmh_name <- paste0(base_name, "_log1p_pie0.95_outCI95_cmh")
            column_patterns <- c(column_patterns, cmh_name)
        }
    }

    # 检查哪些列确实存在
    valid_columns <- column_patterns[column_patterns %in% all_cols]

    # 检查是否所有列都存在相同的genes
    non_na_counts <- sapply(valid_columns, function(col) {
        sum(!is.na(coordObj@meta.data[[col]]))
    })

    # All columns should have the same number of non-NA values, otherwise there might be issues
    if (length(unique(non_na_counts)) > 1) {
        message("[geneSCOPE] !!! Warning: Different clustering results contain different numbers of genes, may cause comparison issues !!!")
    }

    return(valid_columns)
}

#' @title 执行健壮的聚类比较
#' @description 统一数据类型并执行聚类比较
#' @param coordObj CoordObj对象
#' @param grid_size 网格尺寸
#' @param include_cmh 是否包含cmh相关列
#' @param unify_types 是否先统一数据类型
#' @param pct_filters 百分比过滤字符串向量
#' @param ... 传递给plotClusterComparison的其他参数
#' @return plotClusterComparison的返回值
robustClusterCompare <- function(coordObj,
                                 grid_size = 10,
                                 include_cmh = FALSE,
                                 unify_types = TRUE,
                                 pct_filters = c("All", "q85.0", "q90.0", "q95.0", "q99.0", "q99.9"),
                                 ...) {
    # 首先统一数据类型（如果需要）
    if (unify_types) {
        coordObj <- unifyClusteringTypes(coordObj, grid_sizes = grid_size)
    }

    # 获取有效的列名
    valid_columns <- getValidClusterColumns(coordObj, grid_size, include_cmh, pct_filters)

    if (length(valid_columns) == 0) {
        stop("找不到有效的聚类列，请检查网格尺寸和meta.data")
    }

    message("[geneSCOPE] Using ", length(valid_columns), " valid clustering columns for comparison")

    # 执行聚类比较
    result <- plotClusterComparison(
        coordObj = coordObj,
        method_cols = valid_columns,
        ...
    )

    return(result)
}

#' @title 自动列出现有聚类列
#' @description 基于既有 meta.data 列名与给定网格尺寸/百分位模式返回真实存在的列
#' @param coordObj CoordObj
#' @param grid_sizes 整数向量 (例如 c(10,30,55))
#' @param pct_filters 百分位标签 (需与 add / clusterGenes 生成规则一致)
#' @param include_cmh 是否包含 *_cmh 变体（仅当列真实存在才返回，不会产生警告）
#' @return data.frame(columns, grid_size, pct, variant)
#' @export
autoListClusterCols <- function(coordObj,
                                grid_sizes = c(10, 30, 55),
                                pct_filters = c("All", "q85.0", "q90.0", "q95.0", "q99.0", "q99.9"),
                                include_cmh = TRUE) {
    stopifnot(!is.null(coordObj@meta.data))
    cols <- colnames(coordObj@meta.data)
    out <- list()
    for (g in grid_sizes) {
        for (pct in pct_filters) {
            base <- paste0(pct, "_res0.1_grid", g)
            variants <- c(
                base,
                paste0(base, "_log1p"),
                paste0(base, "_log1p_pie0.95"),
                paste0(base, "_log1p_pie0.95_outCI95")
            )
            if (include_cmh) {
                variants <- c(variants, paste0(base, "_log1p_pie0.95_outCI95_cmh"))
            }
            exist <- variants[variants %in% cols]
            if (length(exist)) {
                out[[length(out) + 1]] <- data.frame(
                    column = exist,
                    grid_size = g,
                    pct = pct,
                    stringsAsFactors = FALSE
                )
            }
        }
    }
    if (!length(out)) {
        return(data.frame())
    }
    df <- do.call(rbind, out)
    df$variant <- sub("^[^_]+_res0.1_grid[0-9]+", "base", df$column)
    df$variant <- sub("^.+_log1p_pie0.95_outCI95_cmh$", "log1p_pie0.95_outCI95_cmh", df$column)
    df$variant <- ifelse(grepl("_log1p_pie0.95_outCI95_cmh$", df$column), "log1p_pie0.95_outCI95_cmh",
        ifelse(grepl("_log1p_pie0.95_outCI95$", df$column), "log1p_pie0.95_outCI95",
            ifelse(grepl("_log1p_pie0.95$", df$column), "log1p_pie0.95",
                ifelse(grepl("_log1p$", df$column), "log1p", "base")
            )
        )
    )
    df
}

#' @title 批量绘制网络与树形网络
#' @description 针对若干聚类结果列批量调用 plotNetworkGenes 与 plotDendroNetwork
#' @param coordObj CoordObj
#' @param cluster_columns 需要绘图的聚类列（默认自动探测）
#' @param grid_size_col 是否根据列名自动解析 grid (TRUE)
#' @param output_dir 根输出目录（NULL 则不保存，仅返回列表）
#' @param width,height,dpi 保存参数
#' @param ... 传递给 plotNetworkGenes 的其他参数（如 fdr_source / FDR_max 等）
#' @return 命名列表：每列含 list(network=ggplot, dendro=list(plot=..., graph=...))
#' @export
batchPlotGeneNetworks <- function(coordObj,
                                  cluster_columns = NULL,
                                  grid_sizes = c(10, 30, 55),
                                  output_dir = NULL,
                                  width = 12, height = 12, dpi = 300,
                                  ...) {
    if (is.null(cluster_columns)) {
        cc_df <- autoListClusterCols(coordObj, grid_sizes = grid_sizes, include_cmh = TRUE)
        cluster_columns <- cc_df$column
    }
    if (!length(cluster_columns)) stop("未找到可用聚类列")
    if (!is.null(output_dir) && !dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    out <- list()
    for (col in cluster_columns) {
        # 解析 grid 尺寸
        gsize <- sub("^.*_grid([0-9]+).*", "\\1", col)
        grid_name <- paste0("grid_lenGrid", gsize)
        if (!(grid_name %in% names(coordObj@grid))) {
            message("[geneSCOPE] Skipping ", col, " (missing grid layer ", grid_name, ")")
            next
        }
        # At least two non-NA values
        non_na <- sum(!is.na(coordObj@meta.data[[col]]))
        if (non_na < 2) {
            message("[geneSCOPE] Skipping ", col, " (non-NA genes < 2)")
            next
        }
        message("[geneSCOPE] Plotting: ", col, " / grid=", grid_name, " (non-NA=", non_na, ")")
        title_net <- paste0("Network - ", col, " (grid ", gsize, ")")
        p_net <- tryCatch(
            plotNetworkGenes(
                coordObj = coordObj,
                grid_name = grid_name,
                cluster_vec = col,
                use_consensus_graph = TRUE, # 使用共识图
                graph_slot_name = "g_consensus",
                title = title_net,
                ...
            ),
            error = function(e) {
                message("[geneSCOPE] !!! Warning: plotNetworkGenes failed: ", col, " : ", e$message, " !!!")
                NULL
            }
        )
        # 树状网络
        title_den <- paste0("Dendro Network - ", col, " (grid ", gsize, ")")
        p_den <- tryCatch(
            plotDendroNetwork(
                coordObj = coordObj,
                grid_name = grid_name,
                cluster_vec = col,
                graph_slot_name = "g_consensus",
                IDelta_col_name = paste0("grid_lenGrid", gsize, "_iDelta_raw"),
                title = title_den,
                show_sign = TRUE
            ),
            error = function(e) {
                message("[geneSCOPE] !!! Warning: plotDendroNetwork failed: ", col, " : ", e$message, " !!!")
                NULL
            }
        )
        out[[col]] <- list(network = p_net, dendro = p_den)

        if (!is.null(output_dir)) {
            subdir <- file.path(output_dir, paste0("grid", gsize))
            if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)
            if (!is.null(p_net)) {
                fn1 <- file.path(subdir, paste0("network_", col, ".png"))
                ggplot2::ggsave(fn1, plot = p_net, width = width, height = height, dpi = dpi)
            }
            if (!is.null(p_den) && !is.null(p_den$plot)) {
                fn2 <- file.path(subdir, paste0("dendro_network_", col, ".png"))
                ggplot2::ggsave(fn2, plot = p_den$plot, width = width, height = height, dpi = dpi)
            }
        }
    }
    invisible(out)
}

# ===== 使用示例（脚本中调用） =====
# res_list <- batchPlotGeneNetworks(P5.coord,
#                                   grid_sizes = c(10,30,55),
#                                   output_dir = "/path/to/output",
#                                   fdr_source = "FDR_beta",
#                                   FDR_max = 0.05)
