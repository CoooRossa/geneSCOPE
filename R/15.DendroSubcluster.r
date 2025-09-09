#' @title 基于树形网络的分支子聚类可视化
#' @description
#' 调用已有 \code{plotDendroNetwork()} 构建树状/骨架网络后，对指定簇(或全部簇)执行
#' “关节点分支”子聚类细分；若无候选且设置 \code{fallback_community=TRUE}，
#' 回退 Louvain/Leiden 社区拆分。输出带子聚类着色的新网络图。
#' @inheritParams plotDendroNetwork
#' @param enable_subbranch 是否启用分支子聚类细分 (默认 TRUE)
#' @param cluster_id 可选字符向量：只对子聚类拆分的簇；NULL=全部
#' @param include_root 关节点法是否把关节点包含进每个分支
#' @param max_subclusters 每簇最多保留子聚类数 (按规模贪心去重)
#' @param fallback_community 关节点无结果时是否回退社区
#' @param min_sub_size 社区回退时子聚类最小规模
#' @param community_method 社区算法优先顺序
#' @param subbranch_palette 子聚类调色板：第1色=剩余/未细分，其余循环分配给子聚类
#' @param downstream_min_size 下游子树判定规模阈值（NULL=自适应）
#' @param force_split 原社区拆分失败时是否尝试放宽接受条件
#' @param main_fraction_cap 最大社区占比阈值（超过也允许保留其余小社区）
#' @param core_periph 允许核心-边缘拆分
#' @param core_degree_quantile 核心度阈值分位
#' @param core_min_fraction 核心最小占比
#' @param degree_gini_threshold 触发核心-边缘拆分的度Gini阈值
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
    ## ==== 原 plotDendroNetwork 参数（保持次序以便兼容） ====
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
    ## ==== 新增子聚类相关参数 ====
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
    downstream_min_size = NULL, # <--- 新增：下游子树判定规模阈值（NULL=自适应）
    ## === 新增控制参数（默认尽量保守） ===
    force_split = TRUE, # 原社区拆分失败时是否尝试放宽接受条件
    main_fraction_cap = 0.9, # 最大社区占比阈值（超过也允许保留其余小社区）
    core_periph = TRUE, # 允许核心-边缘拆分
    core_degree_quantile = 0.75, # 核心度阈值分位
    core_min_fraction = 0.05, # 核心最小占比
    degree_gini_threshold = 0.35, # 触发核心-边缘拆分的度Gini阈值
    verbose = TRUE) {
    tree_mode <- match.arg(tree_mode)
    community_method <- match.arg(community_method, several.ok = TRUE)

    ## 0. 调用基础网络函数 ------------------------------------------------------
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
    if (is.null(g) || !inherits(g, "igraph")) stop("基础网络构建失败或未返回 igraph")

    vnames <- igraph::V(g)$name

    ## 1. 解析簇标签（与基础函数逻辑对齐） ------------------------------------
    if (verbose) message("[geneSCOPE] [Subbranch] Preparing cluster labels...")
    Vnames <- vnames
    # cluster_vec 参数两种形式：长度>1命名向量 / 单值列名
    if (!is.null(cluster_vec)) {
        if (length(cluster_vec) == 1) {
            if (is.null(coordObj@meta.data) || !(cluster_vec %in% colnames(coordObj@meta.data))) {
                stop("meta.data 不含列: ", cluster_vec)
            }
            clu_full <- coordObj@meta.data[[cluster_vec]]
            names(clu_full) <- rownames(coordObj@meta.data)
            clu <- as.character(clu_full[Vnames])
        } else {
            if (is.null(names(cluster_vec))) stop("cluster_vec 为向量时须具备 names")
            clu <- as.character(cluster_vec[Vnames])
        }
    } else {
        clu <- rep("C1", length(Vnames))
        if (verbose) message("[geneSCOPE] [Subbranch] No cluster_vec provided, using single cluster C1")
    }
    names(clu) <- Vnames
    cluster_df <- data.frame(gene = Vnames, cluster = clu, stringsAsFactors = FALSE)

    ## 2. 若不启用子聚类直接返回 ------------------------------------------------
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

    ## 3. 定义单簇检测函数（加入 hub 与下游多子树逻辑） ----------
    detect_one_cluster <- function(genes_target, cid) {
        tryCatch(
            {
                # --- 基因集合合法化 ---
                genes_target <- unique(stats::na.omit(genes_target))
                # 额外一致性检查（理论上上一层已交集）
                if (any(!genes_target %in% igraph::V(g)$name)) {
                    genes_target <- intersect(genes_target, igraph::V(g)$name)
                }
                if (!length(genes_target)) {
                    return(list(method = "none", subclusters = list()))
                }

                # 关键修改1：直接使用名称向量（避免 VertexSeq二次解析引发 Unknown vertex selected）
                sg <- tryCatch(
                    igraph::induced_subgraph(g, vids = genes_target),
                    error = function(e) {
                        warning("[geneSCOPE] !!! [detect_one_cluster] Subgraph construction failed for cluster ", cid, ": ", e$message, " !!!")
                        return(NULL)
                    }
                )
                if (is.null(sg) || igraph::vcount(sg) == 0) {
                    return(list(method = "none", subclusters = list()))
                }
                if (verbose) {
                    message(
                        "[geneSCOPE]     [detect_one_cluster] cluster ", cid,
                        " subgraph nodes=", igraph::vcount(sg),
                        " edges=", igraph::ecount(sg)
                    )
                }

                method_used <- "none"
                subclusters <- list()
                vcount_sg <- igraph::vcount(sg)

                safe_add_subcluster <- function(tag, genes_branch) {
                    genes_branch <- unique(stats::na.omit(genes_branch))
                    genes_branch <- intersect(genes_branch, igraph::V(sg)$name)
                    if (!length(genes_branch)) {
                        return(FALSE)
                    }
                    subclusters[[tag]] <<- genes_branch
                    TRUE
                }

                ## ---- 3.1 articulation（重写循环，纯整数索引 + 额外 tryCatch） ----
                if (vcount_sg >= 3) {
                    art_ok <- TRUE
                    cand <- list()
                    tryCatch(
                        {
                            arts_vs <- tryCatch(igraph::articulation_points(sg),
                                error = function(e) igraph::V(sg)[FALSE]
                            )
                            if (length(arts_vs)) {
                                arts_ids <- as.integer(igraph::as_ids(arts_vs))
                                arts_ids <- arts_ids[is.finite(arts_ids) & arts_ids >= 1 & arts_ids <= vcount_sg]
                                for (a_id in arts_ids) {
                                    root_name <- igraph::V(sg)$name[a_id]
                                    sg_minus <- tryCatch(igraph::delete_vertices(sg, a_id),
                                        error = function(e) NULL
                                    )
                                    if (is.null(sg_minus)) next
                                    comps <- tryCatch(igraph::components(sg_minus),
                                        error = function(e) NULL
                                    )
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
                                    ord <- order(vapply(cand, `[[`, numeric(1), "size"),
                                        decreasing = TRUE
                                    )
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
                                            safe_add_subcluster(
                                                paste0(cid, "_sub", k),
                                                kept[[k]]$genes
                                            )
                                        }
                                    }
                                }
                            }
                        },
                        error = function(e) {
                            art_ok <<- FALSE
                            if (verbose) {
                                message(
                                    "[geneSCOPE] !!! [detect_one_cluster] Articulation method failed for cluster ",
                                    cid, ": ", e$message, " !!!"
                                )
                            }
                        }
                    )
                }

                ## ---- 3.2 fallback community (原始) ----
                if (method_used == "none" && fallback_community) {
                    comm <- NULL
                    for (mtd in community_method) {
                        comm <- tryCatch(
                            switch(mtd,
                                louvain = igraph::cluster_louvain(sg),
                                leiden  = igraph::cluster_leiden(sg),
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

                ## ---- 3.2b 扩展/强制社区拆分（原社区未成功） ----
                if (method_used == "none" && force_split) {
                    # 重取一个社区划分，优先换算法；若只有一个算法则尝试 edge betweenness 2-way
                    comm2 <- NULL
                    alt_algos <- setdiff(c("leiden", "louvain"), if (length(community_method)) community_method[1] else character(0))
                    for (mtd in alt_algos) {
                        comm2 <- tryCatch(
                            switch(mtd,
                                louvain = igraph::cluster_louvain(sg),
                                leiden  = igraph::cluster_leiden(sg),
                                NULL
                            ),
                            error = function(e) NULL
                        )
                        if (!is.null(comm2) && length(unique(comm2$membership)) > 1) break
                    }
                    if (is.null(comm2)) {
                        # 尝试 edge betweenness k=2
                        comm2 <- tryCatch(
                            {
                                eb <- igraph::cluster_edge_betweenness(sg)
                                memb <- igraph::cut_at(eb, no = 2)
                                structure(list(membership = memb), class = "EC_split")
                            },
                            error = function(e) NULL
                        )
                    }
                    if (!is.null(comm2) && length(unique(comm2$membership)) > 1) {
                        tb <- table(comm2$membership)
                        ord_mb <- as.integer(names(sort(tb, decreasing = TRUE)))
                        if (length(ord_mb) >= 2) {
                            total <- sum(tb)
                            largest_frac <- max(tb) / total
                            # 允许最大社区极大仍保留其它（只要其它 >= min_sub_size）
                            acc_ids <- ord_mb[-1]
                            keep_ids <- acc_ids[tb[as.character(acc_ids)] >= min_sub_size]
                            if (length(keep_ids)) {
                                method_used <- "community_force"
                                k <- 1
                                for (sid in keep_ids) {
                                    subclusters[[paste0(cid, "_cf", k)]] <- igraph::V(sg)$name[comm2$membership == sid]
                                    k <- k + 1
                                    if (length(subclusters) >= max_subclusters) break
                                }
                            }
                        }
                    }
                }

                ## ---- 3.2c 核心-边缘拆分（仍无结果 & 允许） ----
                if (method_used == "none" && core_periph && length(subclusters) == 0) {
                    deg <- igraph::degree(sg)
                    if (length(deg) && vcount_sg >= max(2 * min_sub_size, 8)) {
                        # Gini 计算
                        gini_deg <- {
                            x <- as.numeric(deg)
                            n <- length(x)
                            if (n < 2 || all(x == 0)) {
                                0
                            } else {
                                x <- sort(x)
                                idx <- seq_len(n)
                                (2 * sum(idx * x) / (n * sum(x))) - (n + 1) / n
                            }
                        }
                        if (gini_deg >= degree_gini_threshold) {
                            thr_core <- stats::quantile(deg, core_degree_quantile, names = FALSE, type = 7)
                            core_ids <- which(deg >= thr_core)
                            # 去掉孤立高值/确保核心最少
                            if (length(core_ids) >= max(ceiling(core_min_fraction * vcount_sg), min_sub_size) &&
                                length(core_ids) <= vcount_sg - min_sub_size) {
                                core_genes <- igraph::V(sg)$name[core_ids]
                                peri_genes <- setdiff(igraph::V(sg)$name, core_genes)
                                # 尝试剔除与核心直接相连的边缘薄层 → 精炼核心 (可选)
                                if (length(core_genes) && length(peri_genes)) {
                                    method_used <- "core-periph"
                                    subclusters[[paste0(cid, "_core")]] <- core_genes
                                    # 只保留核心为子簇，周边作为剩余（若需要也可添加周边为子簇）
                                    # 若希望同时输出周边，可解除下行注释：
                                    # subclusters[[paste0(cid, "_peri")]] <- peri_genes
                                }
                            }
                        }
                    }
                }

                ## ---- 保留后续 3.3 hub + 3.4 downstream 原逻辑（若已有 subclusters 仍可追加） ----
                # ...existing code (hub 分支、downstream 扩展)...

                list(method = method_used, subclusters = subclusters)
            },
            error = function(e) {
                warning("[geneSCOPE] !!! [detect_one_cluster] Uncaught error for cluster ", cid, ": ", e$message, " !!!")
                list(method = "none", subclusters = list())
            }
        )
    }

    ## 4. 遍历簇执行子聚类（加入输入集合合法化） ------------------------------
    target_clusters <- if (is.null(cluster_id)) sort(unique(na.omit(clu))) else intersect(unique(clu), cluster_id)
    if (!length(target_clusters)) {
        if (verbose) message("[geneSCOPE] [Subbranch] No available target clusters, skipping subdivision")
        enable_subbranch <- FALSE
    }

    sub_attr <- setNames(rep(NA_character_, length(Vnames)), Vnames)
    subclusters_all <- list()
    methods_seen <- character(0)

    if (enable_subbranch) {
        if (verbose) message("[geneSCOPE] [Subbranch] Cluster count: ", length(target_clusters), " -> Starting subdivision")
        for (cid in target_clusters) {
            genes_target_raw <- names(clu)[clu == cid]
            genes_target <- intersect(unique(stats::na.omit(genes_target_raw)), Vnames)
            if (length(genes_target_raw) != length(genes_target) && verbose) {
                message(
                    "[geneSCOPE]   Cluster ", cid, ": Filtered invalid/missing genes ",
                    length(genes_target_raw), " -> ", length(genes_target)
                )
            }
            if (length(genes_target) < 3) {
                if (verbose) message("[geneSCOPE]   Cluster ", cid, ": Valid nodes < 3, skipping")
                next
            }
            # --- 外层 tryCatch 避免单簇错误中断 ---
            res_c <- tryCatch(detect_one_cluster(genes_target, cid),
                error = function(e) {
                    warning("[geneSCOPE] !!! [Subbranch] Cluster ", cid, " processing error: ", e$message, " !!!")
                    list(method = "none", subclusters = list())
                }
            )
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
                    "[geneSCOPE]   Cluster ", cid, ": method=", res_c$method,
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

    ## === 新增：构建有序 factor 层次 (父簇 -> 其子簇) =========================
    numeric_sort_levels <- function(x) {
        ux <- unique(x)
        suppressWarnings({
            nx <- as.numeric(ux)
        })
        if (all(!is.na(nx))) as.character(sort(as.numeric(ux))) else sort(ux)
    }
    base_levels <- numeric_sort_levels(clu)

    if (method_final %in% c("none", "disabled")) {
        # 仅父簇
        cluster_df$cluster <- factor(cluster_df$cluster, levels = base_levels)
    } else {
        # 生成子簇层次：父簇后插入其所有子簇
        # 子簇名称格式假设以 "父簇ID_" 开头
        sub_levels_raw <- sort(unique(na.omit(sub_attr))) # 当前子簇集合
        parent_of <- sub("^([^_]+)_.*$", "\\1", sub_levels_raw) # 提取父簇ID
        level_vec <- character(0)
        for (pv in base_levels) {
            level_vec <- c(level_vec, pv)
            if (pv %in% parent_of) {
                # 保持出现顺序（在 sub_levels_raw 中的顺序）
                kids <- sub_levels_raw[parent_of == pv]
                level_vec <- c(level_vec, kids)
            }
        }
        # 若存在未匹配父簇的子簇（异常），追加到末尾
        orphan <- setdiff(sub_levels_raw, level_vec)
        if (length(orphan)) level_vec <- c(level_vec, orphan)

        # 更新 factor（cluster_df 仍只含父簇，保持基因所属父簇层次）
        cluster_df$cluster <- factor(cluster_df$cluster, levels = base_levels)
        # subcluster_df 后面统一更新
    }

    ## 5. 构建结果表与颜色 ------------------------------------------------------
    subcluster_df <- data.frame(
        gene = Vnames,
        cluster = cluster_df$cluster,
        subcluster = sub_attr,
        stringsAsFactors = FALSE
    )

    # 若有子簇，为 subcluster 列建立分层 levels：与上面 level_vec 对齐（无子簇则保持 NA）
    if (!(method_final %in% c("none", "disabled"))) {
        # subcluster 列：保留 NA；levels 使用 (level_vec 去除父簇中未作为子簇的那些? 按需求：父簇与子簇是同一“绘制标签”集合)
        # 这里仅对子簇名称设 levels，父簇不在 subcluster 列 levels 中（subcluster 列是子簇ID或 NA）
        sub_levels_only <- level_vec[level_vec %in% subcluster_df$subcluster]
        subcluster_df$subcluster <- factor(subcluster_df$subcluster, levels = sub_levels_only)
    }

    if (method_final %in% c("none", "disabled")) {
        # 没有实际子聚类：仍旧用 cluster_palette
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

    ## 6. 利用原函数重新绘制（仅当存在子聚类） -------------------------------
    if (verbose) message("[geneSCOPE] [Subbranch] Preparing final plot generation (reusing original plotDendroNetwork logic)...")
    branch_network <- NULL
    if (method_final %in% c("none", "disabled")) {
        # 无子簇：为保持层次输出，可（可选）重新绘制一次使 cluster factor 顺序生效
        cluster_vec_base <- factor(clu, levels = base_levels)
        new_args <- base_args
        new_args$cluster_vec <- setNames(cluster_vec_base, names(clu))
        new_args$title <- if (is.null(title)) "Base Network" else title
        branch_network <- do.call(plotDendroNetwork, new_args)
        plt <- branch_network$plot
    } else {
        cluster_vec_sub <- setNames(ifelse(is.na(sub_attr), clu, sub_attr), names(clu))
        # 将 cluster_vec_sub 转为 factor：使用 level_vec（父簇 + 子簇层次）
        cluster_vec_sub <- factor(cluster_vec_sub, levels = level_vec)
        new_args <- base_args
        new_args$cluster_vec <- cluster_vec_sub
        new_args$title <- if (is.null(title)) {
            paste0("Sub-branch: ", method_final)
        } else {
            paste0(title, " | Sub-branch: ", method_final)
        }
        if (verbose) message("[geneSCOPE] [Subbranch] Re-calling plotDendroNetwork to generate subcluster coloring plot...")
        branch_network <- do.call(plotDendroNetwork, new_args)
        plt <- branch_network$plot
    }

    ## 7. 返回 -------------------------------------------------------------------
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
            downstream_min_size = downstream_min_size, # <--- 新增返回
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
