#' Internal helpers for cross-sample Lee's L comparisons
#' @keywords internal

.is_scope_object <- function(x) inherits(x, "scope_object")

.leesl_compare_pick_column <- function(df, candidates) {
    hit <- intersect(candidates, names(df))
    if (length(hit)) hit[1] else NA_character_
}

.leesl_compare_standardize_pairs <- function(pairs) {
    if (is.null(pairs)) return(NULL)
    if (is.matrix(pairs)) pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
    if (!is.data.frame(pairs)) stop("pairs must be a data.frame or matrix")

    i_col <- .leesl_compare_pick_column(pairs, c("gene_i", "gene1", "gene_1", "g1", "geneA", "gene_a"))
    j_col <- .leesl_compare_pick_column(pairs, c("gene_j", "gene2", "gene_2", "g2", "geneB", "gene_b"))
    if (is.na(i_col) || is.na(j_col)) {
        if (ncol(pairs) >= 2) {
            i_col <- names(pairs)[1]
            j_col <- names(pairs)[2]
        } else {
            stop("pairs must have at least two columns for gene identifiers")
        }
    }

    gene_i <- as.character(pairs[[i_col]])
    gene_j <- as.character(pairs[[j_col]])
    keep <- !is.na(gene_i) & !is.na(gene_j) & nzchar(gene_i) & nzchar(gene_j)
    if (!any(keep)) return(data.frame())
    gene_i <- gene_i[keep]
    gene_j <- gene_j[keep]

    swap <- gene_i > gene_j
    if (any(swap)) {
        tmp <- gene_i[swap]
        gene_i[swap] <- gene_j[swap]
        gene_j[swap] <- tmp
    }

    out <- data.frame(
        gene_i = gene_i,
        gene_j = gene_j,
        stringsAsFactors = FALSE
    )
    out$pair_id <- paste(out$gene_i, out$gene_j, sep = "|")
    out <- out[!duplicated(out$pair_id), , drop = FALSE]
    rownames(out) <- NULL
    out
}

.leesl_compare_expand_pairs <- function(df, pairs_df, sample_id) {
    if (is.null(pairs_df) || !is.data.frame(pairs_df) || !nrow(pairs_df)) return(df)
    if (!is.data.frame(df) || !nrow(df)) {
        df <- data.frame()
    }
    merged <- merge(pairs_df, df, by = "pair_id", all.x = TRUE, sort = FALSE, suffixes = c("_ref", ""))
    if ("gene_i_ref" %in% names(merged)) {
        merged$gene_i <- ifelse(is.na(merged$gene_i), merged$gene_i_ref, merged$gene_i)
        merged$gene_j <- ifelse(is.na(merged$gene_j), merged$gene_j_ref, merged$gene_j)
        merged$gene_i_ref <- NULL
        merged$gene_j_ref <- NULL
    }
    merged$sample_id <- sample_id
    if (!"gene_missing" %in% names(merged)) merged$gene_missing <- FALSE
    missing_row <- is.na(merged$L) & is.na(merged$Pearson_r) & is.na(merged$delta_Lr)
    merged$gene_missing[missing_row] <- TRUE
    merged
}

.leesl_compare_attach_sample_info <- function(df, sample_info) {
    if (is.null(sample_info) || !is.data.frame(sample_info) || !nrow(sample_info)) return(df)
    if (!"sample_id" %in% names(sample_info)) {
        stop("sample_info must contain a sample_id column")
    }
    keep_cols <- unique(c("sample_id", "group", "batch", "donor",
        "edge_density", "components", "modularity", "modularity_Q",
        "mean_degree", "sd_degree", "hub_ratio", "sig_edge_frac"))
    extra <- intersect(names(sample_info), keep_cols)
    sample_info_sub <- unique(sample_info[, extra, drop = FALSE])
    add_cols <- setdiff(names(sample_info_sub), names(df))
    if (!length(add_cols)) return(df)
    sample_info_sub <- sample_info_sub[, c("sample_id", add_cols), drop = FALSE]
    merge(df, sample_info_sub, by = "sample_id", all.x = TRUE, sort = FALSE)
}

.leesl_compare_standardize_pair_table <- function(df, sample_id = NULL, sample_info = NULL) {
    if (!is.data.frame(df)) stop("pair table must be a data.frame")
    df <- as.data.frame(df, stringsAsFactors = FALSE)

    if (!is.null(sample_id)) {
        df$sample_id <- as.character(sample_id)
    } else if (!"sample_id" %in% names(df)) {
        stop("pair table must include sample_id or be provided via sample_id argument")
    }

    gene_i_col <- .leesl_compare_pick_column(df, c("gene_i", "gene1", "gene_1", "g1", "geneA", "gene_a"))
    gene_j_col <- .leesl_compare_pick_column(df, c("gene_j", "gene2", "gene_2", "g2", "geneB", "gene_b"))
    if (is.na(gene_i_col) || is.na(gene_j_col)) stop("pair table must include gene identifiers")

    L_col <- .leesl_compare_pick_column(df, c("L", "LeesL", "LeeL", "Lee_L", "lee_l", "leeL", "lee_L"))
    r_col <- .leesl_compare_pick_column(df, c("Pearson_r", "Pear", "pear", "r", "pearson_r"))
    delta_col <- .leesl_compare_pick_column(df, c("delta_Lr", "Delta", "delta", "Delta_Lr", "delta_lr"))

    out <- data.frame(
        sample_id = as.character(df$sample_id),
        gene_i = as.character(df[[gene_i_col]]),
        gene_j = as.character(df[[gene_j_col]]),
        stringsAsFactors = FALSE
    )

    out$L <- if (!is.na(L_col)) as.numeric(df[[L_col]]) else NA_real_
    out$Pearson_r <- if (!is.na(r_col)) as.numeric(df[[r_col]]) else NA_real_
    if (!is.na(delta_col)) {
        out$delta_Lr <- as.numeric(df[[delta_col]])
    } else if (!is.na(L_col) && !is.na(r_col)) {
        out$delta_Lr <- out$L - out$Pearson_r
    } else {
        out$delta_Lr <- NA_real_
    }

    # Uncertainty fields (L)
    map_numeric <- function(name, candidates) {
        col <- .leesl_compare_pick_column(df, candidates)
        if (!is.na(col)) out[[name]] <<- as.numeric(df[[col]])
    }
    map_numeric("L_se", c("L_se", "se_L", "L_sem", "LeeL_se"))
    map_numeric("L_var", c("L_var", "var_L"))
    map_numeric("L_sd", c("L_sd", "sd_L"))
    map_numeric("L_ci_lo", c("L_ci_lo", "ci_lo", "lo95", "L_lo95"))
    map_numeric("L_ci_hi", c("L_ci_hi", "ci_hi", "hi95", "L_hi95"))
    map_numeric("L_z", c("L_z", "Z", "z", "Z_L"))
    map_numeric("L_p", c("L_p", "P", "p", "pval", "p_value", "p_value_L"))
    map_numeric("L_expected", c("L_expected", "EZ", "L0"))

    # Uncertainty fields (delta)
    map_numeric("delta_se", c("delta_se", "Delta_se", "delta_Lr_se"))
    map_numeric("delta_var", c("delta_var", "Delta_var"))
    map_numeric("delta_sd", c("delta_sd", "Delta_sd"))
    map_numeric("delta_ci_lo", c("delta_ci_lo", "Delta_ci_lo"))
    map_numeric("delta_ci_hi", c("delta_ci_hi", "Delta_ci_hi"))
    map_numeric("delta_z", c("delta_z", "Delta_z"))
    map_numeric("delta_p", c("delta_p", "Delta_p"))

    # Uncertainty fields (Pearson r)
    map_numeric("r_se", c("Pearson_r_se", "r_se", "Pear_se"))
    map_numeric("r_n", c("Pearson_r_n", "r_n", "n"))

    # Support fields
    map_numeric("support_i_n", c("support_i_n", "gene1_support_n", "gene_i_support_n", "gene1_n"))
    map_numeric("support_j_n", c("support_j_n", "gene2_support_n", "gene_j_support_n", "gene2_n"))
    map_numeric("support_i_prop", c("support_i_prop", "gene1_support_prop", "gene_i_support_prop", "gene1_expr_support"))
    map_numeric("support_j_prop", c("support_j_prop", "gene2_support_prop", "gene_j_support_prop", "gene2_expr_support"))
    map_numeric("support_i_pct", c("support_i_pct", "gene1_expr_pct", "gene1_pct", "gene_i_pct"))
    map_numeric("support_j_pct", c("support_j_pct", "gene2_expr_pct", "gene2_pct", "gene_j_pct"))
    map_numeric("overlap_n", c("overlap_n", "pair_overlap_n", "pair_support_n", "overlap"))
    map_numeric("overlap_ratio", c("overlap_ratio", "overlap_pct", "pair_overlap_ratio", "overlap_prop"))

    # Carry grouping fields if present
    for (nm in c("group", "batch", "donor")) {
        if (nm %in% names(df)) out[[nm]] <- df[[nm]]
    }

    out <- .leesl_compare_attach_sample_info(out, sample_info)

    # Standardize pair ordering and id
    swap <- out$gene_i > out$gene_j
    if (any(swap, na.rm = TRUE)) {
        tmp <- out$gene_i[swap]
        out$gene_i[swap] <- out$gene_j[swap]
        out$gene_j[swap] <- tmp
    }
    out$pair_id <- paste(out$gene_i, out$gene_j, sep = "|")
    out
}

.leesl_compare_normalize_support_pct <- function(df,
                                                 support_pct_range = c("auto", "0-100", "0-1"),
                                                 caller = ".leesl_compare_normalize_support_pct") {
    support_pct_range <- match.arg(support_pct_range)
    prop_cols <- intersect(c("support_i_prop", "support_j_prop"), names(df))
    pct_cols <- intersect(c("support_i_pct", "support_j_pct"), names(df))
    cols <- if (length(prop_cols)) prop_cols else pct_cols
    attr(df, "support_pct_internal_scale") <- "0-1"
    attr(df, "support_pct_legacy_scale") <- "0-100"
    attr(df, "support_pct_requested_scale") <- support_pct_range
    attr(df, "support_pct_normalized_scale") <- "missing"
    attr(df, "support_pct_auto_normalized") <- FALSE
    attr(df, "support_pct_invalid_values") <- FALSE

    if (!length(cols)) {
        attr(df, "support_pct_input_scale") <- NA_character_
        attr(df, "support_pct_scale_used") <- "missing"
        return(df)
    }

    vals_all <- unlist(df[, cols, drop = FALSE], use.names = FALSE)
    vals_observed <- vals_all[!is.na(vals_all)]
    if (!length(vals_observed)) {
        attr(df, "support_pct_input_scale") <- NA_character_
        attr(df, "support_pct_scale_used") <- "missing"
        return(df)
    }

    normalized <- .normalize_support_pct_values(
        vals_all,
        support_pct_scale = support_pct_range,
        caller = caller
    )

    normalized_prop <- matrix(
        normalized$values_prop,
        nrow = nrow(df),
        ncol = length(cols),
        byrow = FALSE,
        dimnames = list(NULL, cols)
    )
    normalized_pct <- matrix(
        normalized$values_pct,
        nrow = nrow(df),
        ncol = length(cols),
        byrow = FALSE,
        dimnames = list(NULL, cols)
    )

    if (length(prop_cols)) {
        for (nm in colnames(normalized_prop)) {
            df[[nm]] <- normalized_prop[, nm]
        }
        if ("support_i_prop" %in% names(df)) df$support_i_pct <- normalized_pct[, "support_i_prop"]
        if ("support_j_prop" %in% names(df)) df$support_j_pct <- normalized_pct[, "support_j_prop"]
    } else {
        if ("support_i_pct" %in% names(df)) df$support_i_pct <- normalized_pct[, "support_i_pct"]
        if ("support_j_pct" %in% names(df)) df$support_j_pct <- normalized_pct[, "support_j_pct"]
        df$support_i_prop <- normalized_prop[, "support_i_pct"]
        df$support_j_prop <- normalized_prop[, "support_j_pct"]
    }

    attr(df, "support_pct_input_scale") <- normalized$input_scale
    attr(df, "support_pct_scale_used") <- normalized$scale_used
    attr(df, "support_pct_normalized_scale") <- normalized$normalized_scale
    attr(df, "support_pct_auto_normalized") <- normalized$auto_normalized
    attr(df, "support_pct_invalid_values") <- normalized$invalid_values
    df
}

.leesl_compare_get_lee_stats <- function(scope_obj, grid_name = NULL, lee_layer = NULL) {
    g_layer <- .select_grid_layer(scope_obj, grid_name)
    if (is.null(grid_name)) {
        grid_name <- names(scope_obj@grid)[
            vapply(scope_obj@grid, identical, logical(1), g_layer)
        ]
    }
    if (is.null(lee_layer)) {
        cand <- character(0)
        if (!is.null(scope_obj@stats[[grid_name]])) cand <- names(scope_obj@stats[[grid_name]])
        cand <- c(cand, names(g_layer))
        cand <- unique(cand[grepl("^LeeStats_", cand)])
        if (length(cand) == 1L) {
            lee_layer <- cand
        } else if ("LeeStats_Xz" %in% cand) {
            lee_layer <- "LeeStats_Xz"
        } else if (length(cand)) {
            lee_layer <- cand[1]
        }
    }

    leeStat <- NULL
    if (!is.null(lee_layer) && !is.null(scope_obj@stats[[grid_name]]) &&
        !is.null(scope_obj@stats[[grid_name]][[lee_layer]])) {
        leeStat <- scope_obj@stats[[grid_name]][[lee_layer]]
    }
    if (is.null(leeStat) && !is.null(lee_layer) && !is.null(g_layer[[lee_layer]])) {
        leeStat <- g_layer[[lee_layer]]
    }
    list(leeStat = leeStat, grid_name = grid_name, lee_layer = lee_layer, g_layer = g_layer)
}

.leesl_compare_pick_expr_matrix <- function(g_layer, pear_level = c("cell", "grid"), scope_obj = NULL) {
    pear_level <- match.arg(pear_level)
    pick_matrix <- function(layer) {
        if (is.null(layer)) return(NULL)
        for (nm in c("counts", "raw_counts", "expr", "X", "data", "logCPM", "Xz")) {
            m <- layer[[nm]]
            if (!is.null(m) && (is.matrix(m) || inherits(m, "dgCMatrix"))) return(m)
        }
        NULL
    }

    expr_mat <- pick_matrix(g_layer)
    if (is.null(expr_mat) && pear_level == "cell" && !is.null(scope_obj)) {
        cell_env <- tryCatch(scope_obj@cell, error = function(e) NULL)
        if (is.null(cell_env)) cell_env <- tryCatch(scope_obj@cells, error = function(e) NULL)
        expr_mat <- pick_matrix(cell_env)
    }
    expr_mat
}

.leesl_compare_build_expr_bin <- function(g_layer, genes_use, pear_level = c("cell", "grid"), scope_obj = NULL) {
    pear_level <- match.arg(pear_level)
    n_grid <- NA_integer_
    source <- NA_character_
    expr_bin <- NULL

    if (!is.null(g_layer$counts) && is.data.frame(g_layer$counts) &&
        all(c("gene", "grid_id") %in% colnames(g_layer$counts))) {
        ct <- g_layer$counts
        keep <- ct$gene %in% genes_use
        ct <- ct[keep, , drop = FALSE]
        if (nrow(ct)) {
            grid_ids <- if (!is.null(g_layer$grid_info$grid_id)) g_layer$grid_info$grid_id else unique(ct$grid_id)
            n_grid <- length(grid_ids)
            i <- match(ct$grid_id, grid_ids)
            j <- match(ct$gene, genes_use)
            ok <- !is.na(i) & !is.na(j)
            if (any(ok)) {
                expr_bin <- Matrix::sparseMatrix(
                    i = i[ok],
                    j = j[ok],
                    x = 1,
                    dims = c(n_grid, length(genes_use))
                )
                source <- "counts_table"
            }
        }
    }

    if (is.null(expr_bin)) {
        expr_mat <- .leesl_compare_pick_expr_matrix(g_layer, pear_level = pear_level, scope_obj = scope_obj)
        if (!is.null(expr_mat)) {
            rn <- rownames(expr_mat)
            cn <- colnames(expr_mat)
            inter_col <- intersect(cn, genes_use)
            inter_row <- intersect(rn, genes_use)
            gene_in_col <- length(inter_col) >= length(inter_row)
            if (!gene_in_col) {
                expr_mat <- if (inherits(expr_mat, "dgCMatrix")) Matrix::t(expr_mat) else t(expr_mat)
                cn <- colnames(expr_mat)
                inter_col <- intersect(cn, genes_use)
            }
            if (length(inter_col)) {
                expr_mat <- expr_mat[, inter_col, drop = FALSE]
                expr_bin <- expr_mat > 0
                n_grid <- nrow(expr_mat)
                genes_use <- inter_col
                source <- "expr_matrix"
            }
        }
    }

    list(expr_bin = expr_bin, n_grid = n_grid, genes_use = genes_use, source = source)
}

.leesl_compare_compute_support_overlap <- function(expr_bin, genes_use, pairs_df, n_grid) {
    support_i_n <- support_j_n <- support_i_prop <- support_j_prop <- support_i_pct <- support_j_pct <- overlap_n <- overlap_ratio <- rep(NA_real_, nrow(pairs_df))
    if (is.null(expr_bin) || !length(genes_use) || is.null(n_grid) || n_grid <= 0) {
        return(list(
            support_i_n = support_i_n,
            support_j_n = support_j_n,
            support_i_prop = support_i_prop,
            support_j_prop = support_j_prop,
            support_i_pct = support_i_pct,
            support_j_pct = support_j_pct,
            overlap_n = overlap_n,
            overlap_ratio = overlap_ratio
        ))
    }

    if (inherits(expr_bin, "dgCMatrix")) {
        support_counts <- Matrix::colSums(expr_bin)
    } else {
        support_counts <- colSums(expr_bin)
    }
    support_counts <- as.numeric(support_counts)
    names(support_counts) <- genes_use

    support_i_n <- support_counts[pairs_df$gene_i]
    support_j_n <- support_counts[pairs_df$gene_j]
    support_i_prop <- support_i_n / n_grid
    support_j_prop <- support_j_n / n_grid
    support_i_pct <- support_i_prop * 100
    support_j_pct <- support_j_prop * 100

    pair_i <- match(pairs_df$gene_i, genes_use)
    pair_j <- match(pairs_df$gene_j, genes_use)

    if (length(genes_use) <= 2000 && nrow(pairs_df) <= 2e5) {
        cross <- crossprod(expr_bin)
        overlap_n <- as.numeric(cross[cbind(pair_i, pair_j)])
    } else {
        overlap_n <- vapply(seq_len(nrow(pairs_df)), function(k) {
            i <- pair_i[k]
            j <- pair_j[k]
            if (is.na(i) || is.na(j)) return(NA_real_)
            if (inherits(expr_bin, "dgCMatrix")) {
                sum(expr_bin[, i] & expr_bin[, j])
            } else {
                sum(expr_bin[, i] & expr_bin[, j])
            }
        }, numeric(1))
    }

    denom <- pmin(support_i_n, support_j_n)
    overlap_ratio <- ifelse(is.na(denom) | denom == 0, NA_real_, overlap_n / denom)

    list(
        support_i_n = support_i_n,
        support_j_n = support_j_n,
        support_i_prop = support_i_prop,
        support_j_prop = support_j_prop,
        support_i_pct = support_i_pct,
        support_j_pct = support_j_pct,
        overlap_n = overlap_n,
        overlap_ratio = overlap_ratio
    )
}

.leesl_compare_extract_from_scope <- function(scope_obj,
                                              sample_id,
                                              grid_name = NULL,
                                              lee_stats_layer = NULL,
                                              pear_level = c("cell", "grid"),
                                              pairs = NULL,
                                              verbose = TRUE) {
    pear_level <- match.arg(pear_level)
    stat <- .leesl_compare_get_lee_stats(scope_obj, grid_name = grid_name, lee_layer = lee_stats_layer)
    leeStat <- stat$leeStat
    g_layer <- stat$g_layer
    grid_name <- stat$grid_name

    Lmat <- .get_lee_matrix(scope_obj, grid_name = grid_name, lee_layer = stat$lee_layer)
    genes <- rownames(Lmat)
    if (is.null(genes)) genes <- colnames(Lmat)
    if (is.null(genes)) stop("Lee's L matrix must have row/column names")

    pairs_df <- .leesl_compare_standardize_pairs(pairs)
    if (is.null(pairs_df) || !nrow(pairs_df)) {
        if (inherits(Lmat, "big.matrix")) {
            if (verbose) message("[compareLeesL] Converting big.matrix to matrix for full pair table; consider providing pairs to limit.")
            Lmat_use <- as.matrix(Lmat)
        } else {
            Lmat_use <- Lmat
        }
        ut <- upper.tri(Lmat_use)
        gene_i <- rep(genes, each = length(genes))[ut]
        gene_j <- rep(genes, length(genes))[ut]
        L_vals <- Lmat_use[ut]
        pairs_df <- data.frame(gene_i = gene_i, gene_j = gene_j, stringsAsFactors = FALSE)
        pairs_df$pair_id <- paste(pairs_df$gene_i, pairs_df$gene_j, sep = "|")
        idx_i <- match(gene_i, genes)
        idx_j <- match(gene_j, genes)
    } else {
        idx_i <- match(pairs_df$gene_i, genes)
        idx_j <- match(pairs_df$gene_j, genes)
    }

    missing_pair <- is.na(idx_i) | is.na(idx_j)
    L_vals <- rep(NA_real_, nrow(pairs_df))
    if (any(!missing_pair)) {
        if (inherits(Lmat, "big.matrix")) {
            L_vals[!missing_pair] <- Lmat[cbind(idx_i[!missing_pair], idx_j[!missing_pair])]
        } else {
            L_vals[!missing_pair] <- Lmat[cbind(idx_i[!missing_pair], idx_j[!missing_pair])]
        }
    }

    r_vals <- rep(NA_real_, nrow(pairs_df))
    rmat <- tryCatch(.get_pearson_matrix(scope_obj, grid_name = grid_name, level = pear_level),
        error = function(e) NULL
    )
    if (!is.null(rmat)) {
        common <- intersect(rownames(rmat), genes)
        if (length(common)) {
            rmat <- rmat[genes, genes, drop = FALSE]
            if (any(!missing_pair)) {
                r_vals[!missing_pair] <- rmat[cbind(idx_i[!missing_pair], idx_j[!missing_pair])]
            }
        }
    }

    Z_vals <- P_vals <- rep(NA_real_, nrow(pairs_df))
    if (!is.null(leeStat) && is.list(leeStat)) {
        if (!is.null(leeStat$Z)) {
            Zmat <- leeStat$Z
            if (any(!missing_pair)) {
                Z_vals[!missing_pair] <- Zmat[cbind(idx_i[!missing_pair], idx_j[!missing_pair])]
            }
        }
        if (!is.null(leeStat$P)) {
            Pmat <- leeStat$P
            if (any(!missing_pair)) {
                P_vals[!missing_pair] <- Pmat[cbind(idx_i[!missing_pair], idx_j[!missing_pair])]
            }
        }
    }

    out <- data.frame(
        sample_id = sample_id,
        gene_i = pairs_df$gene_i,
        gene_j = pairs_df$gene_j,
        pair_id = pairs_df$pair_id,
        L = L_vals,
        Pearson_r = r_vals,
        delta_Lr = ifelse(is.na(L_vals) | is.na(r_vals), NA_real_, L_vals - r_vals),
        L_z = Z_vals,
        L_p = P_vals,
        gene_missing = missing_pair,
        stringsAsFactors = FALSE
    )

    # Support and overlap
    genes_use <- unique(c(out$gene_i, out$gene_j))
    expr_info <- .leesl_compare_build_expr_bin(g_layer, genes_use, pear_level = pear_level, scope_obj = scope_obj)
    support <- .leesl_compare_compute_support_overlap(expr_info$expr_bin, expr_info$genes_use, out, expr_info$n_grid)
    out$support_i_n <- support$support_i_n
    out$support_j_n <- support$support_j_n
    out$support_i_prop <- support$support_i_prop
    out$support_j_prop <- support$support_j_prop
    out$support_i_pct <- support$support_i_pct
    out$support_j_pct <- support$support_j_pct
    out$overlap_n <- support$overlap_n
    out$overlap_ratio <- support$overlap_ratio

    # Sample-level metadata for uncertainty
    n_grid <- expr_info$n_grid
    if (is.na(n_grid) && !is.null(g_layer$grid_info)) n_grid <- nrow(g_layer$grid_info)
    W <- g_layer$W
    S0 <- if (!is.null(W)) sum(W) else NA_real_
    EZ <- if (!is.na(n_grid) && n_grid > 1) -1 / (n_grid - 1) else NA_real_
    Var <- if (!is.na(n_grid) && n_grid > 3 && !is.na(S0) && S0 > 0) {
        (n_grid^2 * (n_grid - 2)) / ((n_grid - 1)^2 * (n_grid - 3) * S0)
    } else {
        NA_real_
    }

    sample_meta <- data.frame(
        sample_id = sample_id,
        n_grid = n_grid,
        S0 = S0,
        EZ = EZ,
        Var = Var,
        support_source = expr_info$source,
        stringsAsFactors = FALSE
    )

    qc <- NULL
    if (!is.null(leeStat) && is.list(leeStat) && !is.null(leeStat$qc)) qc <- leeStat$qc
    sample_qc <- data.frame(sample_id = sample_id, stringsAsFactors = FALSE)
    if (is.list(qc)) {
        for (nm in c("edge_density", "components", "modularity_Q", "mean_degree", "sd_degree", "hub_ratio", "sig_edge_frac")) {
            if (!is.null(qc[[nm]])) sample_qc[[nm]] <- as.numeric(qc[[nm]])
        }
    }

    list(pairs = out, sample_meta = sample_meta, sample_qc = sample_qc)
}

.leesl_compare_infer_uncertainty <- function(df, sample_meta = NULL, epsilon = 1e-6) {
    if (!is.data.frame(df) || !nrow(df)) return(df)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    n <- nrow(df)

    # Attach sample-level EZ/Var if present
    if (!is.null(sample_meta) && is.data.frame(sample_meta)) {
        if ("sample_id" %in% names(sample_meta)) {
            add_cols <- setdiff(names(sample_meta), names(df))
            if (length(add_cols)) {
                sample_meta_sub <- sample_meta[, c("sample_id", add_cols), drop = FALSE]
                df <- merge(df, sample_meta_sub, by = "sample_id", all.x = TRUE, sort = FALSE)
            }
        }
    }

    L <- df$L
    L_expected <- if (!is.null(df$L_expected)) df$L_expected else df$EZ

    se <- rep(NA_real_, n)
    src <- rep(NA_character_, n)

    fill <- function(cond, value, name) {
        idx <- which(is.na(se) & cond)
        if (length(idx)) {
            se[idx] <<- value[idx]
            src[idx] <<- name
        }
    }

    if (!is.null(df$L_se)) fill(!is.na(df$L_se), df$L_se, "se")
    if (!is.null(df$L_var)) fill(!is.na(df$L_var), sqrt(df$L_var), "var")
    if (!is.null(df$L_sd)) fill(!is.na(df$L_sd), df$L_sd, "sd")

    ci_width <- rep(NA_real_, n)
    if (!is.null(df$L_ci_lo) && !is.null(df$L_ci_hi)) {
        ci_width <- df$L_ci_hi - df$L_ci_lo
        se_ci <- ci_width / (2 * 1.96)
        fill(!is.na(se_ci), se_ci, "ci95")
    }

    if (!is.null(df$L_z)) {
        z <- df$L_z
        z_ok <- !is.na(z) & abs(z) > 0
        L0 <- ifelse(is.na(L_expected), 0, L_expected)
        se_z <- abs((L - L0) / z)
        fill(z_ok, se_z, "z")
    }

    if (!is.null(df$L_p)) {
        p <- df$L_p
        z <- stats::qnorm(1 - p / 2)
        z_ok <- !is.na(z) & is.finite(z) & abs(z) > 0
        L0 <- ifelse(is.na(L_expected), 0, L_expected)
        se_p <- abs((L - L0) / z)
        fill(z_ok, se_p, "p")
    }

    if (!is.null(df$Var)) {
        fill(!is.na(df$Var) & df$Var > 0, sqrt(df$Var), "analytic_var")
    }

    df$L_se <- se
    df$L_se_source <- src
    df$L_ci_width <- ci_width

    # Pearson r uncertainty (for delta fallback)
    r_se <- if (!is.null(df$r_se)) df$r_se else rep(NA_real_, n)
    if (all(is.na(r_se)) && !is.null(df$r_n) && !is.null(df$Pearson_r)) {
        n_use <- df$r_n
        r <- df$Pearson_r
        ok <- !is.na(n_use) & n_use > 3 & !is.na(r)
        if (any(ok)) {
            se_z <- 1 / sqrt(n_use[ok] - 3)
            r_se[ok] <- se_z * (1 - r[ok]^2)
        }
    }
    df$r_se <- r_se

    # Delta uncertainty
    delta <- df$delta_Lr
    dse <- rep(NA_real_, n)
    dsrc <- rep(NA_character_, n)

    fill_delta <- function(cond, value, name) {
        idx <- which(is.na(dse) & cond)
        if (length(idx)) {
            dse[idx] <<- value[idx]
            dsrc[idx] <<- name
        }
    }

    if (!is.null(df$delta_se)) fill_delta(!is.na(df$delta_se), df$delta_se, "se")
    if (!is.null(df$delta_var)) fill_delta(!is.na(df$delta_var), sqrt(df$delta_var), "var")
    if (!is.null(df$delta_sd)) fill_delta(!is.na(df$delta_sd), df$delta_sd, "sd")
    if (!is.null(df$delta_ci_lo) && !is.null(df$delta_ci_hi)) {
        dwidth <- df$delta_ci_hi - df$delta_ci_lo
        fill_delta(!is.na(dwidth), dwidth / (2 * 1.96), "ci95")
    }
    if (!is.null(df$delta_z)) {
        z <- df$delta_z
        z_ok <- !is.na(z) & abs(z) > 0
        dse_z <- abs(delta / z)
        fill_delta(z_ok, dse_z, "z")
    }
    if (!is.null(df$delta_p)) {
        p <- df$delta_p
        z <- stats::qnorm(1 - p / 2)
        z_ok <- !is.na(z) & is.finite(z) & abs(z) > 0
        dse_p <- abs(delta / z)
        fill_delta(z_ok, dse_p, "p")
    }

    if (all(is.na(dse)) && !all(is.na(df$L_se))) {
        dse <- df$L_se
        dsrc <- ifelse(is.na(df$L_se), NA_character_, "L_se_only")
    } else if (!all(is.na(df$L_se)) && !all(is.na(df$r_se))) {
        comb <- sqrt(df$L_se^2 + df$r_se^2)
        fill_delta(!is.na(comb), comb, "L_r_se")
    }

    df$delta_se <- dse
    df$delta_se_source <- dsrc

    df
}

.leesl_compare_apply_gate <- function(df,
                                      min_support_n = 5,
                                      min_support_prop = 0,
                                      min_overlap_n = 3,
                                      min_overlap_ratio = 0.1,
                                      max_rel_se = Inf,
                                      max_ci_width = Inf,
                                      require_uncertainty = TRUE,
                                      epsilon = 1e-6) {
    if (!is.data.frame(df) || !nrow(df)) return(df)

    if (!"gene_missing" %in% names(df)) df$gene_missing <- FALSE
    if (!"support_i_n" %in% names(df)) df$support_i_n <- NA_real_
    if (!"support_j_n" %in% names(df)) df$support_j_n <- NA_real_
    if (!"support_i_prop" %in% names(df) && "support_i_pct" %in% names(df)) df$support_i_prop <- df$support_i_pct / 100
    if (!"support_j_prop" %in% names(df) && "support_j_pct" %in% names(df)) df$support_j_prop <- df$support_j_pct / 100
    if (!"support_i_prop" %in% names(df)) df$support_i_prop <- NA_real_
    if (!"support_j_prop" %in% names(df)) df$support_j_prop <- NA_real_
    if (!"support_i_pct" %in% names(df)) df$support_i_pct <- NA_real_
    if (!"support_j_pct" %in% names(df)) df$support_j_pct <- NA_real_
    if (!"overlap_n" %in% names(df)) df$overlap_n <- NA_real_
    if (!"overlap_ratio" %in% names(df)) df$overlap_ratio <- NA_real_
    if (!"L_ci_width" %in% names(df)) df$L_ci_width <- NA_real_
    if (!"L_se" %in% names(df)) df$L_se <- NA_real_
    if (!"delta_se" %in% names(df)) df$delta_se <- NA_real_

    support_missing <- is.na(df$support_i_n) | is.na(df$support_j_n)
    support_pct_missing <- is.na(df$support_i_prop) | is.na(df$support_j_prop)
    support_pct_invalid <- (!support_pct_missing) &
        ((df$support_i_prop < 0) | (df$support_i_prop > 1) |
            (df$support_j_prop < 0) | (df$support_j_prop > 1))
    overlap_missing <- is.na(df$overlap_n) | is.na(df$overlap_ratio)

    support_ok <- (!support_missing) & (df$support_i_n >= min_support_n) & (df$support_j_n >= min_support_n)
    if (isTRUE(min_support_prop > 0)) {
        support_pct_ok <- (!support_pct_missing) & (!support_pct_invalid) &
            (df$support_i_prop >= min_support_prop) &
            (df$support_j_prop >= min_support_prop)
        support_ok <- support_ok & support_pct_ok
    }
    overlap_ok <- (!overlap_missing) & (df$overlap_n >= min_overlap_n) & (df$overlap_ratio >= min_overlap_ratio)

    rel_se <- ifelse(is.na(df$L_se), NA_real_, df$L_se / (abs(df$L) + epsilon))
    stability_ok <- TRUE
    if (is.finite(max_rel_se)) stability_ok <- stability_ok & (is.na(rel_se) | rel_se <= max_rel_se)
    if (is.finite(max_ci_width)) stability_ok <- stability_ok & (is.na(df$L_ci_width) | df$L_ci_width <= max_ci_width)

    uncertainty_ok_L <- if (!require_uncertainty) TRUE else !is.na(df$L_se)
    uncertainty_ok_delta <- if (!require_uncertainty) TRUE else !is.na(df$delta_se)

    base_ok <- support_ok & overlap_ok & stability_ok & !df$gene_missing

    df$gate_L <- base_ok & uncertainty_ok_L
    df$gate_delta <- base_ok & uncertainty_ok_delta

    reason_base <- rep(NA_character_, nrow(df))
    reason_base[support_missing] <- "missing_support"
    reason_base[!support_missing & !support_ok] <- "low_support"
    if (isTRUE(min_support_prop > 0)) {
        reason_base[support_pct_missing] <- ifelse(
            is.na(reason_base[support_pct_missing]),
            "missing_support_pct",
            paste(reason_base[support_pct_missing], "missing_support_pct", sep = ";")
        )
        reason_base[support_pct_invalid] <- ifelse(
            is.na(reason_base[support_pct_invalid]),
            "invalid_support_pct",
            paste(reason_base[support_pct_invalid], "invalid_support_pct", sep = ";")
        )
        bad_support_pct <- !support_pct_missing &
            !support_pct_invalid &
            ((df$support_i_prop < min_support_prop) | (df$support_j_prop < min_support_prop))
        reason_base[bad_support_pct] <- ifelse(
            is.na(reason_base[bad_support_pct]),
            "low_support_pct",
            paste(reason_base[bad_support_pct], "low_support_pct", sep = ";")
        )
    }
    reason_base[overlap_missing] <- ifelse(is.na(reason_base[overlap_missing]), "missing_overlap", paste(reason_base[overlap_missing], "missing_overlap", sep = ";"))
    reason_base[!overlap_missing & !overlap_ok] <- ifelse(is.na(reason_base[!overlap_missing & !overlap_ok]), "low_overlap", paste(reason_base[!overlap_missing & !overlap_ok], "low_overlap", sep = ";"))
    reason_base[df$gene_missing] <- ifelse(is.na(reason_base[df$gene_missing]), "missing_gene", paste(reason_base[df$gene_missing], "missing_gene", sep = ";"))
    if (is.finite(max_rel_se)) {
        bad_rel <- !is.na(rel_se) & rel_se > max_rel_se
        reason_base[bad_rel] <- ifelse(is.na(reason_base[bad_rel]), "high_rel_se", paste(reason_base[bad_rel], "high_rel_se", sep = ";"))
    }
    if (is.finite(max_ci_width)) {
        bad_ci <- !is.na(df$L_ci_width) & df$L_ci_width > max_ci_width
        reason_base[bad_ci] <- ifelse(is.na(reason_base[bad_ci]), "wide_ci", paste(reason_base[bad_ci], "wide_ci", sep = ";"))
    }

    reason_L <- reason_base
    if (require_uncertainty) {
        miss_L <- is.na(df$L_se)
        reason_L[miss_L] <- ifelse(is.na(reason_L[miss_L]), "missing_uncertainty", paste(reason_L[miss_L], "missing_uncertainty", sep = ";"))
    }

    reason_delta <- reason_base
    if (require_uncertainty) {
        miss_delta <- is.na(df$delta_se)
        reason_delta[miss_delta] <- ifelse(is.na(reason_delta[miss_delta]), "missing_uncertainty", paste(reason_delta[miss_delta], "missing_uncertainty", sep = ";"))
    }

    df$gate_reason_L <- ifelse(df$gate_L, NA_character_, reason_L)
    df$gate_reason_delta <- ifelse(df$gate_delta, NA_character_, reason_delta)
    df
}

.leesl_compare_build_background <- function(sample_qc, sample_info = NULL) {
    if (!is.data.frame(sample_qc) || !nrow(sample_qc)) {
        sample_qc <- data.frame(sample_id = character(0), stringsAsFactors = FALSE)
    }
    sample_qc <- .leesl_compare_attach_sample_info(sample_qc, sample_info)
    if (!"modularity" %in% names(sample_qc) && "modularity_Q" %in% names(sample_qc)) {
        sample_qc$modularity <- sample_qc$modularity_Q
    }

    metrics <- c("edge_density", "components", "modularity", "mean_degree", "sd_degree", "hub_ratio", "sig_edge_frac")
    present <- intersect(metrics, names(sample_qc))

    score <- rep(NA_real_, nrow(sample_qc))
    if (length(present)) {
        mat <- as.matrix(sample_qc[, present, drop = FALSE])
        # invert metrics where lower is better
        inv <- colnames(mat) %in% c("components", "sd_degree")
        mat[, inv] <- -mat[, inv]
        z <- apply(mat, 2, function(v) {
            if (all(is.na(v))) return(rep(NA_real_, length(v)))
            s <- stats::sd(v, na.rm = TRUE)
            if (is.na(s) || s == 0) return(rep(0, length(v)))
            (v - mean(v, na.rm = TRUE)) / s
        })
        if (is.vector(z)) z <- matrix(z, ncol = 1)
        score <- rowMeans(z, na.rm = TRUE)
    }

    sample_qc$background_score <- score
    label <- rep(NA_character_, length(score))
    good <- !is.na(score)
    if (sum(good) == 1) {
        label[good] <- "single"
    } else if (sum(good) == 2) {
        med <- stats::median(score[good])
        label[good] <- ifelse(score[good] <= med, "low", "high")
    } else if (sum(good) >= 3) {
        qs <- stats::quantile(score[good], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
        qs <- unique(qs)
        if (length(qs) >= 2) {
            label[good] <- cut(score[good], breaks = qs, include.lowest = TRUE, labels = FALSE)
            label[good] <- c("low", "mid", "high")[label[good]]
        }
    }
    sample_qc$background_label <- label
    sample_qc
}

.leesl_compare_pairwise <- function(df, sample_info = NULL, background = NULL, p_adj_mode = c("BH", "none")) {
    p_adj_mode <- match.arg(p_adj_mode)
    if (!is.data.frame(df) || !nrow(df)) return(data.frame())

    samples <- unique(df$sample_id)
    if (length(samples) < 2) return(data.frame())

    out_list <- list()
    combs <- combn(samples, 2, simplify = FALSE)

    for (pair in combs) {
        s1 <- pair[1]
        s2 <- pair[2]
        d1 <- df[df$sample_id == s1, , drop = FALSE]
        d2 <- df[df$sample_id == s2, , drop = FALSE]
        merged <- merge(d1, d2, by = "pair_id", suffixes = c("_1", "_2"), all = FALSE, sort = FALSE)
        if (!nrow(merged)) next
        merged$sample_id_1 <- s1
        merged$sample_id_2 <- s2
        merged$gene_i <- merged$gene_i_1
        merged$gene_j <- merged$gene_j_1

        merged$comparable_L <- merged$gate_L_1 & merged$gate_L_2
        merged$comparable_delta <- merged$gate_delta_1 & merged$gate_delta_2

        merged$diff_L <- merged$L_2 - merged$L_1
        merged$diff_delta <- merged$delta_Lr_2 - merged$delta_Lr_1

        merged$se_diff_L <- sqrt(merged$L_se_1^2 + merged$L_se_2^2)
        merged$se_diff_delta <- sqrt(merged$delta_se_1^2 + merged$delta_se_2^2)

        bad_se_L <- !is.na(merged$se_diff_L) & merged$se_diff_L <= 0
        bad_se_delta <- !is.na(merged$se_diff_delta) & merged$se_diff_delta <= 0
        merged$z_diff_L <- ifelse(merged$comparable_L & !bad_se_L, merged$diff_L / merged$se_diff_L, NA_real_)
        merged$z_diff_delta <- ifelse(merged$comparable_delta & !bad_se_delta, merged$diff_delta / merged$se_diff_delta, NA_real_)
        merged$diff_issue_L <- ifelse(bad_se_L, "nonpositive_se", NA_character_)
        merged$diff_issue_delta <- ifelse(bad_se_delta, "nonpositive_se", NA_character_)

        merged$p_diff_L <- ifelse(is.na(merged$z_diff_L), NA_real_, 2 * stats::pnorm(-abs(merged$z_diff_L)))
        merged$p_diff_delta <- ifelse(is.na(merged$z_diff_delta), NA_real_, 2 * stats::pnorm(-abs(merged$z_diff_delta)))

        if (p_adj_mode != "none") {
            merged$p_adj_L <- stats::p.adjust(merged$p_diff_L, method = p_adj_mode)
            merged$p_adj_delta <- stats::p.adjust(merged$p_diff_delta, method = p_adj_mode)
        } else {
            merged$p_adj_L <- NA_real_
            merged$p_adj_delta <- NA_real_
        }

        # Gate reasons
        merged$gate_reason_L <- ifelse(merged$comparable_L, NA_character_,
            paste0("s1:", merged$gate_reason_L_1, "|s2:", merged$gate_reason_L_2)
        )
        merged$gate_reason_delta <- ifelse(merged$comparable_delta, NA_character_,
            paste0("s1:", merged$gate_reason_delta_1, "|s2:", merged$gate_reason_delta_2)
        )

        # Attach sample info and background
        if (!is.null(sample_info) && is.data.frame(sample_info) && nrow(sample_info)) {
            info <- unique(sample_info)
            colnames(info) <- paste0(colnames(info), "_1")
            merged <- merge(merged, info, by.x = "sample_id_1", by.y = "sample_id_1", all.x = TRUE, sort = FALSE)
            info2 <- unique(sample_info)
            colnames(info2) <- paste0(colnames(info2), "_2")
            merged <- merge(merged, info2, by.x = "sample_id_2", by.y = "sample_id_2", all.x = TRUE, sort = FALSE)
        }
        if (!is.null(background) && is.data.frame(background) && nrow(background)) {
            bg1 <- background
            colnames(bg1) <- paste0(colnames(bg1), "_1")
            merged <- merge(merged, bg1, by.x = "sample_id_1", by.y = "sample_id_1", all.x = TRUE, sort = FALSE)
            bg2 <- background
            colnames(bg2) <- paste0(colnames(bg2), "_2")
            merged <- merge(merged, bg2, by.x = "sample_id_2", by.y = "sample_id_2", all.x = TRUE, sort = FALSE)
        }

        out_list[[length(out_list) + 1]] <- merged
    }

    if (!length(out_list)) return(data.frame())
    do.call(rbind, out_list)
}

.leesl_compare_group <- function(df,
                                 group_col = "group",
                                 min_samples_per_group = 1,
                                 p_adj_mode = c("BH", "none")) {
    p_adj_mode <- match.arg(p_adj_mode)
    if (!is.data.frame(df) || !nrow(df) || !group_col %in% names(df)) return(data.frame())

    df$group <- df[[group_col]]
    df <- df[!is.na(df$group) & nzchar(df$group), , drop = FALSE]
    if (!nrow(df)) return(data.frame())

    groups <- unique(df$group)
    if (length(groups) < 2) return(data.frame())

    # Aggregate within group using precision weights
    agg <- function(metric, se_col, gate_col) {
        by_key <- interaction(df$group, df$pair_id, drop = TRUE)
        rows <- lapply(split(df, by_key), function(sub) {
            sub <- sub[sub[[gate_col]], , drop = FALSE]
            n_use <- length(unique(sub$sample_id))
            if (!n_use || n_use < min_samples_per_group) return(NULL)
            val <- sub[[metric]]
            se <- sub[[se_col]]
            w <- ifelse(is.na(se) | se <= 0, NA_real_, 1 / se^2)
            use_w <- !is.na(w)
            if (any(use_w)) {
                mean_val <- sum(val[use_w] * w[use_w]) / sum(w[use_w])
                se_val <- sqrt(1 / sum(w[use_w]))
                w_source <- "precision"
            } else {
                mean_val <- mean(val, na.rm = TRUE)
                se_val <- NA_real_
                w_source <- "unweighted"
            }
            data.frame(
                group = sub$group[1],
                pair_id = sub$pair_id[1],
                gene_i = sub$gene_i[1],
                gene_j = sub$gene_j[1],
                mean_val = mean_val,
                se_val = se_val,
                n_samples = n_use,
                weight_source = w_source,
                stringsAsFactors = FALSE
            )
        })
        rows <- Filter(Negate(is.null), rows)
        if (!length(rows)) return(data.frame())
        do.call(rbind, rows)
    }

    L_agg <- agg("L", "L_se", "gate_L")
    delta_agg <- agg("delta_Lr", "delta_se", "gate_delta")

    if (!nrow(L_agg) && !nrow(delta_agg)) return(data.frame())

    out_list <- list()
    combs <- combn(groups, 2, simplify = FALSE)

    for (pair in combs) {
        g1 <- pair[1]
        g2 <- pair[2]
        L1 <- L_agg[L_agg$group == g1, , drop = FALSE]
        L2 <- L_agg[L_agg$group == g2, , drop = FALSE]
        merged <- merge(L1, L2, by = "pair_id", suffixes = c("_1", "_2"), all = FALSE, sort = FALSE)
        if (nrow(merged)) {
            merged$group_1 <- g1
            merged$group_2 <- g2
            merged$diff_L <- merged$mean_val_2 - merged$mean_val_1
            merged$se_diff_L <- sqrt(merged$se_val_1^2 + merged$se_val_2^2)
            merged$z_diff_L <- merged$diff_L / merged$se_diff_L
            merged$p_diff_L <- 2 * stats::pnorm(-abs(merged$z_diff_L))
            if (p_adj_mode != "none") {
                merged$p_adj_L <- stats::p.adjust(merged$p_diff_L, method = p_adj_mode)
            } else {
                merged$p_adj_L <- NA_real_
            }
        }

        D1 <- delta_agg[delta_agg$group == g1, , drop = FALSE]
        D2 <- delta_agg[delta_agg$group == g2, , drop = FALSE]
        merged_delta <- merge(D1, D2, by = "pair_id", suffixes = c("_1", "_2"), all = FALSE, sort = FALSE)
        if (nrow(merged_delta)) {
            merged_delta$group_1 <- g1
            merged_delta$group_2 <- g2
            merged_delta$diff_delta <- merged_delta$mean_val_2 - merged_delta$mean_val_1
            merged_delta$se_diff_delta <- sqrt(merged_delta$se_val_1^2 + merged_delta$se_val_2^2)
            merged_delta$z_diff_delta <- merged_delta$diff_delta / merged_delta$se_diff_delta
            merged_delta$p_diff_delta <- 2 * stats::pnorm(-abs(merged_delta$z_diff_delta))
            if (p_adj_mode != "none") {
                merged_delta$p_adj_delta <- stats::p.adjust(merged_delta$p_diff_delta, method = p_adj_mode)
            } else {
                merged_delta$p_adj_delta <- NA_real_
            }
        }

        if (nrow(merged) || nrow(merged_delta)) {
            out_list[[length(out_list) + 1]] <- list(L = merged, delta = merged_delta)
        }
    }

    if (!length(out_list)) return(data.frame())
    # Keep L and delta tables separate but return combined with type tag
    rows <- list()
    for (entry in out_list) {
        if (nrow(entry$L)) {
            tmp <- entry$L
            tmp$metric <- "L"
            rows[[length(rows) + 1]] <- tmp
        }
        if (nrow(entry$delta)) {
            tmp <- entry$delta
            tmp$metric <- "delta_Lr"
            rows[[length(rows) + 1]] <- tmp
        }
    }
    if (!length(rows)) return(data.frame())
    do.call(rbind, rows)
}
