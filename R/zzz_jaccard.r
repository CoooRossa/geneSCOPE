

#' Compute gene–gene Jaccard similarity and build a heatmap object
#'
#' Binarizes each gene's spatial footprint (counts or Xz > 0), computes the
#' gene–gene Jaccard matrix, and constructs a ComplexHeatmap object. No files
#' are written; results are returned to the caller.
#'
#' @param scope_obj geneSCOPE object containing the grid layers.
#' @param grid_name Grid layer name (default "grid30").
#' @param expr_layer Expression layer to binarize, one of c("Xz", "counts").
#' @param min_bins_per_gene Minimum bins per gene to retain (default 8).
#' @param top_n_genes Optional cap on number of genes (by bin presence); NULL keeps all.
#' @param genes_of_interest Optional character vector to subset genes.
#' @param output_dir Ignored (kept for backward compatibility; no files saved).
#' @param heatmap_bg Background color used when drawing the heatmap.
#' @param cluster_col Optional column in meta.data providing module labels.
#' @param keep_cluster_na Whether to keep genes with NA cluster labels.
#' @param split_heatmap_by_cluster Whether to split rows/cols by cluster labels.
#'
#' @return List with elements: matrix (Jaccard), genes, clusters, grid,
#'   min_bins, top_n, heatmap (ComplexHeatmap object), palette (color function).
#' @export
computeGeneGeneJaccardHeatmap <- function(scope_obj,
                                          grid_name = "grid30",
                                          expr_layer = c("Xz", "counts"),
                                          min_bins_per_gene = 8L,
                                          top_n_genes = 60L,
                                          genes_of_interest = NULL,
                                          output_dir = NULL,
                                          heatmap_bg = "#c0c0c0",
                                          cluster_col = NULL,
                                          keep_cluster_na = TRUE,
                                          split_heatmap_by_cluster = FALSE) {
  expr_layer <- match.arg(expr_layer)
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required to draw the Jaccard heatmap.")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required to build the heatmap palette.")
  }

  cluster_lookup <- NULL
  if (!is.null(cluster_col)) {
    meta_df <- scope_obj@meta.data
    if (is.null(meta_df) || !is.data.frame(meta_df)) {
      warning("scope_obj@meta.data is missing or not a data.frame; ignore cluster_col.")
    } else if (!(cluster_col %in% colnames(meta_df))) {
      warning("cluster column ", cluster_col, " is not present in meta.data; ignore this parameter.")
    } else if (is.null(rownames(meta_df))) {
      warning("meta.data has no rownames, cannot map gene → cluster; ignore cluster_col.")
    } else {
      cluster_lookup <- as.character(meta_df[[cluster_col]])
      names(cluster_lookup) <- rownames(meta_df)
    }
  }

  grid_layer <- scope_obj@grid[[grid_name]]
  if (is.null(grid_layer)) {
    warning("Grid layer ", grid_name, " is missing; cannot compute Jaccard.")
    return(invisible(NULL))
  }

  gene_bin_mat <- NULL
  if (identical(expr_layer, "counts")) {
    if (is.null(grid_layer$counts)) {
      warning("Grid layer ", grid_name, " has no counts; cannot compute Jaccard.")
      return(invisible(NULL))
    }
    counts_df <- as.data.frame(grid_layer$counts)
    if (!nrow(counts_df)) {
      warning("Grid layer ", grid_name, " has no gene × bin records.")
      return(invisible(NULL))
    }

    cols_expected <- c("gene", "grid_id", "count")
    if (!all(cols_expected %in% colnames(counts_df))) {
      warning("counts is missing required columns: ", paste(setdiff(cols_expected, colnames(counts_df)), collapse = ", "))
      return(invisible(NULL))
    }

    counts_df$gene <- as.character(counts_df$gene)
    counts_df$grid_id <- as.character(counts_df$grid_id)
    counts_df$count <- as.numeric(as.character(counts_df$count))
    counts_df <- counts_df[!is.na(counts_df$gene) & !is.na(counts_df$grid_id) & counts_df$count > 0, , drop = FALSE]
    if (!nrow(counts_df)) {
      warning("No gene × bin records remain after filtering NA/zero.")
      return(invisible(NULL))
    }

    counts_df <- stats::aggregate(count ~ gene + grid_id, data = counts_df, FUN = sum)
    genes <- sort(unique(counts_df$gene))
    bins <- sort(unique(counts_df$grid_id))

    gene_idx <- match(counts_df$gene, genes)
    bin_idx <- match(counts_df$grid_id, bins)
    gene_bin_mat <- Matrix::sparseMatrix(
      i = gene_idx,
      j = bin_idx,
      x = as.numeric(counts_df$count > 0),
      dims = c(length(genes), length(bins)),
      dimnames = list(genes, bins)
    )
  } else {
    if (is.null(grid_layer$Xz)) {
      warning("Grid layer ", grid_name, " has no Xz; cannot compute Jaccard.")
      return(invisible(NULL))
    }
    xz_mat <- grid_layer$Xz
    xz_mat <- as.matrix(xz_mat)
    if (!nrow(xz_mat) || !ncol(xz_mat)) {
      warning(grid_name, " Xz layer is empty; cannot compute Jaccard.")
      return(invisible(NULL))
    }
    if (is.null(colnames(xz_mat)) || is.null(rownames(xz_mat))) {
      warning("Xz layer lacks row/col names; cannot compute Jaccard.")
      return(invisible(NULL))
    }
    bin_mask <- xz_mat > 0
    bin_mask[is.na(bin_mask)] <- FALSE
    gene_bin_mat <- Matrix::Matrix(t(bin_mask) * 1, sparse = TRUE)
    rownames(gene_bin_mat) <- colnames(xz_mat)
    colnames(gene_bin_mat) <- rownames(xz_mat)
  }

  if (is.null(gene_bin_mat) || !nrow(gene_bin_mat)) {
    warning("No usable gene × bin information; cannot compute Jaccard.")
    return(invisible(NULL))
  }

  if (!is.null(genes_of_interest)) {
    genes_keep <- intersect(unique(genes_of_interest), rownames(gene_bin_mat))
    if (!length(genes_keep)) {
      warning("genes_of_interest not found in current grid.")
      return(invisible(NULL))
    }
    gene_bin_mat <- gene_bin_mat[genes_keep, , drop = FALSE]
  }

  gene_clusters <- NULL
  if (!is.null(cluster_lookup)) {
    gene_names_all <- rownames(gene_bin_mat)
    idx <- match(gene_names_all, names(cluster_lookup))
    found <- !is.na(idx)
    if (!any(found)) {
      warning("cluster column ", cluster_col, " has no overlapping genes with current grid.")
      return(invisible(NULL))
    }
    if (!keep_cluster_na) {
      found <- found & !is.na(cluster_lookup[idx])
      if (!any(found)) {
        warning("All overlapping genes have NA cluster labels; cannot continue Jaccard.")
        return(invisible(NULL))
      }
    }
    gene_bin_mat <- gene_bin_mat[found, , drop = FALSE]
    idx <- idx[found]
    gene_clusters <- cluster_lookup[idx]
    names(gene_clusters) <- rownames(gene_bin_mat)
  }

  if (!nrow(gene_bin_mat)) {
    warning("No genes remain after filtering.")
    return(invisible(NULL))
  }

  bin_presence <- Matrix::rowSums(gene_bin_mat)
  keep_idx <- which(bin_presence >= min_bins_per_gene)
  if (!length(keep_idx)) {
    warning("No genes meet the minimum bin count of ", min_bins_per_gene, ".")
    return(invisible(NULL))
  }

  gene_bin_mat <- gene_bin_mat[keep_idx, , drop = FALSE]
  bin_presence <- bin_presence[keep_idx]
  if (!is.null(gene_clusters)) {
    gene_clusters <- gene_clusters[rownames(gene_bin_mat)]
  }

  if (!is.null(top_n_genes) && is.finite(top_n_genes) && top_n_genes > 0 &&
      nrow(gene_bin_mat) > top_n_genes) {
    ord <- order(bin_presence, decreasing = TRUE)
    keep_ord <- ord[seq_len(top_n_genes)]
    gene_bin_mat <- gene_bin_mat[keep_ord, , drop = FALSE]
    bin_presence <- bin_presence[keep_ord]
    if (!is.null(gene_clusters)) {
      gene_clusters <- gene_clusters[keep_ord]
    }
  }

  if (!is.null(gene_clusters) && split_heatmap_by_cluster) {
    cluster_order <- order(gene_clusters, rownames(gene_bin_mat))
    gene_bin_mat <- gene_bin_mat[cluster_order, , drop = FALSE]
    bin_presence <- bin_presence[cluster_order]
    gene_clusters <- gene_clusters[cluster_order]
  }

  gene_names <- rownames(gene_bin_mat)
  inter_mat <- Matrix::tcrossprod(gene_bin_mat)
  inter_dense <- as.matrix(inter_mat)
  union_dense <- outer(bin_presence, bin_presence, "+") - inter_dense
  jac_mat <- inter_dense
  positive_union <- union_dense > 0
  jac_mat[positive_union] <- jac_mat[positive_union] / union_dense[positive_union]
  jac_mat[!positive_union] <- NA_real_
  diag(jac_mat) <- 1
  rownames(jac_mat) <- gene_names
  colnames(jac_mat) <- gene_names

  col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#113A70", "#f7f7f7", "#B11226"))
  row_split_factor <- column_split_factor <- NULL
  if (!is.null(gene_clusters) && split_heatmap_by_cluster) {
    cluster_labels_plot <- ifelse(is.na(gene_clusters), "NA", as.character(gene_clusters))
    cluster_levels <- unique(cluster_labels_plot)
    row_split_factor <- factor(cluster_labels_plot, levels = cluster_levels)
    column_split_factor <- row_split_factor
  }
  column_title_txt <- NULL
  ht <- ComplexHeatmap::Heatmap(
    jac_mat,
    name = "Jaccard",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 6),
    column_title = column_title_txt,
    column_title_gp = grid::gpar(fontsize = 10),
    row_split = row_split_factor,
    column_split = column_split_factor
  )

  invisible(list(
    matrix = jac_mat,
    genes = gene_names,
    clusters = gene_clusters,
    grid = grid_name,
    min_bins = min_bins_per_gene,
    top_n = top_n_genes,
    heatmap = ht,
    palette = col_fun
  ))
}

#' Compute module-level similarity from density layers
#'
#' Computes similarity between module density profiles using Jaccard, Pearson,
#' Spearman, or cosine. Optionally z-score normalizes densities and filters
#' modules with insufficient grid coverage. Returns the similarity matrix and a
#' ComplexHeatmap object; no files are written.
#'
#' @param scope_obj geneSCOPE object containing density layers.
#' @param grid_name Grid layer name (default "grid30").
#' @param density_cols Optional explicit density columns to use.
#' @param density_pattern Regex to select density columns (default "_density$").
#' @param method Similarity metric: "jaccard", "pearson", "spearman", "cosine".
#' @param density_threshold Threshold for binarization (Jaccard) or presence.
#' @param normalize Whether to z-score densities before similarity (default TRUE).
#' @param min_grids_per_module Minimum grids required per module.
#' @param top_n_modules Optional cap on number of modules; NULL keeps all.
#' @param output_dir Ignored (kept for backward compatibility; no files saved).
#' @param heatmap_bg Background color used when drawing the heatmap.
#'
#' @return List with: matrix (similarity), modules (names), grid, method,
#'   threshold, min_grids, heatmap (ComplexHeatmap object), palette (color
#'   function), labels (clean module labels).
#' @export
computeModuleDensitySimilarity <- function(scope_obj,
                                           grid_name = "grid30",
                                           density_cols = NULL,
                                           density_pattern = "_density$",
                                           method = c("jaccard", "pearson", "spearman", "cosine"),
                                           density_threshold = 0,
                                           normalize = TRUE,
                                           min_grids_per_module = 5L,
                                           top_n_modules = NULL,
                                           output_dir = NULL,
                                           heatmap_bg = "#c0c0c0") {
  method <- match.arg(method)
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required to draw the module similarity heatmap.")
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required to build the heatmap palette.")
  }
  dens_tbl <- scope_obj@density[[grid_name]]
  if (is.null(dens_tbl) || !nrow(dens_tbl)) {
    warning("density layer is missing or empty: ", grid_name)
    return(invisible(NULL))
  }

  if (is.null(rownames(dens_tbl))) {
    grid_info <- scope_obj@grid[[grid_name]]$grid_info
    if (!is.null(grid_info) && "grid_id" %in% colnames(grid_info)) {
      rownames(dens_tbl) <- as.character(grid_info$grid_id)
    }
  }

  if (is.null(density_cols)) {
    density_cols <- colnames(dens_tbl)
    if (!is.null(density_pattern)) {
      density_cols <- grep(density_pattern, density_cols, value = TRUE)
    }
  }

  density_cols <- intersect(density_cols, colnames(dens_tbl))
  if (!length(density_cols)) {
    warning("No matching density columns found.")
    return(invisible(NULL))
  }

  num_mask <- vapply(dens_tbl[density_cols], is.numeric, logical(1))
  if (!all(num_mask)) {
    warning("Non-numeric density columns ignored: ", paste(density_cols[!num_mask], collapse = ", "))
  }
  density_cols <- density_cols[num_mask]
  if (!length(density_cols)) {
    warning("No numeric density columns available for similarity computation.")
    return(invisible(NULL))
  }

  dens_mat <- as.matrix(dens_tbl[, density_cols, drop = FALSE])
  dens_mat[is.na(dens_mat)] <- 0
  if (normalize) {
    dens_mat <- scale(dens_mat, center = TRUE, scale = TRUE)
    dens_mat[is.na(dens_mat)] <- 0
  }

  module_names <- colnames(dens_mat)
  module_labels <- sub("^.*?_([0-9]+)_density$", "\\1", module_names)
  invalid_label <- module_labels == module_names | module_labels == ""
  module_labels[invalid_label] <- module_names[invalid_label]

  if (method == "jaccard") {
    bin_mat <- dens_mat > density_threshold
    bin_mat[is.na(bin_mat)] <- FALSE
    module_grid_mat <- Matrix::Matrix(t(bin_mat) * 1, sparse = TRUE)
    module_presence <- Matrix::rowSums(module_grid_mat)
  } else {
    module_grid_mat <- as.matrix(t(dens_mat))
    module_presence <- rowSums(module_grid_mat > density_threshold)
  }

  keep_idx <- which(module_presence >= min_grids_per_module)
  if (!length(keep_idx)) {
    warning("No modules meet the minimum grid count of ", min_grids_per_module, ".")
    return(invisible(NULL))
  }

  module_grid_mat <- module_grid_mat[keep_idx, , drop = FALSE]
  module_presence <- module_presence[keep_idx]
  module_names <- module_names[keep_idx]

  if (!is.null(top_n_modules) && is.finite(top_n_modules) && top_n_modules > 0 &&
      length(module_names) > top_n_modules) {
    ord <- order(module_presence, decreasing = TRUE)
    keep_ord <- ord[seq_len(top_n_modules)]
    module_grid_mat <- module_grid_mat[keep_ord, , drop = FALSE]
    module_presence <- module_presence[keep_ord]
    module_names <- module_names[keep_ord]
    module_labels <- module_labels[keep_ord]
  }

  sim_mat <- switch(method,
    jaccard = {
      inter_mat <- Matrix::tcrossprod(module_grid_mat)
      inter_dense <- as.matrix(inter_mat)
      union_dense <- outer(module_presence, module_presence, "+") - inter_dense
      positive_union <- union_dense > 0
      jac <- inter_dense
      jac[positive_union] <- jac[positive_union] / union_dense[positive_union]
      jac[!positive_union] <- NA_real_
      diag(jac) <- 1
      jac
    },
    pearson = stats::cor(t(module_grid_mat), method = "pearson"),
    spearman = stats::cor(t(module_grid_mat), method = "spearman"),
    cosine = {
      dense_mat <- module_grid_mat
      num <- dense_mat %*% t(dense_mat)
      norms <- sqrt(diag(num))
      denom <- norms %o% norms
      denom[denom == 0] <- NA_real_
      cos_mat <- num / denom
      diag(cos_mat) <- 1
      cos_mat
    }
  )

  rownames(sim_mat) <- module_names
  colnames(sim_mat) <- module_names

  if (method == "jaccard" || method == "cosine") {
    col_fun <- circlize::colorRamp2(c(0, 0.5, 1), c("#113A70", "#f7f7f7", "#B11226"))
  } else {
    col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#113A70", "#f7f7f7", "#B11226"))
  }

  heatmap_name <- paste0("ModuleSimilarity_", method)
  ht <- ComplexHeatmap::Heatmap(
    sim_mat,
    name = heatmap_name,
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_labels = module_labels,
    column_labels = module_labels,
    column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    column_names_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    heatmap_legend_param = list(title = ifelse(method == "jaccard", "Jaccard", method))
  )

  invisible(list(
    matrix = sim_mat,
    modules = module_names,
    grid = grid_name,
    method = method,
    threshold = density_threshold,
    min_grids = min_grids_per_module,
    heatmap = ht,
    palette = col_fun,
    labels = module_labels
  ))
}

#' Per-gene within vs between module Jaccard test
#'
#' For each gene, computes mean Jaccard to same-module genes and to other
#' modules, then tests the within–between difference with a paired Wilcoxon
#' signed-rank test or sign-flip permutation test.
#'
#' @param jac_mat Gene–gene Jaccard matrix (rownames/colnames = genes).
#' @param gene_modules Named vector of module labels for genes.
#' @param alternative Alternative hypothesis for the test.
#' @param method "wilcoxon" or "permutation".
#' @param n_perm Number of permutations for permutation test.
#' @param seed Optional RNG seed.
#' @param min_within Minimum within-module partners per gene.
#' @param min_between Minimum between-module partners per gene.
#'
#' @return List with gene_summary (per-gene means/counts/diff), test (htest or
#'   permutation result), alternative, method.
#' @export
computeGeneModuleWithinBetweenTest <- function(jac_mat,
                                               gene_modules,
                                               alternative = c("greater", "two.sided", "less"),
                                               method = c("wilcoxon", "permutation"),
                                               n_perm = 10000L,
                                               seed = NULL,
                                               min_within = 1L,
                                               min_between = 1L) {
  alternative <- match.arg(alternative)
  method <- match.arg(method)

  if (is.null(jac_mat) || !is.matrix(jac_mat)) {
    stop("jac_mat must be a numeric matrix with gene names as dimnames.")
  }
  if (is.null(rownames(jac_mat)) || is.null(colnames(jac_mat))) {
    stop("jac_mat must have rownames/colnames for genes.")
  }
  if (is.null(gene_modules)) {
    stop("gene_modules (named vector of module labels) is required.")
  }
  if (is.null(names(gene_modules))) {
    stop("gene_modules must be named by gene.")
  }

  gene_names <- intersect(rownames(jac_mat), names(gene_modules))
  if (!length(gene_names)) {
    stop("No overlap between jac_mat genes and gene_modules names.")
  }

  jac_mat <- jac_mat[gene_names, gene_names, drop = FALSE]
  gene_modules <- gene_modules[gene_names]

  # Drop genes with NA module labels
  keep_mod <- !is.na(gene_modules)
  jac_mat <- jac_mat[keep_mod, keep_mod, drop = FALSE]
  gene_modules <- gene_modules[keep_mod]
  gene_names <- names(gene_modules)

  if (!length(gene_names)) {
    stop("All genes have NA module labels after filtering.")
  }

  per_gene <- lapply(seq_along(gene_names), function(i) {
    g <- gene_names[[i]]
    mod <- gene_modules[[i]]
    row_vals <- jac_mat[g, ]

    within_mask <- gene_modules == mod & names(gene_modules) != g
    between_mask <- gene_modules != mod

    within_vals <- row_vals[within_mask]
    between_vals <- row_vals[between_mask]

    within_vals <- within_vals[!is.na(within_vals)]
    between_vals <- between_vals[!is.na(between_vals)]

    if (length(within_vals) < min_within || length(between_vals) < min_between) {
      return(NULL)
    }

    data.frame(
      gene = g,
      module = mod,
      within_mean = mean(within_vals),
      between_mean = mean(between_vals),
      within_n = length(within_vals),
      between_n = length(between_vals),
      diff = mean(within_vals) - mean(between_vals),
      stringsAsFactors = FALSE
    )
  })

  per_gene_df <- do.call(rbind, per_gene)
  if (is.null(per_gene_df) || !nrow(per_gene_df)) {
    warning("No genes retained for within/between summary (check min_within/min_between).")
    return(invisible(NULL))
  }

  delta <- per_gene_df$diff
  delta <- delta[is.finite(delta)]

  if (!length(delta)) {
    warning("No finite within-between differences to test.")
    return(invisible(list(gene_summary = per_gene_df, test = NULL)))
  }

  test_res <- switch(method,
    wilcoxon = {
      stats::wilcox.test(delta, mu = 0, alternative = alternative, paired = FALSE, exact = FALSE)
    },
    permutation = {
      if (!is.null(seed)) {
        set.seed(seed)
      }
      if (!is.finite(n_perm) || n_perm < 1) {
        stop("n_perm must be a positive integer for permutation test.")
      }
      obs <- mean(delta)
      perm_stats <- replicate(n_perm, {
        mean(delta * sample(c(-1, 1), length(delta), replace = TRUE))
      })
      p_val <- switch(alternative,
        greater = (sum(perm_stats >= obs) + 1) / (n_perm + 1),
        less = (sum(perm_stats <= obs) + 1) / (n_perm + 1),
        two.sided = (sum(abs(perm_stats) >= abs(obs)) + 1) / (n_perm + 1)
      )
      list(
        statistic = obs,
        p.value = p_val,
        method = "paired sign-flip permutation on per-gene within-minus-between means",
        alternative = alternative,
        n_perm = n_perm
      )
    }
  )

  invisible(list(
    gene_summary = per_gene_df,
    test = test_res,
    alternative = alternative,
    method = method
  ))
}

#' Run per-gene within vs between Jaccard workflow on a SCOPE object
#'
#' Wrapper that computes the gene–gene Jaccard matrix on a given grid/layer and
#' then runs the per-gene within vs between test. No files are written.
#'
#' @param scope_obj geneSCOPE object.
#' @param dataset_name Label used in messages.
#' @param grid_name Grid layer name (default "grid30").
#' @param cluster_col Column in meta.data for module labels.
#' @param expr_layer Expression layer to binarize ("Xz" or "counts").
#' @param output_dir Ignored (backward compatibility; no files saved).
#' @param alternative Alternative hypothesis for the test.
#' @param method "wilcoxon" or "permutation".
#' @param n_perm Number of permutations when method = "permutation".
#' @param seed RNG seed.
#'
#' @return List with `jaccard` (matrix, heatmap, clusters) and `within_between`
#'   (test results).
#' @export
runWithinBetweenTest <- function(scope_obj,
                                 dataset_name,
                                 grid_name = "grid30",
                                 cluster_col = "q95_res0.1_grid30_log1p_freq0.95",
                                 expr_layer = "Xz",
                                 output_dir = NULL,
                                 alternative = "greater",
                                 method = "wilcoxon",
                                 n_perm = 10000L,
                                 seed = 123) {
  jac_res <- computeGeneGeneJaccardHeatmap(
    scope_obj = scope_obj,
    grid_name = grid_name,
    expr_layer = expr_layer,
    cluster_col = cluster_col,
    keep_cluster_na = FALSE,   # drop genes without module labels
    top_n_genes = NULL,        # use all labeled genes, no top-N truncation
    split_heatmap_by_cluster = TRUE,
    output_dir = output_dir
  )

  if (is.null(jac_res) || is.null(jac_res$matrix) || is.null(jac_res$clusters)) {
    warning("Skipping ", dataset_name, ": missing Jaccard matrix or clusters.")
    return(invisible(NULL))
  }

  test_res <- computeGeneModuleWithinBetweenTest(
    jac_mat = jac_res$matrix,
    gene_modules = jac_res$clusters,
    alternative = alternative,
    method = method,
    n_perm = n_perm,
    seed = seed
  )

  if (is.null(test_res)) {
    warning("Test failed for ", dataset_name)
    return(invisible(NULL))
  }

  message("✓ ", dataset_name, " within>between test (", method, ", ", alternative, ") p=",
          if (!is.null(test_res$test$p.value)) signif(test_res$test$p.value, 3) else NA)

  invisible(list(
    jaccard = jac_res,
    within_between = test_res
  ))
}
