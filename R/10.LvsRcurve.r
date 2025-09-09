#' @title Simplified and robust Lee's L vs Pearson r bootstrap curve (ensure CI output)
#' @description Fit LOESS from (Lee's L, r) upper triangle point set and generate 95% confidence bands through residual bootstrap.
#' @param coordObj CoordObj
#' @param grid_name Grid layer name
#' @param level "grid" or "cell"
#' @param lee_stats_layer Lee statistics layer
#' @param span LOESS span
#' @param B Bootstrap iterations
#' @param deg 0/1
#' @param ncore Threads
#' @param length_out Number of grid points
#' @param downsample Downsample (ratio or integer)
#' @param n_strata Number of strata (too large will be truncated to unique r values-1)
#' @param k_max LOESS maximum neighbors
#' @param jitter_eps Jitter for r
#' @param ci_method "percentile","basic","bc"
#' @param ci_adjust "none"|"analytic" (only for percentile)
#' @param min_rel_width floor minimum relative width (0 to disable)
#' @param widen_span floor width smoothing span
#' @param curve_name Output name
#' @return Modified CoordObj
#' @export
addLRcurve <- function(coordObj,
                       grid_name,
                       level = c("grid", "cell"),
                       lee_stats_layer = "LeeStats_Xz",
                       span = 0.45,
                       B = 1000,
                       deg = 1,
                       ncore = max(1, parallel::detectCores() - 1),
                       length_out = 1000,
                       downsample = 1,
                       n_strata = 50,
                       k_max = Inf,
                       jitter_eps = 0,
                       ci_method = c("percentile", "basic", "bc"),
                       ci_adjust = c("none", "analytic"),
                       min_rel_width = 0,
                       widen_span = 0.1,
                       curve_name = "LR_curve2") {
  ci_method <- match.arg(ci_method)
  ci_adjust <- match.arg(ci_adjust)
  level <- match.arg(level)
  if (min_rel_width < 0) stop("min_rel_width cannot be negative")
  if (B < 20) message("[geneSCOPE::addLRcurve] !!! Warning: B < 20 may be unstable !!!")

  # 1. Extract matrices
  g_layer <- .selectGridLayer(coordObj, grid_name)
  grid_name <- names(coordObj@grid)[vapply(coordObj@grid, identical, logical(1), g_layer)]
  Lmat <- .getLeeMatrix(coordObj, grid_name, lee_layer = lee_stats_layer)
  rmat <- .getPearsonMatrix(coordObj, grid_name, level = ifelse(level == "grid", "grid", "cell"))

  common <- intersect(rownames(Lmat), rownames(rmat))
  if (length(common) < 2) stop("Insufficient common genes")
  Lmat <- Lmat[common, common, drop = FALSE]
  rmat <- rmat[common, common, drop = FALSE]
  ut <- upper.tri(Lmat, diag = FALSE)
  Lv <- Lmat[ut]
  rv <- rmat[ut]

  # 2. Clean / Downsample
  ok <- is.finite(Lv) & is.finite(rv)
  Lv <- Lv[ok]
  rv <- rv[ok]
  if (!length(Lv)) stop("No valid points")
  if (is.numeric(downsample) && downsample < 1) {
    keep <- sample.int(length(Lv), max(1L, floor(downsample * length(Lv))))
    Lv <- Lv[keep]
    rv <- rv[keep]
  } else if (is.numeric(downsample) && downsample >= 1 && length(Lv) > downsample) {
    keep <- sample.int(length(Lv), downsample)
    Lv <- Lv[keep]
    rv <- rv[keep]
  }
  if (jitter_eps > 0) rv <- jitter(rv, factor = jitter_eps)

  # 3. Stratify (by r quantiles) — fallback if failed
  uniq_r <- sort(unique(rv))
  if (length(uniq_r) < 3) stop("r values too discrete")
  n_strata_eff <- min(n_strata, length(uniq_r) - 1)
  probs <- seq(0, 1, length.out = n_strata_eff + 1)
  brks <- unique(quantile(rv, probs, na.rm = TRUE))
  if (length(brks) < 2) {
    # Fallback: uniform slices
    brks <- seq(min(rv), max(rv), length.out = n_strata_eff + 1)
    brks <- unique(brks)
  }
  if (length(brks) < 2) stop("Cannot establish strata")
  strat <- cut(rv, breaks = brks, include.lowest = TRUE, labels = FALSE)
  ok2 <- !is.na(strat)
  rv <- rv[ok2]
  Lv <- Lv[ok2]
  strat <- strat[ok2]

  # 4. xgrid
  xr <- range(rv)
  if (diff(xr) <= 0) stop("r has no span")
  xgrid <- seq(xr[1], xr[2], length.out = length_out)

  # Directly call unified interface (no longer check old version)
  keep_boot <- TRUE
  adjust_mode <- if (ci_method == "percentile" && ci_adjust == "analytic") 1L else 0L
  res <- loess_residual_bootstrap(
    x = rv, y = Lv, strat = as.integer(strat),
    grid = xgrid,
    B = as.integer(B),
    span = span,
    deg = as.integer(deg),
    n_threads = as.integer(max(1, ncore)),
    k_max = if (is.finite(k_max)) as.integer(k_max) else -1L,
    keep_boot = keep_boot,
    adjust_mode = adjust_mode,
    ci_type = 0L,
    level = 0.95
  )

  fit <- res$fit
  lo <- res$lo
  hi <- res$hi
  if (ci_method == "basic") {
    lo <- res$lo_basic
    hi <- res$hi_basic
  } else if (ci_method == "bc") {
    lo <- res$lo_bc
    hi <- res$hi_bc
  }
  # 6. floor widen
  if (min_rel_width > 0) {
    rng_fit <- diff(range(fit, finite = TRUE))
    if (!is.finite(rng_fit) || rng_fit <= 0) rng_fit <- 1
    width <- hi - lo
    target_w <- min_rel_width * rng_fit
    width <- pmax(width, target_w)
    # Smooth width (LOESS)
    if (is.finite(widen_span) && widen_span > 0 && length(width) > 10) {
      sm <- tryCatch(loess(width ~ xgrid, span = widen_span), error = function(e) NULL)
      if (!is.null(sm)) {
        wp <- tryCatch(predict(sm, xgrid), error = function(e) NULL)
        if (!is.null(wp) && all(is.finite(wp))) {
          width <- pmax(wp, target_w)
        }
      }
    }
    center <- (lo + hi) / 2
    lo <- center - width / 2
    hi <- center + width / 2
  }

  # ---- New: Local residual scale adaptive lower bound (no new parameters) ----------------------------
  # Purpose: Avoid "points more dispersed but CI narrower"; force width >= 2*1.96*local MAD
  {
    if (length(Lv) > 20) {
      # Fit values interpolated to original r
      fit_at_rv <- tryCatch(approx(xgrid, fit, xout = rv, rule = 2, ties = "ordered")$y,
        error = function(e) rep(mean(fit), length(rv))
      )
      res_raw <- Lv - fit_at_rv
      ord_r <- order(rv)
      rv_sorted <- rv[ord_r]
      res_sorted <- res_raw[ord_r]
      n_pts <- length(rv_sorted)
      # Fixed K (internal constant, no new parameters)
      K <- min(200L, n_pts)
      if (K > 5) {
        # Nearest neighbor index function (bidirectional expansion)
        get_knn_idx <- function(x0) {
          pos <- findInterval(x0, rv_sorted)
          l <- pos
          r <- pos + 1
          out <- integer(0)
          while (length(out) < K && (l >= 1 || r <= n_pts)) {
            dl <- if (l >= 1) abs(rv_sorted[l] - x0) else Inf
            dr <- if (r <= n_pts) abs(rv_sorted[r] - x0) else Inf
            if (dl <= dr) {
              if (l >= 1) {
                out <- c(out, l)
                l <- l - 1
              }
            } else {
              if (r <= n_pts) {
                out <- c(out, r)
                r <- r + 1
              }
            }
          }
          out
        }
        # Calculate local MAD
        loc_mad <- vapply(xgrid, function(xx) {
          idx <- get_knn_idx(xx)
          if (!length(idx)) {
            return(NA_real_)
          }
          rr <- res_sorted[idx]
          med <- stats::median(rr)
          1.4826 * stats::median(abs(rr - med))
        }, numeric(1))
        # Global benchmark (prevent all NA)
        med_mad <- stats::median(loc_mad[is.finite(loc_mad) & loc_mad > 0], na.rm = TRUE)
        if (is.finite(med_mad) && med_mad > 0) {
          width <- hi - lo
          # Required minimum width
          w_need <- 2 * 1.96 * loc_mad
          # If local MAD is NA or 0, no restriction
          bad <- !is.finite(w_need) | w_need <= 0
          if (any(!bad)) {
            width_new <- width
            width_new[!bad] <- pmax(width[!bad], w_need[!bad])
            if (!identical(width_new, width)) {
              center <- (lo + hi) / 2
              lo <- center - width_new / 2
              hi <- center + width_new / 2
              # Record diagnostics (add to meta)
              local_mad_diag <- list(
                K = K,
                mad_global_med = med_mad,
                frac_expanded = mean(width_new > width)
              )
            }
          }
        }
      }
    }
  }

  # ---- New: CI edge smoothing (post-processing only, does not change algorithm) ---------------------------
  {
    if (length(lo) >= 15 && all(is.finite(xgrid))) {
      span_s <- 0.06 + 10 / length(lo) # Adaptive small span, smaller with longer length
      lo_s <- tryCatch(
        predict(loess(lo ~ xgrid,
          span = span_s, degree = 1,
          surface = "direct", family = "gaussian"
        )),
        error = function(e) lo
      )
      hi_s <- tryCatch(
        predict(loess(hi ~ xgrid,
          span = span_s, degree = 1,
          surface = "direct", family = "gaussian"
        )),
        error = function(e) hi
      )
      # Use smoothed width as main, keep original center to avoid overall drift
      center_old <- (lo + hi) / 2
      width_new <- pmax(hi_s - lo_s, 0)
      lo_new <- center_old - width_new / 2
      hi_new <- center_old + width_new / 2
      # Prevent numerical error causing reversal
      swap_idx <- which(hi_new < lo_new)
      if (length(swap_idx)) {
        tmp <- lo_new[swap_idx]
        lo_new[swap_idx] <- hi_new[swap_idx]
        hi_new[swap_idx] <- tmp
      }
      # If NA produced, fallback
      if (all(is.finite(lo_new)) && all(is.finite(hi_new))) {
        edge_smooth_info <- list(
          edge_smooth = TRUE,
          edge_smooth_span = span_s,
          edge_smooth_expanded = mean(width_new > (hi - lo))
        )
        lo <- lo_new
        hi <- hi_new
      } else {
        edge_smooth_info <- list(edge_smooth = FALSE)
      }
    } else {
      edge_smooth_info <- list(edge_smooth = FALSE)
    }
  }

  # ---- CI edge smoothing completed, enter length check ----
  ## ---- 6. Write back to @stats -------------------------------------------------------
  n_grid <- length(xgrid)
  if (!all(
    length(fit) == n_grid,
    length(lo) == n_grid,
    length(hi) == n_grid
  )) {
    stop(sprintf(
      "Internal length mismatch: fit=%d lo=%d hi=%d xgrid=%d",
      length(fit), length(lo), length(hi), n_grid
    ))
  }

  if (is.null(coordObj@stats)) coordObj@stats <- list()
  if (is.null(coordObj@stats[[grid_name]])) {
    coordObj@stats[[grid_name]] <- list()
  }
  if (is.null(coordObj@stats[[grid_name]][[lee_stats_layer]])) {
    coordObj@stats[[grid_name]][[lee_stats_layer]] <- list()
  }

  ## ---- meta update ----
  meta_obj <- list(
    B = res$B,
    span = span,
    deg = deg,
    ci_method = ci_method,
    ci_adjust = ci_adjust,
    min_rel_width = min_rel_width,
    n_strata = n_strata,
    k_max = k_max,
    edf = res$edf,
    sigma2_raw = res$sigma2_raw,
    sigma2_edf = res$sigma2_edf,
    resid_global_mad = res$resid_mad,
    adjust_mode = adjust_mode,
    note = "Generated by addLRcurve2 with residual-MAD floor + edge smoothing",
    local_mad_diag = if (exists("local_mad_diag")) local_mad_diag else NULL,
    edge_smooth = if (exists("edge_smooth_info")) edge_smooth_info else NULL
  )
  meta_col <- rep(list(meta_obj), length(fit))

  coordObj@stats[[grid_name]][[lee_stats_layer]][[curve_name]] <-
    data.frame(Pear = xgrid, fit = fit, lo95 = lo, hi95 = hi, meta = I(meta_col))

  invisible(coordObj)
}

#' @title Scatter plot of Lee's L versus Pearson correlation
#'
#' @description
#'   Builds a scatter plot comparing Lee's spatial correlation statistic
#'   (L) with the conventional Pearson \(r\) for every gene pair in a
#'   selected grid layer or at the cell level.  Optionally flips axes and
#'   annotates the largest positive and negative L–rplotLvsR differences
#'   (\code{Delta}).
#'
#' @param coordObj    A \code{CoordObj} containing both Lee's L matrices
#'                    (in \code{LeeStats_Xz}) and Pearson correlation
#'                    matrices.
#' @param grid_name   Character. Grid sub-layer name.
#' @param pear_level  Character, \code{"cell"} or \code{"grid"}; which
#'                    Pearson matrix to use.
#' @param delta_top_n Integer. How many extreme \code{Delta} pairs to label
#'                    on the plot.
#' @param flip        Logical. If \code{TRUE}, put Pearson on the x-axis and
#'                    Lee's L on the y-axis.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{getTopDeltaL}}
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline scale_x_continuous scale_y_continuous labs theme_minimal theme element_text element_rect element_line element_blank
#' @importFrom ggrepel geom_text_repel
#' @importFrom dplyr mutate bind_rows slice_max slice_min distinct
#' @importFrom scales label_number
#' @export
plotLvsR <- function(coordObj,
                     grid_name,
                     pear_level = c("cell", "grid"),
                     delta_top_n = 10,
                     flip = TRUE) {
  pear_level <- match.arg(pear_level)

  ## ---- 1. Lee's L and Pearson r ------------------------------------------------
  L_mat <- .getLeeMatrix(coordObj, grid_name, lee_layer = "LeeStats_Xz")
  r_mat <- .getPearsonMatrix(coordObj, grid_name,
    level = pear_level
  )

  common <- intersect(rownames(L_mat), rownames(r_mat))
  if (length(common) < 2) stop("Insufficient common genes for plotting.")
  L_mat <- L_mat[common, common]
  r_mat <- r_mat[common, common]
  diag(L_mat) <- NA
  diag(r_mat) <- NA

  ## ---- 2. Convert to long table + Δ ---------------------------------------------------------
  genes <- colnames(L_mat)
  ut <- upper.tri(L_mat, diag = FALSE)
  df_long <- data.frame(
    gene1 = rep(genes, each = length(genes))[ut],
    gene2 = rep(genes, length(genes))[ut],
    LeesL = L_mat[ut],
    Pear  = r_mat[ut]
  ) |>
    dplyr::mutate(Delta = LeesL - Pear)

  ## ---- 3. Label points (extreme Δ) ----------------------------------------------------
  df_label <- dplyr::bind_rows(
    dplyr::slice_max(df_long, Delta, n = delta_top_n, with_ties = FALSE),
    dplyr::slice_min(df_long, Delta, n = delta_top_n, with_ties = FALSE)
  ) |>
    dplyr::distinct(gene1, gene2, .keep_all = TRUE) |>
    dplyr::mutate(label = sprintf("%s–%s\nL=%.3f", gene1, gene2, LeesL))

  ## ---- 4. Plot ----------------------------------------------------------------
  if (!flip) {
    p <- ggplot2::ggplot(
      df_long,
      ggplot2::aes(x = LeesL, y = Pear)
    )
    xlab <- "Lee’s L"
    ylab <- "Pearson correlation"
    ttl <- sprintf(
      "Lee’s L vs Pearson  (%s, %s)",
      sub("grid_lenGrid", "Grid ", grid_name),
      ifelse(pear_level == "cell", "cell level", "grid level")
    )
  } else {
    p <- ggplot2::ggplot(
      df_long,
      ggplot2::aes(x = Pear, y = LeesL)
    )
    xlab <- "Pearson correlation"
    ylab <- "Lee’s L"
    ttl <- sprintf(
      "Pearson vs Lee’s L  (%s, %s)",
      sub("grid_lenGrid", "Grid ", grid_name),
      ifelse(pear_level == "cell", "cell level", "grid level")
    )
  }

  p <- p +
    ggplot2::geom_point(alpha = 0.3, size = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", size = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", size = 0.4) +
    ggplot2::scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
    ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.01))

  if (delta_top_n > 0 && nrow(df_label) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = df_label,
      aes(label = label),
      size = 3,
      box.padding = 0.25,
      min.segment.length = 0,
      max.overlaps = Inf
    )
  }

  p + ggplot2::labs(title = ttl, x = xlab, y = ylab) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      panel.title      = ggplot2::element_text(hjust = .5, size = 10),
      panel.border     = ggplot2::element_rect(colour = "black", fill = NA, size = .5),
      panel.background = ggplot2::element_rect(fill = "#c0c0c0", colour = NA),
      panel.grid.major = ggplot2::element_line(colour = "white"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title       = ggplot2::element_text(size = 9),
      axis.text        = ggplot2::element_text(size = 8),
      axis.ticks       = ggplot2::element_blank()
    )
}


#' @title Retrieve top gene pairs by Delta (Lee's L - Pearson r) with simplified permutation test
#' @description Only perform fixed number of permutations on selected top gene pairs; supports optional non-negative truncation of Pearson r to form different definitions of Delta statistic.
#'
#' @details
#' Compared to Seurat::FindMarkers p-values "appear more continuous" due to:
#' \itemize{
#'   \item FindMarkers typically uses parametric/asymptotic tests (e.g., wilcox/t/MAST/DESeq2/LR), directly computing p-values from analytical distributions → theoretically continuous (floating-point display can have multiple decimal places).
#'   \item This function uses a finite number of permutations (perms), p-value resolution is limited by the grid of 1/(perms+1); even with beta / mid
#'         formulas only providing smooth correction, it remains discrete; the uniform scheme adds random noise to break ties, but the expectation remains
#'         equal to the original discrete grid.
#'   \item To achieve a more continuous appearance similar to FindMarkers, you can:
#'     \enumerate{
#'       \item Increase perms (reduce step size);
#'       \item Choose pval_mode='uniform' to randomize levels;
#'       \item Or add an analytical test mode (not yet implemented), for example, establish a joint null model of Lee's L and r, then
#'             use approximate normal/asymptotic distribution to infer Delta (requires additional theoretical support).
#'     }
#'   \item Discrete p-values are not inherently wrong, often statistically more conservative.
#' }
#'
#' @param clamp_mode Character: c("none","ref_only","both") defines Delta statistic:
#'   \describe{
#'     \item{none}{Delta = L - r (original, retains positive and negative correlation information)}
#'     \item{ref_only}{Reference value uses r* = max(r,0), permutation retains original r → p-values conservative}
#'     \item{both}{Current implementation still equivalent to ref_only (reference truncation only), marked as Delta_clamp_Ronly}
#'   }
#' @param p_adj_mode Multiple correction mode (see below).
#' @param pval_mode Character: c("beta","mid","uniform"):
#'   \describe{
#'     \item{beta}{(k+1)/(N+2) Jeffreys smoothing (default)}
#'     \item{mid}{(k+0.5)/(N+1) mid-p, slightly less conservative but still discrete}
#'     \item{uniform}{(k+U)/(N+1) random jitter to alleviate ties (expectation unchanged)}
#'   }
#' @section Definition of Universe:
#'   universe = total_universe = all candidate gene pairs after filtering by pear_range and L_range (nrow(df)).
#'   The top set is just the extreme pairs selected from the universe. In *_universe mode, although only the top is permuted,
#'   multiple correction uses the universe as the total number of hypothesis tests (more conservative).
#' @return data.frame containing gene1,gene2,LeesL,Pear,Delta,raw_p,FDR,stat_type,pval_mode
#' @details
#' (Supplement) Effect of permutation count on "significance threshold":
#' Actual α threshold remains unchanged; increasing perms only reduces Monte Carlo sampling error and refines p-value step size (≈1/(N+1)).
#' Output column mc_se≈sqrt(p*(1-p)/N), and Clopper–Pearson interval
#' can be used to assess whether more perms are needed.
#' @export
getTopDeltaL <- function(coordObj,
                          grid_name,
                          pear_level = c("cell", "grid"),
                          pear_range = c(-1, 1),
                          L_range = c(-1, 1),
                          top_n = 10,
                          direction = c("largest", "smallest", "both"),
                          do_perm = TRUE,
                          perms = 1000,
                          block_side = 8,
                          use_blocks = TRUE,
                          ncores = 1,
                          clamp_mode = c("none", "ref_only", "both"),
                          p_adj_mode = c("BH", "BY", "BH_universe", "BY_universe", "bonferroni"),
                          mem_limit_GB = 2,
                          pval_mode = c("beta", "mid", "uniform")) {
  pear_level <- match.arg(pear_level)
  direction <- match.arg(direction)
  p_adj_mode <- match.arg(p_adj_mode)
  clamp_mode <- match.arg(clamp_mode)
  pval_mode <- match.arg(pval_mode)

  # ---- Core count clipping ----
  os_type <- tryCatch(detectOS(), error = function(e) "linux")
  phys <- parallel::detectCores(FALSE)
  logi <- parallel::detectCores(TRUE)
  safe_cores <- switch(os_type,
    windows = min(4, phys),
    linux   = max(1, min(phys - 2, ceiling(logi * 0.6))),
    macos   = max(1, min(phys - 1, ceiling(logi * 0.5))),
    max(1, phys - 2)
  )
  if (ncores > safe_cores) {
    message("[geneSCOPE::addLRcurve] Adjusting ncores: requested=", ncores, " safe=", safe_cores)
    ncores <- safe_cores
  }
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)
  Sys.setenv(OMP_NUM_THREADS = ncores)

  # ---- Extract matrices ----
  L_mat <- .getLeeMatrix(coordObj, grid_name, lee_layer = "LeeStats_Xz")
  r_mat <- .getPearsonMatrix(coordObj, grid_name, level = pear_level)
  common <- intersect(rownames(L_mat), rownames(r_mat))
  if (length(common) < 2) stop("Insufficient common genes")
  L_mat <- L_mat[common, common]
  r_mat <- r_mat[common, common]
  diag(L_mat) <- NA
  diag(r_mat) <- NA
  ut <- upper.tri(L_mat)
  LeesL_vec <- L_mat[ut]
  Pear_vec <- r_mat[ut]

  # ---- Delta ----
  internal_clamp_mode <- clamp_mode
  if (clamp_mode == "both") {
    message("[geneSCOPE::addLRcurve] clamp_mode='both' currently equivalent to ref_only (reference truncation only)")
    internal_clamp_mode <- "ref_only"
  }
  Pear_for_delta <- if (internal_clamp_mode == "ref_only") pmax(Pear_vec, 0) else Pear_vec

  df <- data.frame(
    gene1 = rep(common, each = length(common))[ut],
    gene2 = rep(common, length(common))[ut],
    LeesL = LeesL_vec,
    Pear  = Pear_vec,
    Delta = LeesL_vec - Pear_for_delta
  )

  # ===== Fix coverage: first use counts long table (grid_id,gene,count) then fallback to matrix =====
  {
    expr_pct_map <- setNames(rep(0, length(common)), common) # Default 0
    g_layer_try <- tryCatch(.selectGridLayer(coordObj, grid_name), error = function(e) NULL)
    coverage_done <- FALSE

    # --- 1. Prefer counts long table ---
    if (!is.null(g_layer_try) && !is.null(g_layer_try$counts)) {
      ct <- g_layer_try$counts
      if (is.data.frame(ct) && all(c("gene", "grid_id") %in% colnames(ct))) {
        total_cells <- if (!is.null(g_layer_try$grid_info)) nrow(g_layer_try$grid_info) else length(unique(ct$grid_id))
        if (total_cells <= 0) total_cells <- NA_real_
        if (inherits(ct, "data.table")) {
          # Only count common genes, uniqueN(grid_id) avoids duplicates
          gene_cells <- ct[gene %in% common, .(cells = uniqueN(grid_id)), by = gene]
        } else if (requireNamespace("data.table", quietly = TRUE)) {
          dct <- data.table::as.data.table(ct)
          gene_cells <- dct[gene %in% common, .(cells = data.table::uniqueN(grid_id)), by = gene]
        } else {
          # base: first split by gene, then deduplicate and count grid_id
          keep <- ct$gene %in% common
          if (any(keep)) {
            spl <- split(ct$grid_id[keep], ct$gene[keep])
            gene_cells <- data.frame(
              gene = names(spl),
              cells = vapply(spl, function(v) length(unique(v)), integer(1)),
              row.names = NULL
            )
          } else {
            gene_cells <- data.frame(gene = character(0), cells = integer(0))
          }
        }
        if (nrow(gene_cells)) {
          expr_pct_map[gene_cells$gene] <- (gene_cells$cells / total_cells) * 100
        }
        coverage_done <- TRUE
      }
    }

    # --- 2. Fallback to matrix (>0 row coverage) ---
    if (!coverage_done) {
      pick_matrix <- function(layer) {
        if (is.null(layer)) {
          return(NULL)
        }
        for (nm in c("counts", "raw_counts", "expr", "X", "data", "logCPM", "Xz")) {
          m <- layer[[nm]]
          # Skip long table data.frame: needs to have dim/be (sparse) matrix
          if (!is.null(m) && (is.matrix(m) || inherits(m, "dgCMatrix"))) {
            return(m)
          }
        }
        NULL
      }
      expr_mat <- pick_matrix(g_layer_try)
      if (is.null(expr_mat) && pear_level == "cell") {
        cell_env <- tryCatch(coordObj@cell, error = function(e) NULL)
        if (is.null(cell_env)) cell_env <- tryCatch(coordObj@cells, error = function(e) NULL)
        expr_mat <- pick_matrix(cell_env)
      }
      if (!is.null(expr_mat)) {
        rn <- rownames(expr_mat)
        cn <- colnames(expr_mat)
        inter_col <- intersect(cn, common)
        inter_row <- intersect(rn, common)
        if (length(inter_col) == 0 && length(inter_row) == 0) {
          # Keep all 0
        } else {
          gene_in_col <- (length(inter_col) >= length(inter_row))
          if (!gene_in_col) {
            expr_mat <- if (inherits(expr_mat, "dgCMatrix")) Matrix::t(expr_mat) else t(expr_mat)
            cn <- colnames(expr_mat)
            inter_col <- intersect(cn, common)
          }
          if (length(inter_col)) {
            nz_counts <- if (inherits(expr_mat, "dgCMatrix")) {
              Matrix::colSums(expr_mat[, inter_col, drop = FALSE] > 0)
            } else {
              colSums(expr_mat[, inter_col, drop = FALSE] > 0)
            }
            expr_pct_map[inter_col] <- (nz_counts / nrow(expr_mat)) * 100
          }
        }
      }
    }

    df$gene1_expr_pct <- round(expr_pct_map[df$gene1], 3)
    df$gene2_expr_pct <- round(expr_pct_map[df$gene2], 3)
  }
  # ===== Coverage calculation ends =====

  # ---- Threshold filtering ----
  df <- df[df$Pear >= pear_range[1] & df$Pear <= pear_range[2] &
    df$LeesL >= L_range[1] & df$LeesL <= L_range[2], ]
  if (!nrow(df)) stop("Thresholds remove all pairs")
  total_universe <- nrow(df)

  # ---- Select top ----
  sel <- switch(direction,
    largest = dplyr::slice_max(df, Delta, n = top_n),
    smallest = dplyr::slice_min(df, Delta, n = top_n),
    both = dplyr::bind_rows(
      dplyr::slice_max(df, Delta, n = top_n),
      dplyr::slice_min(df, Delta, n = top_n)
    )
  )
  sel <- dplyr::distinct(sel, gene1, gene2, .keep_all = TRUE)
  rownames(sel) <- NULL
  if (!nrow(sel) || !do_perm) {
    sel$stat_type <- switch(clamp_mode,
      none     = "Delta",
      ref_only = "Delta_refClamp",
      both     = "Delta_clamp_Ronly"
    )
    return(sel)
  }

  # ---- Permutation preparation ----
  g_layer <- .selectGridLayer(coordObj, grid_name)
  Xz <- g_layer$Xz
  W <- g_layer$W
  grid_info <- g_layer$grid_info
  if (is.null(Xz) || is.null(W)) stop("Grid layer missing Xz/W")
  genes_top <- unique(c(sel$gene1, sel$gene2))
  gene_map <- match(genes_top, colnames(Xz))
  if (any(is.na(gene_map))) stop("Selected genes not found in Xz")
  gene_pairs <- cbind(
    match(sel$gene1, genes_top) - 1L,
    match(sel$gene2, genes_top) - 1L
  )
  delta_ref <- sel$Delta

  block_id <- if (use_blocks) {
    bx <- (grid_info$gx - 1L) %/% block_side
    by <- (grid_info$gy - 1L) %/% block_side
    max_by <- max(by)
    bx * (max_by + 1L) + by + 1L
  } else {
    seq_len(nrow(Xz))
  }
  split_rows <- split(seq_along(block_id), block_id)

  # ---- CSR weights ----
  W_rows <- lapply(seq_len(nrow(W)), function(i) {
    nz <- which(W[i, ] != 0)
    if (length(nz)) list(indices = nz - 1L, values = as.numeric(W[i, nz])) else list(indices = integer(0), values = numeric(0))
  })
  W_row_lengths <- vapply(W_rows, function(x) length(x$indices), integer(1))
  W_row_ptr <- c(0L, cumsum(W_row_lengths))
  W_indices <- unlist(lapply(W_rows, `[[`, "indices"), use.names = FALSE)
  W_values <- unlist(lapply(W_rows, `[[`, "values"), use.names = FALSE)
  if (!length(W_indices)) stop("Weight matrix has no non-zero entries")
  Xz_sub <- Xz[, gene_map, drop = FALSE]
  n_cells <- nrow(Xz_sub)

  # ---- batch size ----
  target_batch <- min(100L, perms)
  max_idx_bytes <- mem_limit_GB * 1024^3 * 0.30
  est_bytes <- function(bs) n_cells * bs * 4
  while (target_batch > 1L && est_bytes(target_batch) > max_idx_bytes) {
    target_batch <- max(1L, floor(target_batch / 2))
  }
  message(
    "[geneSCOPE::addLRcurve] Planned batch_size=", target_batch,
    " (estimated idx_mat ", sprintf("%.2f", est_bytes(target_batch) / 1024^2),
    " MB / limit ", sprintf("%.2f", max_idx_bytes / 1024^2), " MB)"
  )

  # ---- Permutation loop ----
  remaining <- perms
  exceed_count <- rep(0L, nrow(sel))
  perm_threads <- ncores
  attempt <- 1
  while (remaining > 0) {
    bsz <- min(target_batch, remaining)
    success <- FALSE
    while (!success && perm_threads >= 1) {
      message(sprintf(
        "[geneSCOPE::addLRcurve] Permutation attempt #%d threads=%d batch=%d remain=%d clamp_mode=%s",
        attempt, perm_threads, bsz, remaining, clamp_mode
      ))
      attempt <- attempt + 1
      idx_mat <- matrix(integer(0), nrow = n_cells, ncol = bsz)
      for (p in seq_len(bsz)) {
        new_order <- sample(length(split_rows))
        idx_mat[, p] <- unlist(split_rows[new_order], use.names = FALSE) - 1L
      }
      res <- tryCatch(
        {
          if (use_blocks) {
            delta_lr_perm_csr_block(
              Xz_sub, W_indices, W_values, W_row_ptr, idx_mat,
              as.integer(block_id) - 1L, gene_pairs, delta_ref,
              perm_threads
            )
          } else {
            delta_lr_perm_csr(
              Xz_sub, W_indices, W_values, W_row_ptr, idx_mat,
              gene_pairs, delta_ref,
              perm_threads
            )
          }
        },
        error = function(e) e
      )
      if (inherits(res, "error")) {
        message("[geneSCOPE::addLRcurve]  Batch failed: ", conditionMessage(res))
        if (perm_threads > 1) {
          perm_threads <- max(1, floor(perm_threads / 2))
          next
        } else {
          stop("Permutation failed at single-thread: ", conditionMessage(res))
        }
      } else {
        exceed_count <- exceed_count + res
        success <- TRUE
      }
    }
    remaining <- remaining - bsz
  }

  # ---- p-values and multiple correction ----
  N <- perms
  k <- exceed_count
  p_values <- switch(pval_mode,
    beta    = (k + 1) / (N + 2),
    mid     = (k + 0.5) / (N + 1),
    uniform = (k + runif(length(k))) / (N + 1)
  )
  p_values[p_values > 1] <- 1
  mc_se <- sqrt(p_values * (1 - p_values) / N)
  p_ci_lo <- qbeta(0.025, k + 1, N - k + 1)
  p_ci_hi <- qbeta(0.975, k + 1, N - k + 1)

  FDR <- switch(p_adj_mode,
    BH          = stats::p.adjust(p_values, "BH"),
    BY          = stats::p.adjust(p_values, "BY"),
    BH_universe = stats::p.adjust(p_values, "BH", n = total_universe),
    BY_universe = stats::p.adjust(p_values, "BY", n = total_universe),
    bonferroni  = stats::p.adjust(p_values, "bonferroni", n = total_universe)
  )

  sel$raw_p <- p_values
  sel$mc_se <- mc_se
  sel$p_ci_lo <- p_ci_lo
  sel$p_ci_hi <- p_ci_hi
  sel$FDR <- FDR
  sel$stat_type <- switch(clamp_mode,
    none     = "Delta",
    ref_only = "Delta_refClamp",
    both     = "Delta_clamp_Ronly"
  )
  sel$pval_mode <- pval_mode
  if (p_adj_mode == "bonferroni") {
    sel <- sel[sel$FDR < 0.05, , drop = FALSE]
    if (!nrow(sel)) message("[geneSCOPE::addLRcurve] No Bonferroni-significant pairs (FDR < 0.05).")
  }

  sel
}
