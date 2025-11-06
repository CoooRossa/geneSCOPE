#' Compute spatial weights for a grid layer
#'
#' Builds a Queen/hex/rook neighbour list for the active grid cells in a
#' \code{scope_object} layer, converts it to the requested output formats, and
#' stores the results in-place. The implementation shares the same behaviour as
#' the legacy `computeWeights()` helper but delegates major steps to modular
#' helpers so they can be reused by higher level tooling (e.g. fuzzy-queen
#' diagnostics).
#'
#' @param scope_obj A \code{scope_object} with populated \code{@grid} slots.
#' @param grid_name Optional character string selecting the grid sub-layer.
#'   When \code{NULL} and only one layer exists the layer is auto-selected.
#' @param style Character weighting style passed to \code{spdep::nb2listw}.
#'   Default `"B"` (binary).
#' @param store_mat Logical; store the binary adjacency matrix in the grid's
#'   \code{$W} slot (default \code{TRUE}).
#' @param zero_policy Logical forwarded to \code{spdep::nb2listw}; when
#'   \code{TRUE} (default) zero-neighbour regions are allowed.
#' @param store_listw Logical; store the \code{listw} object (default
#'   \code{TRUE}).
#' @param verbose Logical toggle for progress messages.
#' @param topology Character flag controlling the neighbourhood lattice.
#'   `"auto"` inspects the grid geometry and platform metadata to decide
#'   between queen/rook/hex variants。`"fuzzy"` 在生成 Queen 邻接的基础上，
#'   额外存储一份行归一化后的矩阵以支持模糊扩散分析。
#'
#' @return Invisibly returns the modified \code{scope_object}.
#'
#' @importFrom spdep nb2listw
#' @importFrom Matrix Diagonal rowSums
#' @importFrom stats median
#' @export
computeWeights <- function(scope_obj,
                           grid_name = NULL,
                           style = "B",
                           store_mat = TRUE,
                           zero_policy = TRUE,
                           store_listw = TRUE,
                           verbose = TRUE,
                           topology = c("auto", "queen", "rook", "hex", "fuzzy")) {
  topology <- match.arg(topology)
  use_fuzzy <- identical(topology, "fuzzy")
  topo_arg <- if (use_fuzzy) "queen" else topology
  if (verbose) message("[geneSCOPE::computeWeights] Computing spatial weights for grid layer")

  guard <- .spw_thread_guard()
  on.exit(guard$restore(), add = TRUE)

  # 1) Select grid layer and basic metadata
  layer_info <- .spw_select_layer(scope_obj, grid_name)
  grid_name <- layer_info$grid_name
  layer <- layer_info$layer
  meta <- .spw_validate_grid_metadata(layer, grid_name)

  # 2) Determine lattice topology (auto or forced)
  topo_info <- .spw_choose_topology(scope_obj, layer, topo_arg, verbose)

  # 3) Build the full neighbourhood structure and relabel to active cells
  nb_bundle <- .spw_build_neighbourhood(meta, topo_info, zero_policy)
  relabel_nb <- nb_bundle$nb
  zero_idx <- nb_bundle$zero_idx

  # 4) Convert to requested outputs
  listw_obj <- NULL
  if (store_listw || store_mat) {
    listw_obj <- spdep::nb2listw(relabel_nb,
      style = style,
      zero.policy = zero_policy
    )
  }

  if (store_mat) {
    W <- geneSCOPE:::listw_B_omp(relabel_nb)
    diag(W) <- 0
    W@x[] <- 1
    dimnames(W) <- list(meta$grid_id, meta$grid_id)
    scope_obj@grid[[grid_name]]$W <- W
  }

  if (store_listw) {
    scope_obj@grid[[grid_name]]$listw <- listw_obj
  }

  # Optionally restore empty neighbour vectors for downstream consumers
  if (zero_policy && length(zero_idx)) {
    for (j in zero_idx) relabel_nb[[j]] <- integer(0)
  }

  if (use_fuzzy) {
    if (!store_mat) {
      stop("Fuzzy queen topology requires store_mat = TRUE to retain the adjacency matrix.")
    }
    W_base <- scope_obj@grid[[grid_name]]$W
    scope_obj@grid[[grid_name]]$fuzzy <- list(
      W_binary = W_base,
      W_row_normalised = .row_normalize_sparse(W_base)
    )
    if (verbose) {
      message("[geneSCOPE::computeWeights] Stored queen weights in grid[['W']] and fuzzy weights under grid[['fuzzy']].")
    }
  }

  if (verbose) message("[geneSCOPE::computeWeights] Spatial weights computation completed")
  invisible(scope_obj)
}

# ---------------------------------------------------------------------------
# Internal helpers for spatial weights
# ---------------------------------------------------------------------------

.spw_thread_guard <- function() {
  thread_config <- configureThreadsFor("openmp_only",
    ncores_requested = parallel::detectCores(),
    restore_after = TRUE
  )
  restore_threads <- attr(thread_config, "restore_function")

  blas_restore <- NULL
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_threads <- tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) NA_integer_)
    RhpcBLASctl::blas_set_num_threads(1L)
    blas_restore <- function() {
      if (!is.na(old_threads)) {
        RhpcBLASctl::blas_set_num_threads(old_threads)
      }
    }
  }

  list(
    restore = function() {
      if (!is.null(blas_restore)) blas_restore()
      if (!is.null(restore_threads)) restore_threads()
    }
  )
}

.spw_select_layer <- function(scope_obj, grid_name, verbose = FALSE) {
  layer <- .selectGridLayer(scope_obj, grid_name)
  if (is.null(grid_name)) {
    matches <- vapply(scope_obj@grid, identical, logical(1), layer)
    grid_name <- names(scope_obj@grid)[which(matches)][1]
  }
  if (is.null(grid_name) || !nzchar(grid_name)) {
    stop("Unable to resolve grid_name for spatial weights.")
  }
  list(layer = layer, grid_name = grid_name)
}

.spw_validate_grid_metadata <- function(layer, grid_name) {
  grid_info <- layer$grid_info
  required_cols <- c("gx", "gy", "grid_id")
  if (!all(required_cols %in% names(grid_info))) {
    stop("grid_info for layer ", grid_name, " must contain columns: gx, gy, grid_id.")
  }
  if (is.null(layer$xbins_eff) || is.null(layer$ybins_eff)) {
    stop("grid layer missing xbins_eff / ybins_eff metadata")
  }
  gx <- grid_info$gx
  gy <- grid_info$gy
  xbins_eff <- layer$xbins_eff
  ybins_eff <- layer$ybins_eff
  if (any(gx < 1 | gx > xbins_eff | gy < 1 | gy > ybins_eff)) {
    stop("gx/gy exceed xbins_eff × ybins_eff range; check coordinates.")
  }
  list(
    grid_info = data.table::as.data.table(grid_info),
    grid_id = grid_info$grid_id,
    gx = gx,
    gy = gy,
    xbins_eff = xbins_eff,
    ybins_eff = ybins_eff
  )
}

.spw_choose_topology <- function(scope_obj, layer, topology, verbose) {
  forced_topology <- !identical(topology, "auto")
  grid_info <- data.table::as.data.table(layer$grid_info)
  has_centers <- all(c("gx", "gy", "center_x", "center_y") %in% names(grid_info))

  platform_vals <- c(
    if (!is.null(scope_obj@meta.data) && nrow(scope_obj@meta.data) && "platform" %in% names(scope_obj@meta.data))
      as.character(scope_obj@meta.data$platform) else NULL,
    if (!is.null(scope_obj@stats$platform)) as.character(scope_obj@stats$platform) else NULL
  )
  is_visium <- any(grepl("visium", platform_vals, ignore.case = TRUE))
  is_xenium <- any(grepl("xenium", platform_vals, ignore.case = TRUE))

  resolve_hex <- function(force = FALSE, msg_prefix = NULL) {
    hex_info <- .spw_infer_hex_variant(grid_info, force = force)
    if (!is.null(msg_prefix) && verbose) {
      message(msg_prefix, hex_info$variant)
    } else if (force && verbose && isFALSE(hex_info$hex_like)) {
      message("[geneSCOPE::computeWeights] Hex topology requested; defaulting to ", hex_info$variant, " due to ambiguous geometry.")
    } else if (verbose && isTRUE(hex_info$hex_like)) {
      message("[geneSCOPE::computeWeights] Auto-detected hex topology: ", hex_info$variant)
    }
    list(topology = hex_info$variant, base = "hex")
  }

  if (forced_topology) {
    if (identical(topology, "hex")) {
      return(resolve_hex(force = TRUE))
    }
    if (verbose) message("[geneSCOPE::computeWeights] Topology forced by argument: ", topology)
    return(list(topology = topology, base = topology))
  }

  if (is_xenium) {
    if (verbose) message("[geneSCOPE::computeWeights] Platform flagged as Xenium; forcing queen adjacency.")
    return(list(topology = "queen", base = "queen"))
  }

  if (is_visium && has_centers) {
    return(resolve_hex(force = TRUE, msg_prefix = "[geneSCOPE::computeWeights] Platform flagged as Visium; selecting hex orientation: "))
  }

  if (has_centers) {
    detected <- .spw_detect_topology(grid_info, verbose)
    return(detected)
  }

  if (verbose) {
    message("[geneSCOPE::computeWeights] Defaulting to queen adjacency (insufficient geometry for auto detection).")
  }
  list(topology = "queen", base = "queen")
}

.spw_detect_topology <- function(grid_info, verbose) {
  hex_info <- .spw_infer_hex_variant(grid_info, force = FALSE)
  if (isTRUE(hex_info$hex_like)) {
    return(list(topology = hex_info$variant, base = "hex"))
  }
  if (verbose) message("[geneSCOPE::computeWeights] Auto-detected topology: queen")
  list(topology = "queen", base = "queen")
}

.spw_infer_hex_variant <- function(grid_info, force = FALSE) {
  gi <- data.table::as.data.table(grid_info)[, .(gx, gy, center_x, center_y)]
  analysis <- .spw_analyse_hex_geometry(gi)
  hex_like <- isTRUE(analysis$hex_like)

  orientation <- if (hex_like || force) {
    if (is.finite(analysis$horiz) && is.finite(analysis$vert)) {
      if (analysis$horiz >= analysis$vert) "row" else "col"
    } else {
      "row"
    }
  } else {
    NA_character_
  }

  if (is.na(orientation)) {
    return(list(hex_like = hex_like, variant = "hex-oddr", orientation = orientation))
  }

  variant <- if (orientation == "row") {
    .spw_hex_row_parity(gi)
  } else {
    .spw_hex_col_parity(gi)
  }

  list(hex_like = hex_like, variant = variant, orientation = orientation)
}

.spw_analyse_hex_geometry <- function(gi) {
  mk <- function(dr, dc) {
    b <- gi[, .(gx_nb = gx - dr, gy_nb = gy - dc, cx_nb = center_x, cy_nb = center_y)]
    m <- merge(gi, b, by.x = c("gx", "gy"), by.y = c("gx_nb", "gy_nb"), all = FALSE)
    if (!nrow(m)) return(NULL)
    data.table::data.table(dx = m$cx_nb - m$center_x, dy = m$cy_nb - m$center_y)
  }
  edges_xy <- data.table::rbindlist(list(
    mk(0, 1), mk(1, 0), mk(1, 1), mk(1, -1)
  ), use.names = TRUE, fill = TRUE)

  if (is.null(edges_xy) || !nrow(edges_xy)) {
    return(list(hex_like = FALSE, horiz = NA_real_, vert = NA_real_, diag_ratio = NA_real_))
  }

  theta <- abs(atan2(edges_xy$dy, edges_xy$dx) * 180 / pi)
  horiz <- mean(pmin(theta, 180 - theta) <= 10, na.rm = TRUE)
  vert <- mean(abs(theta - 90) <= 10, na.rm = TRUE)
  hex_hits <- mean(abs(theta - 60) <= 10 | abs(theta - 120) <= 10, na.rm = TRUE)
  hex_hits <- if (is.finite(hex_hits)) hex_hits else 0

  neighbor_pairs <- function(dr, dc) {
    b2 <- gi[, .(gx_nb = gx - dr, gy_nb = gy - dc, cx_nb = center_x, cy_nb = center_y)]
    m2 <- merge(gi, b2, by.x = c("gx", "gy"), by.y = c("gx_nb", "gy_nb"), all = FALSE)
    if (!nrow(m2)) return(numeric())
    sqrt((m2$cx_nb - m2$center_x)^2 + (m2$cy_nb - m2$center_y)^2)
  }
  d_ax <- c(neighbor_pairs(0, 1), neighbor_pairs(1, 0))
  d_diag <- neighbor_pairs(1, 1)
  med_ax <- if (length(d_ax)) stats::median(d_ax, na.rm = TRUE) else NA_real_
  med_diag <- if (length(d_diag)) stats::median(d_diag, na.rm = TRUE) else NA_real_
  diag_ratio <- if (is.finite(med_ax) && med_ax > 0 && is.finite(med_diag)) med_diag / med_ax else NA_real_
  ratio_ok <- is.finite(diag_ratio) && abs(diag_ratio - 1) <= 0.2

  list(
    hex_like = (hex_hits >= 0.4) && max(horiz, vert, na.rm = TRUE) < 0.8 && ratio_ok,
    horiz = if (is.finite(horiz)) horiz else 0,
    vert = if (is.finite(vert)) vert else 0,
    diag_ratio = if (is.finite(diag_ratio)) diag_ratio else NA_real_
  )
}

.spw_hex_row_parity <- function(gi) {
  row_diffs <- gi[order(gy, center_x), diff(center_x), by = gy]$V1
  row_diffs <- row_diffs[is.finite(row_diffs)]
  if (!length(row_diffs)) return("hex-oddr")
  cell_width <- stats::median(abs(row_diffs), na.rm = TRUE)
  if (!is.finite(cell_width) || cell_width == 0) return("hex-oddr")
  half_width <- cell_width / 2

  row_shift <- gi[, .(min_x = min(center_x)), by = gy]
  global_min <- min(row_shift$min_x)
  row_shift[, shift := min_x - global_min]
  row_shift[, ratio := shift / half_width]
  row_shift[, nearest := round(ratio)]
  row_shift[, diff := abs(ratio - nearest)]
  row_shift[diff > 0.35, nearest := NA_integer_]
  row_shift[, flag := nearest %% 2]

  odd_mean <- mean(row_shift[gy %% 2 == 1L, flag], na.rm = TRUE)
  even_mean <- mean(row_shift[gy %% 2 == 0L, flag], na.rm = TRUE)

  if (is.nan(odd_mean) && is.nan(even_mean)) return("hex-oddr")
  if (!is.nan(odd_mean) && (is.nan(even_mean) || odd_mean >= 0.5)) return("hex-oddr")
  if (!is.nan(even_mean) && even_mean >= 0.5) return("hex-evenr")

  if (!is.nan(odd_mean) && odd_mean > 0.25) return("hex-oddr")
  "hex-evenr"
}

.spw_hex_col_parity <- function(gi) {
  col_diffs <- gi[order(gx, center_y), diff(center_y), by = gx]$V1
  col_diffs <- col_diffs[is.finite(col_diffs)]
  if (!length(col_diffs)) return("hex-oddq")
  cell_height <- stats::median(abs(col_diffs), na.rm = TRUE)
  if (!is.finite(cell_height) || cell_height == 0) return("hex-oddq")
  half_height <- cell_height / 2

  col_shift <- gi[, .(min_y = min(center_y)), by = gx]
  global_min <- min(col_shift$min_y)
  col_shift[, shift := min_y - global_min]
  col_shift[, ratio := shift / half_height]
  col_shift[, nearest := round(ratio)]
  col_shift[, diff := abs(ratio - nearest)]
  col_shift[diff > 0.35, nearest := NA_integer_]
  col_shift[, flag := nearest %% 2]

  odd_mean <- mean(col_shift[gx %% 2 == 1L, flag], na.rm = TRUE)
  even_mean <- mean(col_shift[gx %% 2 == 0L, flag], na.rm = TRUE)

  if (is.nan(odd_mean) && is.nan(even_mean)) return("hex-oddq")
  if (!is.nan(odd_mean) && (is.nan(even_mean) || odd_mean >= 0.5)) return("hex-oddq")
  if (!is.nan(even_mean) && even_mean >= 0.5) return("hex-evenq")

  if (!is.nan(odd_mean) && odd_mean > 0.25) return("hex-oddq")
  "hex-evenq"
}

.spw_build_neighbourhood <- function(meta, topo_info, zero_policy) {
  topology <- topo_info$topology
  xbins_eff <- meta$xbins_eff
  ybins_eff <- meta$ybins_eff
  gx <- meta$gx
  gy <- meta$gy

  if (topology %in% c("queen", "rook")) {
    full_nb <- geneSCOPE:::grid_nb_omp(ybins_eff, xbins_eff, queen = identical(topology, "queen"))
    cell_id <- (gy - 1L) * xbins_eff + gx
  } else if (topology == "hex-oddr") {
    full_nb <- geneSCOPE:::grid_nb_hex_omp(ybins_eff, xbins_eff, oddr = TRUE)
    cell_id <- (gy - 1L) * xbins_eff + gx
  } else if (topology == "hex-evenr") {
    full_nb <- geneSCOPE:::grid_nb_hex_omp(ybins_eff, xbins_eff, oddr = FALSE)
    cell_id <- (gy - 1L) * xbins_eff + gx
  } else if (topology == "hex-oddq") {
    full_nb <- geneSCOPE:::grid_nb_hexq_omp(ybins_eff, xbins_eff, oddq = TRUE)
    cell_id <- (gx - 1L) * ybins_eff + gy
  } else if (topology == "hex-evenq") {
    full_nb <- geneSCOPE:::grid_nb_hexq_omp(ybins_eff, xbins_eff, oddq = FALSE)
    cell_id <- (gx - 1L) * ybins_eff + gy
  } else {
    stop("Unknown topology: ", topology)
  }

  sub_nb <- full_nb[cell_id]
  relabel_nb <- lapply(sub_nb, function(v) {
    idx <- match(v, cell_id)
    idx[!is.na(idx)]
  })

  attr(relabel_nb, "class") <- "nb"
  attr(relabel_nb, "region.id") <- meta$grid_id
  attr(relabel_nb, "queen") <- identical(topology, "queen")
  attr(relabel_nb, "topology") <- topology

  zero_idx <- integer(0)
  if (zero_policy) {
    zero_idx <- which(vapply(relabel_nb, length, integer(1)) == 0L)
    if (length(zero_idx) > 0L) {
      for (j in zero_idx) {
        relabel_nb[[j]] <- j
      }
    }
  }

  list(nb = relabel_nb, zero_idx = zero_idx)
}

# ---------------------------------------------------------------------------
# Fuzzy Queen utilities
# ---------------------------------------------------------------------------

#' Row-normalise a sparse adjacency matrix.
#'
#' @param W \code{dgCMatrix} adjacency matrix.
#' @return A row-stochastic \code{dgCMatrix}.
#' @keywords internal
.row_normalize_sparse <- function(W) {
  rs <- Matrix::rowSums(W)
  inv_rs <- rep(0, length(rs))
  inv_rs[rs > 0] <- 1 / rs[rs > 0]
  as(Matrix::Diagonal(x = inv_rs) %*% W, "dgCMatrix")
}

#' Diffuse a binary pattern through a spatial adjacency matrix.
#'
#' @param x Numeric vector (0/1) of length \eqn{n}.
#' @param W \code{dgCMatrix} adjacency (optionally row-normalised).
#' @param alpha Diffusion weight per step.
#' @param steps Number of diffusion steps.
#' @param clamp Logical; clamp values into \eqn{[0, 1]} after each step.
#' @return Numeric vector after diffusion.
#' @keywords internal
.diffuse_binary_pattern <- function(x, W, alpha = 0.5, steps = 1L, clamp = TRUE) {
  v <- as.numeric(x)
  if (steps <= 0L || alpha <= 0) return(v)
  for (s in seq_len(steps)) {
    v <- v + alpha * as.numeric(W %*% v)
    if (clamp) v <- pmin(1, pmax(0, v))
  }
  v
}

#' Diffuse multiple binary patterns simultaneously.
#'
#' @param X Matrix (n_cells × n_patterns) of 0/1 values.
#' @param W \code{dgCMatrix} adjacency (optionally row-normalised).
#' @param alpha Diffusion weight per step.
#' @param steps Number of diffusion steps.
#' @return Numeric matrix after diffusion.
#' @keywords internal
.fuzzy_queen_diffuse_mat <- function(X, W, alpha = 0.5, steps = 1L) {
  V <- as.matrix(X)
  if (steps <= 0L || alpha <= 0) return(V)
  for (s in seq_len(steps)) {
    V <- V + alpha * as.matrix(W %*% V)
    V <- pmin(1, pmax(0, V))
  }
  V
}

#' Fuzzy-queen Jaccard similarity between two binary patterns
#'
#' @param x,y Numeric/logical vectors (0/1) of identical length.
#' @param W Sparse \code{dgCMatrix} queen adjacency.
#' @param alpha Diffusion weight per step (0–1).
#' @param steps Integer number of diffusion steps (≥ 0).
#' @param row_normalize Logical; row-normalise \code{W} before diffusion.
#' @return Scalar numeric fuzzy Jaccard similarity in \eqn{[0, 1]} (or \code{NA}
#'   when both patterns are empty).
#' @export
fuzzy_queen_jaccard <- function(x, y,
                                W,
                                alpha = 0.5,
                                steps = 1L,
                                row_normalize = TRUE) {
  if (length(x) != length(y)) stop("x and y must have identical length.")
  n <- length(x)
  if (!inherits(W, "dgCMatrix")) stop("W must be a dgCMatrix.")
  if (nrow(W) != n || ncol(W) != n) stop("Dimensions of W must match length of x and y.")
  if (!is.numeric(alpha) || alpha < 0) stop("alpha must be non-negative.")
  if (!is.numeric(steps) || steps < 0) stop("steps must be a non-negative integer.")

  x_vec <- as.numeric(x > 0)
  y_vec <- as.numeric(y > 0)

  W_use <- if (isTRUE(row_normalize)) .row_normalize_sparse(W) else W

  x_fuzzy <- .diffuse_binary_pattern(x_vec, W_use, alpha = alpha, steps = steps)
  y_fuzzy <- .diffuse_binary_pattern(y_vec, W_use, alpha = alpha, steps = steps)

  num <- sum(pmin(x_fuzzy, y_fuzzy))
  den <- sum(pmax(x_fuzzy, y_fuzzy))
  if (den == 0) return(NA_real_)
  as.numeric(num / den)
}

#' Diffuse multiple binary patterns with fuzzy-queen smoothing
#'
#' @param X Matrix/data.frame whose rows correspond to cells (length equals
#'   \code{nrow(W)}). Values should be 0/1.
#' @param W Sparse \code{dgCMatrix} adjacency.
#' @param alpha Diffusion weight per step.
#' @param steps Number of diffusion steps.
#' @param row_normalize Logical; row-normalise \code{W} prior to diffusion.
#' @return Numeric matrix with the same dimensions as \code{X}.
#' @export
fuzzy_queen_diffuse_mat <- function(X, W, alpha = 0.5, steps = 1L, row_normalize = TRUE) {
  if (!inherits(W, "dgCMatrix")) stop("W must be a dgCMatrix.")
  V <- as.matrix(X)
  if (nrow(W) != nrow(V)) stop("Row count of X must match dimension of W.")
  W_use <- if (isTRUE(row_normalize)) .row_normalize_sparse(W) else W
  .fuzzy_queen_diffuse_mat(V, W_use, alpha = alpha, steps = steps)
}
