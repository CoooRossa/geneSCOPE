#' @title Build CoordObj from Xenium Output with (Multi‑)Scale Grid Aggregation
#' @description
#'   Streams Xenium *transcripts.parquet* once, applies molecule‑level QC, and **creates several
#'   grid layers at different `lenGrid` resolutions in a single pass**.  Segmentation polygons and
#'   centroids are read once and attached unchanged.  Molecule table is discarded to save RAM.
#'
#'   For each `lenGrid[i]`, a layer named `grid_lenGrid<value>` is added under `coordObj@grid`,
#'   containing `grid_info`, `counts`, and helper metadata.
#'
#'   Additional filters include gene diversity thresholds (`min_gene_types`, `max_gene_types`)
#'   and minimum segmentation points (`min_seg_points`) per grid cell.
#'
#' @param xenium_dir   Path to Xenium output folder.
#' @param lenGrid      Numeric **vector**. Each value is a square grid side length (µm).
#' @param filtergenes, max_dist_mol_nuc, filtermolecule, filterqv, coord_file, ncores, verbose
#'                     Same meaning as in v3 (single‑resolution version).
#' @param min_gene_types Minimum number of unique gene types required per grid cell (default 1).
#' @param max_gene_types Maximum number of unique gene types allowed per grid cell (default Inf).
#' @param min_seg_points Minimum number of segmentation points required per grid cell (default 0).
#' @param seg_type     Character. One of "cell", "nucleus", or "both" (default "cell").
#' @param flip_y Logical. Mirror‑flip Y axis if TRUE (default FALSE) so grid aligns
#'        with Xenium Explorer (origin top‑left).
#'
#' @return An S4 `CoordObj` with one grid sub‑layer per `lenGrid` value, plus `centroids` and
#'         `segmentation` in `coord`.
#'
#' @import arrow
#' @import data.table
#' @import sf
#' @import parallel
#' @export

createCoordObj<- function(xenium_dir,
                                          lenGrid,
                                          seg_type          = c("cell", "nucleus", "both"),
                                          filtergenes       = NULL,
                                          max_dist_mol_nuc  = NULL,
                                          filtermolecule    = TRUE,
                                          filterqv          = NULL,
                                          coord_file        = NULL,
                                          ncores            = 1L,
                                          chunk_pts         = 5e5,
                                          min_gene_types    = 1,
                                          max_gene_types    = Inf,
                                          min_seg_points    = 0,
                                          verbose           = TRUE,
                                          flip_y            = FALSE) {
  seg_type <- match.arg(seg_type)
  stopifnot(is.logical(flip_y) && length(flip_y) == 1)
  stopifnot(dir.exists(xenium_dir))
  stopifnot(is.numeric(lenGrid) && length(lenGrid) >= 1 && all(lenGrid > 0))
  options(arrow.num_threads = ncores)

  library(arrow)
  library(data.table)
  library(sf)
  library(parallel)

# ---- safe S4 registration helpers ---------------------------------------
safeSetClass <- function(class, slots, ...) {
  if (is.null(methods::getClassDef(class, where = topenv()))) {
    methods::setClass(class, slots = slots, ...)
  }
}

safeSetMethod <- function(generic, signature, definition, ...) {
  if (is.null(methods::selectMethod(generic, signature, optional = TRUE))) {
    methods::setMethod(generic, signature = signature, definition = definition, ...)
  }
}

  tf_parq   <- file.path(xenium_dir, "transcripts.parquet")
  cell_parq <- file.path(xenium_dir, "cells.parquet")
  if (!file.exists(tf_parq) || !file.exists(cell_parq)) stop("transcripts.parquet or cells.parquet missing")

  # ---- locate segmentation parquet(s) ------------------------------------
  seg_files <- list()
  if (seg_type %in% c("cell", "both")) {
    cell_file <- if (file.exists(file.path(xenium_dir, "cell_boundaries.parquet"))) {
      "cell_boundaries.parquet"
    } else if (file.exists(file.path(xenium_dir, "segmentation_boundaries.parquet"))) {
      "segmentation_boundaries.parquet"
    } else {
      stop("cell segmentation file missing")
    }
    seg_files$cell <- file.path(xenium_dir, cell_file)
  }
  if (seg_type %in% c("nucleus", "both")) {
    nuc_file <- file.path(xenium_dir, "nucleus_boundaries.parquet")
    if (!file.exists(nuc_file)) stop("nucleus_boundaries.parquet missing")
    seg_files$nucleus <- nuc_file
  }

  # ---- centroids ---------------------------------------------------------
  centroids_dt <- as.data.table(read_parquet(cell_parq))[, .(cell = cell_id, x = x_centroid, y = y_centroid)]

  # ---- build dataset with QC filters (dplyr verbs compile to Arrow C++) ----
  ds <- open_dataset(tf_parq)
  if (!is.null(filtergenes))      ds <- ds |> dplyr::filter(feature_name %in% filtergenes)
  if (!is.null(max_dist_mol_nuc)) ds <- ds |> dplyr::filter(nucleus_distance <= max_dist_mol_nuc)
  if (filtermolecule)             ds <- ds |> dplyr::filter(!grepl("^(Unassigned|NegControl|Background)", feature_name))
  if (!is.null(filterqv))         ds <- ds |> dplyr::filter(qv >= filterqv)

  # ---- capture global y_max after bounds_global creation -----------------
  bounds_global <- ds |> summarise(xmin = min(x_location), xmax = max(x_location),
                                   ymin = min(y_location), ymax = max(y_location)) |> collect()
  y_max <- as.numeric(bounds_global$ymax)

  # ---- mirror centroids if flip_y ----------------------------------------
  if (flip_y) centroids_dt[, y := y_max - y]

  # ---- optional user polygon --------------------------------------------
  user_poly <- NULL
  if (!is.null(coord_file)) {
    ext <- tolower(tools::file_ext(coord_file))
    sep_char <- if (ext == "csv") "," else if (ext %in% c("tsv", "txt")) "\t" else "[,\t ]"
    poly_tb  <- utils::read.table(coord_file, header = TRUE, sep = sep_char)
    names(poly_tb) <- tolower(names(poly_tb))
    stopifnot(all(c("x", "y") %in% names(poly_tb)))
    if (flip_y) poly_tb$y <- y_max - poly_tb$y
    user_poly <- st_sfc(st_polygon(list(as.matrix(poly_tb[, c("x", "y")]))))
  }

  # ---- polygon‑based clipping of centroids & segmentation ----------------
  if (!is.null(user_poly)) {
    # centroids inside ROI
    if (nrow(centroids_dt)) {
      centroids_sf <- st_as_sf(centroids_dt, coords = c("x", "y"), crs = NA, remove = FALSE)
      keep_cent    <- lengths(st_within(centroids_sf, user_poly)) > 0
      centroids_dt <- centroids_dt[keep_cent]
    }
    keep_cells <- unique(centroids_dt$cell)
    if (length(keep_cells) == 0) warning("ROI polygon removed all centroids; downstream layers empty")
  } else {
    keep_cells <- NULL
  }

  # ---- helper: clip points to polygon in parallel ------------------------
  clip_poly_parallel <- function(dt_pts) {
    if (is.null(user_poly) || !nrow(dt_pts)) return(dt_pts)
    if (!is.data.table(dt_pts)) dt_pts <- as.data.table(dt_pts)

    n_tot <- nrow(dt_pts)
    idx_split <- split(seq_len(n_tot), ceiling(seq_len(n_tot) / chunk_pts))

    worker <- function(ix) {
      sub <- dt_pts[ix]
      keep <- lengths(st_within(
        st_as_sf(sub, coords = c("x", "y"), crs = NA, remove = FALSE),
        user_poly)) > 0
      sub[keep]
    }

    res <- if (ncores > 1 && .Platform$OS.type != "windows") {
      mclapply(idx_split, worker, mc.cores = ncores)
    } else if (ncores > 1) {
      cl <- makeCluster(ncores)
      clusterEvalQ(cl, library(sf))
      clusterExport(cl, varlist = c("dt_pts", "user_poly", "worker"), envir = environment())
      out <- parLapply(cl, idx_split, worker)
      stopCluster(cl)
      out
    } else {
      lapply(idx_split, worker)
    }
    rbindlist(res)
  }

 ## optional pre‑collect molecules for polygon ---------------------------
  mol_small <- NULL
  if (!is.null(user_poly)) {
    mol_small <- ds |> select(x = x_location, y = y_location, feature_name) |> collect() |> as.data.table()
    if (flip_y) mol_small[, y := y_max - y]
    mol_small <- clip_poly_parallel(mol_small)
  }
  coord_obj <- new("CoordObj", coord = list(centroids = centroids_dt), grid = list(), meta.data = data.frame())

  ## read + convert segmentation layers ----------------------------------
  convert_seg <- function(file, tag) {
     seg_dt <- as.data.table(read_parquet(file))[
              , .(cell = cell_id, x = vertex_x, y = vertex_y, label_id)]
    if (flip_y) seg_dt[, y := y_max - y]
    if (!is.null(keep_cells)) seg_dt <- seg_dt[cell %in% keep_cells]

    coord_obj@coord[[paste0("segmentation_", tag)]] <<- seg_dt
    split_idx <- split(seq_len(nrow(seg_dt)), seg_dt$label_id)
    build_poly <- function(ix) {
      df <- seg_dt[ix, .(x,y)]; if (nrow(df) < 3) return(NULL); st_polygon(list(as.matrix(df))) }
    plist <- if (ncores > 1 && .Platform$OS.type != "windows") {
      mclapply(split_idx, build_poly, mc.cores = ncores)
    } else if (ncores > 1) {
      cl <- makeCluster(ncores); clusterEvalQ(cl, library(sf));
      clusterExport(cl, varlist = c("seg_dt", "build_poly"), envir = environment())
      out <- parLapply(cl, split_idx, build_poly); stopCluster(cl); out
    } else lapply(split_idx, build_poly)
    coord_obj@coord[[paste0("segmentation_polys_", tag)]] <<- st_sfc(plist[!sapply(plist, is.null)])
  }
  lapply(names(seg_files), \(tag) convert_seg(seg_files[[tag]], tag))

  ## grid layer builder ----------------------------------------------------
  build_grid <- function(lg) {
  if (verbose) message("• Building grid lenGrid = ", lg)

  # --- aligned origin so that breaks are multiples of lg -----------------
  # bounds_global already computed above
  x0 <- floor(bounds_global$xmin / lg) * lg
  y0 <- floor(bounds_global$ymin / lg) * lg

  if (!is.null(user_poly)) {
    # ---- user-defined polygon分支 ----------------------------------------
    mol_dt <- copy(mol_small) |>
              mutate(gx  = (x - x0) %/% lg,
                     gy  = ((if (flip_y) (y_max - y) else y) - y0) %/% lg,
                     grid_id = paste0("g", gx, "_", gy))

    cnt <- mol_dt |>
      as.data.table() |>                     # 转 data.table
      (\(dt) dt[, .(count = .N), by = .(grid_id, feature_name)])()   # 匿名函数里用 [...]
    setnames(cnt, "feature_name", "gene")
  } else {
    # ---- 无 polygon，直接 Arrow C++ 聚合 ---------------------------------
    cnt <- ds |>
      mutate(gx  = (x_location - x0) %/% lg,
             gy  = ((if (flip_y) (y_max - y_location) else y_location) - y0) %/% lg,
             grid_id = paste0("g", gx, "_", gy)) |>
      group_by(grid_id, feature_name) |>
      summarise(count = n(), .groups = "drop") |>
      collect() |>
      as.data.table()
    setnames(cnt, "feature_name", "gene")
  }

  # ---- filter by gene diversity ---------------------------------------
  gene_div <- cnt[, uniqueN(gene), by = grid_id]
  keep_grid <- gene_div[V1 >= min_gene_types & V1 <= max_gene_types, grid_id]
  cnt <- cnt[grid_id %in% keep_grid]

  if (min_seg_points > 0) {
    seg_layers <- names(coord_obj@coord)[grepl("^segmentation_", names(coord_obj@coord)) &
                                         !grepl("_polys_", names(coord_obj@coord))]
    seg_dt_all <- rbindlist(coord_obj@coord[seg_layers], use.names = TRUE, fill = TRUE)
    head(seg_dt_all)  # for debugging
    seg_dt_all[, `:=`(gx = (x - x0) %/% lg,
                      gy = ((if (flip_y) (y_max - y) else y) - y0) %/% lg)]
    seg_dt_all[, grid_id := paste0("g", gx, "_", gy)]
    seg_cnt <- seg_dt_all[, .N, by = grid_id]
    valid_seg <- seg_cnt[N >= min_seg_points, grid_id]
    keep_grid <- intersect(keep_grid, valid_seg)
    cnt <- cnt[grid_id %in% keep_grid]
  }

  # ---- 计算全局 bounds ----------------------------------------------------
  # bounds_global already computed above

    x_breaks <- seq(x0, bounds_global$xmax, by = lg)
    if (tail(x_breaks, 1) < bounds_global$xmax) x_breaks <- c(x_breaks, bounds_global$xmax)
    y_breaks <- seq(y0, bounds_global$ymax, by = lg)
    if (tail(y_breaks, 1) < bounds_global$ymax) y_breaks <- c(y_breaks, bounds_global$ymax)
    xbins_eff <- length(x_breaks) - 1L; ybins_eff <- length(y_breaks) - 1L

    # ---- build grid_dt --------------------------------------------------
    grid_dt <- CJ(gx = 0:(xbins_eff - 1L), gy = 0:(ybins_eff - 1L))

    ## ① 先写入边界列——这些列彼此独立
    grid_dt[, `:=`(
      xmin = x_breaks[gx + 1L],
      xmax = x_breaks[gx + 2L],
      ymin = y_breaks[gy + 1L],
      ymax = y_breaks[gy + 2L]
    )]

    ## ② 再用已存在的列计算派生量
    grid_dt[, `:=`(
      center_x = (xmin + xmax) / 2,
      center_y = (ymin + ymax) / 2,
      width    = xmax - xmin,
      height   = ymax - ymin,
      grid_id  = paste0("g", gx, "_", gy),
      idx      = .I
    )]

    grid_dt <- grid_dt[grid_id %in% keep_grid]

    # remove edge grids whose area is not exactly lg*lg
    grid_dt <- grid_dt[abs(width - lg) < 1e-6 & abs(height - lg) < 1e-6]
    # ---- optional Y‑axis mirror flip -----------------------------------
    if (flip_y) {
      y_max <- as.numeric(bounds_global$ymax)
      grid_dt[, `:=`(ymin_tmp = ymin, ymax_tmp = ymax)]
      grid_dt[, `:=`(ymin = y_max - ymax_tmp,
                     ymax = y_max - ymin_tmp)]
      grid_dt[, center_y := (ymin + ymax) / 2]
      grid_dt[, c("ymin_tmp", "ymax_tmp") := NULL]
    }
    keep_grid <- intersect(keep_grid, grid_dt$grid_id)
    cnt <- cnt[grid_id %in% keep_grid]

    coord_obj@grid[[paste0("grid_lenGrid", lg)]] <<- list(
      grid_info = grid_dt,
      counts    = cnt,
      lenGrid   = lg,
      xbins_eff = xbins_eff,
      ybins_eff = ybins_eff)
}

  lapply(sort(unique(lenGrid)), build_grid)
  coord_obj
}

  # Ensure only one build_grid definition remains (remove any below this point)