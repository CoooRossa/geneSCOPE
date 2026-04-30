#' Detect Darwin Native Spatial Default
#' @description
#' Native grid/listw helpers can terminate the R process on some Darwin builds.
#' Keep Lee's L and other numerical native paths available, but route spatial
#' neighbourhood/listw construction through R unless explicitly opted back in.
#' @return Logical flag.
#' @keywords internal
.spw_darwin_native_spatial_disabled <- function() {
  sysname <- tolower(Sys.info()[["sysname"]] %||% "")
  if (!identical(sysname, "darwin")) return(FALSE)
  if (exists(".native_all_explicitly_enabled", mode = "function") &&
      .native_all_explicitly_enabled()) return(FALSE)
  if (isTRUE(getOption("geneSCOPE.allow_darwin_native_spatial", FALSE))) return(FALSE)
  env_val <- tolower(Sys.getenv("GENESCOPE_ALLOW_DARWIN_NATIVE_SPATIAL", ""))
  !(env_val %in% c("1", "true", "yes", "on", "force", "unsafe"))
}

.spw_darwin_native_spatial_reason <- function() {
  "darwin_native_spatial_safe_default"
}

#' Spw Should Enable Parallel Nb
#' @description
#' Internal helper for `.spw_should_enable_parallel_nb`.
#' @param ncores Number of cores/threads to use.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.spw_should_enable_parallel_nb <- function(ncores, verbose,
                                           parent = "computeWeights",
                                           step = "S03") {
  if (ncores < 2L) return(FALSE)
  if (.spw_darwin_native_spatial_disabled()) {
    if (verbose) {
      .log_info(parent, step, "Darwin native spatial OpenMP builder disabled by safe default; using R fallback where available.", verbose)
    }
    return(FALSE)
  }
  if (.native_backend_disabled("geneSCOPE.disable_native_grid_nb")) {
    if (verbose) {
      .log_info(parent, step, "Native OpenMP grid builder disabled by option; using serial/R fallback path.", verbose)
    }
    return(FALSE)
  }
  if (tolower(Sys.getenv("GENESCOPE_FORCE_R_GRID_NB", "")) %in% c("1", "true", "yes", "on")) {
    return(FALSE)
  }
  env_val <- Sys.getenv("GENESCOPE_ENABLE_OMP_GRID_NB", "")
  if (!nzchar(env_val)) return(TRUE)
  env_lower <- tolower(env_val)
  requested <- env_lower %in% c("1", "true", "yes", "on", "force", "unsafe")
  if (!requested) return(FALSE)
  force_enable <- env_lower %in% c("force", "unsafe")
  sysname <- Sys.info()[["sysname"]]
  if (is.null(sysname)) sysname <- ""
  sysname <- tolower(sysname)
  if (identical(sysname, "darwin") && force_enable && verbose) {
    .log_info(parent, step, "OpenMP grid builder forced on macOS (GENESCOPE_ENABLE_OMP_GRID_NB=force). Expect instability.", verbose)
  }
  TRUE
}

#' Spw Build Neighbourhood
#' @description
#' Internal helper for `.spw_build_neighbourhood`.
#' @param meta Parameter value.
#' @param topo_info Parameter value.
#' @param ncores Number of cores/threads to use.
#' @param verbose Logical; whether to emit progress messages.
#' @return Return value used internally.
#' @keywords internal
.spw_build_neighbourhood <- function(meta, topo_info, ncores = 1L, verbose = FALSE,
                                     parent = "computeWeights", step = "S03") {
  topology <- topo_info$topology
  dims <- .spw_validate_grid_shape(meta$xbins_eff, meta$ybins_eff, context = "[geneSCOPE::.compute_weights]")
  xbins_eff <- dims$xbins_eff
  ybins_eff <- dims$ybins_eff
  n_cells <- dims$n_cells
  gx <- meta$gx
  gy <- meta$gy
  if (length(gx) != length(gy)) {
    stop("Grid metadata gx/gy lengths do not match (", length(gx), " vs ", length(gy), ").")
  }
  if (any(!is.finite(gx)) || any(!is.finite(gy))) {
    stop("Grid metadata gx/gy must be finite.")
  }

  use_parallel <- .spw_should_enable_parallel_nb(ncores, verbose,
    parent = parent, step = step
  )

  if (.spw_darwin_native_spatial_disabled() && !topology %in% c("queen", "rook")) {
    stop(
      "Native ", topology, " spatial builder is disabled by Darwin spatial safe default, ",
      "and no R fallback is available for this topology. Set ",
      "options(geneSCOPE.allow_darwin_native_spatial=TRUE) in a clean R session ",
      "only if you accept native spatial helper risk.",
      call. = FALSE
    )
  }

  if (use_parallel && verbose) {
    .log_info(parent, step, "Using native OpenMP grid builder by default (C++ first, R fallback).", verbose)
  } else if (!use_parallel && verbose && ncores > 1L) {
    .log_info(parent, step, "Using serial grid builder (OpenMP gate disabled).", verbose)
  }

  if (.spw_trace_enabled()) {
    gx_range <- range(meta$gx, na.rm = TRUE)
    gy_range <- range(meta$gy, na.rm = TRUE)
    .spw_trace(sprintf(
      "grid[%s] builder: topology=%s dims=%d x %d (cells=%d) active_points=%d xbins_integerish=%s ybins_integerish=%s grid_length=%s gx_range=[%s,%s] gy_range=[%s,%s] backend=%s ncores=%d",
      if (!is.null(meta$grid_name)) meta$grid_name else "<unknown>",
      topology,
      ybins_eff,
      xbins_eff,
      n_cells,
      length(meta$gx),
      if (!is.null(meta$xbins_integerish)) meta$xbins_integerish else NA,
      if (!is.null(meta$ybins_integerish)) meta$ybins_integerish else NA,
      if (!is.null(meta$grid_length)) format(meta$grid_length, digits = 12, trim = TRUE) else NA_character_,
      format(gx_range[1], digits = 8, trim = TRUE),
      format(gx_range[2], digits = 8, trim = TRUE),
      format(gy_range[1], digits = 8, trim = TRUE),
      format(gy_range[2], digits = 8, trim = TRUE),
      if (use_parallel) "omp" else "serial",
      as.integer(ncores)
    ), parent = parent, step = step, verbose = verbose)
  }

  choose_builder <- function(rect_fn, hexr_fn, hexq_fn) {
    if (topology %in% c("queen", "rook")) {
      if (use_parallel) rect_fn$omp else rect_fn$serial
    } else if (topology %in% c("hex-oddr", "hex-evenr")) {
      if (use_parallel) hexr_fn$omp else hexr_fn$serial
    } else if (topology %in% c("hex-oddq", "hex-evenq")) {
      if (use_parallel) hexq_fn$omp else hexq_fn$serial
    } else {
      stop("Unknown topology: ", topology)
    }
  }

  builder <- choose_builder(
    rect_fn = list(
      omp = function() .grid_nb_safe(
        ybins_eff,
        xbins_eff,
        topology = topology,
        queen = identical(topology, "queen"),
        use_parallel = use_parallel,
        verbose = verbose,
        parent = parent,
        step = step
      ),
      serial = function() .grid_nb_safe(
        ybins_eff,
        xbins_eff,
        topology = topology,
        queen = identical(topology, "queen"),
        use_parallel = FALSE,
        verbose = verbose,
        parent = parent,
        step = step
      )
    ),
    hexr_fn = list(
      omp = function() list(
        nb = .grid_nb_hex_omp(ybins_eff, xbins_eff, oddr = identical(topology, "hex-oddr")),
        backend = if (use_parallel) "omp" else "serial"
      ),
      serial = function() list(
        nb = .grid_nb_hex(ybins_eff, xbins_eff, oddr = identical(topology, "hex-oddr")),
        backend = "serial"
      )
    ),
    hexq_fn = list(
      omp = function() list(
        nb = .grid_nb_hexq_omp(ybins_eff, xbins_eff, oddq = identical(topology, "hex-oddq")),
        backend = if (use_parallel) "omp" else "serial"
      ),
      serial = function() list(
        nb = .grid_nb_hexq(ybins_eff, xbins_eff, oddq = identical(topology, "hex-oddq")),
        backend = "serial"
      )
    )
  )

  builder_res <- builder()
  full_nb <- builder_res$nb
  backend_used <- builder_res$backend
  if (topology %in% c("queen", "rook", "hex-oddr", "hex-evenr")) {
    cell_id <- (gy - 1L) * xbins_eff + gx
  } else {
    cell_id <- (gx - 1L) * ybins_eff + gy
  }

  sub_nb <- full_nb[cell_id]
  idx_map <- rep.int(NA_integer_, length(full_nb))
  idx_map[cell_id] <- seq_along(cell_id)

  map_fn <- function(v) {
    idx_map[v][!is.na(idx_map[v])]
  }

  relabel_nb <- if (ncores > 1L) {
    if (verbose) {
      .log_info(parent, step, paste0("Relabelling neighbours with ", ncores, " processes"), verbose)
    }
    mclapply(sub_nb, map_fn, mc.cores = ncores)
  } else {
    lapply(sub_nb, map_fn)
  }

  attr(relabel_nb, "class") <- "nb"
  attr(relabel_nb, "region.id") <- meta$grid_id
  attr(relabel_nb, "queen") <- identical(topology, "queen")
  attr(relabel_nb, "topology") <- topology

  list(nb = relabel_nb, backend = backend_used)
}

# =============================================================================
# Round 5 Fix: R Fallbacks for Native Spatial Functions
# =============================================================================

#' Grid Neighbourhood Builder (R fallback for rectangular grids)
#' @param nrow Number of rows in grid
#' @param ncol Number of columns in grid
#' @param queen Logical; if TRUE use queen contiguity
#' @return nb object
#' @keywords internal
.grid_nb_r <- function(nrow, ncol, queen = TRUE) {
    n_cells <- nrow * ncol
    nb <- vector("list", n_cells)
    for (r in seq_len(nrow)) {
        for (c in seq_len(ncol)) {
            idx <- (r - 1L) * ncol + c
            neighbors <- integer(0)
            if (r > 1L) neighbors <- c(neighbors, (r - 2L) * ncol + c)
            if (r < nrow) neighbors <- c(neighbors, r * ncol + c)
            if (c > 1L) neighbors <- c(neighbors, (r - 1L) * ncol + c - 1L)
            if (c < ncol) neighbors <- c(neighbors, (r - 1L) * ncol + c + 1L)
            if (queen) {
                if (r > 1L && c > 1L) neighbors <- c(neighbors, (r - 2L) * ncol + c - 1L)
                if (r > 1L && c < ncol) neighbors <- c(neighbors, (r - 2L) * ncol + c + 1L)
                if (r < nrow && c > 1L) neighbors <- c(neighbors, r * ncol + c - 1L)
                if (r < nrow && c < ncol) neighbors <- c(neighbors, r * ncol + c + 1L)
            }
            nb[[idx]] <- as.integer(sort(neighbors))
        }
    }
    attr(nb, "class") <- "nb"
    attr(nb, "region.id") <- as.character(seq_len(n_cells))
    attr(nb, "queen") <- queen
    attr(nb, "topology") <- if (queen) "queen" else "rook"
    nb
}

#' grid_nb_omp - R fallback
#' @keywords internal
grid_nb_omp <- function(nrow, ncol, queen = TRUE) {
    if (.native_all_explicitly_enabled()) {
        return(.Call(`_geneSCOPE_grid_nb_omp`, nrow, ncol, queen))
    }
    if (identical(tolower(Sys.info()[["sysname"]]), "darwin")) {
        warning("grid_nb_omp: Darwin spatial safe default; using R fallback.", call. = FALSE)
    }
    .grid_nb_r(nrow, ncol, queen)
}

#' grid_nb_serial - R fallback
#' @keywords internal
grid_nb_serial <- function(nrow, ncol, queen = TRUE) {
    if (.native_all_explicitly_enabled()) {
        return(.Call(`_geneSCOPE_grid_nb_serial`, nrow, ncol, queen))
    }
    .grid_nb_r(nrow, ncol, queen)
}

#' grid_nb - R fallback
#' @keywords internal
grid_nb <- function(nrow, ncol, queen = TRUE) {
    if (.native_all_explicitly_enabled()) {
        return(.Call(`_geneSCOPE_grid_nb`, nrow, ncol, queen))
    }
    .grid_nb_r(nrow, ncol, queen)
}

#' grid_nb_hex - R fallback
#' @keywords internal
grid_nb_hex <- function(nrow, ncol, oddr = TRUE) {
    if (.native_all_explicitly_enabled()) {
        return(.Call(`_geneSCOPE_grid_nb_hex`, nrow, ncol, oddr))
    }
    .grid_nb_hex_r(nrow, ncol, oddr)
}

#' grid_nb_hex_omp - R fallback
#' @keywords internal
grid_nb_hex_omp <- function(nrow, ncol, oddr = TRUE) {
    if (.native_all_explicitly_enabled()) {
        return(.Call(`_geneSCOPE_grid_nb_hex_omp`, nrow, ncol, oddr))
    }
    .grid_nb_hex_r(nrow, ncol, oddr)
}

#' grid_nb_hexq - R fallback
#' @keywords internal
grid_nb_hexq <- function(nrow, ncol, oddq = TRUE) {
    if (.native_all_explicitly_enabled()) {
        return(.Call(`_geneSCOPE_grid_nb_hexq`, nrow, ncol, oddq))
    }
    .grid_nb_hexq_r(nrow, ncol, oddq)
}

#' grid_nb_hexq_omp - R fallback
#' @keywords internal
grid_nb_hexq_omp <- function(nrow, ncol, oddq = TRUE) {
    if (.native_all_explicitly_enabled()) {
        return(.Call(`_geneSCOPE_grid_nb_hexq_omp`, nrow, ncol, oddq))
    }
    .grid_nb_hexq_r(nrow, ncol, oddq)
}

#' Hex grid R implementation (odd-r/even-r)
#' @keywords internal
.grid_nb_hex_r <- function(nrow, ncol, oddr = TRUE) {
    n_cells <- nrow * ncol
    nb <- vector("list", n_cells)
    for (r in seq_len(nrow)) {
        for (c in seq_len(ncol)) {
            idx <- (r - 1L) * ncol + c
            is_odd_row <- (r %% 2L) == 1L
            neighbors <- integer(0)
            if (oddr) {
                if (is_odd_row) {
                    offsets <- matrix(c(-1L,-1L, -1L,0L, 0L,-1L, 0L,1L, 1L,-1L, 1L,0L), ncol=2, byrow=TRUE)
                } else {
                    offsets <- matrix(c(-1L,0L, -1L,1L, 0L,-1L, 0L,1L, 1L,0L, 1L,1L), ncol=2, byrow=TRUE)
                }
            } else {
                if (!is_odd_row) {
                    offsets <- matrix(c(-1L,-1L, -1L,0L, 0L,-1L, 0L,1L, 1L,-1L, 1L,0L), ncol=2, byrow=TRUE)
                } else {
                    offsets <- matrix(c(-1L,0L, -1L,1L, 0L,-1L, 0L,1L, 1L,0L, 1L,1L), ncol=2, byrow=TRUE)
                }
            }
            for (i in seq_len(nrow(offsets))) {
                nr <- r + offsets[i,1L]; nc <- c + offsets[i,2L]
                if (nr >= 1L && nr <= nrow && nc >= 1L && nc <= ncol) {
                    neighbors <- c(neighbors, (nr - 1L) * ncol + nc)
                }
            }
            nb[[idx]] <- as.integer(sort(neighbors))
        }
    }
    attr(nb, "class") <- "nb"
    attr(nb, "region.id") <- as.character(seq_len(n_cells))
    attr(nb, "topology") <- if (oddr) "hex-oddr" else "hex-evenr"
    nb
}

#' Hex grid R implementation (odd-q/even-q)
#' @keywords internal
.grid_nb_hexq_r <- function(nrow, ncol, oddq = TRUE) {
    n_cells <- nrow * ncol
    nb <- vector("list", n_cells)
    for (r in seq_len(nrow)) {
        for (c in seq_len(ncol)) {
            idx <- (r - 1L) * ncol + c
            is_odd_col <- (c %% 2L) == 1L
            neighbors <- integer(0)
            if (oddq) {
                if (is_odd_col) {
                    offsets <- matrix(c(-1L,0L, -1L,1L, 0L,-1L, 0L,1L, 1L,0L, 1L,1L), ncol=2, byrow=TRUE)
                } else {
                    offsets <- matrix(c(-1L,-1L, -1L,0L, 0L,-1L, 0L,1L, 1L,-1L, 1L,0L), ncol=2, byrow=TRUE)
                }
            } else {
                if (!is_odd_col) {
                    offsets <- matrix(c(-1L,0L, -1L,1L, 0L,-1L, 0L,1L, 1L,0L, 1L,1L), ncol=2, byrow=TRUE)
                } else {
                    offsets <- matrix(c(-1L,-1L, -1L,0L, 0L,-1L, 0L,1L, 1L,-1L, 1L,0L), ncol=2, byrow=TRUE)
                }
            }
            for (i in seq_len(nrow(offsets))) {
                nr <- r + offsets[i,1L]; nc <- c + offsets[i,2L]
                if (nr >= 1L && nr <= nrow && nc >= 1L && nc <= ncol) {
                    neighbors <- c(neighbors, (nr - 1L) * ncol + nc)
                }
            }
            nb[[idx]] <- as.integer(sort(neighbors))
        }
    }
    attr(nb, "class") <- "nb"
    attr(nb, "region.id") <- as.character(seq_len(n_cells))
    attr(nb, "topology") <- if (oddq) "hex-oddq" else "hex-evenq"
    nb
}

#' listw_B_omp - R fallback
#' @keywords internal
listw_B_omp <- function(nb) {
    if (.native_all_explicitly_enabled()) {
        return(.Call(`_geneSCOPE_listw_B_omp`, nb))
    }
    if (identical(tolower(Sys.info()[["sysname"]]), "darwin")) {
        warning("listw_B_omp: Darwin spatial safe default; using R fallback.", call. = FALSE)
    }
    .listw_b_r(nb)
}

#' Listw B R implementation
#' @keywords internal
.listw_b_r <- function(nb) {
    n <- length(nb)
    if (n == 0L) return(Matrix::sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(0L,0L)))
    lens <- lengths(nb)
    i <- rep.int(seq_len(n), lens)
    j <- unlist(nb, use.names=FALSE)
    Matrix::sparseMatrix(i=i, j=j, x=rep.int(1.0, length(j)), dims=c(n,n))
}

#' nb2mat - R implementation
#' @keywords internal
nb2mat <- function(nb) {
    if (.native_all_explicitly_enabled()) {
        return(.Call(`_geneSCOPE_nb2mat`, nb))
    }
    as.matrix(.listw_b_r(nb))
}

#' grid_weights_kernel_rect_omp - R fallback
#' @keywords internal
grid_weights_kernel_rect_omp <- function(nrow, ncol, gx, gy, queen=TRUE, radius=2L, kernel="gaussian", sigma=1.0) {
    if (.native_all_explicitly_enabled()) {
        return(.Call(`_geneSCOPE_grid_weights_kernel_rect_omp`, nrow, ncol, gx, gy, queen, radius, kernel, sigma))
    }
    if (identical(tolower(Sys.info()[["sysname"]]), "darwin")) {
        warning("grid_weights_kernel_rect_omp: Darwin spatial safe default; using R fallback.", call. = FALSE)
    }
    .grid_weights_kernel_rect_r(nrow, ncol, gx, gy, queen, radius, kernel, sigma)
}

#' grid_weights_kernel_hexr_omp - R fallback
#' @keywords internal
grid_weights_kernel_hexr_omp <- function(nrow, ncol, gx, gy, oddr=TRUE, radius=2L, kernel="gaussian", sigma=1.0) {
    if (.native_all_explicitly_enabled()) {
        return(.Call(`_geneSCOPE_grid_weights_kernel_hexr_omp`, nrow, ncol, gx, gy, oddr, radius, kernel, sigma))
    }
    if (identical(tolower(Sys.info()[["sysname"]]), "darwin")) {
        warning("grid_weights_kernel_hexr_omp: Darwin spatial safe default; using R fallback.", call. = FALSE)
    }
    .grid_weights_kernel_hexr_r(nrow, ncol, gx, gy, oddr, radius, kernel, sigma)
}

#' grid_weights_kernel_hexq_omp - R fallback
#' @keywords internal
grid_weights_kernel_hexq_omp <- function(nrow, ncol, gx, gy, oddq=TRUE, radius=2L, kernel="gaussian", sigma=1.0) {
    if (.native_all_explicitly_enabled()) {
        return(.Call(`_geneSCOPE_grid_weights_kernel_hexq_omp`, nrow, ncol, gx, gy, oddq, radius, kernel, sigma))
    }
    if (identical(tolower(Sys.info()[["sysname"]]), "darwin")) {
        warning("grid_weights_kernel_hexq_omp: Darwin spatial safe default; using R fallback.", call. = FALSE)
    }
    .grid_weights_kernel_hexq_r(nrow, ncol, gx, gy, oddq, radius, kernel, sigma)
}
