#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
is_cli <- any(grepl("^--file=", commandArgs(trailingOnly = FALSE)))

parse_args <- function(args) {
  out <- list()
  if (!length(args)) return(out)
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    item <- sub("^--", "", arg)
    kv <- strsplit(item, "=", fixed = TRUE)[[1]]
    key <- kv[1]
    value <- if (length(kv) > 1L) kv[2] else "TRUE"
    out[[key]] <- value
  }
  out
}

`%||%` <- function(x, y) if (!is.null(x)) x else y

flag_true <- function(x) {
  if (is.null(x)) return(FALSE)
  tolower(as.character(x)) %in% c("1", "true", "yes", "y", "on")
}

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript repro_sigpipe_createSCOPE.R --xenium_dir=/path/to/outs \\\n",
    "    --grid_length=30 --ncores=8 --chunk_size=5e5 \\\n",
    "    --parallel_backend=auto --batch_size=16 --roi_path=/path/to/roi.csv \\\n",
    "    --self_check=true --debug_parallel=true --debug_parallel_dir=/tmp\n\n",
    "Notes:\n",
    "  - Use --key=value syntax.\n",
    "  - Only runs the ROI clip stages (centroids + prefetch molecules).\n",
    "  - Use --self_check=true for synthetic serial/fork/psock parity.\n",
    sep = ""
  )
  invisible(1L)
}

opts <- parse_args(args)

xenium_dir <- opts$xenium_dir %||% opts$data_dir
if (is.null(xenium_dir)) {
  status <- usage()
  if (is_cli) quit(status = status)
} else {
  grid_length <- as.numeric(opts$grid_length %||% 30)
  if (!is.finite(grid_length) || grid_length <= 0) stop("grid_length must be positive numeric")

  ncores <- suppressWarnings(as.integer(opts$ncores %||% parallel::detectCores()))
  if (!is.finite(ncores) || ncores < 1L) ncores <- 1L

  chunk_size <- as.numeric(opts$chunk_size %||% 5e5)
  if (!is.finite(chunk_size) || chunk_size < 1) stop("chunk_size must be positive numeric")

  batch_size <- suppressWarnings(as.integer(opts$batch_size %||% opts$batch_chunks %||% NA))
  if (!is.na(batch_size)) {
    if (!is.finite(batch_size) || batch_size < 1L) stop("batch_size must be positive integer")
    options(geneSCOPE.psock_batch_chunks = batch_size)
  }

  parallel_backend <- tolower(as.character(opts$parallel_backend %||% "auto"))
  valid_backends <- c("auto", "fork", "psock", "serial")
  if (!(parallel_backend %in% valid_backends)) {
    stop("parallel_backend must be one of: ", paste(valid_backends, collapse = ", "))
  }

  seg_type <- tolower(as.character(opts$seg_type %||% "cell"))
  seg_type <- match.arg(seg_type, c("cell", "nucleus", "both", "none"))

  roi_path <- opts$roi_path %||% NULL
  self_check <- flag_true(opts$self_check)

  debug_parallel <- flag_true(opts$debug_parallel)
  if (debug_parallel) {
    options(geneSCOPE.debug_parallel = TRUE)
    if (!is.null(opts$debug_parallel_dir)) {
      options(geneSCOPE.debug_parallel_dir = opts$debug_parallel_dir)
    }
  }

  if (!requireNamespace("geneSCOPE", quietly = TRUE)) {
    if (requireNamespace("devtools", quietly = TRUE)) {
      devtools::load_all(".")
    } else {
      stop("geneSCOPE is not installed. Install the package or run from repo with devtools available.")
    }
  } else {
    suppressPackageStartupMessages(library(geneSCOPE))
  }

  set.seed(1)

  message("[repro] xenium_dir=", xenium_dir)
  message("[repro] grid_length=", grid_length, " ncores=", ncores, " chunk_size=", chunk_size)
  message("[repro] parallel_backend=", parallel_backend, " seg_type=", seg_type, " self_check=", self_check)
  if (!is.na(batch_size)) message("[repro] batch_size=", batch_size)
  if (!is.null(roi_path)) message("[repro] roi_path=", roi_path)

  run_self_check <- function() {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("self_check requires the sf package.")
    }
    message("[repro] self_check: comparing serial/fork/psock ROI clip outputs")
    pts <- data.frame(x = runif(200), y = runif(200))
    poly <- sf::st_polygon(list(matrix(
      c(0.25, 0.25, 0.75, 0.25, 0.75, 0.75, 0.25, 0.75, 0.25, 0.25),
      ncol = 2, byrow = TRUE
    )))
    region <- sf::st_sfc(poly)
    backends <- c("serial", "psock")
    if (.Platform$OS.type != "windows") backends <- c(backends, "fork")
    self_batch <- if (!is.na(batch_size)) batch_size else 4L
    res_list <- lapply(backends, function(backend) {
      geneSCOPE:::.clip_points_to_region(
        pts,
        region,
        chunk_size = 50,
        workers = 2L,
        parallel_backend = backend,
        psock_batch_chunks = self_batch
      )
    })
    names(res_list) <- backends
    normalize_xy <- function(dt) {
      if (nrow(dt) == 0L) return(dt[, c("x", "y"), drop = FALSE])
      dt <- dt[, c("x", "y"), drop = FALSE]
      dt[order(dt$x, dt$y), , drop = FALSE]
    }
    ref <- normalize_xy(res_list[[1]])
    for (backend in backends[-1]) {
      comp <- normalize_xy(res_list[[backend]])
      if (!isTRUE(all.equal(ref, comp, tolerance = 1e-8))) {
        stop("self_check failed: backend mismatch for ", backend)
      }
    }
    message("[repro] self_check ok: ", paste(backends, collapse = ", "))
  }

  run_clip <- function() {
    state <- geneSCOPE:::.initialize_dataset_builder_state(
      input_dir = xenium_dir,
      grid_length = grid_length,
      seg_type = seg_type,
      filter_genes = NULL,
      max_dist_mol_nuc = 25,
      filtermolecule = TRUE,
      exclude_prefix = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword", "SystemControl", "Negative"),
      filterqv = 20,
      coord_file = roi_path,
      ncores = ncores,
      chunk_pts = chunk_size,
      parallel_backend = parallel_backend,
      min_gene_types = 1,
      max_gene_types = Inf,
      min_seg_points = 1,
      keep_partial_grid = FALSE,
      verbose = TRUE,
      flip_y = TRUE,
      data_type = "xenium"
    )

    state <- geneSCOPE:::.configure_worker_threads(state)
    if (!is.null(state$thread_context$restore_fn)) {
      on.exit(state$thread_context$restore_fn(), add = TRUE)
    }
    state <- geneSCOPE:::.limit_thread_environment(state)
    if (!is.null(state$env_restore)) on.exit(state$env_restore(), add = TRUE)

    state <- geneSCOPE:::.load_dataset_dependencies(state)
    state <- geneSCOPE:::.resolve_input_paths(state)
    state <- geneSCOPE:::.load_centroid_table(state)
    state <- geneSCOPE:::.open_transcript_dataset(state)
    state <- geneSCOPE:::.apply_dataset_filters(state)
    state <- geneSCOPE:::.prepare_roi_geometry(state)

    message("[repro] clipping centroids")
    state <- geneSCOPE:::.clip_points_within_roi(state)

    message("[repro] prefetch + clip ROI molecules")
    state <- geneSCOPE:::.prefetch_roi_molecules(state)

    cells_kept <- if (!is.null(state$objects$centroids)) nrow(state$objects$centroids) else 0L
    mol_kept <- if (!is.null(state$objects$mol_small)) nrow(state$objects$mol_small) else 0L
    message("[repro] done: cells_retained=", cells_kept, " mol_retained=", mol_kept)
  }

  if (self_check) run_self_check()
  geneSCOPE:::.with_data_table_attached(run_clip)
}
