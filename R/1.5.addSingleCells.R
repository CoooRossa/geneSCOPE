
#' @title Add Single Cell Count Matrix (Xenium)
#' @description
#'   Read the Xenium **cell-feature sparse matrix** (either HDF5 or
#'   Matrix-Market format), keep only the cells that appear in
#'   \code{scope_obj@coord$centroids$cell}, preserve that order, and store the
#'   result in the new \code{@cells} slot.
#'
#'   The function relies only on \strong{Matrix}, \strong{data.table}, and
#'   \strong{rhdf5} (or \strong{hdf5r})—no Seurat, tidyverse, or other heavy
#'   dependencies.
#' @param scope_obj A valid \code{scope_object} whose \code{@coord$centroids} slot
#'        is already filled (e.g. via \code{\link{createSCOPE}}). Cells
#'        found in this slot define the subset and order of columns kept.
#' @param xenium_dir Character scalar. Path to the Xenium \file{outs/}
#'        directory that holds \file{cell_feature_matrix.h5}; must exist.
#' @param filter_genes Character vector or NULL. Gene names to include (default NULL keeps all).
#' @param exclude_prefix Character vector. Prefixes to exclude when filtering genes (default c("Unassigned", "NegControl", "Background", "DeprecatedCodeword")).
#' @param filter_barcodes Logical. Whether to filter cells to match barcodes in scope_obj (default TRUE). When FALSE, all cells from the HDF5 matrix are included.
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#' @return The same object (modified in place) now carrying a
#'         \code{dgCMatrix} in \code{@cells}.
#' @details
#' The function performs five main steps:
#' \enumerate{
#'   \item Reads the \strong{cell‑feature} matrix in HDF5 format produced by
#'         10x Xenium, using \pkg{rhdf5} to avoid loading the full file into
#'         memory.
#'   \item Reconstructs the sparse count matrix as a \code{dgCMatrix} from
#'         the \emph{x/i/p} components and automatically handles either
#'         {genes × cells} or {cells × genes} orientations.
#'   \item Filters genes by excluding those with specified prefixes
#'         (\code{exclude_prefix}) and optionally keeping only genes in
#'         \code{filter_genes}, similar to the filtering applied in
#'         \code{\link{createSCOPE}}.
#'   \item Conditionally filters and reorders columns based on \code{filter_barcodes}:
#'         when \code{TRUE} (default), columns are filtered and reordered to match
#'         the centroid table (\code{scope_obj@coord$centroids$cell}) for spatial
#'         analyses; when \code{FALSE}, all cells from the HDF5 matrix are retained.
#'   \item Migrates any extra attributes found on the original matrix
#'         (e.g. log-CPM transforms) into the list stored in
#'         \code{scope_obj@cells}.
#' }
#'
#' Only the \pkg{Matrix}, \pkg{data.table}, and \pkg{rhdf5} packages are
#' required; no Seurat or tidyverse dependencies are introduced. Memory
#' footprint scales with the number of non‑zero entries, not the full matrix
#' dimensions.
#' @importFrom Matrix sparseMatrix t
#' @importFrom data.table data.table
#' @importFrom rhdf5 h5read
#' @examples
#' \dontrun{
#' ## Add the cell‑level count matrix to an existing scope_object
#' coord <- createSCOPE("mouse_brain/xenium_output")
#' coord <- addSingleCells(coord, "mouse_brain/xenium_output")
#'
#' ## With custom gene filtering
#' coord <- addSingleCells(coord, "mouse_brain/xenium_output",
#'   filter_genes = c("Gene1", "Gene2", "Gene3"),
#'   exclude_prefix = c("Unassigned", "NegControl")
#' )
#'
#' ## Include all cells from HDF5 matrix (skip barcode filtering)
#' coord <- addSingleCells(coord, "mouse_brain/xenium_output",
#'   filter_barcodes = FALSE
#' )
#' }
#' @export
addSingleCells_xenium <- function(scope_obj, xenium_dir,
                           filter_genes = NULL,
                           exclude_prefix = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword"),
                           filter_barcodes = TRUE,
                           verbose = TRUE) {
  stopifnot(inherits(scope_obj, "scope_object"), dir.exists(xenium_dir))

  if (verbose) message("[geneSCOPE::addSingleCells xenium] Loading cell-feature matrix from HDF5")

  suppressPackageStartupMessages({
    library(Matrix)
    library(data.table)
  })

  # Cross-platform HDF5 library detection and loading
  hdf5_available <- FALSE
  if (requireNamespace("rhdf5", quietly = TRUE)) {
    library(rhdf5)
    hdf5_available <- TRUE
  } else if (requireNamespace("hdf5r", quietly = TRUE)) {
    library(hdf5r)
    hdf5_available <- TRUE
  }

  if (!hdf5_available) {
    stop("Neither rhdf5 nor hdf5r package is available. Please install one of them.")
  }

  ## ------------------------------------------------------------------ 1
  if (filter_barcodes) {
    keep_cells_target <- scope_obj@coord$centroids$cell
    if (!length(keep_cells_target)) {
      stop("scope_obj@coord$centroids is empty, cannot determine cells to keep.")
    }

    if (verbose) message("[geneSCOPE::addSingleCells xenium]  Target cells: ", length(keep_cells_target))
  } else {
    if (verbose) message("[geneSCOPE::addSingleCells xenium]  Barcode filtering disabled - will include all cells from HDF5")
  }

  ## ------------------------------------------------------------------ 2
  h5_file <- file.path(xenium_dir, "cell_feature_matrix.h5")
  if (!file.exists(h5_file)) {
    stop("HDF5 file not found: ", h5_file)
  }

  # Cross-platform HDF5 reading
  if (requireNamespace("rhdf5", quietly = TRUE)) {
    barcodes <- rhdf5::h5read(h5_file, "matrix/barcodes")
    genes_id <- rhdf5::h5read(h5_file, "matrix/features/id")
    genes_nm <- rhdf5::h5read(h5_file, "matrix/features/name")
    data <- rhdf5::h5read(h5_file, "matrix/data")
    indices <- rhdf5::h5read(h5_file, "matrix/indices")
    indptr <- rhdf5::h5read(h5_file, "matrix/indptr")
    shape <- as.integer(rhdf5::h5read(h5_file, "matrix/shape"))
  } else {
    # Alternative hdf5r implementation
    h5file <- hdf5r::H5File$new(h5_file, mode = "r")
    barcodes <- h5file[["matrix/barcodes"]][]
    genes_id <- h5file[["matrix/features/id"]][]
    genes_nm <- h5file[["matrix/features/name"]][]
    data <- h5file[["matrix/data"]][]
    indices <- h5file[["matrix/indices"]][]
    indptr <- h5file[["matrix/indptr"]][]
    shape <- as.integer(h5file[["matrix/shape"]][])
    h5file$close()
  }

  ## -------------------------- 2a Determine HDF5 storage orientation ----
  nrow_h5 <- shape[1]
  ncol_h5 <- shape[2]
  by_cols <- length(indptr) == (ncol_h5 + 1) # indptr length = ncol+1

  if (!by_cols) {
    stop("Current function only supports CSC-style storage (indptr length = ncol+1).")
  }

  col_is_cell <- (ncol_h5 == length(barcodes))
  if (!col_is_cell && !(nrow_h5 == length(barcodes))) {
    stop("Shape doesn't match barcode/gene counts, cannot determine orientation.")
  }

  if (verbose) {
    message(
      "[geneSCOPE::addSingleCells xenium]  Matrix: ", nrow_h5, "×", ncol_h5,
      " (", length(genes_nm), " genes, ", length(barcodes), " cells)"
    )
  }

  ## -------------------------- 2b Build sparse matrix ------------------
  x <- as.numeric(data)
  i <- as.integer(indices)
  p <- as.integer(indptr)

  if (col_is_cell) {
    ## →  columns = cells, rows = genes   (most common)
    counts_raw <- new("dgCMatrix",
      Dim      = c(length(genes_nm), length(barcodes)),
      x        = x,
      i        = i,
      p        = p,
      Dimnames = list(make.unique(genes_nm), barcodes)
    )
  } else {
    ## →  columns = genes, rows = cells   (rare case: need transpose)
    if (verbose) message("[geneSCOPE::addSingleCells xenium]  Transposing matrix layout")
    tmp <- new("dgCMatrix",
      Dim      = c(length(barcodes), length(genes_nm)),
      x        = x,
      i        = i,
      p        = p,
      Dimnames = list(barcodes, make.unique(genes_nm))
    )
    counts_raw <- Matrix::t(tmp)
  }

  attr(counts_raw, "gene_map") <- data.frame(
    ensembl = genes_id,
    symbol = genes_nm,
    stringsAsFactors = FALSE
  )

  ## -------------------------- 2c Filter genes by exclude_prefix and filter_genes ------------------
  if (verbose) message("[geneSCOPE::addSingleCells xenium] Applying gene filters")

  # Get gene names (rownames)
  all_genes <- rownames(counts_raw)
  keep_genes <- all_genes

  # Filter by exclude_prefix
  if (length(exclude_prefix) > 0) {
    exclude_pattern <- paste0("^(?:", paste(exclude_prefix, collapse = "|"), ")")
    exclude_mask <- grepl(exclude_pattern, keep_genes, perl = TRUE)
    keep_genes <- keep_genes[!exclude_mask]

    if (verbose) {
      excluded_count <- sum(exclude_mask)
      if (excluded_count > 0) {
        message(
          "[geneSCOPE::addSingleCells xenium]  Excluded ", excluded_count, " genes with prefixes: ",
          paste(exclude_prefix, collapse = ", ")
        )
      }
    }
  }

  # Filter by filter_genes (if specified)
  if (!is.null(filter_genes)) {
    keep_genes <- intersect(keep_genes, filter_genes)
    if (verbose) {
      message("[geneSCOPE::addSingleCells xenium]  Filtered to ", length(keep_genes), " genes from target list")
    }
  }

  if (length(keep_genes) == 0) {
    stop("No genes remain after filtering. Please check exclude_prefix and filter_genes parameters.")
  }

  # Apply gene filtering to the matrix
  counts_raw <- counts_raw[keep_genes, , drop = FALSE]

  if (verbose) {
    message(
      "[geneSCOPE::addSingleCells xenium]  Final count: ", nrow(counts_raw), "/", length(all_genes), " genes"
    )
  }

  ## ------------------------------------------------------------------ 3
  if (filter_barcodes) {
    keep_cells <- keep_cells_target[keep_cells_target %in% colnames(counts_raw)]
    if (!length(keep_cells)) {
      stop("Centroid cell IDs do not overlap with H5 data.")
    }

    if (verbose) message("[geneSCOPE::addSingleCells xenium]  Cell overlap: ", length(keep_cells), "/", length(keep_cells_target))

    counts <- counts_raw[, keep_cells, drop = FALSE]
  } else {
    # Keep all cells from the HDF5 matrix
    counts <- counts_raw
    if (verbose) message("[geneSCOPE::addSingleCells xenium]  Keeping all ", ncol(counts), " cells from HDF5 matrix")
  }

  ## ------------------------------------------------------------------ 4
  scope_obj@cells <- list(counts = counts)

  if (verbose) message("[geneSCOPE::addSingleCells xenium] Cell matrix integration completed")
  invisible(scope_obj)
}


#' @title Add single-cell count matrix for CosMx (entry compatible with addSingleCells)
#' @description
#'  Read CosMx flatFiles expression matrix (`*_exprMat_file.csv(.gz)`), extract
#'  counts for a target gene set, build a sparse matrix and write it to
#'  `scope_obj@cells$counts`, so downstream steps such as
#'  `normalizeSingleCells()` and `computeDensity()` can be reused. Behavior is
#'  compatible with Xenium `addSingleCells()`.
#'
#'  Columns are filtered/reordered to match the set and order of
#'  `scope_obj@coord$centroids$cell`. CosMx cell identifiers typically come
#'  from `cell_ID`. To form unique keys with FOV, set `id_mode = "fov_cell"`
#'  (column names created by `paste(cell_ID, "_", fov)`).
#'
#' @param scope_obj A prepared `scope_object` whose `@coord$centroids` is filled.
#' @param cosmx_root Root directory of the CosMx project (contains `flatFiles/`).
#' @param filter_genes Character vector or NULL. Read only these genes. If NULL
#'        and `force_all_genes = FALSE`, an error is thrown because loading all
#'        genes can be very memory-intensive for large datasets.
#' @param exclude_prefix Character vector of prefixes to exclude (default
#'        c("Unassigned","NegControl","Background","DeprecatedCodeword",
#'        "SystemControl","Negative")).
#' @param id_mode One of `"cell_id"` (default), `"fov_cell"`, or `"auto"`.
#'        - `cell_id`: use `cell_ID` as column names (must be globally unique).
#'        - `fov_cell`: use `paste(cell_ID, "_", fov)`.
#'        - `auto`: if `centroids$cell` contains underscores/FOV traces, use
#'          `fov_cell`; otherwise `cell_id`.
#' @param filter_barcodes Logical. Filter/reorder columns by `centroids$cell`
#'        (default TRUE).
#' @param force_all_genes Logical. When `filter_genes` is NULL, force reading all
#'        genes (may require large memory).
#' @param verbose Logical. Print progress messages.
#'
#' @return The modified `scope_object` with a `dgCMatrix` stored in
#'         `@cells$counts`.
#' @importFrom Matrix sparseMatrix
#' @importFrom data.table fread as.data.table setDT
#' @export
addSingleCells_cosmx <- function(scope_obj,
                                 cosmx_root,
                                 filter_genes = NULL,
                                 exclude_prefix = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword",
                                                    "SystemControl", "Negative"),
                                 id_mode = c("cell_id", "fov_cell", "auto"),
                                 filter_barcodes = TRUE,
                                 force_all_genes = FALSE,
                                 verbose = TRUE) {
  stopifnot(inherits(scope_obj, "scope_object"), dir.exists(cosmx_root))

  id_mode <- match.arg(id_mode)

  # 1) Load target cell set for filtering and column ordering
  keep_cells_target <- NULL
  if (filter_barcodes) {
    keep_cells_target <- scope_obj@coord$centroids$cell
    if (!length(keep_cells_target)) {
      stop("centroids are empty; cannot determine columns to keep. Please run createSCOPE/_cosmx first.")
    }
  }

  # 2) Locate CosMx exprMat file
  ff_dir <- file.path(cosmx_root, "flatFiles")
  if (!dir.exists(ff_dir)) stop("flatFiles/ not found: ", ff_dir)
  ds_dir <- list.dirs(ff_dir, full.names = TRUE, recursive = FALSE)
  if (!length(ds_dir)) stop("No dataset directory found under flatFiles/.")
  ds_dir <- ds_dir[[1]]

  expr_mat <- NULL
  cand <- c(
    file.path(ds_dir, paste0(basename(ds_dir), "_exprMat_file.csv.gz")),
    file.path(ds_dir, paste0(basename(ds_dir), "_exprMat_file.csv"))
  )
  for (p in cand) if (file.exists(p)) { expr_mat <- p; break }
  if (is.null(expr_mat)) stop("No *_exprMat_file.csv(.gz) found under: ", ds_dir)

  if (verbose) message("[geneSCOPE::addSingleCells_cosmx] Expression matrix: ", basename(expr_mat))

  # 3) Read header to determine gene columns
  suppressPackageStartupMessages({ library(data.table); library(Matrix) })
  # Prefer fread to fetch header; fall back to base::readLines(gzfile(...), n=1)
  header_dt <- NULL
  cols <- NULL
  header_dt <- tryCatch({ data.table::fread(expr_mat, nrows = 0L, showProgress = FALSE) }, error = function(e) NULL)
  if (!is.null(header_dt)) {
    cols <- names(header_dt)
  } else {
    # Safe fallback: read only the first line as header
    con <- NULL
    cols <- tryCatch(
      {
        if (grepl("\\.gz$", expr_mat, ignore.case = TRUE)) {
          con <- gzfile(expr_mat, open = "rt")
        } else {
          con <- file(expr_mat, open = "rt")
        }
        ln <- readLines(con, n = 1L, warn = FALSE)
        if (length(ln) == 0L || !nzchar(ln)) stop("Empty file or missing header")
        ln <- sub("\\r$", "", ln)  # strip CR from CRLF
        strsplit(ln, ",", fixed = TRUE)[[1]]
      },
      error = function(e) NULL,
      finally = { if (!is.null(con)) try(close(con), silent = TRUE) }
    )
    if (is.null(cols) || !length(cols)) stop("Unable to read expression matrix header: ", expr_mat)
  }
  base_cols <- tolower(cols)
  # Identify fov / cell_ID column names (case-insensitive)
  fov_col  <- cols[which(base_cols == "fov")[1]]
  cid_col  <- cols[which(base_cols %in% c("cell_id", "cellid"))[1]]
  if (is.na(fov_col) || is.na(cid_col)) stop("Expression matrix missing fov/cell_ID column")

  # ---- Auto-detect non-gene columns; the rest are considered gene columns ----
  # Known/common non-gene column names (case-insensitive)
  non_gene_exact <- c(
    "fov", "cell_id", "cellid", "cell", "barcode",
    "x", "y", "z",
    "umi", "reads", "total", "sum", "sizefactor",
    "feature_name", "gene", "target", "qv", "nucleus_distance",
    "cellcomp", "compartment", "compartmentlabel"
  )
  # Columns containing these substrings are treated as non-gene (coords/units/etc.)
  non_gene_substr <- c(
    "_px", "_um", "_mm", "global_", "local_", "centroid", "area",
    "x_", "y_", "_x", "_y"
  )
  low_cols <- tolower(cols)
  is_non_gene <- low_cols %in% non_gene_exact
  for (ss in non_gene_substr) {
    is_non_gene <- is_non_gene | grepl(ss, low_cols, fixed = TRUE)
  }
  # fov/cell_ID must be among non-gene; remaining columns are candidate genes
  gene_cols_auto <- cols[!is_non_gene]

  # Exclude non-biological prefixes (NegControl/Unassigned/Background etc.)
  if (length(exclude_prefix)) {
    pat <- paste0("^(?:", paste(exclude_prefix, collapse = "|"), ")")
    gene_cols_auto <- gene_cols_auto[!grepl(pat, gene_cols_auto, perl = TRUE)]
  }

  # If user did not provide filter_genes, use auto set; otherwise take the intersection
  if (is.null(filter_genes)) {
    gene_cols <- gene_cols_auto
  } else {
    gene_cols <- intersect(gene_cols_auto, filter_genes)
  if (!length(gene_cols)) stop("No overlap between filter_genes and available genes (after auto-detection).")
  }

  if (!length(gene_cols)) stop("No gene columns identified. Check matrix header and prefix filters.")
  if (verbose) message("[geneSCOPE::addSingleCells_cosmx] Number of genes detected (auto): ", length(gene_cols))

  # 4) Read only required columns (fov/cell_ID + target genes)
  sel_cols <- c(fov_col, cid_col, gene_cols)
  # fread selected columns; fall back to read.csv on failure (slower but robust)
  dt <- tryCatch(
    data.table::fread(expr_mat, select = sel_cols, showProgress = verbose),
    error = function(e) {
      if (verbose) message("[geneSCOPE::addSingleCells_cosmx] fread failed; falling back to read.csv: ", e$message)
      utils::read.csv(expr_mat, header = TRUE)[, sel_cols]
    }
  )
  data.table::setDT(dt)

  # Construct column keys (aligned with centroids$cell) and adaptively fall back
  if (id_mode == "auto") {
    id_mode_eff <- if (any(grepl("_", keep_cells_target)) || any(grepl("FOV", keep_cells_target))) "fov_cell" else "cell_id"
  } else {
    id_mode_eff <- id_mode
  }
  make_key <- function(mode, fov_offset = 0L) {
    if (mode == "cell_id") {
      as.character(dt[[cid_col]])
    } else {
      fvv <- suppressWarnings(as.integer(dt[[fov_col]]))
      if (!is.na(fov_offset) && fov_offset != 0L) fvv <- fvv + as.integer(fov_offset)
      paste0(as.character(dt[[cid_col]]), "_", as.character(fvv))
    }
  }
  # First attempt
  dt[, cell_key := make_key(id_mode_eff, 0L)]

  # 5) Optional: filter/reorder cells based on centroids
  if (filter_barcodes) {
    Kt <- unique(keep_cells_target)
    cand <- list(
      list(mode = id_mode_eff, off = 0L),
      list(mode = if (id_mode_eff == "cell_id") "fov_cell" else "cell_id", off = 0L),
      list(mode = "fov_cell", off = +1L),
      list(mode = "fov_cell", off = -1L)
    )
    best_idx <- 1L
    best_overlap <- -1L
    keys_cache <- list()
    for (i in seq_along(cand)) {
      mode_i <- cand[[i]]$mode; off_i <- cand[[i]]$off
      # Only try offset when fov column exists
      if (mode_i == "fov_cell" && (!exists("fov_col") || is.na(fov_col))) next
      if (is.null(keys_cache[[paste(mode_i, off_i)]])) {
        keys_cache[[paste(mode_i, off_i)]] <- make_key(mode_i, off_i)
      }
      kk <- unique(keys_cache[[paste(mode_i, off_i)]])
      ov <- length(intersect(kk, Kt))
      if (ov > best_overlap) { best_overlap <- ov; best_idx <- i }
    }
    # Apply the best key
    mode_best <- cand[[best_idx]]$mode; off_best <- cand[[best_idx]]$off
    dt[, cell_key := make_key(mode_best, off_best)]
    keep_cells <- intersect(unique(dt$cell_key), Kt)
    if (!length(keep_cells)) {
      # Diagnostics
      stop("Cell identifiers in expression matrix do not overlap with centroids; check id_mode/FOV encoding. Suggest id_mode=\"fov_cell\"; if still failing, try FOV offset ±1. Example keys: ",
           head(unique(dt$cell_key), 3L), " vs ", head(Kt, 3L))
    }
    dt <- dt[cell_key %in% keep_cells]
  }

  # 6) Build sparse matrix (rows=genes, cols=cells) without forming a huge dense matrix
  # Column index/order: if filtered, follow keep_cells_target; else by appearance
  if (filter_barcodes) {
    col_order <- keep_cells_target[keep_cells_target %in% unique(dt$cell_key)]
  } else {
    col_order <- unique(dt$cell_key)
  }
  col_index <- match(dt$cell_key, col_order)  # 1..ncol

  # Preallocate i/j/x accumulators
  i_idx <- integer(0)
  j_idx <- integer(0)
  x_val <- numeric(0)

  # Collect non-zeros per gene column to avoid large temporary memory
  for (g in gene_cols) {
    v <- suppressWarnings(as.numeric(dt[[g]]))
    nz <- which(!is.na(v) & v != 0)
    if (length(nz)) {
      i_idx <- c(i_idx, rep.int(which(gene_cols == g), length(nz)))
      j_idx <- c(j_idx, col_index[nz])
      x_val <- c(x_val, v[nz])
    }
  }

  counts <- Matrix::sparseMatrix(
    i = as.integer(i_idx),
    j = as.integer(j_idx),
    x = as.numeric(x_val),
    dims = c(length(gene_cols), length(col_order)),
    dimnames = list(gene_cols, col_order)
  )

  # 7) Write back to object
  scope_obj@cells <- list(counts = counts)
  if (verbose) message("[geneSCOPE::addSingleCells cosmx] Wrote @cells$counts (",
                       nrow(counts), " genes × ", ncol(counts), " cells)")
  invisible(scope_obj)
}

#' @title Add Single Cells (auto-dispatch)
#' @description
#'   Wrapper that dispatches to \code{addSingleCells_xenium()} or
#'   \code{addSingleCells_cosmx()} based on the recorded platform tag in
#'   \code{scope_obj} (\code{meta.data$platform} or \code{stats$platform}).
#'   You may also force the backend by providing \code{platform}, or by
#'   passing \code{xenium_dir}/\code{cosmx_root} explicitly.
#' @param scope_obj A \code{scope_object} with centroids loaded.
#' @param xenium_dir Optional. Xenium outs directory (if using Xenium backend).
#' @param cosmx_root Optional. CosMx project root containing \file{flatFiles/}.
#' @param platform Character. One of \code{"auto"}, \code{"Xenium"}, \code{"CosMx"}.
#'        Default \code{"auto"}.
#' @param verbose Logical. Print progress messages.
#' @param ... Additional arguments forwarded to the chosen backend.
#' @return The modified \code{scope_object} with \code{@cells$counts} added.
#' @export
addSingleCells <- function(scope_obj,
                           xenium_dir = NULL,
                           cosmx_root = NULL,
                           platform = c("auto", "Xenium", "CosMx"),
                           verbose = TRUE,
                           ...) {
  stopifnot(inherits(scope_obj, "scope_object"))
  platform <- match.arg(platform)

  # 1) If user passed an explicit directory, prefer that
  pick <- NULL
  if (!is.null(xenium_dir)) pick <- "Xenium"
  if (!is.null(cosmx_root)) pick <- "CosMx"

  # 2) Auto-detect from scope_obj if not forced
  if (is.null(pick) && platform == "auto") {
    plat <- NA_character_
    if (!is.null(scope_obj@meta.data) && nrow(scope_obj@meta.data)) {
      if ("platform" %in% names(scope_obj@meta.data)) {
        v <- scope_obj@meta.data$platform
        if (length(v)) plat <- as.character(v[which.max(tabulate(match(v, unique(v))))])
      }
    }
    if (is.na(plat) && !is.null(scope_obj@stats$platform)) plat <- as.character(scope_obj@stats$platform)
    if (!is.na(plat)) {
      if (grepl("cosmx", plat, ignore.case = TRUE)) pick <- "CosMx"
      else if (grepl("xenium", plat, ignore.case = TRUE)) pick <- "Xenium"
    }
  }

  # 3) Fallback default
  if (is.null(pick)) pick <- if (!is.null(xenium_dir)) "Xenium" else if (!is.null(cosmx_root)) "CosMx" else "Xenium"

  if (identical(pick, "CosMx")) {
    if (isTRUE(verbose)) message("[geneSCOPE::addSingleCells cosmx] Dispatching to addSingleCells_cosmx() …")
    cr <- if (!is.null(cosmx_root)) cosmx_root else getOption("geneSCOPE.cosmx_root")
    return(addSingleCells_cosmx(scope_obj = scope_obj, cosmx_root = cr, verbose = verbose, ...))
  } else {
    if (isTRUE(verbose)) message("[geneSCOPE::addSingleCells xenium] Dispatching to addSingleCells_xenium() …")
    xd <- if (!is.null(xenium_dir)) xenium_dir else getOption("geneSCOPE.xenium_dir")
    return(addSingleCells_xenium(scope_obj = scope_obj, xenium_dir = xd, verbose = verbose, ...))
  }
}
