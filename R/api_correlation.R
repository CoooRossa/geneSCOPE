#' Compute gene-gene correlations.
#' @description
#' Computes a gene × gene correlation matrix from either single-cell expression
#' (`level = "cell"`) or grid-level expression (`level = "grid"`), stores results
#' into the `scope_object`, and returns the updated object invisibly. The
#' currently supported correlation family is Pearson only; callers who need a
#' different correlation family should compute it externally rather than assume
#' a built-in fallback.
#' @param scope_obj A \code{scope_object} containing expression layers.
#' @param level Correlation level (`cell` or `grid`).
#' @param grid_name Character; name of the grid layer to operate on. If \code{NULL}
#'   and only one grid layer exists, it is auto-selected. If multiple layers exist,
#'   the name must be specified.
#' @param layer Character; name of the expression layer to use (e.g., \code{"logCPM"}, \code{"Xz"}).
#' @param method Requested correlation method. The only supported public value
#'   is `pearson`. Any other value errors explicitly instead of silently
#'   reusing the Pearson backend.
#' @param backend Deprecated compatibility argument. The runtime now uses the
#'   package policy `C++ first, R fallback`; `python` still errors explicitly
#'   because no Python correlation backend is shipped. Inputs must provide at
#'   least two observations
#'   (`n_cells >= 2`); non-finite values error explicitly, and extremely large
#'   finite values (`> 1e300`) trigger a precision warning.
#' @param blocksize Block size passed to the C++ backend (tuning knob).
#' @param ncores Integer; number of threads to use for parallel computation. Default varies by function. Set to 1 for serial execution.
#' @param chunk_size Chunk size used when using file-backed computation.
#' @param memory_limit_gb Numeric; memory threshold in gigabytes. When matrix size exceeds this threshold, chunked or file-backed computation is triggered. Default is 16 GB.
#' @param use_bigmemory Logical; whether to use file-backed matrices for large computations. Default is \code{FALSE}. When \code{TRUE}, temporary files are created in \code{backing_path}.
#' @param backing_path Character; directory for file-backed matrix artifacts. Default is \code{tempdir()}.
#' @param force_compute Whether to force recomputation even if cached outputs exist.
#' @param store_layer Slot name for storing the correlation matrix. When
#' `level = "cell"` and `store_layer` is not supplied, `"_cell"` is appended.
#' @param compute_fdr Whether to compute and store an FDR-adjusted p-value matrix.
#' @param fdr_method Multiple-testing adjustment method passed to `p.adjust()`.
#' @param verbose Logical; whether to emit compact public progress messages.
#'   Default is \code{TRUE}. Set to \code{FALSE} for silent operation, or use
#'   \code{options(geneSCOPE.verbose = FALSE)} for global control.
#' @return The updated `scope_object` (invisibly).
#' @examples
#' \dontrun{
#' scope_obj <- createSCOPE(data_dir = "/path/to/output_dir",　grid_length = 30)
#' scope_obj <- addSingleCells(scope_obj, data_dir = "/path/to/output_dir")
#' scope_obj <- normalizeSingleCells(scope_obj, input_layer = "counts", output_layer = "logCPM", scale_factor = 1e4)
#' scope_obj <- computeCorrelation(scope_obj, level = "cell", layer = "logCPM", method = "pearson", ncores = 16) 
#' scope_obj <- computeCorrelation(scope_obj, level = "grid", grid_name = "grid30", layer = "Xz", method = "pearson", ncores = 16)
#' }
#' @seealso `addSingleCells()`, `normalizeSingleCells()`, `computeL()`
#' @export
computeCorrelation <- function(
    scope_obj,
    level = c("cell", "grid"),
    grid_name = NULL,
    layer = "logCPM",
    method = "pearson",
    backend = "cpp",
    blocksize = 2000,
    ncores = 16,
    chunk_size = 1000,
    memory_limit_gb = 16,
    use_bigmemory = TRUE,
    backing_path = tempdir(),
    force_compute = FALSE,
    store_layer = "pearson_cor",
    compute_fdr = TRUE,
    fdr_method = "BH",
    verbose = TRUE) {
    level <- match.arg(level)
    backend <- if (missing(backend)) {
        "cpp"
    } else {
        .normalize_core_backend(backend, arg = "backend", allow_auto = TRUE)
    }
    verbose <- .resolve_runtime_verbose(verbose)
    if (missing(store_layer) && level == "cell") {
        store_layer <- paste0(store_layer, "_cell")
    }
    .with_log_session("computeCorrelation", verbose = verbose, {
        .log_start("computeCorrelation", verbose = verbose)
        
        .log_key_values("computeCorrelation", list(
            level = level,
            grid = if (!is.null(grid_name)) as.character(grid_name)[1] else NA_character_,
            layer = layer,
            method = method,
            backend = backend,
            store_layer = store_layer
        ))

        result <- .compute_correlation(
            scope_obj = scope_obj,
            level = level,
            grid_name = grid_name,
            layer = layer,
            method = method,
            backend = backend,
            blocksize = blocksize,
            ncores = ncores,
            chunk_size = chunk_size,
            memory_limit_gb = memory_limit_gb,
            use_bigmemory = use_bigmemory,
            backing_path = backing_path,
            force_compute = force_compute,
            store_layer = store_layer,
            compute_fdr = compute_fdr,
            fdr_method = fdr_method,
            verbose = verbose
        )

        cor_holder <- if (identical(level, "cell")) {
            result@stats[["cell"]]
        } else if (!is.null(grid_name)) {
            result@stats[[grid_name]]
        } else {
            NULL
        }
        cor_mat <- if (!is.null(cor_holder)) cor_holder[[store_layer]] else NULL
        fdr_mat <- if (!is.null(cor_holder)) cor_holder[[paste0(store_layer, "_FDR")]] else NULL

        .log_key_values("computeCorrelation", list(
            level = level,
            store_layer = store_layer,
            dims = if (!is.null(cor_mat)) paste(dim(cor_mat), collapse = "x") else NA_character_,
            finite_pairs = if (!is.null(cor_mat)) sum(is.finite(as.numeric(cor_mat))) else NA_integer_,
            sig_fdr_005 = if (!is.null(fdr_mat)) sum(fdr_mat < 0.05, na.rm = TRUE) else NA_integer_
        ))

        .log_done("computeCorrelation", verbose = verbose)
        invisible(result)
    })
}
