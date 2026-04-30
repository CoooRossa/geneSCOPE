#' Compute spatial weights for a grid layer (public wrapper).
#' @description
#' Builds spatial neighbourhood structures and stores adjacency/weights into the
#' specified grid layer of a `scope_object`.
#' @details
#' On macOS/Darwin, the native OpenMP grid builder is disabled by default to avoid
#' thread-safety issues. The R serial fallback is used instead, which is slower
#' but more stable. Performance overhead: the serial builder is typically 10-15%
#' slower than the parallel native builder. To force-enable the native OpenMP
#' builder on Darwin (not recommended), set
#' `GENESCOPE_ENABLE_OMP_GRID_NB=force` in the environment before loading the
#' package.
#' @param scope_obj A \code{scope_object} with grid layers.
#' @param grid_name Character; name of the grid layer to operate on. If \code{NULL}
#'   and only one grid layer exists, it is auto-selected. If multiple layers exist,
#'   the name must be specified.
#' @param style Style for `spdep::nb2listw()` (default `"B"`).
#'   Kernel weights (main-2 parity) are enabled by setting `style = "kernel_gaussian"`
#'   or `style = "kernel_flat"`. In that case, use `kernel_radius`/`kernel_sigma`
#'   to control the kernel shape. For backward compatibility, `max_order` and
#'   `min_neighbors` are only used as fallbacks when the kernel-specific
#'   parameters are not supplied.
#' @param store_mat Store the adjacency matrix in the grid layer when TRUE.
#' @param zero_policy Whether zero-neighbour cells are allowed.
#' @param store_listw Store the `listw` object when TRUE.
#' @param verbose Logical; whether to emit compact public progress messages.
#'   Default is \code{TRUE}. Set to \code{FALSE} for silent operation, or use
#'   \code{options(geneSCOPE.verbose = FALSE)} for global control.
#' @param topology Topology choice (`auto`, `queen`, `rook`, `hex`, `fuzzy_queen`, `fuzzy_hex`).
#' @param nb_order Optional higher-order neighbour expansion.
#' @param min_neighbors Minimum neighbour count when expanding (non-kernel styles).
#' @param max_order Maximum neighbour order (non-kernel styles).
#' @param kernel_radius Kernel radius for `style="kernel_*"` (overrides `max_order`).
#' @param kernel_sigma Kernel sigma for `style="kernel_*"` (overrides `min_neighbors`).
#' @param ncores Number of threads to use.
#' @param backend Deprecated compatibility argument. The runtime now uses the
#'   package policy `C++ first, R fallback`; `python` still errors explicitly
#'   because no Python spatial-weights backend is shipped.
#' @return The modified `scope_object` (invisibly).
#' @examples
#' \dontrun{
#' scope_obj <- computeWeights(scope_obj, grid_name = "grid30", topology = "queen")
#' }
#' @seealso `computeL()`, `fuzzy_queen_jaccard()`
#' @export
computeWeights <- function(
    scope_obj,
    grid_name = NULL,
    style = "B",
    store_mat = TRUE,
    zero_policy = TRUE,
    store_listw = TRUE,
    verbose = TRUE,
    topology = c("auto", "queen", "rook", "hex", "fuzzy_queen", "fuzzy_hex"),
    nb_order = 1L,
    min_neighbors = NULL,
    max_order = NULL,
    kernel_radius = NULL,
    kernel_sigma = NULL,
    ncores = detectCores(logical = TRUE),
    backend = "cpp") {
    backend <- .normalize_core_backend(backend, arg = "backend", allow_auto = TRUE)
    verbose <- .resolve_runtime_verbose(verbose)
    .with_log_session("computeWeights", verbose = verbose, {
        .log_start("computeWeights", verbose = verbose)
        
        resolved_grid_name <- if (!is.null(grid_name) && nzchar(grid_name)) {
            as.character(grid_name)[1]
        } else if (length(scope_obj@grid) == 1L) {
            names(scope_obj@grid)[1]
        } else {
            NULL
        }

        .log_key_values("computeWeights", list(
            grid = resolved_grid_name,
            style = style,
            topology = as.character(topology)[1],
            backend = backend
        ))

        result <- .compute_weights(
            scope_obj = scope_obj,
            grid_name = grid_name,
            style = style,
            store_mat = store_mat,
            zero_policy = zero_policy,
            store_listw = store_listw,
            verbose = verbose,
            topology = topology,
            nb_order = nb_order,
            min_neighbors = min_neighbors,
            max_order = max_order,
            kernel_radius = kernel_radius,
            kernel_sigma = kernel_sigma,
            ncores = ncores,
            backend = backend
        )

        grid_name_after <- resolved_grid_name
        if (is.null(grid_name_after) && length(result@grid) == 1L) {
            grid_name_after <- names(result@grid)[1]
        }
        grid_layer <- if (!is.null(grid_name_after) && grid_name_after %in% names(result@grid)) result@grid[[grid_name_after]] else NULL
        W <- if (!is.null(grid_layer)) grid_layer$W else NULL
        weights_meta <- if (!is.null(grid_layer) && is.list(grid_layer$weights_meta)) grid_layer$weights_meta else NULL

        .log_key_values("computeWeights", list(
            grid = grid_name_after,
            dims = if (!is.null(W)) paste(dim(W), collapse = "x") else NA_character_,
            nnz = if (!is.null(W)) Matrix::nnzero(W) else NA_integer_,
            style = .log_first_scalar(weights_meta$weight_style, attr(W, "weight_style"), style),
            backend_selected = .log_first_scalar(weights_meta$backend_selected, attr(W, "backend"))
        ))

        .log_done("computeWeights", verbose = verbose)
        invisible(result)
    })
}
