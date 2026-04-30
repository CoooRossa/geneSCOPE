#' Validate core numeric compute matrices.
#' @description
#' Shared preflight for native/R fallback compute paths. It rejects single-cell
#' or single-spot inputs for statistics that require at least two observations,
#' rejects non-finite numeric values, and warns on extremely large finite values
#' that can lose precision in cross-products.
#' @param x Matrix-like object.
#' @param caller Function label used in messages.
#' @param observation_axis Whether observations are stored in rows or columns.
#' @param extreme_threshold Finite magnitude that triggers a precision warning.
#' @return `x`, invisibly.
#' @keywords internal
.validate_core_numeric_input <- function(x,
                                         caller,
                                         observation_axis = c("rows", "columns"),
                                         extreme_threshold = 1e300) {
    observation_axis <- match.arg(observation_axis)
    dims <- dim(x)
    if (is.null(dims) || length(dims) != 2L) {
        stop(caller, ": input must be a matrix-like object.", call. = FALSE)
    }

    n_cells <- if (identical(observation_axis, "rows")) dims[[1]] else dims[[2]]
    if (!is.finite(n_cells) || n_cells < 2L) {
        stop("n_cells must be >= 2", call. = FALSE)
    }

    values <- if (inherits(x, "dgCMatrix")) {
        x@x
    } else if (is.matrix(x) || is.data.frame(x)) {
        as.numeric(x)
    } else {
        as.numeric(as.matrix(x))
    }

    if (length(values) && any(!is.finite(values))) {
        stop(caller, ": input contains non-finite values.", call. = FALSE)
    }
    if (length(values) && any(abs(values) > extreme_threshold, na.rm = TRUE)) {
        warning(
            caller,
            ": input contains extremely large finite values (>1e300); cross-products may lose precision. ",
            "Consider scaling or winsorizing before spatial statistics.",
            call. = FALSE,
            immediate. = TRUE
        )
    }

    invisible(x)
}
