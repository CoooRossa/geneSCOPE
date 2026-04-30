#' Normalize core compute backend names.
#' @description
#' Internal helper shared by compute paths that expose an explicit backend
#' policy. User-facing names are `auto`, `cpp`, `python`, and `r`; legacy
#' aliases remain accepted for compatibility.
#' @param backend Backend requested by the caller.
#' @param arg Argument name used in error messages.
#' @param allow_auto Whether `auto` is an accepted value.
#' @return Canonical lower-case backend label.
#' @keywords internal
.normalize_core_backend <- function(backend,
                                    arg = "backend",
                                    allow_auto = TRUE) {
    backend <- as.character(backend)
    if (length(backend) != 1L || is.na(backend) || !nzchar(backend)) {
        stop("`", arg, "` must be a single non-empty string.", call. = FALSE)
    }

    normalized <- tolower(backend)
    normalized <- switch(normalized,
        "native" = "cpp",
        "c++" = "cpp",
        "rcpp" = "cpp",
        "reference" = "r",
        normalized
    )

    choices <- c(if (isTRUE(allow_auto)) "auto", "cpp", "python", "r")
    if (!normalized %in% choices) {
        stop(
            "Unsupported `", arg, "` value '", backend,
            "'. Use one of: ", paste(choices, collapse = ", "),
            ". Legacy alias 'native' maps to 'cpp'.",
            call. = FALSE
        )
    }
    normalized
}

#' Core backend policy descriptor.
#' @description
#' Internal helper that records the intended order without implying that every
#' tier is implemented for every statistic.
#' @param requested Canonical requested backend.
#' @param python_supported Whether a real Python implementation is available.
#' @param cpp_supported Whether a C++ implementation is available.
#' @param r_supported Whether an R reference implementation is available.
#' @return A list suitable for structured reports.
#' @keywords internal
.core_backend_policy <- function(requested,
                                 python_supported = FALSE,
                                 cpp_supported = TRUE,
                                 r_supported = TRUE) {
    requested <- .normalize_core_backend(requested, allow_auto = TRUE)
    list(
        requested = requested,
        auto_order = c("cpp", "r"),
        supported = list(
            cpp = isTRUE(cpp_supported),
            python = isTRUE(python_supported),
            r = isTRUE(r_supported)
        ),
        aliases = list(native = "cpp", R = "r")
    )
}

.core_backend_python_unavailable_message <- function(function_name) {
    paste0(
        "[geneSCOPE::", function_name,
        "] backend='python' is not implemented for this statistic. ",
        "Current geneSCOPE core statistics use C++ first and fall back to the ",
        "R reference implementation on regular backend errors."
    )
}
