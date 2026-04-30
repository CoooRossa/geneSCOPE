#' OpenMP build/runtime info
#' @description
#' Returns a list describing whether the package was compiled with OpenMP, and
#' basic runtime capacity fields. On Darwin/macOS, native OpenMP diagnostics
#' obey the spatial native safe default and return a structured disabled
#' payload unless diagnostics are explicitly allowed.
#' @return Named list.
#' @keywords internal
.native_openmp_info <- function() {
    if (isTRUE(getOption("geneSCOPE.disable_native_openmp_diagnostics", FALSE))) {
        return(list(
            compiled_with_openmp = FALSE,
            omp_max_threads = NA_integer_,
            omp_num_procs = NA_integer_,
            diagnostic_mode = "disabled_by_option"
        ))
    }
    if (exists(".spw_darwin_native_spatial_disabled", mode = "function") &&
        isTRUE(.spw_darwin_native_spatial_disabled()) &&
        !isTRUE(getOption("geneSCOPE.allow_darwin_native_openmp_diagnostics", FALSE))) {
        return(list(
            compiled_with_openmp = FALSE,
            omp_max_threads = NA_integer_,
            omp_num_procs = NA_integer_,
            diagnostic_mode = "disabled_by_darwin_spatial_safe_default"
        ))
    }
    .Call(`_geneSCOPERebuild_native_openmp_info`)
}

#' Set OpenMP threads (best-effort)
#' @description
#' Calls omp_set_num_threads(n) when compiled with OpenMP; otherwise returns a
#' structured NA payload.
#' @param n_threads Integer requested thread count.
#' @return Named list with requested/max threads and compile status.
#' @keywords internal
.native_openmp_set_threads <- function(n_threads = 1L) {
    if (isTRUE(getOption("geneSCOPE.disable_native_openmp_diagnostics", FALSE))) {
        return(list(
            requested_threads = as.integer(n_threads),
            compiled_with_openmp = FALSE,
            omp_max_threads = NA_integer_,
            omp_num_procs = NA_integer_,
            diagnostic_mode = "disabled_by_option"
        ))
    }
    if (exists(".spw_darwin_native_spatial_disabled", mode = "function") &&
        isTRUE(.spw_darwin_native_spatial_disabled()) &&
        !isTRUE(getOption("geneSCOPE.allow_darwin_native_openmp_diagnostics", FALSE))) {
        return(list(
            requested_threads = as.integer(n_threads),
            compiled_with_openmp = FALSE,
            omp_max_threads = NA_integer_,
            omp_num_procs = NA_integer_,
            diagnostic_mode = "disabled_by_darwin_spatial_safe_default"
        ))
    }
    .Call(`_geneSCOPERebuild_native_openmp_set_threads`, as.integer(n_threads))
}
