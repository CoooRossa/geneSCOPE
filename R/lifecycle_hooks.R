
#' On Attach
#' @description
#' Internal helper for `.on_attach`.
#' @param libname Parameter value.
#' @param pkgname Parameter value.
#' @return Return value used internally.
#' @keywords internal
.on_attach <- function(libname, pkgname) {
  packageStartupMessage("[geneSCOPE] Gene Spatial Correlation Of Pairwise Expression")
  packageStartupMessage("[geneSCOPE] Version: ", packageDescription(pkgname)$Version)

  # System capabilities
  ncores <- detectCores()
  packageStartupMessage("[geneSCOPE] Detected ", ncores, " CPU cores")

  # BLAS notice
  packageStartupMessage("[geneSCOPE] BLAS threading disabled to prevent OpenMP conflicts")
  
  # Round 4: Darwin native safety - check if native is disabled
  is_darwin <- Sys.info()["sysname"] == "Darwin"
  native_disabled <- isTRUE(getOption("geneSCOPE.disable_native_all", FALSE))
  
  if (is_darwin && native_disabled) {
    packageStartupMessage(
      "[geneSCOPE] Darwin/macOS detected: Native C++ backends DISABLED by default. ",
      "Using R fallback for all native paths to prevent crashes. ",
      "To re-enable native C++, set options(geneSCOPE.disable_native_all=FALSE)."
    )
  } else if (is_darwin) {
    packageStartupMessage(
      "[geneSCOPE] Darwin/macOS detected: Native C++ backends ENABLED (experimental). ",
      "If you experience crashes, set options(geneSCOPE.disable_native_all=TRUE) before loading geneSCOPE."
    )
  } else {
    packageStartupMessage(
      "[geneSCOPE] Native C++ backends are enabled by default (C++ first, R fallback). ",
      "If a native backend causes a fatal process error, set ",
      "options(geneSCOPE.disable_native_all=TRUE) before loading geneSCOPE."
    )
  }

  # Migrate old option namespace notice
  if (any(nzchar(Sys.getenv(c("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))))) {
    # no-op; environment already set in .on_load
  }
}

#' On Load
#' @description
#' Internal helper for `.on_load`.
#' @param libname Parameter value.
#' @param pkgname Parameter value.
#' @return Return value used internally.
#' @keywords internal
.on_load <- function(libname, pkgname) {
  # 1) Force BLAS to single-thread to avoid OpenMP conflicts
  .force_single_thread_blas()

  # Round 4: Darwin native safety - detect platform and set defaults
  is_darwin <- Sys.info()["sysname"] == "Darwin"
  
  # 2) Set package options (migrate from old FG2CLI.* if present)
  op <- options()
  # New option names under geneSCOPE.*
  # Round 4 Fix: On Darwin, disable all native C++ paths by default to prevent crashes
  # Users can override by setting options(geneSCOPE.disable_native_all=FALSE) before loading
  native_all_default <- if (!is.null(getOption("geneSCOPE.disable_native_all", NULL))) {
    isTRUE(getOption("geneSCOPE.disable_native_all"))
  } else {
    is_darwin
  }
  op_gs <- list(
    geneSCOPE.use_64bit = TRUE,
    geneSCOPE.max_matrix_size = 2^31 - 1,
    geneSCOPE.chunk_size = 1000,
    geneSCOPE.verbose = TRUE,
    geneSCOPE.blas_threads = 1,
    geneSCOPE.enable_64bit_check = FALSE,
    geneSCOPE.lee_l.backend = if (is_darwin) "R" else "cpp",  # Darwin default to R backend
    geneSCOPE.disable_native_all = native_all_default,  # Darwin default TRUE, others FALSE
    geneSCOPE.disable_native_grid_nb = native_all_default,
    geneSCOPE.disable_native_listw_builder = native_all_default,
    geneSCOPE.disable_native_kernel_weights = native_all_default,
    geneSCOPE.disable_native_lee_l_backend = native_all_default,
    geneSCOPE.disable_native_correlation_backend = native_all_default,
    geneSCOPE.disable_native_consensus = native_all_default,
    geneSCOPE.disable_native_permutation = native_all_default,  # New: prevent lee_perm_block crash
    geneSCOPE.disable_native_openmp_diagnostics = FALSE,
    geneSCOPE.allow_darwin_native_openmp_diagnostics = FALSE,
    geneSCOPE.allow_darwin_native_spatial = FALSE
  )
  # Backward-compat: if old options are set, prefer their values
  old_keys <- c("use_64bit", "max_matrix_size", "chunk_size", "verbose", "blas_threads")
  for (k in old_keys) {
    old_name <- paste0("FG2CLI.", k)
    new_name <- paste0("geneSCOPE.", k)
    if (!is.null(getOption(old_name, NULL))) {
      op_gs[[new_name]] <- getOption(old_name)
    }
  }
  toset <- !(names(op_gs) %in% names(op))
  if (any(toset)) options(op_gs[toset])

  # 3) Set BLAS-related environment variables
  .set_blas_environment()

  # 4) Optional 64-bit support probe (disabled by default; opt-in via option/env)
  run_64bit_check <- isTRUE(getOption("geneSCOPE.enable_64bit_check")) ||
    identical(Sys.getenv("GENESCOPE_ENABLE_64BIT_TEST", ""), "1")
  skip_64bit_check <- identical(Sys.getenv("GENESCOPE_SKIP_64BIT_TEST", ""), "1")
  if (run_64bit_check && !skip_64bit_check) {
    tryCatch(
      {
        if (exists(".test_64bit_support", mode = "function")) {
          test_result <- .test_64bit_support()
          if (!isTRUE(test_result)) {
            packageStartupMessage(
              "[geneSCOPE] !!! Warning: 64-bit matrix indexing may not be fully supported. ",
              "Large matrices may fail. Consider recompiling with ARMA_64BIT_WORD=1 !!!"
            )
          }
        }
      },
      error = function(e) { }
    )
  }

  invisible()
}

#' On Unload
#' @description
#' Internal helper for `.on_unload`.
#' @param libpath Parameter value.
#' @return Return value used internally.
#' @keywords internal
.on_unload <- function(libpath) {
  invisible()
}
