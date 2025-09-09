#' Package Loading and Environment Setup
#'
#' @name geneSCOPE-package
#' @title Package Loading and Environment Setup
#' @description
#'   Package initialization functions and utility helpers for cross-platform
#'   compatibility and safe thread management.
#' @keywords internal
NULL

#' @title OS Detection
#' @description Detect the operating system
#' @return Character string: "windows", "macos", or "linux"
#' @export
detectOS <- function() {
  if (.Platform$OS.type == "windows") {
    "windows"
  } else if (Sys.info()["sysname"] == "Darwin") {
    "macos"
  } else {
    "linux"
  }
}

#' @title Get Safe Thread Count
#' @description Determines safe thread count based on platform and available cores
#' @param requested_cores Integer. Requested number of cores
#' @return Integer. Safe number of cores to use
#' @keywords internal
#' @export
getSafeThreadCount <- function(requested_cores) {
  os_type <- detectOS()
  max_cores <- parallel::detectCores()

  # Platform-specific limits
  safe_max <- switch(os_type,
    windows = min(16, max_cores),
    macos = min(20, max_cores),
    linux = min(24, max_cores)
  )

  # Leave at least one core for the OS
  system_reserve <- if (max_cores <= 4) 1 else min(2, max_cores %/% 8)
  safe_max <- max(1, safe_max - system_reserve)

  min(requested_cores, safe_max)
}

#' @title Fix Random Number Generation
#' @description Ensures reproducible random number generation across platforms
#' @param seed Integer. Random seed to set
#' @keywords internal
#' @export
fixRNG <- function(seed = 42) {
  set.seed(seed)
  # Also set OpenMP seed if available
  if (Sys.getenv("OMP_NUM_THREADS") != "") {
    Sys.setenv(OMP_RANDOM_SEED = as.character(seed))
  }
}

#' @title Select Grid Layer
#' @description Helper function to select a grid layer from coordObj
#' @param coordObj A CoordObj
#' @param grid_name Optional grid layer name
#' @return The selected grid layer
#' @keywords internal
#' @export
selectGridLayer <- function(coordObj, grid_name = NULL) {
  if (is.null(grid_name)) {
    if (length(coordObj@grid) == 1) {
      return(coordObj@grid[[1]])
    } else {
      stop("Multiple grid layers found. Please specify grid_name.")
    }
  } else {
    if (grid_name %in% names(coordObj@grid)) {
      return(coordObj@grid[[grid_name]])
    } else {
      stop("Grid layer '", grid_name, "' not found.")
    }
  }
}

# Internal compatibility aliases for C++ function name compatibility
.leeL_cpp_cache <- function(...) leeL_cpp_cache(...)
.leeL_cpp_cols <- function(...) leeL_cpp_cols(...)
.lee_perm_block_cpp <- function(...) lee_perm_block_cpp(...)

# Package initialization and thread management
#' @importFrom RhpcBLASctl blas_set_num_threads omp_set_num_threads get_num_procs
#' @importFrom parallel detectCores
#' @export
.onLoad <- function(libname, pkgname) {
  # 1) Force BLAS to single-thread to avoid OpenMP conflicts
  .forceSingleThreadBLAS()

  # 2) Set package options (migrate from old FG2CLI.* if present)
  op <- options()
  # New option names under geneSCOPE.*
  op_gs <- list(
    geneSCOPE.use_64bit = TRUE,
    geneSCOPE.max_matrix_size = 2^31 - 1,
    geneSCOPE.chunk_size = 1000,
    geneSCOPE.verbose = TRUE,
    geneSCOPE.blas_threads = 1
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
  .setBlasEnvironment()

  # 4) Check for 64-bit support via wrapper if available
  tryCatch(
    {
      if (exists("test_64bit_support", mode = "function")) {
        test_result <- test_64bit_support()
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

  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("[geneSCOPE] Gene Spatial Correlation Of Pairwise Expression")
  packageStartupMessage("[geneSCOPE] Version: ", utils::packageDescription(pkgname)$Version)

  # System capabilities
  ncores <- parallel::detectCores()
  packageStartupMessage("[geneSCOPE] Detected ", ncores, " CPU cores")

  # BLAS notice
  packageStartupMessage("[geneSCOPE] BLAS threading disabled to prevent OpenMP conflicts")

  # Migrate old option namespace notice
  if (any(nzchar(Sys.getenv(c("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))))) {
    # no-op; environment already set in .onLoad
  }
}

# Force BLAS single-threading
.forceSingleThreadBLAS <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    tryCatch(
      {
        RhpcBLASctl::blas_set_num_threads(1)
        RhpcBLASctl::omp_set_num_threads(1)
      },
      error = function(e) { }
    )
  }
  .setBlasEnvironment()
  Sys.setenv(
    OPENBLAS_NUM_THREADS = "1",
    GOTO_NUM_THREADS = "1",
    OMP_NUM_THREADS = "1"
  )
}

# Set BLAS-related environment variables
.setBlasEnvironment <- function() {
  blas_vars <- c(
    "OPENBLAS_NUM_THREADS" = "1",
    "MKL_NUM_THREADS" = "1",
    "VECLIB_MAXIMUM_THREADS" = "1",
    "NUMEXPR_NUM_THREADS" = "1",
    "BLIS_NUM_THREADS" = "1",
    "GOTO_NUM_THREADS" = "1",
    "ATLAS_NUM_THREADS" = "1"
  )

  for (var in names(blas_vars)) {
    do.call(Sys.setenv, setNames(list(blas_vars[[var]]), var))
  }
}

.onUnload <- function(libpath) {
  invisible()
}

# Thread management system
.init_thread_manager <- function() {
  original_settings <- list(
    blas_threads = if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) 1)
    } else {
      1
    },
    env_vars = sapply(c(
      "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
      "VECLIB_MAXIMUM_THREADS", "BLAS_NUM_THREADS"
    ), Sys.getenv, unset = NA_character_)
  )

  assign(".original_thread_settings", original_settings, envir = .GlobalEnv)
  .set_conservative_threads()
}

.set_conservative_threads <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    tryCatch(
      {
        RhpcBLASctl::blas_set_num_threads(1)
      },
      error = function(e) { }
    )
  }
  thread_vars <- c(
    "OMP_NUM_THREADS" = "1",
    "OPENBLAS_NUM_THREADS" = "1",
    "MKL_NUM_THREADS" = "1",
    "VECLIB_MAXIMUM_THREADS" = "1",
    "BLAS_NUM_THREADS" = "1"
  )
  do.call(Sys.setenv, as.list(thread_vars))
  if (requireNamespace("data.table", quietly = TRUE)) {
    tryCatch(
      {
        data.table::setDTthreads(1)
      },
      error = function(e) { }
    )
  }
}

#' @title Smart thread allocation for different computational tasks
#' @description Determines optimal thread counts for different types of parallel operations
#' @param ncores_requested Integer. Requested number of cores
#' @param task_type Character. One of "openmp_only","blas_heavy","mixed","io_bound"
#' @param memory_gb Numeric. Available memory in GB
#' @return List with optimal thread allocations
#' @export
getOptimalThreadConfig <- function(ncores_requested = parallel::detectCores(),
                                   task_type = "mixed",
                                   memory_gb = NULL) {
  max_cores <- parallel::detectCores()
  os_type <- detectOS()

  if (is.null(memory_gb)) {
    memory_gb <- tryCatch(
      {
        if (os_type == "linux") {
          mem_info <- readLines("/proc/meminfo")
          mem_total <- grep("MemTotal", mem_info, value = TRUE)
          as.numeric(gsub("[^0-9]", "", mem_total)) / 1024 / 1024
        } else if (os_type == "macos") {
          mem_bytes <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
          mem_bytes / 1024^3
        } else {
          16
        }
      },
      error = function(e) 16
    )
  }

  platform_limit <- switch(os_type,
    windows = min(16, max_cores),
    macos = min(20, max_cores),
    linux = min(24, max_cores)
  )

  memory_limit <- switch(task_type,
    "openmp_only" = max(1, floor(memory_gb / 1.5)),
    "blas_heavy" = max(1, floor(memory_gb / 2)),
    "mixed" = max(1, floor(memory_gb / 2.5)),
    "io_bound" = max(1, floor(memory_gb / 1))
  )

  thread_config <- switch(task_type,
    "openmp_only" = list(openmp_threads = min(ncores_requested, platform_limit, memory_limit), blas_threads = 1, r_threads = 1),
    "blas_heavy" = list(openmp_threads = 1, blas_threads = min(ncores_requested, platform_limit, memory_limit), r_threads = 1),
    "mixed" = {
      total_safe <- min(ncores_requested, platform_limit, memory_limit)
      openmp_share <- max(1, ceiling(total_safe * 0.7))
      blas_share <- max(1, min(4, total_safe - openmp_share + 1))
      list(openmp_threads = openmp_share, blas_threads = blas_share, r_threads = 1)
    },
    "io_bound" = list(openmp_threads = min(8, max_cores), blas_threads = 1, r_threads = min(ncores_requested, 12))
  )

  total_threads <- thread_config$openmp_threads + thread_config$blas_threads
  max_allowed <- max_cores - 1
  if (total_threads > max_allowed) {
    scale_factor <- max_allowed / total_threads
    thread_config$openmp_threads <- max(1, floor(thread_config$openmp_threads * scale_factor))
    thread_config$blas_threads <- max(1, min(4, floor(thread_config$blas_threads * scale_factor)))
  }

  thread_config$total_used <- thread_config$openmp_threads + thread_config$blas_threads
  thread_config$efficiency <- thread_config$total_used / max_cores
  thread_config$task_type <- task_type
  return(thread_config)
}

#' @title Configure threads for specific computational task
#' @description Sets up optimal thread configuration for a specific task type
#' @inheritParams getOptimalThreadConfig
#' @param restore_after Logical. Restore previous settings after execution
#' @return Previous thread configuration (invisibly)
#' @export
configureThreadsFor <- function(task_type = "mixed",
                                ncores_requested = parallel::detectCores(),
                                restore_after = TRUE) {
  current_config <- list(
    blas_threads = if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) 1)
    } else {
      1
    },
    env_vars = sapply(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"), Sys.getenv, unset = NA_character_)
  )

  optimal <- getOptimalThreadConfig(ncores_requested, task_type)

  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    tryCatch(
      {
        RhpcBLASctl::blas_set_num_threads(optimal$blas_threads)
      },
      error = function(e) {
        message("[geneSCOPE] !!! Warning: Could not set BLAS threads: ", e$message, " !!!")
      }
    )
  }

  do.call(Sys.setenv, list(OMP_NUM_THREADS = as.character(optimal$openmp_threads)))
  thread_env <- list(
    OPENBLAS_NUM_THREADS = as.character(optimal$blas_threads),
    MKL_NUM_THREADS = as.character(optimal$blas_threads)
  )
  do.call(Sys.setenv, thread_env)

  if (restore_after) {
    attr(optimal, "restore_function") <- function() {
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
        tryCatch(
          {
            RhpcBLASctl::blas_set_num_threads(current_config$blas_threads)
          },
          error = function(e) { }
        )
      }
      for (var in names(current_config$env_vars)) {
        if (is.na(current_config$env_vars[[var]])) {
          Sys.unsetenv(var)
        } else {
          do.call(Sys.setenv, setNames(list(current_config$env_vars[[var]]), var))
        }
      }
    }
  }

  invisible(optimal)
}

#' @title Execute function with specific thread configuration
#' @description Temporarily sets thread configuration for a function execution
#' @param expr Expression to execute
#' @param task_type Character. Type of computational task
#' @param ncores_requested Integer. Requested cores
#' @return Result of expression execution
#' @export
withThreadConfig <- function(expr, task_type = "mixed", ncores_requested = parallel::detectCores()) {
  config <- configureThreadsFor(task_type, ncores_requested, restore_after = TRUE)
  on.exit({
    restore_fn <- attr(config, "restore_function")
    if (!is.null(restore_fn)) restore_fn()
  })

  eval(expr, envir = parent.frame())
}
