# Internal helpers for runtime dependency governance.
#' Ensure Arrow Available
#' @description
#' Internal helper for `.ensure_arrow_available`.
#' @param context String label used for tracing/debug messages.
#' @return Return value used internally.
#' @keywords internal
.ensure_arrow_available <- function(context = "[geneSCOPE::.create_scope]") {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop(
      paste0(
        context, " Missing runtime dependency 'arrow'. ",
        "Install it via install.packages('arrow') or supply parquet outputs from an environment where Arrow is available."
      ),
      call. = FALSE
    )
  }
}

#' With Data Table Attached
#' @description
#' Internal helper for `.with_data_table_attached`.
#' @param expr Parameter value.
#' @param context String label used for tracing/debug messages.
#' @return Return value used internally.
#' @keywords internal
.with_data_table_attached <- function(expr, context = "[geneSCOPE::runtime]") {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(
      paste0(
        context, " Missing runtime dependency 'data.table'. ",
        "Install it via install.packages('data.table') to enable ':=' semantics."
      ),
      call. = FALSE
    )
  }
  if (!is.function(expr)) {
    stop(context, " .with_data_table_attached requires a function argument", call. = FALSE)
  }

  already_attached <- "package:data.table" %in% search()
  if (!already_attached) {
    # This is an intentional design choice for data.table NSE support
    # The package is in Imports, but we need it in search path for := syntax
    suppressPackageStartupMessages(library(data.table, warn.conflicts = FALSE))
  }
  on.exit({
    if (!already_attached && "package:data.table" %in% search()) {
      detach("package:data.table", unload = FALSE, character.only = TRUE)
    }
  }, add = TRUE)

  expr()
}

#' Thread Guard Trace Enabled
#' @description
#' Internal helper for `.thread_guard_trace_enabled`.
#' @return Return value used internally.
#' @keywords internal
.thread_guard_trace_enabled <- local({
  .cache <- NULL
  function() {
    if (is.null(.cache)) {
      flag <- Sys.getenv("GENESCOPE_THREAD_GUARD_TRACE", "")
      .cache <<- tolower(flag) %in% c("1", "true", "yes", "on")
    }
    .cache
  }
})

#' Thread Guard Trace
#' @description
#' Internal helper for `.thread_guard_trace`.
#' @param context String label used for tracing/debug messages.
#' @param before Parameter value.
#' @param after Parameter value.
#' @param arrow_threads Parameter value.
#' @param builder_threads Parameter value.
#' @param requested_ncores Parameter value.
#' @param other_threads Parameter value.
#' @return Return value used internally.
#' @keywords internal
.thread_guard_trace <- function(context,
                                before,
                                after,
                                arrow_threads = NA_integer_,
                                builder_threads = NA_integer_,
                                requested_ncores = NA_integer_,
                                other_threads = NA_integer_) {
  if (!.thread_guard_trace_enabled()) return(invisible(NULL))

  fmt <- function(x) {
    if (is.null(x) || (length(x) && is.na(x))) return("NA")
    if (identical(x, "")) return("unset")
    paste0(x)
  }

  message(
    "[thread_guard] ",
    context,
    " | requested_ncores=", fmt(requested_ncores),
    " builder_threads=", fmt(builder_threads),
    if (!is.na(other_threads)) paste0(" other_threads=", fmt(other_threads)) else "",
    " | arrow_threads=", fmt(arrow_threads),
    " | before{OMP=", fmt(before$omp_num_threads),
    " OPENBLAS=", fmt(before$openblas_threads),
    " MKL=", fmt(before$mkl_threads),
    " VECLIB=", fmt(before$veclib_threads),
    " BLIS=", fmt(before$blis_threads),
    " NUMEXPR=", fmt(before$numexpr_threads),
    " RCPP_PARALLEL=", fmt(before$rcpp_parallel_threads),
    "} after{OMP=", fmt(after$omp_num_threads),
    " OPENBLAS=", fmt(after$openblas_threads),
    " MKL=", fmt(after$mkl_threads),
    " VECLIB=", fmt(after$veclib_threads),
    " BLIS=", fmt(after$blis_threads),
    " NUMEXPR=", fmt(after$numexpr_threads),
    " RCPP_PARALLEL=", fmt(after$rcpp_parallel_threads),
    "}"
  )
}

#' Deterministic Arrow Threads
#' @description
#' Internal helper for `.deterministic_arrow_threads`.
#' @param ncores_requested Parameter value.
#' @return Return value used internally.
#' @keywords internal
.deterministic_arrow_threads <- function(ncores_requested) {
  cores_detected <- suppressWarnings(detectCores(logical = TRUE))
  cores_detected <- if (length(cores_detected) && is.finite(cores_detected)) cores_detected else 2L
  requested <- suppressWarnings(as.integer(ncores_requested))
  requested <- if (length(requested) && is.finite(requested) && requested > 0L) requested else cores_detected
  max(2L, min(requested, cores_detected))
}

#' Apply Thread Config
#' @description
#' Internal helper for `.apply_thread_config`.
#' @param config Parameter value.
#' @return Return value used internally.
#' @keywords internal
.apply_thread_config <- function(config) {
  # Configure OpenMP threads
  Sys.setenv(OMP_NUM_THREADS = as.character(config$openmp_threads))
  if (!is.null(config$openmp_threads) && is.finite(config$openmp_threads)) {
    try(.native_openmp_set_threads(as.integer(config$openmp_threads)), silent = TRUE)
  }
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = as.character(config$r_threads))

  # Ensure BLAS remains single-threaded
  .force_blas_single_thread()

  # Configure Arrow threads when available
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow_threads <- config[["arrow_threads"]]
    arrow_threads <- .deterministic_arrow_threads(if (is.null(arrow_threads)) config[["ncores_requested"]] else arrow_threads)
    tryCatch(
      {
        if (exists("set_cpu_count", where = asNamespace("arrow"))) {
          arrow::set_cpu_count(arrow_threads)
        }
        if (exists("set_io_thread_count", where = asNamespace("arrow"))) {
          arrow::set_io_thread_count(arrow_threads)
        }
      },
      error = function(e) NULL
    )
  }
}

#' Configure for compute-intensive operations
#' @description
#' Internal helper for `.configure_compute_intensive`.
#' @param ncores_requested Parameter value.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.configure_compute_intensive <- function(ncores_requested) {
  max_cores <- max(1L, detectCores())
  safe_cores <- max(1L, min(ncores_requested, max_cores))

  list(
    openmp_threads = safe_cores,
    r_threads = 1, # keep R single-threaded
    blas_threads = 1 # force BLAS single-thread
  )
}

#' Configure for I/O bound operations
#' @description
#' Internal helper for `.configure_io_bound`.
#' @param ncores_requested Parameter value.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.configure_io_bound <- function(ncores_requested) {
  max_cores <- max(1L, detectCores())
  safe_cores <- max(1L, min(ncores_requested, max_cores))

  list(
    openmp_threads = safe_cores, # use all requested cores for OpenMP
    r_threads = safe_cores, # parallel R-level tasks can also use requested cores
    blas_threads = 1
  )
}

#' Configure for mixed operations
#' @description
#' Internal helper for `.configure_mixed`.
#' @param ncores_requested Parameter value.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.configure_mixed <- function(ncores_requested) {
  max_cores <- max(1L, detectCores())
  safe_cores <- max(1L, min(ncores_requested, max_cores))
  openmp_share <- safe_cores

  list(
    openmp_threads = openmp_share,
    r_threads = max(1, safe_cores - openmp_share),
    blas_threads = 1
  )
}

#' Configure for OpenMP-only operations
#' @description
#' Internal helper for `.configure_open_mp_only`.
#' @param ncores_requested Parameter value.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.configure_open_mp_only <- function(ncores_requested) {
  max_cores <- max(1L, detectCores())
  safe_cores <- max(1L, min(ncores_requested, max_cores))

  list(
    openmp_threads = safe_cores,
    r_threads = 1,
    blas_threads = 1 # still force single-thread BLAS
  )
}

#' Force Blas Single Thread
#' @description
#' Internal helper for `.force_blas_single_thread`.
#' @return Return value used internally.
#' @keywords internal
.force_blas_single_thread <- function() {
  # 1. RhpcBLASctl method
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    tryCatch(
      {
        RhpcBLASctl::blas_set_num_threads(1)
      },
      error = function(e) NULL
    )
  }

  # 2. Environment variable method (more reliable)
  blas_vars <- c(
    "OPENBLAS_NUM_THREADS" = "1",
    "MKL_NUM_THREADS" = "1",
    "VECLIB_MAXIMUM_THREADS" = "1",
    "NUMEXPR_NUM_THREADS" = "1",
    "BLIS_NUM_THREADS" = "1",
    "GOTO_NUM_THREADS" = "1",
    "ATLAS_NUM_THREADS" = "1",
    "LAPACK_NUM_THREADS" = "1"
  )

  for (var in names(blas_vars)) {
    do.call(Sys.setenv, setNames(list(blas_vars[[var]]), var))
  }

  # 3. Set data.table threads to 1
  if (requireNamespace("data.table", quietly = TRUE)) {
    tryCatch(
      {
        setDTthreads(1)
      },
      error = function(e) NULL
    )
  }
}

#' Fix Rng
#' @description
#' Internal helper for `.fix_rng`.
#' @param seed Random seed.
#' @return Return value used internally.
#' @keywords internal
.fix_rng <- function(seed = 42) {
  set.seed(seed)
  # Also set OpenMP seed if available
  if (Sys.getenv("OMP_NUM_THREADS") != "") {
    Sys.setenv(OMP_RANDOM_SEED = as.character(seed))
  }
}

#' Detect Os
#' @description
#' Internal helper for `.detect_os`.
#' @keywords internal
#' @return Return value used internally.
.detect_os <- function() {
  if (.Platform$OS.type == "windows") {
    return("windows")
  } else if (Sys.info()["sysname"] == "Darwin") {
    return("macos")
  } else {
    return("linux")
  }
}

#' Detect Container
#' @description
#' Internal helper for `.detect_container`.
#' @keywords internal
#' @return Return value used internally.
.detect_container <- function() {
  if (.Platform$OS.type == "windows") return(FALSE)
  if (!identical(Sys.info()[["sysname"]], "Linux")) return(FALSE)
  if (file.exists("/.dockerenv")) return(TRUE)
  if (nzchar(Sys.getenv("container", ""))) return(TRUE)

  cgroup <- tryCatch(readLines("/proc/1/cgroup", warn = FALSE), error = function(e) character())
  if (length(cgroup) == 0L) return(FALSE)
  any(grepl("docker|kubepods|containerd|podman|lxc", cgroup))
}

#' Get Optimal Thread Config
#' @description
#' Internal helper for `.get_optimal_thread_config`.
#' @param ncores_requested Parameter value.
#' @param task_type Parameter value.
#' @param memory_gb Parameter value.
#' @return Return value used internally.
#' @keywords internal
.get_optimal_thread_config <- function(ncores_requested = detectCores(),
                                    task_type = "mixed",
                                    memory_gb = NULL) {
  max_cores <- detectCores()

  safe_cores <- max(1L, min(ncores_requested, max_cores))

  thread_config <- switch(task_type,
    "openmp_only" = list(openmp_threads = safe_cores, blas_threads = 1, r_threads = 1),
    "blas_heavy" = list(openmp_threads = 1, blas_threads = safe_cores, r_threads = 1),
    "mixed" = {
      openmp_share <- max(1, ceiling(safe_cores * 0.7))
      r_share <- max(1, safe_cores - openmp_share)
      list(openmp_threads = openmp_share, blas_threads = 1, r_threads = r_share)
    },
    "io_bound" = list(openmp_threads = max(1, min(8, safe_cores)), blas_threads = 1, r_threads = safe_cores)
  )

  thread_config$total_used <- safe_cores
  thread_config$efficiency <- safe_cores / max_cores
  thread_config$task_type <- task_type
  return(thread_config)
}

#' Get Safe Thread Count
#' @description
#' Internal helper for `.get_safe_thread_count`.
#' @param max_requested Numeric threshold.
#' @param default Parameter value.
#' @return Return value used internally.
#' @keywords internal
.get_safe_thread_count <- function(max_requested = detectCores(),
                                default = NULL) {
  total_cores <- detectCores()
  safe_limit <- max(1L, total_cores - 1L)
  requested <- if (!is.null(default)) default else max_requested
  req_val <- suppressWarnings(as.integer(requested)[1])
  if (!length(req_val) || is.na(req_val)) req_val <- safe_limit
  return(min(req_val, safe_limit, na.rm = TRUE))
}

#' Configure Threads For
#' @description
#' Internal helper for `.configure_threads_for`.
#' @param operation_type Parameter value.
#' @param ncores_requested Parameter value.
#' @param restore_after Parameter value.
#' @return A list with configuration/metadata used internally.
#' @keywords internal
.configure_threads_for <- function(operation_type = c("compute_intensive", "io_bound", "mixed", "openmp_only"),
                                 ncores_requested = 1,
                                 restore_after = FALSE) {
  operation_type <- match.arg(operation_type)

  # Save current state
  old_state <- .save_thread_state()

  # Force BLAS single thread to avoid OpenBLAS warnings
  .force_blas_single_thread()

  # Configure based on operation type
  config <- switch(operation_type,
    compute_intensive = .configure_compute_intensive(ncores_requested),
    io_bound = .configure_io_bound(ncores_requested),
    mixed = .configure_mixed(ncores_requested),
    openmp_only = .configure_open_mp_only(ncores_requested)
  )

  # Apply configuration
  .apply_thread_config(config)

  # Prepare return value
  result <- list(
    operation_type = operation_type,
    openmp_threads = config$openmp_threads,
    r_threads = config$r_threads,
    blas_threads = 1 # always 1
  )

  if (restore_after) {
    restore_function <- function() {
      .restore_thread_state(old_state)
    }
    attr(result, "restore_function") <- restore_function
  }

  return(result)
}

#' With Thread Config
#' @description
#' Internal helper for `.with_thread_config`.
#' @param expr Parameter value.
#' @param task_type Parameter value.
#' @param ncores_requested Parameter value.
#' @return Return value used internally.
#' @keywords internal
.with_thread_config <- function(expr, task_type = "mixed", ncores_requested = detectCores()) {
  config <- .configure_threads_for(task_type, ncores_requested, restore_after = TRUE)
  on.exit({
    restore_fn <- attr(config, "restore_function")
    if (!is.null(restore_fn)) restore_fn()
  })

  eval(expr, envir = parent.frame())
}

#' Configure Threads For V2
#' @description
#' Internal helper for `.configure_threads_for_v2`.
#' @param ... Additional arguments (currently unused).
#' @keywords internal
#' @return A list with configuration/metadata used internally.
.configure_threads_for_v2 <- function(...) {
  .configure_threads_for(...)
}

#' With Thread Config V2
#' @description
#' Internal helper for `.with_thread_config_v2`.
#' @param ... Additional arguments (currently unused).
#' @keywords internal
#' @return Return value used internally.
.with_thread_config_v2 <- function(...) {
  .with_thread_config(...)
}

#' Get Safe Thread Count V2
#' @description
#' Internal helper for `.get_safe_thread_count_v2`.
#' @param ... Additional arguments (currently unused).
#' @keywords internal
#' @return Return value used internally.
.get_safe_thread_count_v2 <- function(...) {
  .get_safe_thread_count(...)
}

#' Get Optimal Thread Config V2
#' @description
#' Internal helper for `.get_optimal_thread_config_v2`.
#' @param ... Additional arguments (currently unused).
#' @keywords internal
#' @return Return value used internally.
.get_optimal_thread_config_v2 <- function(...) {
  .get_optimal_thread_config(...)
}

#' Detect Os V2
#' @description
#' Internal helper for `.detect_os_v2`.
#' @param ... Additional arguments (currently unused).
#' @keywords internal
#' @return Return value used internally.
.detect_os_v2 <- function(...) {
  .detect_os(...)
}

#' Fix Rng V2
#' @description
#' Internal helper for `.fix_rng_v2`.
#' @param ... Additional arguments (currently unused).
#' @keywords internal
#' @return Return value used internally.
.fix_rng_v2 <- function(...) {
  .fix_rng(...)
}

#' Get System Memory Gb
#' @description
#' Internal helper for `.get_system_memory_gb`.
#' @return Return value used internally.
#' @keywords internal
.get_system_memory_gb <- function() {
  os_type <- .detect_os()
  tryCatch({
    if (os_type == "linux") {
      mem_info <- readLines("/proc/meminfo")
      mem_total <- grep("MemTotal", mem_info, value = TRUE)
      mem_kb <- as.numeric(gsub("[^0-9]", "", mem_total))
      mem_kb / 1024 / 1024
    } else if (os_type == "macos") {
      mem_bytes <- as.numeric(system("sysctl -n hw.memsize", intern = TRUE))
      mem_bytes / 1024^3
    } else {
      # Windows or fallback; conservatively assume 32GB if unknown
      32
    }
  }, error = function(e) 32)
}

#' Restore thread state
#' @description
#' Internal helper for `.restore_thread_state`.
#' @param old_state Parameter value.
#' @return Return value used internally.
#' @keywords internal
.restore_thread_state <- function(old_state) {
  # Restore environment variables
  for (var in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")) {
    key <- switch(var,
      "OMP_NUM_THREADS" = "omp_num_threads",
      "OPENBLAS_NUM_THREADS" = "openblas_threads",
      "MKL_NUM_THREADS" = "mkl_threads"
    )

    if (is.na(old_state[[key]])) {
      Sys.unsetenv(var)
    } else {
      do.call(Sys.setenv, setNames(list(old_state[[key]]), var))
    }
  }

  # Keep BLAS single-threaded to avoid conflicts
  .force_blas_single_thread()
}

#' Save Thread State
#' @description
#' Internal helper for `.save_thread_state`.
#' @return Return value used internally.
#' @keywords internal
.save_thread_state <- function() {
  list(
    omp_num_threads = Sys.getenv("OMP_NUM_THREADS", unset = NA),
    openblas_threads = Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA),
    mkl_threads = Sys.getenv("MKL_NUM_THREADS", unset = NA),
    veclib_threads = Sys.getenv("VECLIB_MAXIMUM_THREADS", unset = NA),
    blis_threads = Sys.getenv("BLIS_NUM_THREADS", unset = NA),
    numexpr_threads = Sys.getenv("NUMEXPR_NUM_THREADS", unset = NA),
    rcpp_parallel_threads = Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA),
    blas_threads = if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) 1)
    } else {
      1
    }
  )
}

#' Logging Enabled
#' @description
#' Internal helper for .log_enabled.
#' @param verbose Logical; whether to emit progress messages.
#' @return Logical flag indicating whether logging is enabled.
#' @keywords internal
.log_enabled <- function(verbose) {
    isTRUE(verbose)
}

#' Resolve Public Runtime Log Verbosity
#' @description
#' Resolves the `verbose` flag for runtime logging. Key logs are always
#' emitted by `.log_key()`, while detailed step logs use this resolved flag.
#' @param verbose Compact/public verbosity flag.
#' @return Logical flag for detailed logs.
#' @keywords internal
.resolve_runtime_verbose <- function(verbose = TRUE) {
    isTRUE(verbose)
}

#' Detailed Runtime Log Enabled
#' @description
#' Internal helper for distinguishing compact public logs from full debug traces.
#' Always returns FALSE now that debug_verbose has been removed.
#' @param verbose Logical verbosity flag (unused, kept for API compatibility).
#' @return Logical flag indicating whether detailed logs are enabled (always FALSE).
#' @keywords internal
.log_detail_enabled <- function(verbose) {
    FALSE
}

#' Compact Runtime Log Enabled
#' @description
#' Internal helper for compact public logging mode.
#' @param verbose Logical verbosity flag.
#' @return Logical flag indicating whether compact logs should be filtered.
#' @keywords internal
.log_compact_enabled <- function(verbose) {
    .log_enabled(verbose) && !.log_detail_enabled(verbose)
}

#' Runtime Log Importance
#' @description
#' Lightweight classifier used to suppress low-value compact logs while keeping
#' warnings, failures, fallback decisions, and inference-critical summaries.
#' @param parent Parent function name.
#' @param step Step id.
#' @param msg Message body.
#' @param kind Log kind (`info`, `backend`, `step`).
#' @return Integer importance score from 1 to 5.
#' @keywords internal
.log_importance_score <- function(parent, step = NULL, msg = "", kind = "info") {
    text <- tolower(paste(parent, step, msg, kind, collapse = " "))
    if (grepl(
        paste(
            c(
                "error", "failed", "failure", "fatal", "unsupported",
                "not supported", "invalid", "missing", "cannot", "disabled",
                "non-finite", "out of bounds", "trust", "signature"
            ),
            collapse = "|"
        ),
        text
    )) {
        return(5L)
    }
    if (grepl(
        paste(
            c(
                "warning", "warn", "fallback", "falling back", "retry",
                "skipped", "native", "permutation", "fdr",
                "q-value", "q-guard", "approximate_q", "guard",
                "bigmemory", "dense", "memory", "c\\+\\+", "cpp",
                "lee_l", "lee's l"
            ),
            collapse = "|"
        ),
        text
    )) {
        return(4L)
    }
    if (identical(kind, "step")) return(3L)
    if (grepl("thread|openmp|blas|ncores|cores|debug|diagnostic", text)) return(2L)
    1L
}

#' Compact Step Visibility
#' @description
#' Keeps compact public logs focused on long-running, user-relevant steps. Fast
#' helper APIs such as `computeDensity()` remain silent unless warnings/errors
#' occur or full debug logging is requested.
#' @param parent Parent function name.
#' @param step Step id.
#' @return Logical flag.
#' @keywords internal
.log_compact_step_visible <- function(parent, step) {
    step_id <- .log_step_id(step)
    visible_steps <- list(
        createSCOPE = c("S04", "S11", "S12", "S12b", "S13", "S14", "S15"),
        addSingleCells = c("S01"),
        normalizeMoleculesInGrid = c("S02", "S04", "S05"),
        normalizeSingleCells = c("S03", "S04"),
        computeWeights = c("S03", "S05"),
        computeL = c("S03", "S05", "S06", "S08"),
        computeCorrelation = c("S04", "S05", "S06"),
        compareLeesL = c("S02", "S03", "S04"),
        clusterGenes = c("S05", "S07", "S08", "S09"),
        computeLvsRCurve = c("S02", "S06")
    )
    allowed <- visible_steps[[as.character(parent)[1]]]
    !is.null(allowed) && step_id %in% allowed
}

#' Compact Step Field Relevance
#' @description
#' Compact step logs are only useful when they show user-facing inputs or
#' important outputs. This helper suppresses normal intermediate ENTER/DONE
#' lines that do not contain such fields.
#' @param extra Additional step text passed by the caller.
#' @param phase Step phase (`enter` or `done`).
#' @return Logical flag.
#' @keywords internal
.log_compact_step_extra_relevant <- function(extra, phase = c("enter", "done")) {
    if (is.null(extra) || !nzchar(extra)) return(FALSE)
    phase <- match.arg(phase)
    text <- tolower(as.character(extra)[1])
    param_keys <- c(
        "grid", "grid_name", "layer", "layer_name", "level", "method",
        "mode", "backend", "style", "topology", "ncores", "threads",
        "perms", "block_size", "mem_limit", "use_bigmemory", "chunk_size",
        "norm_layer", "resolution", "pct_min", "q_span", "threshold",
        "approximate_q", "sctransform", "prefer", "detected"
    )
    result_keys <- c(
        "dims", "nnz", "sig", "assigned", "cluster", "cluster_name",
        "cells", "molecules", "grids", "grid_layers", "n_genes",
        "n_grids", "retained", "stored", "layer", "mem_mode",
        "threads_final", "fdr_main_method", "backend_selected", "ok",
        "outputs", "count", "counts"
    )
    keys <- if (identical(phase, "enter")) param_keys else c(param_keys, result_keys)
    grepl(paste0("(^|[ _-])(", paste(keys, collapse = "|"), ")(=|[ _-]|$)"), text, perl = TRUE)
}

#' Logging Timestamp
#' @description
#' Internal helper for .log_ts.
#' @return Character timestamp formatted for logs.
#' @keywords internal
.log_ts <- function() {
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

#' Normalize Step Id
#' @description
#' Internal helper to normalize step ids.
#' @param step Step id (e.g., "S01" or 1).
#' @return Step id in Sxx form.
#' @keywords internal
.log_step_id <- function(step) {
    if (is.null(step)) return("S??")
    step_chr <- as.character(step)
    if (grepl("^S[0-9]{2}$", step_chr)) return(step_chr)
    if (grepl("^[0-9]+$", step_chr)) return(sprintf("S%02d", as.integer(step_chr)))
    step_chr
}

#' Logging Prefix
#' @description
#' Internal helper for .log_prefix.
#' @param parent Parent function name.
#' @param step Step id.
#' @return String prefix for structured logs.
#' @keywords internal
.log_prefix <- function(parent, step) {
    paste0("[", .log_ts(), "][", parent, "][", .log_step_id(step), "]")
}

.genescope_log_state <- local({
    env <- new.env(parent = emptyenv())
    env$stack <- list()
    env
})

#' Sanitize Runtime Log Message
#' @description
#' Redacts unstable absolute prefixes from user-facing logs.
#' @param msg Message string.
#' @return Sanitized scalar character string.
#' @keywords internal
.sanitize_log_message <- function(msg) {
    out <- paste(as.character(msg), collapse = "")
    temp_prefix <- normalizePath(tempdir(), winslash = "/", mustWork = FALSE)
    home_prefix <- normalizePath(path.expand("~"), winslash = "/", mustWork = FALSE)
    prefix_pattern <- function(path) {
        if (is.na(path) || !nzchar(path)) return(NULL)
        escaped <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", path, perl = TRUE)
        gsub("/", "/+", escaped, fixed = TRUE)
    }
    if (!is.na(temp_prefix) && nzchar(temp_prefix)) {
        out <- gsub(prefix_pattern(temp_prefix), "<tmp>", out, perl = TRUE)
    }
    out <- gsub("/\\S*Rtmp\\S*", "<tmp>", out, perl = TRUE)
    if (!is.na(home_prefix) && nzchar(home_prefix)) {
        out <- gsub(prefix_pattern(home_prefix), "~", out, perl = TRUE)
    }
    out
}

#' Get Active Log Session
#' @description
#' Returns the current buffered log session, if any.
#' @return Environment or `NULL`.
#' @keywords internal
.active_log_session <- function() {
    stack <- .genescope_log_state$stack
    if (!length(stack)) return(NULL)
    stack[[length(stack)]]
}

#' Begin Runtime Log Session
#' @description
#' Starts a log session used to track wrapper nesting.
#' @param parent Parent function label.
#' @return Session environment.
#' @keywords internal
.begin_log_session <- function(parent = NA_character_) {
    session <- new.env(parent = emptyenv())
    session$parent <- as.character(parent)[1]
    session$lines <- character()
    session$closed <- FALSE
    .genescope_log_state$stack <- c(.genescope_log_state$stack, list(session))
    session
}

#' Emit Runtime Log Line
#' @description
#' Emits a formatted runtime log line immediately.
#' @param line Pre-formatted log line.
#' @return Invisible NULL.
#' @keywords internal
.emit_log_line <- function(line) {
    line <- .sanitize_log_message(line)
    message(line)
    invisible(NULL)
}

#' Flush Buffered Log Session
#' @description
#' Flushes the current session to its parent session or to `message()`.
#' @param session Session environment.
#' @return Invisible NULL.
#' @keywords internal
.flush_log_session <- function(session = .active_log_session()) {
    if (is.null(session) || isTRUE(session$closed)) return(invisible(NULL))
    lines <- session$lines
    session$lines <- character()
    if (!length(lines)) return(invisible(NULL))

    stack <- .genescope_log_state$stack
    idx <- which(vapply(stack, identical, logical(1), session))
    if (length(idx) && idx[[1]] > 1L) {
        parent_session <- stack[[idx[[1]] - 1L]]
        parent_session$lines <- c(parent_session$lines, lines)
        return(invisible(NULL))
    }

    for (line in lines) {
        message(line)
    }
    invisible(NULL)
}

#' End Buffered Log Session
#' @description
#' Flushes and removes a buffered log session.
#' @param session Session environment.
#' @return Invisible NULL.
#' @keywords internal
.end_log_session <- function(session = .active_log_session()) {
    if (is.null(session) || isTRUE(session$closed)) return(invisible(NULL))
    .flush_log_session(session)
    stack <- .genescope_log_state$stack
    idx <- which(vapply(stack, identical, logical(1), session))
    if (length(idx)) {
        .genescope_log_state$stack <- stack[-idx[[1]]]
    }
    session$closed <- TRUE
    invisible(NULL)
}

#' Run Expression Within Runtime Log Session
#' @description
#' Tracks log session nesting while emitting runtime logs immediately.
#' @param parent Parent function label.
#' @param verbose Logical flag carried for interface symmetry.
#' @param expr Expression to evaluate.
#' @return Result of `expr`.
#' @keywords internal
.with_log_session <- function(parent, verbose = TRUE, expr) {
    session <- .begin_log_session(parent)
    on.exit(.end_log_session(session), add = TRUE)
    tryCatch(
        force(expr),
        error = function(e) {
            .log_key(parent, conditionMessage(e), level = "ERROR")
            stop(e)
        }
    )
}

#' Null-Coalescing Helper
#' @description
#' Internal helper returning the first non-null/non-empty scalar candidate.
#' @param ... Candidate values.
#' @return First usable scalar or `NULL`.
#' @keywords internal
.log_first_scalar <- function(...) {
    vals <- list(...)
    for (val in vals) {
        if (is.null(val) || !length(val)) next
        if (is.list(val) && !is.data.frame(val)) next
        scalar <- tryCatch(as.character(val)[1], error = function(e) NA_character_)
        if (!is.na(scalar) && nzchar(scalar)) return(scalar)
    }
    NULL
}

#' Log Key
#' @description
#' Emit a compact, user-facing key log that appears in both compact and
#' detailed modes.
#' @param parent Parent function name.
#' @param msg Message body.
#' @param level Log level label, usually `KEY`, `WARN`, or `ERROR`.
#' @return Invisible NULL.
#' @keywords internal
.log_key <- function(parent, msg, level = "KEY") {
    .emit_log_line(paste0(.log_prefix(parent, "KEY"), " ", level, ": ", msg))
    invisible(NULL)
}

#' Log Key-Value Summary
#' @description
#' Internal helper for compact summary lines.
#' @param parent Parent function name.
#' @param fields Named list of scalar summary fields.
#' @param level Log level label.
#' @return Invisible NULL.
#' @keywords internal
.log_key_values <- function(parent, fields, level = "KEY") {
    keep <- vapply(fields, function(x) {
        !(is.null(x) || !length(x) || (length(x) == 1L && is.na(x)))
    }, logical(1))
    fields <- fields[keep]
    msg <- if (length(fields)) {
        paste(
            paste0(names(fields), "=", vapply(fields, function(x) as.character(x)[1], character(1))),
            collapse = " "
        )
    } else {
        "no_summary_fields=TRUE"
    }
    .log_key(parent, msg, level = level)
}

#' Log Info
#' @description
#' Internal helper for .log_info.
#' @param parent Parent function name.
#' @param step Step id.
#' @param msg Message body.
#' @param verbose Logical; whether to emit progress messages.
#' @return Invisible NULL.
#' @keywords internal
.log_info <- function(parent, step, msg, verbose) {
    if (!.log_enabled(verbose)) return(invisible(NULL))
    if (.log_compact_enabled(verbose) &&
        .log_importance_score(parent, step, msg, kind = "info") < 4L) {
        return(invisible(NULL))
    }
    .emit_log_line(paste0(.log_prefix(parent, step), " INFO: ", msg))
    invisible(NULL)
}

#' Log Step
#' @description
#' Internal helper for .log_step.
#' @param parent Parent function name.
#' @param step Step id.
#' @param title Step title.
#' @param verbose Logical; whether to emit progress messages.
#' @return List with enter() and done() closures.
#' @keywords internal
.log_step <- function(parent, step, title, verbose) {
    start_time <- Sys.time()
    step_id <- .log_step_id(step)

    enter <- function(extra = NULL) {
        if (!.log_enabled(verbose)) return(invisible(NULL))
        if (.log_compact_enabled(verbose)) {
            if (!.log_compact_step_visible(parent, step_id) ||
                !.log_compact_step_extra_relevant(extra, phase = "enter")) {
                return(invisible(NULL))
            }
        }
        msg <- title
        if (!is.null(extra) && nzchar(extra)) msg <- paste0(msg, " ", extra)
        .emit_log_line(paste0(.log_prefix(parent, step_id), " ENTER: ", msg))
        invisible(NULL)
    }

    done <- function(extra = NULL) {
        if (!.log_enabled(verbose)) return(invisible(NULL))
        if (.log_compact_enabled(verbose)) {
            if (!.log_compact_step_visible(parent, step_id) ||
                !.log_compact_step_extra_relevant(extra, phase = "done")) {
                return(invisible(NULL))
            }
        }
        elapsed_ms <- as.integer(round(as.numeric(difftime(Sys.time(), start_time, units = "secs")) * 1000))
        msg <- paste0(title, " elapsed=", elapsed_ms, "ms")
        if (!is.null(extra) && nzchar(extra)) msg <- paste0(msg, " ", extra)
        .emit_log_line(paste0(.log_prefix(parent, step_id), " DONE: ", msg))
        invisible(NULL)
    }

    list(enter = enter, done = done, start_time = start_time, step_id = step_id)
}

#' Log Backend
#' @description
#' Internal helper for .log_backend.
#' @param parent Parent function name.
#' @param step Step id.
#' @param key Backend key name.
#' @param value Backend value.
#' @param reason Optional reason string.
#' @param verbose Logical; whether to emit progress messages.
#' @return Invisible NULL.
#' @keywords internal
.log_backend <- function(parent, step, key, value, reason = NULL, verbose) {
    if (!.log_enabled(verbose)) return(invisible(NULL))
    msg <- paste0(key, "=", value)
    if (!is.null(reason) && nzchar(reason)) msg <- paste0(msg, " reason=", reason)
    if (.log_compact_enabled(verbose) &&
        .log_importance_score(parent, step, msg, kind = "backend") < 4L &&
        !.log_compact_step_extra_relevant(msg, phase = "done")) {
        return(invisible(NULL))
    }
    .emit_log_line(paste0(.log_prefix(parent, step), " BACKEND: ", msg))
    invisible(NULL)
}

#' Log Function Start
#' @description
#' Emit a START marker for function entry. Visible in both compact and debug modes when verbose=TRUE.
#' @param parent Parent function name.
#' @param msg Optional additional message.
#' @param verbose Logical; whether to emit the log. Defaults to getOption("geneSCOPE.verbose", TRUE).
#' @return Invisible NULL.
#' @keywords internal
.log_start <- function(parent, msg = NULL, verbose = getOption("geneSCOPE.verbose", TRUE)) {
    if (!isTRUE(verbose)) return(invisible(NULL))
    full_msg <- if (!is.null(msg) && nzchar(msg)) paste0("START: ", msg) else "START"
    .emit_log_line(paste0(.log_prefix(parent, "START"), " ", full_msg))
    invisible(NULL)
}

#' Log Function Done
#' @description
#' Emit a DONE marker for function exit. Visible in both compact and debug modes when verbose=TRUE.
#' @param parent Parent function name.
#' @param msg Optional additional message.
#' @param verbose Logical; whether to emit the log. Defaults to getOption("geneSCOPE.verbose", TRUE).
#' @return Invisible NULL.
#' @keywords internal
.log_done <- function(parent, msg = NULL, verbose = getOption("geneSCOPE.verbose", TRUE)) {
    if (!isTRUE(verbose)) return(invisible(NULL))
    full_msg <- if (!is.null(msg) && nzchar(msg)) paste0("DONE: ", msg) else "DONE"
    .emit_log_line(paste0(.log_prefix(parent, "DONE"), " ", full_msg))
    invisible(NULL)
}

#' Log Warning
#' @description
#' Convenience wrapper for emitting WARN level logs.
#' @param parent Parent function name.
#' @param msg Warning message.
#' @return Invisible NULL.
#' @keywords internal
.log_warn <- function(parent, msg) {
    .log_key(parent, msg, level = "WARN")
}

#' Log Error
#' @description
#' Convenience wrapper for emitting ERROR level logs.
#' @param parent Parent function name.
#' @param msg Error message.
#' @return Invisible NULL.
#' @keywords internal
.log_error <- function(parent, msg) {
    .log_key(parent, msg, level = "ERROR")
}
