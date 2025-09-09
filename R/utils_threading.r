#' Threading and BLAS Management Utilities
#' 
#' @description
#' Utilities to safely manage threading in the presence of both OpenMP 
#' and BLAS libraries, preventing conflicts that cause OpenBLAS warnings.

#' @title Configure threads for different operation types
#' @description
#' Intelligently configures OpenMP and BLAS threading based on operation type
#' and system capabilities, with automatic restoration.
#' @param operation_type Character. One of:
#'   - "compute_intensive": Heavy computation (Lee's L, correlations)
#'   - "io_bound": File I/O and data processing  
#'   - "mixed": Mixed workloads
#'   - "openmp_only": Pure OpenMP, disable BLAS threading
#' @param ncores_requested Integer. Requested number of cores
#' @param restore_after Logical. If TRUE, return restoration function
#' @return List with thread configuration and optional restore function
#' @export
configureThreadsFor <- function(operation_type = c("compute_intensive", "io_bound", "mixed", "openmp_only"),
                                ncores_requested = 1,
                                restore_after = FALSE) {
  operation_type <- match.arg(operation_type)
  
  # 保存当前状态
  old_state <- .saveThreadState()
  
  # 强制BLAS单线程以避免OpenBLAS警告
  .forceBLASSingleThread()
  
  # 根据操作类型配置
  config <- switch(operation_type,
    compute_intensive = .configureComputeIntensive(ncores_requested),
    io_bound = .configureIOBound(ncores_requested), 
    mixed = .configureMixed(ncores_requested),
    openmp_only = .configureOpenMPOnly(ncores_requested)
  )
  
  # 应用配置
  .applyThreadConfig(config)
  
  # 准备返回值
  result <- list(
    operation_type = operation_type,
    openmp_threads = config$openmp_threads,
    r_threads = config$r_threads,
    blas_threads = 1  # 始终为1
  )
  
  if (restore_after) {
    restore_function <- function() {
      .restoreThreadState(old_state)
    }
    attr(result, "restore_function") <- restore_function
  }
  
  return(result)
}

#' @title Force BLAS to single thread
#' @description Aggressively sets all known BLAS libraries to single thread
.forceBLASSingleThread <- function() {
  # 1. RhpcBLASctl方法
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    tryCatch({
      RhpcBLASctl::blas_set_num_threads(1)
    }, error = function(e) NULL)
  }
  
  # 2. 环境变量方法（更可靠）
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
  
  # 3. 设置data.table线程为1
  if (requireNamespace("data.table", quietly = TRUE)) {
    tryCatch({
      data.table::setDTthreads(1)
    }, error = function(e) NULL)
  }
}

#' @title Save current thread state
.saveThreadState <- function() {
  list(
    omp_num_threads = Sys.getenv("OMP_NUM_THREADS", unset = NA),
    openblas_threads = Sys.getenv("OPENBLAS_NUM_THREADS", unset = NA),
    mkl_threads = Sys.getenv("MKL_NUM_THREADS", unset = NA),
    blas_threads = if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      tryCatch(RhpcBLASctl::blas_get_num_procs(), error = function(e) 1)
    } else {
      1
    }
  )
}

#' @title Restore thread state
.restoreThreadState <- function(old_state) {
  # 恢复环境变量
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
  
  # 但是保持BLAS为单线程，防止冲突
  .forceBLASSingleThread()
}

#' @title Configure for compute-intensive operations
.configureComputeIntensive <- function(ncores_requested) {
  max_cores <- parallel::detectCores()
  safe_cores <- min(ncores_requested, max_cores - 1, 16)  # 保守限制
  
  list(
    openmp_threads = safe_cores,
    r_threads = 1,  # R本身单线程
    blas_threads = 1  # BLAS强制单线程
  )
}

#' @title Configure for I/O bound operations  
.configureIOBound <- function(ncores_requested) {
  max_cores <- parallel::detectCores()
  safe_cores <- min(ncores_requested, max_cores - 1, 8)
  
  list(
    openmp_threads = 2,  # 少量OpenMP
    r_threads = safe_cores,  # 更多R线程用于I/O
    blas_threads = 1
  )
}

#' @title Configure for mixed operations
.configureMixed <- function(ncores_requested) {
  max_cores <- parallel::detectCores()
  safe_cores <- min(ncores_requested, max_cores - 1, 12)
  
  list(
    openmp_threads = max(1, safe_cores %/% 2),
    r_threads = max(1, safe_cores %/% 2), 
    blas_threads = 1
  )
}

#' @title Configure for OpenMP-only operations
.configureOpenMPOnly <- function(ncores_requested) {
  max_cores <- parallel::detectCores()
  safe_cores <- min(ncores_requested, max_cores - 1, 16)
  
  list(
    openmp_threads = safe_cores,
    r_threads = 1,
    blas_threads = 1  # 依然强制单线程
  )
}

#' @title Apply thread configuration
.applyThreadConfig <- function(config) {
  # 设置OpenMP线程
  Sys.setenv(OMP_NUM_THREADS = as.character(config$openmp_threads))
  
  # 确保BLAS保持单线程
  .forceBLASSingleThread()
  
  # 设置Arrow线程（如果可用）
  if (requireNamespace("arrow", quietly = TRUE)) {
    tryCatch({
      if (exists("set_cpu_count", where = asNamespace("arrow"))) {
        arrow::set_cpu_count(config$r_threads)
      }
      if (exists("set_io_thread_count", where = asNamespace("arrow"))) {
        arrow::set_io_thread_count(config$r_threads) 
      }
    }, error = function(e) NULL)
  }
}

#' @title Detect operating system
#' @description Simple OS detection utility
#' @return Character string: "windows", "macos", or "linux"
#' @export
detectOS <- function() {
  if (.Platform$OS.type == "windows") {
    return("windows")
  } else if (Sys.info()["sysname"] == "Darwin") {
    return("macos")
  } else {
    return("linux")
  }
}

#' @title Get safe thread count
#' @description Returns a conservative thread count that should work on most systems
#' @param max_requested Maximum threads requested
#' @return Safe number of threads to use
#' @export  
getSafeThreadCount <- function(max_requested = parallel::detectCores()) {
  total_cores <- parallel::detectCores()
  os_type <- detectOS()
  
  # 系统特定的安全限制
  safe_limit <- switch(os_type,
    windows = min(8, total_cores - 1),
    macos = min(12, total_cores - 1), 
    linux = min(16, total_cores - 1)
  )
  
  return(min(max_requested, safe_limit, na.rm = TRUE))
}
