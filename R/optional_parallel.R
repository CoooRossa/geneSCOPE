#' @title Parallel Processing Optional Module
#' @description Adapter layer for parallel computing with future and future.apply.
#' @name optional-parallel
#' @keywords internal
NULL

#' Check if future is available
#' @return Logical
#' @keywords internal
.future_available <- function() {
  requireNamespace("future", quietly = TRUE)
}

#' Check if future.apply is available
#' @return Logical
#' @keywords internal
.future_apply_available <- function() {
  requireNamespace("future.apply", quietly = TRUE)
}

#' Require future or stop
#' @return Invisible TRUE
#' @keywords internal
.require_future_or_stop <- function() {
  if (!.future_available()) {
    stop(
      "Parallel processing requires future package.\n",
      "Install with:\n",
      "  install.packages('future')"
    )
  }
  invisible(TRUE)
}

#' Require future.apply or stop
#' @return Invisible TRUE
#' @keywords internal
.require_future_apply_or_stop <- function() {
  if (!.future_apply_available()) {
    stop(
      "Parallel apply operations require future.apply package.\n",
      "Install with:\n",
      "  install.packages('future.apply')"
    )
  }
  invisible(TRUE)
}

#' Set parallel plan (adapter)
#' @description Wrapper around future::plan with sensible defaults
#' @param strategy Parallel strategy ("sequential", "multisession", "multicore", "cluster")
#' @param workers Number of workers (default NULL for auto)
#' @return Invisible NULL
#' @keywords internal
.set_parallel_plan <- function(strategy = c("sequential", "multisession", 
                                              "multicore", "cluster"),
                                workers = NULL) {
  .require_future_or_stop()
  strategy <- match.arg(strategy)
  
  if (is.null(workers)) {
    future::plan(strategy)
  } else {
    future::plan(strategy, workers = workers)
  }
  invisible(NULL)
}

#' Parallel lapply (adapter)
#' @description Wrapper around future.apply::future_lapply
#' @param X List or vector
#' @param FUN Function to apply
#' @param ... Additional arguments passed to FUN
#' @param future.scheduling Scheduling strategy
#' @return List of results
#' @keywords internal
.parallel_lapply <- function(X, FUN, ..., future.scheduling = 1.0) {
  .require_future_apply_or_stop()
  future.apply::future_lapply(X, FUN, ..., future.scheduling = future.scheduling)
}

#' Parallel sapply (adapter)
#' @description Wrapper around future.apply::future_sapply
#' @param X List or vector
#' @param FUN Function to apply
#' @param ... Additional arguments passed to FUN
#' @param simplify Simplify result (default TRUE)
#' @param future.scheduling Scheduling strategy
#' @return Vector or matrix of results
#' @keywords internal
.parallel_sapply <- function(X, FUN, ..., simplify = TRUE, 
                               future.scheduling = 1.0) {
  .require_future_apply_or_stop()
  future.apply::future_sapply(X, FUN, ..., simplify = simplify, 
                               future.scheduling = future.scheduling)
}

#' Check if running in parallel mode
#' @description Check if future is configured for parallel execution
#' @return Logical
#' @keywords internal
.is_parallel_mode <- function() {
  if (!.future_available()) return(FALSE)
  
  plan <- future::plan()
  !identical(class(plan), c("sequential", "strategy", "future"))
}
