#' @title Python Bridge Optional Module
#' @description Adapter layer for Python integration via reticulate.
#'   Required for Xenium transcript fallback conversion and pyarrow support.
#' @name optional-python-bridge
#' @keywords internal
NULL

#' Check if reticulate is available
#' @return Logical
#' @keywords internal
.reticulate_available <- function() {
  requireNamespace("reticulate", quietly = TRUE)
}

#' Check if pyarrow is available in Python environment
#' @return Logical
#' @keywords internal
.pyarrow_available <- function() {
  if (!.reticulate_available()) return(FALSE)
  tryCatch({
    reticulate::py_module_available("pyarrow")
  }, error = function(e) FALSE)
}

#' Require reticulate or stop
#' @return Invisible TRUE
#' @keywords internal
.require_reticulate_or_stop <- function() {
  if (!.reticulate_available()) {
    stop(
      "Python integration requires reticulate package.\n",
      "Install with:\n",
      "  install.packages('reticulate')\n",
      "Then set up Python environment:\n",
      "  reticulate::install_miniconda()\n",
      "  reticulate::py_install('pyarrow')"
    )
  }
  invisible(TRUE)
}

#' Require pyarrow or stop
#' @return Invisible TRUE
#' @keywords internal
.require_pyarrow_or_stop <- function() {
  .require_reticulate_or_stop()
  
  if (!.pyarrow_available()) {
    stop(
      "pyarrow Python module required for Xenium transcript conversion.\n",
      "Install with:\n",
      "  reticulate::py_install('pyarrow')"
    )
  }
  invisible(TRUE)
}

#' Import Python module (adapter)
#' @description Wrapper around reticulate::import with error handling
#' @param module Module name
#' @param ... Additional arguments
#' @return Python module object
#' @keywords internal
.import_python <- function(module, ...) {
  .require_reticulate_or_stop()
  reticulate::import(module, ...)
}

#' Check if Python module is available (adapter)
#' @description Wrapper around reticulate::py_module_available
#' @param module Module name
#' @return Logical
#' @keywords internal
.py_module_available <- function(module) {
  if (!.reticulate_available()) return(FALSE)
  reticulate::py_module_available(module)
}

#' Use specific Python environment (adapter)
#' @description Set Python environment for reticulate
#' @param env_type Environment type ("conda", "virtualenv", "python")
#' @param env_name Environment name or path
#' @return Invisible NULL
#' @keywords internal
.use_python_env <- function(env_type = c("conda", "virtualenv", "python"), 
                             env_name) {
  .require_reticulate_or_stop()
  env_type <- match.arg(env_type)
  
  if (env_type == "conda") {
    reticulate::use_condaenv(env_name, required = TRUE)
  } else if (env_type == "virtualenv") {
    reticulate::use_virtualenv(env_name, required = TRUE)
  } else {
    reticulate::use_python(env_name, required = TRUE)
  }
  invisible(NULL)
}

#' Convert Python object to R (adapter)
#' @description Wrapper around reticulate::py_to_r
#' @param py_obj Python object
#' @return R object
#' @keywords internal
.py_to_r <- function(py_obj) {
  if (!.reticulate_available()) {
    stop("reticulate package required")
  }
  reticulate::py_to_r(py_obj)
}

#' Convert R object to Python (adapter)
#' @description Wrapper around reticulate::r_to_py
#' @param r_obj R object
#' @param convert Whether to convert on access
#' @return Python object
#' @keywords internal
.r_to_py <- function(r_obj, convert = TRUE) {
  if (!.reticulate_available()) {
    stop("reticulate package required")
  }
  reticulate::r_to_py(r_obj, convert = convert)
}
