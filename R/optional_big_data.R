#' @title Big Data Processing Optional Module
#' @description Adapter layer for large-scale data processing.
#'   Requires bigmemory, hdf5r, or rhdf5 packages.
#' @name optional-big-data
#' @keywords internal
NULL

#' Check if bigmemory is available
#' @return Logical
#' @keywords internal
.bigmemory_available <- function() {
  requireNamespace("bigmemory", quietly = TRUE)
}

#' Check if hdf5r is available
#' @return Logical
#' @keywords internal
.hdf5r_available <- function() {
  requireNamespace("hdf5r", quietly = TRUE)
}

#' Check if rhdf5 is available
#' @return Logical
#' @keywords internal
.rhdf5_available <- function() {
  requireNamespace("rhdf5", quietly = TRUE)
}

#' Require bigmemory or stop
#' @return Invisible TRUE
#' @keywords internal
.require_bigmemory_or_stop <- function() {
  if (!.bigmemory_available()) {
    stop(
      "Big matrix operations require bigmemory package.\n",
      "Install with:\n",
      "  install.packages('bigmemory')"
    )
  }
  invisible(TRUE)
}

#' Require HDF5 support or stop
#' @return Invisible TRUE
#' @keywords internal
.require_hdf5_or_stop <- function() {
  if (!.hdf5r_available() && !.rhdf5_available()) {
    stop(
      "HDF5 file operations require hdf5r or rhdf5 package.\n",
      "Install with:\n",
      "  install.packages('hdf5r')\n",
      "  # or\n",
      "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
      "  BiocManager::install('rhdf5')"
    )
  }
  invisible(TRUE)
}

#' Create big.matrix (adapter)
#' @description Wrapper around bigmemory::big.matrix
#' @param nrow Number of rows
#' @param ncol Number of columns
#' @param type Type of matrix (default "double")
#' @param backingfile Optional backing file path
#' @param ... Additional arguments
#' @return big.matrix object
#' @keywords internal
.create_big_matrix <- function(nrow, ncol, type = "double", 
                                 backingfile = NULL, ...) {
  .require_bigmemory_or_stop()
  
  if (is.null(backingfile)) {
    bigmemory::big.matrix(nrow = nrow, ncol = ncol, type = type, ...)
  } else {
    bigmemory::filebacked.big.matrix(
      nrow = nrow, ncol = ncol, type = type,
      backingfile = backingfile, ...
    )
  }
}

#' Convert big.matrix to regular matrix (adapter)
#' @description Wrapper around bigmemory::as.matrix
#' @param bigmat big.matrix object
#' @return Regular R matrix
#' @keywords internal
.big_matrix_to_matrix <- function(bigmat) {
  if (!.bigmemory_available()) {
    stop("bigmemory package required")
  }
  bigmemory::as.matrix(bigmat)
}

#' Read HDF5 file (adapter)
#' @description Read data from HDF5 file using available package
#' @param file Path to HDF5 file
#' @param name Dataset name within file
#' @return Matrix or array
#' @keywords internal
.read_hdf5 <- function(file, name) {
  .require_hdf5_or_stop()
  
  if (.hdf5r_available()) {
    h5 <- hdf5r::H5File$new(file, mode = "r")
    on.exit(h5$close(), add = TRUE)
    h5[[name]]$read()
  } else if (.rhdf5_available()) {
    rhdf5::h5read(file, name)
  } else {
    stop("No HDF5 package available")
  }
}

#' Write HDF5 file (adapter)
#' @description Write data to HDF5 file using available package
#' @param data Data to write
#' @param file Path to HDF5 file
#' @param name Dataset name within file
#' @return Invisible NULL
#' @keywords internal
.write_hdf5 <- function(data, file, name) {
  .require_hdf5_or_stop()
  
  if (.hdf5r_available()) {
    h5 <- hdf5r::H5File$new(file, mode = "a")
    on.exit(h5$close(), add = TRUE)
    h5[[name]] <- data
  } else if (.rhdf5_available()) {
    rhdf5::h5write(data, file, name)
  }
  invisible(NULL)
}
