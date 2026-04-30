#' @title Imaging Optional Module
#' @description Adapter layer for image processing (JPEG, PNG).
#' @name optional-imaging
#' @keywords internal
NULL

#' Check if jpeg is available
#' @return Logical
#' @keywords internal
.jpeg_available <- function() {
  requireNamespace("jpeg", quietly = TRUE)
}

#' Check if png is available
#' @return Logical
#' @keywords internal
.png_available <- function() {
  requireNamespace("png", quietly = TRUE)
}

#' Require jpeg or stop
#' @return Invisible TRUE
#' @keywords internal
.require_jpeg_or_stop <- function() {
  if (!.jpeg_available()) {
    stop(
      "JPEG image processing requires jpeg package.\n",
      "Install with:\n",
      "  install.packages('jpeg')"
    )
  }
  invisible(TRUE)
}

#' Require png or stop
#' @return Invisible TRUE
#' @keywords internal
.require_png_or_stop <- function() {
  if (!.png_available()) {
    stop(
      "PNG image processing requires png package.\n",
      "Install with:\n",
      "  install.packages('png')"
    )
  }
  invisible(TRUE)
}

#' Read JPEG image (adapter)
#' @description Wrapper around jpeg::readJPEG
#' @param source Path to JPEG file
#' @param native Whether to return native raster (default FALSE)
#' @return Raster array
#' @keywords internal
.read_jpeg <- function(source, native = FALSE) {
  .require_jpeg_or_stop()
  jpeg::readJPEG(source, native = native)
}

#' Write JPEG image (adapter)
#' @description Wrapper around jpeg::writeJPEG
#' @param image Image data array
#' @param target Output file path
#' @param quality JPEG quality (1-100)
#' @return Invisible NULL
#' @keywords internal
.write_jpeg <- function(image, target, quality = 75) {
  .require_jpeg_or_stop()
  jpeg::writeJPEG(image, target = target, quality = quality)
  invisible(NULL)
}

#' Read PNG image (adapter)
#' @description Wrapper around png::readPNG
#' @param source Path to PNG file
#' @param native Whether to return native raster (default FALSE)
#' @param info Whether to return image info (default FALSE)
#' @return Raster array
#' @keywords internal
.read_png <- function(source, native = FALSE, info = FALSE) {
  .require_png_or_stop()
  png::readPNG(source, native = native, info = info)
}

#' Write PNG image (adapter)
#' @description Wrapper around png::writePNG
#' @param image Image data array
#' @param target Output file path (default NULL for console)
#' @return Invisible NULL
#' @keywords internal
.write_png <- function(image, target = NULL) {
  .require_png_or_stop()
  png::writePNG(image, target = target)
  invisible(NULL)
}
