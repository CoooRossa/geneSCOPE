#' @title Spatial Geometry Optional Module
#' @description Adapter layer for spatial geometry operations using sf package.
#'   This module provides spatial predicates and geometry operations.
#' @name optional-spatial-geometry
#' @keywords internal
NULL

#' Check if sf is available
#' @return Logical
#' @keywords internal
.sf_available <- function() {
  requireNamespace("sf", quietly = TRUE)
}

#' Require sf or stop with informative error
#' @return Invisible TRUE if available
#' @keywords internal
.require_sf_or_stop <- function() {
  if (!.sf_available()) {
    stop(
      "Spatial geometry operations require sf package.\n",
      "Install with:\n",
      "  install.packages('sf')\n",
      "Note: sf requires system dependencies (GDAL, PROJ, GEOS)."
    )
  }
  invisible(TRUE)
}

#' Create sf object (adapter)
#' @description Wrapper around sf::st_as_sf
#' @param x Data frame or matrix with coordinates
#' @param ... Additional arguments passed to st_as_sf
#' @return sf object
#' @keywords internal
.create_sf_object <- function(x, ...) {
  .require_sf_or_stop()
  sf::st_as_sf(x, ...)
}

#' Create simple feature geometry collection (adapter)
#' @description Wrapper around sf::st_sfc
#' @param ... Geometry objects
#' @param crs Coordinate reference system
#' @return sfc object
#' @keywords internal
.create_geometry_collection <- function(..., crs = NA) {
  .require_sf_or_stop()
  sf::st_sfc(..., crs = crs)
}

#' Create polygon geometry (adapter)
#' @description Wrapper around sf::st_polygon
#' @param x List of matrices defining polygon rings
#' @return Polygon geometry
#' @keywords internal
.create_polygon <- function(x) {
  .require_sf_or_stop()
  sf::st_polygon(x)
}

#' Get geometry from sf object (adapter)
#' @description Wrapper around sf::st_geometry
#' @param x sf object
#' @return sfc geometry list column
#' @keywords internal
.get_geometry <- function(x) {
  .require_sf_or_stop()
  sf::st_geometry(x)
}

#' Get bounding box (adapter)
#' @description Wrapper around sf::st_bbox
#' @param obj sf object or extent
#' @return Named numeric vector with bounding box
#' @keywords internal
.get_bbox <- function(obj) {
  .require_sf_or_stop()
  sf::st_bbox(obj)
}

#' Get/set coordinate reference system (adapter)
#' @description Wrapper around sf::st_crs
#' @param x sf object
#' @param value CRS value to set (optional)
#' @return CRS object or modified sf object
#' @keywords internal
.get_or_set_crs <- function(x, value = NULL) {
  .require_sf_or_stop()
  if (is.null(value)) {
    sf::st_crs(x)
  } else {
    sf::st_crs(x) <- value
    x
  }
}

#' Union geometries (adapter)
#' @description Wrapper around sf::st_union
#' @param x sf object
#' @param y Optional second sf object
#' @return Unioned geometry
#' @keywords internal
.union_geometries <- function(x, y = NULL) {
  .require_sf_or_stop()
  if (is.null(y)) {
    sf::st_union(x)
  } else {
    sf::st_union(x, y)
  }
}

#' Spatial predicate: within (adapter)
#' @description Test if geometries are within other geometries
#' @param x sf object
#' @param y sf object
#' @return Logical matrix
#' @keywords internal
.spatial_within <- function(x, y) {
  .require_sf_or_stop()
  sf::st_within(x, y)
}

#' Make geometry valid (adapter)
#' @description Wrapper around sf::st_make_valid
#' @param x sf object
#' @return sf object with valid geometries
#' @keywords internal
.make_geometry_valid <- function(x) {
  .require_sf_or_stop()
  sf::st_make_valid(x)
}

#' Convert to sf object (adapter)
#' @description Generic conversion to sf
#' @param x Object to convert
#' @param ... Additional arguments
#' @return sf object
#' @keywords internal
.as_sf <- function(x, ...) {
  .require_sf_or_stop()
  sf::st_as_sf(x, ...)
}
