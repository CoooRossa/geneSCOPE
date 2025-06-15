## ---- S4 class & methods ----
if (!methods::isClass("CoordObj")) {
  methods::setClass(
    "CoordObj",
    slots = c(
      coord     = "list",       # centroids, segmentation, etc.
      grid      = "list",       # one or more grid layers
      meta.data = "data.frame", # misc. metadata
      cells     = "list" 
    ))
}

if (is.null(methods::selectMethod("initialize", "CoordObj", optional = TRUE))) {
  methods::setMethod(
    "initialize", "CoordObj",
    function(.Object, coord = list(), grid = list(), meta.data = data.frame()) {
      .Object@coord     <- coord
      .Object@grid      <- grid
      .Object@meta.data <- meta.data
      .Object@cells     <- cells
      .Object
    }
  )
}

#' @title Upgrade `CoordObj` Cells Slot
#' @description Converts the `cells` slot from a sparse matrix to a list format.
#' @param obj A `CoordObj` object with the `cells` slot as a sparse matrix.
#' @return The modified `CoordObj` with the `cells` slot upgraded to a list format.
#' @export
upgradeCellsSlot <- function(obj) {
  if (!inherits(obj@cells, "dgCMatrix")) return(obj)  # 已是新格式

  lst <- list(counts = obj@cells)                    # 原矩阵命名为 counts
  atts <- attributes(obj@cells)

  # 把 attr 中附带的矩阵也转移进列表
  extra <- atts[setdiff(names(atts), c("dim", "dimnames", "i", "p", "x", "factors"))]
  extra <- extra[vapply(extra, inherits, logical(1), "dgCMatrix")]
  lst <- c(lst, extra)

  obj@cells <- lst
  obj
}


#' @useDynLib FG2CLI, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

