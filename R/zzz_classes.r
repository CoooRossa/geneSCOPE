## ---- S4 class & methods ----
if (!methods::isClass("CoordObj")) {
  methods::setClass("CoordObj",
                    slots = c(
                      coord     = "list",
                      grid      = "list",
                      meta.data = "data.frame"
                    ))
}

if (is.null(methods::selectMethod("initialize", "CoordObj", optional = TRUE))) {
  methods::setMethod(
    "initialize", "CoordObj",
    function(.Object, coord = list(), grid = list(), meta.data = data.frame()) {
      .Object@coord     <- coord
      .Object@grid      <- grid
      .Object@meta.data <- meta.data
      .Object
    }
  )
}