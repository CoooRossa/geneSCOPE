## ---- S4 class & methods ----
if (!methods::isClass("scope_object")) {
  methods::setClass(
    "scope_object",
    slots = c(
      coord     = "list", # centroids, segmentation, …
      grid      = "list", # grid layers
      meta.data = "data.frame", # miscellaneous gene / grid meta
      cells     = "list",
      stats     = "list", # statistics layers
      density   = "list" # density layers
    )
  )
}

if (is.null(methods::selectMethod("initialize", "scope_object", optional = TRUE))) {
  methods::setMethod(
    "initialize", "scope_object",
    function(.Object,
             coord = list(),
             grid = list(),
             meta.data = data.frame(),
             cells = list(),
             stats = list(),
             density = list()) {
      ## slot assignment ----
      .Object@coord <- coord
      .Object@grid <- grid
      .Object@meta.data <- meta.data
      .Object@cells <- cells
      .Object@stats <- stats
      .Object@density <- density

      ## basic validation – avoid downstream function errors ----
      if (nrow(.Object@meta.data) > 0 && is.null(rownames(.Object@meta.data))) {
        warning(
          "meta.data has no row-names; downstream functions assume gene names. ",
          "Consider setting row-names now."
        )
      }
      .Object
    }
  )

  ## concise object printout (PDF suggestion) ----
  methods::setMethod(
    "show", "scope_object",
    function(object) {
      cat("scope_object with:\n")
      cat("  coord:", length(object@coord), "elements\n")
      cat("  grid:", length(object@grid), "layers\n")
      cat("  cells:", length(object@cells), "matrices\n")
      cat("  stats:", length(object@stats), "result sets\n")
      cat("  density:", length(object@density), "density tables\n")
      cat("  meta.data:", nrow(object@meta.data), "genes x", ncol(object@meta.data), "features\n")
    }
  )

  ## summary() – quick overview ----
  methods::setMethod(
    "summary", signature(object = "scope_object"),
    function(object, ...) {
      cat("=== scope_object Summary ===\n\n")

      # Coordinate data
      if (length(object@coord) > 0) {
        cat("Coordinate data:\n")
        for (nm in names(object@coord)) {
          if (is.data.frame(object@coord[[nm]])) {
            cat(" ", nm, ":", nrow(object@coord[[nm]]), "rows\n")
          } else {
            cat(" ", nm, ":", length(object@coord[[nm]]), "elements\n")
          }
        }
        cat("\n")
      }

      # Grid layers
      if (length(object@grid) > 0) {
        cat("Grid layers:\n")
        for (nm in names(object@grid)) {
          g <- object@grid[[nm]]
          if ("grid_info" %in% names(g)) {
            cat(" ", nm, ":", nrow(g$grid_info), "grid cells\n")
          }
        }
        cat("\n")
      }

      # Cell data
      if (length(object@cells) > 0) {
        cat("Cell matrices:\n")
        for (nm in names(object@cells)) {
          mat <- object@cells[[nm]]
          if (is.matrix(mat) || inherits(mat, "Matrix")) {
            cat(" ", nm, ":", nrow(mat), "x", ncol(mat), "\n")
          }
        }
        cat("\n")
      }

      # Stats
      if (length(object@stats) > 0) {
        cat("Analysis results:\n")
        for (nm in names(object@stats)) {
          cat(" ", nm, ":", length(object@stats[[nm]]), "result sets\n")
        }
        cat("\n")
      }

      # Meta data
      if (nrow(object@meta.data) > 0) {
        cat("Gene metadata:", nrow(object@meta.data), "genes\n")
        if (ncol(object@meta.data) > 0) {
          cat("  Columns:", paste(colnames(object@meta.data), collapse = ", "), "\n")
        }
      }
    }
  )
}
