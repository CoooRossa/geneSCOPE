#' Attach histology assets to a scope object.
#' @description
#' Adds Visium-style histology images and scalefactors to an existing grid layer
#' within a `scope_object`.
#' @param scope_obj A `scope_object` to annotate.
#' @param grid_name Grid layer name to update.
#' @param png_path Path to the histology PNG.
#' @param json_path Optional Visium scalefactors JSON; used when `scalefactors` is not provided.
#' @param level Image resolution level (`lowres` or `hires`).
#' @param crop_bbox_px Optional pixel-space bounding box to crop the PNG before attaching.
#' @param roi_bbox Optional ROI bounding box metadata stored alongside the image.
#' @param coord_type Coordinate system for scaling (`visium` or `manual`).
#' @param scalefactors Optional list of scalefactors to use instead of the JSON file.
#' @param y_origin Coordinate origin policy (`auto`, `top-left`, `bottom-left`).
#' @return Updated `scope_object` with histology metadata.
#' @examples
#' \dontrun{
#' scope_obj <- createSCOPE(data_dir = "/path/to/output_dir")
#' scope_obj <- addScopeHistology(scope_obj, grid_name = "grid55", png_path = "/path/to/spatial/tissue_hires_image.png", json_path = "/path/to/spatial/scalefactors_json.json")
#' }
#' @seealso `plotGridBoundary()`, `createSCOPE_visium()`
#' @export
addScopeHistology <- function(
    scope_obj,
    grid_name,
    png_path,
    json_path = NULL,
    level = c("lowres", "hires"),
    crop_bbox_px = NULL,
    roi_bbox = NULL,
    coord_type = c("visium", "manual"),
    scalefactors = NULL,
    y_origin = c("auto", "top-left", "bottom-left")) {
    if (missing(level)) {
        level <- NULL
    }
    .add_scope_histology(
        scope_obj = scope_obj,
        grid_name = grid_name,
        png_path = png_path,
        json_path = json_path,
        level = level,
        crop_bbox_px = crop_bbox_px,
        roi_bbox = roi_bbox,
        coord_type = coord_type,
        scalefactors = scalefactors,
        y_origin = y_origin
    )
}

#' Add single-cell spatial data.
#' @description
#' Attaches platform-specific single-cell layers (CosMx or Xenium) to an existing
#' `scope_object`.
#' @param scope_obj A `scope_object` to extend.
#' @param data_dir Optional data directory (legacy alias for `xenium_dir` / `cosmx_root`).
#' @param xenium_dir Optional Xenium run directory.
#' @param cosmx_root Optional CosMx project root.
#' @param platform Deprecated platform hint (`auto`, `Xenium`, or `CosMx`).
#'   Use `prefer=` instead. Explicit use emits a deprecation warning and maps
#'   to `prefer=` when unambiguous.
#' @param prefer Preferred single-cell loader (`auto`, `xenium`, or `cosmx`).
#' @param verbose Logical; whether to emit compact public progress messages.
#'   Default is \code{TRUE}. Set to \code{FALSE} for silent operation, or use
#'   \code{options(geneSCOPE.verbose = FALSE)} for global control.
#' @param ... Additional arguments passed to the platform-specific loader.
#' @return Updated `scope_object` with single-cell layers.
#' @examples
#' \dontrun{
#' scope_obj <- createSCOPE(data_dir = "/path/to/output_dir", grid_length = 30)
#' scope_obj <- addSingleCells(scope_obj, data_dir = "/path/to/xenium_run")
#' }
#' @seealso `createSCOPE_xenium()`, `createSCOPE_cosmx()`
#' @export
addSingleCells <- function(
    scope_obj,
    data_dir = NULL,
    xenium_dir = NULL,
    cosmx_root = NULL,
    platform = c("auto", "Xenium", "CosMx"),
    prefer = c("auto", "xenium", "cosmx"),
    verbose = TRUE,
    ...) {
    verbose <- .resolve_runtime_verbose(verbose)
    platform_missing <- missing(platform)
    prefer_missing <- missing(prefer)
    prefer <- match.arg(prefer)
    if (!platform_missing) {
        platform_value <- match.arg(platform)
        platform_prefer <- switch(tolower(platform_value),
            auto = "auto",
            xenium = "xenium",
            cosmx = "cosmx"
        )
        .Deprecated(
            "prefer",
            package = "geneSCOPE",
            msg = "addSingleCells(platform=) is deprecated; use addSingleCells(prefer=) instead."
        )
        if (prefer_missing) {
            prefer <- platform_prefer
        } else if (!identical(platform_prefer, prefer)) {
            stop("addSingleCells(platform=) conflicts with prefer=. Use only prefer= to select the single-cell loader.", call. = FALSE)
        }
    }
    
    # Handle data_dir parameter (legacy alias)
    if (!is.null(data_dir)) {
        if (is.null(xenium_dir) && is.null(cosmx_root)) {
            guess <- .guess_singlecell_platform(data_dir)
            if (!is.null(guess$platform)) {
                if (identical(guess$platform, "Xenium")) {
                    xenium_dir <- guess$path
                } else if (identical(guess$platform, "CosMx")) {
                    cosmx_root <- guess$path
                }
            } else {
                # Fall back to prefer parameter
                if (identical(prefer, "cosmx")) {
                    cosmx_root <- data_dir
                } else if (identical(prefer, "xenium")) {
                    xenium_dir <- data_dir
                }
            }
        }
    }
    
    .add_single_cells(
        scope_obj = scope_obj,
        xenium_dir = xenium_dir,
        cosmx_root = cosmx_root,
        prefer = prefer,
        verbose = verbose,
        ...
    )
}

#' Create a scope object from platform outputs.
#' @description
#' Auto-detects Xenium/CosMx/Visium inputs under `data_dir` (or uses an explicit
#' `prefer` override) and dispatches to the matching constructor. The Visium
#' constructor uses the Seurat spatial builder with SCTransform by default.
#' Non-SCT Visium construction is no longer a supported public path and fails
#' fast with an explicit error. HDF5-backed Visium inputs require `hdf5r`
#' before `Seurat::Load10X_Spatial()` is called.
#' @param data_dir Root directory containing Xenium/CosMx/Visium outputs.
#' @param prefer Override platform detection (`auto`, `xenium`, `cosmx`, `visium`).
#' @param grid_length Grid size in microns for Xenium/CosMx (required for these platforms).
#' @param seg_type Segmentation type for Xenium (`cell`, `nucleus`, or `both`).
#' @param roi_path Optional path to ROI coordinate file.
#' @param num_workers Number of parallel workers for data processing.
#' @param filter_genes Optional character vector of genes to filter out.
#' @param max_distance_to_nucleus Maximum distance (in microns) from molecule to nucleus.
#' @param filter_molecules Whether to apply molecule quality filtering.
#' @param exclude_prefixes Character vector of gene prefixes to exclude.
#' @param min_quality Minimum quality score threshold for molecules.
#' @param min_segmentation_points Minimum points per segmentation unit.
#' @param parallel_backend Parallel backend choice (`auto`, `fork`, `psock`, `serial`).
#' @param verbose Logical; whether to emit compact public progress messages.
#'   Default is \code{TRUE}. Set to \code{FALSE} for silent operation, or use
#'   \code{options(geneSCOPE.verbose = FALSE)} for global control.
#' @param sctransform Run SCTransform when using the Visium Seurat builder.
#'   Must be `TRUE` for Visium inputs in this release.
#' @param ... Additional arguments forwarded to the platform constructor.
#' @return A constructed `scope_object`.
#' @examples
#' \dontrun{
#' scope_obj <- createSCOPE("/path/to/output_dir", grid_length = 30)
#' scope_obj <- createSCOPE("/path/to/visium_run")
#' }
#' @seealso `createSCOPE_visium()`, `createSCOPE_xenium()`, `createSCOPE_cosmx()`
#' @export
createSCOPE <- function(
    data_dir = NULL,
    prefer = c("auto", "xenium", "cosmx", "visium"),
    grid_length = NULL,
    seg_type = c("cell", "nucleus", "both"),
    roi_path = NULL,
    num_workers = 1L,
    filter_genes = NULL,
    max_distance_to_nucleus = 25,
    filter_molecules = TRUE,
    exclude_prefixes = c("Unassigned", "NegControl", "Background", "DeprecatedCodeword",
                         "SystemControl", "Negative"),
    min_quality = 20,
    min_segmentation_points = 0,
    parallel_backend = c("auto", "fork", "psock", "serial"),
    verbose = TRUE,
    sctransform = TRUE,
    ...) {
    dots <- list(...)
    prefer_missing <- missing(prefer)
    seg_type_missing <- missing(seg_type)
    parallel_backend_missing <- missing(parallel_backend)
    
    # Handle deprecated platform= parameter
    if (!is.null(dots$platform)) {
        platform <- tolower(as.character(dots$platform)[1])
        if (!platform %in% c("auto", "xenium", "cosmx", "visium")) {
            stop("createSCOPE(platform=) must be one of 'auto', 'xenium', 'cosmx', or 'visium'. Use prefer= instead.", call. = FALSE)
        }
        .Deprecated(
            "prefer",
            package = "geneSCOPE",
            msg = "createSCOPE(platform=) is deprecated and no longer silently ignored; use createSCOPE(prefer=) instead."
        )
        if (prefer_missing) {
            prefer <- platform
        } else {
            prefer_requested <- match.arg(prefer)
            if (!identical(platform, prefer_requested)) {
                stop("createSCOPE(platform=) conflicts with prefer=. Use only prefer= to select the platform.", call. = FALSE)
            }
        }
        dots$platform <- NULL
    }
    
    prefer <- match.arg(prefer)
    if (!seg_type_missing) {
        seg_type <- match.arg(seg_type)
    } else {
        seg_type <- "cell"
    }
    if (!parallel_backend_missing) {
        parallel_backend <- match.arg(parallel_backend)
    } else {
        parallel_backend <- "auto"
    }
    
    verbose <- .resolve_runtime_verbose(verbose)
    .with_log_session("createSCOPE", verbose = verbose, {
        .log_start("createSCOPE", verbose = verbose)
        
        detected <- NA_character_
        if (!is.null(data_dir) && dir.exists(data_dir)) {
            heuristic_detected <- if (
                dir.exists(file.path(data_dir, "filtered_feature_bc_matrix")) ||
                dir.exists(file.path(data_dir, "raw_feature_bc_matrix")) ||
                file.exists(file.path(data_dir, "spatial", "scalefactors_json.json"))
            ) {
                "visium"
            } else if (
                dir.exists(file.path(data_dir, "flatFiles")) ||
                dir.exists(file.path(data_dir, "DecodedFiles"))
            ) {
                "cosmx"
            } else if (
                file.exists(file.path(data_dir, "transcripts.parquet")) ||
                file.exists(file.path(data_dir, "transcripts.parquet.gz")) ||
                file.exists(file.path(data_dir, "cell_feature_matrix.h5")) ||
                dir.exists(file.path(data_dir, "analysis"))
            ) {
                "xenium"
            } else {
                NA_character_
            }
            detected <- if (!identical(prefer, "auto")) prefer else heuristic_detected
        }

        .log_key_values("createSCOPE", list(
            data_dir = if (!is.null(data_dir)) basename(path.expand(data_dir)) else NA_character_,
            prefer = prefer,
            detected = detected,
            sctransform = sctransform,
            grid_length = grid_length,
            seg_type = seg_type,
            num_workers = num_workers
        ))

        # Merge explicit parameters with dots
        explicit_args <- list(
            grid_length = grid_length,
            seg_type = seg_type,
            roi_path = roi_path,
            num_workers = num_workers,
            filter_genes = filter_genes,
            max_distance_to_nucleus = max_distance_to_nucleus,
            filter_molecules = filter_molecules,
            exclude_prefixes = exclude_prefixes,
            min_quality = min_quality,
            min_segmentation_points = min_segmentation_points,
            parallel_backend = parallel_backend
        )
        
        # Remove NULL values from explicit_args to let internal defaults take over
        explicit_args <- explicit_args[!sapply(explicit_args, is.null)]
        
        result <- do.call(.create_scope, c(list(
            data_dir = data_dir,
            prefer = prefer,
            verbose = verbose,
            sctransform = sctransform
        ), explicit_args, dots))

        .log_key_values("createSCOPE", list(
            grid_layers = if (!is.null(result@grid)) length(result@grid) else 0L,
            cell_layers = if (!is.null(result@cells)) length(result@cells) else 0L,
            coord_layers = if (!is.null(result@coord)) length(result@coord) else 0L
        ))

        .log_done("createSCOPE", verbose = verbose)
        result
    })
}
