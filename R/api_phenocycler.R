#' Create a marker-family policy for PhenoCycler bundles.
#' @description
#' Builds an explicit policy object that classifies PhenoCycler features into
#' punctate, diffuse-region-like, or boundary-associated families without
#' changing any legacy geneSCOPE defaults.
#' @param punctate_sources Character vector of source types treated as punctate.
#' @param diffuse_sources Character vector of source types treated as diffuse.
#' @param boundary_sources Character vector of source types treated as
#'   boundary-associated.
#' @param feature_family_overrides Optional named character vector mapping
#'   `feature_name -> family`.
#' @param source_family_overrides Optional named character vector mapping
#'   `source_type -> family`.
#' @param version Policy version string stored in manifests.
#' @return A named list with class `genescope_phenocycler_marker_policy`.
#' @export
makePhenoCyclerMarkerPolicy <- function(
    punctate_sources = c("punctate_event", "connected_component_event"),
    diffuse_sources = c("support_grid_point", "support_superpixel"),
    boundary_sources = c("support_contour_anchor"),
    feature_family_overrides = NULL,
    source_family_overrides = NULL,
    version = "phenocycler_marker_policy_v1"
) {
    policy <- list(
        punctate_sources = unique(as.character(punctate_sources)),
        diffuse_sources = unique(as.character(diffuse_sources)),
        boundary_sources = unique(as.character(boundary_sources)),
        feature_family_overrides = .pc_named_chr(feature_family_overrides),
        source_family_overrides = .pc_named_chr(source_family_overrides),
        defaults = c(
            punctate = "event",
            diffuse_region_like = "grid_field",
            boundary_associated = "weighted_anchor"
        ),
        version = as.character(version)[1]
    )
    class(policy) <- c("genescope_phenocycler_marker_policy", "list")
    policy
}

#' Create a diffuse-marker representation config for PhenoCycler bundles.
#' @description
#' Defines how diffuse and boundary-associated markers should be digitised into
#' grid fields, weighted anchors, and proxy support regions.
#' @param grid_lengths Numeric vector of grid lengths to materialise.
#' @param field_weight_col Column used as the weight when aggregating diffuse
#'   signal to field grids.
#' @param anchor_weight_col Column used as the weight when generating anchors.
#' @param area_col Column used to derive local scale and proxy regions.
#' @param region_mode Strategy used for support-region output.
#' @param include_boundary_in_dual Whether boundary-associated supports are
#'   included in the dual promotion layer.
#' @param min_weight Minimum retained weight.
#' @param version Config version string stored in manifests.
#' @return A named list with class `genescope_phenocycler_diffuse_config`.
#' @export
makePhenoCyclerDiffuseConfig <- function(
    grid_lengths = c(20, 40),
    field_weight_col = "support_score",
    anchor_weight_col = "support_score",
    area_col = "support_area_um2",
    region_mode = c("bbox_proxy", "none"),
    include_boundary_in_dual = TRUE,
    min_weight = 0,
    version = "phenocycler_diffuse_config_v1"
) {
    region_mode <- match.arg(region_mode)
    cfg <- list(
        grid_lengths = sort(unique(as.numeric(grid_lengths))),
        field_weight_col = as.character(field_weight_col)[1],
        anchor_weight_col = as.character(anchor_weight_col)[1],
        area_col = as.character(area_col)[1],
        region_mode = region_mode,
        include_boundary_in_dual = isTRUE(include_boundary_in_dual),
        min_weight = as.numeric(min_weight)[1],
        version = as.character(version)[1]
    )
    class(cfg) <- c("genescope_phenocycler_diffuse_config", "list")
    cfg
}

#' Create a geneSCOPE adapter config for PhenoCycler bundles.
#' @description
#' Defines how PhenoCycler representations are attached to a `scope_object` and
#' how, if explicitly requested, they are promoted into standard grid layers for
#' legacy geneSCOPE downstream functions.
#' @param grid_lengths Numeric vector of grid lengths to materialise.
#' @param promote Default promotion behaviour (`none`, `field`, `punctate`,
#'   `boundary`, or `dual`).
#' @param field_layer_prefix Prefix for promoted diffuse field layers.
#' @param punctate_layer_prefix Prefix for promoted punctate event layers.
#' @param boundary_layer_prefix Prefix for promoted boundary-support layers.
#' @param dual_layer_prefix Prefix for promoted dual layers.
#' @param namespace_dual_features Whether dual layers namespace feature names by
#'   family.
#' @param attach_density_tables Whether to attach wide field tables under
#'   `scope_obj@density`.
#' @param version Config version string stored in manifests.
#' @return A named list with class `genescope_phenocycler_adapter_config`.
#' @export
makePhenoCyclerAdapterConfig <- function(
    grid_lengths = c(20, 40),
    promote = c("none", "field", "punctate", "boundary", "dual"),
    field_layer_prefix = "pc_field",
    punctate_layer_prefix = "pc_punctate",
    boundary_layer_prefix = "pc_boundary",
    dual_layer_prefix = "pc_dual",
    namespace_dual_features = TRUE,
    attach_density_tables = TRUE,
    version = "phenocycler_adapter_config_v1"
) {
    promote <- match.arg(promote)
    cfg <- list(
        grid_lengths = sort(unique(as.numeric(grid_lengths))),
        promote = promote,
        field_layer_prefix = as.character(field_layer_prefix)[1],
        punctate_layer_prefix = as.character(punctate_layer_prefix)[1],
        boundary_layer_prefix = as.character(boundary_layer_prefix)[1],
        dual_layer_prefix = as.character(dual_layer_prefix)[1],
        namespace_dual_features = isTRUE(namespace_dual_features),
        attach_density_tables = isTRUE(attach_density_tables),
        version = as.character(version)[1]
    )
    class(cfg) <- c("genescope_phenocycler_adapter_config", "list")
    cfg
}

#' Ingest an organised PhenoCycler bundle.
#' @description
#' Reads the organised PhenoCycler bundle produced by the upstream QPTIFF
#' workflow, applies a marker-family policy, and materialises coordinate-first
#' plus field-first representations without touching any legacy geneSCOPE
#' ingest path.
#' @param bundle_dir Path to the organised PhenoCycler bundle.
#' @param marker_policy Marker-family policy object.
#' @param diffuse_config Diffuse representation config.
#' @param adapter_config geneSCOPE adapter config.
#' @param verbose Emit progress messages when TRUE.
#' @return A named list with class `genescope_phenocycler_bundle`.
#' @export
ingestPhenoCyclerBundle <- function(
    bundle_dir,
    marker_policy = makePhenoCyclerMarkerPolicy(),
    diffuse_config = makePhenoCyclerDiffuseConfig(),
    adapter_config = makePhenoCyclerAdapterConfig(),
    verbose = TRUE
) {
    bundle_dir <- normalizePath(bundle_dir, mustWork = TRUE)
    .pc_require_pkg("arrow", "ingestPhenoCyclerBundle()")
    .pc_require_pkg("jsonlite", "ingestPhenoCyclerBundle()")

    spec <- .pc_bundle_spec(bundle_dir)
    .pc_validate_bundle_spec(spec)

    events <- .pc_read_table(spec$protein_events_path)
    supports <- .pc_read_table(spec$signal_supports_path)
    cells <- .pc_read_optional_table(spec$cells_path, spec$cells_csv_path)
    cell_boundaries <- .pc_read_optional_table(spec$cell_boundaries_path)
    nucleus_boundaries <- .pc_read_optional_table(spec$nucleus_boundaries_path)
    data_semantics <- .pc_read_optional_json(spec$data_semantics_path)
    dataset_manifest <- .pc_read_optional_json(spec$dataset_manifest_path)

    events <- .pc_prepare_events(events, marker_policy)
    supports <- .pc_prepare_supports(supports, marker_policy, diffuse_config)
    feature_meta <- .pc_build_feature_metadata(events, supports, marker_policy)
    weighted_anchors <- .pc_build_weighted_anchors(supports)
    support_regions <- .pc_build_support_regions(supports, diffuse_config)
    bounds <- .pc_compute_bundle_bounds(events, supports, cells, cell_boundaries, nucleus_boundaries)
    grid_lengths <- sort(unique(c(diffuse_config$grid_lengths, adapter_config$grid_lengths)))
    grid_layers <- .pc_build_grid_layers(
        events = events,
        supports = supports,
        bounds = bounds,
        grid_lengths = grid_lengths,
        adapter_config = adapter_config,
        diffuse_config = diffuse_config
    )

    manifest <- list(
        bundle_dir = bundle_dir,
        protein_events_path = spec$protein_events_path,
        signal_supports_path = spec$signal_supports_path,
        cells_path = .pc_first_existing(c(spec$cells_path, spec$cells_csv_path)),
        cell_boundaries_path = .pc_first_existing(spec$cell_boundaries_path),
        nucleus_boundaries_path = .pc_first_existing(spec$nucleus_boundaries_path),
        data_semantics_path = .pc_first_existing(spec$data_semantics_path),
        dataset_manifest_path = .pc_first_existing(spec$dataset_manifest_path),
        marker_policy_version = marker_policy$version,
        diffuse_config_version = diffuse_config$version,
        adapter_config_version = adapter_config$version,
        event_rows = nrow(events),
        support_rows = nrow(supports),
        feature_count = nrow(feature_meta),
        grid_lengths = grid_lengths,
        bounds = bounds
    )

    if (verbose) {
        message(
            "[geneSCOPE::ingestPhenoCyclerBundle] bundle=",
            basename(bundle_dir),
            " events=", nrow(events),
            " supports=", nrow(supports),
            " features=", nrow(feature_meta),
            " grids=", paste(grid_lengths, collapse = ",")
        )
    }

    out <- list(
        manifest = manifest,
        marker_policy = marker_policy,
        diffuse_config = diffuse_config,
        adapter_config = adapter_config,
        inputs = list(
            protein_events = events,
            signal_supports = supports,
            cells = cells,
            cell_boundaries = cell_boundaries,
            nucleus_boundaries = nucleus_boundaries
        ),
        metadata = list(
            data_semantics = data_semantics,
            dataset_manifest = dataset_manifest,
            feature_metadata = feature_meta
        ),
        representations = list(
            weighted_anchors = weighted_anchors,
            support_regions = support_regions,
            grid_layers = grid_layers
        )
    )
    class(out) <- c("genescope_phenocycler_bundle", "list")
    out
}

#' Create a scope object from an organised PhenoCycler bundle.
#' @description
#' Builds a new `scope_object` using only the additive PhenoCycler module. This
#' does not modify any legacy geneSCOPE constructor or dispatch behaviour.
#' @param bundle_dir Path to the organised PhenoCycler bundle.
#' @param marker_policy Marker-family policy object.
#' @param diffuse_config Diffuse representation config.
#' @param adapter_config geneSCOPE adapter config.
#' @param promote Which PhenoCycler representation, if any, should be promoted
#'   into standard `@grid` layers immediately.
#' @param verbose Emit progress messages when TRUE.
#' @return A `scope_object`.
#' @export
createSCOPE_phenocycler <- function(
    bundle_dir,
    marker_policy = makePhenoCyclerMarkerPolicy(),
    diffuse_config = makePhenoCyclerDiffuseConfig(),
    adapter_config = makePhenoCyclerAdapterConfig(),
    promote = adapter_config$promote,
    verbose = TRUE
) {
    scope_obj <- methods::new(
        "scope_object",
        coord = list(),
        grid = list(),
        meta.data = data.frame(row.names = character()),
        cells = list(),
        stats = list(platform = "phenocycler", modalities = list()),
        density = list()
    )
    addPhenoCycler(
        scope_obj = scope_obj,
        bundle_dir = bundle_dir,
        marker_policy = marker_policy,
        diffuse_config = diffuse_config,
        adapter_config = adapter_config,
        promote = promote,
        verbose = verbose
    )
}

#' Attach PhenoCycler data to an existing scope object.
#' @description
#' Adds a PhenoCycler modality sidecar to an existing `scope_object`. New layers
#' are opt-in and default to staying under modality-specific namespaces.
#' @param scope_obj A `scope_object`.
#' @param bundle_dir Path to the organised PhenoCycler bundle.
#' @param bundle Optional already-ingested bundle object from
#'   `ingestPhenoCyclerBundle()`.
#' @param marker_policy Marker-family policy object.
#' @param diffuse_config Diffuse representation config.
#' @param adapter_config geneSCOPE adapter config.
#' @param promote Which representation, if any, to promote into standard
#'   `@grid` layers.
#' @param overwrite_existing_modality Whether to overwrite an existing
#'   PhenoCycler sidecar.
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object`.
#' @export
addPhenoCycler <- function(
    scope_obj,
    bundle_dir = NULL,
    bundle = NULL,
    marker_policy = makePhenoCyclerMarkerPolicy(),
    diffuse_config = makePhenoCyclerDiffuseConfig(),
    adapter_config = makePhenoCyclerAdapterConfig(),
    promote = adapter_config$promote,
    overwrite_existing_modality = FALSE,
    verbose = TRUE
) {
    if (!inherits(scope_obj, "scope_object")) {
        stop("`scope_obj` must be a scope_object.")
    }
    if (is.null(bundle)) {
        if (is.null(bundle_dir)) stop("Provide either `bundle_dir` or `bundle`.")
        bundle <- ingestPhenoCyclerBundle(
            bundle_dir = bundle_dir,
            marker_policy = marker_policy,
            diffuse_config = diffuse_config,
            adapter_config = adapter_config,
            verbose = verbose
        )
    }
    if (!is.list(scope_obj@stats$modalities)) {
        scope_obj@stats$modalities <- list()
    }
    if (!overwrite_existing_modality && !is.null(scope_obj@stats$modalities$phenocycler)) {
        stop("scope_obj already has a `phenocycler` modality sidecar. Set `overwrite_existing_modality=TRUE` to replace it.")
    }

    scope_obj@stats$modalities$phenocycler <- list(
        manifest = bundle$manifest,
        marker_policy = bundle$marker_policy,
        diffuse_config = bundle$diffuse_config,
        adapter_config = bundle$adapter_config,
        feature_metadata = bundle$metadata$feature_metadata,
        data_semantics = bundle$metadata$data_semantics,
        dataset_manifest = bundle$metadata$dataset_manifest,
        grid_layers = bundle$representations$grid_layers,
        support_regions = bundle$representations$support_regions
    )

    scope_obj@coord$phenocycler_events <- bundle$inputs$protein_events
    scope_obj@coord$phenocycler_supports <- bundle$inputs$signal_supports
    scope_obj@coord$phenocycler_weighted_anchors <- bundle$representations$weighted_anchors
    scope_obj@coord$phenocycler_cell_centroids <- bundle$inputs$cells
    scope_obj@coord$phenocycler_segmentation_cell <- bundle$inputs$cell_boundaries
    scope_obj@coord$phenocycler_segmentation_nucleus <- bundle$inputs$nucleus_boundaries

    if (isTRUE(bundle$adapter_config$attach_density_tables)) {
        density_tables <- .pc_density_tables_from_grid_layers(bundle$representations$grid_layers)
        for (nm in names(density_tables)) {
            scope_obj@density[[nm]] <- density_tables[[nm]]
        }
    }

    promote <- match.arg(as.character(promote)[1], c("none", "field", "punctate", "boundary", "dual"))
    if (!identical(promote, "none")) {
        scope_obj <- promotePhenoCyclerToGridLayer(
            scope_obj = scope_obj,
            source = promote,
            grid_lengths = bundle$adapter_config$grid_lengths,
            overwrite = FALSE,
            verbose = verbose
        )
    }

    if (verbose) {
        message(
            "[geneSCOPE::addPhenoCycler] attached modality with ",
            nrow(bundle$metadata$feature_metadata),
            " features; promote=", promote
        )
    }
    scope_obj
}

#' Promote a PhenoCycler representation into standard grid layers.
#' @description
#' Copies a selected PhenoCycler representation into the standard `@grid`
#' namespace so existing geneSCOPE grid-based functions can be used explicitly.
#' This bridge is opt-in and never changes legacy defaults.
#' @param scope_obj A `scope_object` with an attached PhenoCycler modality.
#' @param source Which representation to promote.
#' @param grid_lengths Optional subset of grid lengths to promote.
#' @param overwrite Whether to overwrite an existing grid layer of the same name.
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object`.
#' @export
promotePhenoCyclerToGridLayer <- function(
    scope_obj,
    source = c("field", "punctate", "boundary", "dual"),
    grid_lengths = NULL,
    overwrite = FALSE,
    verbose = TRUE
) {
    if (!inherits(scope_obj, "scope_object")) {
        stop("`scope_obj` must be a scope_object.")
    }
    source <- match.arg(source)
    modality <- scope_obj@stats$modalities$phenocycler
    if (is.null(modality)) {
        stop("No `phenocycler` modality is attached to this scope_object.")
    }

    grid_layers <- modality$grid_layers
    lengths_available <- as.numeric(sub("^grid", "", names(grid_layers)))
    if (is.null(grid_lengths)) {
        grid_lengths <- lengths_available
    }
    grid_lengths <- sort(unique(as.numeric(grid_lengths)))
    missing_lengths <- setdiff(grid_lengths, lengths_available)
    if (length(missing_lengths)) {
        stop("Requested grid lengths not available in phenocycler modality: ", paste(missing_lengths, collapse = ", "))
    }

    for (lg in grid_lengths) {
        key <- paste0("grid", .pc_fmt_number(lg))
        layer_bundle <- grid_layers[[key]]
        if (is.null(layer_bundle[[source]])) {
            stop("No `", source, "` layer available for ", key, ".")
        }
        target_name <- .pc_promoted_layer_name(modality$adapter_config, source, lg)
        if (!overwrite && !is.null(scope_obj@grid[[target_name]])) {
            stop("Grid layer `", target_name, "` already exists. Set `overwrite=TRUE` to replace it.")
        }
        feature_meta <- .pc_promoted_feature_metadata(modality$feature_metadata, source, modality$adapter_config)
        promoted_layer <- layer_bundle[[source]]
        density_table <- NULL
        if (isTRUE(modality$adapter_config$attach_density_tables)) {
            density_table <- layer_bundle$density[[source]]
        }
        prepared <- .pc_prepare_promoted_feature_namespace(
            scope_obj = scope_obj,
            promoted_layer = promoted_layer,
            density_table = density_table,
            feature_metadata = feature_meta,
            namespace_prefix = "pc::"
        )
        promoted_layer <- prepared$promoted_layer
        density_table <- prepared$density_table
        feature_meta <- prepared$feature_metadata
        if (!is.null(density_table)) {
            promoted_layer$densityDF <- density_table
            scope_obj@density[[target_name]] <- density_table
            scope_obj@density[[paste0(target_name, "_density")]] <- density_table
        }
        scope_obj@grid[[target_name]] <- promoted_layer
        scope_obj <- .pc_bind_feature_metadata(scope_obj, feature_meta)
        if (verbose) {
            message(
                "[geneSCOPE::promotePhenoCyclerToGridLayer] ",
                source, " -> @grid$", target_name,
                " (", nrow(layer_bundle[[source]]$grid_info), " tiles",
                if (length(prepared$collisions)) paste0("; namespaced ", length(prepared$collisions), " colliding feature(s)") else "",
                ")"
            )
        }
    }
    scope_obj
}

.pc_bundle_spec <- function(bundle_dir) {
    list(
        bundle_dir = bundle_dir,
        protein_events_path = file.path(bundle_dir, "molecules", "protein_events.parquet"),
        signal_supports_path = file.path(bundle_dir, "molecules", "signal_supports.parquet"),
        cells_path = file.path(bundle_dir, "cells", "cells.parquet"),
        cells_csv_path = file.path(bundle_dir, "cells", "cells.csv.gz"),
        cell_boundaries_path = file.path(bundle_dir, "segmentation", "cell_boundaries.parquet"),
        nucleus_boundaries_path = file.path(bundle_dir, "segmentation", "nucleus_boundaries.parquet"),
        data_semantics_path = file.path(bundle_dir, "metadata", "data_semantics.json"),
        dataset_manifest_path = file.path(bundle_dir, "metadata", "dataset_manifest.json")
    )
}

.pc_validate_bundle_spec <- function(spec) {
    required <- c("protein_events_path", "signal_supports_path")
    missing <- required[!vapply(spec[required], file.exists, logical(1))]
    if (length(missing)) {
        stop(
            "PhenoCycler bundle is missing required files: ",
            paste(vapply(missing, function(k) basename(spec[[k]]), character(1)), collapse = ", ")
        )
    }
    invisible(spec)
}

.pc_prepare_events <- function(events, marker_policy) {
    if (!nrow(events)) return(events)
    dt <- as.data.table(events)
    dt[, marker_family := .pc_classify_family(feature_name, source_type, marker_policy)]
    dt[, representation_default := marker_policy$defaults[marker_family]]
    dt
}

.pc_prepare_supports <- function(supports, marker_policy, diffuse_config) {
    if (!nrow(supports)) return(supports)
    dt <- as.data.table(supports)
    dt[, marker_family := .pc_classify_family(feature_name, source_type, marker_policy)]
    dt[, representation_default := marker_policy$defaults[marker_family]]
    dt[, support_weight := .pc_extract_numeric(dt, diffuse_config$field_weight_col, default = 1)]
    dt[!is.finite(support_weight), support_weight := 0]
    dt[support_weight < diffuse_config$min_weight, support_weight := 0]
    dt[, local_scale := .pc_support_scale(dt, diffuse_config$area_col)]
    dt
}

.pc_build_feature_metadata <- function(events, supports, marker_policy) {
    events_dt <- if (!is.null(events) && nrow(events)) {
        unique(as.data.table(events)[, .(
            feature_name,
            marker_family,
            representation_default,
            source_observed = as.character(source_type)
        )])
    } else {
        data.table(feature_name = character(), marker_family = character(), representation_default = character(), source_observed = character())
    }
    supports_dt <- if (!is.null(supports) && nrow(supports)) {
        unique(as.data.table(supports)[, .(
            feature_name,
            marker_family,
            representation_default,
            source_observed = as.character(source_type)
        )])
    } else {
        data.table(feature_name = character(), marker_family = character(), representation_default = character(), source_observed = character())
    }
    feature_dt <- rbindlist(list(events_dt, supports_dt), use.names = TRUE, fill = TRUE)
    if (!nrow(feature_dt)) {
        return(data.frame(
            feature_name = character(),
            marker_family = character(),
            representation_default = character(),
            source_observed = character(),
            row.names = character()
        ))
    }
    out <- feature_dt[, .(
        marker_family = marker_family[[1]],
        representation_default = representation_default[[1]],
        source_observed = paste(sort(unique(source_observed)), collapse = "|")
    ), by = feature_name]
    out <- out[order(feature_name)]
    out_df <- as.data.frame(out)
    rownames(out_df) <- out_df$feature_name
    out_df
}

.pc_build_weighted_anchors <- function(supports) {
    if (is.null(supports) || !nrow(supports)) {
        return(data.frame(
            anchor_id = character(),
            feature_name = character(),
            x_location = numeric(),
            y_location = numeric(),
            weight = numeric(),
            local_scale = numeric(),
            marker_family = character(),
            representation_type = character(),
            source_type = character(),
            stringsAsFactors = FALSE
        ))
    }
    dt <- as.data.table(supports)
    out <- dt[, .(
        anchor_id = as.character(.pc_val_or_seq(support_id, "anchor")),
        feature_name = as.character(feature_name),
        x_location = as.numeric(x_location),
        y_location = as.numeric(y_location),
        weight = as.numeric(support_weight),
        local_scale = as.numeric(local_scale),
        marker_family = as.character(marker_family),
        representation_type = fifelse(
            marker_family == "boundary_associated",
            "boundary_weighted_anchor",
            "diffuse_weighted_anchor"
        ),
        source_type = as.character(source_type)
    )]
    as.data.frame(out)
}

.pc_build_support_regions <- function(supports, diffuse_config) {
    if (is.null(supports) || !nrow(supports) || identical(diffuse_config$region_mode, "none")) {
        return(data.frame(
            region_id = character(),
            feature_name = character(),
            x_center = numeric(),
            y_center = numeric(),
            support_area = numeric(),
            support_bbox_xmin = numeric(),
            support_bbox_xmax = numeric(),
            support_bbox_ymin = numeric(),
            support_bbox_ymax = numeric(),
            marker_family = character(),
            representation_type = character(),
            stringsAsFactors = FALSE
        ))
    }
    dt <- as.data.table(supports)
    scale_use <- dt$local_scale
    scale_use[!is.finite(scale_use)] <- 0
    out <- dt[, .(
        region_id = as.character(.pc_val_or_seq(support_id, "region")),
        feature_name = as.character(feature_name),
        x_center = as.numeric(x_location),
        y_center = as.numeric(y_location),
        support_area = .pc_extract_numeric(.SD, diffuse_config$area_col, default = NA_real_),
        support_bbox_xmin = as.numeric(x_location) - scale_use,
        support_bbox_xmax = as.numeric(x_location) + scale_use,
        support_bbox_ymin = as.numeric(y_location) - scale_use,
        support_bbox_ymax = as.numeric(y_location) + scale_use,
        marker_family = as.character(marker_family),
        representation_type = "support_region_bbox_proxy"
    )]
    as.data.frame(out)
}

.pc_compute_bundle_bounds <- function(events, supports, cells, cell_boundaries, nucleus_boundaries) {
    xs <- c(
        .pc_pull_num(events, "x_location"),
        .pc_pull_num(supports, "x_location"),
        .pc_pull_num(cells, "x_centroid"),
        .pc_pull_num(cell_boundaries, "vertex_x"),
        .pc_pull_num(nucleus_boundaries, "vertex_x")
    )
    ys <- c(
        .pc_pull_num(events, "y_location"),
        .pc_pull_num(supports, "y_location"),
        .pc_pull_num(cells, "y_centroid"),
        .pc_pull_num(cell_boundaries, "vertex_y"),
        .pc_pull_num(nucleus_boundaries, "vertex_y")
    )
    xs <- xs[is.finite(xs)]
    ys <- ys[is.finite(ys)]
    if (!length(xs) || !length(ys)) {
        stop("Could not determine bundle bounds from PhenoCycler inputs.")
    }
    list(
        xmin = min(xs),
        xmax = max(xs),
        ymin = min(ys),
        ymax = max(ys)
    )
}

.pc_build_grid_layers <- function(events,
                                  supports,
                                  bounds,
                                  grid_lengths,
                                  adapter_config,
                                  diffuse_config) {
    out <- vector("list", length(grid_lengths))
    names(out) <- paste0("grid", vapply(grid_lengths, .pc_fmt_number, character(1)))
    punctate_dt <- as.data.table(events)
    support_dt <- as.data.table(supports)

    for (lg in grid_lengths) {
        grid_info <- .pc_build_grid_info(bounds, lg)
        punctate_counts <- .pc_aggregate_to_grid(
            dt = punctate_dt,
            grid_info = grid_info,
            grid_length = lg,
            x_col = "x_location",
            y_col = "y_location",
            feature_col = "feature_name",
            value_col = NULL
        )
        field_supports <- support_dt[marker_family == "diffuse_region_like"]
        boundary_supports <- support_dt[marker_family == "boundary_associated"]
        field_counts <- .pc_aggregate_to_grid(
            dt = field_supports,
            grid_info = grid_info,
            grid_length = lg,
            x_col = "x_location",
            y_col = "y_location",
            feature_col = "feature_name",
            value_col = "support_weight"
        )
        boundary_counts <- .pc_aggregate_to_grid(
            dt = boundary_supports,
            grid_info = grid_info,
            grid_length = lg,
            x_col = "x_location",
            y_col = "y_location",
            feature_col = "feature_name",
            value_col = "support_weight"
        )
        dual_counts <- .pc_build_dual_counts(
            punctate_counts = punctate_counts,
            field_counts = field_counts,
            boundary_counts = boundary_counts,
            namespace_dual_features = isTRUE(adapter_config$namespace_dual_features),
            include_boundary = isTRUE(diffuse_config$include_boundary_in_dual)
        )
        out[[paste0("grid", .pc_fmt_number(lg))]] <- list(
            punctate = .pc_make_grid_layer(grid_info, punctate_counts, lg, "punctate_event_layer"),
            field = .pc_make_grid_layer(grid_info, field_counts, lg, "diffuse_grid_field"),
            boundary = .pc_make_grid_layer(grid_info, boundary_counts, lg, "boundary_support_field"),
            dual = .pc_make_grid_layer(grid_info, dual_counts, lg, "phenocycler_dual_layer"),
            density = list(
                punctate = .pc_counts_to_density_table(grid_info, punctate_counts),
                field = .pc_counts_to_density_table(grid_info, field_counts),
                boundary = .pc_counts_to_density_table(grid_info, boundary_counts),
                dual = .pc_counts_to_density_table(grid_info, dual_counts)
            )
        )
    }
    out
}

.pc_make_grid_layer <- function(grid_info, counts_dt, grid_length, representation_type) {
    counts_dt <- as.data.table(counts_dt)
    counts_dt <- counts_dt[is.finite(count) & count != 0]
    list(
        grid_info = as.data.frame(grid_info),
        counts = as.data.frame(counts_dt),
        grid_length = as.numeric(grid_length),
        xbins_eff = max(grid_info$gx),
        ybins_eff = max(grid_info$gy),
        representation_type = representation_type,
        source_platform = "phenocycler"
    )
}

.pc_counts_to_density_table <- function(grid_info, counts_dt) {
    grid_key <- as.data.table(grid_info)[, .(grid_id)]
    counts_dt <- as.data.table(counts_dt)
    if (!nrow(counts_dt)) {
        out <- as.data.frame(grid_key)
        rownames(out) <- out$grid_id
        return(out)
    }
    wide <- data.table::dcast(
        counts_dt,
        grid_id ~ gene,
        value.var = "count",
        fill = 0
    )
    out <- merge(grid_key, wide, by = "grid_id", all.x = TRUE, sort = FALSE)
    out[is.na(out)] <- 0
    out <- as.data.frame(out)
    rownames(out) <- out$grid_id
    out
}

.pc_density_tables_from_grid_layers <- function(grid_layers) {
    out <- list()
    for (key in names(grid_layers)) {
        layer_bundle <- grid_layers[[key]]
        for (src in names(layer_bundle$density)) {
            out[[paste0("phenocycler_", src, "_", key)]] <- layer_bundle$density[[src]]
        }
    }
    out
}

.pc_prepare_promoted_feature_namespace <- function(scope_obj,
                                                   promoted_layer,
                                                   density_table = NULL,
                                                   feature_metadata = NULL,
                                                   namespace_prefix = "pc::") {
    if (is.null(feature_metadata) || !nrow(feature_metadata)) {
        return(list(
            promoted_layer = promoted_layer,
            density_table = density_table,
            feature_metadata = feature_metadata,
            collisions = character(),
            rename_map = character()
        ))
    }

    feature_metadata <- as.data.frame(feature_metadata)
    if (is.null(rownames(feature_metadata))) {
        rownames(feature_metadata) <- .pc_coalesce(feature_metadata$feature_name, seq_len(nrow(feature_metadata)))
    }

    existing_rows <- if (nrow(scope_obj@meta.data)) rownames(scope_obj@meta.data) else character()
    collisions <- intersect(rownames(feature_metadata), existing_rows)
    if (!length(collisions)) {
        return(list(
            promoted_layer = promoted_layer,
            density_table = density_table,
            feature_metadata = feature_metadata,
            collisions = character(),
            rename_map = character()
        ))
    }

    rename_map <- stats::setNames(paste0(namespace_prefix, collisions), collisions)

    if (!is.null(promoted_layer$counts) && nrow(promoted_layer$counts) && "gene" %in% names(promoted_layer$counts)) {
        idx <- promoted_layer$counts$gene %in% names(rename_map)
        if (any(idx)) {
            promoted_layer$counts$gene[idx] <- unname(rename_map[promoted_layer$counts$gene[idx]])
        }
    }

    if (!is.null(density_table)) {
        cols <- intersect(names(rename_map), colnames(density_table))
        if (length(cols)) {
            colnames(density_table)[match(cols, colnames(density_table))] <- unname(rename_map[cols])
        }
    }

    old_rows <- rownames(feature_metadata)
    row_idx <- old_rows %in% names(rename_map)
    if (any(row_idx)) {
        rownames(feature_metadata)[row_idx] <- unname(rename_map[old_rows[row_idx]])
    }
    if ("feature_name" %in% names(feature_metadata)) {
        feat_idx <- feature_metadata$feature_name %in% names(rename_map)
        if (any(feat_idx)) {
            feature_metadata$feature_name[feat_idx] <- unname(rename_map[feature_metadata$feature_name[feat_idx]])
        }
    }

    list(
        promoted_layer = promoted_layer,
        density_table = density_table,
        feature_metadata = feature_metadata,
        collisions = collisions,
        rename_map = rename_map
    )
}

.pc_build_dual_counts <- function(punctate_counts,
                                  field_counts,
                                  boundary_counts,
                                  namespace_dual_features = TRUE,
                                  include_boundary = TRUE) {
    parts <- list()
    if (!is.null(punctate_counts) && nrow(punctate_counts)) {
        dt <- copy(as.data.table(punctate_counts))
        if (namespace_dual_features) dt[, gene := paste0("punctate::", gene)]
        parts[[length(parts) + 1L]] <- dt
    }
    if (!is.null(field_counts) && nrow(field_counts)) {
        dt <- copy(as.data.table(field_counts))
        if (namespace_dual_features) dt[, gene := paste0("diffuse::", gene)]
        parts[[length(parts) + 1L]] <- dt
    }
    if (isTRUE(include_boundary) && !is.null(boundary_counts) && nrow(boundary_counts)) {
        dt <- copy(as.data.table(boundary_counts))
        if (namespace_dual_features) dt[, gene := paste0("boundary::", gene)]
        parts[[length(parts) + 1L]] <- dt
    }
    if (!length(parts)) {
        return(data.table(grid_id = character(), gene = character(), count = numeric()))
    }
    merged <- rbindlist(parts, use.names = TRUE, fill = TRUE)
    merged[, .(count = sum(as.numeric(count), na.rm = TRUE)), by = .(grid_id, gene)]
}

.pc_build_grid_info <- function(bounds, grid_length) {
    grid_length <- as.numeric(grid_length)[1]
    x0 <- floor(bounds$xmin / grid_length) * grid_length
    x1 <- ceiling(bounds$xmax / grid_length) * grid_length
    y0 <- floor(bounds$ymin / grid_length) * grid_length
    y1 <- ceiling(bounds$ymax / grid_length) * grid_length
    x_starts <- seq(x0, x1 - grid_length, by = grid_length)
    y_starts <- seq(y0, y1 - grid_length, by = grid_length)
    if (!length(x_starts) || !length(y_starts)) {
        stop("Unable to construct grid info for grid_length=", grid_length)
    }
    grid_dt <- CJ(gx = seq_along(x_starts), gy = seq_along(y_starts), sorted = TRUE)
    grid_dt[, xmin := x_starts[gx]]
    grid_dt[, xmax := xmin + grid_length]
    grid_dt[, ymin := y_starts[gy]]
    grid_dt[, ymax := ymin + grid_length]
    grid_dt[, center_x := xmin + grid_length / 2]
    grid_dt[, center_y := ymin + grid_length / 2]
    grid_dt[, width := grid_length]
    grid_dt[, height := grid_length]
    grid_dt[, idx := .I]
    grid_dt[, grid_id := paste0("g", gx, "_", gy)]
    data.table::setcolorder(
        grid_dt,
        c("grid_id", "gx", "gy", "xmin", "xmax", "ymin", "ymax", "center_x", "center_y", "width", "height", "idx")
    )
    grid_dt[]
}

.pc_aggregate_to_grid <- function(dt,
                                  grid_info,
                                  grid_length,
                                  x_col,
                                  y_col,
                                  feature_col,
                                  value_col = NULL) {
    if (is.null(dt) || !nrow(dt)) {
        return(data.table(grid_id = character(), gene = character(), count = numeric()))
    }
    dt <- as.data.table(dt)
    if (!all(c(x_col, y_col, feature_col) %in% names(dt))) {
        return(data.table(grid_id = character(), gene = character(), count = numeric()))
    }
    grid_info <- as.data.table(grid_info)
    x0 <- min(grid_info$xmin)
    y0 <- min(grid_info$ymin)
    xbins_eff <- max(grid_info$gx)
    ybins_eff <- max(grid_info$gy)
    work <- copy(dt)
    work[, x__ := as.numeric(get(x_col))]
    work[, y__ := as.numeric(get(y_col))]
    work[, feature__ := as.character(get(feature_col))]
    if (is.null(value_col)) {
        work[, value__ := 1]
    } else {
        work[, value__ := .pc_extract_numeric(.SD, value_col, default = 0)]
    }
    work <- work[is.finite(x__) & is.finite(y__) & is.finite(value__) & nzchar(feature__)]
    if (!nrow(work)) {
        return(data.table(grid_id = character(), gene = character(), count = numeric()))
    }
    work[, gx := floor((x__ - x0) / grid_length) + 1L]
    work[, gy := floor((y__ - y0) / grid_length) + 1L]
    work <- work[gx >= 1L & gx <= xbins_eff & gy >= 1L & gy <= ybins_eff]
    if (!nrow(work)) {
        return(data.table(grid_id = character(), gene = character(), count = numeric()))
    }
    work[, grid_id := paste0("g", gx, "_", gy)]
    work[, .(count = sum(value__, na.rm = TRUE)), by = .(grid_id, gene = feature__)]
}

.pc_bind_feature_metadata <- function(scope_obj, feature_metadata) {
    if (is.null(feature_metadata) || !nrow(feature_metadata)) return(scope_obj)
    feature_metadata <- as.data.frame(feature_metadata)
    if (is.null(rownames(feature_metadata))) {
        rownames(feature_metadata) <- .pc_coalesce(feature_metadata$feature_name, seq_len(nrow(feature_metadata)))
    }
    if (!nrow(scope_obj@meta.data)) {
        scope_obj@meta.data <- feature_metadata
        return(scope_obj)
    }
    all_rows <- union(rownames(scope_obj@meta.data), rownames(feature_metadata))
    existing <- scope_obj@meta.data
    missing_cols <- setdiff(names(feature_metadata), names(existing))
    if (length(missing_cols)) {
        for (nm in missing_cols) existing[[nm]] <- NA
    }
    missing_cols_new <- setdiff(names(existing), names(feature_metadata))
    if (length(missing_cols_new)) {
        for (nm in missing_cols_new) feature_metadata[[nm]] <- NA
    }
    existing <- existing[, names(feature_metadata), drop = FALSE]
    feature_metadata <- feature_metadata[, names(existing), drop = FALSE]
    merged <- rbind(existing[setdiff(rownames(existing), rownames(feature_metadata)), , drop = FALSE], feature_metadata)
    merged <- merged[match(all_rows, rownames(merged)), , drop = FALSE]
    scope_obj@meta.data <- merged
    scope_obj
}

.pc_promoted_feature_metadata <- function(feature_metadata, source, adapter_config) {
    if (is.null(feature_metadata) || !nrow(feature_metadata)) return(feature_metadata)
    df <- as.data.frame(feature_metadata)
    if (identical(source, "dual")) {
        parts <- list()
        punctate <- df[df$marker_family == "punctate", , drop = FALSE]
        if (nrow(punctate)) {
            rownames(punctate) <- paste0(if (adapter_config$namespace_dual_features) "punctate::" else "", punctate$feature_name)
            punctate$feature_name <- rownames(punctate)
            punctate$representation_default <- "dual_punctate"
            parts[[length(parts) + 1L]] <- punctate
        }
        diffuse <- df[df$marker_family == "diffuse_region_like", , drop = FALSE]
        if (nrow(diffuse)) {
            rownames(diffuse) <- paste0(if (adapter_config$namespace_dual_features) "diffuse::" else "", diffuse$feature_name)
            diffuse$feature_name <- rownames(diffuse)
            diffuse$representation_default <- "dual_diffuse"
            parts[[length(parts) + 1L]] <- diffuse
        }
        boundary <- df[df$marker_family == "boundary_associated", , drop = FALSE]
        if (nrow(boundary)) {
            rownames(boundary) <- paste0(if (adapter_config$namespace_dual_features) "boundary::" else "", boundary$feature_name)
            boundary$feature_name <- rownames(boundary)
            boundary$representation_default <- "dual_boundary"
            parts[[length(parts) + 1L]] <- boundary
        }
        if (!length(parts)) return(df[0, , drop = FALSE])
        out <- do.call(rbind, parts)
        return(out)
    }
    family_keep <- switch(
        source,
        punctate = "punctate",
        field = "diffuse_region_like",
        boundary = "boundary_associated",
        stop("Unknown promoted source: ", source)
    )
    df[df$marker_family == family_keep, , drop = FALSE]
}

.pc_promoted_layer_name <- function(adapter_config, source, grid_length) {
    suffix <- .pc_fmt_number(grid_length)
    prefix <- switch(
        source,
        field = adapter_config$field_layer_prefix,
        punctate = adapter_config$punctate_layer_prefix,
        boundary = adapter_config$boundary_layer_prefix,
        dual = adapter_config$dual_layer_prefix,
        stop("Unknown source: ", source)
    )
    paste0(prefix, suffix)
}

.pc_classify_family <- function(feature_name, source_type, marker_policy) {
    feature_name <- as.character(feature_name)
    source_type <- as.character(source_type)
    out <- rep("diffuse_region_like", length(source_type))
    src_override <- marker_policy$source_family_overrides
    feat_override <- marker_policy$feature_family_overrides
    if (length(source_type)) {
        out[source_type %in% marker_policy$punctate_sources] <- "punctate"
        out[source_type %in% marker_policy$diffuse_sources] <- "diffuse_region_like"
        out[source_type %in% marker_policy$boundary_sources] <- "boundary_associated"
    }
    if (length(src_override)) {
        idx <- match(source_type, names(src_override))
        keep <- !is.na(idx)
        out[keep] <- unname(src_override[idx[keep]])
    }
    if (length(feat_override)) {
        idx <- match(feature_name, names(feat_override))
        keep <- !is.na(idx)
        out[keep] <- unname(feat_override[idx[keep]])
    }
    out
}

.pc_extract_numeric <- function(df, column, default = NA_real_) {
    if (!nzchar(column) || !column %in% names(df)) {
        return(rep(default, nrow(df)))
    }
    suppressWarnings(as.numeric(df[[column]]))
}

.pc_support_scale <- function(df, area_col) {
    area <- .pc_extract_numeric(df, area_col, default = NA_real_)
    out <- rep(NA_real_, length(area))
    keep <- is.finite(area) & area > 0
    out[keep] <- sqrt(area[keep] / pi)
    out
}

.pc_pull_num <- function(df, col) {
    if (is.null(df) || !is.data.frame(df) || !col %in% names(df)) return(numeric())
    suppressWarnings(as.numeric(df[[col]]))
}

.pc_val_or_seq <- function(x, prefix) {
    if (is.null(x)) return(paste0(prefix, "_", seq_len(0)))
    x_chr <- as.character(x)
    miss <- is.na(x_chr) | !nzchar(x_chr)
    if (any(miss)) x_chr[miss] <- paste0(prefix, "_", which(miss))
    x_chr
}

.pc_first_existing <- function(paths) {
    paths <- paths[file.exists(paths)]
    if (!length(paths)) return(NA_character_)
    normalizePath(paths[[1]], mustWork = TRUE)
}

.pc_read_table <- function(path) {
    as.data.frame(arrow::read_parquet(path))
}

.pc_read_optional_table <- function(path, fallback_csv = NULL) {
    if (!is.null(path) && file.exists(path)) {
        return(.pc_read_table(path))
    }
    if (!is.null(fallback_csv) && file.exists(fallback_csv)) {
        return(as.data.frame(data.table::fread(fallback_csv)))
    }
    data.frame()
}

.pc_read_optional_json <- function(path) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    jsonlite::read_json(path, simplifyVector = TRUE)
}

.pc_require_pkg <- function(pkg, fn) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package `", pkg, "` is required for ", fn, call. = FALSE)
    }
    invisible(TRUE)
}

.pc_named_chr <- function(x) {
    if (is.null(x)) return(character())
    x <- as.character(x)
    if (is.null(names(x))) {
        stop("Override vectors must be named.", call. = FALSE)
    }
    x
}

.pc_fmt_number <- function(x) {
    out <- format(as.numeric(x)[1], scientific = FALSE, trim = TRUE)
    if (grepl("\\.", out, fixed = FALSE)) {
        out <- sub("0+$", "", out)
        out <- sub("\\.$", "", out)
    }
    out
}

.pc_write_example_bundle <- function(out_dir) {
    .pc_require_pkg("arrow", ".pc_write_example_bundle()")
    .pc_require_pkg("jsonlite", ".pc_write_example_bundle()")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir, "molecules"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir, "cells"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir, "segmentation"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir, "metadata"), recursive = TRUE, showWarnings = FALSE)

    protein_events <- data.table(
        event_id = paste0("e", 1:8),
        feature_name = c("PDGFRA", "PDGFRA", "VIM", "VIM", "SLC1A3", "SLC1A3", "CLDN5", "CLDN5"),
        x_location = c(12, 18, 52, 58, 95, 102, 140, 148),
        y_location = c(12, 20, 18, 24, 20, 30, 18, 26),
        z_location = 0,
        source_type = c("punctate_event", "connected_component_event", "punctate_event", "connected_component_event", "punctate_event", "punctate_event", "connected_component_event", "connected_component_event"),
        intensity_raw = c(120, 110, 90, 92, 60, 65, 130, 128),
        intensity_norm = c(1.2, 1.1, 0.9, 0.92, 0.6, 0.65, 1.3, 1.28),
        event_area_um2 = c(1, 2, 1, 2, 1, 1, 2, 2),
        qv_or_score = 30,
        cell_id = c("c1", "c1", "c2", "c2", "c3", "c3", "c4", "c4"),
        assignment_confidence = "approximate",
        assignment_status = "interior",
        overlaps_nucleus = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
        tile_id = "tile_001",
        fov_name = "fov_1",
        segmentation_version = "route_b_example",
        event_call_method = "synthetic",
        resolution_level_used = 1,
        valid_pixel_fraction_local = 1,
        is_recovered_region = FALSE
    )
    signal_supports <- data.table(
        support_id = paste0("s", 1:8),
        feature_name = c("AQP4", "AQP4", "AQP4", "COL4A1", "COL4A1", "GFAP", "GFAP", "GFAP"),
        x_location = c(16, 60, 104, 35, 132, 22, 76, 126),
        y_location = c(44, 48, 46, 84, 82, 128, 130, 124),
        source_type = c("support_superpixel", "support_superpixel", "support_superpixel", "support_contour_anchor", "support_contour_anchor", "support_grid_point", "support_grid_point", "support_grid_point"),
        support_score = c(0.9, 0.8, 0.85, 0.7, 0.75, 0.4, 0.45, 0.42),
        support_area_um2 = c(50, 60, 55, 20, 20, 25, 30, 28),
        cell_id = c("c1", "c2", "c3", "c1", "c4", "c1", "c2", "c4"),
        assignment_confidence = c("approximate", "approximate", "approximate", "boundary-associated", "boundary-associated", "approximate", "approximate", "approximate"),
        assignment_status = c("interior_support", "interior_support", "interior_support", "boundary", "boundary", "interior_support", "interior_support", "interior_support"),
        tile_id = "tile_001",
        fov_name = "fov_1",
        segmentation_version = "route_b_example",
        support_call_method = "synthetic",
        resolution_level_used = 1,
        valid_pixel_fraction_local = 1,
        is_recovered_region = FALSE,
        assignment_semantics = c("support_summary", "support_summary", "support_summary", "boundary_support", "boundary_support", "support_summary", "support_summary", "support_summary"),
        assignment_ambiguity = c("low", "low", "low", "medium", "medium", "low", "low", "low")
    )
    cells <- data.table(
        cell_id = paste0("c", 1:4),
        x_centroid = c(15, 55, 100, 145),
        y_centroid = c(20, 20, 24, 22),
        protein_event_counts = c(2, 2, 2, 2),
        approximate_event_count = c(2, 2, 2, 2),
        interior_support_count = c(2, 2, 1, 2),
        interior_support_score_sum = c(1.3, 1.25, 0.85, 1.17),
        boundary_support_count = c(1, 0, 0, 1),
        boundary_support_score_sum = c(0.7, 0, 0, 0.75),
        total_intensity = c(230, 182, 125, 258),
        cell_area = c(180, 190, 175, 200),
        nucleus_area = c(55, 58, 52, 60),
        dominant_marker = c("PDGFRA", "VIM", "SLC1A3", "CLDN5"),
        segmentation_version = "route_b_example",
        assignment_policy_version = "example_v1"
    )
    mk_poly <- function(id, cx, cy, half) {
        data.table(
            cell_id = id,
            vertex_x = c(cx - half, cx + half, cx + half, cx - half),
            vertex_y = c(cy - half, cy - half, cy + half, cy + half),
            vertex_index = 1:4,
            fov = "fov_1",
            label_id = id
        )
    }
    cell_boundaries <- rbindlist(list(
        mk_poly("c1", 15, 20, 10),
        mk_poly("c2", 55, 20, 10),
        mk_poly("c3", 100, 24, 10),
        mk_poly("c4", 145, 22, 10)
    ))
    nucleus_boundaries <- rbindlist(list(
        mk_poly("c1", 15, 20, 4),
        mk_poly("c2", 55, 20, 4),
        mk_poly("c3", 100, 24, 4),
        mk_poly("c4", 145, 22, 4)
    ))

    arrow::write_parquet(protein_events, file.path(out_dir, "molecules", "protein_events.parquet"))
    arrow::write_parquet(signal_supports, file.path(out_dir, "molecules", "signal_supports.parquet"))
    arrow::write_parquet(cells, file.path(out_dir, "cells", "cells.parquet"))
    arrow::write_parquet(cell_boundaries, file.path(out_dir, "segmentation", "cell_boundaries.parquet"))
    arrow::write_parquet(nucleus_boundaries, file.path(out_dir, "segmentation", "nucleus_boundaries.parquet"))
    jsonlite::write_json(
        list(
            primary_analytical_layers = list(
                protein_events = "discrete punctate/connected-component events",
                signal_supports = "diffuse or boundary support summaries"
            ),
            count_like_layers = list(protein_events = TRUE, signal_supports = FALSE)
        ),
        file.path(out_dir, "metadata", "data_semantics.json"),
        auto_unbox = TRUE,
        pretty = TRUE
    )
    jsonlite::write_json(
        list(dataset_id = "phenocycler_minimal_bundle", version = "v1"),
        file.path(out_dir, "metadata", "dataset_manifest.json"),
        auto_unbox = TRUE,
        pretty = TRUE
    )
    invisible(out_dir)
}

.pc_coalesce <- function(x, y) {
    if (is.null(x)) y else x
}
