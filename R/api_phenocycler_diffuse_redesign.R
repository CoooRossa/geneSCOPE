#' Create a diffuse-marker redesign config for PhenoCycler bundles.
#' @description
#' Defines family-specific candidate representations for the explicit diffuse
#' redesign entrypoints. This config is additive and does not alter any legacy
#' or v1 PhenoCycler defaults.
#' @param grid_lengths Numeric vector of grid lengths to materialise.
#' @param region_shape Proxy shape used for region-first support output.
#' @param overlap_mode Rasterisation mode used for diffuse support regions.
#' @param anchor_target_area_um2 Approximate support area represented by one
#'   anchor in the multi-anchor view.
#' @param max_anchors_per_support Maximum anchors emitted per support region.
#' @param boundary_segment_vertices Number of local cell-boundary vertices to
#'   retain per boundary-associated support.
#' @param boundary_anchor_count Number of sampled anchors emitted per boundary
#'   segment.
#' @param feature_route_overrides Optional named character vector mapping
#'   `feature_name -> primary_route`.
#' @param compare_markers Character vector of markers prioritised in validation.
#' @param min_weight Minimum retained support weight.
#' @param version Config version string stored in manifests.
#' @return A named list with class
#'   `genescope_phenocycler_diffuse_redesign_config`.
#' @export
makePhenoCyclerDiffuseRedesignConfig <- function(
    grid_lengths = c(20, 40),
    region_shape = c("square_proxy"),
    overlap_mode = c("area_overlap"),
    anchor_target_area_um2 = 2500,
    max_anchors_per_support = 9,
    boundary_segment_vertices = 7,
    boundary_anchor_count = 5,
    feature_route_overrides = NULL,
    compare_markers = c("AQP4", "TMEM119", "F4/80", "TUBB3", "p-MYPT"),
    min_weight = 0,
    version = "phenocycler_diffuse_redesign_config_v1"
) {
    region_shape <- match.arg(region_shape)
    overlap_mode <- match.arg(overlap_mode)
    cfg <- list(
        grid_lengths = sort(unique(as.numeric(grid_lengths))),
        region_shape = region_shape,
        overlap_mode = overlap_mode,
        anchor_target_area_um2 = as.numeric(anchor_target_area_um2)[1],
        max_anchors_per_support = as.integer(max_anchors_per_support)[1],
        boundary_segment_vertices = as.integer(boundary_segment_vertices)[1],
        boundary_anchor_count = as.integer(boundary_anchor_count)[1],
        feature_route_overrides = .pc_named_chr(feature_route_overrides),
        compare_markers = unique(as.character(compare_markers)),
        min_weight = as.numeric(min_weight)[1],
        family_primary_routes = c(
            punctate = "event",
            diffuse_region_like = "D3_grid_field_overlap",
            boundary_associated = "D4_boundary_segment_support"
        ),
        candidate_routes = c(
            D1 = "support_region_proxy",
            D2 = "multi_anchor_weighted_support",
            D3 = "grid_field_overlap",
            D4 = "boundary_segment_support",
            D5 = "family_specific_hybrid"
        ),
        version = as.character(version)[1]
    )
    class(cfg) <- c("genescope_phenocycler_diffuse_redesign_config", "list")
    cfg
}

#' Create a redesign adapter config for PhenoCycler diffuse representations.
#' @description
#' Defines how redesigned representations are attached to a `scope_object` and
#' optionally promoted into explicit grid namespaces.
#' @param grid_lengths Numeric vector of grid lengths to materialise.
#' @param promote Default promotion behaviour (`none`, `field`, `boundary`,
#'   `hybrid`, or `punctate`).
#' @param field_layer_prefix Prefix for redesigned diffuse field layers.
#' @param boundary_layer_prefix Prefix for redesigned boundary field layers.
#' @param hybrid_layer_prefix Prefix for redesigned hybrid layers.
#' @param punctate_layer_prefix Prefix for redesigned punctate layers.
#' @param namespace_hybrid_features Whether hybrid layers namespace features by
#'   family.
#' @param attach_density_tables Whether to attach wide density tables.
#' @param version Config version string stored in manifests.
#' @return A named list with class
#'   `genescope_phenocycler_redesign_adapter_config`.
#' @export
makePhenoCyclerRedesignAdapterConfig <- function(
    grid_lengths = c(20, 40),
    promote = c("none", "field", "boundary", "hybrid", "punctate"),
    field_layer_prefix = "pcdr_field",
    boundary_layer_prefix = "pcdr_boundary",
    hybrid_layer_prefix = "pcdr_hybrid",
    punctate_layer_prefix = "pcdr_punctate",
    namespace_hybrid_features = TRUE,
    attach_density_tables = TRUE,
    version = "phenocycler_redesign_adapter_config_v1"
) {
    promote <- match.arg(promote)
    cfg <- list(
        grid_lengths = sort(unique(as.numeric(grid_lengths))),
        promote = promote,
        field_layer_prefix = as.character(field_layer_prefix)[1],
        boundary_layer_prefix = as.character(boundary_layer_prefix)[1],
        hybrid_layer_prefix = as.character(hybrid_layer_prefix)[1],
        punctate_layer_prefix = as.character(punctate_layer_prefix)[1],
        namespace_hybrid_features = isTRUE(namespace_hybrid_features),
        attach_density_tables = isTRUE(attach_density_tables),
        version = as.character(version)[1]
    )
    class(cfg) <- c("genescope_phenocycler_redesign_adapter_config", "list")
    cfg
}

#' Ingest an organised PhenoCycler bundle with redesigned diffuse semantics.
#' @description
#' Materialises explicit D1-D5 candidate layers without changing the existing
#' `ingestPhenoCyclerBundle()` surface or defaults.
#' @param bundle_dir Path to the organised PhenoCycler bundle.
#' @param marker_policy Marker-family policy object.
#' @param redesign_config Diffuse redesign config.
#' @param adapter_config Redesign adapter config.
#' @param verbose Emit progress messages when TRUE.
#' @return A named list with class `genescope_phenocycler_redesign_bundle`.
#' @export
ingestPhenoCyclerBundleRedesign <- function(
    bundle_dir,
    marker_policy = makePhenoCyclerMarkerPolicy(),
    redesign_config = makePhenoCyclerDiffuseRedesignConfig(),
    adapter_config = makePhenoCyclerRedesignAdapterConfig(),
    verbose = TRUE
) {
    bundle_dir <- normalizePath(bundle_dir, mustWork = TRUE)
    .pc_require_pkg("arrow", "ingestPhenoCyclerBundleRedesign()")
    legacy_bundle <- ingestPhenoCyclerBundle(
        bundle_dir = bundle_dir,
        marker_policy = marker_policy,
        diffuse_config = makePhenoCyclerDiffuseConfig(
            grid_lengths = redesign_config$grid_lengths,
            min_weight = redesign_config$min_weight
        ),
        adapter_config = makePhenoCyclerAdapterConfig(
            grid_lengths = adapter_config$grid_lengths,
            promote = "none"
        ),
        verbose = FALSE
    )

    supports <- as.data.table(legacy_bundle$inputs$signal_supports)
    events <- as.data.table(legacy_bundle$inputs$protein_events)
    cells <- legacy_bundle$inputs$cells
    cell_boundaries <- legacy_bundle$inputs$cell_boundaries
    nucleus_boundaries <- legacy_bundle$inputs$nucleus_boundaries
    feature_meta <- .pcdr_feature_metadata(
        legacy_bundle$metadata$feature_metadata,
        redesign_config
    )

    d1 <- .pcdr_build_support_region_proxies(supports, redesign_config)
    d2 <- .pcdr_build_multi_anchors(d1$regions, redesign_config)
    d4 <- .pcdr_build_boundary_segments(
        supports = supports,
        cell_boundaries = cell_boundaries,
        redesign_config = redesign_config
    )

    bounds <- legacy_bundle$manifest$bounds
    grid_lengths <- sort(unique(c(redesign_config$grid_lengths, adapter_config$grid_lengths)))
    d3_grid_layers <- .pcdr_build_grid_layers(
        events = events,
        region_summary = d1$regions,
        boundary_anchors = d4$anchors,
        bounds = bounds,
        grid_lengths = grid_lengths,
        adapter_config = adapter_config
    )

    family_routes <- .pcdr_family_route_table(feature_meta, redesign_config)

    manifest <- list(
        bundle_dir = bundle_dir,
        legacy_marker_policy_version = legacy_bundle$marker_policy$version,
        redesign_config_version = redesign_config$version,
        redesign_adapter_config_version = adapter_config$version,
        event_rows = nrow(events),
        support_rows = nrow(supports),
        feature_count = nrow(feature_meta),
        grid_lengths = grid_lengths,
        bounds = bounds,
        d1_region_rows = nrow(d1$regions),
        d2_anchor_rows = nrow(d2),
        d4_segment_rows = nrow(d4$segments),
        selected_compare_markers = redesign_config$compare_markers,
        candidate_routes = unname(redesign_config$candidate_routes),
        family_primary_routes = unname(redesign_config$family_primary_routes),
        recommended_default_scheme = list(
            punctate = "legacy_event",
            diffuse = "D3_grid_field_overlap",
            boundary = "D4_boundary_segment_support",
            system = "D5_family_specific_hybrid"
        )
    )

    if (verbose) {
        message(
            "[geneSCOPE::ingestPhenoCyclerBundleRedesign] bundle=",
            basename(bundle_dir),
            " regions=", nrow(d1$regions),
            " anchors=", nrow(d2),
            " boundary_segments=", nrow(d4$segments),
            " grids=", paste(grid_lengths, collapse = ",")
        )
    }

    out <- list(
        manifest = manifest,
        marker_policy = marker_policy,
        redesign_config = redesign_config,
        adapter_config = adapter_config,
        inputs = list(
            protein_events = as.data.frame(events),
            signal_supports = as.data.frame(supports),
            cells = cells,
            cell_boundaries = cell_boundaries,
            nucleus_boundaries = nucleus_boundaries
        ),
        metadata = list(
            data_semantics = legacy_bundle$metadata$data_semantics,
            dataset_manifest = legacy_bundle$metadata$dataset_manifest,
            feature_metadata = feature_meta,
            family_routes = family_routes
        ),
        legacy = list(
            manifest = legacy_bundle$manifest,
            representations = legacy_bundle$representations
        ),
        representations = list(
            D1_support_regions = d1,
            D2_multi_anchors = as.data.frame(d2),
            D3_grid_layers = d3_grid_layers,
            D4_boundary_segments = d4,
            D5_hybrid = list(
                family_routes = family_routes,
                primary_diffuse = "D3_grid_field_overlap",
                primary_boundary = "D4_boundary_segment_support",
                primary_punctate = "legacy_event"
            )
        )
    )
    class(out) <- c("genescope_phenocycler_redesign_bundle", "list")
    out
}

#' Create a scope object from an organised PhenoCycler bundle using redesign
#' entrypoints.
#' @description
#' This constructor is explicit opt-in and does not modify the v1
#' `createSCOPE_phenocycler()` path.
#' @param bundle_dir Path to the organised PhenoCycler bundle.
#' @param marker_policy Marker-family policy object.
#' @param redesign_config Diffuse redesign config.
#' @param adapter_config Redesign adapter config.
#' @param promote Which redesigned representation, if any, should be promoted
#'   into `@grid`.
#' @param verbose Emit progress messages when TRUE.
#' @return A `scope_object`.
#' @export
createSCOPE_phenocycler_redesign <- function(
    bundle_dir,
    marker_policy = makePhenoCyclerMarkerPolicy(),
    redesign_config = makePhenoCyclerDiffuseRedesignConfig(),
    adapter_config = makePhenoCyclerRedesignAdapterConfig(),
    promote = adapter_config$promote,
    verbose = TRUE
) {
    scope_obj <- methods::new(
        "scope_object",
        coord = list(),
        grid = list(),
        meta.data = data.frame(row.names = character()),
        cells = list(),
        stats = list(platform = "phenocycler_redesign", modalities = list()),
        density = list()
    )
    addPhenoCyclerRedesign(
        scope_obj = scope_obj,
        bundle_dir = bundle_dir,
        marker_policy = marker_policy,
        redesign_config = redesign_config,
        adapter_config = adapter_config,
        promote = promote,
        verbose = verbose
    )
}

#' Attach redesigned PhenoCycler representations to an existing scope object.
#' @description
#' Adds redesigned diffuse representations to a sidecar namespace and keeps
#' promotion explicit.
#' @param scope_obj A `scope_object`.
#' @param bundle_dir Path to the organised PhenoCycler bundle.
#' @param redesign_bundle Optional already-ingested redesign bundle.
#' @param marker_policy Marker-family policy object.
#' @param redesign_config Diffuse redesign config.
#' @param adapter_config Redesign adapter config.
#' @param promote Which redesigned representation, if any, to promote.
#' @param overwrite_existing_modality Whether to overwrite an existing redesign
#'   modality sidecar.
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object`.
#' @export
addPhenoCyclerRedesign <- function(
    scope_obj,
    bundle_dir = NULL,
    redesign_bundle = NULL,
    marker_policy = makePhenoCyclerMarkerPolicy(),
    redesign_config = makePhenoCyclerDiffuseRedesignConfig(),
    adapter_config = makePhenoCyclerRedesignAdapterConfig(),
    promote = adapter_config$promote,
    overwrite_existing_modality = FALSE,
    verbose = TRUE
) {
    if (!inherits(scope_obj, "scope_object")) {
        stop("`scope_obj` must be a scope_object.")
    }
    if (is.null(redesign_bundle)) {
        if (is.null(bundle_dir)) stop("Provide either `bundle_dir` or `redesign_bundle`.")
        redesign_bundle <- ingestPhenoCyclerBundleRedesign(
            bundle_dir = bundle_dir,
            marker_policy = marker_policy,
            redesign_config = redesign_config,
            adapter_config = adapter_config,
            verbose = verbose
        )
    }
    if (!is.list(scope_obj@stats$modalities)) {
        scope_obj@stats$modalities <- list()
    }
    if (!overwrite_existing_modality &&
        !is.null(scope_obj@stats$modalities$phenocycler_redesign)) {
        stop("scope_obj already has a `phenocycler_redesign` modality sidecar. Set `overwrite_existing_modality=TRUE` to replace it.")
    }

    scope_obj@stats$modalities$phenocycler_redesign <- list(
        manifest = redesign_bundle$manifest,
        marker_policy = redesign_bundle$marker_policy,
        redesign_config = redesign_bundle$redesign_config,
        adapter_config = redesign_bundle$adapter_config,
        feature_metadata = redesign_bundle$metadata$feature_metadata,
        family_routes = redesign_bundle$metadata$family_routes,
        grid_layers = redesign_bundle$representations$D3_grid_layers,
        D1_support_regions = redesign_bundle$representations$D1_support_regions,
        D4_boundary_segments = redesign_bundle$representations$D4_boundary_segments
    )

    scope_obj@coord$phenocycler_redesign_events <- redesign_bundle$inputs$protein_events
    scope_obj@coord$phenocycler_redesign_supports <- redesign_bundle$inputs$signal_supports
    scope_obj@coord$phenocycler_redesign_d1_regions <- redesign_bundle$representations$D1_support_regions$regions
    scope_obj@coord$phenocycler_redesign_d1_region_vertices <- redesign_bundle$representations$D1_support_regions$vertices
    scope_obj@coord$phenocycler_redesign_d2_multi_anchors <- redesign_bundle$representations$D2_multi_anchors
    scope_obj@coord$phenocycler_redesign_d4_segments <- redesign_bundle$representations$D4_boundary_segments$segments
    scope_obj@coord$phenocycler_redesign_d4_segment_vertices <- redesign_bundle$representations$D4_boundary_segments$vertices
    scope_obj@coord$phenocycler_redesign_d4_segment_anchors <- redesign_bundle$representations$D4_boundary_segments$anchors
    scope_obj@coord$phenocycler_redesign_cell_centroids <- redesign_bundle$inputs$cells
    scope_obj@coord$phenocycler_redesign_segmentation_cell <- redesign_bundle$inputs$cell_boundaries
    scope_obj@coord$phenocycler_redesign_segmentation_nucleus <- redesign_bundle$inputs$nucleus_boundaries

    if (isTRUE(redesign_bundle$adapter_config$attach_density_tables)) {
        density_tables <- .pcdr_density_tables_from_grid_layers(
            redesign_bundle$representations$D3_grid_layers
        )
        for (nm in names(density_tables)) {
            scope_obj@density[[nm]] <- density_tables[[nm]]
        }
    }

    promote <- match.arg(
        as.character(promote)[1],
        c("none", "field", "boundary", "hybrid", "punctate")
    )
    if (!identical(promote, "none")) {
        scope_obj <- promotePhenoCyclerRedesignToGridLayer(
            scope_obj = scope_obj,
            source = promote,
            grid_lengths = redesign_bundle$adapter_config$grid_lengths,
            overwrite = FALSE,
            verbose = verbose
        )
    }

    if (verbose) {
        message(
            "[geneSCOPE::addPhenoCyclerRedesign] attached redesign modality with ",
            nrow(redesign_bundle$metadata$feature_metadata),
            " features; promote=", promote
        )
    }
    scope_obj
}

#' Promote redesigned PhenoCycler representations into `@grid`.
#' @description
#' Copies an explicit redesigned representation into the standard `@grid`
#' namespace so downstream grid functions can consume it intentionally.
#' @param scope_obj A `scope_object`.
#' @param source Which redesigned representation to promote.
#' @param grid_lengths Optional subset of grid lengths to promote.
#' @param overwrite Whether to overwrite an existing grid layer.
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object`.
#' @export
promotePhenoCyclerRedesignToGridLayer <- function(
    scope_obj,
    source = c("field", "boundary", "hybrid", "punctate"),
    grid_lengths = NULL,
    overwrite = FALSE,
    verbose = TRUE
) {
    if (!inherits(scope_obj, "scope_object")) {
        stop("`scope_obj` must be a scope_object.")
    }
    source <- match.arg(source)
    modality <- scope_obj@stats$modalities$phenocycler_redesign
    if (is.null(modality)) {
        stop("No `phenocycler_redesign` modality is attached to this scope_object.")
    }

    grid_layers <- modality$grid_layers
    lengths_available <- as.numeric(sub("^grid", "", names(grid_layers)))
    if (is.null(grid_lengths)) {
        grid_lengths <- lengths_available
    }
    grid_lengths <- sort(unique(as.numeric(grid_lengths)))
    missing_lengths <- setdiff(grid_lengths, lengths_available)
    if (length(missing_lengths)) {
        stop(
            "Requested grid lengths not available in phenocycler redesign modality: ",
            paste(missing_lengths, collapse = ", ")
        )
    }

    for (lg in grid_lengths) {
        key <- paste0("grid", .pc_fmt_number(lg))
        layer_bundle <- grid_layers[[key]]
        if (is.null(layer_bundle[[source]])) {
            stop("No `", source, "` layer available for ", key, ".")
        }
        target_name <- .pcdr_promoted_layer_name(modality$adapter_config, source, lg)
        if (!overwrite && !is.null(scope_obj@grid[[target_name]])) {
            stop("Grid layer `", target_name, "` already exists. Set `overwrite=TRUE` to replace it.")
        }
        promoted_layer <- layer_bundle[[source]]
        density_table <- NULL
        if (isTRUE(modality$adapter_config$attach_density_tables)) {
            density_table <- layer_bundle$density[[source]]
        }
        feature_meta <- .pcdr_promoted_feature_metadata(
            modality$feature_metadata,
            source,
            modality$adapter_config
        )
        prepared <- .pc_prepare_promoted_feature_namespace(
            scope_obj = scope_obj,
            promoted_layer = promoted_layer,
            density_table = density_table,
            feature_metadata = feature_meta,
            namespace_prefix = "pcdr::"
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
                "[geneSCOPE::promotePhenoCyclerRedesignToGridLayer] ",
                source, " -> @grid$", target_name,
                " (", nrow(layer_bundle[[source]]$grid_info), " tiles",
                if (length(prepared$collisions)) paste0("; namespaced ", length(prepared$collisions), " colliding feature(s)") else "",
                ")"
            )
        }
    }
    scope_obj
}

#' Compare legacy and redesigned diffuse PhenoCycler representations.
#' @description
#' Builds marker-level and ROI-level comparison artefacts without modifying any
#' legacy tables or overlays.
#' @param bundle_dir Path to the organised PhenoCycler bundle.
#' @param marker_policy Marker-family policy object.
#' @param redesign_config Diffuse redesign config.
#' @param adapter_config Redesign adapter config.
#' @param selected_markers Optional marker subset for validation artefacts.
#' @param output_dir Directory where comparison artefacts should be written.
#' @param verbose Emit progress messages when TRUE.
#' @return A named list describing the comparison and written artefacts.
#' @export
comparePhenoCyclerDiffuseRepresentations <- function(
    bundle_dir,
    marker_policy = makePhenoCyclerMarkerPolicy(),
    redesign_config = makePhenoCyclerDiffuseRedesignConfig(),
    adapter_config = makePhenoCyclerRedesignAdapterConfig(),
    selected_markers = redesign_config$compare_markers,
    output_dir = tempdir(),
    verbose = TRUE
) {
    .pc_require_pkg("jsonlite", "comparePhenoCyclerDiffuseRepresentations()")
    .pc_require_pkg("ggplot2", "comparePhenoCyclerDiffuseRepresentations()")
    bundle_dir <- normalizePath(bundle_dir, mustWork = TRUE)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    legacy_bundle <- ingestPhenoCyclerBundle(
        bundle_dir = bundle_dir,
        marker_policy = marker_policy,
        diffuse_config = makePhenoCyclerDiffuseConfig(
            grid_lengths = redesign_config$grid_lengths,
            min_weight = redesign_config$min_weight
        ),
        adapter_config = makePhenoCyclerAdapterConfig(
            grid_lengths = adapter_config$grid_lengths,
            promote = "none"
        ),
        verbose = FALSE
    )
    redesign_bundle <- ingestPhenoCyclerBundleRedesign(
        bundle_dir = bundle_dir,
        marker_policy = marker_policy,
        redesign_config = redesign_config,
        adapter_config = adapter_config,
        verbose = FALSE
    )

    selected_markers <- unique(as.character(selected_markers))
    marker_level <- .pcdr_marker_level_validation(
        bundle_dir = bundle_dir,
        legacy_bundle = legacy_bundle,
        redesign_bundle = redesign_bundle,
        selected_markers = selected_markers
    )
    overlay_compare <- .pcdr_overlay_compare_manifest(
        bundle_dir = bundle_dir,
        legacy_bundle = legacy_bundle,
        redesign_bundle = redesign_bundle,
        selected_markers = selected_markers,
        output_dir = output_dir
    )
    candidate_compare <- .pcdr_candidate_compare_table(marker_level)
    case_studies <- .pcdr_case_study_table(marker_level, candidate_compare)
    modeling_readiness <- .pcdr_modeling_readiness_table()

    marker_path <- file.path(output_dir, "marker_level_validation.tsv")
    overlay_path <- file.path(output_dir, "overlay_compare_manifest.tsv")
    candidate_path <- file.path(output_dir, "representation_candidate_compare.tsv")
    case_study_path <- file.path(output_dir, "selected_marker_case_studies.tsv")
    modeling_path <- file.path(output_dir, "modeling_readiness_summary.tsv")
    summary_path <- file.path(output_dir, "representation_compare_summary.json")
    data.table::fwrite(marker_level, marker_path, sep = "\t")
    data.table::fwrite(overlay_compare$manifest, overlay_path, sep = "\t")
    data.table::fwrite(candidate_compare, candidate_path, sep = "\t")
    data.table::fwrite(case_studies, case_study_path, sep = "\t")
    data.table::fwrite(modeling_readiness, modeling_path, sep = "\t")
    jsonlite::write_json(
        list(
            bundle_dir = bundle_dir,
            selected_markers = selected_markers,
            marker_rows = nrow(marker_level),
            overlay_rows = nrow(overlay_compare$manifest),
            candidate_rows = nrow(candidate_compare),
            case_study_rows = nrow(case_studies),
            modeling_rows = nrow(modeling_readiness),
            plot_dir = normalizePath(overlay_compare$plot_dir, mustWork = TRUE)
        ),
        summary_path,
        auto_unbox = TRUE,
        pretty = TRUE
    )

    if (verbose) {
        message(
            "[geneSCOPE::comparePhenoCyclerDiffuseRepresentations] markers=",
            nrow(marker_level),
            " overlays=", nrow(overlay_compare$manifest),
            " output_dir=", normalizePath(output_dir, mustWork = TRUE)
        )
    }

    list(
        marker_level = marker_level,
        overlay_compare = overlay_compare$manifest,
        candidate_compare = candidate_compare,
        case_studies = case_studies,
        modeling_readiness = modeling_readiness,
        paths = list(
            marker_level_tsv = marker_path,
            overlay_compare_tsv = overlay_path,
            candidate_compare_tsv = candidate_path,
            case_studies_tsv = case_study_path,
            modeling_readiness_tsv = modeling_path,
            summary_json = summary_path,
            plot_dir = overlay_compare$plot_dir
        )
    )
}

.pcdr_feature_metadata <- function(feature_metadata, redesign_config) {
    df <- as.data.frame(feature_metadata)
    if (!nrow(df)) return(df)
    df$primary_route <- .pcdr_primary_route_for_feature(
        feature_name = df$feature_name,
        marker_family = df$marker_family,
        redesign_config = redesign_config
    )
    df$modeling_primary_layer <- ifelse(
        df$marker_family == "boundary_associated",
        "D4_boundary_segment_support",
        ifelse(df$marker_family == "punctate", "event", "D3_grid_field_overlap")
    )
    rownames(df) <- df$feature_name
    df
}

.pcdr_primary_route_for_feature <- function(feature_name, marker_family, redesign_config) {
    out <- unname(redesign_config$family_primary_routes[as.character(marker_family)])
    out[is.na(out) | !nzchar(out)] <- "D3_grid_field_overlap"
    overrides <- redesign_config$feature_route_overrides
    if (length(overrides)) {
        idx <- match(as.character(feature_name), names(overrides))
        keep <- !is.na(idx)
        out[keep] <- unname(overrides[idx[keep]])
    }
    out
}

.pcdr_family_route_table <- function(feature_metadata, redesign_config) {
    if (is.null(feature_metadata) || !nrow(feature_metadata)) {
        return(data.frame(
            feature_name = character(),
            marker_family = character(),
            primary_route = character(),
            modeling_primary_layer = character(),
            stringsAsFactors = FALSE
        ))
    }
    out <- as.data.frame(feature_metadata)[, c(
        "feature_name",
        "marker_family",
        "primary_route",
        "modeling_primary_layer"
    ), drop = FALSE]
    out$hybrid_route <- ifelse(
        out$marker_family == "punctate",
        "legacy_event",
        ifelse(
            out$marker_family == "boundary_associated",
            "D4_boundary_segment_support",
            "D3_grid_field_overlap"
        )
    )
    out
}

.pcdr_build_support_region_proxies <- function(supports, redesign_config) {
    if (is.null(supports) || !nrow(supports)) {
        empty_regions <- data.frame(
            region_id = character(),
            support_id = character(),
            feature_name = character(),
            marker_family = character(),
            source_type = character(),
            x_center = numeric(),
            y_center = numeric(),
            support_weight = numeric(),
            local_scale = numeric(),
            support_area_um2 = numeric(),
            area_proxy_um2 = numeric(),
            proxy_half_side = numeric(),
            bbox_xmin = numeric(),
            bbox_xmax = numeric(),
            bbox_ymin = numeric(),
            bbox_ymax = numeric(),
            representation_candidate = character(),
            representation_role = character(),
            stringsAsFactors = FALSE
        )
        empty_vertices <- data.frame(
            region_id = character(),
            support_id = character(),
            feature_name = character(),
            marker_family = character(),
            vertex_order = integer(),
            vertex_x = numeric(),
            vertex_y = numeric(),
            stringsAsFactors = FALSE
        )
        return(list(regions = empty_regions, vertices = empty_vertices))
    }
    dt <- as.data.table(supports)
    dt <- dt[support_weight > redesign_config$min_weight]
    if (!nrow(dt)) {
        return(.pcdr_build_support_region_proxies(dt, redesign_config))
    }
    dt[, region_id := as.character(.pc_val_or_seq(support_id, "pcdr_region"))]
    dt[, proxy_half_side := fifelse(
        is.finite(support_area_um2) & support_area_um2 > 0,
        sqrt(support_area_um2) / 2,
        fifelse(is.finite(local_scale) & local_scale > 0, local_scale, 0.5)
    )]
    dt[, area_proxy_um2 := (proxy_half_side * 2) ^ 2]
    dt[, bbox_xmin := x_location - proxy_half_side]
    dt[, bbox_xmax := x_location + proxy_half_side]
    dt[, bbox_ymin := y_location - proxy_half_side]
    dt[, bbox_ymax := y_location + proxy_half_side]
    dt[, representation_candidate := fifelse(
        marker_family == "boundary_associated",
        "D4_boundary_proxy_seed",
        "D1_support_region_proxy"
    )]
    dt[, representation_role := fifelse(
        marker_family == "diffuse_region_like",
        "candidate_primary",
        "candidate_auxiliary"
    )]
    regions <- dt[, .(
        region_id,
        support_id = as.character(support_id),
        feature_name = as.character(feature_name),
        marker_family = as.character(marker_family),
        source_type = as.character(source_type),
        x_center = as.numeric(x_location),
        y_center = as.numeric(y_location),
        support_weight = as.numeric(support_weight),
        local_scale = as.numeric(local_scale),
        support_area_um2 = as.numeric(support_area_um2),
        area_proxy_um2 = as.numeric(area_proxy_um2),
        proxy_half_side = as.numeric(proxy_half_side),
        bbox_xmin = as.numeric(bbox_xmin),
        bbox_xmax = as.numeric(bbox_xmax),
        bbox_ymin = as.numeric(bbox_ymin),
        bbox_ymax = as.numeric(bbox_ymax),
        representation_candidate = as.character(representation_candidate),
        representation_role = as.character(representation_role)
    )]
    vertices <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
        rr <- regions[i]
        data.table(
            region_id = rr$region_id,
            support_id = rr$support_id,
            feature_name = rr$feature_name,
            marker_family = rr$marker_family,
            vertex_order = 1:5,
            vertex_x = c(rr$bbox_xmin, rr$bbox_xmax, rr$bbox_xmax, rr$bbox_xmin, rr$bbox_xmin),
            vertex_y = c(rr$bbox_ymin, rr$bbox_ymin, rr$bbox_ymax, rr$bbox_ymax, rr$bbox_ymin)
        )
    }), use.names = TRUE, fill = TRUE)
    list(
        regions = as.data.frame(regions),
        vertices = as.data.frame(vertices)
    )
}

.pcdr_anchor_offsets <- function(n) {
    n <- max(1L, min(as.integer(n), 9L))
    base <- matrix(
        c(
            0, 0,
            -0.65, 0,
            0.65, 0,
            0, -0.65,
            0, 0.65,
            -0.55, -0.55,
            0.55, -0.55,
            -0.55, 0.55,
            0.55, 0.55
        ),
        byrow = TRUE,
        ncol = 2
    )
    base[seq_len(n), , drop = FALSE]
}

.pcdr_build_multi_anchors <- function(region_summary, redesign_config) {
    if (is.null(region_summary) || !nrow(region_summary)) {
        return(data.frame(
            anchor_id = character(),
            region_id = character(),
            support_id = character(),
            feature_name = character(),
            marker_family = character(),
            source_type = character(),
            x_location = numeric(),
            y_location = numeric(),
            weight = numeric(),
            local_scale = numeric(),
            anchor_rank = integer(),
            anchor_count = integer(),
            representation_type = character(),
            stringsAsFactors = FALSE
        ))
    }
    dt <- as.data.table(region_summary)
    dt <- dt[marker_family == "diffuse_region_like"]
    if (!nrow(dt)) {
        return(data.frame(
            anchor_id = character(),
            region_id = character(),
            support_id = character(),
            feature_name = character(),
            marker_family = character(),
            source_type = character(),
            x_location = numeric(),
            y_location = numeric(),
            weight = numeric(),
            local_scale = numeric(),
            anchor_rank = integer(),
            anchor_count = integer(),
            representation_type = character(),
            stringsAsFactors = FALSE
        ))
    }
    dt[, anchor_count := pmax(
        1L,
        pmin(
            redesign_config$max_anchors_per_support,
            ceiling(area_proxy_um2 / redesign_config$anchor_target_area_um2)
        )
    )]
    out <- rbindlist(lapply(seq_len(nrow(dt)), function(i) {
        rr <- dt[i]
        offsets <- .pcdr_anchor_offsets(rr$anchor_count)
        data.table(
            anchor_id = paste0(rr$region_id, "_a", seq_len(nrow(offsets))),
            region_id = rr$region_id,
            support_id = rr$support_id,
            feature_name = rr$feature_name,
            marker_family = rr$marker_family,
            source_type = rr$source_type,
            x_location = rr$x_center + offsets[, 1] * rr$proxy_half_side,
            y_location = rr$y_center + offsets[, 2] * rr$proxy_half_side,
            weight = rr$support_weight / nrow(offsets),
            local_scale = rr$proxy_half_side,
            anchor_rank = seq_len(nrow(offsets)),
            anchor_count = nrow(offsets),
            representation_type = "D2_multi_anchor_weighted_support"
        )
    }), use.names = TRUE, fill = TRUE)
    as.data.frame(out)
}

.pcdr_build_boundary_segments <- function(supports, cell_boundaries, redesign_config) {
    empty_segments <- data.frame(
        segment_id = character(),
        support_id = character(),
        feature_name = character(),
        marker_family = character(),
        source_type = character(),
        cell_id = character(),
        x_center = numeric(),
        y_center = numeric(),
        support_weight = numeric(),
        segment_length = numeric(),
        matched_cell_boundary = logical(),
        representation_type = character(),
        stringsAsFactors = FALSE
    )
    empty_vertices <- data.frame(
        segment_id = character(),
        support_id = character(),
        feature_name = character(),
        marker_family = character(),
        vertex_order = integer(),
        vertex_x = numeric(),
        vertex_y = numeric(),
        stringsAsFactors = FALSE
    )
    empty_anchors <- data.frame(
        anchor_id = character(),
        segment_id = character(),
        support_id = character(),
        feature_name = character(),
        marker_family = character(),
        source_type = character(),
        x_location = numeric(),
        y_location = numeric(),
        weight = numeric(),
        local_scale = numeric(),
        representation_type = character(),
        stringsAsFactors = FALSE
    )
    if (is.null(supports) || !nrow(supports)) {
        return(list(segments = empty_segments, vertices = empty_vertices, anchors = empty_anchors))
    }
    supports_dt <- as.data.table(supports)
    supports_dt <- supports_dt[marker_family == "boundary_associated" & support_weight > redesign_config$min_weight]
    if (!nrow(supports_dt)) {
        return(list(segments = empty_segments, vertices = empty_vertices, anchors = empty_anchors))
    }

    boundary_map <- list()
    if (!is.null(cell_boundaries) && nrow(cell_boundaries)) {
        cb <- as.data.table(cell_boundaries)
        cb <- cb[cell_id %in% unique(supports_dt$cell_id)]
        if (nrow(cb)) {
            cb <- cb[order(cell_id, vertex_index)]
            boundary_map <- split(cb[, .(vertex_x, vertex_y)], cb$cell_id)
        }
    }
    half_window <- max(1L, floor(redesign_config$boundary_segment_vertices / 2))
    out <- lapply(seq_len(nrow(supports_dt)), function(i) {
        rr <- supports_dt[i]
        poly <- boundary_map[[as.character(rr$cell_id)]]
        seg <- NULL
        matched <- FALSE
        if (!is.null(poly) && nrow(poly)) {
            px <- as.numeric(poly$vertex_x)
            py <- as.numeric(poly$vertex_y)
            idx0 <- which.min((px - rr$x_location) ^ 2 + (py - rr$y_location) ^ 2)
            ids <- ((idx0 - half_window - 1):(idx0 + half_window - 1)) %% nrow(poly) + 1L
            seg <- data.table(
                vertex_x = px[ids],
                vertex_y = py[ids]
            )
            matched <- TRUE
        } else {
            span <- max(rr$local_scale, sqrt(.pc_coalesce(rr$support_area_um2, 1)) / 2, 1)
            seg <- data.table(
                vertex_x = c(rr$x_location - span, rr$x_location + span),
                vertex_y = c(rr$y_location, rr$y_location)
            )
        }
        samples <- .pcdr_sample_path(seg, redesign_config$boundary_anchor_count)
        seg_len <- .pcdr_path_length(seg)
        list(
            segment = data.table(
                segment_id = paste0("segment_", rr$support_id),
                support_id = as.character(rr$support_id),
                feature_name = as.character(rr$feature_name),
                marker_family = as.character(rr$marker_family),
                source_type = as.character(rr$source_type),
                cell_id = as.character(rr$cell_id),
                x_center = as.numeric(rr$x_location),
                y_center = as.numeric(rr$y_location),
                support_weight = as.numeric(rr$support_weight),
                segment_length = as.numeric(seg_len),
                matched_cell_boundary = matched,
                representation_type = "D4_boundary_segment_support"
            ),
            vertices = data.table(
                segment_id = paste0("segment_", rr$support_id),
                support_id = as.character(rr$support_id),
                feature_name = as.character(rr$feature_name),
                marker_family = as.character(rr$marker_family),
                vertex_order = seq_len(nrow(seg)),
                vertex_x = as.numeric(seg$vertex_x),
                vertex_y = as.numeric(seg$vertex_y)
            ),
            anchors = data.table(
                anchor_id = paste0("segment_", rr$support_id, "_b", seq_len(nrow(samples))),
                segment_id = paste0("segment_", rr$support_id),
                support_id = as.character(rr$support_id),
                feature_name = as.character(rr$feature_name),
                marker_family = as.character(rr$marker_family),
                source_type = as.character(rr$source_type),
                x_location = as.numeric(samples$vertex_x),
                y_location = as.numeric(samples$vertex_y),
                weight = as.numeric(rr$support_weight) / nrow(samples),
                local_scale = max(rr$local_scale, 1),
                representation_type = "D4_boundary_segment_anchor"
            )
        )
    })
    list(
        segments = as.data.frame(rbindlist(lapply(out, `[[`, "segment"), use.names = TRUE, fill = TRUE)),
        vertices = as.data.frame(rbindlist(lapply(out, `[[`, "vertices"), use.names = TRUE, fill = TRUE)),
        anchors = as.data.frame(rbindlist(lapply(out, `[[`, "anchors"), use.names = TRUE, fill = TRUE))
    )
}

.pcdr_path_length <- function(seg) {
    if (is.null(seg) || nrow(seg) < 2) return(0)
    dx <- diff(as.numeric(seg$vertex_x))
    dy <- diff(as.numeric(seg$vertex_y))
    sum(sqrt(dx ^ 2 + dy ^ 2))
}

.pcdr_sample_path <- function(seg, n) {
    if (is.null(seg) || !nrow(seg)) {
        return(data.table(vertex_x = numeric(), vertex_y = numeric()))
    }
    if (nrow(seg) == 1L || n <= 1L) {
        return(data.table(vertex_x = seg$vertex_x[[1]], vertex_y = seg$vertex_y[[1]]))
    }
    dx <- diff(as.numeric(seg$vertex_x))
    dy <- diff(as.numeric(seg$vertex_y))
    lens <- sqrt(dx ^ 2 + dy ^ 2)
    cumlen <- c(0, cumsum(lens))
    total <- tail(cumlen, 1)
    if (!is.finite(total) || total <= 0) {
        return(data.table(
            vertex_x = rep(seg$vertex_x[[1]], n),
            vertex_y = rep(seg$vertex_y[[1]], n)
        ))
    }
    targets <- seq(0, total, length.out = max(2L, as.integer(n)))
    out <- vector("list", length(targets))
    for (i in seq_along(targets)) {
        tt <- targets[[i]]
        idx <- max(which(cumlen <= tt))
        idx <- min(idx, length(lens))
        frac <- if (lens[[idx]] > 0) (tt - cumlen[[idx]]) / lens[[idx]] else 0
        x0 <- seg$vertex_x[[idx]]
        y0 <- seg$vertex_y[[idx]]
        x1 <- seg$vertex_x[[idx + 1L]]
        y1 <- seg$vertex_y[[idx + 1L]]
        out[[i]] <- data.table(
            vertex_x = x0 + frac * (x1 - x0),
            vertex_y = y0 + frac * (y1 - y0)
        )
    }
    rbindlist(out, use.names = TRUE, fill = TRUE)
}

.pcdr_build_grid_layers <- function(events,
                                    region_summary,
                                    boundary_anchors,
                                    bounds,
                                    grid_lengths,
                                    adapter_config) {
    out <- vector("list", length(grid_lengths))
    names(out) <- paste0("grid", vapply(grid_lengths, .pc_fmt_number, character(1)))
    punctate_dt <- as.data.table(events)
    region_dt <- as.data.table(region_summary)
    boundary_dt <- as.data.table(boundary_anchors)

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
        diffuse_counts <- .pcdr_region_overlap_to_grid(
            region_dt[marker_family == "diffuse_region_like"],
            grid_info,
            lg
        )
        boundary_counts <- .pc_aggregate_to_grid(
            dt = boundary_dt,
            grid_info = grid_info,
            grid_length = lg,
            x_col = "x_location",
            y_col = "y_location",
            feature_col = "feature_name",
            value_col = "weight"
        )
        hybrid_counts <- .pcdr_build_hybrid_counts(
            punctate_counts = punctate_counts,
            diffuse_counts = diffuse_counts,
            boundary_counts = boundary_counts,
            namespace_hybrid_features = isTRUE(adapter_config$namespace_hybrid_features)
        )
        out[[paste0("grid", .pc_fmt_number(lg))]] <- list(
            punctate = .pc_make_grid_layer(grid_info, punctate_counts, lg, "legacy_punctate_event_layer"),
            field = .pc_make_grid_layer(grid_info, diffuse_counts, lg, "D3_grid_field_overlap"),
            boundary = .pc_make_grid_layer(grid_info, boundary_counts, lg, "D4_boundary_segment_field"),
            hybrid = .pc_make_grid_layer(grid_info, hybrid_counts, lg, "D5_family_specific_hybrid"),
            density = list(
                punctate = .pc_counts_to_density_table(grid_info, punctate_counts),
                field = .pc_counts_to_density_table(grid_info, diffuse_counts),
                boundary = .pc_counts_to_density_table(grid_info, boundary_counts),
                hybrid = .pc_counts_to_density_table(grid_info, hybrid_counts)
            )
        )
    }
    out
}

.pcdr_region_overlap_to_grid <- function(region_dt, grid_info, grid_length) {
    if (is.null(region_dt) || !nrow(region_dt)) {
        return(data.table(grid_id = character(), gene = character(), count = numeric()))
    }
    region_dt <- as.data.table(region_dt)
    grid_info <- as.data.table(grid_info)
    x0 <- min(grid_info$xmin)
    y0 <- min(grid_info$ymin)
    xbins_eff <- max(grid_info$gx)
    ybins_eff <- max(grid_info$gy)

    work <- copy(region_dt)
    work <- work[is.finite(bbox_xmin) & is.finite(bbox_xmax) & is.finite(bbox_ymin) &
        is.finite(bbox_ymax) & is.finite(area_proxy_um2) & area_proxy_um2 > 0]
    if (!nrow(work)) {
        return(data.table(grid_id = character(), gene = character(), count = numeric()))
    }
    work[, gx_min := pmax(1L, floor((bbox_xmin - x0) / grid_length) + 1L)]
    work[, gx_max := pmin(xbins_eff, floor(((bbox_xmax - 1e-9) - x0) / grid_length) + 1L)]
    work[, gy_min := pmax(1L, floor((bbox_ymin - y0) / grid_length) + 1L)]
    work[, gy_max := pmin(ybins_eff, floor(((bbox_ymax - 1e-9) - y0) / grid_length) + 1L)]
    work <- work[gx_min <= gx_max & gy_min <= gy_max]
    if (!nrow(work)) {
        return(data.table(grid_id = character(), gene = character(), count = numeric()))
    }

    single <- work[gx_min == gx_max & gy_min == gy_max, .(
        grid_id = paste0("g", gx_min, "_", gy_min),
        gene = as.character(feature_name),
        count = as.numeric(support_weight)
    )]
    multi <- work[gx_min != gx_max | gy_min != gy_max]
    multi_parts <- if (nrow(multi)) {
        rbindlist(lapply(seq_len(nrow(multi)), function(i) {
            rr <- multi[i]
            cells <- CJ(gx = rr$gx_min:rr$gx_max, gy = rr$gy_min:rr$gy_max, sorted = TRUE)
            cells[, xmin := x0 + (gx - 1L) * grid_length]
            cells[, xmax := xmin + grid_length]
            cells[, ymin := y0 + (gy - 1L) * grid_length]
            cells[, ymax := ymin + grid_length]
            cells[, overlap_x := pmax(0, pmin(rr$bbox_xmax, xmax) - pmax(rr$bbox_xmin, xmin))]
            cells[, overlap_y := pmax(0, pmin(rr$bbox_ymax, ymax) - pmax(rr$bbox_ymin, ymin))]
            cells[, overlap_area := overlap_x * overlap_y]
            cells <- cells[overlap_area > 0]
            if (!nrow(cells)) return(NULL)
            cells[, .(
                grid_id = paste0("g", gx, "_", gy),
                gene = rr$feature_name,
                count = rr$support_weight * overlap_area / rr$area_proxy_um2
            )]
        }), use.names = TRUE, fill = TRUE)
    } else {
        data.table(grid_id = character(), gene = character(), count = numeric())
    }
    combined <- rbindlist(list(single, multi_parts), use.names = TRUE, fill = TRUE)
    combined[, .(count = sum(as.numeric(count), na.rm = TRUE)), by = .(grid_id, gene)]
}

.pcdr_build_hybrid_counts <- function(punctate_counts,
                                      diffuse_counts,
                                      boundary_counts,
                                      namespace_hybrid_features = TRUE) {
    parts <- list()
    if (!is.null(punctate_counts) && nrow(punctate_counts)) {
        dt <- copy(as.data.table(punctate_counts))
        if (namespace_hybrid_features) dt[, gene := paste0("punctate::", gene)]
        parts[[length(parts) + 1L]] <- dt
    }
    if (!is.null(diffuse_counts) && nrow(diffuse_counts)) {
        dt <- copy(as.data.table(diffuse_counts))
        if (namespace_hybrid_features) dt[, gene := paste0("diffuse::", gene)]
        parts[[length(parts) + 1L]] <- dt
    }
    if (!is.null(boundary_counts) && nrow(boundary_counts)) {
        dt <- copy(as.data.table(boundary_counts))
        if (namespace_hybrid_features) dt[, gene := paste0("boundary::", gene)]
        parts[[length(parts) + 1L]] <- dt
    }
    if (!length(parts)) {
        return(data.table(grid_id = character(), gene = character(), count = numeric()))
    }
    merged <- rbindlist(parts, use.names = TRUE, fill = TRUE)
    merged[, .(count = sum(as.numeric(count), na.rm = TRUE)), by = .(grid_id, gene)]
}

.pcdr_density_tables_from_grid_layers <- function(grid_layers) {
    out <- list()
    for (key in names(grid_layers)) {
        layer_bundle <- grid_layers[[key]]
        for (src in names(layer_bundle$density)) {
            out[[paste0("phenocycler_redesign_", src, "_", key)]] <- layer_bundle$density[[src]]
        }
    }
    out
}

.pcdr_promoted_layer_name <- function(adapter_config, source, grid_length) {
    suffix <- .pc_fmt_number(grid_length)
    prefix <- switch(
        source,
        field = adapter_config$field_layer_prefix,
        boundary = adapter_config$boundary_layer_prefix,
        hybrid = adapter_config$hybrid_layer_prefix,
        punctate = adapter_config$punctate_layer_prefix,
        stop("Unknown source: ", source)
    )
    paste0(prefix, suffix)
}

.pcdr_promoted_feature_metadata <- function(feature_metadata, source, adapter_config) {
    if (is.null(feature_metadata) || !nrow(feature_metadata)) return(feature_metadata)
    df <- as.data.frame(feature_metadata)
    if (identical(source, "hybrid")) {
        parts <- list()
        for (fam in c("punctate", "diffuse_region_like", "boundary_associated")) {
            sub_df <- df[df$marker_family == fam, , drop = FALSE]
            if (!nrow(sub_df)) next
            prefix <- switch(
                fam,
                punctate = "punctate::",
                diffuse_region_like = "diffuse::",
                boundary_associated = "boundary::"
            )
            if (isTRUE(adapter_config$namespace_hybrid_features)) {
                rownames(sub_df) <- paste0(prefix, sub_df$feature_name)
                sub_df$feature_name <- rownames(sub_df)
            }
            sub_df$representation_default <- paste0("hybrid_", fam)
            parts[[length(parts) + 1L]] <- sub_df
        }
        if (!length(parts)) return(df[0, , drop = FALSE])
        return(do.call(rbind, parts))
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

.pcdr_marker_level_validation <- function(bundle_dir,
                                          legacy_bundle,
                                          redesign_bundle,
                                          selected_markers) {
    supports <- as.data.table(legacy_bundle$inputs$signal_supports)
    legacy_regions <- as.data.table(legacy_bundle$representations$support_regions)
    legacy_anchors <- as.data.table(legacy_bundle$representations$weighted_anchors)
    redesign_regions <- as.data.table(redesign_bundle$representations$D1_support_regions$regions)
    redesign_anchors <- as.data.table(redesign_bundle$representations$D2_multi_anchors)
    boundary_segments <- as.data.table(redesign_bundle$representations$D4_boundary_segments$segments)
    feature_meta <- as.data.table(redesign_bundle$metadata$feature_metadata)

    overlay_review_path <- file.path(bundle_dir, "qc", "protein_signal_overlay_review.tsv")
    overlay_priority_path <- file.path(bundle_dir, "qc", "protein_signal_overlay_priority_list.tsv")
    overlay_review <- if (file.exists(overlay_review_path)) {
        fread(overlay_review_path)
    } else {
        data.table(
            marker_name = character(),
            alignment_visual_score = numeric(),
            offset_risk_score = numeric(),
            structure_match_score = numeric(),
            representation_mismatch_risk = numeric()
        )
    }
    overlay_priority <- if (file.exists(overlay_priority_path)) {
        fread(overlay_priority_path)
    } else {
        data.table(
            marker_name = character(),
            priority_for_followup = character()
        )
    }

    legacy_field20 <- as.data.table(legacy_bundle$representations$grid_layers$grid20$field$counts)
    redesign_field20 <- as.data.table(redesign_bundle$representations$D3_grid_layers$grid20$field$counts)
    legacy_boundary20 <- as.data.table(legacy_bundle$representations$grid_layers$grid20$boundary$counts)
    redesign_boundary20 <- as.data.table(redesign_bundle$representations$D3_grid_layers$grid20$boundary$counts)

    features_keep <- feature_meta[
        marker_family %in% c("diffuse_region_like", "boundary_associated") |
            feature_name %in% selected_markers
    ]
    out <- rbindlist(lapply(seq_len(nrow(features_keep)), function(i) {
        ff <- features_keep[i]
        feature <- ff$feature_name
        family <- ff$marker_family
        old_s <- supports[feature_name == feature]
        old_r <- legacy_regions[feature_name == feature]
        old_a <- legacy_anchors[feature_name == feature]
        new_r <- redesign_regions[feature_name == feature]
        new_a <- redesign_anchors[feature_name == feature]
        new_b <- boundary_segments[feature_name == feature]
        review <- overlay_review[marker_name == feature]
        priority <- overlay_priority[marker_name == feature]
        old_field <- if (family == "boundary_associated") legacy_boundary20[gene == feature] else legacy_field20[gene == feature]
        new_field <- if (family == "boundary_associated") redesign_boundary20[gene == feature] else redesign_field20[gene == feature]
        data.table(
            marker_name = feature,
            marker_family = family,
            signal_type = .pcdr_scalar_or_blank(old_s$source_type),
            priority_for_followup = .pcdr_scalar_or_blank(priority$priority_for_followup),
            legacy_support_rows = nrow(old_s),
            legacy_weighted_anchor_rows = nrow(old_a),
            legacy_support_region_rows = nrow(old_r),
            redesign_d1_region_rows = nrow(new_r),
            redesign_d2_anchor_rows = nrow(new_a),
            redesign_d4_segment_rows = nrow(new_b),
            legacy_weight_sum = sum(as.numeric(old_s$support_weight), na.rm = TRUE),
            redesign_d2_weight_sum = sum(as.numeric(new_a$weight), na.rm = TRUE),
            redesign_d4_weight_sum = sum(as.numeric(new_b$support_weight), na.rm = TRUE),
            legacy_area_sum = sum(as.numeric(old_s$support_area_um2), na.rm = TRUE),
            redesign_d1_area_sum = sum(as.numeric(new_r$area_proxy_um2), na.rm = TRUE),
            legacy_grid20_nonzero_tiles = uniqueN(old_field$grid_id),
            redesign_grid20_nonzero_tiles = uniqueN(new_field$grid_id),
            legacy_grid20_mass = sum(as.numeric(old_field$count), na.rm = TRUE),
            redesign_grid20_mass = sum(as.numeric(new_field$count), na.rm = TRUE),
            review_rows = nrow(review),
            review_alignment_median = .pcdr_num_median(review$alignment_visual_score),
            review_offset_median = .pcdr_num_median(review$offset_risk_score),
            review_structure_median = .pcdr_num_median(review$structure_match_score),
            representation_mismatch_median = .pcdr_num_median(review$representation_mismatch_risk),
            primary_route = ff$primary_route,
            modeling_primary_layer = ff$modeling_primary_layer,
            validation_status = .pcdr_marker_validation_status(
                marker_family = family,
                review = review,
                legacy_rows = nrow(old_s),
                redesign_regions = nrow(new_r),
                redesign_anchors = nrow(new_a),
                redesign_segments = nrow(new_b)
            )
        )
    }), use.names = TRUE, fill = TRUE)
    out[order(-xtfrm(priority_for_followup), marker_family, marker_name)]
}

.pcdr_candidate_compare_table <- function(marker_level) {
    dt <- as.data.table(marker_level)
    if (!nrow(dt)) {
        return(data.frame(
            marker_name = character(),
            marker_family = character(),
            candidate_id = character(),
            primary_status = character(),
            usage_role = character(),
            keep_decision = character(),
            rationale = character(),
            stringsAsFactors = FALSE
        ))
    }
    rows <- unlist(lapply(seq_len(nrow(dt)), function(i) {
        rr <- dt[i]
        fam <- rr$marker_family
        list(
            data.table(
                marker_name = rr$marker_name,
                marker_family = fam,
                candidate_id = "legacy_single_anchor",
                primary_status = "legacy_reference_only",
                usage_role = "reference",
                keep_decision = "retained_for_compare_only",
                rationale = "Retained only as the pre-redesign reference layer."
            ),
            data.table(
                marker_name = rr$marker_name,
                marker_family = fam,
                candidate_id = "D1_support_region_proxy",
                primary_status = if (fam == "diffuse_region_like") "candidate" else "auxiliary",
                usage_role = "region_first",
                keep_decision = "keep",
                rationale = "Preserves support extent and avoids brightest-pixel misreading."
            ),
            data.table(
                marker_name = rr$marker_name,
                marker_family = fam,
                candidate_id = "D2_multi_anchor_weighted_support",
                primary_status = if (fam == "diffuse_region_like") "candidate" else "auxiliary",
                usage_role = "coordinate_first_modeling",
                keep_decision = "keep",
                rationale = "Provides coordinate-first density reconstruction without collapsing to a single point."
            ),
            data.table(
                marker_name = rr$marker_name,
                marker_family = fam,
                candidate_id = "D3_grid_field_overlap",
                primary_status = if (fam == "diffuse_region_like") "recommended_primary" else "candidate",
                usage_role = "field_first",
                keep_decision = "keep",
                rationale = "Best match to diffuse field semantics and geneSCOPE grid workflows."
            ),
            data.table(
                marker_name = rr$marker_name,
                marker_family = fam,
                candidate_id = "D4_boundary_segment_support",
                primary_status = if (fam == "boundary_associated") "recommended_primary" else "not_applicable",
                usage_role = "boundary_first",
                keep_decision = if (fam == "boundary_associated") "keep" else "not_applicable",
                rationale = "Uses local cell-boundary geometry instead of treating contour anchors as free interior points."
            ),
            data.table(
                marker_name = rr$marker_name,
                marker_family = fam,
                candidate_id = "D5_family_specific_hybrid",
                primary_status = "system_default_recommendation",
                usage_role = "hybrid_adapter",
                keep_decision = "keep",
                rationale = "Routes punctate to events, diffuse to D3, and boundary markers to D4."
            )
        )
    }), recursive = FALSE)
    rbindlist(rows, use.names = TRUE, fill = TRUE)
}

.pcdr_case_study_table <- function(marker_level, candidate_compare) {
    marker_dt <- as.data.table(marker_level)
    candidate_dt <- as.data.table(candidate_compare)
    if (!nrow(marker_dt)) {
        return(data.frame(
            marker_name = character(),
            marker_family = character(),
            validation_status = character(),
            priority_for_followup = character(),
            recommended_primary = character(),
            system_recommendation = character(),
            auxiliary_views = character(),
            main_rationale = character(),
            review_signal = character(),
            stringsAsFactors = FALSE
        ))
    }

    recommended_primary <- candidate_dt[
        primary_status == "recommended_primary",
        .(recommended_primary = paste(candidate_id, collapse = "; ")),
        by = .(marker_name)
    ]
    system_rec <- candidate_dt[
        primary_status == "system_default_recommendation",
        .(system_recommendation = paste(candidate_id, collapse = "; ")),
        by = .(marker_name)
    ]
    auxiliary <- candidate_dt[
        keep_decision == "keep" &
            !(primary_status %in% c("recommended_primary", "system_default_recommendation")) &
            candidate_id != "legacy_single_anchor" &
            candidate_id != "D4_boundary_segment_support" |
            (keep_decision == "keep" &
                candidate_id == "D4_boundary_segment_support" &
                marker_family == "diffuse_region_like"),
        .(auxiliary_views = paste(candidate_id, collapse = "; ")),
        by = .(marker_name)
    ]

    out <- merge(marker_dt, recommended_primary, by = "marker_name", all.x = TRUE)
    out <- merge(out, system_rec, by = "marker_name", all.x = TRUE)
    out <- merge(out, auxiliary, by = "marker_name", all.x = TRUE)
    out[, main_rationale := fifelse(
        marker_family == "boundary_associated",
        "Contour and membrane-like supports are better expressed as boundary-tied segments than as free interior points.",
        "Diffuse supports describe extended bright territory, so a single centroid is too easy to misread as the brightest-point coordinate."
    )]
    out[, review_signal := fifelse(
        is.finite(representation_mismatch_median),
        paste0(
            "mismatch_median=", .pc_fmt_number(representation_mismatch_median),
            "; offset_median=", .pc_fmt_number(review_offset_median),
            "; structure_median=", .pc_fmt_number(review_structure_median)
        ),
        "no_overlay_scores_available"
    )]
    out[, .(
        marker_name,
        marker_family,
        validation_status,
        priority_for_followup,
        recommended_primary = .pc_coalesce(recommended_primary, ""),
        system_recommendation = .pc_coalesce(system_recommendation, ""),
        auxiliary_views = .pc_coalesce(auxiliary_views, ""),
        main_rationale,
        review_signal
    )]
}

.pcdr_modeling_readiness_table <- function() {
    data.frame(
        marker_family_scope = c(
            "diffuse_region_like",
            "diffuse_region_like",
            "diffuse_region_like",
            "diffuse_region_like",
            "boundary_associated",
            "mixed_family"
        ),
        analysis_type = c(
            "density_gradient",
            "KDE",
            "graph_model",
            "PDE_field_model",
            "boundary_graph_or_contour_model",
            "geneSCOPE_multi_family_grid_analysis"
        ),
        preferred_representation = c(
            "D3_grid_field_overlap",
            "D2_multi_anchor_weighted_support",
            "D2_multi_anchor_weighted_support",
            "D3_grid_field_overlap",
            "D4_boundary_segment_support",
            "D5_family_specific_hybrid"
        ),
        secondary_representation = c(
            "D1_support_region_proxy",
            "D3_grid_field_overlap",
            "D1_support_region_proxy",
            "D1_support_region_proxy",
            "D3_grid_field_overlap",
            "family-specific primaries"
        ),
        genescope_consumption = c(
            "pcdr_field* or phenocycler_redesign_field_grid*",
            "phenocycler_redesign_d2_multi_anchors",
            "phenocycler_redesign_d2_multi_anchors",
            "pcdr_field* or phenocycler_redesign_field_grid*",
            "phenocycler_redesign_d4_segments and phenocycler_redesign_d4_segment_anchors",
            "pcdr_hybrid* or phenocycler_redesign_hybrid_grid*"
        ),
        rationale = c(
            "Regular grids preserve diffuse mass continuously and match gradient-based downstream methods.",
            "Weighted anchors retain coordinate-first flexibility for kernel reconstruction without collapsing each support to one point.",
            "Anchors plus local scale are the lightest bridge into neighborhood and graph constructions.",
            "Field-first layers are the most natural input to PDE and continuous-state models.",
            "Boundary-associated supports should remain tied to contour geometry rather than interior proxies.",
            "Hybrid routing is the safest cross-family entry when punctate, diffuse, and boundary channels must coexist in one geneSCOPE layer."
        ),
        stringsAsFactors = FALSE
    )
}

.pcdr_overlay_compare_manifest <- function(bundle_dir,
                                           legacy_bundle,
                                           redesign_bundle,
                                           selected_markers,
                                           output_dir) {
    .pc_require_pkg("jsonlite", ".pcdr_overlay_compare_manifest()")
    manifest_path <- file.path(bundle_dir, "metadata", "protein_signal_overlay_manifest.json")
    validation_manifest_path <- file.path(bundle_dir, "metadata", "full_dataset_validation_manifest.json")
    overlay_manifest <- if (file.exists(manifest_path)) {
        jsonlite::read_json(manifest_path, simplifyVector = TRUE)
    } else {
        NULL
    }
    validation_manifest <- if (file.exists(validation_manifest_path)) {
        jsonlite::read_json(validation_manifest_path, simplifyVector = TRUE)
    } else {
        NULL
    }
    legacy_supports <- as.data.table(legacy_bundle$inputs$signal_supports)
    redesign_regions <- as.data.table(redesign_bundle$representations$D1_support_regions$regions)
    redesign_anchors <- as.data.table(redesign_bundle$representations$D2_multi_anchors)
    boundary_segments <- as.data.table(redesign_bundle$representations$D4_boundary_segments$segments)
    boundary_vertices <- as.data.table(redesign_bundle$representations$D4_boundary_segments$vertices)
    region_vertices <- as.data.table(redesign_bundle$representations$D1_support_regions$vertices)
    plot_dir <- file.path(output_dir, "overlay_compare_plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    overlay_rows <- list()
    row_id <- 0L
    for (marker in selected_markers) {
        rois <- .pcdr_marker_roi_table(validation_manifest, marker)
        if (!nrow(rois)) next
        marker_folder <- .pcdr_overlay_folder(overlay_manifest, marker)
        for (i in seq_len(nrow(rois))) {
            roi <- rois[i, , drop = FALSE]
            xmin <- roi$x_px_level1
            xmax <- xmin + roi$width_px_level1
            ymin <- roi$y_px_level1
            ymax <- ymin + roi$height_px_level1
            old_roi <- legacy_supports[
                feature_name == marker &
                    x_location >= xmin & x_location <= xmax &
                    y_location >= ymin & y_location <= ymax
            ]
            region_roi <- redesign_regions[
                feature_name == marker &
                    bbox_xmax >= xmin & bbox_xmin <= xmax &
                    bbox_ymax >= ymin & bbox_ymin <= ymax
            ]
            anchor_roi <- redesign_anchors[
                feature_name == marker &
                    x_location >= xmin & x_location <= xmax &
                    y_location >= ymin & y_location <= ymax
            ]
            boundary_roi <- boundary_segments[
                feature_name == marker &
                    x_center >= xmin & x_center <= xmax &
                    y_center >= ymin & y_center <= ymax
            ]
            plot_path <- file.path(plot_dir, paste0(.pcdr_safe_name(marker), "_", roi$roi_id, "_representation_compare.png"))
            .pcdr_write_roi_compare_plot(
                plot_path = plot_path,
                marker = marker,
                roi = roi,
                old_roi = old_roi,
                region_roi = region_roi,
                region_vertices = region_vertices[feature_name == marker & region_id %in% region_roi$region_id],
                anchor_roi = anchor_roi,
                boundary_vertices = boundary_vertices[feature_name == marker & segment_id %in% boundary_roi$segment_id]
            )
            row_id <- row_id + 1L
            overlay_rows[[row_id]] <- data.table(
                marker_name = marker,
                roi_id = roi$roi_id,
                stratum = roi$stratum,
                x_px_level1 = roi$x_px_level1,
                y_px_level1 = roi$y_px_level1,
                width_px_level1 = roi$width_px_level1,
                height_px_level1 = roi$height_px_level1,
                legacy_support_count_roi = nrow(old_roi),
                redesign_region_count_roi = nrow(region_roi),
                redesign_anchor_count_roi = nrow(anchor_roi),
                redesign_boundary_segment_count_roi = nrow(boundary_roi),
                raw_channel_path = .pcdr_overlay_asset(marker_folder, roi$roi_id, "raw"),
                legacy_overlay_path = .pcdr_overlay_asset(marker_folder, roi$roi_id, "signal"),
                routeB_overlay_path = .pcdr_overlay_asset(marker_folder, roi$roi_id, "routeB"),
                redesign_compare_plot_path = normalizePath(plot_path, mustWork = TRUE)
            )
        }
    }
    manifest_dt <- if (length(overlay_rows)) {
        rbindlist(overlay_rows, use.names = TRUE, fill = TRUE)
    } else {
        data.table(
            marker_name = character(),
            roi_id = character(),
            stratum = character(),
            x_px_level1 = numeric(),
            y_px_level1 = numeric(),
            width_px_level1 = numeric(),
            height_px_level1 = numeric(),
            legacy_support_count_roi = integer(),
            redesign_region_count_roi = integer(),
            redesign_anchor_count_roi = integer(),
            redesign_boundary_segment_count_roi = integer(),
            raw_channel_path = character(),
            legacy_overlay_path = character(),
            routeB_overlay_path = character(),
            redesign_compare_plot_path = character()
        )
    }
    list(manifest = manifest_dt, plot_dir = plot_dir)
}

.pcdr_marker_roi_table <- function(validation_manifest, marker) {
    if (is.null(validation_manifest) || is.null(validation_manifest$signal_feature_review_rois)) {
        return(data.frame())
    }
    rois <- validation_manifest$signal_feature_review_rois[[marker]]
    if (is.null(rois)) return(data.frame())
    as.data.frame(rois)
}

.pcdr_overlay_folder <- function(overlay_manifest, marker) {
    if (is.null(overlay_manifest) || is.null(overlay_manifest$participating_markers_manual_overlay)) {
        return(NA_character_)
    }
    tab <- as.data.table(overlay_manifest$participating_markers_manual_overlay)
    hit <- tab[marker_name == marker, overlay_folder]
    if (!length(hit)) return(NA_character_)
    as.character(hit[[1]])
}

.pcdr_overlay_asset <- function(marker_folder, roi_id, kind = c("raw", "signal", "routeB")) {
    kind <- match.arg(kind)
    if (is.na(marker_folder) || !nzchar(marker_folder)) return(NA_character_)
    suffix <- switch(
        kind,
        raw = "_raw.png",
        signal = "_signal_overlay.png",
        routeB = "_signal_routeB_overlay.png"
    )
    path <- file.path(marker_folder, paste0(roi_id, suffix))
    if (!file.exists(path)) return(NA_character_)
    normalizePath(path, mustWork = TRUE)
}

.pcdr_write_roi_compare_plot <- function(plot_path,
                                         marker,
                                         roi,
                                         old_roi,
                                         region_roi,
                                         region_vertices,
                                         anchor_roi,
                                         boundary_vertices) {
    old_roi <- as.data.table(old_roi)
    region_vertices <- as.data.table(region_vertices)
    anchor_roi <- as.data.table(anchor_roi)
    boundary_vertices <- as.data.table(boundary_vertices)
    g <- ggplot2::ggplot() +
        ggplot2::coord_fixed(
            xlim = c(roi$x_px_level1, roi$x_px_level1 + roi$width_px_level1),
            ylim = c(roi$y_px_level1 + roi$height_px_level1, roi$y_px_level1)
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = paste0(marker, " / ", roi$roi_id),
            subtitle = "Legacy support centers vs redesigned region/anchor support"
        )
    if (nrow(region_vertices)) {
        g <- g + ggplot2::geom_polygon(
            data = region_vertices,
            mapping = ggplot2::aes(x = vertex_x, y = vertex_y, group = region_id),
            fill = "#5CA4A9",
            alpha = 0.20,
            colour = "#5CA4A9",
            linewidth = 0.2
        )
    }
    if (nrow(boundary_vertices)) {
        g <- g + ggplot2::geom_path(
            data = boundary_vertices,
            mapping = ggplot2::aes(x = vertex_x, y = vertex_y, group = segment_id),
            colour = "#D1495B",
            linewidth = 0.45
        )
    }
    if (nrow(anchor_roi)) {
        g <- g + ggplot2::geom_point(
            data = anchor_roi,
            mapping = ggplot2::aes(x = x_location, y = y_location),
            colour = "#1F7A8C",
            size = 0.45,
            alpha = 0.65
        )
    }
    if (nrow(old_roi)) {
        g <- g + ggplot2::geom_point(
            data = old_roi,
            mapping = ggplot2::aes(x = x_location, y = y_location),
            colour = "#EDAE49",
            size = 0.8,
            alpha = 0.85
        )
    }
    ggplot2::ggsave(
        filename = plot_path,
        plot = g,
        width = 5,
        height = 5,
        dpi = 160
    )
    invisible(plot_path)
}

.pcdr_marker_validation_status <- function(marker_family,
                                           review,
                                           legacy_rows,
                                           redesign_regions,
                                           redesign_anchors,
                                           redesign_segments) {
    mismatch_med <- .pcdr_num_median(review$representation_mismatch_risk)
    if (identical(marker_family, "boundary_associated")) {
        return(if (redesign_segments > 0) "PASS" else "CONDITIONAL_PASS")
    }
    if (identical(marker_family, "diffuse_region_like")) {
        if (legacy_rows > 0 && redesign_regions > 0 && redesign_anchors > legacy_rows) {
            return(if (is.finite(mismatch_med) && mismatch_med >= 4) "PASS" else "CONDITIONAL_PASS")
        }
        return("BLOCK")
    }
    "PASS"
}

.pcdr_scalar_or_blank <- function(x) {
    if (is.null(x) || !length(x)) return("")
    as.character(x[[1]])
}

.pcdr_num_median <- function(x) {
    if (is.null(x) || !length(x)) return(NA_real_)
    stats::median(as.numeric(x), na.rm = TRUE)
}

.pcdr_safe_name <- function(x) {
    gsub("[^A-Za-z0-9]+", "_", as.character(x))
}

#' Build the explicit diffuse-redesign bundle.
#' @description
#' Alias for `ingestPhenoCyclerBundleRedesign()` that preserves the
#' representation-first naming used by the专项 validation layer.
#' @inheritParams ingestPhenoCyclerBundleRedesign
#' @return A named list with class `genescope_phenocycler_redesign_bundle`.
#' @export
buildPhenoCyclerDiffuseRedesign <- function(
    bundle_dir,
    marker_policy = makePhenoCyclerMarkerPolicy(),
    redesign_config = makePhenoCyclerDiffuseRedesignConfig(),
    adapter_config = makePhenoCyclerRedesignAdapterConfig(),
    verbose = TRUE
) {
    .Deprecated("ingestPhenoCyclerBundleRedesign")
    ingestPhenoCyclerBundleRedesign(
        bundle_dir = bundle_dir,
        marker_policy = marker_policy,
        redesign_config = redesign_config,
        adapter_config = adapter_config,
        verbose = verbose
    )
}

#' Create a scope object from an organised PhenoCycler bundle using the
#' diffuse-redesign alias.
#' @description
#' Alias for `createSCOPE_phenocycler_redesign()` that preserves the
#' diffuse-redesign naming contract.
#' @param bundle_dir Path to the organised PhenoCycler bundle.
#' @param marker_policy Marker-family policy object.
#' @param redesign_config Diffuse redesign config.
#' @param adapter_config Redesign adapter config.
#' @param promote Which redesigned representation, if any, should be promoted
#'   into `@grid`.
#' @param verbose Emit progress messages when TRUE.
#' @return A `scope_object`.
#' @export
createSCOPE_phenocycler_diffuse_redesign <- function(
    bundle_dir,
    marker_policy = makePhenoCyclerMarkerPolicy(),
    redesign_config = makePhenoCyclerDiffuseRedesignConfig(),
    adapter_config = makePhenoCyclerRedesignAdapterConfig(),
    promote = adapter_config$promote,
    verbose = TRUE
) {
    .Deprecated("createSCOPE_phenocycler_redesign")
    createSCOPE_phenocycler_redesign(
        bundle_dir = bundle_dir,
        marker_policy = marker_policy,
        redesign_config = redesign_config,
        adapter_config = adapter_config,
        promote = promote,
        verbose = verbose
    )
}

#' Attach redesigned PhenoCycler representations through the
#' diffuse-redesign alias.
#' @description
#' Alias for `addPhenoCyclerRedesign()` that preserves the diffuse-redesign
#' naming contract while keeping the implementation in the current redesign
#' module.
#' @param scope_obj A `scope_object`.
#' @param bundle_dir Path to the organised PhenoCycler bundle.
#' @param redesign_bundle Optional already-ingested redesign bundle.
#' @param marker_policy Marker-family policy object.
#' @param redesign_config Diffuse redesign config.
#' @param adapter_config Redesign adapter config.
#' @param promote Which redesigned representation, if any, to promote.
#' @param overwrite_existing_modality Whether to overwrite an existing redesign
#'   modality sidecar.
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object`.
#' @export
addPhenoCyclerDiffuseRedesign <- function(
    scope_obj,
    bundle_dir = NULL,
    redesign_bundle = NULL,
    marker_policy = makePhenoCyclerMarkerPolicy(),
    redesign_config = makePhenoCyclerDiffuseRedesignConfig(),
    adapter_config = makePhenoCyclerRedesignAdapterConfig(),
    promote = adapter_config$promote,
    overwrite_existing_modality = FALSE,
    verbose = TRUE
) {
    .Deprecated("addPhenoCyclerRedesign")
    addPhenoCyclerRedesign(
        scope_obj = scope_obj,
        bundle_dir = bundle_dir,
        redesign_bundle = redesign_bundle,
        marker_policy = marker_policy,
        redesign_config = redesign_config,
        adapter_config = adapter_config,
        promote = promote,
        overwrite_existing_modality = overwrite_existing_modality,
        verbose = verbose
    )
}

#' Promote redesigned PhenoCycler representations into `@grid` through the
#' diffuse-redesign alias.
#' @description
#' Alias for `promotePhenoCyclerRedesignToGridLayer()` that preserves the
#' diffuse-redesign naming contract.
#' @param scope_obj A `scope_object`.
#' @param source Which redesigned representation to promote.
#' @param grid_lengths Optional subset of grid lengths to promote.
#' @param overwrite Whether to overwrite an existing grid layer.
#' @param verbose Emit progress messages when TRUE.
#' @return The modified `scope_object`.
#' @export
promotePhenoCyclerDiffuseRedesignToGridLayer <- function(
    scope_obj,
    source = c("field", "boundary", "hybrid", "punctate"),
    grid_lengths = NULL,
    overwrite = FALSE,
    verbose = TRUE
) {
    .Deprecated("promotePhenoCyclerRedesignToGridLayer")
    promotePhenoCyclerRedesignToGridLayer(
        scope_obj = scope_obj,
        source = source,
        grid_lengths = grid_lengths,
        overwrite = overwrite,
        verbose = verbose
    )
}

.pcdr_write_example_bundle <- function(out_dir) {
    .pc_require_pkg("arrow", ".pcdr_write_example_bundle()")
    .pc_require_pkg("jsonlite", ".pcdr_write_example_bundle()")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir, "molecules"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir, "cells"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir, "segmentation"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_dir, "metadata"), recursive = TRUE, showWarnings = FALSE)

    protein_events <- data.table(
        event_id = paste0("evt_", 1:8),
        feature_name = c("PDGFRA", "PDGFRA", "CD68", "CD68", "SATB2", "SATB2", "CLDN5", "CLDN5"),
        x_location = c(18, 24, 62, 68, 108, 114, 154, 160),
        y_location = c(18, 26, 22, 28, 20, 30, 18, 26),
        z_location = 0,
        source_type = c(
            "punctate_event", "connected_component_event",
            "punctate_event", "punctate_event",
            "connected_component_event", "connected_component_event",
            "connected_component_event", "connected_component_event"
        ),
        intensity_raw = c(120, 110, 90, 95, 150, 144, 80, 76),
        intensity_norm = c(1.2, 1.1, 0.9, 0.95, 1.5, 1.44, 0.8, 0.76),
        event_area_um2 = c(1, 2, 1, 1, 3, 3, 2, 2),
        qv_or_score = 30,
        cell_id = c("c1", "c1", "c2", "c2", "c3", "c3", "c4", "c4"),
        assignment_confidence = "exact",
        assignment_status = "interior",
        overlaps_nucleus = c(TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
        tile_id = "tile_001",
        fov_name = "fov_redesign",
        segmentation_version = "redesign_example_v1",
        event_call_method = "synthetic",
        resolution_level_used = 1,
        valid_pixel_fraction_local = 1,
        is_recovered_region = FALSE
    )
    signal_supports <- data.table(
        support_id = c(
            paste0("AQP4_sup_", 1:3),
            paste0("TUBB3_sup_", 1:3),
            paste0("pMYPT_sup_", 1:4),
            paste0("TMEM119_sup_", 1:6),
            paste0("F480_sup_", 1:6)
        ),
        feature_name = c(
            rep("AQP4", 3),
            rep("TUBB3", 3),
            rep("p-MYPT", 4),
            rep("TMEM119", 6),
            rep("F4/80", 6)
        ),
        x_location = c(
            24, 78, 126,
            32, 88, 138,
            20, 40, 60, 80,
            18, 22, 26, 150, 154, 158,
            16, 20, 24, 148, 152, 156
        ),
        y_location = c(
            64, 70, 68,
            112, 118, 114,
            150, 150, 150, 150,
            38, 42, 46, 42, 46, 50,
            92, 96, 100, 96, 100, 104
        ),
        source_type = c(
            rep("support_superpixel", 6),
            rep("support_grid_point", 4),
            rep("support_contour_anchor", 12)
        ),
        support_score = c(
            75, 62, 68,
            24, 28, 26,
            41, 42, 43, 40,
            52, 54, 53, 56, 55, 57,
            26, 28, 27, 30, 29, 31
        ),
        support_area_um2 = c(
            62000, 58000, 61000,
            49000, 51000, 47000,
            4218, 4218, 4218, 4218,
            rep(37.07055, 6),
            rep(37.07055, 6)
        ),
        cell_id = c(
            "c1", "c2", "c3",
            "c1", "c2", "c4",
            "c1", "c2", "c3", "c4",
            "c1", "c1", "c1", "c4", "c4", "c4",
            "c2", "c2", "c2", "c3", "c3", "c3"
        ),
        assignment_confidence = c(
            0.45, 0.30, 0.30,
            0.45, 0.45, 0.30,
            rep(0.45, 4),
            0.25, 0.20, 0.25, 0.25, 0.20, 0.25,
            0.25, 0.25, 0.20, 0.25, 0.25, 0.20
        ),
        assignment_status = c(
            "segmentation_sensitive_interior",
            "inside_cell_boundary_band",
            "unassigned",
            "segmentation_sensitive_interior",
            "inside_cell_boundary_band",
            "segmentation_sensitive_interior",
            rep("segmentation_sensitive_interior", 4),
            rep("boundary_associated", 6),
            rep("boundary_associated", 6)
        ),
        tile_id = "tile_001",
        fov_name = "fov_redesign",
        segmentation_version = "redesign_example_v1",
        support_call_method = c(
            rep("superpixel_centroid", 6),
            rep("grid_sampling", 4),
            rep("contour_projection", 12)
        ),
        resolution_level_used = 1,
        valid_pixel_fraction_local = 1,
        is_recovered_region = FALSE,
        assignment_semantics = c(
            rep("inside_cell", 10),
            rep("boundary_associated", 12)
        ),
        assignment_ambiguity = c(
            rep("medium", 6),
            rep("low", 4),
            rep("medium", 12)
        )
    )
    cells <- data.table(
        cell_id = paste0("c", 1:4),
        x_centroid = c(20, 70, 120, 160),
        y_centroid = c(20, 70, 120, 70),
        protein_event_counts = c(2, 2, 2, 2),
        approximate_event_count = c(2, 2, 2, 2),
        interior_support_count = c(3, 4, 4, 4),
        interior_support_score_sum = c(140, 132, 111, 138),
        boundary_support_count = c(3, 3, 3, 3),
        boundary_support_score_sum = c(159, 81, 84, 168),
        total_intensity = c(230, 182, 225, 258),
        cell_area = c(220, 210, 240, 230),
        nucleus_area = c(60, 58, 64, 61),
        dominant_marker = c("PDGFRA", "CD68", "SATB2", "CLDN5"),
        segmentation_version = "redesign_example_v1",
        assignment_policy_version = "assignment_policy_redesign_v1"
    )
    mk_poly <- function(id, cx, cy, half) {
        data.table(
            cell_id = id,
            vertex_x = c(cx - half, cx + half, cx + half, cx - half),
            vertex_y = c(cy - half, cy - half, cy + half, cy + half),
            vertex_index = 1:4,
            fov = "fov_redesign",
            label_id = id
        )
    }
    cell_boundaries <- rbindlist(list(
        mk_poly("c1", 20, 20, 12),
        mk_poly("c2", 70, 70, 12),
        mk_poly("c3", 120, 120, 12),
        mk_poly("c4", 160, 70, 12)
    ))
    nucleus_boundaries <- rbindlist(list(
        mk_poly("c1", 20, 20, 5),
        mk_poly("c2", 70, 70, 5),
        mk_poly("c3", 120, 120, 5),
        mk_poly("c4", 160, 70, 5)
    ))

    arrow::write_parquet(protein_events, file.path(out_dir, "molecules", "protein_events.parquet"))
    arrow::write_parquet(signal_supports, file.path(out_dir, "molecules", "signal_supports.parquet"))
    arrow::write_parquet(cells, file.path(out_dir, "cells", "cells.parquet"))
    arrow::write_parquet(cell_boundaries, file.path(out_dir, "segmentation", "cell_boundaries.parquet"))
    arrow::write_parquet(nucleus_boundaries, file.path(out_dir, "segmentation", "nucleus_boundaries.parquet"))
    jsonlite::write_json(
        list(
            dataset_id = "phenocycler_diffuse_redesign_bundle",
            version = "v1",
            purpose = "minimal validation bundle for family-aware diffuse redesign"
        ),
        file.path(out_dir, "metadata", "dataset_manifest.json"),
        auto_unbox = TRUE,
        pretty = TRUE
    )
    jsonlite::write_json(
        list(
            primary_analytical_layers = list(
                protein_events = "punctate or connected-component events",
                signal_supports = "diffuse or boundary support summaries"
            ),
            recommended_primary_views = list(
                punctate = "event_layer",
                diffuse = "grid_field",
                boundary = "boundary_weighted_anchor"
            )
        ),
        file.path(out_dir, "metadata", "data_semantics.json"),
        auto_unbox = TRUE,
        pretty = TRUE
    )
    invisible(out_dir)
}
