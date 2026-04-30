#' Run the PhenoCycler v1 additive-contract guard.
#' @description
#' Captures the current contract of the previously introduced additive
#' PhenoCycler v1 module and compares it with a pinned baseline. This guard is
#' separate from `runBackwardCompatGuard()` so the diffuse-redesign work can
#' prove it did not change current PhenoCycler v1 defaults.
#' @param package_root Path to the geneSCOPE source tree. Defaults to the
#'   nearest parent directory containing `DESCRIPTION`.
#' @param baseline_path Path to the PhenoCycler v1 baseline JSON. Defaults to
#'   `inst/extdata/backward_compat/phenocycler_contract_baseline.json` under
#'   `package_root`.
#' @param output_dir Directory where comparison artefacts should be written.
#' @param verbose Emit progress messages when TRUE.
#' @return A named list describing the current contract, the comparison result,
#'   and written artefact paths.
#' @export
runPhenoCyclerContractGuard <- function(
    package_root = NULL,
    baseline_path = NULL,
    output_dir = tempdir(),
    verbose = TRUE
) {
    .pc_require_pkg("jsonlite", "runPhenoCyclerContractGuard()")
    .pc_require_pkg("digest", "runPhenoCyclerContractGuard()")
    package_root <- .pc_guess_package_root(package_root)
    baseline_path <- .pc_coalesce(baseline_path, file.path(
        package_root,
        "inst", "extdata", "backward_compat", "phenocycler_contract_baseline.json"
    ))
    if (!file.exists(baseline_path)) {
        stop("PhenoCycler contract baseline not found: ", baseline_path)
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    baseline <- jsonlite::read_json(baseline_path, simplifyVector = FALSE)
    current <- .pcp_capture_contract(package_root = package_root, verbose = verbose)
    comparison <- .pcp_compare_contract(baseline, current)

    comparison_path <- file.path(output_dir, "phenocycler_old_vs_new_exact_comparison.tsv")
    current_path <- file.path(output_dir, "current_phenocycler_contract.json")
    summary_path <- file.path(output_dir, "phenocycler_contract_summary.json")

    data.table::fwrite(comparison$rows, comparison_path, sep = "\t")
    jsonlite::write_json(current, current_path, auto_unbox = TRUE, pretty = TRUE)
    jsonlite::write_json(
        list(
            pass = comparison$pass,
            summary = comparison$summary,
            baseline_path = normalizePath(baseline_path, mustWork = TRUE),
            current_path = normalizePath(current_path, mustWork = TRUE),
            comparison_path = normalizePath(comparison_path, mustWork = TRUE)
        ),
        summary_path,
        auto_unbox = TRUE,
        pretty = TRUE
    )

    if (verbose) {
        message(
            "[geneSCOPE::runPhenoCyclerContractGuard] pass=",
            comparison$pass,
            " rows=", nrow(comparison$rows),
            " output_dir=", normalizePath(output_dir, mustWork = TRUE)
        )
    }

    list(
        pass = comparison$pass,
        summary = comparison$summary,
        rows = comparison$rows,
        baseline = baseline,
        current = current,
        artefacts = list(
            comparison_tsv = comparison_path,
            current_json = current_path,
            summary_json = summary_path
        )
    )
}

.pcp_capture_contract <- function(package_root, verbose = TRUE) {
    exports <- c(
        "makePhenoCyclerMarkerPolicy",
        "makePhenoCyclerDiffuseConfig",
        "makePhenoCyclerAdapterConfig",
        "ingestPhenoCyclerBundle",
        "createSCOPE_phenocycler",
        "addPhenoCycler",
        "promotePhenoCyclerToGridLayer"
    )
    export_formals <- stats::setNames(
        lapply(exports, function(nm) .pc_canon_formals(getExportedValue("geneSCOPE", nm))),
        exports
    )
    real_bundle_dir <- .pcp_guess_real_bundle_dir(package_root)
    list(
        captured_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
        package_root = normalizePath(package_root, mustWork = TRUE),
        contract = list(
            exports = exports,
            export_formals = export_formals,
            config_defaults = .pcp_digest_config_defaults(),
            minimal_fixture = .pcp_capture_minimal_fixture(verbose = verbose),
            real_bundle_summary = .pcp_capture_real_bundle_summary(real_bundle_dir, verbose = verbose)
        )
    )
}

.pcp_guess_real_bundle_dir <- function(package_root) {
    candidates <- c(
        file.path(dirname(package_root), "251021_Yoneyama-kun_Pcdh15", "genescope_coordinate_first"),
        file.path(dirname(dirname(package_root)), "251021_Yoneyama-kun_Pcdh15", "genescope_coordinate_first")
    )
    hits <- candidates[file.exists(file.path(candidates, "molecules", "protein_events.parquet"))]
    if (!length(hits)) return(NULL)
    normalizePath(hits[[1]], mustWork = TRUE)
}

.pcp_digest_config_defaults <- function() {
    list(
        marker_policy = digest::digest(makePhenoCyclerMarkerPolicy(), algo = "sha256"),
        diffuse_config = digest::digest(makePhenoCyclerDiffuseConfig(), algo = "sha256"),
        adapter_config = digest::digest(makePhenoCyclerAdapterConfig(), algo = "sha256")
    )
}

.pcp_capture_minimal_fixture <- function(verbose = TRUE) {
    ns <- asNamespace("geneSCOPE")
    fixture <- file.path(tempdir(), "genescope_phenocycler_contract_fixture")
    get(".pc_write_example_bundle", envir = ns)(fixture)
    bundle <- geneSCOPE::ingestPhenoCyclerBundle(fixture, verbose = FALSE)
    obj_default <- geneSCOPE::createSCOPE_phenocycler(fixture, verbose = FALSE)
    obj_field <- geneSCOPE::createSCOPE_phenocycler(fixture, promote = "field", verbose = FALSE)
    obj_promoted <- geneSCOPE::promotePhenoCyclerToGridLayer(obj_default, source = "field", overwrite = FALSE, verbose = FALSE)

    if (verbose) {
        message(
            "[geneSCOPE::.pcp_capture_minimal_fixture] fixture=",
            normalizePath(fixture, mustWork = TRUE)
        )
    }

    list(
        manifest = list(
            event_rows = bundle$manifest$event_rows,
            support_rows = bundle$manifest$support_rows,
            feature_count = bundle$manifest$feature_count,
            grid_lengths = unname(as.numeric(bundle$manifest$grid_lengths)),
            bounds = bundle$manifest$bounds
        ),
        digests = list(
            protein_events = .pcp_digest_table(bundle$inputs$protein_events, c("event_id")),
            signal_supports = .pcp_digest_table(bundle$inputs$signal_supports, c("support_id")),
            feature_metadata = .pcp_digest_table(bundle$metadata$feature_metadata, c("feature_name")),
            weighted_anchors = .pcp_digest_table(bundle$representations$weighted_anchors, c("anchor_id")),
            support_regions = .pcp_digest_table(bundle$representations$support_regions, c("region_id")),
            grid20_punctate = .pcp_digest_grid_layer(bundle$representations$grid_layers$grid20$punctate),
            grid20_field = .pcp_digest_grid_layer(bundle$representations$grid_layers$grid20$field),
            grid20_boundary = .pcp_digest_grid_layer(bundle$representations$grid_layers$grid20$boundary),
            grid20_dual = .pcp_digest_grid_layer(bundle$representations$grid_layers$grid20$dual),
            grid40_punctate = .pcp_digest_grid_layer(bundle$representations$grid_layers$grid40$punctate),
            grid40_field = .pcp_digest_grid_layer(bundle$representations$grid_layers$grid40$field),
            grid40_boundary = .pcp_digest_grid_layer(bundle$representations$grid_layers$grid40$boundary),
            grid40_dual = .pcp_digest_grid_layer(bundle$representations$grid_layers$grid40$dual)
        ),
        scope_default = .pcp_digest_scope_state(obj_default),
        scope_field = .pcp_digest_scope_state(obj_field),
        scope_promoted = .pcp_digest_scope_state(obj_promoted)
    )
}

.pcp_capture_real_bundle_summary <- function(real_bundle_dir, verbose = TRUE) {
    if (is.null(real_bundle_dir)) {
        return(list(found = FALSE))
    }
    bundle <- geneSCOPE::ingestPhenoCyclerBundle(real_bundle_dir, verbose = FALSE)
    supports <- as.data.table(bundle$inputs$signal_supports)
    src_files <- c(
        protein_events = file.path(real_bundle_dir, "molecules", "protein_events.parquet"),
        signal_supports = file.path(real_bundle_dir, "molecules", "signal_supports.parquet")
    )
    hashes <- stats::setNames(
        lapply(src_files, function(path) unname(tools::md5sum(path))),
        names(src_files)
    )
    if (verbose) {
        message(
            "[geneSCOPE::.pcp_capture_real_bundle_summary] bundle=",
            real_bundle_dir
        )
    }
    list(
        found = TRUE,
        bundle_dir = real_bundle_dir,
        source_hashes = hashes,
        event_rows = bundle$manifest$event_rows,
        support_rows = bundle$manifest$support_rows,
        feature_rows = nrow(bundle$metadata$feature_metadata),
        marker_family_counts = supports[, .N, by = marker_family][order(marker_family)],
        weighted_anchor_rows = nrow(bundle$representations$weighted_anchors),
        support_region_rows = nrow(bundle$representations$support_regions),
        grid_lengths = unname(as.numeric(bundle$manifest$grid_lengths)),
        feature_metadata_digest = .pcp_digest_table(bundle$metadata$feature_metadata, c("feature_name"))
    )
}

.pcp_digest_scope_state <- function(obj) {
    list(
        coord_names = sort(names(obj@coord)),
        grid_names = sort(names(obj@grid)),
        density_names = sort(names(obj@density)),
        meta_data = .pcp_digest_table(obj@meta.data, c("feature_name")),
        coord_digests = stats::setNames(
            lapply(names(obj@coord), function(nm) .pcp_digest_table(obj@coord[[nm]], c("feature_name", "cell_id", "x_location", "y_location"))),
            names(obj@coord)
        ),
        grid_digests = stats::setNames(
            lapply(names(obj@grid), function(nm) .pcp_digest_grid_layer(obj@grid[[nm]])),
            names(obj@grid)
        ),
        density_digests = stats::setNames(
            lapply(names(obj@density), function(nm) .pcp_digest_table(obj@density[[nm]], c("grid_id"))),
            names(obj@density)
        )
    )
}

.pcp_compare_contract <- function(baseline, current) {
    b <- baseline$contract
    c <- current$contract
    rows <- list()

    append_row <- function(section, item, baseline_value, current_value, status, note = "") {
        rows[[length(rows) + 1L]] <<- data.table(
            section = section,
            item = item,
            baseline = .pc_scalar_string(baseline_value),
            current = .pc_scalar_string(current_value),
            status = status,
            note = note
        )
    }

    for (nm in b$exports) {
        append_row(
            "exports",
            nm,
            TRUE,
            nm %in% unlist(c$exports, use.names = FALSE),
            if (nm %in% unlist(c$exports, use.names = FALSE)) "PASS" else "FAIL",
            if (!(nm %in% unlist(c$exports, use.names = FALSE))) "PhenoCycler v1 export missing" else ""
        )
    }

    for (nm in names(b$export_formals)) {
        status <- if (identical(b$export_formals[[nm]], c$export_formals[[nm]])) "PASS" else "FAIL"
        append_row(
            "export_formals",
            nm,
            b$export_formals[[nm]],
            c$export_formals[[nm]],
            status,
            if (status == "FAIL") "PhenoCycler v1 export signature drift" else ""
        )
    }

    for (nm in names(b$config_defaults)) {
        status <- if (identical(b$config_defaults[[nm]], c$config_defaults[[nm]])) "PASS" else "FAIL"
        append_row(
            "config_defaults",
            nm,
            b$config_defaults[[nm]],
            c$config_defaults[[nm]],
            status,
            if (status == "FAIL") "PhenoCycler v1 default config drift" else ""
        )
    }

    .pcp_append_nested_rows(rows = rows, section = "minimal_manifest", baseline = b$minimal_fixture$manifest, current = c$minimal_fixture$manifest)
    .pcp_append_nested_rows(rows = rows, section = "minimal_digests", baseline = b$minimal_fixture$digests, current = c$minimal_fixture$digests)
    .pcp_append_nested_rows(rows = rows, section = "scope_default", baseline = b$minimal_fixture$scope_default, current = c$minimal_fixture$scope_default)
    .pcp_append_nested_rows(rows = rows, section = "scope_field", baseline = b$minimal_fixture$scope_field, current = c$minimal_fixture$scope_field)
    .pcp_append_nested_rows(rows = rows, section = "scope_promoted", baseline = b$minimal_fixture$scope_promoted, current = c$minimal_fixture$scope_promoted)
    .pcp_append_nested_rows(rows = rows, section = "real_bundle_summary", baseline = b$real_bundle_summary, current = c$real_bundle_summary)

    rows_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
    fail_rows <- rows_dt[status != "PASS"]
    list(
        pass = !nrow(fail_rows),
        summary = rows_dt[, .N, by = status][order(status)],
        rows = rows_dt
    )
}

.pcp_append_nested_rows <- function(rows, section, baseline, current, prefix = NULL) {
    append_row <- function(item, baseline_value, current_value, status, note = "") {
        rows[[length(rows) + 1L]] <<- data.table(
            section = section,
            item = item,
            baseline = .pc_scalar_string(baseline_value),
            current = .pc_scalar_string(current_value),
            status = status,
            note = note
        )
    }

    b_names <- names(baseline)
    c_names <- names(current)
    all_names <- union(b_names, c_names)
    for (nm in all_names) {
        item <- if (is.null(prefix)) nm else paste(prefix, nm, sep = ".")
        b_val <- baseline[[nm]]
        c_val <- current[[nm]]
        if (is.list(b_val) && is.list(c_val) && !is.data.frame(b_val) && !is.data.frame(c_val)) {
            rows <- .pcp_append_nested_rows(rows, section, b_val, c_val, prefix = item)
        } else {
            status <- if (identical(b_val, c_val)) "PASS" else "FAIL"
            note <- if (status == "FAIL") paste(section, "drift") else ""
            append_row(item, b_val, c_val, status, note)
        }
    }
    rows
}

.pcp_digest_table <- function(x, order_cols = character()) {
    if (is.null(x)) return(NA_character_)
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    if (!nrow(df)) {
        return(digest::digest(list(names = names(df), rows = 0L), algo = "sha256"))
    }
    df[] <- lapply(df, function(col) {
        if (is.factor(col)) as.character(col) else col
    })
    if (length(order_cols)) {
        order_cols <- intersect(order_cols, names(df))
    }
    if (!length(order_cols)) {
        order_cols <- names(df)
    }
    ord <- do.call(order, c(df[order_cols], na.last = TRUE))
    df <- df[ord, , drop = FALSE]
    rownames(df) <- NULL
    digest::digest(df, algo = "sha256")
}

.pcp_digest_grid_layer <- function(layer) {
    if (is.null(layer)) return(NA_character_)
    payload <- list(
        grid_info = .pcp_digest_table(layer$grid_info, c("grid_id")),
        counts = .pcp_digest_table(layer$counts, c("grid_id", "gene")),
        grid_length = as.numeric(layer$grid_length),
        xbins_eff = layer$xbins_eff,
        ybins_eff = layer$ybins_eff,
        representation_type = layer$representation_type,
        source_platform = layer$source_platform
    )
    digest::digest(payload, algo = "sha256")
}
