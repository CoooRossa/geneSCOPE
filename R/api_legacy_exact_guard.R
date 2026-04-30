#' Run the strict legacy exactness guard.
#' @description
#' Captures a stricter legacy contract than `runBackwardCompatGuard()` by
#' digesting actual grid payloads, sparse weights, normalised matrices, and a
#' wider set of frozen legacy files. This guard exists to prove exact numerical
#' legacy preservation during additive redesign work.
#' @param package_root Path to the geneSCOPE source tree. Defaults to the nearest
#'   parent directory containing `DESCRIPTION`.
#' @param baseline_path Path to the strict baseline JSON. Defaults to
#'   `inst/extdata/backward_compat/legacy_exact_contract_baseline.json` under
#'   `package_root`.
#' @param output_dir Directory where comparison artefacts should be written.
#' @param verbose Emit progress messages when TRUE.
#' @return A named list describing the current strict contract, the comparison
#'   result, and written artefact paths.
#' @export
runLegacyExactGuard <- function(
    package_root = NULL,
    baseline_path = NULL,
    output_dir = tempdir(),
    verbose = TRUE
) {
    .pc_require_pkg("jsonlite", "runLegacyExactGuard()")
    .pc_require_pkg("digest", "runLegacyExactGuard()")

    package_root <- .pc_guess_package_root(package_root)
    baseline_path <- .pc_coalesce(baseline_path, file.path(
        package_root,
        "inst", "extdata", "backward_compat", "legacy_exact_contract_baseline.json"
    ))
    if (!file.exists(baseline_path)) {
        stop("Legacy exact baseline not found: ", baseline_path)
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    baseline <- jsonlite::read_json(baseline_path, simplifyVector = FALSE)
    current <- .pcleg_capture_exact_contract(package_root = package_root, verbose = verbose)
    comparison <- .pcleg_compare_exact_contract(baseline, current)

    comparison_path <- file.path(output_dir, "legacy_exact_old_vs_new_comparison.tsv")
    current_path <- file.path(output_dir, "current_legacy_exact_contract.json")
    summary_path <- file.path(output_dir, "legacy_exact_summary.json")

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
            "[geneSCOPE::runLegacyExactGuard] pass=",
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

.pcleg_capture_exact_contract <- function(package_root, verbose = TRUE) {
    exports <- sort(getNamespaceExports("geneSCOPE"))
    export_formals <- stats::setNames(
        lapply(exports, function(nm) .pc_canon_formals(getExportedValue("geneSCOPE", nm))),
        exports
    )
    scope_slots <- methods::slotNames("scope_object")
    frozen_files <- c(
        "R/internal_leesL.R",
        "R/api_leesL.R",
        "R/internal_metrics_tests.R",
        "R/api_metrics_tests.R",
        "R/api_module_quality.R",
        "R/self_call_shims.R",
        "R/zzz_classes.R",
        "R/gridselect_external.R",
        "R/lifecycle_hooks.R",
        "R/api_data_construction.R",
        "R/internal_ingest.R",
        "R/api_normalization.R",
        "R/api_spatial_weights.R",
        "R/internal_spatial_weights.R",
        "R/api_plotting.R"
    )
    frozen_file_hashes <- stats::setNames(
        lapply(frozen_files, function(rel) {
            path <- file.path(package_root, rel)
            if (!file.exists(path)) return(NA_character_)
            unname(tools::md5sum(path))
        }),
        frozen_files
    )
    runtime_digests <- .pcleg_capture_runtime_digests(verbose = verbose)
    list(
        captured_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
        package_root = normalizePath(package_root, mustWork = TRUE),
        contract = list(
            exports = exports,
            export_formals = export_formals,
            scope_slots = scope_slots,
            frozen_file_hashes = frozen_file_hashes,
            runtime_digests = runtime_digests
        )
    )
}

.pcleg_capture_runtime_digests <- function(verbose = TRUE) {
    fixture <- file.path(tempdir(), "genescope_legacy_exact_guard_fixture")
    .pc_build_legacy_guard_fixture(fixture)
    ns <- asNamespace("geneSCOPE")
    create_xen <- get("createSCOPE_xenium", envir = ns)
    .pc_with_legacy_guard_native_fallbacks(ns, {
        obj1 <- create_xen(
            xenium_dir = fixture,
            grid_length = 20,
            seg_type = "both",
            num_workers = 1,
            verbose = FALSE,
            flip_y = FALSE,
            filter_molecules = FALSE,
            min_quality = NULL
        )
        obj2 <- createSCOPE(
            data_dir = fixture,
            prefer = "xenium",
            grid_length = 20,
            seg_type = "both",
            num_workers = 1,
            verbose = FALSE,
            flip_y = FALSE,
            filter_molecules = FALSE,
            min_quality = NULL
        )
        obj1n <- normalizeMoleculesInGrid(obj1, grid_name = "grid20", verbose = FALSE)
        obj1w <- computeWeights(obj1n, grid_name = "grid20", verbose = FALSE, ncores = 1)
        list(
            createSCOPE_xenium = digest::digest(.pcleg_scope_snapshot(obj1), algo = "sha256"),
            createSCOPE_dispatch = digest::digest(.pcleg_scope_snapshot(obj2), algo = "sha256"),
            normalizeMoleculesInGrid = digest::digest(.pcleg_scope_snapshot(obj1n), algo = "sha256"),
            computeWeights = digest::digest(.pcleg_scope_snapshot(obj1w), algo = "sha256")
        )
    })
}

.pcleg_scope_snapshot <- function(obj) {
    list(
        coord = .pcleg_digestable(obj@coord),
        grid = .pcleg_digestable(obj@grid),
        meta_data = .pcleg_digestable(obj@meta.data),
        cells = .pcleg_digestable(obj@cells),
        density = .pcleg_digestable(obj@density),
        stats = .pcleg_digestable(obj@stats)
    )
}

.pcleg_digestable <- function(x) {
    if (is.null(x)) return(NULL)
    if (methods::is(x, "Matrix")) {
        return(.pcleg_canon_matrix(x))
    }
    if (is.data.frame(x)) {
        return(.pcleg_canon_data_frame(x))
    }
    if (is.matrix(x)) {
        return(list(
            kind = "matrix",
            dim = dim(x),
            dimnames = dimnames(x),
            values = as.vector(unclass(x))
        ))
    }
    if (is.list(x)) {
        nms <- names(x)
        if (is.null(nms)) nms <- as.character(seq_along(x))
        return(stats::setNames(lapply(x, .pcleg_digestable), nms))
    }
    if (is.factor(x)) {
        return(as.character(x))
    }
    if (inherits(x, "Date") || inherits(x, "POSIXt")) {
        return(as.character(x))
    }
    if (is.atomic(x)) {
        return(unclass(x))
    }
    if (methods::is(x, "S4")) {
        slots <- methods::slotNames(x)
        return(stats::setNames(lapply(slots, function(nm) .pcleg_digestable(methods::slot(x, nm))), slots))
    }
    as.character(x)
}

.pcleg_canon_data_frame <- function(df) {
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    df[] <- lapply(df, function(col) {
        if (is.factor(col)) {
            as.character(col)
        } else if (inherits(col, "POSIXt") || inherits(col, "Date")) {
            as.character(col)
        } else {
            col
        }
    })
    if (nrow(df)) {
        order_cols <- intersect(
            c(
                "grid_id", "gene", "feature_name", "cell_id", "event_id",
                "support_id", "anchor_id", "vertex_index", "gx", "gy",
                "x_location", "y_location", "x_centroid", "y_centroid"
            ),
            names(df)
        )
        if (!length(order_cols)) {
            order_cols <- names(df)
        }
        ord <- do.call(order, c(lapply(df[order_cols], .pcleg_orderable), list(na.last = TRUE)))
        df <- df[ord, , drop = FALSE]
        rownames(df) <- NULL
    }
    list(
        kind = "data.frame",
        names = names(df),
        rows = nrow(df),
        columns = lapply(df, .pcleg_digestable)
    )
}

.pcleg_orderable <- function(x) {
    if (is.list(x)) return(vapply(x, function(v) paste(capture.output(str(v)), collapse = " "), character(1)))
    if (is.factor(x)) return(as.character(x))
    if (inherits(x, "POSIXt") || inherits(x, "Date")) return(as.character(x))
    x
}

.pcleg_canon_matrix <- function(x) {
    trip <- Matrix::summary(x)
    ord <- if (nrow(trip)) order(trip$i, trip$j, na.last = TRUE) else integer()
    if (length(ord)) trip <- trip[ord, , drop = FALSE]
    list(
        kind = "Matrix",
        class = class(x),
        dim = dim(x),
        dimnames = dimnames(x),
        triplet = list(
            i = if (nrow(trip)) trip$i else integer(),
            j = if (nrow(trip)) trip$j else integer(),
            x = if (nrow(trip)) trip$x else numeric()
        )
    )
}

.pcleg_compare_exact_contract <- function(baseline, current) {
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

    baseline_exports <- unlist(b$exports, use.names = FALSE)
    current_exports <- unlist(c$exports, use.names = FALSE)
    missing_exports <- setdiff(baseline_exports, current_exports)
    if (length(missing_exports)) {
        for (nm in missing_exports) {
            append_row("exports", nm, TRUE, FALSE, "FAIL", "legacy export missing from strict current namespace")
        }
    } else {
        append_row("exports", "baseline_subset_present", length(baseline_exports), length(current_exports), "PASS", "all strict legacy exports remain available")
    }

    for (nm in names(b$export_formals)) {
        status <- if (identical(b$export_formals[[nm]], c$export_formals[[nm]])) "PASS" else "FAIL"
        append_row(
            "export_formals",
            nm,
            b$export_formals[[nm]],
            c$export_formals[[nm]],
            status,
            if (status == "FAIL") "legacy signature drift under strict guard" else ""
        )
    }

    append_row(
        "scope_slots",
        "scope_object",
        paste(unlist(b$scope_slots, use.names = FALSE), collapse = ","),
        paste(unlist(c$scope_slots, use.names = FALSE), collapse = ","),
        if (identical(unlist(b$scope_slots, use.names = FALSE), unlist(c$scope_slots, use.names = FALSE))) "PASS" else "FAIL",
        if (!identical(unlist(b$scope_slots, use.names = FALSE), unlist(c$scope_slots, use.names = FALSE))) "scope_object slot contract changed" else ""
    )

    for (rel in names(b$frozen_file_hashes)) {
        status <- if (identical(b$frozen_file_hashes[[rel]], c$frozen_file_hashes[[rel]])) "PASS" else "FAIL"
        append_row(
            "frozen_file_hashes",
            rel,
            b$frozen_file_hashes[[rel]],
            c$frozen_file_hashes[[rel]],
            status,
            if (status == "FAIL") "strict frozen legacy file changed" else ""
        )
    }

    for (nm in names(b$runtime_digests)) {
        status <- if (identical(b$runtime_digests[[nm]], c$runtime_digests[[nm]])) "PASS" else "FAIL"
        append_row(
            "runtime_digests",
            nm,
            b$runtime_digests[[nm]],
            c$runtime_digests[[nm]],
            status,
            if (status == "FAIL") "strict legacy runtime digest drift" else ""
        )
    }

    rows_dt <- rbindlist(rows, use.names = TRUE, fill = TRUE)
    list(
        pass = !any(rows_dt$status == "FAIL"),
        summary = rows_dt[, .N, by = status][order(status)],
        rows = rows_dt
    )
}
