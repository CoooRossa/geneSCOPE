#' Run the strict backward-compatibility guard.
#' @description
#' Captures the current legacy contract of the geneSCOPE source tree and compares
#' it with a pinned pre-change baseline. The guard is intentionally explicit
#' about the legacy baseline boundary: runtime digests for legacy constructors
#' and `computeWeights()` are captured under a controlled R fallback harness,
#' while historical Lee's L runtime digests remain omitted from this guard and
#' are covered by the dedicated native/reference backend tests.
#' @param package_root Path to the geneSCOPE source tree. Defaults to the nearest
#'   parent directory containing `DESCRIPTION`.
#' @param baseline_path Path to the pre-change baseline JSON. Defaults to
#'   `inst/extdata/backward_compat/legacy_contract_baseline.json` under
#'   `package_root`.
#' @param output_dir Directory where comparison artefacts should be written.
#' @param verbose Emit progress messages when TRUE.
#' @return A named list describing the current contract, the comparison result,
#'   and written artefact paths.
#' @export
runBackwardCompatGuard <- function(
    package_root = NULL,
    baseline_path = NULL,
    output_dir = tempdir(),
    verbose = TRUE
) {
    .pc_require_pkg("jsonlite", "runBackwardCompatGuard()")
    .pc_require_pkg("arrow", "runBackwardCompatGuard()")
    .pc_require_pkg("digest", "runBackwardCompatGuard()")

    package_root <- .pc_guess_package_root(package_root)
    baseline_path <- .pc_coalesce(baseline_path, file.path(
        package_root,
        "inst", "extdata", "backward_compat", "legacy_contract_baseline.json"
    ))
    if (!file.exists(baseline_path)) {
        stop("Backward-compat baseline not found: ", baseline_path)
    }
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    baseline <- jsonlite::read_json(baseline_path, simplifyVector = FALSE)
    current <- .pc_capture_legacy_contract(package_root = package_root, verbose = verbose)
    comparison <- .pc_compare_legacy_contract(baseline, current)

    comparison_path <- file.path(output_dir, "old_vs_new_exact_comparison.tsv")
    current_path <- file.path(output_dir, "current_legacy_contract.json")
    summary_path <- file.path(output_dir, "backward_compat_summary.json")

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
            "[geneSCOPE::runBackwardCompatGuard] pass=",
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

.pc_capture_legacy_contract <- function(package_root, verbose = TRUE) {
    ns <- asNamespace("geneSCOPE")
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
        "R/lifecycle_hooks.R"
    )
    frozen_file_hashes <- stats::setNames(
        lapply(frozen_files, function(rel) {
            path <- file.path(package_root, rel)
            if (!file.exists(path)) return(NA_character_)
            unname(tools::md5sum(path))
        }),
        frozen_files
    )
    runtime_digests <- .pc_capture_legacy_runtime_digests(verbose = verbose)
    list(
        captured_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
        package_root = normalizePath(package_root, mustWork = TRUE),
        contract = list(
            exports = exports,
            export_formals = export_formals,
            scope_slots = scope_slots,
            frozen_file_hashes = frozen_file_hashes,
            runtime_digests = runtime_digests,
            omitted_runtime = list(
                computeL = "Runtime digest intentionally omitted from the legacy compatibility guard because the historical baseline avoided invoking Lee's L native runtime. Native/reference Lee's L closure is verified by dedicated tests outside this guard."
            )
        )
    )
}

.pc_capture_legacy_runtime_digests <- function(verbose = TRUE) {
    fixture <- file.path(tempdir(), "genescope_legacy_guard_runtime_fixture")
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
            createSCOPE_xenium = digest::digest(.pc_legacy_object_snapshot(obj1), algo = "sha256"),
            createSCOPE_dispatch = digest::digest(.pc_legacy_object_snapshot(obj2), algo = "sha256"),
            normalizeMoleculesInGrid = digest::digest(.pc_legacy_object_snapshot(obj1n), algo = "sha256"),
            computeWeights = digest::digest(.pc_legacy_object_snapshot(obj1w), algo = "sha256")
        )
    })
}

.pc_compare_legacy_contract <- function(baseline, current) {
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
        for (nm in missing_exports) append_row("exports", nm, TRUE, FALSE, "FAIL", "legacy export missing from current namespace")
    } else {
        append_row("exports", "baseline_subset_present", length(baseline_exports), length(current_exports), "PASS", "all legacy exports remain available; additive exports are allowed")
    }

    baseline_formals <- b$export_formals
    current_formals <- c$export_formals
    for (nm in names(baseline_formals)) {
        status <- if (identical(baseline_formals[[nm]], current_formals[[nm]])) "PASS" else "FAIL"
        append_row("export_formals", nm, baseline_formals[[nm]], current_formals[[nm]], status, if (status == "FAIL") "legacy function signature drift" else "")
    }

    baseline_slots <- unlist(b$scope_slots, use.names = FALSE)
    current_slots <- unlist(c$scope_slots, use.names = FALSE)
    append_row(
        "scope_slots",
        "scope_object",
        paste(baseline_slots, collapse = ","),
        paste(current_slots, collapse = ","),
        if (identical(baseline_slots, current_slots)) "PASS" else "FAIL",
        if (!identical(baseline_slots, current_slots)) "scope_object slot contract changed" else ""
    )

    baseline_hashes <- b$frozen_file_hashes
    current_hashes <- c$frozen_file_hashes
    for (rel in names(baseline_hashes)) {
        status <- if (identical(baseline_hashes[[rel]], current_hashes[[rel]])) "PASS" else "FAIL"
        append_row("frozen_file_hashes", rel, baseline_hashes[[rel]], current_hashes[[rel]], status, if (status == "FAIL") "frozen legacy file changed" else "")
    }

    baseline_runtime <- b$runtime_digests
    current_runtime <- c$runtime_digests
    for (nm in names(baseline_runtime)) {
        status <- if (identical(baseline_runtime[[nm]], current_runtime[[nm]])) "PASS" else "FAIL"
        append_row("runtime_digests", nm, baseline_runtime[[nm]], current_runtime[[nm]], status, if (status == "FAIL") "legacy runtime digest drift under controlled fallback harness" else "")
    }

    comparison_rows <- rbindlist(rows, use.names = TRUE, fill = TRUE)
    pass <- !any(comparison_rows$status == "FAIL")
    summary <- comparison_rows[, .N, by = status][order(status)]
    list(pass = pass, rows = comparison_rows, summary = summary)
}

.pc_with_legacy_guard_native_fallbacks <- function(ns, expr) {
    expr_sub <- substitute(expr)
    patch <- list(
        ".native_openmp_info" = function() NULL,
        ".native_openmp_set_threads" = function(n_threads = 1L) {
            list(
                compiled_with_openmp = NA,
                requested_threads = as.integer(n_threads),
                omp_max_threads = as.integer(n_threads),
                omp_num_procs = as.integer(n_threads)
            )
        },
        ".grid_nb_safe" = function(nrow,
                                    ncol,
                                    topology = c("queen", "rook"),
                                    queen = TRUE,
                                    use_parallel = FALSE,
                                    verbose = FALSE,
                                    parent = "computeWeights",
                                    step = "S03") {
            list(nb = .grid_nb_r(nrow, ncol, queen = queen), backend = "R")
        },
        ".listw_b_omp" = function(nb) {
            n <- length(nb)
            ii <- rep.int(seq_len(n), lengths(nb))
            jj <- if (length(nb)) unlist(nb, use.names = FALSE) else integer()
            Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, n))
        }
    )
    orig <- lapply(names(patch), function(nm) get(nm, envir = ns, inherits = FALSE))
    names(orig) <- names(patch)
    for (nm in names(patch)) unlockBinding(nm, ns)
    on.exit({
        for (nm in names(orig)) {
            if (bindingIsLocked(nm, ns)) unlockBinding(nm, ns)
            assign(nm, orig[[nm]], envir = ns)
            lockBinding(nm, ns)
        }
    }, add = TRUE)
    for (nm in names(patch)) assign(nm, patch[[nm]], envir = ns)
    for (nm in rev(names(patch))) lockBinding(nm, ns)
    eval(expr_sub, envir = parent.frame())
}

.pc_build_legacy_guard_fixture <- function(fixture_dir) {
    .pc_require_pkg("arrow", ".pc_build_legacy_guard_fixture()")
    dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)
    transcripts <- data.table(
        x_location = c(10, 12, 11, 30, 32, 31, 50, 52, 51, 70, 72, 71, 15, 35, 55, 75),
        y_location = c(10, 12, 11, 10, 12, 11, 10, 12, 11, 10, 12, 11, 30, 30, 30, 30),
        feature_name = c("G1", "G1", "G2", "G1", "G2", "G2", "G1", "G3", "G3", "G2", "G3", "G3", "G1", "G2", "G1", "G3"),
        qv = 30,
        nucleus_distance = 5
    )
    cells <- data.table(
        cell_id = c("c1", "c2", "c3", "c4"),
        x_centroid = c(11, 31, 51, 71),
        y_centroid = c(11, 11, 11, 11)
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
        mk_poly("c1", 11, 11, 6),
        mk_poly("c2", 31, 11, 6),
        mk_poly("c3", 51, 11, 6),
        mk_poly("c4", 71, 11, 6)
    ))
    nucleus_boundaries <- rbindlist(list(
        mk_poly("c1", 11, 11, 3),
        mk_poly("c2", 31, 11, 3),
        mk_poly("c3", 51, 11, 3),
        mk_poly("c4", 71, 11, 3)
    ))
    arrow::write_parquet(transcripts, file.path(fixture_dir, "transcripts.parquet"))
    arrow::write_parquet(cells, file.path(fixture_dir, "cells.parquet"))
    arrow::write_parquet(cell_boundaries, file.path(fixture_dir, "cell_boundaries.parquet"))
    arrow::write_parquet(nucleus_boundaries, file.path(fixture_dir, "nucleus_boundaries.parquet"))
    invisible(fixture_dir)
}

.pc_legacy_object_snapshot <- function(obj) {
    list(
        coord = list(
            centroids = obj@coord$centroids,
            segmentation_cell = obj@coord$segmentation_cell,
            segmentation_nucleus = obj@coord$segmentation_nucleus
        ),
        grid = lapply(obj@grid, function(g) {
            list(
                grid_info = g$grid_info,
                counts = g$counts,
                grid_length = g$grid_length,
                xbins_eff = g$xbins_eff,
                ybins_eff = g$ybins_eff,
                has_W = !is.null(g$W),
                has_Xz = !is.null(g$Xz)
            )
        }),
        meta.data = obj@meta.data,
        stats = list(platform = obj@stats$platform)
    )
}

.pc_canon_formals <- function(fun) {
    lapply(formals(fun), function(x) paste(deparse(x), collapse = " "))
}

.pc_guess_package_root <- function(package_root = NULL) {
    if (!is.null(package_root)) {
        return(normalizePath(package_root, mustWork = TRUE))
    }
    cur <- normalizePath(getwd(), mustWork = TRUE)
    for (cand in c(cur, dirname(cur), dirname(dirname(cur)))) {
        if (file.exists(file.path(cand, "DESCRIPTION"))) {
            return(cand)
        }
    }
    stop("Could not infer `package_root`; please supply the geneSCOPE source directory.")
}

.pc_scalar_string <- function(x) {
    if (is.null(x)) return("")
    if (length(x) == 1L && !is.list(x)) return(as.character(x))
    paste(capture.output(str(x, give.attr = FALSE, vec.len = 8)), collapse = " ")
}
