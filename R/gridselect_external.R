.gridselect_coalesce <- function(x, default) {
    if (is.null(x)) default else x
}

.gridselect_bool <- function(x, default = FALSE) {
    if (is.null(x)) return(default)
    if (is.logical(x) && length(x) == 1) return(isTRUE(x))
    if (is.numeric(x) && length(x) == 1) return(!is.na(x) && x != 0)
    if (is.character(x) && length(x) == 1) {
        if (x %in% c("0", "FALSE", "false", "F", "f")) return(FALSE)
        if (x %in% c("1", "TRUE", "true", "T", "t")) return(TRUE)
    }
    isTRUE(as.logical(x)[1])
}

.gridselect_num <- function(x, default = NA_real_) {
    if (is.null(x)) return(default)
    suppressWarnings(as.numeric(x))[1]
}

.gridselect_path_is_auto <- function(x) {
    if (is.null(x) || !length(x)) return(TRUE)
    x <- as.character(x)[1]
    if (is.na(x)) return(TRUE)
    x <- trimws(x)
    if (!nzchar(x)) return(TRUE)
    tolower(x) %in% c("auto", "default")
}

.gridselect_is_abs_path <- function(path) {
    if (!is.character(path) || length(path) != 1 || !nzchar(path)) return(FALSE)
    path <- as.character(path)[1]
    if (startsWith(path, "~")) return(TRUE)
    if (.Platform$OS.type == "windows") {
        return(grepl("^[A-Za-z]:[/\\\\]", path) || startsWith(path, "\\\\") || startsWith(path, "//"))
    }
    startsWith(path, "/")
}

.gridselect_can_create_under <- function(path) {
    if (!is.character(path) || length(path) != 1 || !nzchar(path)) return(FALSE)
    path <- normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
    parent <- dirname(path)
    last <- ""
    while (!dir.exists(parent) && !identical(last, parent)) {
        last <- parent
        parent <- dirname(parent)
    }
    if (!dir.exists(parent)) return(FALSE)
    file.access(parent, 2) == 0
}

.gridselect_guess_checkout_root_from_repo_pin <- function(repo_pin) {
    if (!is.character(repo_pin) || length(repo_pin) != 1 || !nzchar(repo_pin)) return(NA_character_)
    if (!file.exists(repo_pin)) return(NA_character_)
    pin_dir <- dirname(normalizePath(repo_pin, winslash = "/", mustWork = TRUE))
    pkg_root <- dirname(pin_dir)
    has_desc <- file.exists(file.path(pkg_root, "DESCRIPTION"))
    has_r <- dir.exists(file.path(pkg_root, "R"))
    has_git <- file.exists(file.path(pkg_root, ".git"))
    if (isTRUE(has_desc) && isTRUE(has_r) && isTRUE(has_git)) return(pkg_root)
    NA_character_
}

.gridselect_default_repo_path <- function(repo_pin = NA_character_) {
    opt <- getOption("geneSCOPE.idelta_gridselect_repo_path", NULL)
    if (!.gridselect_path_is_auto(opt)) {
        opt <- as.character(opt)[1]
        opt <- trimws(opt)
        if (nzchar(opt)) return(opt)
    }

    env <- Sys.getenv("GENESCOPE_IDELTA_GRIDSELECT_REPO_PATH", unset = "")
    if (nzchar(env)) return(env)

    candidates <- character()

    checkout_root <- .gridselect_guess_checkout_root_from_repo_pin(repo_pin)
    if (is.character(checkout_root) && length(checkout_root) == 1 && nzchar(checkout_root) && dir.exists(checkout_root)) {
        candidates <- c(candidates, file.path(dirname(checkout_root), "idelta-gridselect"))
    }

    cache_root <- tryCatch(tools::R_user_dir("geneSCOPE", which = "cache"), error = function(e) NA_character_)
    if (is.character(cache_root) && length(cache_root) == 1 && !is.na(cache_root) && nzchar(cache_root)) {
        candidates <- c(candidates, file.path(cache_root, "idelta-gridselect"))
    }

    candidates <- c(candidates, file.path(tempdir(), "geneSCOPE", "idelta-gridselect"))

    for (cand in candidates) {
        if (.gridselect_can_create_under(cand)) return(cand)
    }
    candidates[[length(candidates)]]
}

.gridselect_resolve_repo_path <- function(repo_path, repo_pin = NA_character_) {
    repo_path <- if (.gridselect_path_is_auto(repo_path)) {
        .gridselect_default_repo_path(repo_pin = repo_pin)
    } else {
        as.character(repo_path)[1]
    }
    if (is.na(repo_path) || !nzchar(repo_path)) {
        repo_path <- .gridselect_default_repo_path(repo_pin = repo_pin)
    }

    repo_path <- trimws(repo_path)

    if (!.gridselect_is_abs_path(repo_path) &&
        is.character(repo_pin) && length(repo_pin) == 1 && nzchar(repo_pin) && file.exists(repo_pin)) {
        pin_dir <- dirname(normalizePath(repo_pin, winslash = "/", mustWork = TRUE))
        repo_path <- file.path(pin_dir, repo_path)
    }

    normalizePath(path.expand(repo_path), winslash = "/", mustWork = FALSE)
}

.gridselect_default_widths_prescreen_v2 <- function() {
    as.integer(c(1, 2, 4, 8, 16, 24, 32, 40, 50, 60, 80, 100, 120, 150, 200, 250, 300))
}

.gridselect_default_widths_anchor_v2 <- function() {
    as.integer(c(1:40, 45, 50, 55, 60, 65, 70, 80, 90, 100, 110, 120, 130, 140, 150, 175, 200, 225, 250, 275, 300))
}

.gridselect_write_widths_file <- function(widths_um, out_path, header = NULL) {
    widths_um <- suppressWarnings(as.integer(widths_um))
    widths_um <- widths_um[is.finite(widths_um)]
    widths_um <- widths_um[widths_um >= 1L]
    widths_um <- sort(unique(widths_um))
    if (length(widths_um) == 0) stop("widths file would be empty: ", out_path)
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    con <- file(out_path, open = "wt")
    on.exit(close(con), add = TRUE)
    if (is.character(header) && length(header) > 0) {
        for (h in header) writeLines(paste0("# ", h), con = con, sep = "\n", useBytes = TRUE)
    }
    writeLines(as.character(widths_um), con = con, sep = "\n", useBytes = TRUE)
    invisible(out_path)
}

.gridselect_wlog <- function(wrapper_log_path, msg, echo = TRUE) {
    dir.create(dirname(wrapper_log_path), recursive = TRUE, showWarnings = FALSE)
    line <- sprintf("[gridselect] %s %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
    if (isTRUE(echo)) cat(line, "\n", sep = "")
    cat(line, "\n", file = wrapper_log_path, append = TRUE, sep = "")
}

.gridselect_is_executable <- function(path) {
    if (!is.character(path) || length(path) != 1 || !nzchar(path)) return(FALSE)
    if (!file.exists(path)) return(FALSE)
    file.access(path, 1) == 0
}

.gridselect_find_repo_file_up <- function(rel_path, max_up = 6) {
    d <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
    for (i in 0:max_up) {
        cand <- file.path(d, rel_path)
        if (file.exists(cand)) return(cand)
        parent <- dirname(d)
        if (identical(parent, d)) break
        d <- parent
    }
    NA_character_
}

.gridselect_kv_read <- function(path) {
    stopifnot(is.character(path), length(path) == 1, nzchar(path))
    if (!file.exists(path)) stop("Missing repo pin file: ", path)
    lines <- readLines(path, warn = FALSE)
    out <- list()
    for (ln in lines) {
        ln <- trimws(ln)
        if (!nzchar(ln) || startsWith(ln, "#")) next
        eq <- regexpr("=", ln, fixed = TRUE)
        if (eq <= 0) next
        key <- trimws(substr(ln, 1, eq - 1))
        val <- trimws(substr(ln, eq + 1, nchar(ln)))
        if (!nzchar(key)) next
        out[[key]] <- val
    }
    out
}

.gridselect_resolve_repo_pin_path <- function(repo_pin = NULL) {
    if (is.character(repo_pin) && length(repo_pin) == 1 && nzchar(repo_pin)) {
        return(normalizePath(repo_pin, winslash = "/", mustWork = TRUE))
    }
    env_pin <- Sys.getenv("GENESCOPE_IDELTA_GRIDSELECT_REPO_PIN", unset = "")
    if (nzchar(env_pin) && file.exists(env_pin)) {
        return(normalizePath(env_pin, winslash = "/", mustWork = TRUE))
    }
    cand <- .gridselect_find_repo_file_up(file.path("docs", "idelta_gridselect_repo_pin.txt"), max_up = 6)
    if (is.character(cand) && length(cand) == 1 && nzchar(cand) && file.exists(cand)) {
        return(normalizePath(cand, winslash = "/", mustWork = TRUE))
    }
    pkg_pin <- system.file("docs", "idelta_gridselect_repo_pin.txt", package = "geneSCOPE")
    if (is.character(pkg_pin) && length(pkg_pin) == 1 && nzchar(pkg_pin) && file.exists(pkg_pin)) {
        return(normalizePath(pkg_pin, winslash = "/", mustWork = TRUE))
    }
    NA_character_
}

.gridselect_load_repo_cfg <- function(repo_pin = NULL) {
    pin_path <- .gridselect_resolve_repo_pin_path(repo_pin)
    cfg <- if (is.character(pin_path) && nzchar(pin_path) && file.exists(pin_path)) .gridselect_kv_read(pin_path) else list()

    list(
        repo_pin = if (is.character(pin_path) && nzchar(pin_path)) pin_path else NA_character_,
        repo_path = .gridselect_coalesce(cfg$REPO_PATH, "AUTO"),
        repo_url = .gridselect_coalesce(cfg$REPO_URL, "https://github.com/CoooRossa/idelta-gridselect-commitonly.git"),
        pin_ref = .gridselect_coalesce(cfg$PIN_REF, ""),
        auto_clone = .gridselect_coalesce(cfg$AUTO_CLONE, "1"),
        auto_pull = .gridselect_coalesce(cfg$AUTO_PULL, "0"),
        build_if_missing = .gridselect_coalesce(cfg$BUILD_IF_MISSING, "0"),
        bin_rel_path = .gridselect_coalesce(cfg$BIN_REL_PATH, "target/release/idelta-gridselect"),
        runner_rel_path = .gridselect_coalesce(cfg$RUNNER_REL_PATH, file.path("integrations", "geneSCOPE", "run_idelta_gridselect.sh"))
    )
}

.gridselect_repo_boot_log_path <- function(out_dir) {
    file.path(out_dir, "logs", "idelta_repo_bootstrap.log")
}

.gridselect_repo_boot_log <- function(path, msg) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    line <- sprintf("[idelta_repo_bootstrap] %s %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
    cat(line, "\n", file = path, append = TRUE, sep = "")
}

.gridselect_git <- function(args, log_path) {
    stopifnot(is.character(args))
    dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
    out <- suppressWarnings(system2("git", args = shQuote(args), stdout = TRUE, stderr = TRUE))
    if (!is.character(out)) out <- character(0)
    if (length(out) > 0) cat(paste0(out, "\n"), file = log_path, append = TRUE, sep = "")
    status <- attr(out, "status")
    if (!is.null(status) && is.finite(suppressWarnings(as.numeric(status))[1]) && suppressWarnings(as.numeric(status))[1] != 0) {
        stop("git failed (see logs): ", log_path)
    }
    invisible(TRUE)
}

.gridselect_ensure_idelta_repo <- function(cfg, out_dir) {
    cfg$repo_path <- .gridselect_resolve_repo_path(cfg$repo_path, repo_pin = cfg$repo_pin)
    log_path <- .gridselect_repo_boot_log_path(out_dir)

    .gridselect_repo_boot_log(log_path, paste0("repo_path=", cfg$repo_path))
    if (is.character(cfg$repo_pin) && nzchar(cfg$repo_pin)) .gridselect_repo_boot_log(log_path, paste0("repo_pin=", cfg$repo_pin))
    if (nzchar(cfg$repo_url)) .gridselect_repo_boot_log(log_path, paste0("repo_url=", cfg$repo_url))
    if (nzchar(cfg$pin_ref)) .gridselect_repo_boot_log(log_path, paste0("pin_ref=", cfg$pin_ref))
    .gridselect_repo_boot_log(
        log_path,
        paste0("auto_clone=", cfg$auto_clone, " auto_pull=", cfg$auto_pull, " build_if_missing=", cfg$build_if_missing)
    )

    repo_present <- dir.exists(cfg$repo_path) && dir.exists(file.path(cfg$repo_path, ".git"))
    if (!repo_present) {
        auto_clone <- .gridselect_bool(cfg$auto_clone, default = FALSE)
        if (!isTRUE(auto_clone)) {
            stop("External repo missing: ", cfg$repo_path, " (set AUTO_CLONE=1 in repo pin file to enable auto-clone).")
        }
        if (!nzchar(cfg$repo_url)) {
            stop("AUTO_CLONE=1 but REPO_URL is empty; set REPO_URL in the repo pin file.")
        }
        if (dir.exists(cfg$repo_path)) {
            # Avoid clobbering non-empty directories.
            kids <- list.files(cfg$repo_path, all.files = TRUE, no.. = TRUE)
            if (length(kids) > 0) stop("Repo path exists but is not a git repo: ", cfg$repo_path, " (directory is not empty).")
        } else {
            dir.create(dirname(cfg$repo_path), recursive = TRUE, showWarnings = FALSE)
        }

        .gridselect_repo_boot_log(log_path, "cloning repo")
        .gridselect_git(c("clone", cfg$repo_url, cfg$repo_path), log_path = log_path)
        repo_present <- dir.exists(file.path(cfg$repo_path, ".git"))
    }

    auto_pull <- .gridselect_bool(cfg$auto_pull, default = FALSE)
    if (isTRUE(auto_pull)) {
        .gridselect_repo_boot_log(log_path, "auto_pull=1: fetch --all (tags disabled)")
        .gridselect_git(c("-C", cfg$repo_path, "fetch", "--all"), log_path = log_path)
    }

    if (nzchar(cfg$pin_ref)) {
        if (!grepl("^[0-9a-f]{7,40}$", cfg$pin_ref, ignore.case = TRUE)) {
            stop("PIN_REF must be a git commit hash (7-40 hex); tag-based pinning is not allowed: ", cfg$pin_ref)
        }
        .gridselect_repo_boot_log(log_path, paste0("checkout pin_ref=", cfg$pin_ref))
        .gridselect_git(c("-C", cfg$repo_path, "checkout", cfg$pin_ref), log_path = log_path)
    } else if (isTRUE(auto_pull)) {
        .gridselect_repo_boot_log(log_path, "auto_pull=1: pull --ff-only")
        .gridselect_git(c("-C", cfg$repo_path, "pull", "--ff-only"), log_path = log_path)
    }

    cfg
}

.gridselect_run_idelta_gridselect_runner <- function(
    out_dir,
    xenium_dir = NULL,
    transcripts_path = NULL,
    roi_csv,
    threads = 32,
    knee_search_min_um = 5,
    descending_delta_um = 2,
    max_width_um = NULL,
    skip_molecules = FALSE,
    skip_tool = FALSE,
    repo_pin = NULL,
    strict = TRUE
) {
    has_xenium <- is.character(xenium_dir) && length(xenium_dir) == 1 && nzchar(xenium_dir)
    has_transcripts <- is.character(transcripts_path) && length(transcripts_path) == 1 && nzchar(transcripts_path)
    if (has_xenium && has_transcripts) stop("Provide only one of `xenium_dir` or `transcripts_path`.")
    if (!has_xenium && !has_transcripts) stop("Provide one of `xenium_dir` or `transcripts_path`.")
    if (!(is.character(roi_csv) && length(roi_csv) == 1 && nzchar(roi_csv))) stop("`roi_csv` is required.")
    if (!(is.character(out_dir) && length(out_dir) == 1 && nzchar(out_dir))) stop("`out_dir` is required.")

    out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    staging_dir <- file.path(out_dir, "staging")
    dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
    logs_dir <- file.path(out_dir, "logs")
    dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

    roi_csv_in <- normalizePath(roi_csv, winslash = "/", mustWork = TRUE)
    roi_vertices_path <- file.path(staging_dir, "roi_vertices.csv")
    if (!file.exists(roi_vertices_path) ||
        file.info(roi_vertices_path)$mtime < file.info(roi_csv_in)$mtime) {
        .gridselect_normalize_roi_vertices(roi_csv_in, roi_vertices_path)
    }

	    # Resolve/validate transcripts inputs early (before repo bootstrap).
	    if (has_xenium) {
	        xenium_dir <- .gridselect_resolve_xenium_dir(xenium_dir)
	        xenium_dir <- normalizePath(xenium_dir, winslash = "/", mustWork = FALSE)
	        if (!dir.exists(xenium_dir)) stop("xenium_dir not found: ", xenium_dir)
	        tx <- .gridselect_guess_xenium_file(xenium_dir, "transcripts")
	        if (!is.character(tx) || length(tx) != 1 || is.na(tx) || !nzchar(tx)) {
	            stop(
	                "Could not find Xenium transcripts under: ", xenium_dir,
	                "\nExpected one of: transcripts.csv.gz, transcripts.csv, transcripts.parquet"
	            )
	        }
	        if (grepl("\\.zarr(\\.zip)?$", tx, ignore.case = TRUE) || (dir.exists(tx) && grepl("\\.zarr$", tx, ignore.case = TRUE))) {
	            stop(
	                "Xenium transcripts appear to be Zarr (", basename(tx), "), which is not supported by the external idelta-gridselect runner.\n",
	                "Please export `transcripts.parquet` or `transcripts.csv.gz` and rerun, or provide `transcripts_path=` pointing to a CSV/Parquet file."
	            )
	        }
	    }

    if (has_transcripts) {
        transcripts_path <- normalizePath(transcripts_path, winslash = "/", mustWork = FALSE)
        if (!(file.exists(transcripts_path) || dir.exists(transcripts_path))) {
            stop("transcripts_path not found: ", transcripts_path)
        }
        if (grepl("\\.zarr(\\.zip)?$", transcripts_path, ignore.case = TRUE) ||
            (dir.exists(transcripts_path) && grepl("\\.zarr$", transcripts_path, ignore.case = TRUE))) {
            stop(
                "transcripts_path points to a Zarr store, which is not supported: ", transcripts_path,
                "\nProvide a `transcripts.parquet` or `transcripts.csv(.gz)` file instead."
            )
        }
    }

    cfg <- .gridselect_load_repo_cfg(repo_pin = repo_pin)
    cfg <- .gridselect_ensure_idelta_repo(cfg = cfg, out_dir = out_dir)

    runner <- normalizePath(file.path(cfg$repo_path, cfg$runner_rel_path), winslash = "/", mustWork = FALSE)
    if (!file.exists(runner)) {
        stop(
            "Runner script not found: ",
            runner,
            "\nExpected the external repo to contain RUNNER_REL_PATH (default: integrations/geneSCOPE/run_idelta_gridselect.sh)."
        )
    }

    threads_num <- suppressWarnings(as.numeric(threads))[1]
    if (!is.finite(threads_num) || threads_num < 1) stop("`threads` must be an integer >= 1.")
    threads_int <- as.integer(round(threads_num))

    knee_search_min_um_num <- suppressWarnings(as.numeric(knee_search_min_um))[1]
    descending_delta_um_num <- suppressWarnings(as.numeric(descending_delta_um))[1]
    if (!is.finite(knee_search_min_um_num) || abs(knee_search_min_um_num - round(knee_search_min_um_num)) > 1e-9) {
        stop("`knee_search_min_um` must be an integer.")
    }
    if (!is.finite(descending_delta_um_num) || abs(descending_delta_um_num - round(descending_delta_um_num)) > 1e-9) {
        stop("`descending_delta_um` must be an integer.")
    }
    knee_search_min_um_int <- as.integer(round(knee_search_min_um_num))
    descending_delta_um_int <- as.integer(round(descending_delta_um_num))
    if (knee_search_min_um_int < 0) stop("`knee_search_min_um` must be >= 0.")
    if (descending_delta_um_int < 0) stop("`descending_delta_um` must be >= 0.")

    max_width_um_num <- suppressWarnings(as.numeric(max_width_um))[1]
    max_width_um_arg <- NULL
    if (is.finite(max_width_um_num)) {
        if (abs(max_width_um_num - round(max_width_um_num)) > 1e-9) stop("`max_width_um` must be an integer (or NULL).")
        max_width_um_arg <- as.integer(round(max_width_um_num))
        if (max_width_um_arg < 1L) stop("`max_width_um` must be >= 1 (or NULL).")
    }

    args <- c(
        "--out", out_dir,
        "--roi_csv", roi_vertices_path,
        "--threads", as.character(threads_int),
        "--knee-search-min-um", as.character(knee_search_min_um_int),
        "--descending-delta-um", as.character(descending_delta_um_int)
    )
    if (!is.null(max_width_um_arg)) args <- c(args, "--max-width-um", as.character(max_width_um_arg))
    if (has_xenium) args <- c(args, "--xenium_dir", xenium_dir)
    if (has_transcripts) args <- c(args, "--transcripts", transcripts_path)
    if (.gridselect_bool(skip_molecules, default = FALSE)) args <- c(args, "--skip_molecules")
    if (.gridselect_bool(skip_tool, default = FALSE)) args <- c(args, "--skip_tool")
    if (is.character(cfg$repo_pin) && nzchar(cfg$repo_pin)) args <- c(args, "--repo_pin", cfg$repo_pin)

    runner_log <- file.path(logs_dir, "run_idelta_gridselect.log")
    out <- suppressWarnings(system2("sh", args = shQuote(c(runner, args)), stdout = TRUE, stderr = runner_log))
    if (!is.character(out)) out <- character(0)
    status <- attr(out, "status")
    if (!is.null(status) && is.finite(suppressWarnings(as.numeric(status))[1]) && suppressWarnings(as.numeric(status))[1] != 0) {
        stop("idelta-gridselect runner failed (see logs): ", runner_log)
    }

    rec_line <- out[grepl("^RECOMMENDED_GRID_UM=", out)]
    recommended_grid_um_runner <- NA_real_
    if (length(rec_line) == 1) {
        recommended_grid_um_runner <- suppressWarnings(as.numeric(sub("^RECOMMENDED_GRID_UM=", "", rec_line)))[1]
    } else if (isTRUE(strict)) {
        warning("Could not parse RECOMMENDED_GRID_UM from runner stdout; see: ", runner_log)
    }

    list(
        runner = runner,
        runner_log = runner_log,
        runner_stdout = out,
        recommended_grid_um_runner = recommended_grid_um_runner,
        repo_path = cfg$repo_path,
        repo_pin = cfg$repo_pin
    )
}

.gridselect_locate_python_helper <- function() {
    # Installed package path (inst/scripts -> <pkg>/scripts)
    p <- system.file("scripts", "gridselect_make_molecules_bbox.py", package = "geneSCOPE")
    if (is.character(p) && nzchar(p) && file.exists(p)) return(p)

    # Repo-local paths (common dev usage).
    cand <- .gridselect_find_repo_file_up(file.path("scripts", "gridselect_make_molecules_bbox.py"))
    if (is.character(cand) && nzchar(cand) && file.exists(cand)) return(cand)

    cand2 <- .gridselect_find_repo_file_up(file.path("inst", "scripts", "gridselect_make_molecules_bbox.py"))
    if (is.character(cand2) && nzchar(cand2) && file.exists(cand2)) return(cand2)

    NA_character_
}

.gridselect_normalize_roi_vertices <- function(roi_csv_in, roi_csv_out) {
    lines <- readLines(roi_csv_in, warn = FALSE)
    lines <- lines[!grepl("^\\s*#", lines)]
    lines <- lines[nzchar(trimws(lines))]
    if (length(lines) < 2) stop("ROI file had no data rows after stripping comment/blank lines: ", roi_csv_in)

    delim <- if (grepl("\t", lines[[1]], fixed = TRUE)) "\t" else ","
    df <- utils::read.table(
        text = paste(lines, collapse = "\n"),
        header = TRUE,
        sep = delim,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        comment.char = ""
    )
    if (nrow(df) < 3) stop("ROI needs >= 3 vertices: ", roi_csv_in)

    headers_low <- tolower(colnames(df))
    idx_x <- match(TRUE, headers_low %in% c("x_um", "x"))
    idx_y <- match(TRUE, headers_low %in% c("y_um", "y"))
    if (is.na(idx_x) || is.na(idx_y)) {
        stop("ROI header must include x_um,y_um (or X,Y). Got: ", paste(colnames(df), collapse = ", "))
    }

    x <- suppressWarnings(as.numeric(df[[idx_x]]))
    y <- suppressWarnings(as.numeric(df[[idx_y]]))
    if (any(!is.finite(x)) || any(!is.finite(y))) stop("ROI vertices contain non-numeric x/y values: ", roi_csv_in)

    out <- data.frame(x_um = x, y_um = y, check.names = FALSE)
    dir.create(dirname(roi_csv_out), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(out, file = roi_csv_out, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
    invisible(roi_csv_out)
}

.gridselect_guess_xenium_file <- function(xenium_dir, base_name) {
    if (!is.character(xenium_dir) || length(xenium_dir) != 1 || !nzchar(xenium_dir)) return(NA_character_)
    d <- normalizePath(xenium_dir, winslash = "/", mustWork = FALSE)
    if (!dir.exists(d)) return(NA_character_)
    cand <- c(
        file.path(d, paste0(base_name, ".csv.gz")),
        file.path(d, paste0(base_name, ".csv")),
        file.path(d, paste0(base_name, ".parquet")),
        file.path(d, paste0(base_name, ".pq")),
        file.path(d, paste0(base_name, ".zarr.zip")),
        file.path(d, paste0(base_name, ".zarr"))
    )
    cand <- cand[file.exists(cand) | dir.exists(cand)]
    if (length(cand) == 0) return(NA_character_)
    normalizePath(cand[[1]], winslash = "/", mustWork = TRUE)
}

.gridselect_resolve_xenium_dir <- function(xenium_dir) {
    if (!is.character(xenium_dir) || length(xenium_dir) != 1 || !nzchar(xenium_dir)) {
        return(xenium_dir)
    }
    d <- normalizePath(xenium_dir, winslash = "/", mustWork = FALSE)
    if (!dir.exists(d)) return(d)
    tx <- .gridselect_guess_xenium_file(d, "transcripts")
    if (is.character(tx) && nzchar(tx) && !is.na(tx)) return(d)
    outs <- file.path(d, "outs")
    if (dir.exists(outs)) {
        tx2 <- .gridselect_guess_xenium_file(outs, "transcripts")
        if (is.character(tx2) && nzchar(tx2) && !is.na(tx2)) {
            return(normalizePath(outs, winslash = "/", mustWork = TRUE))
        }
    }
    d
}

.gridselect_bbox_from_xy_table <- function(
    path,
    x_candidates = c("x_um", "x", "x_centroid", "x_location"),
    y_candidates = c("y_um", "y", "y_centroid", "y_location"),
    strict = TRUE
) {
    if (!is.character(path) || length(path) != 1 || !nzchar(path)) stop("bbox source path is required.")
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)

    if (grepl("\\.zarr(\\.zip)?$", path, ignore.case = TRUE) ||
        (dir.exists(path) && grepl("\\.zarr$", path, ignore.case = TRUE))) {
        stop(
            "Zarr transcripts/cells are not supported for bbox auto-detection: ", path,
            "\nProvide `roi_csv=` explicitly, or provide `bbox_csv=` pointing to a CSV/TSV (or Parquet; requires package 'arrow')."
        )
    }

    if (grepl("\\.(parquet|pq)$", path, ignore.case = TRUE)) {
        if (!requireNamespace("arrow", quietly = TRUE)) {
            stop("bbox source is Parquet but package 'arrow' is required: ", path)
        }
        ds <- tryCatch(arrow::open_dataset(path), error = function(e) e)
        if (inherits(ds, "error")) {
            stop("Failed to open Parquet for bbox scan: ", path, " (", conditionMessage(ds), ")")
        }

        cols <- tryCatch(names(ds), error = function(e) character())
        if (!length(cols)) stop("Failed to read Parquet schema for bbox scan: ", path)
        cols_low <- tolower(cols)
        x_idx <- match(TRUE, cols_low %in% tolower(x_candidates))
        y_idx <- match(TRUE, cols_low %in% tolower(y_candidates))
        if (is.na(x_idx) || is.na(y_idx)) {
            stop(
                "Parquet bbox source is missing x/y columns.\nExpected x in: ",
                paste(x_candidates, collapse = ", "),
                "\nExpected y in: ",
                paste(y_candidates, collapse = ", "),
                "\nAvailable: ",
                paste(cols, collapse = ", "),
                "\nFile: ",
                path
            )
        }
        x_col <- cols[[x_idx]]
        y_col <- cols[[y_idx]]

        reader <- tryCatch({
            if (inherits(ds, "Dataset")) {
                scan_builder <- ds$NewScan()
                scan_builder$Project(c(x_col, y_col))
                scan_builder$UseThreads(TRUE)
                scan_builder$BatchSize(500000L)
                scan_builder$Finish()$ToRecordBatchReader()
            } else {
                arrow::as_record_batch_reader(ds)
            }
        }, error = function(e) e)
        if (inherits(reader, "error")) {
            stop("Failed to create Arrow reader for bbox scan: ", path, " (", conditionMessage(reader), ")")
        }

        xmin <- Inf
        xmax <- -Inf
        ymin <- Inf
        ymax <- -Inf
        n <- 0L

        repeat {
            batch <- tryCatch(reader$read_next_batch(), error = function(e) e)
            if (inherits(batch, "error")) {
                stop("Parquet bbox scan failed: ", path, " (", conditionMessage(batch), ")")
            }
            if (is.null(batch)) break
            tbl <- tryCatch(as.data.frame(batch), error = function(e) NULL)
            if (is.null(tbl) || nrow(tbl) == 0) next

            x <- suppressWarnings(as.numeric(tbl[[x_col]]))
            y <- suppressWarnings(as.numeric(tbl[[y_col]]))
            ok <- is.finite(x) & is.finite(y)
            if (!any(ok)) next

            xmin_b <- suppressWarnings(min(x[ok], na.rm = TRUE))
            xmax_b <- suppressWarnings(max(x[ok], na.rm = TRUE))
            ymin_b <- suppressWarnings(min(y[ok], na.rm = TRUE))
            ymax_b <- suppressWarnings(max(y[ok], na.rm = TRUE))
            if (is.finite(xmin_b)) xmin <- min(xmin, xmin_b)
            if (is.finite(xmax_b)) xmax <- max(xmax, xmax_b)
            if (is.finite(ymin_b)) ymin <- min(ymin, ymin_b)
            if (is.finite(ymax_b)) ymax <- max(ymax, ymax_b)
            n <- n + sum(ok)
        }

        if (n == 0L || any(!is.finite(c(xmin, xmax, ymin, ymax)))) {
            stop("Parquet bbox source had no numeric x/y rows: ", path)
        }

        return(list(
            min_x = xmin,
            max_x = xmax,
            min_y = ymin,
            max_y = ymax,
            n = n,
            x_col = tolower(x_col),
            y_col = tolower(y_col),
            source = path
        ))
    }

    py <- Sys.which("python3")
    py_ok <- nzchar(py)
    if (py_ok) {
        x_cand <- paste(tolower(x_candidates), collapse = ",")
        y_cand <- paste(tolower(y_candidates), collapse = ",")
        py_code <- paste(
            "import sys,csv,gzip,math",
            "path=sys.argv[1]",
            "x_cands=[c for c in sys.argv[2].split(',') if c]",
            "y_cands=[c for c in sys.argv[3].split(',') if c]",
            "def open_text(p):",
            "    return gzip.open(p,'rt',newline='') if p.lower().endswith('.gz') else open(p,'r',newline='')",
            "def iter_data_lines(fp):",
            "    for ln in fp:",
            "        if not ln.strip():",
            "            continue",
            "        if ln.lstrip().startswith('#'):",
            "            continue",
            "        yield ln",
            "sample_lines=[]",
            "with open_text(path) as f:",
            "    for ln in iter_data_lines(f):",
            "        sample_lines.append(ln)",
            "        if len(sample_lines)>=50:",
            "            break",
            "if not sample_lines:",
            "    sys.stderr.write('empty after stripping comments\\n'); sys.exit(2)",
            "sample=''.join(sample_lines)",
            "try:",
            "    dialect=csv.Sniffer().sniff(sample,delimiters=[',','\\t'])",
            "except csv.Error:",
            "    class D(csv.Dialect):",
            "        delimiter=','; quotechar='\"'; doublequote=True; skipinitialspace=False; lineterminator='\\n'; quoting=csv.QUOTE_MINIMAL",
            "    dialect=D()",
            "with open_text(path) as f:",
            "    r=csv.reader(iter_data_lines(f),dialect=dialect)",
            "    hdr=next(r)",
            "    hdr_low=[h.strip().strip('\"').lower() for h in hdr]",
            "    x_idx=None; y_idx=None",
            "    for c in x_cands:",
            "        if c in hdr_low:",
            "            x_idx=hdr_low.index(c); break",
            "    for c in y_cands:",
            "        if c in hdr_low:",
            "            y_idx=hdr_low.index(c); break",
            "    if x_idx is None or y_idx is None:",
            "        sys.stderr.write('missing x/y columns; header=' + ','.join(hdr_low) + '\\n'); sys.exit(3)",
            "    min_x=math.inf; max_x=-math.inf; min_y=math.inf; max_y=-math.inf; n=0",
            "    for row in r:",
            "        if len(row)<=max(x_idx,y_idx):",
            "            continue",
            "        try:",
            "            x=float(row[x_idx]); y=float(row[y_idx])",
            "        except Exception:",
            "            continue",
            "        if math.isfinite(x) and math.isfinite(y):",
            "            if x<min_x: min_x=x",
            "            if x>max_x: max_x=x",
            "            if y<min_y: min_y=y",
            "            if y>max_y: max_y=y",
            "            n+=1",
            "    if n==0 or not (math.isfinite(min_x) and math.isfinite(min_y) and math.isfinite(max_x) and math.isfinite(max_y)):",
            "        sys.stderr.write('no numeric coordinate rows\\n'); sys.exit(4)",
            "    print(f'{min_x}\\t{max_x}\\t{min_y}\\t{max_y}\\t{n}\\t{hdr_low[x_idx]}\\t{hdr_low[y_idx]}')",
            sep = "\n"
        )
        err_log <- tempfile(pattern = "gridselect_bbox_", fileext = ".log")
        on.exit(if (file.exists(err_log)) unlink(err_log), add = TRUE)
        out <- suppressWarnings(system2(
            py,
            args = shQuote(c("-c", py_code, path, x_cand, y_cand)),
            stdout = TRUE,
            stderr = err_log
        ))
        status <- attr(out, "status")
        if (is.null(status) || suppressWarnings(as.numeric(status))[1] == 0) {
            parts <- strsplit(out[[1]], "\t", fixed = TRUE)[[1]]
            if (length(parts) >= 4) {
                min_x <- suppressWarnings(as.numeric(parts[[1]]))[1]
                max_x <- suppressWarnings(as.numeric(parts[[2]]))[1]
                min_y <- suppressWarnings(as.numeric(parts[[3]]))[1]
                max_y <- suppressWarnings(as.numeric(parts[[4]]))[1]
                n <- if (length(parts) >= 5) suppressWarnings(as.integer(parts[[5]]))[1] else NA_integer_
                x_col <- if (length(parts) >= 6) parts[[6]] else NA_character_
                y_col <- if (length(parts) >= 7) parts[[7]] else NA_character_
                if (all(is.finite(c(min_x, max_x, min_y, max_y))) && min_x < max_x && min_y < max_y) {
                    return(list(
                        min_x = min_x,
                        max_x = max_x,
                        min_y = min_y,
                        max_y = max_y,
                        n = n,
                        x_col = x_col,
                        y_col = y_col,
                        source = path
                    ))
                }
            }
        } else if (isTRUE(strict)) {
            err_tail <- tryCatch(readLines(err_log, warn = FALSE), error = function(e) character(0))
            err_tail <- tail(err_tail, 50)
            warning("bbox scan via python3 failed for: ", path, "\n", paste(err_tail, collapse = "\n"))
        }
    }

    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("python3 bbox scan failed/unavailable and package 'data.table' is not available.")
    }
    hdr <- tryCatch(
        data.table::fread(path, nrows = 0L, showProgress = FALSE, comment.char = "#"),
        error = function(e) NULL
    )
    if (is.null(hdr)) stop("Failed to read header for bbox source: ", path)
    headers <- names(hdr)
    headers_low <- tolower(headers)
    x_idx <- match(TRUE, headers_low %in% tolower(x_candidates))
    y_idx <- match(TRUE, headers_low %in% tolower(y_candidates))
    if (is.na(x_idx) || is.na(y_idx)) stop("bbox source missing x/y columns: ", path)
    x_col <- headers[[x_idx]]
    y_col <- headers[[y_idx]]
    dt <- data.table::fread(path, select = c(x_col, y_col), showProgress = FALSE, comment.char = "#")
    x <- suppressWarnings(as.numeric(dt[[x_col]]))
    y <- suppressWarnings(as.numeric(dt[[y_col]]))
    keep <- is.finite(x) & is.finite(y)
    if (!any(keep)) stop("bbox source had no numeric x/y rows: ", path)
    list(
        min_x = min(x[keep]),
        max_x = max(x[keep]),
        min_y = min(y[keep]),
        max_y = max(y[keep]),
        n = sum(keep),
        x_col = tolower(x_col),
        y_col = tolower(y_col),
        source = path
    )
}

.gridselect_write_roi_bbox_vertices <- function(min_x, max_x, min_y, max_y, out_path) {
    if (!all(is.finite(c(min_x, max_x, min_y, max_y)))) stop("bbox must be finite.")
    if (!(min_x < max_x && min_y < max_y)) stop("bbox min/max is invalid (min_x < max_x and min_y < max_y are required).")
    out <- data.frame(
        x_um = c(min_x, max_x, max_x, min_x),
        y_um = c(min_y, min_y, max_y, max_y),
        check.names = FALSE
    )
    dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.table(out, file = out_path, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)
    invisible(out_path)
}

.gridselect_resolve_roi_csv_for_tool_run <- function(
    roi_csv = NULL,
    xenium_dir = NULL,
    transcripts_path = NULL,
    staging_dir,
    wrapper_log_path = NULL,
    bbox_csv = NULL,
    bbox_pad_um = 0,
    strict = TRUE
) {
    if (is.character(roi_csv) && length(roi_csv) == 1 && nzchar(roi_csv)) {
        roi_csv <- normalizePath(roi_csv, winslash = "/", mustWork = TRUE)
        return(roi_csv)
    }

    bbox_source <- NULL
    bbox_source_label <- NA_character_

    if (is.character(bbox_csv) && length(bbox_csv) == 1 && nzchar(bbox_csv)) {
        bbox_source <- normalizePath(bbox_csv, winslash = "/", mustWork = TRUE)
        bbox_source_label <- "bbox_csv"
    } else if (is.character(xenium_dir) && length(xenium_dir) == 1 && nzchar(xenium_dir)) {
        cells_path <- .gridselect_guess_xenium_file(xenium_dir, "cells")
        if (is.character(cells_path) && nzchar(cells_path) && file.exists(cells_path)) {
            bbox_source <- cells_path
            bbox_source_label <- "xenium_dir:cells"
        } else {
            tx_path <- .gridselect_guess_xenium_file(xenium_dir, "transcripts")
            if (is.character(tx_path) && nzchar(tx_path) && file.exists(tx_path)) {
                bbox_source <- tx_path
                bbox_source_label <- "xenium_dir:transcripts"
            }
        }
    } else if (is.character(transcripts_path) && length(transcripts_path) == 1 && nzchar(transcripts_path)) {
        bbox_source <- normalizePath(transcripts_path, winslash = "/", mustWork = TRUE)
        bbox_source_label <- "transcripts_path"
    }

    if (is.null(bbox_source) || !nzchar(bbox_source) || !file.exists(bbox_source)) {
        stop("When running without `roi_csv`, provide `bbox_csv=` (any CSV/TSV with x/y columns) or `xenium_dir=`/`transcripts_path=`.")
    }

    bbox_pad_um_num <- suppressWarnings(as.numeric(bbox_pad_um))[1]
    if (!is.finite(bbox_pad_um_num) || bbox_pad_um_num < 0) stop("`bbox_pad_um` must be a finite number >= 0.")

    bbox <- .gridselect_bbox_from_xy_table(path = bbox_source, strict = strict)
    min_x <- bbox$min_x - bbox_pad_um_num
    max_x <- bbox$max_x + bbox_pad_um_num
    min_y <- bbox$min_y - bbox_pad_um_num
    max_y <- bbox$max_y + bbox_pad_um_num

    dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
    auto_roi_path <- file.path(staging_dir, "roi_auto_bbox.csv")
    .gridselect_write_roi_bbox_vertices(min_x, max_x, min_y, max_y, out_path = auto_roi_path)

    if (is.character(wrapper_log_path) && nzchar(wrapper_log_path)) {
        .gridselect_wlog(
            wrapper_log_path,
            paste0(
                "auto_roi_bbox source=", bbox_source_label,
                " file=", bbox_source,
                " cols=", bbox$x_col, "/", bbox$y_col,
                " n=", ifelse(is.finite(bbox$n), bbox$n, "NA"),
                " pad_um=", bbox_pad_um_num,
                " x=[", sprintf("%.2f", min_x), ",", sprintf("%.2f", max_x),
                "] y=[", sprintf("%.2f", min_y), ",", sprintf("%.2f", max_y), "]"
            ),
            echo = FALSE
        )
    }

    auto_roi_path
}

.gridselect_autofind_bin <- function(search_roots = character(0), max_depth = 6) {
    # Fast explicit candidates first.
    candidates <- c(
        file.path(getwd(), "target", "release", "idelta-gridselect"),
        file.path(getwd(), "..", "idelta-gridselect", "target", "release", "idelta-gridselect")
    )
    for (c in candidates) {
        if (.gridselect_is_executable(c)) return(normalizePath(c, winslash = "/", mustWork = TRUE))
    }

    roots <- character(0)
    if (length(search_roots) > 0) roots <- c(roots, search_roots)
    default_root <- file.path(path.expand("~"), "Documents", "Code")
    if (dir.exists(default_root)) roots <- c(roots, default_root)

    roots <- unique(roots)
    roots <- roots[nzchar(roots)]
    for (root in roots) {
        root <- normalizePath(root, winslash = "/", mustWork = FALSE)
        if (!dir.exists(root)) next

        # BFS over directories up to max_depth; check for */idelta-gridselect/target/release/idelta-gridselect
        queue_path <- root
        queue_depth <- 0L
        head <- 1L
        while (head <= length(queue_path)) {
            d <- queue_path[[head]]
            depth <- queue_depth[[head]]
            head <- head + 1L

            cand_a <- file.path(d, "idelta-gridselect", "target", "release", "idelta-gridselect")
            if (.gridselect_is_executable(cand_a)) return(normalizePath(cand_a, winslash = "/", mustWork = TRUE))
            cand_b <- file.path(d, "target", "release", "idelta-gridselect")
            if (.gridselect_is_executable(cand_b)) return(normalizePath(cand_b, winslash = "/", mustWork = TRUE))

            if (depth >= max_depth) next
            children <- list.files(d, full.names = TRUE, recursive = FALSE, all.files = FALSE, include.dirs = TRUE)
            if (length(children) == 0) next
            info <- suppressWarnings(file.info(children))
            dirs <- rownames(info)[isTRUE(info$isdir)]
            if (length(dirs) == 0) next

            base <- basename(dirs)
            keep <- !(base %in% c(".git", "target", "node_modules", "__pycache__")) & !startsWith(base, ".")
            dirs <- dirs[keep]
            if (length(dirs) == 0) next

            queue_path <- c(queue_path, dirs)
            queue_depth <- c(queue_depth, rep.int(depth + 1L, length(dirs)))
        }
    }

    NA_character_
}

.gridselect_resolve_tool <- function(
    out_dir,
    tool_bin = NULL,
    docker_image = NULL,
    apptainer_sif = NULL,
    auto_find = TRUE,
    search_roots = character(0),
    resolve_path,
    wrapper_log_path
) {
    # Priority rules: args first, then env vars, then auto-find.
    tool_bin_arg <- tool_bin
    if (is.character(tool_bin_arg) && nzchar(tool_bin_arg)) {
        tool_bin_arg <- normalizePath(tool_bin_arg, winslash = "/", mustWork = FALSE)
    } else {
        tool_bin_arg <- NULL
    }

    env_bin <- Sys.getenv("IDELTA_GRIDSELECT_BIN", unset = "")
    if (!nzchar(env_bin)) env_bin <- NULL

    docker_arg <- docker_image
    if (!(is.character(docker_arg) && nzchar(docker_arg))) docker_arg <- NULL
    env_docker <- Sys.getenv("IDELTA_GRIDSELECT_DOCKER_IMAGE", unset = "")
    if (!nzchar(env_docker)) env_docker <- NULL

    sif_arg <- apptainer_sif
    if (is.character(sif_arg) && nzchar(sif_arg)) {
        sif_arg <- normalizePath(sif_arg, winslash = "/", mustWork = FALSE)
    } else {
        sif_arg <- NULL
    }
    env_sif <- Sys.getenv("IDELTA_GRIDSELECT_APPTAINER_SIF", unset = "")
    if (nzchar(env_sif)) env_sif <- normalizePath(env_sif, winslash = "/", mustWork = FALSE) else env_sif <- NULL

    resolved <- list(mode = NA_character_, source = NA_character_, command = NA_character_, base_args = character(0), tool_bin = NA_character_)

    if (!is.null(tool_bin_arg)) {
        if (!.gridselect_is_executable(tool_bin_arg)) stop("--tool_bin was provided but is not executable: ", tool_bin_arg)
        resolved$mode <- "bin"
        resolved$source <- "arg:tool_bin"
        resolved$command <- tool_bin_arg
        resolved$tool_bin <- tool_bin_arg
    } else if (!is.null(env_bin)) {
        if (!.gridselect_is_executable(env_bin)) stop("IDELTA_GRIDSELECT_BIN is set but is not executable: ", env_bin)
        resolved$mode <- "bin"
        resolved$source <- "env:IDELTA_GRIDSELECT_BIN"
        resolved$command <- env_bin
        resolved$tool_bin <- env_bin
    } else if (nzchar(Sys.which("idelta-gridselect"))) {
        path_bin <- Sys.which("idelta-gridselect")
        if (.gridselect_is_executable(path_bin)) {
            resolved$mode <- "bin"
            resolved$source <- "path:idelta-gridselect"
            resolved$command <- path_bin
            resolved$tool_bin <- path_bin
        }
    } else if (!is.null(docker_arg)) {
        resolved$mode <- "docker"
        resolved$source <- "arg:docker_image"
        resolved$command <- "docker"
        resolved$base_args <- c(
            "run",
            "--rm",
            "-v", paste0(normalizePath(out_dir, winslash = "/", mustWork = TRUE), ":/out"),
            "-w", "/out",
            "--entrypoint", "idelta-gridselect",
            docker_arg
        )
    } else if (!is.null(env_docker)) {
        resolved$mode <- "docker"
        resolved$source <- "env:IDELTA_GRIDSELECT_DOCKER_IMAGE"
        resolved$command <- "docker"
        resolved$base_args <- c(
            "run",
            "--rm",
            "-v", paste0(normalizePath(out_dir, winslash = "/", mustWork = TRUE), ":/out"),
            "-w", "/out",
            "--entrypoint", "idelta-gridselect",
            env_docker
        )
    } else if (!is.null(sif_arg)) {
        if (!file.exists(sif_arg)) stop("--apptainer_sif was provided but file not found: ", sif_arg)
        resolved$mode <- "apptainer"
        resolved$source <- "arg:apptainer_sif"
        resolved$command <- if (nzchar(Sys.which("apptainer"))) "apptainer" else if (nzchar(Sys.which("singularity"))) "singularity" else NA_character_
        if (is.na(resolved$command)) stop("apptainer/singularity not found but apptainer_sif was provided.")
        resolved$base_args <- c("run", "--bind", paste0(normalizePath(out_dir, winslash = "/", mustWork = TRUE), ":/out"), sif_arg)
    } else if (!is.null(env_sif)) {
        if (!file.exists(env_sif)) stop("IDELTA_GRIDSELECT_APPTAINER_SIF is set but file not found: ", env_sif)
        resolved$mode <- "apptainer"
        resolved$source <- "env:IDELTA_GRIDSELECT_APPTAINER_SIF"
        resolved$command <- if (nzchar(Sys.which("apptainer"))) "apptainer" else if (nzchar(Sys.which("singularity"))) "singularity" else NA_character_
        if (is.na(resolved$command)) stop("apptainer/singularity not found but IDELTA_GRIDSELECT_APPTAINER_SIF is set.")
        resolved$base_args <- c("run", "--bind", paste0(normalizePath(out_dir, winslash = "/", mustWork = TRUE), ":/out"), env_sif)
    } else if (isTRUE(auto_find)) {
        found <- .gridselect_autofind_bin(search_roots = search_roots, max_depth = 6)
        if (is.character(found) && nzchar(found) && .gridselect_is_executable(found)) {
            resolved$mode <- "bin"
            resolved$source <- "auto_find"
            resolved$command <- found
            resolved$tool_bin <- found
        }
    }

    if (is.na(resolved$mode) || !nzchar(resolved$mode)) {
        stop("No external tool configured. Provide `tool_bin=` (preferred) or set IDELTA_GRIDSELECT_BIN, or use docker/apptainer options.")
    }

    dir.create(dirname(resolve_path), recursive = TRUE, showWarnings = FALSE)
    lines <- c(
        paste0("resolved_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste0("mode=", resolved$mode),
        paste0("source=", resolved$source)
    )
    if (identical(resolved$mode, "bin")) {
        lines <- c(lines, paste0("tool_bin=", resolved$tool_bin))
        .gridselect_wlog(wrapper_log_path, paste0("RESOLVED_TOOL_BIN=", resolved$tool_bin), echo = TRUE)
    } else if (identical(resolved$mode, "docker")) {
        lines <- c(lines, paste0("docker_cmd=", paste(c(resolved$command, resolved$base_args), collapse = " ")))
        .gridselect_wlog(wrapper_log_path, "execution_mode=docker", echo = TRUE)
    } else if (identical(resolved$mode, "apptainer")) {
        lines <- c(lines, paste0("apptainer_cmd=", paste(c(resolved$command, resolved$base_args), collapse = " ")))
        .gridselect_wlog(wrapper_log_path, "execution_mode=apptainer", echo = TRUE)
    }
    if (length(search_roots) > 0) lines <- c(lines, paste0("search_roots=", paste(search_roots, collapse = " ")))
    writeLines(lines, con = resolve_path, useBytes = TRUE)

    resolved
}

.gridselect_make_molecules_bbox <- function(
    roi_vertices_csv,
    transcripts_path,
    out_tsv,
    wrapper_log_path,
    python_bin = "python3"
) {
    py_script <- .gridselect_locate_python_helper()
    if (!nzchar(py_script) || !file.exists(py_script)) {
        stop("Could not locate python helper `gridselect_make_molecules_bbox.py`. If using the repo, run from the repo root or reinstall geneSCOPE.")
    }

    args <- c(
        py_script,
        "--roi_csv", roi_vertices_csv,
        "--transcripts_gz", transcripts_path,
        "--out_tsv", out_tsv
    )
    mols_log_path <- file.path(dirname(wrapper_log_path), "make_molecules_bbox.log")
    .gridselect_wlog(
        wrapper_log_path,
        paste0("make_molecules_bbox python=", python_bin, " log=", mols_log_path),
        echo = TRUE
    )
    if (file.exists(mols_log_path)) file.remove(mols_log_path)
    python_cmd <- python_bin
    if (grepl("[[:space:]]", python_cmd)) python_cmd <- shQuote(python_cmd)
    status <- suppressWarnings(system2(python_cmd, args = shQuote(args), stdout = mols_log_path, stderr = mols_log_path))
    if (!is.numeric(status) || status != 0) {
        stop("bbox molecules generation failed (see logs): ", mols_log_path)
    }
    invisible(out_tsv)
}

.gridselect_default_exclude_prefixes <- function() {
    c(
        "Unassigned",
        "NegControl",
        "Background",
        "DeprecatedCodeword",
        "SystemControl",
        "Negative",
        "BlankCodeword",
        "Blank"
    )
}

.gridselect_roi_bbox_from_vertices <- function(roi_vertices_csv) {
    if (!is.character(roi_vertices_csv) || length(roi_vertices_csv) != 1 || !nzchar(roi_vertices_csv)) {
        stop("roi_vertices_csv must be a single path.")
    }
    roi_vertices_csv <- normalizePath(roi_vertices_csv, winslash = "/", mustWork = TRUE)
    df <- tryCatch(
        utils::read.table(
            roi_vertices_csv,
            header = TRUE,
            sep = ",",
            stringsAsFactors = FALSE,
            check.names = FALSE,
            comment.char = ""
        ),
        error = function(e) NULL
    )
    if (is.null(df) || !nrow(df)) stop("ROI vertices file is empty: ", roi_vertices_csv)

    cn_low <- tolower(names(df))
    x_col <- names(df)[match(TRUE, cn_low %in% c("x_um", "x"))]
    y_col <- names(df)[match(TRUE, cn_low %in% c("y_um", "y"))]
    if (!length(x_col) || is.na(x_col) || !length(y_col) || is.na(y_col)) {
        stop("ROI vertices must have columns x_um/y_um (or x/y): ", roi_vertices_csv)
    }

    x <- suppressWarnings(as.numeric(df[[x_col]]))
    y <- suppressWarnings(as.numeric(df[[y_col]]))
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L) stop("ROI vertices must contain >= 3 finite points: ", roi_vertices_csv)

    list(
        min_x = min(x[ok]),
        max_x = max(x[ok]),
        min_y = min(y[ok]),
        max_y = max(y[ok])
    )
}

.gridselect_make_molecules_bbox_parquet_arrow <- function(
    roi_vertices_csv,
    transcripts_path,
    out_tsv,
    wrapper_log_path,
    exclude_prefixes = .gridselect_default_exclude_prefixes(),
    batch_size = 500000L
) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("Package 'arrow' is required to stream Parquet transcripts for molecules bbox extraction.")
    }
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required to write molecules bbox TSV efficiently.")
    }
    roi_vertices_csv <- normalizePath(roi_vertices_csv, winslash = "/", mustWork = TRUE)
    transcripts_path <- normalizePath(transcripts_path, winslash = "/", mustWork = TRUE)
    out_tsv <- normalizePath(out_tsv, winslash = "/", mustWork = FALSE)
    dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

    bb <- .gridselect_roi_bbox_from_vertices(roi_vertices_csv)
    min_x <- bb$min_x
    max_x <- bb$max_x
    min_y <- bb$min_y
    max_y <- bb$max_y
    if (!all(is.finite(c(min_x, max_x, min_y, max_y))) || !(min_x < max_x && min_y < max_y)) {
        stop("ROI bbox is invalid: ", roi_vertices_csv)
    }

    ds <- tryCatch(arrow::open_dataset(transcripts_path), error = function(e) e)
    if (inherits(ds, "error")) {
        stop("Failed to open Parquet transcripts: ", transcripts_path, " (", conditionMessage(ds), ")")
    }
    cols <- tryCatch(names(ds), error = function(e) character())
    if (!length(cols)) stop("Parquet transcripts have no columns: ", transcripts_path)
    cols_low <- tolower(cols)

    pick_col <- function(cands) {
        idx <- match(TRUE, cols_low %in% tolower(cands))
        if (is.na(idx)) NA_character_ else cols[[idx]]
    }
    gene_col <- pick_col(c("feature_name", "gene", "gene_name", "target", "target_name"))
    x_col <- pick_col(c("x_location", "x_um", "x"))
    y_col <- pick_col(c("y_location", "y_um", "y"))
    if (is.na(gene_col) || is.na(x_col) || is.na(y_col)) {
        stop(
            "Parquet transcripts missing required columns.\nExpected gene in: feature_name/gene/gene_name/target/target_name\nExpected x in: x_location/x_um/x\nExpected y in: y_location/y_um/y\nAvailable: ",
            paste(cols, collapse = ", "),
            "\nFile: ",
            transcripts_path
        )
    }

    reader <- tryCatch({
        if (inherits(ds, "Dataset")) {
            scan_builder <- ds$NewScan()
            scan_builder$Project(c(gene_col, x_col, y_col))
            scan_builder$UseThreads(TRUE)
            scan_builder$BatchSize(as.integer(batch_size))
            scan_builder$Finish()$ToRecordBatchReader()
        } else {
            arrow::as_record_batch_reader(ds)
        }
    }, error = function(e) e)
    if (inherits(reader, "error")) {
        stop("Failed to create Arrow reader for Parquet transcripts: ", transcripts_path, " (", conditionMessage(reader), ")")
    }

    mols_log_path <- file.path(dirname(wrapper_log_path), "make_molecules_bbox.log")
    .gridselect_wlog(
        wrapper_log_path,
        paste0("make_molecules_bbox parquet_arrow=", transcripts_path, " log=", mols_log_path),
        echo = TRUE
    )
    if (file.exists(mols_log_path)) file.remove(mols_log_path)

    con <- file(out_tsv, open = "wt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    writeLines("gene\tx_um\ty_um", con = con, sep = "\n", useBytes = TRUE)
    close(con)

    n_in <- 0L
    n_out <- 0L
    n_excl_prefix <- 0L

    drop_prefix <- NULL
    prefixes <- as.character(exclude_prefixes)
    prefixes <- prefixes[!is.na(prefixes) & nzchar(prefixes)]
    if (length(prefixes)) drop_prefix <- prefixes

    repeat {
        batch <- tryCatch(reader$read_next_batch(), error = function(e) e)
        if (inherits(batch, "error")) {
            stop("Parquet transcripts scan failed: ", transcripts_path, " (", conditionMessage(batch), ")")
        }
        if (is.null(batch)) break
        tbl <- tryCatch(as.data.frame(batch), error = function(e) NULL)
        if (is.null(tbl) || nrow(tbl) == 0) next
        n_in <- n_in + nrow(tbl)

        x <- suppressWarnings(as.numeric(tbl[[x_col]]))
        y <- suppressWarnings(as.numeric(tbl[[y_col]]))
        gene <- as.character(tbl[[gene_col]])
        ok <- is.finite(x) & is.finite(y) & !is.na(gene) & nzchar(gene)
        if (!any(ok)) next
        ok <- ok & x >= min_x & x <= max_x & y >= min_y & y <= max_y
        if (!any(ok)) next

        if (!is.null(drop_prefix)) {
            drop <- rep(FALSE, length(gene))
            for (pref in drop_prefix) drop <- drop | startsWith(gene, pref)
            n_excl_prefix <- n_excl_prefix + sum(ok & drop)
            ok <- ok & !drop
            if (!any(ok)) next
        }

        dt_out <- data.table::data.table(
            gene = gene[ok],
            x_um = x[ok],
            y_um = y[ok]
        )
        data.table::fwrite(
            dt_out,
            file = out_tsv,
            sep = "\t",
            append = TRUE,
            col.names = FALSE
        )
        n_out <- n_out + nrow(dt_out)
    }

    msg <- paste0(
        "[gridselect_make_molecules_bbox] wrote ",
        format(n_out, big.mark = ","),
        " / ",
        format(n_in, big.mark = ","),
        " transcripts rows to ",
        out_tsv,
        " (excluded_by_prefix=",
        format(n_excl_prefix, big.mark = ","),
        ")"
    )
    cat(msg, "\n", file = mols_log_path, append = TRUE, sep = "")
    invisible(out_tsv)
}

.gridselect_tool_map_out_dir_path <- function(host_path, out_dir, mode) {
    if (identical(mode, "bin")) {
        return(normalizePath(host_path, winslash = "/", mustWork = FALSE))
    }
    out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
    host_path <- normalizePath(host_path, winslash = "/", mustWork = FALSE)
    if (!startsWith(host_path, out_dir)) {
        stop("Internal path mapping error: file is not under out_dir.\nfile=", host_path, "\nout_dir=", out_dir)
    }
    rel <- substr(host_path, nchar(out_dir) + 1L, nchar(host_path))
    rel <- sub("^/", "", rel)
    if (!nzchar(rel)) return("/out")
    rel <- gsub("\\\\", "/", rel)
    paste0("/out/", rel)
}

.gridselect_run_idelta_gridselect_tool <- function(
    tool,
    out_dir,
    transcripts_path,
    roi_vertices_path,
    threads,
    knee_search_min_um,
    descending_delta_um,
    max_width_um = NULL,
    skip_molecules = FALSE,
    skip_tool = FALSE,
    wrapper_log_path,
    strict = TRUE
) {
    stopifnot(is.list(tool), is.character(out_dir), length(out_dir) == 1)
    out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
    transcripts_path <- normalizePath(transcripts_path, winslash = "/", mustWork = TRUE)
    roi_vertices_path <- normalizePath(roi_vertices_path, winslash = "/", mustWork = TRUE)
    staging_dir <- file.path(out_dir, "staging")
    logs_dir <- file.path(out_dir, "logs")
    tool_out_dir <- file.path(out_dir, "tool")
    dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(tool_out_dir, recursive = TRUE, showWarnings = FALSE)

    mols_bbox_path <- file.path(staging_dir, "molecules_bbox.tsv")
    knee_summary_path <- file.path(tool_out_dir, "knee_summary.tsv")
    tool_log <- file.path(logs_dir, "idelta_gridselect.log")

    is_parquet <- grepl("\\.(parquet|pq)$", transcripts_path, ignore.case = TRUE)
    is_csv <- grepl("\\.csv(\\.gz)?$", transcripts_path, ignore.case = TRUE)
    if (!is_parquet && !is_csv) {
        stop("Unsupported transcripts file type: ", transcripts_path, " (expected .csv(.gz) or .parquet/.pq)")
    }

    skip_molecules <- .gridselect_bool(skip_molecules, default = FALSE)
    skip_tool <- .gridselect_bool(skip_tool, default = FALSE)

    if (isTRUE(skip_molecules)) {
        if (!file.exists(mols_bbox_path)) stop("--skip_molecules was set but missing: ", mols_bbox_path)
    } else {
        if (is_csv) {
            .gridselect_make_molecules_bbox(
                roi_vertices_csv = roi_vertices_path,
                transcripts_path = transcripts_path,
                out_tsv = mols_bbox_path,
                wrapper_log_path = wrapper_log_path,
                python_bin = "python3"
            )
        } else if (is_parquet) {
            used_tool_subcmd <- FALSE
            if (identical(tool$mode, "bin") && is.character(tool$command) && nzchar(tool$command)) {
                mols_log_path <- file.path(logs_dir, "make_molecules_bbox.log")
                cmd <- tool$command
                args <- c(
                    tool$base_args,
                    "make-molecules-bbox",
                    "--roi-csv", roi_vertices_path,
                    "--transcripts", transcripts_path,
                    "--out-tsv", mols_bbox_path
                )
                status <- suppressWarnings(system2(cmd, args = shQuote(args), stdout = mols_log_path, stderr = mols_log_path))
                if (is.numeric(status) && status == 0) {
                    used_tool_subcmd <- TRUE
                }
            }
            if (!isTRUE(used_tool_subcmd)) {
                if (requireNamespace("arrow", quietly = TRUE)) {
                    .gridselect_make_molecules_bbox_parquet_arrow(
                        roi_vertices_csv = roi_vertices_path,
                        transcripts_path = transcripts_path,
                        out_tsv = mols_bbox_path,
                        wrapper_log_path = wrapper_log_path
                    )
                } else {
                    stop(
                        "Parquet transcripts require either:\n",
                        "  - idelta-gridselect with subcommand `make-molecules-bbox` (preferred; tool_bin/PATH), or\n",
                        "  - R package 'arrow' (host-side streaming prefilter).\n",
                        "Current tool mode=", tool$mode, " source=", tool$source,
                        "\ntranscripts=",
                        transcripts_path
                    )
                }
            }
        }
    }

    if (isTRUE(skip_tool)) {
        if (!file.exists(knee_summary_path)) stop("--skip_tool was set but missing: ", knee_summary_path)
        return(list(mode = tool$mode, source = tool$source, tool_log = tool_log, molecules_bbox = mols_bbox_path))
    }

    threads_num <- suppressWarnings(as.numeric(threads))[1]
    if (!is.finite(threads_num) || threads_num < 1) stop("`threads` must be an integer >= 1.")
    threads_int <- as.integer(round(threads_num))

    knee_search_min_um_num <- suppressWarnings(as.numeric(knee_search_min_um))[1]
    descending_delta_um_num <- suppressWarnings(as.numeric(descending_delta_um))[1]
    if (!is.finite(knee_search_min_um_num) || abs(knee_search_min_um_num - round(knee_search_min_um_num)) > 1e-9) {
        stop("`knee_search_min_um` must be an integer.")
    }
    if (!is.finite(descending_delta_um_num) || abs(descending_delta_um_num - round(descending_delta_um_num)) > 1e-9) {
        stop("`descending_delta_um` must be an integer.")
    }
    knee_search_min_um_int <- as.integer(round(knee_search_min_um_num))
    descending_delta_um_int <- as.integer(round(descending_delta_um_num))

    max_width_um_num <- suppressWarnings(as.numeric(max_width_um))[1]
    max_width_um_arg <- NULL
    if (is.finite(max_width_um_num)) {
        if (abs(max_width_um_num - round(max_width_um_num)) > 1e-9) stop("`max_width_um` must be an integer (or NULL).")
        max_width_um_arg <- as.integer(round(max_width_um_num))
        if (max_width_um_arg < 1L) stop("`max_width_um` must be >= 1 (or NULL).")
    }

    widths_prescreen_host <- NULL
    widths_anchor_host <- NULL
    if (!is.null(max_width_um_arg)) {
        w_pre <- .gridselect_default_widths_prescreen_v2()
        w_pre <- w_pre[w_pre <= max_width_um_arg]
        w_anc <- .gridselect_default_widths_anchor_v2()
        w_anc <- w_anc[w_anc <= max_width_um_arg]
        if (length(w_pre) < 3L) stop("max_width_um too small: prescreen widths would have < 3 entries.")
        if (length(w_anc) < 6L) stop("max_width_um too small: anchor widths would have < 6 entries (min_points=6).")

        widths_prescreen_host <- file.path(staging_dir, paste0("widths_prescreen_v2_max", max_width_um_arg, ".txt"))
        widths_anchor_host <- file.path(staging_dir, paste0("widths_anchor_v2_max", max_width_um_arg, ".txt"))
        .gridselect_write_widths_file(
            w_pre,
            widths_prescreen_host,
            header = c("Generated by geneSCOPE gridSelect()", paste0("base=widths_prescreen_v2 max_width_um=", max_width_um_arg))
        )
        .gridselect_write_widths_file(
            w_anc,
            widths_anchor_host,
            header = c("Generated by geneSCOPE gridSelect()", paste0("base=widths_anchor_v2 max_width_um=", max_width_um_arg))
        )
    }

    mode <- tool$mode
    mols_bbox_arg <- .gridselect_tool_map_out_dir_path(mols_bbox_path, out_dir = out_dir, mode = mode)
    roi_arg <- .gridselect_tool_map_out_dir_path(roi_vertices_path, out_dir = out_dir, mode = mode)
    tool_out_arg <- .gridselect_tool_map_out_dir_path(tool_out_dir, out_dir = out_dir, mode = mode)
    widths_pre_arg <- if (!is.null(widths_prescreen_host)) .gridselect_tool_map_out_dir_path(widths_prescreen_host, out_dir = out_dir, mode = mode) else NULL
    widths_anc_arg <- if (!is.null(widths_anchor_host)) .gridselect_tool_map_out_dir_path(widths_anchor_host, out_dir = out_dir, mode = mode) else NULL

    cmd <- tool$command
    args <- c(
        tool$base_args,
        "recommend-grid",
        "--molecules", mols_bbox_arg,
        "--roi", roi_arg,
        "--out", tool_out_arg,
        "--w0-um", "1",
        "--tile-um", "250",
        "--prescreen-fraction", "0.25",
        "--prescreen-reps", "25",
        "--seed", "1",
        "--knee-window-um", "20",
        "--knee-search-min-um", as.character(knee_search_min_um_int),
        "--descending-delta-um", as.character(descending_delta_um_int),
        "--threads", as.character(threads_int),
        "--b-min", "50",
        "--n-min", "100",
        "--min-points", "6",
        "--min-informative-frac-ge2", "0.6"
    )
    if (!is.null(widths_pre_arg) && !is.null(widths_anc_arg)) {
        args <- c(args, "--widths-prescreen", widths_pre_arg, "--widths-anchor", widths_anc_arg)
    }
    if (file.exists(tool_log)) file.remove(tool_log)
    status <- suppressWarnings(system2(cmd, args = shQuote(args), stdout = tool_log, stderr = tool_log))
    if (!is.numeric(status) || status != 0) {
        stop("idelta-gridselect failed (see logs): ", tool_log)
    }
    if (!file.exists(knee_summary_path) && isTRUE(strict)) {
        warning("Tool run succeeded but knee_summary.tsv missing under: ", tool_out_dir)
    }
    list(mode = mode, source = tool$source, tool_log = tool_log, molecules_bbox = mols_bbox_path)
}

compute_grid_idelta_qc <- function(
    tool_out_dir,
    grid_um_eval = NULL,
    epsilon = 1,
    q_bins = 5,
    ratio_high = 2,
    loss_min = 0.2
) {
    tool_out_dir <- normalizePath(tool_out_dir, winslash = "/", mustWork = TRUE)
    knee_summary_path <- file.path(tool_out_dir, "knee_summary.tsv")
    final_curves_path <- file.path(tool_out_dir, "final_curves.tsv.gz")

    if (!file.exists(knee_summary_path)) stop("Missing required tool output: ", knee_summary_path)
    if (!file.exists(final_curves_path)) {
        warning("Missing required tool output for QC: ", final_curves_path)
    }

    to_num <- function(x) suppressWarnings(as.numeric(x))
    to_int <- function(x) suppressWarnings(as.integer(round(as.numeric(x))))

    epsilon <- to_num(epsilon)[1]
    if (!is.finite(epsilon) || epsilon < 0) stop("`epsilon` must be a non-negative number.")
    q_bins <- to_int(q_bins)[1]
    if (!is.finite(q_bins) || q_bins < 2) stop("`q_bins` must be an integer >= 2.")
    ratio_high <- to_num(ratio_high)[1]
    if (!is.finite(ratio_high) || ratio_high <= 1) stop("`ratio_high` must be > 1.")
    loss_min <- to_num(loss_min)[1]
    if (!is.finite(loss_min) || loss_min < 0 || loss_min > 1) stop("`loss_min` must be in [0, 1].")

    knee_summary <- utils::read.delim(knee_summary_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    if (!("gene" %in% colnames(knee_summary))) stop("knee_summary.tsv missing required column 'gene': ", knee_summary_path)
    summary_row <- knee_summary[knee_summary$gene == "__SUMMARY__", , drop = FALSE]
    if (nrow(summary_row) != 1) stop("Expected exactly one __SUMMARY__ row in: ", knee_summary_path)

    recommended_grid_um <- to_num(if ("recommended_grid_um" %in% colnames(summary_row)) summary_row$recommended_grid_um else NA)[1]
    if (!is.finite(recommended_grid_um)) recommended_grid_um <- NA_real_

    grid_um_eval_num <- to_num(grid_um_eval)[1]
    if (!is.finite(grid_um_eval_num)) {
        if (is.finite(recommended_grid_um)) {
            grid_um_eval_num <- recommended_grid_um
        } else {
            grid_um_eval_num <- NA_real_
        }
    }
    if (!is.finite(grid_um_eval_num) || grid_um_eval_num <= 0) {
        stop("Could not resolve grid_um_eval (provide `grid_um_eval=` or ensure knee_summary.tsv has recommended_grid_um).")
    }

    gene_rows <- knee_summary[knee_summary$gene != "__SUMMARY__", , drop = FALSE]
    if (nrow(gene_rows) == 0) stop("knee_summary.tsv had no per-gene rows: ", knee_summary_path)

    gene_status_missing <- !("gene_status" %in% colnames(gene_rows))
    if (gene_status_missing) warning("knee_summary.tsv missing gene_status; treating all non-__SUMMARY__ rows as informative for QC.")
    informative_rows <- if (!gene_status_missing) gene_rows[gene_rows$gene_status == "informative", , drop = FALSE] else gene_rows
    informative_rows <- informative_rows[nzchar(informative_rows$gene), , drop = FALSE]

    n_genes_total <- nrow(gene_rows)
    n_genes_informative <- nrow(informative_rows)

    if (n_genes_informative == 0) {
        qc_summary <- data.frame(
            grid_um_eval = grid_um_eval_num,
            recommended_grid_um = recommended_grid_um,
            epsilon = epsilon,
            q_bins = q_bins,
            ratio_high = ratio_high,
            loss_min = loss_min,
            n_genes_total = n_genes_total,
            n_genes_informative = 0,
            idelta_mean = NA_real_,
            idelta_sd = NA_real_,
            idelta_cv = NA_real_,
            idelta_median = NA_real_,
            idelta_q1 = NA_real_,
            idelta_q3 = NA_real_,
            idelta_iqr = NA_real_,
            idelta_robust_cv = NA_real_,
            frac_near_zero = NA_real_,
            n_flag_near_zero = 0L,
            n_flag_scale_mismatch = 0L,
            n_flag_loss = 0L,
            n_unsuitable_any = 0L,
            frac_flag_near_zero = NA_real_,
            frac_flag_scale_mismatch = NA_real_,
            frac_flag_loss = NA_real_,
            frac_unsuitable_any = NA_real_,
            stringsAsFactors = FALSE
        )
        gene_bins <- data.frame(
            gene = character(0),
            grid_um_eval = numeric(0),
            width_um_selected = numeric(0),
            idelta_at_grid = numeric(0),
            tier_q5 = character(0),
            knee_um = numeric(0),
            gene_status = character(0),
            stringsAsFactors = FALSE
        )
        return(list(qc_summary = qc_summary, gene_bins = gene_bins))
    }

    informative_rows$knee_um <- to_num(if ("knee_um" %in% colnames(informative_rows)) informative_rows$knee_um else NA)
    informative_genes <- sort(unique(as.character(informative_rows$gene)))

    final_curves <- if (file.exists(final_curves_path)) {
        utils::read.delim(gzfile(final_curves_path), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    } else {
        data.frame()
    }
    if (!("gene" %in% colnames(final_curves)) || !("width_um" %in% colnames(final_curves))) {
        stop("final_curves.tsv.gz missing required columns gene/width_um: ", final_curves_path)
    }

    # Restrict to informative genes for QC (deterministic ordering).
    fc <- final_curves[final_curves$gene %in% informative_genes, , drop = FALSE]
    fc$width_um <- to_num(fc$width_um)
    fc <- fc[is.finite(fc$width_um), , drop = FALSE]
    fc$gene <- as.character(fc$gene)

    pick_row_per_gene <- function(df, target_w) {
        # Deterministic: within each gene, choose the row with minimal |width-target|,
        # breaking ties by smaller width.
        df <- df[order(abs(df$width_um - target_w), df$width_um), , drop = FALSE]
        df[1, , drop = FALSE]
    }

    # Selected per-gene row at the evaluation grid.
    split_fc <- split(fc, fc$gene)
    picked_list <- lapply(informative_genes, function(g) {
        dfg <- split_fc[[g]]
        if (is.null(dfg) || nrow(dfg) == 0) return(NULL)
        pick_row_per_gene(dfg, grid_um_eval_num)
    })
    picked <- do.call(rbind, picked_list)
    if (is.null(picked) || nrow(picked) == 0) stop("No final_curves rows matched informative genes.")

    picked$width_um_selected <- to_num(picked$width_um)
    picked$idelta_smooth_at_grid <- to_num(if ("idelta_smooth" %in% colnames(picked)) picked$idelta_smooth else NA)
    picked$idelta_raw_at_grid <- to_num(if ("idelta_raw" %in% colnames(picked)) picked$idelta_raw else NA)
    picked$idelta_at_grid <- ifelse(is.finite(picked$idelta_smooth_at_grid), picked$idelta_smooth_at_grid, picked$idelta_raw_at_grid)

    # Join knees/status onto picked table.
    gene_meta <- informative_rows[, intersect(c("gene", "knee_um", "gene_status"), colnames(informative_rows)), drop = FALSE]
    gene_meta$gene <- as.character(gene_meta$gene)
    picked$gene <- as.character(picked$gene)
    picked <- merge(picked, gene_meta, by = "gene", all.x = TRUE, sort = FALSE)

    # Max Idelta per gene across widths (for loss flag).
    if ("idelta_smooth" %in% colnames(fc)) {
        y_for_max <- to_num(fc$idelta_smooth)
    } else if ("idelta_raw" %in% colnames(fc)) {
        y_for_max <- to_num(fc$idelta_raw)
    } else {
        y_for_max <- rep(NA_real_, nrow(fc))
    }
    fc_max_tbl <- data.frame(gene = fc$gene, y = y_for_max, status = if ("status" %in% colnames(fc)) fc$status else NA_character_)
    fc_max_tbl <- fc_max_tbl[is.finite(fc_max_tbl$y), , drop = FALSE]
    if (nrow(fc_max_tbl) > 0) {
        max_by_gene <- tapply(fc_max_tbl$y, fc_max_tbl$gene, max, na.rm = TRUE)
        max_df <- data.frame(gene = names(max_by_gene), idelta_max_gene = as.numeric(max_by_gene), stringsAsFactors = FALSE)
        picked <- merge(picked, max_df, by = "gene", all.x = TRUE, sort = FALSE)
    } else {
        picked$idelta_max_gene <- NA_real_
    }
    picked$loss_ratio <- ifelse(
        is.finite(picked$idelta_at_grid) & is.finite(picked$idelta_max_gene) & picked$idelta_max_gene > 0,
        picked$idelta_at_grid / picked$idelta_max_gene,
        NA_real_
    )

    # Flags.
    picked$flag_near_zero <- is.finite(picked$idelta_at_grid) & (picked$idelta_at_grid < epsilon)

    # Degenerate/0-1 occupancy regime diagnostics (optional, transparent).
    have_deg_cols <- all(c("sum_n_n1", "n_bins_ge2", "N_total") %in% colnames(picked))
    if (have_deg_cols) {
        sum_n_n1 <- to_num(picked$sum_n_n1)
        n_bins_ge2 <- to_num(picked$n_bins_ge2)
        n_total <- to_num(picked$N_total)
        picked$flag_deg_01 <- (sum_n_n1 == 0 & is.finite(n_total) & n_total >= 100) | (is.finite(n_bins_ge2) & n_bins_ge2 < 50)
    } else {
        picked$flag_deg_01 <- NA
    }

    # Scale mismatch reference: use the evaluation grid (grid_um_eval), not the tool-recommended grid.
    # This keeps QC fully tied to the user-chosen grid size (or the default resolved grid_um_eval).
    g <- grid_um_eval_num
    picked$flag_scale_mismatch <- ifelse(
        is.finite(picked$knee_um) & is.finite(g) & picked$knee_um > 0 & g > 0,
        abs(log(picked$knee_um / g)) > log(ratio_high),
        NA
    )
    picked$flag_loss <- is.finite(picked$loss_ratio) & (picked$loss_ratio < loss_min)
    picked$unsuitable_any <- (picked$flag_near_zero %in% TRUE) | (picked$flag_scale_mismatch %in% TRUE) | (picked$flag_loss %in% TRUE)

    # Tiering (quantile bins) on idelta_at_grid (informative genes).
    tier_levels <- paste0("Q", seq_len(q_bins))
    tier_for <- rep(NA_character_, nrow(picked))
    ok <- is.finite(picked$idelta_at_grid)
    if (sum(ok) >= q_bins) {
        x <- picked$idelta_at_grid[ok]
        qs <- as.numeric(stats::quantile(x, probs = seq(0, 1, length.out = q_bins + 1), na.rm = TRUE, type = 7))
        if (length(unique(qs)) == length(qs)) {
            tier_for[ok] <- as.character(cut(x, breaks = qs, include.lowest = TRUE, labels = tier_levels))
        } else {
            # Fallback for heavy ties: rank-based bins (deterministic after sorting by gene).
            r <- rank(x, ties.method = "average", na.last = "keep")
            idx <- pmin(q_bins, pmax(1L, floor((r - 1) / length(r) * q_bins) + 1L))
            tier_for[ok] <- tier_levels[idx]
        }
    }
    picked$tier_q5 <- factor(tier_for, levels = tier_levels, ordered = TRUE)

    # Final deterministic output columns.
    picked$grid_um_eval <- grid_um_eval_num
    picked <- picked[order(picked$gene), , drop = FALSE]

    gene_bins <- data.frame(
        gene = picked$gene,
        grid_um_eval = picked$grid_um_eval,
        width_um_selected = picked$width_um_selected,
        idelta_at_grid = picked$idelta_at_grid,
        tier_q5 = as.character(picked$tier_q5),
        knee_um = to_num(picked$knee_um),
        gene_status = if ("gene_status" %in% colnames(picked)) as.character(picked$gene_status) else "informative",
        flag_near_zero = as.logical(picked$flag_near_zero),
        flag_deg_01 = picked$flag_deg_01,
        flag_scale_mismatch = picked$flag_scale_mismatch,
        flag_loss = as.logical(picked$flag_loss),
        unsuitable_any = as.logical(picked$unsuitable_any),
        loss_ratio = picked$loss_ratio,
        idelta_max_gene = picked$idelta_max_gene,
        status_at_grid = if ("status" %in% colnames(picked)) as.character(picked$status) else NA_character_,
        N_total_at_grid = if ("N_total" %in% colnames(picked)) to_num(picked$N_total) else NA_real_,
        Q_bins_at_grid = if ("Q_bins" %in% colnames(picked)) to_num(picked$Q_bins) else NA_real_,
        n_bins_ge2_at_grid = if ("n_bins_ge2" %in% colnames(picked)) to_num(picked$n_bins_ge2) else NA_real_,
        frac_bins_ge2_at_grid = if ("frac_bins_ge2" %in% colnames(picked)) to_num(picked$frac_bins_ge2) else NA_real_,
        max_bin_count_at_grid = if ("max_bin_count" %in% colnames(picked)) to_num(picked$max_bin_count) else NA_real_,
        sum_n_n1_at_grid = if ("sum_n_n1" %in% colnames(picked)) to_num(picked$sum_n_n1) else NA_real_,
        idelta_raw_at_grid = picked$idelta_raw_at_grid,
        idelta_smooth_at_grid = picked$idelta_smooth_at_grid,
        stringsAsFactors = FALSE
    )

    v <- gene_bins$idelta_at_grid
    n_inf <- nrow(gene_bins)
    v_ok <- v[is.finite(v)]
    q1 <- if (length(v_ok) > 0) as.numeric(stats::quantile(v_ok, 0.25, na.rm = TRUE, type = 7)) else NA_real_
    q3 <- if (length(v_ok) > 0) as.numeric(stats::quantile(v_ok, 0.75, na.rm = TRUE, type = 7)) else NA_real_
    med <- if (length(v_ok) > 0) stats::median(v_ok, na.rm = TRUE) else NA_real_
    iqr <- if (is.finite(q1) && is.finite(q3)) (q3 - q1) else NA_real_

    id_mean <- if (length(v_ok) > 0) mean(v_ok, na.rm = TRUE) else NA_real_
    id_sd <- if (length(v_ok) > 1) stats::sd(v_ok, na.rm = TRUE) else NA_real_
    id_cv <- if (is.finite(id_mean) && id_mean != 0) id_sd / id_mean else NA_real_
    id_robust_cv <- if (is.finite(med) && med != 0 && is.finite(iqr)) iqr / med else NA_real_

    n_flag_near_zero <- sum(gene_bins$flag_near_zero %in% TRUE, na.rm = TRUE)
    n_flag_scale_mismatch <- sum(gene_bins$flag_scale_mismatch %in% TRUE, na.rm = TRUE)
    n_flag_loss <- sum(gene_bins$flag_loss %in% TRUE, na.rm = TRUE)
    n_unsuitable_any <- sum(gene_bins$unsuitable_any %in% TRUE, na.rm = TRUE)

    tier_counts <- table(factor(gene_bins$tier_q5, levels = tier_levels))
    tier_counts <- as.integer(tier_counts)
    names(tier_counts) <- tier_levels

    qc_summary <- data.frame(
        grid_um_eval = grid_um_eval_num,
        recommended_grid_um = recommended_grid_um,
        epsilon = epsilon,
        q_bins = q_bins,
        ratio_high = ratio_high,
        loss_min = loss_min,
        n_genes_total = n_genes_total,
        n_genes_informative = n_genes_informative,
        n_idelta_finite = length(v_ok),
        idelta_mean = id_mean,
        idelta_sd = id_sd,
        idelta_cv = id_cv,
        idelta_median = med,
        idelta_q1 = q1,
        idelta_q3 = q3,
        idelta_iqr = iqr,
        idelta_robust_cv = id_robust_cv,
        frac_near_zero = n_flag_near_zero / n_inf,
        n_flag_near_zero = n_flag_near_zero,
        n_flag_scale_mismatch = n_flag_scale_mismatch,
        n_flag_loss = n_flag_loss,
        n_unsuitable_any = n_unsuitable_any,
        frac_flag_near_zero = n_flag_near_zero / n_inf,
        frac_flag_scale_mismatch = n_flag_scale_mismatch / n_inf,
        frac_flag_loss = n_flag_loss / n_inf,
        frac_unsuitable_any = n_unsuitable_any / n_inf,
        stringsAsFactors = FALSE
    )
    # Add remaining tier columns deterministically.
    for (i in seq_len(q_bins)) {
        nm <- paste0("tier_Q", i)
        qc_summary[[nm]] <- tier_counts[[paste0("Q", i)]]
    }

    list(qc_summary = qc_summary, gene_bins = gene_bins)
}

.gridselect_consume_tool_outputs <- function(
    tool_out_dir,
    plots_out_dir,
    log_out_path = NULL,
    best_grid_style = TRUE,
    strict = TRUE,
    smooth_method = NULL
) {
    stopifnot(is.character(tool_out_dir), length(tool_out_dir) == 1, nzchar(tool_out_dir))
    stopifnot(is.character(plots_out_dir), length(plots_out_dir) == 1, nzchar(plots_out_dir))
    stopifnot(is.null(log_out_path) || (is.character(log_out_path) && length(log_out_path) == 1))

    tool_out_dir <- normalizePath(tool_out_dir, mustWork = TRUE)
    plots_out_dir <- normalizePath(plots_out_dir, mustWork = FALSE)
    dir.create(plots_out_dir, recursive = TRUE, showWarnings = FALSE)

    knee_summary_path <- file.path(tool_out_dir, "knee_summary.tsv")
    avg_curve_path <- file.path(tool_out_dir, "avg_curve.tsv")
    final_curves_path <- file.path(tool_out_dir, "final_curves.tsv.gz")

    if (!file.exists(knee_summary_path)) stop("Missing required tool output: ", knee_summary_path)

    to_num <- function(x) suppressWarnings(as.numeric(x))
    get_col <- function(df, name) if (name %in% colnames(df)) df[[name]] else rep(NA, nrow(df))

    knee_summary <- utils::read.delim(knee_summary_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    if (!("gene" %in% colnames(knee_summary))) stop("knee_summary.tsv missing required column 'gene': ", knee_summary_path)

    summary_row <- knee_summary[knee_summary$gene == "__SUMMARY__", , drop = FALSE]
    if (nrow(summary_row) != 1) stop("Expected exactly one summary row (gene='__SUMMARY__') in: ", knee_summary_path)

    recommended_grid_um <- to_num(get_col(summary_row, "recommended_grid_um"))[1]
    overall_knee_um <- to_num(get_col(summary_row, "overall_knee_um"))[1]
    prescreen_global_knee_um <- to_num(get_col(summary_row, "prescreen_global_knee_um"))[1]

    n_genes_total <- to_num(get_col(summary_row, "n_genes_total"))[1]
    n_genes_informative <- to_num(get_col(summary_row, "n_genes_informative"))[1]
    n_genes_excluded_degenerate <- to_num(get_col(summary_row, "n_genes_excluded_degenerate"))[1]

    recommendation_status <- as.character(get_col(summary_row, "recommendation_status"))[1]
    recommendation_note <- as.character(get_col(summary_row, "recommendation_note"))[1]

    b_min <- to_num(get_col(summary_row, "b_min"))[1]
    n_min <- to_num(get_col(summary_row, "n_min"))[1]
    min_points <- to_num(get_col(summary_row, "min_points"))[1]
    min_informative_frac_ge2 <- to_num(get_col(summary_row, "min_informative_frac_ge2"))[1]

    knee_mean <- to_num(get_col(summary_row, "knee_mean"))[1]
    knee_sd <- to_num(get_col(summary_row, "knee_sd"))[1]
    knee_min <- to_num(get_col(summary_row, "knee_min"))[1]
    knee_max <- to_num(get_col(summary_row, "knee_max"))[1]

    knee_of_mean_curve_um <- to_num(get_col(summary_row, "knee_of_mean_curve_um"))[1]
    overall_knee_median_um <- to_num(get_col(summary_row, "overall_knee_median_um"))[1]
    overall_knee_trimmed_mean_um <- to_num(get_col(summary_row, "overall_knee_trimmed_mean_um"))[1]
    overall_knee_sd_um <- to_num(get_col(summary_row, "overall_knee_sd_um"))[1]
    overall_knee_iqr_um <- to_num(get_col(summary_row, "overall_knee_iqr_um"))[1]
    knee_search_min_um <- to_num(get_col(summary_row, "knee_search_min_um"))[1]
    descending_delta_um <- to_num(get_col(summary_row, "descending_delta_um"))[1]
    widths_anchor_version <- as.character(get_col(summary_row, "widths_anchor_version"))[1]

    # Backward-compatible fallbacks:
    # - Older idelta-gridselect versions wrote `overall_knee_um` as the knee of the mean curve.
    # - Newer versions split this into `overall_knee_um` (mean of gene knees) and `knee_of_mean_curve_um`.
    gene_rows_all <- knee_summary[knee_summary$gene != "__SUMMARY__", , drop = FALSE]
    informative_rows <- if ("gene_status" %in% colnames(gene_rows_all)) {
        gene_rows_all[gene_rows_all$gene_status == "informative", , drop = FALSE]
    } else {
        gene_rows_all
    }
    knee_vals_informative <- to_num(get_col(informative_rows, "knee_um"))
    knee_vals_informative <- knee_vals_informative[is.finite(knee_vals_informative)]
    knee_sorted <- sort(knee_vals_informative)

    upper_median <- function(x_sorted) {
        if (length(x_sorted) == 0) return(NA_real_)
        x_sorted[floor(length(x_sorted) / 2) + 1]
    }
    trimmed_mean_u32 <- function(x_sorted, trim = 0.10) {
        if (length(x_sorted) == 0) return(NA_real_)
        n <- length(x_sorted)
        k <- floor(max(0, min(0.49, trim)) * n)
        if (2 * k >= n) k <- floor(n / 2)
        slice <- x_sorted[(k + 1):(n - k)]
        if (length(slice) == 0) return(NA_real_)
        mean(slice)
    }
    iqr_u32 <- function(x_sorted) {
        if (length(x_sorted) == 0) return(NA_real_)
        n <- length(x_sorted)
        q1 <- x_sorted[floor(n / 4) + 1]
        q3 <- x_sorted[floor((3 * n) / 4) + 1]
        q3 - q1
    }

    if (!is.finite(overall_knee_median_um)) overall_knee_median_um <- upper_median(knee_sorted)
    if (!is.finite(overall_knee_trimmed_mean_um)) overall_knee_trimmed_mean_um <- trimmed_mean_u32(knee_sorted, 0.10)
    if (!is.finite(overall_knee_sd_um) && length(knee_sorted) > 0) overall_knee_sd_um <- if (length(knee_sorted) >= 2) stats::sd(knee_sorted) else 0
    if (!is.finite(overall_knee_iqr_um)) overall_knee_iqr_um <- iqr_u32(knee_sorted)

    has_knee_of_mean_curve <- "knee_of_mean_curve_um" %in% colnames(knee_summary)
    if (!has_knee_of_mean_curve) {
        # Legacy schema: interpret old `overall_knee_um` as knee-of-mean, and use the mean/median of gene-wise knees instead.
        if (is.finite(overall_knee_um) && !is.finite(knee_of_mean_curve_um)) knee_of_mean_curve_um <- overall_knee_um
        if (is.finite(knee_mean)) {
            overall_knee_um <- knee_mean
        } else if (length(knee_sorted) > 0) {
            overall_knee_um <- mean(knee_sorted)
        }
        if (is.finite(overall_knee_median_um)) {
            recommended_grid_um <- overall_knee_median_um
            recommendation_status <- "wrapper_policy_median_gene_knees_legacy_tool"
            recommendation_note <- "legacy knee_summary.tsv schema (missing knee_of_mean_curve_um); wrapper set recommended_grid_um=median of informative gene knees; rerun idelta-gridselect>=v0.1.2 for fully auditable summary columns"
        }
        if (isTRUE(strict)) {
            warning("Legacy idelta-gridselect knee_summary.tsv schema detected; rerun with idelta-gridselect>=v0.1.2 to refresh outputs (denser widths + richer summary).")
        }
    } else {
        if (!is.finite(overall_knee_um) && is.finite(knee_mean)) overall_knee_um <- knee_mean
    }

    smooth_method_norm <- NULL
    if (!is.null(smooth_method) && is.character(smooth_method) && length(smooth_method) >= 1) {
        smooth_method_norm <- tolower(trimws(smooth_method[[1]]))
        if (!nzchar(smooth_method_norm)) smooth_method_norm <- NULL
        if (!is.null(smooth_method_norm) && smooth_method_norm %in% c("null", "none", "tool", "current")) {
            smooth_method_norm <- NULL
        }
    }
    smooth_method_label <- if (is.null(smooth_method_norm)) {
        "tool_default (idelta_smooth from idelta-gridselect)"
    } else if (identical(smooth_method_norm, "loess")) {
        "loess (R stats::loess on idelta_raw; plots only)"
    } else {
        stop("Unsupported smooth_method: ", smooth_method, " (use NULL or 'loess')")
    }

    fmt2 <- function(x) ifelse(is.finite(x), sprintf("%.2f", x), "NA")
    fmt_int <- function(x) ifelse(is.finite(x), as.character(as.integer(round(x))), "NA")

    degenerate_thresholds_line <- paste0(
        "degenerate_thresholds: ",
        "minimum n_bins_ge2 (bins with >=2 molecules) per gene/width (b_min)=", fmt_int(b_min),
        "; degenerate (0/1 occupancy) if sum_n_n1==0 and N_total>=n_min (n_min)=", fmt_int(n_min),
        "; minimum valid widths per gene for knee estimation (min_points)=", fmt_int(min_points),
        "; at recommended grid, minimum fraction of informative genes with n_bins_ge2>=b_min (min_informative_frac_ge2)=", fmt2(min_informative_frac_ge2)
    )

    lines <- c(
        sprintf("overall_knee_um (mean of gene knees): %s", fmt2(overall_knee_um)),
        sprintf("overall_knee_median_um: %s", fmt2(overall_knee_median_um)),
        sprintf("overall_knee_trimmed_mean_um (trim=0.10): %s", fmt2(overall_knee_trimmed_mean_um)),
        sprintf("overall_knee_sd_um: %s", fmt2(overall_knee_sd_um)),
        sprintf("overall_knee_iqr_um: %s", fmt2(overall_knee_iqr_um)),
        sprintf("n_genes_total: %s", fmt_int(n_genes_total)),
        sprintf("n_genes_informative: %s", fmt_int(n_genes_informative)),
        degenerate_thresholds_line
    )

    cat(paste0(lines, collapse = "\n"), "\n", sep = "")
    if (!is.null(log_out_path) && nzchar(log_out_path)) {
        dir.create(dirname(log_out_path), recursive = TRUE, showWarnings = FALSE)
        writeLines(lines, con = log_out_path, useBytes = TRUE)
    }

    if (isTRUE(strict)) {
        if (is.finite(n_genes_total) && n_genes_total > 0 && is.finite(n_genes_excluded_degenerate)) {
            frac_excl <- n_genes_excluded_degenerate / n_genes_total
            if (is.finite(frac_excl) && frac_excl > 0.5) {
                warning(
                    sprintf(
                        "High degenerate exclusion fraction: %.1f%% (%d / %d). Small-width Idelta may be in a 0/1-occupancy regime.",
                        100 * frac_excl,
                        as.integer(n_genes_excluded_degenerate),
                        as.integer(n_genes_total)
                    )
                )
            }
        }
        if (is.finite(n_genes_informative) && n_genes_informative == 0) {
            warning("No informative genes after degenerate filtering; recommendation may be a deterministic fallback.")
        }
        if (nzchar(recommendation_status) && grepl("fallback", recommendation_status, fixed = TRUE)) {
            warning("Recommendation used a fallback: ", recommendation_status)
        }
    }

    if (isTRUE(best_grid_style)) {
        have_plot_pkgs <- requireNamespace("ggplot2", quietly = TRUE) &&
            requireNamespace("dplyr", quietly = TRUE) &&
            requireNamespace("tibble", quietly = TRUE)
        if (!isTRUE(have_plot_pkgs)) {
            if (isTRUE(strict)) warning("Missing ggplot2/dplyr/tibble; skipping best_grid-style plots.")
        } else {
            red_color <- "#E41A1C"
            blue_color <- "#377EB8"

            # ---- Read diagnostics tables ----
            final_curves <- NULL
            if (file.exists(final_curves_path)) {
                final_curves <- utils::read.delim(
                    gzfile(final_curves_path),
                    sep = "\t",
                    stringsAsFactors = FALSE,
                    check.names = FALSE
                )
            } else if (isTRUE(strict)) {
                warning("final_curves.tsv.gz not found; annotated plots may be incomplete: ", final_curves_path)
            }

            avg_curve <- NULL
            if (file.exists(avg_curve_path)) {
                avg_curve <- utils::read.delim(
                    avg_curve_path,
                    sep = "\t",
                    stringsAsFactors = FALSE,
                    check.names = FALSE
                )
                if (!all(c("width_um", "idelta_mean") %in% colnames(avg_curve))) {
                    if (isTRUE(strict)) warning("avg_curve.tsv missing expected columns; will recompute mean curve from final_curves: ", avg_curve_path)
                    avg_curve <- NULL
                }
            }

            # ---- Informative gene knees ----
            gene_rows <- knee_summary[knee_summary$gene != "__SUMMARY__", , drop = FALSE]
            if (!("knee_um" %in% colnames(gene_rows))) {
                if (isTRUE(strict)) warning("knee_summary.tsv missing knee_um; skipping plots.")
            } else {
                if ("gene_status" %in% colnames(gene_rows)) {
                    informative <- gene_rows[gene_rows$gene_status == "informative", , drop = FALSE]
                } else {
                    informative <- gene_rows
                }
                informative$knee_um <- to_num(informative$knee_um)
                informative <- informative[is.finite(informative$knee_um) & nzchar(informative$gene), , drop = FALSE]

                knee_vals <- informative$knee_um

                # Optional plot-only smoothing override: compute loess-smoothed curves from idelta_raw.
                final_plot <- NULL
                if (identical(smooth_method_norm, "loess")) {
                    if (!is.null(final_curves) &&
                        all(c("gene", "width_um", "idelta_raw", "status") %in% colnames(final_curves)) &&
                        nrow(informative) > 0) {
                        loess_span <- 0.1
                        final_plot <- final_curves |>
                            dplyr::filter(
                                gene %in% informative$gene,
                                status == "ok",
                                is.finite(width_um),
                                is.finite(idelta_raw)
                            ) |>
                            dplyr::group_by(gene) |>
                            dplyr::group_modify(~ {
                                df <- .x |>
                                    dplyr::arrange(width_um) |>
                                    dplyr::select(gene, width_um, idelta_raw)
                                x <- df$width_um
                                y <- df$idelta_raw
                                if (length(unique(x)) < 6) {
                                    df$idelta_smooth <- pmax(y, 0)
                                    return(df)
                                }
                                fit <- try(stats::loess(y ~ x, span = loess_span), silent = TRUE)
                                if (inherits(fit, "try-error")) {
                                    df$idelta_smooth <- pmax(stats::approx(x, y, xout = x, rule = 2)$y, 0)
                                    return(df)
                                }
                                preds <- suppressWarnings(stats::predict(fit, newdata = data.frame(x = x)))
                                df$idelta_smooth <- pmax(as.numeric(preds), 0)
                                df
                            }) |>
                            dplyr::ungroup()

                        if (nrow(final_plot) > 0) {
                            avg_curve <- final_plot |>
                                dplyr::group_by(width_um) |>
                                dplyr::summarise(idelta_mean = mean(idelta_smooth, na.rm = TRUE), .groups = "drop") |>
                                dplyr::arrange(width_um)
                        } else {
                            final_plot <- NULL
                            if (isTRUE(strict)) warning("smooth_method='loess' requested but no finite idelta_raw points were available; falling back to tool idelta_smooth.")
                        }
                    } else if (isTRUE(strict)) {
                        warning("smooth_method='loess' requested but final_curves.tsv.gz missing required columns; falling back to tool idelta_smooth.")
                    }
                }

                # ---- Mean curve fallback ----
                if (is.null(avg_curve) && !is.null(final_curves) &&
                    all(c("width_um", "idelta_smooth", "status") %in% colnames(final_curves))) {
                    avg_curve <- final_curves |>
                        dplyr::filter(status == "ok", is.finite(idelta_smooth)) |>
                        dplyr::group_by(width_um) |>
                        dplyr::summarise(idelta_mean = mean(idelta_smooth, na.rm = TRUE), .groups = "drop") |>
                        dplyr::arrange(width_um)
                }

                # Knee range ribbon (min-max across informative knees).
                x_rng <- if (length(knee_vals) > 0) {
                    c(min(knee_vals), max(knee_vals))
                } else if (!is.null(avg_curve) && nrow(avg_curve) > 0 && any(is.finite(avg_curve$width_um))) {
                    range(avg_curve$width_um, na.rm = TRUE)
                } else if (!is.null(final_curves) && any(is.finite(final_curves$width_um))) {
                    range(final_curves$width_um, na.rm = TRUE)
                } else {
                    c(0, 0)
                }

                y_max <- 100
                if (!is.null(avg_curve) && any(is.finite(avg_curve$idelta_mean))) {
                    y_max <- min(100, max(avg_curve$idelta_mean, na.rm = TRUE))
                } else if (!is.null(final_curves) && any(is.finite(final_curves$idelta_smooth))) {
                    y_max <- min(100, max(final_curves$idelta_smooth, na.rm = TRUE))
                }
                knee_range <- data.frame(x = x_rng, ymin = 0, ymax = y_max)

                # ---- 1) Mean curve with knee ribbon ----
                p_main <- if (!is.null(avg_curve) && nrow(avg_curve) > 0) {
                    ggplot2::ggplot(avg_curve, ggplot2::aes(x = width_um, y = idelta_mean)) +
                        ggplot2::geom_ribbon(
                            data = knee_range,
                            ggplot2::aes(x = x, ymin = ymin - 5, ymax = ymax + 5),
                            fill = blue_color,
                            alpha = 0.5,
                            linewidth = 0,
                            color = blue_color,
                            inherit.aes = FALSE
                        ) +
                        ggplot2::geom_line(color = "black", linewidth = 1) +
                        {
                            if (is.finite(overall_knee_um)) {
                                ggplot2::geom_vline(
                                    xintercept = overall_knee_um,
                                    color = red_color,
                                    linewidth = 1,
                                    linetype = "solid"
                                )
                            }
                        } +
                        ggplot2::scale_x_reverse(
                            breaks = seq(
                                min(avg_curve$width_um, na.rm = TRUE),
                                max(avg_curve$width_um, na.rm = TRUE),
                                by = 5
                            )
                        ) +
                        ggplot2::labs(
                            title = "scope - Knee Point Analysis",
                            x = "Grid size (umm)",
                            y = "Mean Curve Idelta"
                        ) +
                        ggplot2::theme_minimal(base_size = 8) +
                        ggplot2::theme(
                            panel.border = ggplot2::element_rect(color = "black", fill = NA),
                            panel.background = ggplot2::element_rect(fill = "#C0C0C0", color = NA),
                            plot.background = ggplot2::element_rect(fill = "#ffffff", color = NA),
                            axis.title = ggplot2::element_text(size = 9),
                            axis.text.y = ggplot2::element_text(size = 8),
                            axis.text.x = ggplot2::element_text(size = 8, angle = 45, hjust = 1),
                            axis.ticks = ggplot2::element_line(color = "black"),
                            plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
                            legend.position = "none",
                            plot.margin = grid::unit(c(1, 1, 1, 1), "cm")
                        ) +
                        ggplot2::scale_y_continuous(
                            limits = c(0, min(100, max(avg_curve$idelta_mean, na.rm = TRUE))),
                            oob = scales::squish,
                            breaks = scales::breaks_width(5)
                        )
                } else {
                    ggplot2::ggplot() +
                        ggplot2::theme_void() +
                        ggplot2::annotate("text", x = 0, y = 0, label = "Mean curve unavailable")
                }

                ggplot2::ggsave(
                    filename = file.path(plots_out_dir, "scope_iDelta_mean_curve_knee_stable.png"),
                    plot = p_main,
                    width = 6,
                    height = 5,
                    dpi = 600
                )

                # ---- 2) All-genes annotated plot ----
                p_all <- NULL
                if (is.null(final_plot) && !is.null(final_curves) &&
                    all(c("gene", "width_um", "idelta_smooth", "status") %in% colnames(final_curves)) &&
                    nrow(informative) > 0) {
                    final_plot <- final_curves |>
                        dplyr::filter(gene %in% informative$gene, status == "ok", is.finite(idelta_smooth)) |>
                        dplyr::select(gene, width_um, idelta_smooth)
                }

                if (!is.null(final_plot) && nrow(final_plot) > 0 && nrow(informative) > 0) {
                    gene_knee_with_y <- final_plot |>
                        dplyr::inner_join(informative[, c("gene", "knee_um")], by = "gene") |>
                        dplyr::group_by(gene) |>
                        dplyr::slice_min(abs(width_um - knee_um), n = 1, with_ties = FALSE) |>
                        dplyr::ungroup() |>
                        dplyr::select(gene, knee_um, knee_y = idelta_smooth)

                    vline_df <- data.frame(type = "Overall Knee (mean of gene knees)", x = overall_knee_um)

                    p_all <- ggplot2::ggplot(final_plot, ggplot2::aes(x = width_um, y = idelta_smooth, group = gene)) +
                        ggplot2::geom_line(alpha = 0.5, linewidth = 0.5, color = "gray60") +
                        ggplot2::geom_ribbon(
                            data = knee_range,
                            ggplot2::aes(x = x, ymin = ymin - 5, ymax = ymax + 5),
                            fill = blue_color,
                            alpha = 0.3,
                            linewidth = 0,
                            color = blue_color,
                            inherit.aes = FALSE
                        ) +
                        {
                            if (!is.null(avg_curve) && nrow(avg_curve) > 0) {
                                ggplot2::geom_line(
                                    data = avg_curve,
                                    ggplot2::aes(x = width_um, y = idelta_mean, group = 1),
                                    color = "black",
                                    linewidth = 1,
                                    inherit.aes = FALSE
                                )
                            }
                        } +
                        ggplot2::geom_point(
                            data = gene_knee_with_y,
                            ggplot2::aes(x = knee_um, y = knee_y, color = "Gene Wise Knee"),
                            size = 1,
                            alpha = 0.5,
                            show.legend = TRUE,
                            inherit.aes = FALSE
                        ) +
                        {
                            if (is.finite(overall_knee_um)) {
                                ggplot2::geom_vline(
                                    data = vline_df,
                                    ggplot2::aes(xintercept = x, color = type),
                                    linetype = "solid",
                                    linewidth = 1,
                                    show.legend = TRUE,
                                    inherit.aes = FALSE
                                )
                            }
                        } +
                        ggplot2::scale_color_manual(
                            name = NULL,
                            breaks = c("Gene Wise Knee", "Overall Knee (mean of gene knees)"),
                            values = c(
                                "Gene Wise Knee" = blue_color,
                                "Overall Knee (mean of gene knees)" = red_color
                            )
                        ) +
                        ggplot2::guides(
                            color = ggplot2::guide_legend(
                                override.aes = list(
                                    linetype = c(0, 1),
                                    shape = c(16, NA)
                                )
                            )
                        ) +
                        ggplot2::scale_x_reverse(
                            breaks = seq(
                                min(final_plot$width_um, na.rm = TRUE),
                                max(final_plot$width_um, na.rm = TRUE),
                                by = 10
                            )
                        ) +
                        ggplot2::labs(
                            title = "Idelta vs Width (scope)",
                            x = "Grid width (umm)",
                            y = "Idelta"
                        ) +
                        ggplot2::theme_minimal(base_size = 8) +
                        ggplot2::theme(
                            panel.border = ggplot2::element_rect(color = "black", fill = NA),
                            panel.background = ggplot2::element_rect(fill = "#c0c0c0", color = NA),
                            plot.background = ggplot2::element_rect(fill = "#ffffff", color = NA),
                            axis.title = ggplot2::element_text(size = 9),
                            axis.text.y = ggplot2::element_text(size = 8),
                            axis.text.x = ggplot2::element_text(size = 8, angle = 45, hjust = 1),
                            axis.ticks = ggplot2::element_line(color = "black"),
                            plot.title = ggplot2::element_text(hjust = 0.5, size = 10),
                            legend.position = c(0.02, 0.98),
                            legend.justification = c(0, 1),
                            legend.background = ggplot2::element_blank(),
                            legend.key = ggplot2::element_blank(),
                            legend.text = ggplot2::element_text(size = 7),
                            plot.margin = grid::unit(c(1, 1, 1, 1), "cm")
                        ) +
                        ggplot2::scale_y_continuous(
                            limits = c(0, min(100, max(final_plot$idelta_smooth, na.rm = TRUE))),
                            oob = scales::squish,
                            breaks = scales::breaks_width(5)
                        )
                } else {
                    p_all <- ggplot2::ggplot() +
                        ggplot2::theme_void() +
                        ggplot2::annotate("text", x = 0, y = 0, label = "Annotated plot unavailable")
                }

                ggplot2::ggsave(
                    filename = file.path(plots_out_dir, "scope_iDelta_all_genes_annotated_knee_stable.png"),
                    plot = p_all,
                    width = 6,
                    height = 5,
                    dpi = 600
                )

                # ---- 3) Gene-wise knee distribution (LQC-style) ----
                theme_lqc <- function() {
                    ggplot2::theme_bw(base_size = 8) +
                        ggplot2::theme(
                            plot.title = ggplot2::element_text(size = 12, colour = "black"),
                            panel.background = ggplot2::element_rect(fill = "#c0c0c0", colour = NA),
                            panel.grid.major.x = ggplot2::element_blank(),
                            panel.grid.minor = ggplot2::element_blank(),
                            plot.background = ggplot2::element_rect(fill = "#ffffff", color = NA),
                            panel.grid.major.y = ggplot2::element_line(colour = "#c0c0c0", linewidth = .2),
                            panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = .5),
                            axis.line = ggplot2::element_line(colour = "black", linewidth = .3),
                            axis.ticks = ggplot2::element_line(colour = "black", linewidth = .3),
                            axis.text = ggplot2::element_text(size = 10, colour = "black"),
                            axis.title = ggplot2::element_text(size = 12, colour = "black"),
                            plot.margin = grid::unit(c(0, 0, 0, 0), "cm")
                        )
                }

                knee_df <- data.frame(gene = informative$gene, knee_um = knee_vals)

                p_knee_violin <- if (nrow(knee_df) > 0) {
                    ggplot2::ggplot(knee_df, ggplot2::aes(x = "Genes", y = knee_um)) +
                        {
                            if (nrow(knee_df) >= 2) {
                                ggplot2::geom_violin(fill = "white", colour = "black", linewidth = .3, trim = FALSE)
                            }
                        } +
                        ggplot2::geom_boxplot(
                            width = .15,
                            outlier.size = .4,
                            outlier.stroke = .2,
                            fill = "white",
                            colour = "black",
                            linewidth = .25
                        ) +
                        {
                            if (nrow(knee_df) < 2) {
                                ggplot2::geom_point(
                                    shape = 21,
                                    size = 1.2,
                                    fill = "white",
                                    colour = "black",
                                    stroke = .25
                                )
                            }
                        } +
                        ggplot2::labs(title = "Gene Wise Knee (scope)", x = NULL, y = "Knee (umm)") +
                        theme_lqc()
                } else {
                    ggplot2::ggplot() +
                        ggplot2::theme_void() +
                        ggplot2::annotate("text", x = 0, y = 0, label = "No informative knees")
                }

                p_knee_hist <- if (nrow(knee_df) > 0) {
                    ggplot2::ggplot(knee_df, ggplot2::aes(x = knee_um)) +
                        ggplot2::geom_histogram(bins = 30, fill = "white", colour = "black", linewidth = .3) +
                        ggplot2::labs(title = "Gene Wise Knee Histogram (scope)", x = "Knee (umm)", y = "Frequency") +
                        theme_lqc()
                } else {
                    ggplot2::ggplot() +
                        ggplot2::theme_void() +
                        ggplot2::annotate("text", x = 0, y = 0, label = "No informative knees")
                }

                ggplot2::ggsave(
                    filename = file.path(plots_out_dir, "scope_geneWiseKnee_UIK_violin.png"),
                    plot = p_knee_violin,
                    width = 3,
                    height = 2.5,
                    dpi = 600
                )
                ggplot2::ggsave(
                    filename = file.path(plots_out_dir, "scope_geneWiseKnee_UIK_histogram.png"),
                    plot = p_knee_hist,
                    width = 3,
                    height = 2.5,
                    dpi = 600
                )
            }

        }
    }

    if (isTRUE(strict) && (!is.finite(n_genes_total) || !is.finite(n_genes_informative)) && file.exists(final_curves_path)) {
        warning("Tool summary is missing gene-count diagnostics; consider updating idelta-gridselect for degenerate-regime reporting.")
    }

    list(
        recommended_grid_um = recommended_grid_um,
        overall_knee_um = overall_knee_um,
        knee_of_mean_curve_um = knee_of_mean_curve_um,
        overall_knee_median_um = overall_knee_median_um,
        n_genes_total = n_genes_total,
        n_genes_informative = n_genes_informative,
        n_genes_excluded_degenerate = n_genes_excluded_degenerate,
        grid_qc_summary_path = NA_character_,
        grid_gene_bins_path = NA_character_,
        tool_out_dir = tool_out_dir,
        plots_out_dir = plots_out_dir
    )
}

#' Recommend an optimal grid size using the external Idelta grid selector
#'
#' One-step wrapper that can (i) run the external `idelta-gridselect recommend-grid`
#' tool on Xenium inputs and (ii) generate best_grid-style plots + summary.
#'
#' If `xenium_dir`/`transcripts_path`/`roi_csv` are not provided, this function
#' instead consumes existing outputs under `file.path(out_dir, "tool")` and
#' regenerates plots/summary under `file.path(out_dir, "plots")`.
#'
#' @param xenium_dir Xenium outs directory containing `transcripts.csv.gz` (or `transcripts.csv` / `transcripts.parquet`).
#' @param transcripts_path Optional explicit transcripts path (overrides `xenium_dir`).
#' @param roi_csv ROI vertices CSV/TSV (comments allowed; header may be `X,Y`). If NULL, a bbox ROI is auto-generated.
#' @param bbox_csv Optional coordinate CSV/TSV used to auto-generate a bbox ROI when `roi_csv` is NULL (e.g., Xenium `cells.csv.gz`).
#' @param bbox_pad_um Optional padding (umm) applied to the auto-generated bbox ROI (default 0).
#' @param out_dir Output directory root; creates `staging/`, `tool/`, `plots/`, `logs/`.
#' @param tool_bin Optional path to `idelta-gridselect` binary (highest priority).
#' @param docker_image Optional docker image to run (next priority).
#' @param apptainer_sif Optional Apptainer/Singularity SIF path to run (next priority).
#' @param auto_find Logical; if TRUE and no tool is provided, search common local build locations.
#' @param threads Number of threads passed to the external tool (`--threads`).
#' @param knee_search_min_um Minimum grid width considered when selecting the knee.
#' @param descending_delta_um Minimum descending-grid delta used by the external selector.
#' @param max_width_um Optional maximum grid width for selector search.
#' @param smooth_method Optional smoothing method passed through to the external selector.
#' @param grid_um_eval Optional grid size for immediate QC evaluation after selection.
#' @param epsilon Near-zero threshold for I-delta QC summaries.
#' @param q_bins Number of quantile tiers used by QC summaries.
#' @param ratio_high Knee scale mismatch threshold used by QC summaries.
#' @param loss_min Loss flag threshold used by QC summaries.
#' @param ... Extra options: `strict`, `best_grid_style`, `search_roots`, `skip_tool`, `skip_plots`.
#'
#' @return A list with the recommended grid size and output paths.
#' @export
gridSelect <- function(
    xenium_dir = NULL,
    transcripts_path = NULL,
    roi_csv = NULL,
    bbox_csv = NULL,
    bbox_pad_um = 0,
    out_dir,
    tool_bin = NULL,
    docker_image = NULL,
    apptainer_sif = NULL,
    auto_find = TRUE,
    threads = NULL,
    knee_search_min_um = 5,
    descending_delta_um = 2,
    max_width_um = NULL,
    smooth_method = NULL,
    grid_um_eval = NULL,
    epsilon = 1,
    q_bins = 5,
    ratio_high = 2,
    loss_min = 0.2,
    ...
) {
    dots <- list(...)
    strict <- .gridselect_bool(.gridselect_coalesce(dots$strict, TRUE), default = TRUE)
    best_grid_style <- .gridselect_bool(.gridselect_coalesce(dots$best_grid_style, TRUE), default = TRUE)
    skip_tool <- .gridselect_bool(.gridselect_coalesce(dots$skip_tool, FALSE), default = FALSE)
    skip_plots <- .gridselect_bool(.gridselect_coalesce(dots$skip_plots, FALSE), default = FALSE)
    search_roots <- .gridselect_coalesce(dots$search_roots, character(0))
    if (is.null(search_roots)) search_roots <- character(0)
    if (!is.character(search_roots)) search_roots <- as.character(search_roots)
    skip_molecules <- .gridselect_bool(.gridselect_coalesce(dots$skip_molecules, FALSE), default = FALSE)
    repo_pin <- .gridselect_coalesce(dots$repo_pin, NULL)
    bbox_csv_eff <- .gridselect_coalesce(
        bbox_csv,
        .gridselect_coalesce(
            dots$bbox_csv,
            .gridselect_coalesce(dots$coord_csv, .gridselect_coalesce(dots$coord_file, NULL))
        )
    )

    if (missing(out_dir) || !is.character(out_dir) || length(out_dir) != 1 || !nzchar(out_dir)) {
        stop("`out_dir` is required.")
    }
    out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)

    staging_dir <- file.path(out_dir, "staging")
    logs_dir <- file.path(out_dir, "logs")
    tool_out_dir <- file.path(out_dir, "tool")
    plots_out_dir <- file.path(out_dir, "plots")

    dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(tool_out_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plots_out_dir, recursive = TRUE, showWarnings = FALSE)

    wrapper_log_path <- file.path(logs_dir, "wrapper.log")
    resolve_path <- file.path(logs_dir, "idelta_gridselect.resolve.txt")
    args_stamp_path <- file.path(logs_dir, "idelta_gridselect.args.txt")
    tool_log_path <- file.path(logs_dir, "idelta_gridselect.log")
    summary_log_path <- .gridselect_coalesce(dots$log_out_path, file.path(logs_dir, "best_grid_summary.txt"))

    t0 <- Sys.time()
    .gridselect_wlog(wrapper_log_path, paste0("out_dir=", out_dir), echo = FALSE)

    if (isTRUE(skip_plots)) best_grid_style <- FALSE
    knee_search_min_um_num <- suppressWarnings(as.numeric(knee_search_min_um))[1]
    descending_delta_um_num <- suppressWarnings(as.numeric(descending_delta_um))[1]
    if (!is.finite(knee_search_min_um_num) || abs(knee_search_min_um_num - round(knee_search_min_um_num)) > 1e-9) {
        stop("`knee_search_min_um` must be an integer.")
    }
    if (!is.finite(descending_delta_um_num) || abs(descending_delta_um_num - round(descending_delta_um_num)) > 1e-9) {
        stop("`descending_delta_um` must be an integer.")
    }
    knee_search_min_um <- as.integer(round(knee_search_min_um_num))
    descending_delta_um <- as.integer(round(descending_delta_um_num))
    if (!is.finite(knee_search_min_um) || knee_search_min_um < 0) {
        stop("`knee_search_min_um` must be a non-negative integer.")
    }
    if (!is.finite(descending_delta_um) || descending_delta_um < 0) {
        stop("`descending_delta_um` must be a non-negative integer.")
    }

    max_width_um_num <- suppressWarnings(as.numeric(max_width_um))[1]
    if (is.finite(max_width_um_num)) {
        if (abs(max_width_um_num - round(max_width_um_num)) > 1e-9) {
            stop("`max_width_um` must be an integer (or NULL).")
        }
        max_width_um <- as.integer(round(max_width_um_num))
        if (!is.finite(max_width_um) || max_width_um < 1L) stop("`max_width_um` must be >= 1 (or NULL).")
        .gridselect_wlog(wrapper_log_path, paste0("max_width_um=", max_width_um), echo = FALSE)
    } else {
        max_width_um <- NULL
    }

    grid_um_eval_num <- suppressWarnings(as.numeric(grid_um_eval))[1]
    if (!is.null(grid_um_eval)) {
        if (!is.finite(grid_um_eval_num)) stop("`grid_um_eval` must be a finite integer (or NULL).")
        if (abs(grid_um_eval_num - round(grid_um_eval_num)) > 1e-9) stop("`grid_um_eval` must be an integer (or NULL).")
        grid_um_eval <- as.integer(round(grid_um_eval_num))
        if (!is.finite(grid_um_eval) || grid_um_eval < 1L) stop("`grid_um_eval` must be >= 1 (or NULL).")
    }

    epsilon_num <- suppressWarnings(as.numeric(epsilon))[1]
    if (!is.finite(epsilon_num) || epsilon_num < 0) stop("`epsilon` must be a finite number >= 0.")
    epsilon <- epsilon_num

    q_bins_num <- suppressWarnings(as.numeric(q_bins))[1]
    if (!is.finite(q_bins_num) || abs(q_bins_num - round(q_bins_num)) > 1e-9) stop("`q_bins` must be an integer >= 2.")
    q_bins <- as.integer(round(q_bins_num))
    if (!is.finite(q_bins) || q_bins < 2L) stop("`q_bins` must be an integer >= 2.")

    ratio_high_num <- suppressWarnings(as.numeric(ratio_high))[1]
    if (!is.finite(ratio_high_num) || ratio_high_num <= 1) stop("`ratio_high` must be a finite number > 1.")
    ratio_high <- ratio_high_num

    loss_min_num <- suppressWarnings(as.numeric(loss_min))[1]
    if (!is.finite(loss_min_num) || loss_min_num < 0 || loss_min_num > 1) stop("`loss_min` must be between 0 and 1.")
    loss_min <- loss_min_num

    # Mode selection:
    # - If xenium_dir/transcripts_path are provided, run the external tool (prefer
    #   a locally available idelta-gridselect binary/docker/apptainer; fallback to repo runner).
    # - Otherwise, consume existing tool outputs under out_dir/tool.
    has_xenium <- is.character(xenium_dir) && length(xenium_dir) == 1 && nzchar(xenium_dir)
    has_transcripts <- is.character(transcripts_path) && length(transcripts_path) == 1 && nzchar(transcripts_path)
    wants_tool_run <- has_xenium || has_transcripts
    if (has_xenium) xenium_dir <- .gridselect_resolve_xenium_dir(xenium_dir)
    tool_log_path_eff <- tool_log_path

    if (wants_tool_run) {
        if (has_xenium && has_transcripts) stop("Provide only one of `xenium_dir` or `transcripts_path`.")

        # Resolve/validate transcripts inputs early so both tool modes (direct tool vs repo runner)
        # behave consistently.
        transcripts_eff <- NULL
        xenium_dir_eff <- NULL
        if (has_xenium) {
            xenium_dir_eff <- normalizePath(xenium_dir, winslash = "/", mustWork = FALSE)
            if (!dir.exists(xenium_dir_eff)) stop("xenium_dir not found: ", xenium_dir_eff)
            tx <- .gridselect_guess_xenium_file(xenium_dir_eff, "transcripts")
            if (!is.character(tx) || length(tx) != 1 || is.na(tx) || !nzchar(tx)) {
                stop(
                    "Could not find Xenium transcripts under: ", xenium_dir_eff,
                    "\nExpected one of: transcripts.csv.gz, transcripts.csv, transcripts.parquet"
                )
            }
            if (grepl("\\.zarr(\\.zip)?$", tx, ignore.case = TRUE) || (dir.exists(tx) && grepl("\\.zarr$", tx, ignore.case = TRUE))) {
                stop(
                    "Xenium transcripts appear to be Zarr (", basename(tx), "), which is not supported.\n",
                    "Please export `transcripts.parquet` or `transcripts.csv.gz` and rerun, or provide `transcripts_path=` pointing to a CSV/Parquet file."
                )
            }
            transcripts_eff <- tx
        } else {
            transcripts_eff <- normalizePath(transcripts_path, winslash = "/", mustWork = FALSE)
            if (!(file.exists(transcripts_eff) || dir.exists(transcripts_eff))) {
                stop("transcripts_path not found: ", transcripts_eff)
            }
            if (grepl("\\.zarr(\\.zip)?$", transcripts_eff, ignore.case = TRUE) ||
                (dir.exists(transcripts_eff) && grepl("\\.zarr$", transcripts_eff, ignore.case = TRUE))) {
                stop(
                    "transcripts_path points to a Zarr store, which is not supported: ", transcripts_eff,
                    "\nProvide a `transcripts.parquet` or `transcripts.csv(.gz)` file instead."
                )
            }
        }

        roi_csv_eff <- .gridselect_resolve_roi_csv_for_tool_run(
            roi_csv = roi_csv,
            xenium_dir = if (has_xenium) xenium_dir_eff else NULL,
            transcripts_path = transcripts_eff,
            staging_dir = staging_dir,
            wrapper_log_path = wrapper_log_path,
            bbox_csv = bbox_csv_eff,
            bbox_pad_um = bbox_pad_um,
            strict = strict
        )

        threads_eff <- if (is.null(threads)) max(1L, as.integer(parallel::detectCores())) else threads

        # Normalize ROI vertices for tool execution (direct tool + repo runner both expect x_um/y_um CSV).
        roi_csv_in <- normalizePath(roi_csv_eff, winslash = "/", mustWork = TRUE)
        roi_vertices_path <- file.path(staging_dir, "roi_vertices.csv")
        if (!file.exists(roi_vertices_path) ||
            file.info(roi_vertices_path)$mtime < file.info(roi_csv_in)$mtime) {
            .gridselect_normalize_roi_vertices(roi_csv_in, roi_vertices_path)
        }

        # Prefer running a locally-available tool (tool_bin/env/PATH/docker/apptainer/auto_find).
        tool_requested <- !is.null(tool_bin) || !is.null(docker_image) || !is.null(apptainer_sif) ||
            nzchar(Sys.getenv("IDELTA_GRIDSELECT_BIN", unset = "")) ||
            nzchar(Sys.getenv("IDELTA_GRIDSELECT_DOCKER_IMAGE", unset = "")) ||
            nzchar(Sys.getenv("IDELTA_GRIDSELECT_APPTAINER_SIF", unset = ""))

        tool_try <- tryCatch(
            .gridselect_resolve_tool(
                out_dir = out_dir,
                tool_bin = tool_bin,
                docker_image = docker_image,
                apptainer_sif = apptainer_sif,
                auto_find = auto_find,
                search_roots = search_roots,
                resolve_path = resolve_path,
                wrapper_log_path = wrapper_log_path
            ),
            error = function(e) e
        )

        if (!inherits(tool_try, "error")) {
            tool <- tool_try
            .gridselect_wlog(
                wrapper_log_path,
                paste0("using direct idelta-gridselect execution mode=", tool$mode, " source=", tool$source),
                echo = FALSE
            )
            tool_out <- .gridselect_run_idelta_gridselect_tool(
                tool = tool,
                out_dir = out_dir,
                transcripts_path = transcripts_eff,
                roi_vertices_path = roi_vertices_path,
                threads = threads_eff,
                knee_search_min_um = knee_search_min_um,
                descending_delta_um = descending_delta_um,
                max_width_um = max_width_um,
                skip_molecules = skip_molecules,
                skip_tool = skip_tool,
                wrapper_log_path = wrapper_log_path,
                strict = strict
            )
            tool_log_path_eff <- tool_out$tool_log
        } else {
            if (isTRUE(tool_requested)) stop(tool_try)
            .gridselect_wlog(
                wrapper_log_path,
                paste0("direct tool not available; falling back to repo runner via repo pin (", conditionMessage(tool_try), ")"),
                echo = FALSE
            )
            tool_run <- .gridselect_run_idelta_gridselect_runner(
                out_dir = out_dir,
                xenium_dir = if (has_xenium) xenium_dir_eff else NULL,
                transcripts_path = if (!has_xenium) transcripts_eff else NULL,
                roi_csv = roi_vertices_path,
                threads = threads_eff,
                knee_search_min_um = knee_search_min_um,
                descending_delta_um = descending_delta_um,
                max_width_um = max_width_um,
                skip_molecules = skip_molecules,
                skip_tool = skip_tool,
                repo_pin = repo_pin,
                strict = strict
            )
            tool_log_path_eff <- tool_run$runner_log

            resolve_lines <- c(
                "MODE=repo_clone_check+runner",
                paste0("REPO_PATH=", tool_run$repo_path),
                paste0("RUNNER=", tool_run$runner),
                paste0("RUNNER_LOG=", tool_run$runner_log),
                paste0(
                    "RECOMMENDED_GRID_UM_RUNNER=",
                    ifelse(is.finite(tool_run$recommended_grid_um_runner), tool_run$recommended_grid_um_runner, "NA")
                )
            )
            writeLines(resolve_lines, con = resolve_path, useBytes = TRUE)
        }
    } else {
        .gridselect_wlog(wrapper_log_path, "no xenium inputs provided; consuming existing tool outputs under out_dir/tool", echo = FALSE)
        if (!file.exists(file.path(tool_out_dir, "knee_summary.tsv"))) {
            stop("Missing tool outputs under out_dir/tool. Provide xenium_dir/transcripts_path to run the tool.")
        }
    }

    # Consume outputs + generate plots/summary.
    res <- .gridselect_consume_tool_outputs(
        tool_out_dir = tool_out_dir,
        plots_out_dir = plots_out_dir,
        log_out_path = summary_log_path,
        best_grid_style = best_grid_style,
        strict = strict,
        smooth_method = smooth_method
    )

    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    .gridselect_wlog(wrapper_log_path, paste0("done seconds=", sprintf("%.1f", dt)), echo = FALSE)

    res$out_dir <- out_dir
    res$logs_dir <- logs_dir
    res$tool_log_path <- tool_log_path_eff
    res$resolve_path <- resolve_path
    res
}

# Backward-compatible alias (legacy API name).
recommendGridExternal <- function(...) gridSelect(...)

# Backward-compatible alias (legacy snake_case name).
recommend_grid_external <- function(...) gridSelect(...)

#' Compute grid-at-choice Idelta QC (dispersion, tiers, unsuitable flags) from external tool outputs
#'
#' This is a post-recommendation diagnostic step: it consumes `knee_summary.tsv`
#' and `final_curves.tsv.gz` under `tool_out_dir` (or `<out_dir>/tool`) and writes
#' QC tables + plots under `plots_out_dir` (or `<out_dir>/plots`).
#'
#' @param out_dir Optional wrapper output root that contains `tool/` and `plots/`.
#' @param xenium_dir Optional Xenium outs directory; when supplied the external tool is run first.
#' @param transcripts_path Optional explicit transcripts path; when supplied the external tool is run first.
#' @param roi_csv Optional ROI vertices CSV/TSV passed to the external tool.
#' @param bbox_csv Optional coordinate CSV/TSV used for bbox ROI generation.
#' @param bbox_pad_um Optional padding applied to generated bbox ROI coordinates.
#' @param threads Number of threads passed to the external tool.
#' @param knee_search_min_um Minimum grid width considered when selecting the knee.
#' @param descending_delta_um Minimum descending-grid delta used by the external selector.
#' @param max_width_um Optional maximum grid width for selector search.
#' @param skip_molecules Skip molecule preprocessing when external tool outputs already exist.
#' @param skip_tool Skip external tool execution and consume existing outputs.
#' @param repo_pin Optional source repository pin recorded in generated QC metadata.
#' @param tool_out_dir Directory containing external tool outputs (`knee_summary.tsv`, `final_curves.tsv.gz`).
#' @param plots_out_dir Output directory for QC tables and plots.
#' @param grid_um_eval Grid size (umm) to evaluate; defaults to `recommended_grid_um` from `knee_summary.tsv`.
#' @param epsilon Near-zero threshold for `idelta_at_grid` (default 1).
#' @param q_bins Number of quantile tiers (default 5).
#' @param ratio_high Knee scale mismatch threshold (default 2.0), comparing `knee_um` to `grid_um_eval`.
#' @param loss_min Loss flag threshold (default 0.2).
#' @param best_grid_style When TRUE, generate best_grid-style QC plots (requires ggplot2).
#' @param strict When TRUE, emit warnings on missing fields or plotting dependencies.
#' @param log_out_path Optional text log path to write the QC summary block.
#' @param flag Optional free-text flag written to the QC summary output.
#'
#' @return A list with `qc_summary`, `gene_bins`, file paths, and plot paths.
#' @export
gridSelectQC <- function(
    out_dir = NULL,
    xenium_dir = NULL,
    transcripts_path = NULL,
    roi_csv = NULL,
    bbox_csv = NULL,
    bbox_pad_um = 0,
    threads = 32,
    knee_search_min_um = 5,
    descending_delta_um = 2,
    max_width_um = NULL,
    skip_molecules = FALSE,
    skip_tool = FALSE,
    repo_pin = NULL,
    tool_out_dir = NULL,
    plots_out_dir = NULL,
    grid_um_eval = NULL,
    epsilon = 1,
    q_bins = 5,
    ratio_high = 2,
    loss_min = 0.2,
    best_grid_style = TRUE,
    strict = TRUE,
    log_out_path = NULL,
    flag = NULL
) {
    # Optional tool run/bootstrap step (delegated to scripts/external_tools/*).
    has_xenium <- is.character(xenium_dir) && length(xenium_dir) == 1 && nzchar(xenium_dir)
    has_transcripts <- is.character(transcripts_path) && length(transcripts_path) == 1 && nzchar(transcripts_path)
    wants_tool_run <- has_xenium || has_transcripts
    if (has_xenium) xenium_dir <- .gridselect_resolve_xenium_dir(xenium_dir)
    recommended_grid_um_runner <- NA_real_
    if (wants_tool_run) {
        if (has_xenium && has_transcripts) stop("Provide only one of `xenium_dir` or `transcripts_path`.")
        if (is.null(out_dir) || !is.character(out_dir) || length(out_dir) != 1 || !nzchar(out_dir)) {
            stop("When running the tool, `out_dir` is required (tool outputs are written under out_dir/tool).")
        }
        out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
        logs_dir <- file.path(out_dir, "logs")
        staging_dir <- file.path(out_dir, "staging")
        dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
        wrapper_log_path <- file.path(logs_dir, "wrapper.log")
        resolve_path <- file.path(logs_dir, "idelta_gridselect.resolve.txt")

        # Resolve transcripts to a concrete file path for direct-tool execution.
        transcripts_eff <- NULL
        xenium_dir_eff <- NULL
        if (has_xenium) {
            xenium_dir_eff <- normalizePath(xenium_dir, winslash = "/", mustWork = FALSE)
            if (!dir.exists(xenium_dir_eff)) stop("xenium_dir not found: ", xenium_dir_eff)
            tx <- .gridselect_guess_xenium_file(xenium_dir_eff, "transcripts")
            if (!is.character(tx) || length(tx) != 1 || is.na(tx) || !nzchar(tx)) {
                stop(
                    "Could not find Xenium transcripts under: ", xenium_dir_eff,
                    "\nExpected one of: transcripts.csv.gz, transcripts.csv, transcripts.parquet"
                )
            }
            if (grepl("\\.zarr(\\.zip)?$", tx, ignore.case = TRUE) || (dir.exists(tx) && grepl("\\.zarr$", tx, ignore.case = TRUE))) {
                stop(
                    "Xenium transcripts appear to be Zarr (", basename(tx), "), which is not supported.\n",
                    "Please export `transcripts.parquet` or `transcripts.csv.gz` and rerun, or provide `transcripts_path=` pointing to a CSV/Parquet file."
                )
            }
            transcripts_eff <- tx
        } else {
            transcripts_eff <- normalizePath(transcripts_path, winslash = "/", mustWork = FALSE)
            if (!(file.exists(transcripts_eff) || dir.exists(transcripts_eff))) stop("transcripts_path not found: ", transcripts_eff)
            if (grepl("\\.zarr(\\.zip)?$", transcripts_eff, ignore.case = TRUE) ||
                (dir.exists(transcripts_eff) && grepl("\\.zarr$", transcripts_eff, ignore.case = TRUE))) {
                stop(
                    "transcripts_path points to a Zarr store, which is not supported: ", transcripts_eff,
                    "\nProvide a `transcripts.parquet` or `transcripts.csv(.gz)` file instead."
                )
            }
        }

        roi_csv_eff <- .gridselect_resolve_roi_csv_for_tool_run(
            roi_csv = roi_csv,
            xenium_dir = if (has_xenium) xenium_dir_eff else NULL,
            transcripts_path = transcripts_eff,
            staging_dir = staging_dir,
            wrapper_log_path = wrapper_log_path,
            bbox_csv = bbox_csv,
            bbox_pad_um = bbox_pad_um,
            strict = strict
        )

        # Normalize ROI vertices to x_um/y_um CSV.
        roi_csv_in <- normalizePath(roi_csv_eff, winslash = "/", mustWork = TRUE)
        roi_vertices_path <- file.path(staging_dir, "roi_vertices.csv")
        if (!file.exists(roi_vertices_path) ||
            file.info(roi_vertices_path)$mtime < file.info(roi_csv_in)$mtime) {
            .gridselect_normalize_roi_vertices(roi_csv_in, roi_vertices_path)
        }

        # Prefer direct tool if available (PATH/env/docker/apptainer); otherwise fallback to repo runner.
        tool_requested <- nzchar(Sys.getenv("IDELTA_GRIDSELECT_BIN", unset = "")) ||
            nzchar(Sys.getenv("IDELTA_GRIDSELECT_DOCKER_IMAGE", unset = "")) ||
            nzchar(Sys.getenv("IDELTA_GRIDSELECT_APPTAINER_SIF", unset = ""))

        tool_try <- tryCatch(
            .gridselect_resolve_tool(
                out_dir = out_dir,
                tool_bin = NULL,
                docker_image = NULL,
                apptainer_sif = NULL,
                auto_find = TRUE,
                search_roots = character(0),
                resolve_path = resolve_path,
                wrapper_log_path = wrapper_log_path
            ),
            error = function(e) e
        )

        if (!inherits(tool_try, "error")) {
            tool <- tool_try
            .gridselect_wlog(
                wrapper_log_path,
                paste0("using direct idelta-gridselect execution mode=", tool$mode, " source=", tool$source),
                echo = FALSE
            )
            .gridselect_run_idelta_gridselect_tool(
                tool = tool,
                out_dir = out_dir,
                transcripts_path = transcripts_eff,
                roi_vertices_path = roi_vertices_path,
                threads = threads,
                knee_search_min_um = knee_search_min_um,
                descending_delta_um = descending_delta_um,
                max_width_um = max_width_um,
                skip_molecules = skip_molecules,
                skip_tool = skip_tool,
                wrapper_log_path = wrapper_log_path,
                strict = strict
            )
        } else {
            if (isTRUE(tool_requested)) stop(tool_try)
            .gridselect_wlog(
                wrapper_log_path,
                paste0("direct tool not available; falling back to repo runner via repo pin (", conditionMessage(tool_try), ")"),
                echo = FALSE
            )
            tool_run <- .gridselect_run_idelta_gridselect_runner(
                out_dir = out_dir,
                xenium_dir = if (has_xenium) xenium_dir_eff else NULL,
                transcripts_path = if (!has_xenium) transcripts_eff else NULL,
                roi_csv = roi_vertices_path,
                threads = threads,
                knee_search_min_um = knee_search_min_um,
                descending_delta_um = descending_delta_um,
                max_width_um = max_width_um,
                skip_molecules = skip_molecules,
                skip_tool = skip_tool,
                repo_pin = repo_pin,
                strict = strict
            )
            recommended_grid_um_runner <- tool_run$recommended_grid_um_runner

            resolve_lines <- c(
                "MODE=repo_clone_check+runner",
                paste0("REPO_PATH=", tool_run$repo_path),
                paste0("RUNNER=", tool_run$runner),
                paste0("RUNNER_LOG=", tool_run$runner_log),
                paste0(
                    "RECOMMENDED_GRID_UM_RUNNER=",
                    ifelse(is.finite(tool_run$recommended_grid_um_runner), tool_run$recommended_grid_um_runner, "NA")
                )
            )
            writeLines(resolve_lines, con = resolve_path, useBytes = TRUE)
        }
    }

    if (!is.null(out_dir)) {
        if (!is.character(out_dir) || length(out_dir) != 1 || !nzchar(out_dir)) stop("`out_dir` must be a single path (or NULL).")
        out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
        if (is.null(tool_out_dir)) tool_out_dir <- file.path(out_dir, "tool")
        if (is.null(plots_out_dir)) plots_out_dir <- file.path(out_dir, "plots")
    }
    if (is.null(tool_out_dir) || !is.character(tool_out_dir) || length(tool_out_dir) != 1 || !nzchar(tool_out_dir)) {
        stop("Provide `tool_out_dir` (or `out_dir`).")
    }
    if (is.null(plots_out_dir) || !is.character(plots_out_dir) || length(plots_out_dir) != 1 || !nzchar(plots_out_dir)) {
        stop("Provide `plots_out_dir` (or `out_dir`).")
    }

    tool_out_dir <- normalizePath(tool_out_dir, winslash = "/", mustWork = TRUE)
    plots_out_dir <- normalizePath(plots_out_dir, winslash = "/", mustWork = FALSE)
    dir.create(plots_out_dir, recursive = TRUE, showWarnings = FALSE)

    qc <- compute_grid_idelta_qc(
        tool_out_dir = tool_out_dir,
        grid_um_eval = grid_um_eval,
        epsilon = epsilon,
        q_bins = q_bins,
        ratio_high = ratio_high,
        loss_min = loss_min
    )

    qc_summary_path <- file.path(plots_out_dir, "grid_qc_summary.tsv")
    gene_bins_path <- file.path(plots_out_dir, "grid_gene_bins.tsv.gz")

    utils::write.table(
        qc$qc_summary,
        file = qc_summary_path,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )
    con <- gzfile(gene_bins_path, open = "wt")
    on.exit(close(con), add = TRUE)
    utils::write.table(
        qc$gene_bins,
        file = con,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )

    to_num <- function(x) suppressWarnings(as.numeric(x))
    qs <- qc$qc_summary[1, , drop = FALSE]
    q_bins_out <- suppressWarnings(as.integer(qs$q_bins))[1]
    if (!is.finite(q_bins_out) || q_bins_out < 2L) q_bins_out <- suppressWarnings(as.integer(q_bins))[1]
    if (!is.finite(q_bins_out) || q_bins_out < 2L) q_bins_out <- 5L

    tier_parts <- character(0)
    for (i in seq_len(q_bins_out)) {
        nm <- paste0("tier_Q", i)
        tier_parts <- c(tier_parts, sprintf("%s=%s", nm, if (nm %in% colnames(qs)) suppressWarnings(as.integer(qs[[nm]]))[1] else NA_integer_))
    }
    out_dir_print <- out_dir
    if (is.null(out_dir_print) || !is.character(out_dir_print) || length(out_dir_print) != 1 || !nzchar(out_dir_print)) {
        if (is.character(tool_out_dir) && length(tool_out_dir) == 1 && nzchar(tool_out_dir) && identical(basename(tool_out_dir), "tool")) {
            out_dir_print <- dirname(tool_out_dir)
        } else if (is.character(plots_out_dir) && length(plots_out_dir) == 1 && nzchar(plots_out_dir) && identical(basename(plots_out_dir), "plots")) {
            out_dir_print <- dirname(plots_out_dir)
        } else {
            out_dir_print <- tool_out_dir
        }
    }
    out_dir_print <- normalizePath(out_dir_print, winslash = "/", mustWork = FALSE)

    qc_lines <- c(
        sprintf("out_dir: %s", out_dir_print),
        sprintf("grid_um_eval: %s", ifelse(is.finite(to_num(qs$grid_um_eval)), to_num(qs$grid_um_eval), "NA")),
        sprintf(
            "idelta_at_grid median/IQR: %s / %s",
            ifelse(is.finite(to_num(qs$idelta_median)), to_num(qs$idelta_median), "NA"),
            ifelse(is.finite(to_num(qs$idelta_iqr)), to_num(qs$idelta_iqr), "NA")
        ),
        sprintf(
            "near_zero_frac (idelta < epsilon=%s): %s",
            ifelse(is.finite(to_num(qs$epsilon)), to_num(qs$epsilon), "NA"),
            ifelse(is.finite(to_num(qs$frac_near_zero)), sprintf("%.3f", to_num(qs$frac_near_zero)), "NA")
        ),
        sprintf("tiers (q_bins=%s): %s", q_bins_out, paste(tier_parts, collapse = " ")),
        sprintf(
            "unsuitable flags (near_zero/scale_mismatch/loss/any): %s/%s/%s/%s",
            ifelse(is.finite(to_num(qs$n_flag_near_zero)), as.integer(qs$n_flag_near_zero), "NA"),
            ifelse(is.finite(to_num(qs$n_flag_scale_mismatch)), as.integer(qs$n_flag_scale_mismatch), "NA"),
            ifelse(is.finite(to_num(qs$n_flag_loss)), as.integer(qs$n_flag_loss), "NA"),
            ifelse(is.finite(to_num(qs$n_unsuitable_any)), as.integer(qs$n_unsuitable_any), "NA")
        )
    )

    # Optional: print gene lists for selected flags/tiers.
    # `flag` may include any of: near_zero/scale_mismatch/loss/any or tier labels Q1..Qk.
    flag_lines <- character(0)
    if (!is.null(flag)) {
        flags <- flag
        if (!is.character(flags)) flags <- as.character(flags)
        flags <- trimws(flags)
        flags <- flags[!is.na(flags) & nzchar(flags)]
        if (length(flags) > 0) {
            seen <- character(0)
            genes <- as.character(qc$gene_bins$gene)
            tier <- if ("tier_q5" %in% colnames(qc$gene_bins)) as.character(qc$gene_bins$tier_q5) else rep(NA_character_, length(genes))
            flag_near_zero <- if ("flag_near_zero" %in% colnames(qc$gene_bins)) qc$gene_bins$flag_near_zero %in% TRUE else rep(FALSE, length(genes))
            flag_scale_mismatch <- if ("flag_scale_mismatch" %in% colnames(qc$gene_bins)) qc$gene_bins$flag_scale_mismatch %in% TRUE else rep(FALSE, length(genes))
            flag_loss <- if ("flag_loss" %in% colnames(qc$gene_bins)) qc$gene_bins$flag_loss %in% TRUE else rep(FALSE, length(genes))
            unsuitable_any <- if ("unsuitable_any" %in% colnames(qc$gene_bins)) qc$gene_bins$unsuitable_any %in% TRUE else rep(FALSE, length(genes))

            quote_join <- function(x) {
                x <- unique(as.character(x))
                x <- x[!is.na(x) & nzchar(x)]
                x <- sort(x)
                if (length(x) == 0) return("")
                paste0("\"", x, "\"", collapse = ",")
            }

            for (f in flags) {
                f_norm <- tolower(f)
                label <- NULL
                mask <- NULL

                if (f_norm %in% c("near_zero", "scale_mismatch", "loss", "any")) {
                    label <- f_norm
                    mask <- switch(
                        f_norm,
                        near_zero = flag_near_zero,
                        scale_mismatch = flag_scale_mismatch,
                        loss = flag_loss,
                        any = unsuitable_any
                    )
                } else if (grepl("^q[0-9]+$", f_norm)) {
                    k <- suppressWarnings(as.integer(sub("^q", "", f_norm)))
                    if (!is.finite(k) || k < 1L) stop("Invalid tier flag: ", f, " (use Q1..Q", q_bins_out, ")")
                    if (is.finite(q_bins_out) && k > q_bins_out) stop("Invalid tier flag: ", f, " (q_bins=", q_bins_out, ")")
                    label <- paste0("Q", k)
                    mask <- !is.na(tier) & tier == label
                } else {
                    stop("Unsupported flag: ", f, " (use near_zero/scale_mismatch/loss/any or Q1..Q", q_bins_out, ")")
                }

                if (label %in% seen) next
                seen <- c(seen, label)
                sel <- if (!is.null(mask)) genes[which(mask)] else character(0)
                flag_lines <- c(flag_lines, paste0(label, ":", quote_join(sel)))
            }
        }
    }

    if (length(flag_lines) > 0) qc_lines <- c(qc_lines, flag_lines)
    cat(paste0(qc_lines, collapse = "\n"), "\n", sep = "")
    if (!is.null(log_out_path) && is.character(log_out_path) && length(log_out_path) == 1 && nzchar(log_out_path)) {
        dir.create(dirname(log_out_path), recursive = TRUE, showWarnings = FALSE)
        writeLines(qc_lines, con = log_out_path, useBytes = TRUE)
    }

    plot_paths <- character(0)
    if (isTRUE(best_grid_style)) {
        have_plot_pkgs <- requireNamespace("ggplot2", quietly = TRUE) &&
            requireNamespace("dplyr", quietly = TRUE) &&
            requireNamespace("tibble", quietly = TRUE)
        if (!isTRUE(have_plot_pkgs)) {
            if (isTRUE(strict)) warning("Missing ggplot2/dplyr/tibble; skipping QC plots.")
        } else {
            red_color <- "#E41A1C"
            blue_color <- "#377EB8"

            theme_lqc_qc <- function() {
                ggplot2::theme_bw(base_size = 8) +
                    ggplot2::theme(
                        plot.title = ggplot2::element_text(size = 12, colour = "black"),
                        panel.background = ggplot2::element_rect(fill = "#c0c0c0", colour = NA),
                        panel.grid.major.x = ggplot2::element_blank(),
                        panel.grid.minor = ggplot2::element_blank(),
                        plot.background = ggplot2::element_rect(fill = "#ffffff", color = NA),
                        panel.grid.major.y = ggplot2::element_line(colour = "#c0c0c0", linewidth = .2),
                        panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = .5),
                        axis.line = ggplot2::element_line(colour = "black", linewidth = .3),
                        axis.ticks = ggplot2::element_line(colour = "black", linewidth = .3),
                        axis.text = ggplot2::element_text(size = 10, colour = "black"),
                        axis.title = ggplot2::element_text(size = 12, colour = "black"),
                        plot.margin = grid::unit(c(0, 0, 0, 0), "cm")
                    )
            }

            qc_bins <- qc$gene_bins
            qc_grid <- suppressWarnings(as.numeric(qc$qc_summary$grid_um_eval))[1]
            if (!is.finite(qc_grid)) qc_grid <- suppressWarnings(as.numeric(qc$qc_summary$recommended_grid_um))[1]

            # 1) Distribution
            qc_bins$idelta_at_grid <- suppressWarnings(as.numeric(qc_bins$idelta_at_grid))
            qc_ok <- qc_bins[is.finite(qc_bins$idelta_at_grid), , drop = FALSE]
            p_dist <- if (nrow(qc_ok) > 0) {
                ggplot2::ggplot(qc_ok, ggplot2::aes(x = idelta_at_grid)) +
                    ggplot2::geom_histogram(
                        ggplot2::aes(y = ggplot2::after_stat(density)),
                        bins = 40,
                        fill = "white",
                        colour = "black",
                        linewidth = .3
                    ) +
                    ggplot2::geom_density(color = blue_color, linewidth = .5) +
                    ggplot2::labs(
                        title = sprintf("Idelta at chosen grid (%.0f umm)", qc_grid),
                        x = "Idelta at grid",
                        y = "Density"
                    ) +
                    theme_lqc_qc()
            } else {
                ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0, y = 0, label = "No finite Idelta at grid")
            }
            p1 <- file.path(plots_out_dir, "scope_gridQC_idelta_distribution.png")
            ggplot2::ggsave(filename = p1, plot = p_dist, width = 3, height = 2.5, dpi = 600)
            plot_paths <- c(plot_paths, p1)

            # 2) Tiers
            tier_levels <- paste0("Q", seq_len(q_bins_out))
            qc_bins$tier_q5 <- factor(as.character(qc_bins$tier_q5), levels = tier_levels, ordered = TRUE)
            tier_tbl <- as.data.frame(table(qc_bins$tier_q5), stringsAsFactors = FALSE)
            colnames(tier_tbl) <- c("tier", "n")
            p_tiers <- if (nrow(tier_tbl) > 0) {
                ggplot2::ggplot(tier_tbl, ggplot2::aes(x = tier, y = n)) +
                    ggplot2::geom_col(fill = "white", colour = "black", linewidth = .3) +
                    ggplot2::labs(
                        title = sprintf("Idelta tiers at grid (%.0f umm)", qc_grid),
                        x = "Tier (quantile bins)",
                        y = "Genes"
                    ) +
                    theme_lqc_qc()
            } else {
                ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0, y = 0, label = "No tier data")
            }
            p2 <- file.path(plots_out_dir, "scope_gridQC_idelta_tiers.png")
            ggplot2::ggsave(filename = p2, plot = p_tiers, width = 3, height = 2.5, dpi = 600)
            plot_paths <- c(plot_paths, p2)

            # 3) Unsuitable flags (not mutually exclusive)
            n_qc <- nrow(qc_bins)
            flags_df <- data.frame(
                flag = c("near_zero", "scale_mismatch", "loss", "any"),
                n = c(
                    sum(qc_bins$flag_near_zero %in% TRUE, na.rm = TRUE),
                    sum(qc_bins$flag_scale_mismatch %in% TRUE, na.rm = TRUE),
                    sum(qc_bins$flag_loss %in% TRUE, na.rm = TRUE),
                    sum(qc_bins$unsuitable_any %in% TRUE, na.rm = TRUE)
                ),
                stringsAsFactors = FALSE
            )
            flags_df$frac <- if (n_qc > 0) flags_df$n / n_qc else NA_real_
            if (!requireNamespace("scales", quietly = TRUE)) {
                if (isTRUE(strict)) warning("Missing scales; using numeric y-axis for flag fractions.")
                p_flags <- ggplot2::ggplot(flags_df, ggplot2::aes(x = flag, y = frac)) +
                    ggplot2::geom_col(fill = "white", colour = "black", linewidth = .3) +
                    ggplot2::scale_y_continuous(limits = c(0, 1)) +
                    ggplot2::labs(
                        title = sprintf("Unsuitable flags at grid (%.0f umm)", qc_grid),
                        x = NULL,
                        y = "Fraction of informative genes"
                    ) +
                    theme_lqc_qc()
            } else {
                p_flags <- ggplot2::ggplot(flags_df, ggplot2::aes(x = flag, y = frac)) +
                    ggplot2::geom_col(fill = "white", colour = "black", linewidth = .3) +
                    ggplot2::scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
                    ggplot2::labs(
                        title = sprintf("Unsuitable flags at grid (%.0f umm)", qc_grid),
                        x = NULL,
                        y = "Fraction of informative genes"
                    ) +
                    theme_lqc_qc()
            }
            p3 <- file.path(plots_out_dir, "scope_gridQC_unsuitable_flags.png")
            ggplot2::ggsave(filename = p3, plot = p_flags, width = 3, height = 2.5, dpi = 600)
            plot_paths <- c(plot_paths, p3)

            # 4) Knee vs Idelta scatter
            qc_scatter <- qc_bins
            qc_scatter$knee_um <- suppressWarnings(as.numeric(qc_scatter$knee_um))
            qc_scatter <- qc_scatter[is.finite(qc_scatter$knee_um) & is.finite(qc_scatter$idelta_at_grid), , drop = FALSE]
            p_scatter <- if (nrow(qc_scatter) > 0) {
                qc_scatter$unsuitable_any <- qc_scatter$unsuitable_any %in% TRUE
                ggplot2::ggplot(qc_scatter, ggplot2::aes(x = knee_um, y = idelta_at_grid)) +
                    ggplot2::geom_point(ggplot2::aes(color = unsuitable_any), size = 1.0, alpha = 0.7) +
                    ggplot2::scale_color_manual(name = "unsuitable_any", values = c("FALSE" = "black", "TRUE" = red_color)) +
                    ggplot2::labs(
                        title = "Knee vs Idelta at chosen grid",
                        x = "knee_um (per-gene)",
                        y = "Idelta at chosen grid"
                    ) +
                    theme_lqc_qc()
            } else {
                ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::annotate("text", x = 0, y = 0, label = "No knee/Idelta pairs")
            }
            p4 <- file.path(plots_out_dir, "scope_gridQC_knee_vs_idelta.png")
            ggplot2::ggsave(filename = p4, plot = p_scatter, width = 6, height = 5, dpi = 600)
            plot_paths <- c(plot_paths, p4)
        }
    }

    list(
        recommended_grid_um_runner = recommended_grid_um_runner,
        qc_summary = qc$qc_summary,
        gene_bins = qc$gene_bins,
        qc_summary_path = qc_summary_path,
        gene_bins_path = gene_bins_path,
        plot_paths = plot_paths,
        tool_out_dir = tool_out_dir,
        plots_out_dir = plots_out_dir
    )
}

# Backward-compatible alias (legacy API name).
computeGridIDeltaQCExternal <- function(...) gridSelectQC(...)
