#' Normalize support percentage values to the package canonical scale.
#' @description
#' Internal helper for APIs that accept support percentage values. geneSCOPE
#' now keeps a canonical internal support scale on `[0,1]` for gating and
#' diagnostics, while still returning legacy percent-scale views (`0-100`) for
#' backward-compatible table columns and reports.
#' @param values Numeric vector of support percentage values.
#' @param support_pct_scale One of `auto`, `0-1`, or `0-100`.
#' @param caller Label used in warnings.
#' @return List containing canonical `values`/`values_prop`, legacy
#'   `values_pct`, and scale diagnostics.
#' @keywords internal
.normalize_support_pct_values <- function(values,
                                          support_pct_scale = c("auto", "0-1", "0-100"),
                                          caller = ".normalize_support_pct_values") {
    support_pct_scale <- match.arg(support_pct_scale)
    values <- as.numeric(values)

    out <- list(
        values = values,
        values_prop = values,
        values_pct = values,
        requested_scale = support_pct_scale,
        input_scale = NA_character_,
        scale_used = "missing",
        normalized_scale = "missing",
        internal_scale = "0-1",
        legacy_scale = "0-100",
        auto_normalized = FALSE,
        invalid_values = FALSE
    )

    observed <- values[!is.na(values)]
    if (!length(observed)) {
        return(out)
    }

    has_negative <- any(observed < 0)
    has_gt_100 <- any(observed > 100 + 1e-8)
    invalid_values <- has_negative || has_gt_100

    detected_scale <- if (has_negative && has_gt_100) {
        "invalid-negative-and->100"
    } else if (has_negative) {
        "invalid-negative"
    } else if (has_gt_100) {
        "invalid->100"
    } else if (all(observed <= 1 + 1e-8)) {
        "0-1"
    } else if (all(observed <= 100 + 1e-8)) {
        "0-100"
    } else {
        "invalid"
    }

    scale_used <- "0-1"
    normalized_scale <- "0-1"
    auto_normalized <- FALSE
    values_prop <- values
    values_pct <- values

    if (invalid_values) {
        warning(
            caller, ": support_pct contains values outside the valid domain; keeping values unchanged and recording the scale as invalid.",
            call. = FALSE
        )
        scale_used <- "unchanged-invalid"
        normalized_scale <- "invalid"
    } else {
        declared_scale <- if (identical(support_pct_scale, "auto")) detected_scale else support_pct_scale
        if (is.na(declared_scale) || !declared_scale %in% c("0-1", "0-100")) {
            declared_scale <- detected_scale
        }
        auto_normalized <- identical(support_pct_scale, "auto") &&
            identical(detected_scale, "0-100")

        if (identical(declared_scale, "0-100")) {
            values_pct <- values
            values_prop <- values / 100
            scale_used <- "0-100->0-1"
            if (identical(support_pct_scale, "auto")) {
                warning(
                    caller, ": support_pct auto-detected 0-100 input and normalized it to the canonical [0,1] scale.",
                    call. = FALSE
                )
            }
        } else if (identical(detected_scale, "0-100")) {
            warning(
                caller, ": support_pct_scale='0-1' but observed values look like 0-100; normalizing them to the canonical [0,1] scale.",
                call. = FALSE
            )
            values_prop <- values / 100
            values_pct <- values
            scale_used <- "0-100->0-1"
        } else {
            values_prop <- values
            values_pct <- values * 100
            scale_used <- "0-1"
        }
    }

    out$values <- values_prop
    out$values_prop <- values_prop
    out$values_pct <- values_pct
    out$input_scale <- detected_scale
    out$scale_used <- scale_used
    out$normalized_scale <- normalized_scale
    out$auto_normalized <- auto_normalized
    out$invalid_values <- invalid_values
    out
}

#' Normalize support threshold values to the canonical `[0,1]` scale.
#' @description
#' Internal helper for threshold arguments such as `min_support_pct`.
#' @param value Numeric scalar threshold.
#' @param support_pct_scale One of `auto`, `0-1`, or `0-100`.
#' @param detected_scale Optional detected input scale from the data path.
#' @param caller Label used in warnings.
#' @return Numeric scalar on the canonical `[0,1]` support scale.
#' @keywords internal
.normalize_support_pct_threshold <- function(value,
                                             support_pct_scale = c("auto", "0-1", "0-100"),
                                             detected_scale = NA_character_,
                                             caller = ".normalize_support_pct_threshold") {
    support_pct_scale <- match.arg(support_pct_scale)
    value <- suppressWarnings(as.numeric(value)[1])
    if (!is.finite(value) || value <= 0) {
        return(0)
    }

    declared_scale <- support_pct_scale
    if (identical(declared_scale, "auto")) {
        if (!is.na(detected_scale) && detected_scale %in% c("0-1", "0-100")) {
            declared_scale <- detected_scale
        } else if (value <= 1 + 1e-8) {
            declared_scale <- "0-1"
        } else {
            declared_scale <- "0-100"
        }
    }

    if (identical(declared_scale, "0-100")) {
        return(value / 100)
    }

    if (value > 1 + 1e-8) {
        warning(
            caller, ": support threshold declared as 0-1 but value exceeds 1; interpreting it as percent-like input and normalizing to [0,1].",
            call. = FALSE
        )
        return(value / 100)
    }

    value
}
