#' Check style-B variance assumptions on a real weight matrix.
#'
#' @description
#' Computes the analytic style-B variance term for a supplied spatial weight
#' matrix and compares it with a permutation-estimated null variance for one
#' observed vector pair, or for a bounded set of feature pairs from a real
#' expression matrix supplied through `X`. This is a diagnostic helper for assessing whether the
#' large-sample, symmetric-binary, zero-diagonal style-B assumptions are
#' reasonable for a specific graph. It is not used by the main `computeL()`
#' inference path.
#'
#' @param W Spatial weight matrix.
#' @param x Numeric observed vector for the first feature. Optional when `X`
#'   is supplied.
#' @param y Numeric observed vector for the second feature. Defaults to `x`.
#' @param X Optional feature matrix with rows matching `nrow(W)` and columns as
#'   genes/features. If supplied, a bounded set of feature pairs is checked.
#' @param max_pairs Maximum number of feature pairs to check when `X` is
#'   supplied.
#' @param n_perm Number of permutations used for the null variance estimate.
#' @param seed Optional random seed.
#' @param center_scale Center and scale `x` and `y` before computing Lee's L.
#' @param warn_ratio_threshold Warn when a permutation/style-B variance ratio
#'   is outside `[1 / warn_ratio_threshold, warn_ratio_threshold]`.
#' @return A list with analytic variance, permutation variance, ratio, observed
#'   L, permutation summary, and assumption flags.
#' @examples
#' \dontrun{
#' out <- styleB_check_realW(scope_obj@grid$grid30$W, x, y, n_perm = 499)
#' out$variance_ratio_perm_over_styleB
#' matrix_out <- styleB_check_realW(scope_obj@grid$grid30$W, X = expr, n_perm = 199)
#' matrix_out$summary
#' }
#' @export
styleB_check_realW <- function(W,
                               x = NULL,
                               y = x,
                               X = NULL,
                               max_pairs = 25L,
                               n_perm = 199L,
                               seed = NULL,
                               center_scale = TRUE,
                               warn_ratio_threshold = 2) {
    if (!inherits(W, "Matrix")) {
        W <- Matrix::Matrix(W, sparse = TRUE)
    }
    W <- suppressWarnings(methods::as(W, "dgCMatrix"))
    n <- nrow(W)
    if (ncol(W) != n) {
        stop("W must be square.", call. = FALSE)
    }
    if (any(!is.finite(W@x))) {
        stop("W must contain only finite values.", call. = FALSE)
    }
    n_perm <- as.integer(n_perm)[1]
    if (!is.finite(n_perm) || n_perm < 2L) {
        stop("n_perm must be an integer >= 2.", call. = FALSE)
    }
    max_pairs <- as.integer(max_pairs)[1]
    if (!is.finite(max_pairs) || max_pairs < 1L) {
        stop("max_pairs must be an integer >= 1.", call. = FALSE)
    }
    warn_ratio_threshold <- as.numeric(warn_ratio_threshold)[1]
    if (!is.finite(warn_ratio_threshold) || warn_ratio_threshold <= 1) {
        stop("warn_ratio_threshold must be > 1.", call. = FALSE)
    }
    if (!is.null(seed)) {
        old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        } else {
            NULL
        }
        on.exit({
            if (is.null(old_seed)) {
                if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
                    rm(".Random.seed", envir = .GlobalEnv)
                }
            } else {
                assign(".Random.seed", old_seed, envir = .GlobalEnv)
            }
        }, add = TRUE)
        set.seed(seed)
    }

    z <- function(v) {
        if (!isTRUE(center_scale)) return(v)
        vv <- v - mean(v)
        s <- stats::sd(vv)
        if (!is.finite(s) || s <= 0) return(rep(0, length(v)))
        vv / s
    }
    lee_one <- function(a, b) {
        s0 <- sum(W)
        denom <- sqrt(sum(a^2) * sum(b^2))
        if (!is.finite(s0) || s0 == 0 || !is.finite(denom) || denom == 0) {
            return(NA_real_)
        }
        out <- as.numeric((n / s0) * (Matrix::t(a) %*% (W %*% b)) / denom)
        if (is.finite(out)) out else NA_real_
    }

    S0 <- sum(W)
    S1 <- suppressWarnings(0.5 * sum((W + Matrix::t(W))^2))
    S2 <- sum((Matrix::rowSums(W) + Matrix::colSums(W))^2)
    EZ <- -1 / (n - 1)
    var_B <- (n^2 * S1 - n * S2 + 3 * S0^2) / (S0^2 * (n^2 - 1)) - EZ^2
    diag_vals <- Matrix::diag(W)
    symmetric <- isTRUE(suppressWarnings(Matrix::isSymmetric(W)))
    binary <- length(W@x) == 0L || all(W@x %in% c(0, 1))
    zero_diag <- all(!is.finite(diag_vals) | diag_vals == 0)
    large_n <- n >= 100L
    assumptions_hold <- symmetric && binary && zero_diag && large_n

    check_pair <- function(x_vec, y_vec, pair_label = NA_character_) {
        x_vec <- as.numeric(x_vec)
        y_vec <- as.numeric(y_vec)
        if (length(x_vec) != n || length(y_vec) != n) {
            stop("x and y must have length nrow(W).", call. = FALSE)
        }
        if (any(!is.finite(x_vec)) || any(!is.finite(y_vec))) {
            stop("x, y, and X must contain only finite values.", call. = FALSE)
        }
        x_vec <- z(x_vec)
        y_vec <- z(y_vec)
        L_obs <- lee_one(x_vec, y_vec)
        perm_L <- numeric(n_perm)
        for (i in seq_len(n_perm)) {
            perm_L[[i]] <- lee_one(x_vec, sample(y_vec, length(y_vec), replace = FALSE))
        }
        var_perm <- stats::var(perm_L, na.rm = TRUE)
        p_two_sided <- (sum(abs(perm_L) >= abs(L_obs), na.rm = TRUE) + 1) / (sum(is.finite(perm_L)) + 1)
        ratio <- var_perm / var_B
        data.frame(
            pair = pair_label,
            L_observed = L_obs,
            Var_B = var_B,
            Var_from_permutation = var_perm,
            variance_ratio_perm_over_styleB = ratio,
            permutation_p_two_sided = p_two_sided,
            permutation_L_mean = mean(perm_L, na.rm = TRUE),
            permutation_L_sd = stats::sd(perm_L, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    }

    if (!assumptions_hold) {
        warning(
            "styleB_check_realW: style-B large-sample 2x interpretation is boundary-limited; inspect assumption flags.",
            call. = FALSE
        )
    }

    assumption_list <- list(
        symmetric = symmetric,
        binary = binary,
        zero_diagonal = zero_diag,
        large_n = large_n,
        applies_styleB_large_sample_rule = assumptions_hold
    )

    if (!is.null(X)) {
        if (!inherits(X, "Matrix")) {
            X <- as.matrix(X)
        }
        if (nrow(X) != n && ncol(X) == n) {
            X <- Matrix::t(X)
        }
        if (nrow(X) != n) {
            stop("X must have nrow(X) == nrow(W), or ncol(X) == nrow(W) for transposable input.", call. = FALSE)
        }
        X <- as.matrix(X)
        storage.mode(X) <- "double"
        finite_cols <- which(colSums(is.finite(X)) == n)
        if (length(finite_cols) < 1L) {
            stop("X must contain at least one fully finite feature column.", call. = FALSE)
        }
        X <- X[, finite_cols, drop = FALSE]
        feature_names <- colnames(X) %||% paste0("feature", seq_len(ncol(X)))
        pairs <- utils::combn(seq_len(ncol(X)), 2L)
        if (!ncol(pairs)) {
            pairs <- rbind(1L, 1L)
        }
        pairs <- pairs[, seq_len(min(base::ncol(pairs), max_pairs)), drop = FALSE]
        pair_results <- vector("list", base::ncol(pairs))
        for (i in seq_len(base::ncol(pairs))) {
            pair_results[[i]] <- check_pair(
                X[, pairs[1L, i]],
                X[, pairs[2L, i]],
                paste(feature_names[pairs[1L, i]], feature_names[pairs[2L, i]], sep = "__")
            )
        }
        pair_results <- do.call(rbind, pair_results)
        ratios <- pair_results$variance_ratio_perm_over_styleB
        out <- list(
            mode = "matrix",
            n = n,
            n_features_checked = length(unique(as.integer(pairs))),
            n_pairs_checked = nrow(pair_results),
            n_perm = n_perm,
            Var_B = var_B,
            pair_results = pair_results,
            summary = data.frame(
                ratio_median = stats::median(ratios, na.rm = TRUE),
                ratio_min = min(ratios, na.rm = TRUE),
                ratio_max = max(ratios, na.rm = TRUE),
                ratio_outside_threshold_fraction = mean(
                    ratios < (1 / warn_ratio_threshold) | ratios > warn_ratio_threshold,
                    na.rm = TRUE
                ),
                stringsAsFactors = FALSE
            ),
            assumptions = assumption_list
        )
        if (isTRUE(out$summary$ratio_outside_threshold_fraction > 0)) {
            warning(
                "styleB_check_realW: empirical permutation variance differs materially from analytic style-B variance for at least one checked feature pair.",
                call. = FALSE
            )
        }
        return(out)
    }

    if (is.null(x)) {
        stop("Either x/y or X must be supplied.", call. = FALSE)
    }
    pair <- check_pair(x, y, "x__y")
    ratio <- pair$variance_ratio_perm_over_styleB[[1]]
    if (is.finite(ratio) && (ratio < (1 / warn_ratio_threshold) || ratio > warn_ratio_threshold)) {
        warning(
            "styleB_check_realW: empirical permutation variance differs materially from analytic style-B variance.",
            call. = FALSE
        )
    }

    list(
        n = n,
        n_perm = n_perm,
        L_observed = pair$L_observed[[1]],
        Var_B = var_B,
        Var_from_permutation = pair$Var_from_permutation[[1]],
        variance_ratio_perm_over_styleB = ratio,
        permutation_p_two_sided = pair$permutation_p_two_sided[[1]],
        permutation_L_mean = pair$permutation_L_mean[[1]],
        permutation_L_sd = pair$permutation_L_sd[[1]],
        assumptions = assumption_list
    )
}
