#' Default Trust Mode
#' @description
#' Normalizes trust mode values for signature verification.
#' @param trust_mode Optional mode override.
#' @param default Default mode when `trust_mode` is `NULL`.
#' @param supported Supported modes.
#' @return Scalar character string.
#' @keywords internal
.trust_chain_default_mode <- function(trust_mode = NULL,
                                      default = "ed25519",
                                      supported = c("ed25519")) {
    if (is.null(trust_mode) || !length(trust_mode) || is.na(trust_mode[[1]]) || !nzchar(trust_mode[[1]])) {
        return(default)
    }
    trust_mode <- tolower(as.character(trust_mode)[1])
    if (!trust_mode %in% supported) {
        stop("Unsupported trust mode: ", trust_mode, call. = FALSE)
    }
    trust_mode
}

#' Resolve Trust Root Path
#' @description
#' Resolves an explicit or configured trust-root path.
#' @param trust_root_path Optional direct path.
#' @param option_name Optional option name.
#' @param envvar Optional environment variable name.
#' @param default Default fallback path.
#' @return Normalized path or `NULL`.
#' @keywords internal
.trust_chain_resolve_root <- function(trust_root_path = NULL,
                                      option_name = NULL,
                                      envvar = NULL,
                                      default = NULL) {
    path <- trust_root_path
    if (is.null(path) && !is.null(option_name) && nzchar(option_name)) {
        path <- getOption(option_name, NULL)
    }
    if (is.null(path) && !is.null(envvar) && nzchar(envvar)) {
        env_path <- Sys.getenv(envvar, unset = "")
        if (nzchar(env_path)) path <- env_path
    }
    if (is.null(path) && !is.null(default) && nzchar(default)) {
        path <- default
    }
    if (is.null(path) || !length(path)) return(NULL)
    path <- as.character(path)[1]
    if (is.na(path) || !nzchar(path)) return(NULL)
    normalizePath(path.expand(path), winslash = "/", mustWork = FALSE)
}

#' Require OpenSSL R Package
#' @description
#' Ed25519 verification depends on the `openssl` R package.
#' @return Invisible TRUE.
#' @keywords internal
.trust_chain_require_openssl <- function() {
    if (!requireNamespace("openssl", quietly = TRUE)) {
        stop("Package 'openssl' is required for Ed25519 trust verification.", call. = FALSE)
    }
    invisible(TRUE)
}

#' Decode Base64 Signature File
#' @description
#' Reads a base64-encoded signature file and returns raw bytes.
#' @param sig_path Signature file path.
#' @return Raw vector.
#' @keywords internal
.decode_trust_signature_file <- function(sig_path) {
    if (is.null(sig_path) || !nzchar(sig_path) || !file.exists(sig_path)) {
        stop("Signature file not found: ", sig_path, call. = FALSE)
    }
    .trust_chain_require_openssl()
    sig_text <- paste(readLines(sig_path, warn = FALSE, encoding = "UTF-8"), collapse = "")
    if (!nzchar(sig_text)) {
        stop("Signature payload is empty: ", sig_path, call. = FALSE)
    }
    sig_raw <- tryCatch(openssl::base64_decode(sig_text), error = function(e) e)
    if (inherits(sig_raw, "error") || !is.raw(sig_raw) || !length(sig_raw)) {
        stop("Signature payload could not be decoded: ", sig_path, call. = FALSE)
    }
    sig_raw
}

#' Read Ed25519 Trust Root
#' @description
#' Loads an Ed25519 public key from a trust-root path or bundled PEM text.
#' @param trust_root_path Optional trust-root path.
#' @param bundled_public_key_pem Optional bundled PEM string.
#' @return List with loaded key and path metadata.
#' @keywords internal
.read_ed25519_trust_root <- function(trust_root_path = NULL,
                                     bundled_public_key_pem = NULL) {
    .trust_chain_require_openssl()
    key_path <- trust_root_path
    key_tmp <- FALSE
    if (is.null(key_path) || !nzchar(key_path)) {
        if (is.null(bundled_public_key_pem) || !nzchar(bundled_public_key_pem)) {
            stop("No trust root path or bundled Ed25519 public key is available.", call. = FALSE)
        }
        key_path <- tempfile("genescope_trust_root_", fileext = ".pem")
        writeLines(bundled_public_key_pem, key_path, useBytes = TRUE)
        key_tmp <- TRUE
    }
    if (!file.exists(key_path)) {
        stop("Trust root not found: ", key_path, call. = FALSE)
    }
    key_obj <- tryCatch(openssl::read_pubkey(key_path), error = function(e) e)
    if (inherits(key_obj, "error")) {
        stop("Failed to read Ed25519 trust root: ", conditionMessage(key_obj), call. = FALSE)
    }
    list(
        key = key_obj,
        path = normalizePath(key_path, winslash = "/", mustWork = FALSE),
        temporary = key_tmp
    )
}

#' Verify Ed25519-Signed File
#' @description
#' Verifies a base64-encoded Ed25519 signature against a file on disk.
#' @param data_path Path to the signed file.
#' @param sig_path Signature file path.
#' @param trust_mode Trust mode (currently `ed25519` only).
#' @param trust_root_path Optional trust-root path.
#' @param bundled_public_key_pem Optional bundled PEM string.
#' @param key_id Key identifier to record in the result.
#' @return Verification result list.
#' @keywords internal
.verify_ed25519_file_signature <- function(data_path,
                                           sig_path,
                                           trust_mode = "ed25519",
                                           trust_root_path = NULL,
                                           bundled_public_key_pem = NULL,
                                           key_id = NA_character_) {
    trust_mode <- .trust_chain_default_mode(trust_mode)
    if (!identical(trust_mode, "ed25519")) {
        stop("Unsupported trust mode for file verification: ", trust_mode, call. = FALSE)
    }
    if (is.null(data_path) || !nzchar(data_path) || !file.exists(data_path)) {
        return(list(ok = FALSE, reason = "missing_signed_file", data_path = data_path))
    }
    if (is.null(sig_path) || !nzchar(sig_path) || !file.exists(sig_path)) {
        return(list(ok = FALSE, reason = "missing_signature", data_path = data_path))
    }

    root <- tryCatch(
        .read_ed25519_trust_root(
            trust_root_path = trust_root_path,
            bundled_public_key_pem = bundled_public_key_pem
        ),
        error = function(e) e
    )
    if (inherits(root, "error")) {
        return(list(ok = FALSE, reason = conditionMessage(root), data_path = data_path))
    }
    on.exit({
        if (isTRUE(root$temporary) && file.exists(root$path)) unlink(root$path)
    }, add = TRUE)

    signature <- tryCatch(.decode_trust_signature_file(sig_path), error = function(e) e)
    if (inherits(signature, "error")) {
        return(list(ok = FALSE, reason = conditionMessage(signature), data_path = data_path))
    }

    payload <- readBin(data_path, what = "raw", n = file.info(data_path)$size[[1]])
    verified <- tryCatch(openssl::ed25519_verify(payload, signature, root$key), error = function(e) e)
    if (inherits(verified, "error")) {
        return(list(ok = FALSE, reason = conditionMessage(verified), data_path = data_path))
    }
    if (!isTRUE(verified)) {
        return(list(ok = FALSE, reason = "signature_verification_failed", data_path = data_path))
    }

    list(
        ok = TRUE,
        reason = "ok",
        data_path = normalizePath(data_path, winslash = "/", mustWork = FALSE),
        signature_path = normalizePath(sig_path, winslash = "/", mustWork = FALSE),
        signature_sha256 = .compute_file_sha256(sig_path),
        trust_root_path = root$path,
        trust_root_sha256 = .compute_file_sha256(root$path),
        signature_key_id = as.character(key_id)[1],
        trust_mode = trust_mode
    )
}
