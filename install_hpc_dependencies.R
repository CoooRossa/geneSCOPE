#!/usr/bin/env Rscript
# HPC-Optimized Dependency Installation for geneSCOPE
# This script handles problematic dependencies in HPC environments

cat("[geneSCOPE] Starting HPC-optimized dependency installation...\n")

# Function to install packages with fallback strategies
install_with_fallback <- function(packages, method = "auto") {
    success <- character(0)
    failed <- character(0)

    for (pkg in packages) {
        cat(sprintf("[geneSCOPE] Installing %s...\n", pkg))

        # Try multiple installation strategies
        installed <- FALSE

        # Strategy 1: Standard CRAN
        if (!installed) {
            tryCatch(
                {
                    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
                        install.packages(pkg,
                            repos = "https://cloud.r-project.org/",
                            dependencies = FALSE, INSTALL_opts = "--no-docs"
                        )
                        if (require(pkg, character.only = TRUE, quietly = TRUE)) {
                            installed <- TRUE
                            cat(sprintf("[geneSCOPE] ✓ %s installed successfully\n", pkg))
                        }
                    } else {
                        installed <- TRUE
                        cat(sprintf("[geneSCOPE] ✓ %s already available\n", pkg))
                    }
                },
                error = function(e) {
                    cat(sprintf("[geneSCOPE] Strategy 1 failed for %s: %s\n", pkg, e$message))
                }
            )
        }

        # Strategy 2: Binary installation (if available)
        if (!installed) {
            tryCatch(
                {
                    install.packages(pkg,
                        repos = "https://cloud.r-project.org/",
                        type = "binary", dependencies = FALSE
                    )
                    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
                        installed <- TRUE
                        cat(sprintf("[geneSCOPE] ✓ %s installed from binary\n", pkg))
                    }
                },
                error = function(e) {
                    cat(sprintf("[geneSCOPE] Binary installation failed for %s: %s\n", pkg, e$message))
                }
            )
        }

        # Strategy 3: Install from source with minimal configuration
        if (!installed) {
            tryCatch(
                {
                    install.packages(pkg,
                        repos = "https://cloud.r-project.org/",
                        type = "source", dependencies = FALSE,
                        configure.args = "--disable-shared",
                        INSTALL_opts = c("--no-docs", "--no-html", "--no-demo")
                    )
                    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
                        installed <- TRUE
                        cat(sprintf("[geneSCOPE] ✓ %s installed from source\n", pkg))
                    }
                },
                error = function(e) {
                    cat(sprintf("[geneSCOPE] Source installation failed for %s: %s\n", pkg, e$message))
                }
            )
        }

        if (installed) {
            success <- c(success, pkg)
        } else {
            failed <- c(failed, pkg)
            cat(sprintf("[geneSCOPE] !!! Failed to install %s !!!\n", pkg))
        }
    }

    return(list(success = success, failed = failed))
}

# Function to disable httpgd-related options
disable_problematic_options <- function() {
    cat("[geneSCOPE] Disabling problematic graphics options...\n")

    # Disable httpgd if it's causing issues
    options(device = "png")
    options(bitmapType = "cairo")

    # Set conservative CRAN mirror
    options(repos = c(CRAN = "https://cloud.r-project.org/"))

    # Disable parallel installation that might cause issues
    options(Ncpus = 1)

    cat("[geneSCOPE] Graphics options configured for HPC\n")
}

# Essential packages in dependency order
essential_packages <- c(
    # Base dependencies
    "Rcpp", "RcppArmadillo", "RcppEigen",

    # Core utilities
    "data.table", "dplyr", "tidyr", "scales",

    # Graphics and plotting (minimal set)
    "ggplot2", "ggrepel", "viridis",

    # Spatial analysis
    "sp", "sf", "spdep",

    # Network analysis
    "igraph",

    # Parallel processing
    "foreach", "future.apply"
)

# Optional packages (install if possible, skip if not)
optional_packages <- c(
    "ggforce", "ggnewscale", "ggraph", "arrow", "RhpcBLASctl"
)

# Main installation process
main <- function() {
    cat("[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] HPC Dependency Installation Started\n")
    cat("[geneSCOPE] ========================================\n")

    # Configure environment
    disable_problematic_options()

    # Check R version
    r_version <- R.Version()
    cat(sprintf(
        "[geneSCOPE] R version: %s.%s.%s\n",
        r_version$major, r_version$minor, r_version$year
    ))

    # Install essential packages
    cat("[geneSCOPE] Installing essential packages...\n")
    essential_result <- install_with_fallback(essential_packages)

    # Install optional packages
    cat("[geneSCOPE] Installing optional packages...\n")
    optional_result <- install_with_fallback(optional_packages)

    # Summary
    total_success <- length(c(essential_result$success, optional_result$success))
    total_failed <- length(c(essential_result$failed, optional_result$failed))

    cat("[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] Installation Summary\n")
    cat("[geneSCOPE] ========================================\n")
    cat(sprintf("[geneSCOPE] Successfully installed: %d packages\n", total_success))
    cat(sprintf("[geneSCOPE] Failed to install: %d packages\n", total_failed))

    if (length(essential_result$failed) > 0) {
        cat("[geneSCOPE] !!! Critical failures (essential packages) !!!\n")
        cat(paste("[geneSCOPE]", essential_result$failed, collapse = "\n"))
        cat("\n")
    }

    if (length(optional_result$failed) > 0) {
        cat("[geneSCOPE] Optional package failures (can be ignored):\n")
        cat(paste("[geneSCOPE]", optional_result$failed, collapse = "\n"))
        cat("\n")
    }

    # Check if we can proceed
    critical_missing <- essential_result$failed
    if (length(critical_missing) > 0) {
        cat("[geneSCOPE] !!! Cannot proceed - critical dependencies missing !!!\n")
        cat("[geneSCOPE] Try the following manual installation commands:\n")
        for (pkg in critical_missing) {
            cat(sprintf("install.packages('%s', repos = 'https://cloud.r-project.org/', dependencies = FALSE)\n", pkg))
        }
        return(FALSE)
    } else {
        cat("[geneSCOPE] ✓ All essential dependencies installed successfully!\n")
        cat("[geneSCOPE] You can now install geneSCOPE with:\n")
        cat("[geneSCOPE] devtools::install_github('CoooRossa/geneSCOPE', dependencies = FALSE)\n")
        return(TRUE)
    }
}

# Run if called as script
if (!interactive()) {
    success <- main()
    if (!success) {
        quit(status = 1)
    }
} else {
    cat("[geneSCOPE] HPC dependency installer loaded. Run main() to start installation.\n")
}
