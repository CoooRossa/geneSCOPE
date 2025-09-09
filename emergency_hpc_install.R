#!/usr/bin/env Rscript
# Emergency HPC Installer - Bypasses httpgd issues
# This script specifically addresses the httpgd dependency problem

cat("[geneSCOPE] Emergency HPC Installer Started\n")
cat("[geneSCOPE] Targeting httpgd dependency issues\n")

# Disable all graphics-related problematic packages
disable_graphics_dependencies <- function() {
    cat("[geneSCOPE] Disabling graphics dependencies...\n")

    # Set environment to avoid httpgd
    Sys.setenv("R_BROWSER" = "false")
    Sys.setenv("DISPLAY" = "")

    # Configure safe graphics options
    options(device = "png")
    options(bitmapType = "cairo")
    options(menu.graphics = FALSE)

    # Disable Suggests dependencies that might pull in httpgd
    options(install.packages.check.source = "no")
    options(install.packages.compile.from.source = "never")

    cat("[geneSCOPE] Graphics dependencies disabled\n")
}

# Install packages without Suggests dependencies
install_without_suggests <- function(packages) {
    cat("[geneSCOPE] Installing packages without Suggests dependencies...\n")

    success <- character(0)
    failed <- character(0)

    for (pkg in packages) {
        cat(sprintf("[geneSCOPE] Installing %s (no Suggests)...\n", pkg))

        tryCatch(
            {
                install.packages(pkg,
                    dependencies = c("Depends", "Imports"), # Exclude Suggests
                    repos = "https://cloud.r-project.org/",
                    type = "binary",
                    INSTALL_opts = c("--no-docs", "--no-html")
                )

                if (require(pkg, character.only = TRUE, quietly = TRUE)) {
                    success <- c(success, pkg)
                    cat(sprintf("[geneSCOPE] ✓ %s installed\n", pkg))
                } else {
                    failed <- c(failed, pkg)
                }
            },
            error = function(e) {
                cat(sprintf("[geneSCOPE] !!! %s failed: %s\n", pkg, e$message))
                failed <<- c(failed, pkg)
            }
        )
    }

    return(list(success = success, failed = failed))
}

# Core packages that definitely don't need httpgd
core_packages <- c(
    "Rcpp",
    "RcppArmadillo",
    "data.table",
    "dplyr",
    "Matrix",
    "foreach"
)

# Graphics packages (install carefully)
safe_graphics <- c(
    "ggplot2",
    "scales",
    "viridis"
)

# Spatial/network packages
spatial_packages <- c(
    "sp",
    "igraph"
)

# Main function
emergency_install <- function() {
    cat("[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] Emergency HPC Installation\n")
    cat("[geneSCOPE] Bypassing httpgd dependencies\n")
    cat("[geneSCOPE] ========================================\n")

    # Step 1: Configure safe environment
    disable_graphics_dependencies()

    # Step 2: Install core packages first
    cat("[geneSCOPE] Installing core packages...\n")
    core_result <- install_without_suggests(core_packages)

    # Step 3: Install graphics packages carefully
    cat("[geneSCOPE] Installing graphics packages...\n")
    graphics_result <- install_without_suggests(safe_graphics)

    # Step 4: Install spatial packages
    cat("[geneSCOPE] Installing spatial packages...\n")
    spatial_result <- install_without_suggests(spatial_packages)

    # Step 5: Try to install devtools minimally
    cat("[geneSCOPE] Installing devtools...\n")
    devtools_result <- install_without_suggests("devtools")

    # Summary
    all_success <- c(
        core_result$success, graphics_result$success,
        spatial_result$success, devtools_result$success
    )
    all_failed <- c(
        core_result$failed, graphics_result$failed,
        spatial_result$failed, devtools_result$failed
    )

    cat("[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] Emergency Installation Summary\n")
    cat("[geneSCOPE] ========================================\n")
    cat(sprintf("[geneSCOPE] Successful: %d packages\n", length(all_success)))
    cat(sprintf("[geneSCOPE] Failed: %d packages\n", length(all_failed)))

    if (length(all_failed) > 0) {
        cat("[geneSCOPE] Failed packages:\n")
        cat(paste("[geneSCOPE]", all_failed, collapse = "\n"))
        cat("\n")
    }

    # Check if we can proceed with geneSCOPE
    essential_check <- c("Rcpp", "data.table", "ggplot2")
    missing_essential <- setdiff(essential_check, all_success)

    if (length(missing_essential) == 0) {
        cat("[geneSCOPE] ✓ Essential packages available - can install geneSCOPE\n")

        # Try installing geneSCOPE
        if ("devtools" %in% all_success) {
            cat("[geneSCOPE] Attempting geneSCOPE installation...\n")
            tryCatch(
                {
                    devtools::install_github("haenolab/geneSCOPE",
                        dependencies = FALSE,
                        upgrade = "never"
                    )

                    if (require("geneSCOPE", quietly = TRUE)) {
                        cat("[geneSCOPE] ✓ geneSCOPE installed successfully!\n")
                        return(TRUE)
                    }
                },
                error = function(e) {
                    cat("[geneSCOPE] geneSCOPE installation failed:", e$message, "\n")
                }
            )
        }
    } else {
        cat("[geneSCOPE] !!! Cannot proceed - missing essential packages:\n")
        cat(paste("[geneSCOPE]", missing_essential, collapse = "\n"))
        cat("\n")
    }

    cat("[geneSCOPE] If this fails, try manual installation:\n")
    cat("[geneSCOPE] 1. Load R module: module load R\n")
    cat("[geneSCOPE] 2. Install packages one by one\n")
    cat("[geneSCOPE] 3. Use: install.packages('pkg', dependencies=FALSE)\n")

    return(FALSE)
}

# Run the emergency installer
if (!interactive()) {
    emergency_install()
} else {
    cat("[geneSCOPE] Emergency installer loaded. Run emergency_install() to start.\n")
}
