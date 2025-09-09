#!/usr/bin/env Rscript
# Simplified geneSCOPE Installation for HPC
# Use this after dependencies are installed with install_hpc_dependencies.R

cat("[geneSCOPE] Starting geneSCOPE installation...\n")

# Configure environment for HPC
configure_hpc_environment <- function() {
    cat("[geneSCOPE] Configuring HPC environment...\n")

    # Disable problematic graphics options
    options(device = "png")
    options(bitmapType = "cairo")
    Sys.setenv("R_BROWSER" = "false")

    # Set conservative installation options
    options(repos = c(CRAN = "https://cloud.r-project.org/"))
    options(timeout = 300)
    options(Ncpus = 1)

    # Disable documentation building during installation
    Sys.setenv("R_INSTALL_TAR" = "internal")

    cat("[geneSCOPE] Environment configured\n")
}

# Check if essential dependencies are available
check_dependencies <- function() {
    cat("[geneSCOPE] Checking essential dependencies...\n")

    essential <- c("Rcpp", "RcppArmadillo", "data.table", "ggplot2", "igraph")
    missing <- character(0)

    for (pkg in essential) {
        if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
            missing <- c(missing, pkg)
        }
    }

    if (length(missing) > 0) {
        cat("[geneSCOPE] !!! Missing essential dependencies !!!\n")
        cat(paste("[geneSCOPE] Missing:", paste(missing, collapse = ", ")), "\n")
        cat("[geneSCOPE] Please run install_hpc_dependencies.R first\n")
        return(FALSE)
    } else {
        cat("[geneSCOPE] ✓ All essential dependencies found\n")
        return(TRUE)
    }
}

# Install devtools if needed
install_devtools <- function() {
    cat("[geneSCOPE] Checking devtools...\n")

    if (!require("devtools", quietly = TRUE)) {
        cat("[geneSCOPE] Installing devtools...\n")
        tryCatch(
            {
                install.packages("devtools",
                    dependencies = FALSE,
                    repos = "https://cloud.r-project.org/"
                )
                if (require("devtools", quietly = TRUE)) {
                    cat("[geneSCOPE] ✓ devtools installed\n")
                    return(TRUE)
                } else {
                    cat("[geneSCOPE] !!! devtools installation failed !!!\n")
                    return(FALSE)
                }
            },
            error = function(e) {
                cat("[geneSCOPE] !!! devtools installation error:", e$message, "!!!\n")
                return(FALSE)
            }
        )
    } else {
        cat("[geneSCOPE] ✓ devtools already available\n")
        return(TRUE)
    }
}

# Install geneSCOPE from GitHub
install_genescope <- function() {
    cat("[geneSCOPE] Installing geneSCOPE from GitHub...\n")

    tryCatch(
        {
            # Install without dependencies (already handled)
            devtools::install_github("haenolab/geneSCOPE",
                dependencies = FALSE,
                upgrade = "never",
                build_opts = c("--no-resave-data", "--no-manual"),
                quiet = FALSE
            )

            # Test if installation succeeded
            if (require("geneSCOPE", quietly = TRUE)) {
                cat("[geneSCOPE] ✓ geneSCOPE installed successfully!\n")
                return(TRUE)
            } else {
                cat("[geneSCOPE] !!! geneSCOPE installation failed - package not loadable !!!\n")
                return(FALSE)
            }
        },
        error = function(e) {
            cat("[geneSCOPE] !!! geneSCOPE installation error:", e$message, "!!!\n")
            return(FALSE)
        }
    )
}

# Verify installation
verify_installation <- function() {
    cat("[geneSCOPE] Verifying installation...\n")

    tryCatch(
        {
            library(geneSCOPE)

            # Check exported functions
            exported_functions <- ls("package:geneSCOPE")
            cat(sprintf("[geneSCOPE] Found %d exported functions\n", length(exported_functions)))

            # Check version
            version <- packageVersion("geneSCOPE")
            cat(sprintf("[geneSCOPE] Package version: %s\n", version))

            # Test a simple function if available
            if ("generateSpatialNetwork" %in% exported_functions) {
                cat("[geneSCOPE] ✓ Core functions accessible\n")
            }

            return(TRUE)
        },
        error = function(e) {
            cat("[geneSCOPE] !!! Verification failed:", e$message, "!!!\n")
            return(FALSE)
        }
    )
}

# Main installation function
main <- function() {
    cat("[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] geneSCOPE HPC Installation Started\n")
    cat("[geneSCOPE] ========================================\n")

    # Step 1: Configure environment
    configure_hpc_environment()

    # Step 2: Check dependencies
    if (!check_dependencies()) {
        cat("[geneSCOPE] !!! Installation aborted - dependencies missing !!!\n")
        return(FALSE)
    }

    # Step 3: Install devtools
    if (!install_devtools()) {
        cat("[geneSCOPE] !!! Installation aborted - devtools failed !!!\n")
        return(FALSE)
    }

    # Step 4: Install geneSCOPE
    if (!install_genescope()) {
        cat("[geneSCOPE] !!! Installation aborted - geneSCOPE installation failed !!!\n")
        return(FALSE)
    }

    # Step 5: Verify installation
    if (!verify_installation()) {
        cat("[geneSCOPE] !!! Installation verification failed !!!\n")
        return(FALSE)
    }

    # Success message
    cat("[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] Installation Completed Successfully!\n")
    cat("[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] You can now use geneSCOPE with:\n")
    cat("[geneSCOPE] library(geneSCOPE)\n")
    cat("[geneSCOPE] help(package = 'geneSCOPE')\n")

    return(TRUE)
}

# Print usage if run interactively
if (interactive()) {
    cat("[geneSCOPE] geneSCOPE HPC installer loaded.\n")
    cat("[geneSCOPE] Run main() to start installation, or use:\n")
    cat("[geneSCOPE] source('install_genescope_hpc.R')\n")
} else {
    # Run automatically if called as script
    success <- main()
    if (!success) {
        quit(status = 1)
    }
}
