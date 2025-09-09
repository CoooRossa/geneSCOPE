#!/usr/bin/env Rscript
# HPC Environment Diagnostics for geneSCOPE
# This script helps identify specific issues in HPC environments

cat("[geneSCOPE] Starting HPC environment diagnostics...\n")

# System information gathering
gather_system_info <- function() {
    cat("\n[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] System Information\n")
    cat("[geneSCOPE] ========================================\n")

    # R version
    r_info <- R.Version()
    cat(sprintf(
        "[geneSCOPE] R version: %s.%s.%s (%s)\n",
        r_info$major, r_info$minor, r_info$year, r_info$platform
    ))

    # Operating system
    cat(sprintf("[geneSCOPE] OS: %s\n", Sys.info()["sysname"]))
    cat(sprintf("[geneSCOPE] Machine: %s\n", Sys.info()["machine"]))

    # Library paths
    cat("[geneSCOPE] Library paths:\n")
    for (i in seq_along(.libPaths())) {
        cat(sprintf("[geneSCOPE]   %d: %s\n", i, .libPaths()[i]))
    }

    # Environment variables
    cat("[geneSCOPE] Key environment variables:\n")
    env_vars <- c("R_LIBS_USER", "R_LIBS_SITE", "OMP_NUM_THREADS", "DISPLAY")
    for (var in env_vars) {
        value <- Sys.getenv(var, unset = NA)
        if (!is.na(value) && value != "") {
            cat(sprintf("[geneSCOPE]   %s = %s\n", var, value))
        }
    }
}

# Check compiler availability
check_compilers <- function() {
    cat("\n[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] Compiler Check\n")
    cat("[geneSCOPE] ========================================\n")

    # Check if Rcpp can compile
    tryCatch(
        {
            library(Rcpp)
            cat("[geneSCOPE] ✓ Rcpp available\n")

            # Try a simple compilation
            simple_cpp <- "int add(int x, int y) { return x + y; }"
            Rcpp::cppFunction(simple_cpp)
            cat("[geneSCOPE] ✓ C++ compilation works\n")
        },
        error = function(e) {
            cat("[geneSCOPE] !!! C++ compilation failed:", e$message, "\n")
        }
    )

    # Check make tools
    make_available <- nzchar(Sys.which("make"))
    cat(sprintf("[geneSCOPE] Make available: %s\n", make_available))

    # Check compiler
    gcc_available <- nzchar(Sys.which("gcc"))
    gpp_available <- nzchar(Sys.which("g++"))
    cat(sprintf("[geneSCOPE] GCC available: %s\n", gcc_available))
    cat(sprintf("[geneSCOPE] G++ available: %s\n", gpp_available))
}

# Test package installation capability
test_package_installation <- function() {
    cat("\n[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] Package Installation Test\n")
    cat("[geneSCOPE] ========================================\n")

    # Test simple package installation
    test_pkg <- "cli" # Small, common package

    cat(sprintf("[geneSCOPE] Testing installation of '%s'...\n", test_pkg))

    tryCatch(
        {
            if (!require(test_pkg, character.only = TRUE, quietly = TRUE)) {
                install.packages(test_pkg,
                    repos = "https://cloud.r-project.org/",
                    dependencies = FALSE, quiet = TRUE
                )

                if (require(test_pkg, character.only = TRUE, quietly = TRUE)) {
                    cat("[geneSCOPE] ✓ Package installation works\n")
                } else {
                    cat("[geneSCOPE] !!! Package installation failed\n")
                }
            } else {
                cat("[geneSCOPE] ✓ Test package already available\n")
            }
        },
        error = function(e) {
            cat("[geneSCOPE] !!! Package installation error:", e$message, "\n")
        }
    )
}

# Check specific problematic packages
check_problematic_packages <- function() {
    cat("\n[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] Problematic Package Check\n")
    cat("[geneSCOPE] ========================================\n")

    problematic <- c("httpgd", "devtools", "ggplot2", "sf", "rgeos")

    for (pkg in problematic) {
        cat(sprintf("[geneSCOPE] Checking %s: ", pkg))

        if (require(pkg, character.only = TRUE, quietly = TRUE)) {
            version <- packageVersion(pkg)
            cat(sprintf("✓ available (v%s)\n", version))
        } else {
            cat("✗ not available\n")

            # Try to install and see what happens
            cat(sprintf("[geneSCOPE] Attempting to install %s...\n", pkg))
            tryCatch(
                {
                    install.packages(pkg,
                        repos = "https://cloud.r-project.org/",
                        dependencies = FALSE, quiet = TRUE
                    )
                    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
                        cat(sprintf("[geneSCOPE] ✓ %s installed successfully\n", pkg))
                    } else {
                        cat(sprintf("[geneSCOPE] !!! %s installation failed\n", pkg))
                    }
                },
                error = function(e) {
                    cat(sprintf("[geneSCOPE] !!! %s error: %s\n", pkg, e$message))
                }
            )
        }
    }
}

# Check internet connectivity and repositories
check_repositories <- function() {
    cat("\n[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] Repository Access Check\n")
    cat("[geneSCOPE] ========================================\n")

    # Current repositories
    repos <- getOption("repos")
    cat("[geneSCOPE] Configured repositories:\n")
    for (name in names(repos)) {
        cat(sprintf("[geneSCOPE]   %s: %s\n", name, repos[name]))
    }

    # Test CRAN access
    cat("[geneSCOPE] Testing CRAN access...\n")
    tryCatch(
        {
            available <- available.packages(repos = "https://cloud.r-project.org/")[1:5, 1:2]
            cat("[geneSCOPE] ✓ CRAN accessible\n")
            cat(sprintf(
                "[geneSCOPE] Sample packages available: %s\n",
                paste(rownames(available)[1:3], collapse = ", ")
            ))
        },
        error = function(e) {
            cat("[geneSCOPE] !!! CRAN access failed:", e$message, "\n")
        }
    )

    # Test GitHub access (for devtools)
    cat("[geneSCOPE] Testing GitHub access...\n")
    tryCatch(
        {
            url <- "https://api.github.com/repos/CoooRossa/geneSCOPE"
            response <- readLines(url, warn = FALSE)
            cat("[geneSCOPE] ✓ GitHub accessible\n")
        },
        error = function(e) {
            cat("[geneSCOPE] !!! GitHub access failed:", e$message, "\n")
        }
    )
}

# Generate recommendations
generate_recommendations <- function() {
    cat("\n[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] Recommendations\n")
    cat("[geneSCOPE] ========================================\n")

    cat("[geneSCOPE] Based on the diagnostics above:\n\n")

    cat("[geneSCOPE] 1. If C++ compilation failed:\n")
    cat("[geneSCOPE]    - Load development tools: module load gcc\n")
    cat("[geneSCOPE]    - Check: which gcc; which g++\n\n")

    cat("[geneSCOPE] 2. If package installation failed:\n")
    cat("[geneSCOPE]    - Use: install.packages('pkg', dependencies=FALSE)\n")
    cat("[geneSCOPE]    - Set: options(repos='https://cloud.r-project.org/')\n\n")

    cat("[geneSCOPE] 3. If httpgd causes issues:\n")
    cat("[geneSCOPE]    - Use: emergency_hpc_install.R\n")
    cat("[geneSCOPE]    - Disable graphics: options(device='png')\n\n")

    cat("[geneSCOPE] 4. If internet access is limited:\n")
    cat("[geneSCOPE]    - Download packages manually\n")
    cat("[geneSCOPE]    - Install from local files\n\n")

    cat("[geneSCOPE] 5. For geneSCOPE installation:\n")
    cat("[geneSCOPE]    - First run: install_hpc_dependencies.R\n")
    cat("[geneSCOPE]    - Then run: install_genescope_hpc.R\n")
}

# Main diagnostic function
run_diagnostics <- function() {
    cat("[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] HPC Environment Diagnostics\n")
    cat("[geneSCOPE] ========================================\n")

    gather_system_info()
    check_compilers()
    test_package_installation()
    check_problematic_packages()
    check_repositories()
    generate_recommendations()

    cat("\n[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] Diagnostics Complete\n")
    cat("[geneSCOPE] ========================================\n")
    cat("[geneSCOPE] Save this output and share with support if needed.\n")
}

# Auto-run if called as script
if (!interactive()) {
    run_diagnostics()
} else {
    cat("[geneSCOPE] HPC diagnostics loaded. Run run_diagnostics() to start.\n")
}
