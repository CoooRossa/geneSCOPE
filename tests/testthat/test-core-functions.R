# Basic Package Functionality Tests for geneSCOPE
# Focus on essential functions that should always work

library(testthat)
library(geneSCOPE)

test_that("detectOS function works correctly", {
    os_result <- detectOS()
    expect_true(os_result %in% c("windows", "macos", "linux", "unknown"))
    expect_type(os_result, "character")
})

test_that("getSafeThreadCount returns valid thread count", {
    thread_count <- getSafeThreadCount(4)
    expect_type(thread_count, "double") # Function returns numeric
    expect_gte(thread_count, 1)
    expect_lte(thread_count, 4)
})

test_that("fixRNG function sets reproducible state", {
    # Test that fixRNG creates reproducible results
    fixRNG(42)
    rand1 <- runif(5)

    fixRNG(42)
    rand2 <- runif(5)

    expect_equal(rand1, rand2)
})

# Test that critical functions exist and are exported
test_that("Core functions are properly exported", {
    core_functions <- c(
        "createSCOPE", "computeWeights",
        "addLeeStats", "detectOS", "fixRNG"
    )

    for (func in core_functions) {
        expect_true(exists(func),
            info = paste("Function", func, "should be exported")
        )
        expect_true(is.function(get(func)),
            info = paste(func, "should be a function")
        )
    }
})

test_that("Package loads without errors", {
    expect_silent(library(geneSCOPE))

    # Check that package has expected metadata
    expect_true("geneSCOPE" %in% loadedNamespaces())

    # Check version exists
    version_info <- packageVersion("geneSCOPE")
    expect_s3_class(version_info, "package_version")
})
