test_that("geneSCOPE package loads correctly", {
  expect_true(require(geneSCOPE, quietly = TRUE))
  expect_true("geneSCOPE" %in% loadedNamespaces())
})

test_that("basic utility functions work", {
  # Test detectOS function
  os <- detectOS()
  expect_type(os, "character")
  expect_true(nchar(os) > 0)
  expect_true(os %in% c("windows", "macos", "linux", "unknown"))

  # Test getSafeThreadCount function with required parameter
  threads <- getSafeThreadCount(4)
  expect_type(threads, "double") # Function returns numeric, not integer
  expect_gte(threads, 1)
  expect_lte(threads, 4)
})
