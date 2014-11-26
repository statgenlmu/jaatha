# Load the test environment
source("../testthat/setup-custom-simulation-interface.R")
source("../testthat/setup-demographic-models.R")

# Check that we can test all aspekts of Jaatha
context("Test Environment")

test_that("seqgen and msms are available", {
  expect_true(test_seqgen)
  expect_true(test_msms)
})