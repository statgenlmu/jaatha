# Load the test environment
source("../testsetup/custom_simulation_interface.R")
source("../testsetup/demographic_models.R")

# Check that we can test all aspekts of Jaatha
context("Test Environment")

test_that("seqgen and msms are available", {
  expect_true(test_seqgen)
  expect_true(test_msms)
})