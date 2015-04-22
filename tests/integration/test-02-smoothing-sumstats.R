context("Smoothing SumStats")

test_that("test.smoothing", {
  smooth_jaatha <- Jaatha.initialSearch(smooth_jaatha, 25, 2)
  pStartPoints <- Jaatha.getLikelihoods(smooth_jaatha, initial_search = TRUE)
  expect_equal(nrow(pStartPoints), 4)
  smooth_jaatha <- Jaatha.refinedSearch(smooth_jaatha, 1, 25, 25, max.steps = 9)
  expect_true(ncol(smooth_jaatha@likelihood.table) == 4)
  expect_true(nrow(smooth_jaatha@likelihood.table) >= 5)
})