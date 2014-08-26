context("Smoothing SumStats")

test_that("test.smoothing", {
  smooth.jaatha <- Jaatha.initialSearch(smooth.jaatha, 25, 2)
  pStartPoints <- Jaatha.getStartingPoints(smooth.jaatha)
  expect_equal(nrow(pStartPoints), 4)
  smooth.jaatha <- Jaatha.refinedSearch(smooth.jaatha, 1, 25, 25, max.steps = 9)
  expect_true(ncol(smooth.jaatha@likelihood.table) == 4)
  expect_true(nrow(smooth.jaatha@likelihood.table) >= 5)
})

test_that("test.smoothing.border", {
  smooth.jaatha <- Jaatha.initialSearch(smooth.border.jaatha, 25, 2)
  pStartPoints <- Jaatha.getStartingPoints(smooth.jaatha)
  expect_equal(nrow(pStartPoints), 4)
})