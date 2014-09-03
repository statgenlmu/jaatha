context("Confidence Intervals")

test_that("Calculation of Confidence Intervals", {
  jaatha <- Jaatha.initialSearch(jaatha.csi, 10, 1)
  jaatha <- Jaatha.refinedSearch(jaatha, 1, 10, max.steps=5)
  jaatha <- Jaatha.confidenceIntervals(jaatha, 0.95, 20, 2)
  estimates <- t(Jaatha.getLikelihoods(jaatha)[1, -(1:2), drop = FALSE])
  conf.ints <- jaatha@conf.ints
  expect_equal(dim(conf.ints), c(2, 2))
  expect_true(all(conf.ints[, 1] <= estimates & estimates <= conf.ints[, 2]))
  
  # Simulate simulation on multiple machines using the subset feature
  logs <- tempfile('j-test-logs')
  Jaatha.confidenceIntervals(jaatha, 0.95, 20, 2, logs, 1:10)
  Jaatha.confidenceIntervals(jaatha, 0.95, 20, 2, logs, 11:20)
  jaatha <- Jaatha.getCIsFromLogs(jaatha, 0.95, logs)
  expect_equal(jaatha@conf.ints, conf.ints)
  unlink(logs, recursive = TRUE)
})