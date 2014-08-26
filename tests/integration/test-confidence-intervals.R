context("Confidence Intervals")

test_that("seqgen and msms are available", {
  expect_true(test_seqgen)
  expect_true(test_msms)
})

test_that("Calculation of Confidence Intervals", {
  jaatha <- Jaatha.initialSearch(jaatha.csi, 20, 1)
  jaatha <- Jaatha.refinedSearch(jaatha, 1, 20, max.steps=10)
  jaatha <- Jaatha.confidenceIntervals(jaatha, 0.95, 10, 1)
  estimates <- t(Jaatha.getLikelihoods(jaatha)[1, -(1:2), drop = FALSE])
  conf.ints <- jaatha@conf.ints
  expect_equal(dim(conf.ints), c(2, 2))
  expect_true(all(conf.ints[, 1] <= estimates & estimates <= conf.ints[, 2]))
})