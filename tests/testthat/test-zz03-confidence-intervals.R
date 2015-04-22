context("Confidence Intervals")


test_that("convertSimDataToSumStats works", {
  skip("Need to clone sumstats for conf. intervals")
  sim_data <- simulateWithinBlock(1, block.test, jaatha.csi)[[1]]
  sim_data$data <- rep(7.5, 6) 
  sum_stats <- list(csi=R6::R6Class("Stat_PoiInd", 
                                    inherit = jaatha:::Stat_Base)$new(csi.obs, 'csi'))
  sum_stats <- convertSimDataToSumStats(sim_data, sum_stats)
  expect_equal(length(sum_stats), length(jaatha.csi@sum.stats))
  expect_false(all(jaatha.csi@sum.stats[[1]]$get_data() == 7.5))
  expect_equal(sum_stats[[1]]$get_data(), rep(7.5, 6))
})


test_that("Calculation of Confidence Intervals", {
  skip("Need to clone sumstats for conf. intervals")
  skip_on_cran()
  jaatha <- Jaatha.initialSearch(jaatha.csi, 10, 2)
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