context("Confidence Intervals calculation")

test_that("convertSimDataToSumStats works", {
  sim_data <- simulateWithinBlock(1, block.test, jaatha.csi)[[1]]
  sim_data$data <- rep(7.5, 6)
  sum_stats <- convertSimDataToSumStats(sim_data, jaatha.csi@sum.stats)
  expect_equal(length(sum_stats), length(jaatha.csi@sum.stats))
  expect_equal(sum_stats[[1]]$get_data(), rep(7.5, 6))
})