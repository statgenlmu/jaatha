context("Confidence Intervals calculation")

test_that("convertSimDataToSumStats works", {
  sim.data <- simulateWithinBlock(10, block.test, jaatha.csi)
  sim.data$poisson.vector <- rep(7.5, 6)
  sum.stats <- convertSimDataToSumStats(sim.data, jaatha.csi@sum.stats)
  expect_true(length(sum.stats) == length(jaatha.csi@sum.stats))
  expect_true(sum.stats$poisson.vector$method == jaatha.csi@sum.stats$poisson.vector$method)
  expect_true(all(sum.stats$poisson.vector$value == 7.5))
  expect_true(all(sum.stats$poisson.vector$value.transformed ==  7.5))
})