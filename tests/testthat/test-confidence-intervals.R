context("Confidence Intervals calculation")

# test_that("convertSimDataToSumStats works", {
#   sim_data <- simulateWithinBlock(1, block.test, jaatha.csi)[[1]]
#   sim_data$data <- rep(7.5, 6) 
#   sum_stats <- list(csi=R6::R6Class("Stat_PoiInd", 
#                                     inherit = jaatha:::Stat_Base)$new(csi.obs, 'csi'))
#   sum_stats <- convertSimDataToSumStats(sim_data, sum_stats)
#   expect_equal(length(sum_stats), length(jaatha.csi@sum.stats))
#   expect_false(all(jaatha.csi@sum.stats[[1]]$get_data() == 7.5))
#   expect_equal(sum_stats[[1]]$get_data(), rep(7.5, 6))
# })