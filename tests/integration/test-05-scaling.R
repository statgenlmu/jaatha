context("Scaling")

test_that("Scaling produces resonable results", {
  set.seed(2378421)
  jaatha <- Jaatha.initialize(sumstat_tt, dm_tt, scaling_factor = 5)
  jaatha <- Jaatha.initialSearch(jaatha, sim = 100, blocks.per.par = 1)
  expect_that(Jaatha.getLikelihoods(jaatha, initial_search = TRUE)[4], 
              is_less_than(7))
  
  jaatha <- Jaatha.refinedSearch(jaatha, 1, 20, max.steps = 20)
  expect_that(Jaatha.getLikelihoods(jaatha)[4], is_less_than(7))
})