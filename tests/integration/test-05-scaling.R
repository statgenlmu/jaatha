context("Scaling")

test_that("Scaling produces resonable results", {
  dm <- dm.createDemographicModel(5:6, 200)
  dm <- dm.addMutation(dm, 1, 10)
  dm <- dm.addSpeciationEvent(dm, time.point = ".5")
  sum_stats <- dm.simSumStats(dm.addSummaryStatistic(dm, 'seg.sites'), c(2))
  
  jaatha <- Jaatha.initialize(sum_stats, dm, scaling.factor = 20)
  expect_equal(dm.getLociNumber(jaatha@opts[['dm']]), 10)
  jaatha <- Jaatha.initialSearch(jaatha, sim = 50, blocks.per.par = 1)
  expect_that(Jaatha.getStartingPoints(jaatha)[2], is_less_than(5))
  
  jaatha <- Jaatha.refinedSearch(jaatha, 1, 20, max.steps = 20)
  expect_true(all(Jaatha.getLikelihoods(jaatha)[,3] < 4))
})