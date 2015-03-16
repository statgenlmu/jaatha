context("coalsimr initialization")

model <- coalsimr::coal_model(10:11, 10, 100) +
  coalsimr::feat_pop_merge(coalsimr::par_range("tau", .1, 2), 2, 1) +
  coalsimr::feat_mutation(coalsimr::par_range("theta", 1, 10)) +
  coalsimr::feat_recombination(1) +
  coalsimr::sumstat_jsfs()

data <- simulate(model + coalsimr::sumstat_seg_sites(), pars=c(1, 5))


test_that("coalsimr initialization works", {  
  # dm + jsfs
  checkType(Jaatha.initialize(data, model), "jaatha")
  expect_error(Jaatha.initialize(model))
  expect_error(Jaatha.initialize(NULL))
  
  # cores
  jaatha <- Jaatha.initialize(data, model)
  expect_equal(jaatha@cores, 1)
  jaatha <- Jaatha.initialize(sumstat_tt, dm_tt, cores = 2)
  expect_equal(jaatha@cores, 2)
  
  # smoothing
  jaatha <- Jaatha.initialize(data, model, smoothing = TRUE)
  expect_equal(length(jaatha@sum_stats), 2)
  expect_that(jaatha@sum_stats$jsfs$get_data(), is_a('data.frame'))
  expect_equal(colnames(jaatha@sum_stats$jsfs$get_data()), c('X1', 'X2', 'sum.stat'))
  expect_that(sum(jaatha@sum_stats$jsfs$get_data()), is_more_than(0))
  expect_that(jaatha@sum_stats$jsfs$get_model(), is_a('character'))
  
  expect_that(jaatha@sum_stats$border_jsfs, is_a('Stat_PoiInd'))
  expect_that(sum(jaatha@sum_stats$border_jsfs$get_data()), is_more_than(0))
})


test_that("PG initialization works with scaling", {
  skip("Temporarily deactived")
  jaatha <- Jaatha.initialize(data, model, scaling_factor = 5)
  expect_equal(jaatha@scaling_factor, 5)
  expect_equal(coalsimr::get_locus_number(jaatha@opts[['dm']]), 2L)
})


test_that("initialization with FPC statistic", {
  jaatha.fpc <- Jaatha.initialize(data, model + 
                                    coalsimr::sumstat_four_gamete('fgc'))
  
  fpc_stat <- jaatha.fpc@sum_stats[['fgc']]
  expect_false(is.null(fpc_stat))
  expect_that(sum(fpc_stat$get_data()), is_more_than(0))
  expect_that(sum(fpc_stat$get_data()), 
              is_less_than(coalsimr:::get_locus_number(dm_tt)+1))
})


# test_that("Initialization of PMC statistic", {
#   dm.pmc <- dm.addSummaryStatistic(dm.fpc, 'pmc')
#   
#   # Without groups
#   jaatha.pmc <- Jaatha.initialize(sum.stats.fpc, dm.pmc)
#   expect_true(sum(jaatha.pmc@sum.stats[["pmc"]]$value) > 0)
#   expect_false(is.null(jaatha.pmc@opts$dm@options[["pmc_breaks_private"]]))
#   expect_false(is.null(jaatha.pmc@opts$dm@options[["pmc_breaks_fixed"]]))  
#   
#   # With groups
#   dm.pmc <- dm.addSummaryStatistic(dm.grp, 'pmc', group=1)
#   dm.pmc <- dm.addSummaryStatistic(dm.pmc, 'pmc', group=3)
#   
#   jaatha.pmc <- Jaatha.initialize(sum.stats.grp, dm.pmc)
#   expect_equal(length(jaatha.pmc@sum.stats), 5)
#   expect_false(is.null(jaatha.pmc@sum.stats$pmc.1))
#   expect_true(sum(jaatha.pmc@sum.stats[["pmc.1"]]$value) > 0)
#   expect_true(is.null(jaatha.pmc@sum.stats$pmc.2))
#   expect_false(is.null(jaatha.pmc@sum.stats$pmc.3))
#   expect_true(sum(jaatha.pmc@sum.stats[["pmc.3"]]$value) > 0)
# })