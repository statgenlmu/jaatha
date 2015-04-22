context("coala initialization")

set.seed(112233)
model <- coala::coal_model(10:11, 10, 100) +
  coala::feat_pop_merge(coala::par_range("tau", .1, 2), 2, 1) +
  coala::feat_mutation(coala::par_range("theta", 1, 10)) +
  coala::feat_recombination(1) +
  coala::sumstat_jsfs()

data <- simulate(model + coala::sumstat_seg_sites(), pars=c(1, 5))


test_that("coala initialization works", {  
  # dm + jsfs
  expect_true(is_jaatha(Jaatha.initialize(data, model)))
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
  expect_equal(coala::get_locus_number(jaatha@opts[['dm']]), 2L)
})


test_that("initialization with FPC statistic works", {
  jaatha.fpc <- Jaatha.initialize(data, model + 
                                    coala::sumstat_four_gamete('fgc'))
  
  fpc_stat <- jaatha.fpc@sum_stats[['fgc']]
  expect_false(is.null(fpc_stat))
  expect_that(sum(fpc_stat$get_data()), is_more_than(0))
  expect_that(sum(fpc_stat$get_data()), 
              is_less_than(coala:::get_locus_number(dm_tt)+1))
})


test_that("initialization with iHH statistic works", {
  if (!requireNamespace("rehh", quietly = TRUE)) skip("rehh not installed")
  model <- model + coala::sumstat_ihh("ihh", position = .5)
  jaatha.fpc <- Jaatha.initialize(data, model)
  
  stat <- jaatha.fpc@sum_stats[['ihh']]
  expect_false(is.null(stat))
  expect_that(sum(stat$get_data()), is_more_than(0))
  expect_that(sum(stat$get_data()), 
              is_less_than(coala:::get_locus_number(dm_tt)+1))
})


test_that("initialization with Omega' statistic works", {
  model <- model + coala::sumstat_omegaprime("op")
  jaatha.fpc <- Jaatha.initialize(data, model)
  
  stat <- jaatha.fpc@sum_stats[['op']]
  expect_false(is.null(stat))
  expect_that(sum(stat$get_data()), is_more_than(0))
  expect_that(sum(stat$get_data()), 
              is_less_than(coala:::get_locus_number(dm_tt)+1))
})
