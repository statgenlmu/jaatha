context("Coala Interface")

test_that("it creates a jaatha model from a coala one", {
  skip_if_not_installed("coala")
  coala_model <- coala::coal_model(10:15, 100) + 
    coala::feat_mutation(coala::par_range("theta", 1, 5)) +
    coala::feat_migration(coala::par_range('m', 1, 5), symmetric = TRUE) +
    coala::sumstat_sfs()
  jaatha_model <- create_jaatha_model(coala_model, test = FALSE)
  
  # Check the par_ranges
  par_ranges <- jaatha_model$get_par_ranges()
  expect_equal(par_ranges$get_par_number(), 2)
  expect_equal(par_ranges$get_par_names(), c("theta", "m"))
  expect_equal(par_ranges$denormalize(c(0, 0)), c(theta=1, m=1))
  expect_equal(par_ranges$denormalize(c(1, 1)), c(theta=5, m=5))
  
  # Check the summary statistics
  sum_stats <- jaatha_model$get_sum_stats()
  expect_that(sum_stats, is_a("list"))
  expect_equal(length(sum_stats), 1)
})


test_that("conversion of coala sumstats works", {
  skip_if_not_installed("coala")
  model <- coala::coal_model(10:15, 100) + 
    coala::feat_mutation(coala::par_range("theta", 1, 5)) +
    coala::feat_migration(coala::par_range('m', 1, 5), symmetric = TRUE)
  
  expect_equal(convert_coala_sumstats(model), list())
  expect_equal(length(convert_coala_sumstats(model + coala::sumstat_sfs())), 1)
  #expect_equal(length(convert_coala_sumstats(model + coala::sumstat_jsfs())), 1)
})


# test_that("initialization of iHH sumstat works", {
#   skip_if_not_installed("rehh")
#   stat <- coala::sumstat_ihh(position = .5)
#   ihh = Stat_Ihh$new(sumstat_tt$seg_sites, dm_tt, stat, c(.1, .5, .9))
#   expect_that(ihh$get_data(), is_a("integer"))
#   expect_that(ihh$get_data(), is_a("integer"))
#   expect_that(sum(ihh$get_data()), is_more_than(0))
#   expect_that(ihh$get_breaks(), is_a("list"))
#   expect_equal(length(ihh$get_breaks()), 1)
#   
#   stat <- coala::sumstat_ihh()
#   ihh = Stat_Ihh$new(sumstat_tt$seg_sites, dm_tt, stat, c(.1, .5, .9))
#   expect_that(ihh$get_data(), is_a("integer"))
#   expect_that(ihh$get_data(), is_a("integer"))
#   expect_that(sum(ihh$get_data()), is_more_than(0))
#   expect_that(ihh$get_breaks(), is_a("list"))
#   expect_equal(length(ihh$get_breaks()), 1)
# })
# 
# 
# test_that("initialization of OmegaPrime sumstat works", {
#   stat <- coala:::sumstat_omegaprime()
#   ihh = Stat_OmegaPrime$new(sumstat_tt$seg_sites, dm_tt, stat, c(.1, .5, .9))
#   expect_that(ihh$get_data(), is_a("integer"))
#   expect_that(ihh$get_data(), is_a("integer"))
#   expect_that(sum(ihh$get_data()), is_more_than(0))
#   expect_that(ihh$get_breaks(), is_a("list"))
#   expect_equal(length(ihh$get_breaks()), 1)
# })