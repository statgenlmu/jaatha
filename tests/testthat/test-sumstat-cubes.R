context("SumStat Cubes")

test_that("initialization of iHH sumstat works", {
  if (!requireNamespace("rehh", quietly = TRUE)) skip("rehh not installed")
  stat <- coala::sumstat_ihh(position = .5)
  ihh = Stat_Ihh$new(sumstat_tt$seg_sites, dm_tt, stat, c(.1, .5, .9))
  expect_that(ihh$get_data(), is_a("integer"))
  expect_that(ihh$get_data(), is_a("integer"))
  expect_that(sum(ihh$get_data()), is_more_than(0))
  expect_that(ihh$get_breaks(), is_a("list"))
  expect_equal(length(ihh$get_breaks()), 1)
})


test_that("initialization of OmegaPrime sumstat works", {
  stat <- coala::sumstat_omegaprime()
  ihh = Stat_OmegaPrime$new(sumstat_tt$seg_sites, dm_tt, stat, c(.1, .5, .9))
  expect_that(ihh$get_data(), is_a("integer"))
  expect_that(ihh$get_data(), is_a("integer"))
  expect_that(sum(ihh$get_data()), is_more_than(0))
  expect_that(ihh$get_breaks(), is_a("list"))
  expect_equal(length(ihh$get_breaks()), 1)
})
