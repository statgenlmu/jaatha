context("SumStat SFS")

model <- coala::coal_model(10, 10) + 
  coala::feat_mutation(coala::par_range("theta", 1, 10)) +
  coala::sumstat_sfs("SFS_01") +
  coala::sumstat_seg_sites()
data <- simulate(model, pars = 5)


test_that("using the SFS works", {
  set.seed(17)
  sfs = Stat_sfs$new(data$seg_sites, model,
                     coala::get_summary_statistics(model)$SFS_01)
  
  expect_that(sum(sfs$get_data()), is_more_than(0))
  expect_equal(sfs$transform(data), sfs$get_data())
})


test_that("initialization with sfs works", {
  capture.output(j <- Jaatha.initialize(data, model))
})