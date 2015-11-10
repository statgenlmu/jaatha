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
  expect_equal(par_ranges$denormalize(c(0, 0)), c(theta = 1, m = 1))
  expect_equal(par_ranges$denormalize(c(1, 1)), c(theta = 5, m = 5))
  
  # Check the summary statistics
  sum_stats <- jaatha_model$get_sum_stats()
  expect_that(sum_stats, is_a("list"))
  expect_equal(length(sum_stats), 1)
})


test_that("conversion of coala sumstats works", {
  skip_if_not_installed("coala")
  model <- coala::coal_model(10:15, 1) + 
    coala::feat_mutation(coala::par_range("theta", 1, 5)) +
    coala::feat_migration(coala::par_range('m', 1, 5), symmetric = TRUE)
  
  expect_equal(convert_coala_sumstats(model), list())
  model <- model + coala::sumstat_sfs()
  expect_equal(length(convert_coala_sumstats(model)), 1)
  model <- model + coala::sumstat_jsfs()
  expect_equal(length(convert_coala_sumstats(model)), 2)
  model <- model + coala::sumstat_ihh()
  stats <- convert_coala_sumstats(model, ihs_breaks = c(.5, .7, .9))
  expect_equal(length(stats), 3)
  
  data <- simulate(model, pars = c(2, 2))
  expect_equal(stats$jsfs$calculate(data, NULL), sum_jsfs(data$jsfs))
  expect_equal(stats$sfs$calculate(data, NULL), data$sfs)
})


test_that("the SFS can be used in Jaatha", {
  skip_if_not_installed("coala")
  coala_model <- coala::coal_model(5, 2) + 
    coala::feat_mutation(coala::par_range("theta", 1, 5)) +
    coala::sumstat_sfs()

  jaatha_model <- create_jaatha_model(coala_model, test = FALSE)
  expect_equal(length(jaatha_model$get_sum_stats()), 1)
  sim_data <- jaatha_model$test()
  value <- jaatha_model$get_sum_stats()[[1]]$calculate(sim_data)
  expect_that(value, is_a("numeric"))
  expect_equal(length(value), 4)
  
  jaatha_model <- create_jaatha_model(coala_model, "sums", test = FALSE)
  jaatha_data <- create_jaatha_data(sim_data, jaatha_model)
  expect_true(is_jaatha_data(jaatha_data))
})


test_that("the JSFS can be used in Jaatha", {
  skip_if_not_installed("coala")
  coala_model <- coala::coal_model(5:6, 2) + 
    coala::feat_mutation(coala::par_range("theta", 1, 5)) +
    coala::feat_migration(coala::par_range("m", 1, 5), symmetric = TRUE) +
    coala::sumstat_jsfs()
  
  jaatha_model <- create_jaatha_model(coala_model, "sums", test = FALSE)
  expect_equal(length(jaatha_model$get_sum_stats()), 1)
  sim_data <- jaatha_model$test()
  value <- jaatha_model$get_sum_stats()[[1]]$calculate(sim_data)
  expect_that(value, is_a("numeric"))
  expect_equal(length(value), 23)
  
  jaatha_model <- create_jaatha_model(coala_model, "sums", test = FALSE)
  jaatha_data <- create_jaatha_data(sim_data, jaatha_model)
  expect_true(is_jaatha_data(jaatha_data))
  
  jaatha_model <- create_jaatha_model(coala_model, "none", test = FALSE)
  expect_equal(length(jaatha_model$get_sum_stats()), 1)
  value <- jaatha_model$get_sum_stats()[[1]]$calculate(sim_data)
  expect_that(value, is_a("numeric"))
  expect_equal(length(value), 6 * 7 - 2)
  

})


test_that("the four gamete stat can be used in Jaatha", {
  skip_if_not_installed("coala")
  
  coala_model <- coala::coal_model(10, 2) + 
    coala::feat_mutation(coala::par_range("theta", 10, 50)) +
    coala::sumstat_four_gamete("fg")
  
  jaatha_model <- create_jaatha_model(coala_model, test = FALSE)
  jaatha_data <- create_jaatha_data(jaatha_model$test(), jaatha_model)
  sim_results <- jaatha_model$simulate(matrix(.5, 1, 1), jaatha_data)
  expect_equal(length(jaatha_model$get_sum_stats()), 1)
  expect_equal(sum(sim_results[[1]]$fg), 2)
})


test_that("MCMF can be used in Jaatha", {
  skip_if_not_installed("coala")
  coala_version <- package_version(packageVersion("coala"))
  if (coala_version$major == 0 & coala_version$minor == 1) {
    skip("Insufficient coala version")
  }
  
  coala_model <- coala::coal_model(10, 2) + 
    coala::feat_mutation(coala::par_range("theta", 10, 50)) +
    coala::sumstat_mcmf("mcmf")
  
  jaatha_model <- create_jaatha_model(coala_model, test = FALSE)
  jaatha_data <- create_jaatha_data(jaatha_model$test(), jaatha_model)
  sim_results <- jaatha_model$simulate(matrix(.5, 1, 1), jaatha_data)
  expect_equal(length(jaatha_model$get_sum_stats()), 1)
  expect_equal(sum(sim_results[[1]]$mcmf), 2)
})


test_that("segsites_list can be used to import data", {
  skip_if_not_installed("coala")
  
  # Create a segsites_list
  seg_sites <- matrix(c(1, 0, 0, 0, 1,
                        1, 1, 0, 1, 0,
                        1, 0, 0, 1, 1,
                        1, 0, 0, 1, 0), 4, 5, byrow = TRUE)
  attr(seg_sites, "positions") <- c(0.1, 0.2, 0.5, 0.7, 0.9)
  segsites_list <- list(seg_sites, seg_sites, seg_sites)
  class(segsites_list) <- "segsites_list"
  
  # Create a coala & jaatha model
  coala_model <- coala::coal_model(c(2, 2), 3, 100) + 
    coala::feat_mutation(coala::par_range("theta", 1, 10)) +
    coala::sumstat_sfs()
  jaatha_model <- create_jaatha_model(coala_model, test = FALSE)
  
  # ------------------------------------------------------------
  # Test import without trios
  # ------------------------------------------------------------
  jaatha_data <- create_jaatha_data(segsites_list, jaatha_model, coala_model)
  expect_equal(jaatha_data$get_values(), list(sfs = c(3, 3, 3)))
  
  # Test that an error is thrown for missing data
  segsites_list2 <- segsites_list
  segsites_list2[2] <- NA
  expect_error(create_jaatha_data(segsites_list2, jaatha_model, coala_model))
  segsites_list2[2] <- NULL
  expect_error(create_jaatha_data(segsites_list2, jaatha_model, coala_model))
  
  # ------------------------------------------------------------
  # Test import with trios
  # ------------------------------------------------------------
  coala_model <- coala::coal_model(c(2, 2)) + 
    coala::locus_trio() +
    coala::feat_mutation(coala::par_range("theta", 1, 10)) +
    coala::sumstat_sfs()
  jaatha_data <- create_jaatha_data(segsites_list, jaatha_model, coala_model, 
                                    trios = list(1:3))
  expect_equal(jaatha_data$get_values(), list(sfs = c(1, 1, 1)))
})
