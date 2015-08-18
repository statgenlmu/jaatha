context("Jaatha Model")

test_that("jaatha model can be initialized", {
  model <- create_test_model()
  expect_true(is_par_ranges(model$get_par_ranges()))
  expect_equal(model$get_par_number(), model$get_par_ranges()$get_par_number())
  expect_equal(model$get_scaling_factor(), 1)
})


test_that("adding summary statistics works", {
  sim_func <- function(x) rpois(10, x)
  par_ranges = matrix(c(0.1, 0.1, 10, 10), 2, 2)
  
  model <- create_jaatha_model(sim_func, par_ranges, list(stat_identity()))
  expect_equal(model$get_sum_stats(), list("id" = stat_identity()))
  
  model <- create_jaatha_model(sim_func, par_ranges, list(stat_identity(),
                                                          stat_sum()))
  expect_equal(model$get_sum_stats(), list("id" = stat_identity(), 
                                           "sum" = stat_sum()))
  
  expect_error(create_jaatha_model(sim_func, par_ranges, list(stat_identity(),
                                                              stat_identity())))
})


test_that("simulation works", {
  model <- create_test_model()
  data <- create_test_data(model)
  
  # One parameter combination
  set.seed(1)
  res <- model$simulate(pars = matrix(c(1, 1), 1), data)
  expect_that(res, is_a("list"))
  expect_equal(length(res), 1)
  expect_equal(length(res[[1]]), length(model$get_sum_stats()) + 2)
  expect_equal(names(res[[1]]), 
               c(names(model$get_sum_stats()), "pars", "pars_normal"))
  expect_equivalent(res[[1]]$pars, c(10, 10))
  expect_equivalent(res[[1]]$pars_normal, c(1, 1))
  
  # Check reproducibility
  set.seed(1)
  res2 <- model$simulate(pars = matrix(c(1, 1), 1), data)
  expect_equal(res, res2)
  
  # Two parameter combinations
  res <- model$simulate(pars = matrix(c(1, 0, 1, 0), 2), data)
  expect_equal(length(res), 2)
  expect_equal(length(res[[1]]), length(model$get_sum_stats()) + 2)
  expect_equal(length(res[[2]]), length(model$get_sum_stats()) + 2)
  expect_equivalent(res[[1]]$pars_normal, c(1, 1))
  expect_equivalent(res[[2]]$pars_normal, c(0, 0))
  
  # Errors
  expect_error(model$simulate(pars = matrix(c(1, -0.1), 1), data))
  expect_error(model$simulate(pars = matrix(c(1, 1.1), 1), data))
  expect_error(model$simulate(pars = matrix(c(1, 1), 1), "data"))
  expect_error(model$simulate(pars = matrix(c(1, 1), 1), data, "blub"))
})