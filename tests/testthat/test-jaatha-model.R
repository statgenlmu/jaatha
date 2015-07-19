context("Jaatha Model")

test_that("jaatha model can be initialized", {
  model <- create_test_model()
})


test_that("adding summary statistics works", {
  sim_func <- function(x) rpois(10, x)
  par_ranges = matrix(c(0.1, 0.1, 10, 10), 2, 2)
  
  model <- create_jaatha_model(sim_func, par_ranges, list(stat_identity()))
  expect_equal(model$get_sum_stats(), list("id"=stat_identity()))
  
  model <- create_jaatha_model(sim_func, par_ranges, list(stat_identity(),
                                                          stat_sum()))
  expect_equal(model$get_sum_stats(), list("id"=stat_identity(), 
                                           "sum" = stat_sum()))
  
  expect_error(create_jaatha_model(sim_func, par_ranges, list(stat_identity(),
                                                              stat_identity())))
})


test_that("simulation works", {
  model <- create_test_model()
  data <- create_test_data(model)
  
  res <- model$simulate(pars = c(1, 1), seed = 1, data)
  expect_that(res, is_a("list"))
  expect_equal(length(res), length(model$get_sum_stats()) + 2)
  expect_equal(names(res), 
               c(names(model$get_sum_stats()), "pars", "pars_normal"))
  expect_equivalent(res$pars, c(10, 10))
  expect_equivalent(res$pars_normal, c(1, 1))
  
  res2 <- model$simulate(pars = c(1, 1), seed = 1, data)
  expect_equal(res, res2)
})
