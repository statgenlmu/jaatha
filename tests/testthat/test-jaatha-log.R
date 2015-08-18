context("Logging")

test_that("log initialization works", {
  model <- create_test_model()
  log <- create_jaatha_log(2, 100, model)
  expect_equivalent(log$get_estimates(1), 
                    data.frame(rep(NA, 100), NA, NA, NA, NA))
  expect_equivalent(log$get_estimates(2), 
                    data.frame(rep(NA, 100), NA, NA, NA, NA))
})


test_that("logging estimates works", {
  model <- create_test_model()
  log <- create_jaatha_log(2, 10, model)
  
  estimate <- list(par = 1:model$get_par_number(), value = 0.1)
  log$log_estimate(1, 1, estimate)

  expect_equivalent(log$get_estimates(1)[1, ], 
                    data.frame(1, 1, estimate$value, 1, 2))
  
  estimate$value <- 0.2
  log$log_estimate(1, 2, estimate)
  expect_equivalent(log$get_estimates(1)[2, ], 
                    data.frame(1, 2, estimate$value, 1, 2))
  
  log$log_estimate(2, 1, estimate)
  expect_equivalent(log$get_estimates(2)[1, ], 
                    data.frame(2, 1, estimate$value, 1, 2))
})


test_that("getting best estimates works", {
  model <- create_test_model()
  log <- create_jaatha_log(2, 10, model)
  for (i in c(4, 3, 1, 5, 10, 2, 7, 9, 8) / 10) {
    log$log_estimate(1, i * 10, list(par = rep(i, model$get_par_number()), 
                                     value = i))
  }
  
  for (i in 11:20 / 10) {
    log$log_estimate(2, i * 10, list(par = rep(i, model$get_par_number()), 
                                     value = i))
  }
  
  expect_equivalent(log$get_best_estimates(1),
                    data.frame(rep = c(2, 1),
                               step = c(20, 10),
                               llh = c(2.0, 1.0),
                               p1 = c(2.0, 1.0),
                               p2 = c(2.0, 1.0)))
  
  expect_equivalent(log$get_best_estimates(2),
                    data.frame(rep = c(2, 2, 1, 1),
                               step = c(20, 19, 10, 9),
                               llh = c(2.0, 1.9, 1.0, 0.9),
                               p1 = c(2.0, 1.9, 1.0, 0.9),
                               p2 = c(2.0, 1.9, 1.0, 0.9)))
})


test_that("getting best estimates works", {
  model <- create_test_model()
  log <- create_jaatha_log(2, 10, model)
  log$log_estimate("final", 1, list(par = rep(1, model$get_par_number()),
                                    value = 1))
  log$log_estimate("final", 2, list(par = rep(2, model$get_par_number()),
                                    value = 2))
  
  expect_equivalent(log$create_results(),
                    list(param = c(2, 2),
                         loglikelihood = 2,
                         converged = FALSE))
  
  log$log_convergence(1)
  expect_equal(log$create_results()$converged, FALSE)
  log$log_convergence(2)
  expect_equal(log$create_results()$converged, TRUE)
  
})