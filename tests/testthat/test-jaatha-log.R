context("Logging")

test_that("log initialization works", {
  model <- create_test_model()
  log <- create_jaatha_log(2, 100, model)
  expect_equivalent(log$get_estimates(), matrix(NA, 200, 5))
})


test_that("logging estimates works", {
  model <- create_test_model()
  log <- create_jaatha_log(2, 10, model)
  
  estimate <- list(par = 1:model$get_par_number(), value = 0.1)
  log$log_estimate(1, 1, estimate)
  expect_equivalent(log$get_estimates()[1, ], 
                    c(1, 1, estimate$value, estimate$par))
  
  log$log_estimate(1, 2, estimate)
  expect_equivalent(log$get_estimates()[1, ], 
                    c(1, 1, estimate$value, estimate$par))
  expect_equivalent(log$get_estimates()[2, ], 
                    c(1, 2, estimate$value, estimate$par))
  
  estimate$value <- 0.2
  log$log_estimate(1, 2, estimate)
  expect_equivalent(log$get_estimates()[2, ], 
                    c(1, 2, estimate$value, estimate$par))
  
  log$log_estimate(2, 1, estimate)
  expect_equivalent(log$get_estimates()[11, ], 
                    c(2, 1, estimate$value, estimate$par))
})
