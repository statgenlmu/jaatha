context("Jaatha estimation")

test_that("main function works", {
  model <- create_test_model()
  data <- create_test_data(model)
  
  results <- jaatha(model, data, repetitions = 2, sim = 10, cores = 1, 
                    max_steps = 15)
  
  expect_is(results, "list")
  expect_true(is.finite(results$loglikelihood))
  expect_true(all(results$param > 1))
  expect_true(is_single_logical(results$converged))
})