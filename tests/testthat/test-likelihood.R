context("Likelihood estimation")

test_that("llh is approximatied for basic statistics", {
  model <- create_test_model()
  data <- create_test_data(model)
  block <- create_block(matrix(c(0, 0, 1, 1), 2))
  sim_data <- model$simulate(block$sample_pars(10), data, 1)
  glms <- fit_glm(model, sim_data)
  
  llh <- approximate_llh(model$get_sum_stats()[[1]], data, c(.5, .5), glms, 1)
  expect_true(is.numeric(llh))
  expect_true(llh <= 0)
  
  llh2 <- approximate_llh(model$get_sum_stats()[[1]], data, c(.5, .5), glms, 2)
  expect_true(is.numeric(llh2))
  expect_true(llh2 <= 0)
  expect_true(llh != llh2)
})


test_that("calcStatLLH works for Stat_PoiSmooth", {
  skip("Smoothing not implemented")
  glm_fit <- fit_glm(smooth_stat, smooth_sim_data)
  
  ll <- calcStatLLH(smooth_stat, glm_fit, c(x=.5, y=.5), scaling_factor=1)
  expect_true(is.numeric(ll))
  expect_true(0 <= exp(ll) | exp(ll) <= 1)
})


test_that("llh is approximatied for complete models", {
  model <- create_test_model()
  data <- create_test_data(model)
  block <- create_block(matrix(c(0, 0, 1, 1), 2))
  sim_data <- model$simulate(block$sample_pars(10), data, 1)
  glms <- fit_glm(model, sim_data)
  
  llh <- approximate_llh(model, data, c(.5, .5), glms)
  expect_true(is.numeric(llh))
  expect_true(llh <= 0)
  
  llh1 <- approximate_llh(model$get_sum_stats()[[1]], data, c(.5, .5), glms, 1)
  llh2 <- approximate_llh(model$get_sum_stats()[[2]], data, c(.5, .5), glms, 1)
  expect_equal(llh, llh1 + llh2)
})


test_that("llh optimization works", {
  model <- create_test_model()
  data <- create_test_data(model)
  block <- create_block(matrix(c(0, 0, 1, 1), 2))
  sim_data <- model$simulate(block$sample_pars(10), data, 1)
  glms <- fit_glm(model, sim_data)
  
  opt_llh <- optimize_llh(block, model, data, glms)
  expect_equal(length(opt_llh$par), 2)
  expect_true(all(opt_llh$par > 0 & opt_llh$par < 1))
  expect_true(opt_llh$value < 0)
})
