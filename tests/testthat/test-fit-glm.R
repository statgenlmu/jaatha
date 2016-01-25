context("GLM Fitting")


test_that("fit_glm works for PoiInd", {
  model <- create_test_model()
  data <- create_test_data(model)
  block <- create_block(matrix(c(0, 0, 1, 1), 2))
  sim_data <- model$simulate(block$sample_pars(10), data, 1)
  
  glms <- fit_glm(model$get_sum_stats()[[1]], sim_data)
  expect_true(is.list(glms))
  expect_equal(10, length(glms))
  
  glms <- fit_glm(model$get_sum_stats()[[2]], sim_data)
  expect_true(is.list(glms))
  expect_equal(1, length(glms))
  
  expect_error(fit_glm(1:10, sim_data))
})


test_that("glm fitting throws an error if coefficients are 0", {
  model <- create_test_model()
  data <- create_test_data(model)
  block <- create_block(matrix(c(0, 0, 1, 1), 2))
  sim_data <- model$simulate(block$sample_pars(5), data, 1)
  sim_data <- lapply(sim_data, function(x) {
    x$id <- rep(5, 10)
    x$sum <- 0
    x
  })
  
  expect_error(fit_glm(model$get_sum_stats()[[1]], sim_data))
  expect_error(fit_glm(model$get_sum_stats()[[2]], sim_data))
})


test_that("fit_glm works for PoiSmooth", {
  skip("Smoothing not implemented")
  glms_fitted_smooth <- fit_glm(smooth_jaatha, smooth_sim_data)
  expect_true(is.list(glms_fitted_smooth))
  expect_equal(1, length(glms_fitted_smooth))
  expect_true(is.list(glms_fitted_smooth$csi))
  expect_equal(1, length(glms_fitted_smooth$csi))
  expect_true("glm" %in% is(glms_fitted_smooth$csi[[1]]))
})


test_that("fit_glm works", {
  model <- create_test_model()
  data <- create_test_data(model)
  block <- create_block(matrix(c(0, 0, 1, 1), 2))
  sim_data <- model$simulate(block$sample_pars(10), data, 1)
  
  glm_fit <- fit_glm(model, sim_data)
  expect_that(glm_fit, is_a("list"))
  expect_equal(length(glm_fit), 2)
  expect_equal(names(glm_fit), c("id", "sum"))
})
