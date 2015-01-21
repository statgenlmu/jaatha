context("GLM Fitting")

test_that("fitGlm works for PoiInd", {
  glms.fitted.csi <- fitGlm(jaatha.csi, sim.data.csi)
  expect_true(is.list(glms.fitted.csi))
  expect_equal(1, length(glms.fitted.csi))
  expect_equal(length(sim.data.csi[[1]]$poisson.vector), 
               length(glms.fitted.csi$poisson.vector))
  expect_true(all(sapply(glms.fitted.csi$csi, is)[1, ] == "glm"))
})

test_that("fitGlm works for PoiSmooth", {
  glms.fitted.smooth <- fitGlm(smooth_jaatha, smooth_sim_data)
  expect_true(is.list(glms.fitted.smooth))
  expect_equal(1, length(glms.fitted.smooth))
  expect_true(is.list(glms.fitted.smooth$csi))
  expect_equal(1, length(glms.fitted.smooth$csi))
  expect_true("glm" %in% is(glms.fitted.smooth$csi[[1]]))
})

test_that("fitGlm works", {
  glm_fit <- fitGlm(jaatha.csi, sim.data.csi)
  expect_that(glm_fit, is_a('list'))
  expect_equal(length(glm_fit), 1)
  expect_false(is.null(glm_fit$csi))
  expect_that(glm_fit$csi, is_a('list'))
})
