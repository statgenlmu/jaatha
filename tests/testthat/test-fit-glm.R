context("GLM Fitting")

test_that("test.convertSimResultsToDataFrame", {
  smooth.df <- convertSimResultsToDataFrame(smooth.sim.data, 
                                            "mat")
  expect_true(is.data.frame(smooth.df))
  expect_equal(5, ncol(smooth.df))
  expect_equal(length(as.vector(smooth.sim.data[[1]]$mat)) * 
                 length(smooth.sim.data), nrow(smooth.df))
  expect_false(is.null(colnames(smooth.df)))
  expect_true(all(colnames(smooth.df) == c("x", "y", "X1", 
                                           "X2", "sum.stat")))
  smooth.df <- convertSimResultsToDataFrame(smooth.sim.data, 
                                            "mat", smooth.border.sum.stats$mat$border.mask)
  expect_true(is.data.frame(smooth.df))
  expect_equal(5, ncol(smooth.df))
  expect_true(all(smooth.df$X1 %in% 1:10))
  expect_true(all(smooth.df$X2 %in% 1:12))
  test.sim.data <- list(list(pars.normal = c(u = 1, v = 2, 
                                             w = 3), ar = array(1, c(2, 3, 4))))
  smooth.df <- convertSimResultsToDataFrame(test.sim.data, 
                                            "ar")
  expect_true(is.data.frame(smooth.df))
  expect_equal(dim(smooth.df), c(24, 7))
  expect_true(all(colnames(smooth.df) == c("u", "v", "w", "X1", 
                                           "X2", "X3", "sum.stat")))
})

test_that("test.fitGlm", {
  glms.fitted.csi <- fitGlm(sim.data.csi, jaatha.csi)
  expect_true(is.list(glms.fitted.csi))
  expect_equal(1, length(glms.fitted.csi))
  expect_equal(length(sim.data.csi[[1]]$poisson.vector), 
               length(glms.fitted.csi$poisson.vector))
  expect_true(all(sapply(glms.fitted.csi$poisson.vector, is)[1, ] == "glm"))
  glms.fitted.smooth <- fitGlm(smooth.sim.data, smooth.jaatha)
  expect_true(is.list(glms.fitted.smooth))
  expect_equal(1, length(glms.fitted.smooth))
  expect_true(is.list(glms.fitted.smooth$mat))
  expect_equal(1, length(glms.fitted.smooth$mat))
  expect_true("glm" %in% is(glms.fitted.smooth$mat[[1]]))
})

test_that("test.fitGlm.Smoothing", {
  glm.fitted <- fitGlm(smooth.sim.data, smooth.jaatha)
  expect_true(is.list(glm.fitted$mat))
  expect_equal(1, length(glm.fitted$mat))
  expect_true("glm" %in% is(glm.fitted$mat[[1]]))
})

test_that("test.fitGlm.Smoothing.border", {
  glm.fitted <- fitGlm(smooth.sim.data, smooth.border.jaatha)
  expect_true(is.list(glm.fitted$mat))
  expect_equal(2, length(glm.fitted$mat))
  expect_true("glm" %in% is(glm.fitted$mat[["smooth"]]))
  expect_true(is.list(glm.fitted$mat$border))
  expect_equal(40, length(glm.fitted$mat$border))
})

