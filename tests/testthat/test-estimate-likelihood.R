context("Likelihood estimation")

test_that("test.estimateLikelihood.PoiSmoothing", {

  glm.fitted <- fitGlm(smooth.sim.data, smooth.jaatha)
  ll <- estimateLogLikelihood(c(x = 3, y = 3), glm.fitted, 
                              smooth.jaatha@sum.stats)
  expect_true(is.numeric(ll))
  expect_true(0 <= exp(ll) | exp(ll) <= 1)
  glm.fitted <- fitGlm(smooth.sim.data, smooth.border.jaatha)
  ll2 <- estimateLogLikelihood(c(x = 3, y = 3), glm.fitted, 
                               smooth.border.jaatha@sum.stats)
  expect_true(is.numeric(ll2))
  expect_true(0 <= exp(ll2) | exp(ll2) <= 1)
  expect_true(ll != ll2)
})

test_that("test.estimateLikelihood.PoiTransformed", {
  sim.data.csi <- jaatha:::simulateWithinBlock(10, block.test, jaatha.csi)
  glm.fitted <- fitGlm(sim.data.csi, jaatha.csi)
  ll <- estimateLogLikelihood(c(x = 1, y = 1), glm.fitted, 
                              jaatha.csi@sum.stats)
  expect_true(is.numeric(ll))
  expect_true(0 <= exp(ll) | exp(ll) <= 1)
})

