context("Likelihood estimation")

test_that("calcStatLLH works for Stat_PoiInd", {
  glm_fitted <- fitGlm(csi.sum.stat, sim.data.csi)
  
  ll <- calcStatLLH(csi.sum.stat, glm_fitted, c(x=.5, y=.5), scaling_factor=1)
  expect_true(is.numeric(ll))
  expect_true(0 <= exp(ll) | exp(ll) <= 1)
  
  ll2 <- calcStatLLH(csi.sum.stat, glm_fitted, c(x=.5, y=.5), scaling_factor=2)
  expect_true(is.numeric(ll))
  expect_true(0 <= exp(ll) | exp(ll) <= 1)
  expect_true(ll != ll2)
})


test_that("calcStatLLH works for Stat_PoiSmooth", {
  glm_fit <- fitGlm(smooth_stat, smooth_sim_data)
  
  ll <- calcStatLLH(smooth_stat, glm_fit, c(x=.5, y=.5), scaling_factor=1)
  expect_true(is.numeric(ll))
  expect_true(0 <= exp(ll) | exp(ll) <= 1)
})


test_that("estimateLogLikelihood works", {
  glm_fit <- fitGlm(jaatha.csi, sim.data.csi)
  ll <- estimateLogLikelihood(c(x=3, y=3), glm_fit, jaatha.csi@sum.stats)
  expect_true(is.numeric(ll))
  expect_true(0 <= exp(ll) | exp(ll) <= 1)
  
  glm_fit <- fitGlm(csi.sum.stat, sim.data.csi)
  ll2 <- calcStatLLH(csi.sum.stat, glm_fit, c(x=3, y=3), scaling_factor = 1)
  expect_equal(ll, ll2)
})
