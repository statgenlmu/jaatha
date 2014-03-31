test.estimateLikelihood.PoiTransformed <- function() {
  glm.fitted <- fitGlm(sim.data.csi, jaatha.csi)
  ll <- estimateLogLikelihood(c(x=1, y=1), glm.fitted, jaatha.csi@sum.stats)
  checkTrue( is.numeric(ll) )
  checkTrue( 0 <= exp(ll) | exp(ll) <= 1 )
}

test.estimateLikelihood.PoiSmoothing <- function() {
  glm.fitted <- fitGlm(smooth.sim.data, smooth.jaatha)
  ll <- estimateLogLikelihood(c(x=3, y=3), glm.fitted, smooth.jaatha@sum.stats)
  checkTrue( is.numeric(ll) )
  checkTrue( 0 <= exp(ll) | exp(ll) <= 1 )

  glm.fitted <- fitGlm(smooth.sim.data, smooth.border.jaatha)
  ll2 <- estimateLogLikelihood(c(x=3, y=3), glm.fitted, smooth.border.jaatha@sum.stats)
  checkTrue( is.numeric(ll2) )
  checkTrue( 0 <= exp(ll2) | exp(ll2) <= 1 )
  checkTrue( ll != ll2 )
}
