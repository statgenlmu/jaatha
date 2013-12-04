test.simLikelihood <- function() {
  lh1 <- simLikelihood(jaatha.csi, 10, c(.5,.5))   
  checkTrue(is.numeric(lh1))
  checkTrue(lh1 != 0)
  lh2 <- simLikelihood(smooth.jaatha, 20, c(.5,.5))
  checkTrue(is.numeric(lh2))
  checkTrue(lh2 != 0)
}
