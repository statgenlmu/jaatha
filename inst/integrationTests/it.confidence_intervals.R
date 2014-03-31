test.confidenceIntervals <- function() {
  jaatha <- Jaatha.confidenceIntervals(jaatha.csi, 0.95, 10, 1)
  estimates <- t(Jaatha.getLikelihoods(jaatha.csi)[1, -(1:2), drop=FALSE])
  conf.ints <- jaatha@conf.ints
  checkEquals( c(2,2), dim(conf.ints))
  checkTrue( all(conf.ints[,1] <= estimates & estimates <= conf.ints[,2]) )  
}
