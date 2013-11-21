test.confInts <- function() {
  jaatha <- Jaatha.confidenceIntervals(jaatha.csi, 0.95, 10, 2)
  estimates <- t(Jaatha.getLikelihoods(jaatha)[1, -(1:2), drop=FALSE])
  conf.ints <- jaatha@conf.ints
  checkTrue( all(conf.ints[,1] <= estimates & estimates <= conf.ints[,2]) )  
}
