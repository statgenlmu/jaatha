test.customSimulator <- function() {
  checkTrue( length(jaatha.csi@starting.positions) == 4 )

  estimates <- Jaatha.getLikelihoods(jaatha.csi)
  checkTrue( all( estimates[,1] != 0 ) )
  checkTrue( all( 0.1 <= estimates[,3] & estimates[,4] <= 10 ) )
}
