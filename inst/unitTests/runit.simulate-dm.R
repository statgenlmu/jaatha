runit.simulateDemographicModel <- function() {
  dm <- dm.createThetaTauModel(c(20, 23), 100)
  jsfs <- matrix(rpois(21*24, 2.5), 21, 24)
  jaatha <- Jaatha.initialize(dm, jsfs=jsfs)
  sim.pars <- matrix(c(1, 5, 2, 7), 2, 2)

  sim <- simulateDemographicModel(jaatha, sim.pars)
  checkTrue( is.matrix(sim) )
  checkTrue( nrow(sim) == 2 )
  checkTrue( sum(sim) > 0 )
}
