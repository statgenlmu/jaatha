sampleFromModel <- function(x, y){
  return( c(rpois(10, x), rpois(10, y)) )
}

sim.func <- function(jaatha, sim.pars) {
  t(apply(sim.pars, 1, 
          function(x) sampleFromModel(x[1], x[2])))
}

test.customSimulator <- function() {
  set.seed(25)
  data.observed <- sampleFromModel(3, 5)
  par.ranges <- matrix(c(0.1, 0.1, 10, 10), 2, 2)
  rownames(par.ranges) <- c('x', 'y')
  colnames(par.ranges) <- c('min', 'max')

  jaatha <- new('Jaatha', sim.func, par.ranges, data.observed) 
  jaatha <- Jaatha.initialSearch(jaatha, 20, 2)
  jaatha <- Jaatha.refinedSearch(jaatha, 1, 20)
  estimates <- t(Jaatha.getLikelihoods(jaatha)[1, -(1:2), drop=FALSE])
  checkTrue( 2 <= estimates[1] & estimates[1] <= 4 )
  checkTrue( 4 <= estimates[2] & estimates[2] <= 6 )

  jaatha <- Jaatha.confidenceIntervals(jaatha, 0.95, 10, 2)

  estimates <- t(Jaatha.getLikelihoods(jaatha)[1, -(1:2), drop=FALSE])
  conf.ints <- jaatha@conf.ints
  checkTrue( all(conf.ints[,1] <= estimates & estimates <= conf.ints[,2]) )  
}
