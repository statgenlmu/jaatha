# --------------------------------------------------------------
# calc_likelihood.R
# A function to calculate the likelihood of specific parsameter 
# combinations using simulations. 
# 
# Authors:  Lisha Naduvilezhath & Paul R. Staab
# Date:     2013-09-04
# Licence:  GPLv3 or later
# --------------------------------------------------------------

## Funtion to calculate the likelihood based on simulations with the
## given parameters.  Order of parameters should be the same as needed
## for the simulate-function (in Simulator.R).
calcLikelihood <- function(jaatha, nSimulations, pars){
  .log2("Called Jaatha.calcLikelihood()")
  pars <- matrix(rep(pars, each=nSimulations), nrow=nSimulations)

  sim.packages <- createSimulationPackages(pars, jaatha@sim.package.size)
  seeds <- generateSeeds(length(sim.packages)+1)

  # Simulate each package, maybe on different cores
  i <- NULL   # To make R CMD check stop complaining
  simSS  <- foreach(i = seq(along = sim.packages), .combine='rbind') %dopar% {
    set.seed(seeds[i])
    sim.pars <- .deNormalize(jaatha, sim.packages[[i]])
    sumStats <- jaatha@simFunc(jaatha, sim.pars)
    return(sumStats)
  }
  set.seed(seeds[length(seeds)])

  # Average the values of each summary statistic
  simSS <- apply(simSS, 2, mean)

  .log2("Calculating Likelihood...")
  simSS[simSS==0] <- 0.5
  logL <- jaatha@sumStats * log(simSS) - simSS - .logfac(jaatha@sumStats)
  .log2("Finished Jaatha.calcLikelihood(). Return:",logL)
  return(logL)
}


## Returns the logarithm of the factorial of k. Recursively implemented.
.logfac <- function(k) {
  return( log(gamma( k+1 )) )
}
