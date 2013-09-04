# --------------------------------------------------------------
# calc_likelihood.R
# A function to calculate the likelihood of specific parsameter 
# combinations using simulations. 
# 
# Authors:  Paul R. Staab & Lisha Mathew 
# Date:     2013-09-04
# Licence:  GPLv3 or later
# --------------------------------------------------------------

## Funtion to calculate the likelihood based on simulations with the
## given parameters.  Order of parameters should be the same as needed
## for the simulate-function (in Simulator.R).
calcLikelihood <- function(jaatha, sim, pars){
  .log2("Called Jaatha.calcLikelihood()")
  pars <- matrix(rep(pars, each=sim), nrow=sim)

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
  logL <- sum(jaatha@sumStats * log(simSS) - simSS - calcLogFactorial(jaatha@sumStats))
  .log2("Finished Jaatha.calcLikelihood(). Return:",logL)
  return(logL)
}


calcLogFactorial <- function(k) {
  if (!isJaathaVariable("logfacs")) setJaathaVariable("logfacs", c(0)) 
  logfacs <- getJaathaVariable("logfacs")

  maxk <- max(k)
  if (maxk > length(logfacs)) {
    l <- length(logfacs) + 1
    logfacs[l:maxk] <- 0
    for (i in l:maxk) {
      logfacs[i] <- logfacs[i-1] + log(i)
    }
    setJaathaVariable("logfacs", logfacs)
  }

  ret <- rep(0, length(k))
  ret[k!=0] <- logfacs[k[k!=0]] 
  return(ret)
}
