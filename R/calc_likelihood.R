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
  sim.pars <- matrix(pars, sim, length(pars), byrow=TRUE)

  sim.data <- runSimulations(sim.pars, jaatha@cores, jaatha)

  logL <- 0
  # Average the values of each summary statistic
  for (sum.stat in jaatha@sum.stats) {
    if (sum.stat$method == "poisson.transformed") {
    values <- t(sapply(sim.data, function(x) sum.stat$transformation(x$jsfs))) 
    simSS <- apply(values, 2, mean)

    simSS[simSS==0] <- 0.5
    sum.stat.value <- sum.stat$transformation(sum.stat$value)
    logL <- logL + sum(sum.stat.value * log(simSS) - simSS - calcLogFactorial(sum.stat.value)) 
    }
  }
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
