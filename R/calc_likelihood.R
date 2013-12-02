# --------------------------------------------------------------
# calc_likelihood.R
# A function to calculate the likelihood of specific parsameter 
# combinations using simulations. 
# 
# Authors:  Paul R. Staab & Lisha Mathew 
# Date:     2013-11-28
# Licence:  GPLv3 or later
# --------------------------------------------------------------

## Funtion to calculate the likelihood based on simulations with the
## given parameters.  Order of parameters should be the same as needed
## for the simulate-function (in Simulator.R).
simLikelihood <- function(jaatha, sim, pars) {
  sim.pars <- matrix(pars, sim, length(pars), byrow=TRUE)
  sim.data <- runSimulations(sim.pars, jaatha@cores, jaatha)

  # Average the values of each summary statistic
  log.cl <- 0

  sum.stats <- jaatha@sum.stats
  for (sum.stat in names(sum.stats)) {
    if (sum.stats[[sum.stat]]$method %in% c("poisson.transformed", "poisson.independent")) {
      values <- t(sapply(sim.data, 
                         function(x) sum.stats[[sum.stat]]$transformation(x[[sum.stat]]))) 
      simSS <- apply(values, 2, mean)

      simSS[simSS==0] <- 0.5
      sum.stat.value <- sum.stats[[sum.stat]]$value.transformed
      log.cl <- log.cl + 
        sum(sum.stat.value * log(simSS) - simSS - calcLogFactorial(sum.stat.value)) 
    }
    else if (sum.stats[[sum.stat]]$method == "poisson.smoothing") {
      sum.stat.value <- as.vector(sum.stats[[sum.stat]]$value)
      sim.mean <- apply(sapply(sim.data, function(x) as.vector(x[[sum.stat]])), 1, mean)
      if (any(sim.mean == 0)) warning("Warning: A summary statistic was always 0,
                                      likelihood will be inaccurate. Try
                                      increasing sim.final") 
      sim.mean[sim.mean == 0] <- .5
      log.cl <- log.cl + 
        sum(sum.stat.value * log(sim.mean) - sim.mean - calcLogFactorial(sum.stat.value)) 
    }
    else stop("Unsupported SumStat method")
  }
  return(log.cl)
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
