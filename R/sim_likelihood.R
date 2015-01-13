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
  sim_pars <- matrix(pars, sim, length(pars), byrow=TRUE)
  sim_data <- runSimulations(sim_pars, jaatha@cores, jaatha)

  # Average the values of each summary statistic
  sum(sapply(jaatha@sum.stats, simLogLLH, sim_data, getScalingFactor(jaatha)))
}

simLogLLH <- function(sum_stat, ...) UseMethod("simLogLLH")
simLogLLH.default <- function(sum_stat, ...) stop('Unkown Summary Statistic')

simLogLLH.Stat_PoiInd <- function(sum_stat, sim_data, scaling_factor = 1) {
  values <- t(sapply(sim_data,
                     function(x) sum_stat$transform(x))) 
  simSS <- apply(values, 2, mean)
  
  simSS[simSS==0] <- 0.5
  if (scaling_factor != 1) simSS <- simSS * scaling_factor
  
  sum(sum_stat$get_data() * log(simSS) - simSS - calcLogFactorial(sum_stat$get_data())) 
}

simLogLLH.Stat_PoiSmooth <- function(sum_stat, sim_data, scaling_factor = 1) {
  sim_results <- sapply(sim_data,
                        function(x) sum_stat$transform(x)$sum.stat)
  sim_mean <- apply(sim_results, 1, mean)
  
  if (any(sim_mean == 0)) {
    warning(paste("A summary statistic was always 0, likelihood will",
                  "be inaccurate. Try increasing sim.final"))
    sim_mean[sim.mean == 0] <- .5
  }
  
  obs <- sum_stat$get_data()$sum.stat
  sum(obs * log(sim_mean) - sim_mean - calcLogFactorial(obs)) 
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
