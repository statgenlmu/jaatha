# --------------------------------------------------------------
# DataProcessor.R
# Various functions for importing and processing data from external 
# sources.
# 
# Authors:  Lisha Naduvilezhath & Paul R. Staab
# Date:     2012-10-05
# Licence:  GPLv3 or later
# --------------------------------------------------------------

## Funtion to calculate the likelihood based on simulations with the
## given parameters.  Order of parameters should be the same as needed
## for the simulate-function (in Simulator.R).
Jaatha.calcLikelihood <- function(jObject, nSimulations, par){
    i <- NULL   # To make R CMD check stop complaining

	.log2("Called Jaatha.calcLikelihood()")
    par <- matrix(rep(par, each=nSimulations), nrow=nSimulations)

    sim.packages <- createSimulationPackages(par, jObject@sim.package.size)
    seeds <- generateSeeds(length(sim.packages)+1)

    # Simulate each package, maybe on different cores
    simSS  <- foreach(i = seq(along = sim.packages), .combine='rbind') %dopar% {
      set.seed(seeds[i])
      sim.pars <- .deNormalize(jObject, sim.packages[[i]])
      sumStats <- dm.simSumStats(jObject@dm, sim.pars, jObject@sum.stats.func)
      return(sumStats)
    }
    set.seed(seeds[length(seeds)])

    # Average the values of each summary statistic
    simSS <- apply(simSS, 2, mean)

	.log2("Calculating Likelihood...")
	logL <- 0
	simSS[simSS==0] <- 0.5
	for (s in 1:jObject@nTotalSumstat){
		logL <- logL + jObject@sumStats[s] * log(simSS[s]) - simSS[s] - .logfac(jObject@sumStats[s])
	}
	.log2("Finished Jaatha.calcLikelihood(). Return:",logL)
	return(logL)
}


## Returns the logarithm of the factorial of k. Recursively implemented.
.logfac <- function(k) {
	if(k<2) return(0)
	if(k>10) return(lgamma(k+1))
	else return(log(k)+.logfac(k-1))
}
