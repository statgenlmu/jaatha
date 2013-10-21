# --------------------------------------------------------------
# simulate_within_block.R
# Jaatha main simulation function and related helpers 
# 
# Authors:  Paul R. Staab & Lisha Mathew
# Date:     2013-09-04
# Licence:  GPLv3 or later
# --------------------------------------------------------------

## This function simulates nSamp random parameter combinations within
## lower and upper bounds for each parameter (boarders values
## needed to be between 0 and 1!!) and with the demographic
## model specified in simulate() with theta=5.  Returns the
## random parameters and the corresponding summary statistics. 
simulateWithinBlock<- function(sim, block, jaatha) {
  # Sample random simulation parameters
  randompar <- aperm(array(runif(jaatha@nPar*sim,
                                 min=block@border[,1],
                                 max=block@border[,2]),
                           dim=c(jaatha@nPar,sim)))
  
  # Add the corners of the block to sim parameters
  randompar <- rbind(randompar, getCorners(block))

  # Create "packages" of parameters combinations for possible parallelization.
  sim.packages <- createSimulationPackages(randompar, jaatha@sim.package.size)
  seeds <- generateSeeds(length(sim.packages)+1)

  i <- NULL # To make R CMD check stop complaining
  # Simulate each package, maybe on different cores
  sum.stats <- foreach(i = seq(along = sim.packages), .combine='rbind') %dopar% {
    set.seed(seeds[i])
    sim.pars <- .deNormalize(jaatha, 
                             sim.packages[[i]])
    sumStats <- jaatha@simFunc(jaatha, sim.pars)
    return(sumStats)
  }

  set.seed(seeds[length(seeds)])

  # Create combined output
  sim.result <- data.frame(cbind(randompar, sum.stats))
  colnames(sim.result) <- c(jaatha@par.names, paste("SS", 1:ncol(sum.stats),
                                                    sep=""))  
  .log2("Finished simulating for this block")
  return(sim.result)
}

createSimulationPackages <- function(random.par, package.size) {
  if (package.size == 0) return(list(random.par))

  sim.packages <- list()
  num.pars <- nrow(random.par)
  i <- 0
  
  while (i < num.pars/package.size) {
    i <- i + 1
    lower <- (i-1)*package.size+1
    upper <-  min(i*package.size, num.pars)
    sim.packages[[i]] <- random.par[lower:upper, , drop=F]
  }

  return(sim.packages)  
}

