# --------------------------------------------------------------
# simulate_within_block.R
# Jaatha main simulation function and related helpers 
# 
# Authors:  Paul R. Staab & Lisha Mathew
# Date:     2013-11-07
# Licence:  GPLv3 or later
# --------------------------------------------------------------

simulateWithinBlock <- function(sim, block, jaatha) {
  # Sample random simulation parameters
  sim.pars <- aperm(array(runif(jaatha@nPar*sim,
                                 min=block@border[,1],
                                 max=block@border[,2]),
                           dim=c(jaatha@nPar,sim)))
  
  # Add the corners of the block to sim parameters
  sim.pars <- rbind(sim.pars, getCorners(block))

  runSimulations(sim.pars, jaatha@cores, jaatha)
}


runSimulations <- function(pars, cores, jaatha) {
  seeds <- generateSeeds(length(pars)+1)

  if (cores == 1) {
    sum.stats <- lapply(1:nrow(pars), 
                        function(i) {
                          set.seed(seeds[i]) 
                          sim.pars <- denormalize(pars[i, ], jaatha)
                          jaatha@simFunc(jaatha, sim.pars)
                        })
  } else {
    sum.stats <- mclapply(1:nrow(pars), 
                          function(i) {
                            set.seed(seeds[i]) 
                            sim.pars <- denormalize(pars[i, ], jaatha)
                            jaatha@simFunc(jaatha, sim.pars)
                          },
                          mc.preschedule=TRUE, mc.cores=cores)
  }

  set.seed(seeds[length(seeds)])
  return(sum.stats)
}
