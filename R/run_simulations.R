# --------------------------------------------------------------
# Authors:  Paul R. Staab
# Date:     2013-11-13
# Licence:  GPLv3 or later
# --------------------------------------------------------------

#' The function that actually executes the simulations from within Jaatha.
#' 
#' @param pars A matrix with parameter values for the simulations. Each row is a
#'         set of model parameters. Should be in jaatha internal 0-1 scaling.
#' @param cores The number of cores to use. Using more than one core requires
#'         the 'multicore' package.
#' @param jaatha The current jaatha object
#' @return A list, where each entry is a list of summary statistics for a 
#'         simulation.
runSimulations <- function(pars, cores, jaatha) {
  checkType(pars, c('mat', 'num'))
  stopifnot(ncol(pars) == jaatha@nPar)
  colnames(pars) <- jaatha@par.names
  seeds <- generateSeeds(length(pars)+1)

  if (cores == 1) {
    sum.stats <- lapply(1:nrow(pars), runSimulation, pars=pars, 
                        seeds=seeds, jaatha=jaatha) 
  } else {
    sum.stats <- mclapply(1:nrow(pars), runSimulation, pars=pars, 
                          seeds=seeds, jaatha=jaatha,
                          mc.preschedule=TRUE, mc.cores=cores)
  }

  set.seed(seeds[length(seeds)])
  return(sum.stats)
}

runSimulation <- function(i, pars, seeds, jaatha) {
  set.seed(seeds[i]) 
  sim.pars <- denormalize(pars[i, ], jaatha)
  sim.results <- jaatha@simFunc(jaatha, sim.pars)
  sim.results$pars <- sim.pars
  sim.results$pars.normal <- pars[i, ]
  return(sim.results)
}
