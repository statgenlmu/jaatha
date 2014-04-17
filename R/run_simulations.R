# --------------------------------------------------------------
# Authors:  Paul R. Staab
# Date:     2013-11-14
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
  stopifnot(ncol(pars) == getParNumber(jaatha))
  stopifnot(all( 0-1e-5 <= pars & pars <= 1 + 1e-5 ))
  colnames(pars) <- getParNames(jaatha)
  seeds <- generateSeeds(length(pars)+1)

  sim.data <- mclapply(1:nrow(pars), runSimulation, pars=pars, 
                       seeds=seeds, jaatha=jaatha,
                       mc.preschedule=TRUE, mc.cores=cores)

  if ( "try-error" %in% unlist(sapply(sim.data, is)) ) 
    stop("Error while simulating. Check your sim.func") 

  set.seed(seeds[length(seeds)])
  return(sim.data)
}

runSimulation <- function(i, pars, seeds, jaatha) {
  set.seed(seeds[i]) 
  sim.pars <- denormalize(pars[i, ], jaatha)
  sim.results <- jaatha@simFunc(sim.pars, jaatha)
  sim.results$pars <- sim.pars
  sim.results$pars.normal <- pars[i, ]
  return(sim.results)
}
