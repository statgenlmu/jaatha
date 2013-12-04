# --------------------------------------------------------------
# Authors:  Paul R. Staab & Lisha Mathew
# Date:     2013-11-13
# Licence:  GPLv3 or later
# --------------------------------------------------------------

#' This function executes simulations for a number of random 
#' parameter values inside a block.
#' 
#' @param sim the number of simulations to execute.
#' @param block the block in which the simulations are executed.
#' @param jaatha The current jaatha object
#' @return A list of simulation results, where is simulation result is a list of
#'         summary statistics.
simulateWithinBlock <- function(sim, block, jaatha) {
  # Sample random simulation parameters
  sim.pars <- aperm(array(runif(getParNumber(jaatha)*sim,
                                 min=block@border[,1],
                                 max=block@border[,2]),
                           dim=c(getParNumber(jaatha),sim)))
  
  # Add the corners of the block to sim parameters
  sim.pars <- rbind(sim.pars, getCorners(block))

  runSimulations(sim.pars, jaatha@cores, jaatha)
}
