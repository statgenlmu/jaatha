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
  seeds <- sampleSeed(length(pars)+1)

  sim.data <- mclapply(1:nrow(pars), runSimulation, pars=pars, 
                       seeds=seeds, jaatha=jaatha,
                       mc.preschedule=TRUE, mc.cores=cores)

  if ( "try-error" %in% unlist(sapply(sim.data, is)) ) {
    cat("Error(s) while simulating:\n")
    sapply(sim.data["try-error" %in% unlist(sapply(sim.data, is))],
           function(x) print(attr(x, 'condition')))
    stop("One or more errors while simulating. Check your sim.func")
  }
  
  set.seed(seeds[length(seeds)])
  sim.data
}

runSimulation <- function(i, pars, seeds, jaatha) {
  # Set the seed & prepare parameters
  set.seed(seeds[i])
  sim_pars <- denormalize(pars[i, ], jaatha)
  
  # Simulate
  sim_results <- jaatha@simFunc(sim_pars, jaatha)
  
  # Calculate Summary Statistics
  sim_sum_stats <- lapply(jaatha@sum_stats, function(sum_stat) {
    sum_stat$transform(sim_results)
  })
  
  # Add the parameter values
  sim_sum_stats$pars <- sim_pars
  sim_sum_stats$pars.normal <- pars[i, ]
  
  sim_sum_stats
}


test_simulation <- function(jaatha, quite=FALSE) {
  pars <- matrix(rep(0.5, nrow(jaatha@par.ranges)), nrow = 1)
  
  time <- system.time(
    a <- runSimulation(1, pars, sampleSeed(1), jaatha)
  )['elapsed']

  if (time > 30) warning('Each simulation takes about ', round(time),
                         's, Jaatha might run for a long time.')
  if (!quite && time <= 30) {
    if (time < 1) message('A simulation takes less than a second')
    else message('A simulation takes about ', round(time), 's')
  }
  
  invisible(NULL)
}
