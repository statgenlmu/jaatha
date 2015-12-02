#' Jaatha: Simulation-based maximum likelihood parameter estimation
#' 
#' Jaatha is a simulation based maximum likelihood parameter estimation
#' methods.
#' 
#' @name jaatha-package
#' @docType package
#' @keywords package
#' @importFrom parallel mclapply
#' @importFrom assertthat assert_that is.string is.count
NULL

#' Simulation based maximum likelihood estimation
#' 
#' @param model The model used for the estimation. 
#'   See \code{\link{create_jaatha_model}}.
#' @param data The data used for the estimation.
#'   See \code{\link{create_jaatha_data}}.
#' @param repetitions The number of independend optimizations that will be
#'   conducted. You should use a value greater than one here, to minimize
#'   the chance that the algorithms is stuck in a local maximum.
#' @param sim The number of simulations conducted for each step.
#' @param max_steps The maximal number of steps, in case Jaatha fails to 
#'   converge.
#' @param init_method Determines how the starting position of each repetition
#'   is chosen. See below for a description of the different options. 
#' @param cores The number of CPU cores that will be used for the simulations.
#'   The relies on the \pkg{parallel} package, and consequenlty only one
#'   core is supported on Windows.
#' @param verbose If \code{TRUE}, information about the optimization algorithm
#'   is printed.
#' @return TBR
#' 
#' @section Algorithm:
#'   TBR
#'   
#' @section Initialization Methods:
#'   TBR
#'   
#' @author Paul Staab and Lisha Mathew
#' @export
jaatha <- function(model, data, 
                   repetitions = 3, 
                   sim = model$get_par_number() * 25, 
                   max_steps = 100, 
                   init_method = c("initial-search", "zoom-in", "middle"),
                   cores = 1,
                   verbose = TRUE) {
  
  # Check parameters
  assert_that(is_jaatha_model(model))
  assert_that(is_jaatha_data(data))
  assert_that(is.count(repetitions))
  assert_that(is.count(sim))
  assert_that(is.count(cores))
  
  # Setup
  sim_cache <- create_sim_cache()
  log <- create_jaatha_log(model, data, repetitions, sim, 
                           max_steps, init_method, verbose)
  
  # Get start positions
  log$log_initialization(init_method[1])
  start_pos <- get_start_pos(model, data, repetitions, sim, init_method, cores,
                             sim_cache = sim_cache)
  block_width <- 0.1
  
  for (rep in 1:repetitions) {
    estimate <- start_pos[rep, ]
    log$log_new_rep(rep, estimate)
    likelihood <- -Inf
    last_lh_increase <- 0
    
    for (step in 1:max_steps) {
      block <- create_block(cbind(estimate - block_width * .5,
                                  estimate + block_width * .5), 
                            cut = TRUE)
      
      local_ml <- estimate_local_ml(block, model, data, sim, cores, sim_cache)
      log$log_estimate(rep, step, local_ml)
      estimate <- local_ml$par
      
      if (local_ml$value > likelihood) {
        likelihood <- local_ml$value
        last_lh_increase <- step
      }
      
      if (step >= last_lh_increase + 10) {
        log$log_convergence(rep)
        break
      }
    }
  }
  
  # get presice llh values for best estimates
  log$log_llh_correction()
  best_values <- log$get_best_estimates(5)
  for (i in 1:nrow(best_values)) {
    llh <- estimate_llh(model, data, as.numeric(best_values[i, -(1:3)]), 
                        100, cores, TRUE)
    log$log_estimate("final", i, llh, best_values[i, 3])
  }
  
  # return the results
  log$create_results()
}
