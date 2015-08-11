#' Jaatha: Simulation-based maximum likelihood parameter estimation
#' 
#' Jaatha is a simulation based maximum likelihood parameter estimation
#' methods.
#' 
#' @name jaatha-package
#' @docType package
#' @keywords package
#' @importFrom parallel mclapply
#' @importFrom assertthat assert_that
NULL

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
  assert_that(is_positive_int(repetitions))
  assert_that(is_positive_int(sim))
  assert_that(is_positive_int(cores))
  
  # Setup
  sim_cache <- create_sim_cache()
  log <- create_jaatha_log(repetitions, max_steps, model, verbose)
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
        log$log_convergence()
        break
      }
    }
  }
}
