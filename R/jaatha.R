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
                   sim = model$get_par_ranges()$get_par_number() * 25, 
                   max_steps = 200, 
                   initialization = c("zoom-in", "initial-search", "middle"),
                   cores = 1) {
  
  # Check parameters
  assert_that(is_jaatha_model(model))
  assert_that(is_jaatha_data(data))
  assert_that(is_positive_int(repetitions))
  assert_that(is_positive_int(sim))
  assert_that(is_positive_int(cores))
  assert_that(any(initialization[1] == c("zoom-in", "initial-search", "middle")))
  
  # Setup
  likelihood_table <- create_likelihood_table(jaatha, max.steps)
  
  
  
}