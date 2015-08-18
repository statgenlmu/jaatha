calc_poisson_llh <- function(data, stat, loglambda, scaling_factor = 1) {
  # Upscale predicted expectation value if we use scaling
  if (scaling_factor != 1) loglambda <- loglambda + log(scaling_factor)
  
  sum(data$get_values(stat) * loglambda - 
        exp(loglambda) - data$get_log_factorial(stat)) 
}


approximate_llh <- function(x, data, param, glm_fitted, ...) {
  "approximates the log-likelihood using the fitted glms"
  UseMethod("approximate_llh")
}

#' @export
approximate_llh.default <- function(x, data, param, glm_fitted, ...) { 
  stop("Unknown Summary Statistic")
}


#' @export
approximate_llh.jaatha_model <- function(x, data, param, glm_fitted) {
  assert_that(is_jaatha_data(data))
  assert_that(is.numeric(param))
  assert_that(is.list(glm_fitted))
  sum(vapply(x$get_sum_stats(), approximate_llh, numeric(1),
             data, param, glm_fitted, x$get_scaling_factor()))
}


#' @importFrom stats predict.glm
#' @export
approximate_llh.jaatha_stat_basic  <- function(x, data, param, glm_fitted, 
                                               scaling_factor) {
  
  loglambda <- sapply(glm_fitted[[x$get_name()]], function(glm_obj) {
    glm_obj$coefficients %*% c(1, param)
  })
  
  # Calculate the Poission log-likelihood
  calc_poisson_llh(data, x, loglambda, scaling_factor)
}


#' @importFrom stats optim
optimize_llh <- function(block, model, data, glms) {
  best_value <- optim(block$get_middle(),
                      function(param) {
                        approximate_llh(model, data, param, glms)
                      },
                      lower = block$get_border()[ , 1, drop = FALSE], 
                      upper = block$get_border()[ , 2, drop = FALSE],
                      method = "L-BFGS-B", 
                      control = list(fnscale = -1))
  
  assert_that(block$includes(best_value$par))
  best_value
}


estimate_local_ml <- function(block, model, data, sim, cores, sim_cache) {
  sim_data <- model$simulate(pars = block$sample_pars(sim, add_corners = TRUE), 
                             data = data,
                             cores = cores)
  
  # Cache simulation & load older simulations within this block
  sim_cache$add(sim_data)
  sim_data <- sim_cache$get_sim_data(block)
  assert_that(length(sim_data) >= sim)
  
  # Fit glms and find maximal likelihood value
  glms <- fit_glm(model, sim_data)
  optimize_llh(block, model, data, glms)
}


#' Estimate the Log-Likelihood for a given parameter combination
#' 
#' This function estimates the Log-likelihood value for a given
#' parameter combination. It conducts a number of simulations for
#' the parameter combination, averages the summary statistics to
#' esimate their expected values, and uses them to calculate the
#' likelihood. For a resonable number of simulation, this is more
#' precise than the glm fitting used in the main algorithm.
#' 
#' @inheritParams jaatha
#' @param parameter The parameter combination for which the loglikelihood
#'          will be estimated.
#' @param sim The number of simulations that will be used for averaging the
#'          expectation values of the summary statistics.
#' @param normalized For internal use. Indicates whether the parameter
#'          combination is normalized to [0, 1]-scale, or on its natural
#'          scale.
#' @export
estimate_llh <- function(model, data, parameter, sim = 100, 
                         cores = 1, normalized = FALSE) {
  
  assert_that(is_jaatha_model(model))
  assert_that(is_jaatha_data(data))
  assert_that(is.numeric(parameter))
  assert_that(is_positive_int(sim))
  assert_that(is_positive_int(cores))
  assert_that(is_single_logical(normalized))
  
  if (!normalized) parameter <- model$get_par_ranges()$normalize(parameter)
  sim_pars <- matrix(parameter, sim, length(parameter), byrow = TRUE)
  sim_data <- model$simulate(sim_pars, data, cores)
  
  llh <- sum(vapply(names(model$get_sum_stats()), function(stat){
    stat_values <- sapply(sim_data, function(x) x[[stat]])
    if (!is.matrix(stat_values)) stat_values <- matrix(stat_values, nrow = 1)
    log_means <- log(rowMeans(stat_values))
    calc_poisson_llh(data, stat, log_means, model$get_scaling_factor())
  }, numeric(1)))
  
  list(param = parameter,
       value = llh)
}
