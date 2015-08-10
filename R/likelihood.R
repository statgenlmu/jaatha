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
  
  # Upscale predicted expectation value if we use scaling
  if (scaling_factor != 1) loglambda <- loglambda + log(scaling_factor)
  
  # Calculate the Poission log-likelihood
  sum(data$get_values(x) * loglambda - 
        exp(loglambda) - data$get_log_factorial(x)) 
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
  
  glms <- fit_glm(model, sim_data)
  optimize_llh(block, model, data, glms)
}
