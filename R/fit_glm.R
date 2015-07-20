fit_glm <- function(x, sim_data, ...) UseMethod("fit_glm")
fit_glm.default <- function(x, sim_data, ...) {
  stop("Unknown Summary Statistic")
}


fit_glm.jaatha_model <- function(x, sim_data, ...) { 
  "Fits a GLM to the simulation results"
  glm_fitted <- list()
  lapply(x$get_sum_stats(), fit_glm, sim_data, ...)
}


fit_glm.jaatha_stat_basic <- function(x, sim_data, ...) {
  "Fits a GLM for each entry of the simulation results"
  Y <- do.call(rbind, lapply(sim_data, function(data) data[[x$get_name()]]))
  X <- do.call(rbind, lapply(sim_data, function(data) data$pars_normal))
  
  glms <- lapply(1:ncol(Y), function(i) {
    glm.fit(X, Y[ , i], family=poisson("log"), control = list(maxit = 200))
  })
  
  sapply(glms, function(x){if (!x$converged) stop("GLM did not converge")})
  glms
}


#' Fits a GLM for a summary statistics of type "poisson.smoothed"
#'
#' @param sim_data Results from simulations
#' @param sum_stat Name of the summary statistics
#' @return A list with one fitted GLM
fit_glm.Stat_PoiSmooth <- function(sum_stat, sim_data) {
  par_names <- names(sim_data[[1]]$pars)
  model <- paste0("sum.stat ~ ",
                  "(", sum_stat$get_model(), ")",  
                  "*(", paste(par_names, collapse="+"), ")") 

  sim_data_df <- do.call(rbind, lapply(sim_data, function(sim_result) {
    pars <- matrix(sim_result$pars.normal, 1,
                   length(sim_result$pars.normal), byrow=TRUE)
    colnames(pars) <- names(sim_result$pars.normal)
    data.frame(pars, sim_result[[sum_stat$get_name()]])
  }))

  smooth_glm  <- glm(model, data=sim_data_df, family=poisson("log"), 
                     model = FALSE, x = FALSE, y = FALSE,
                     control = list(maxit = 200))
  if (!smooth_glm$converged) stop("GLM did not converge")
  
  list(smooth=smooth_glm)
}
