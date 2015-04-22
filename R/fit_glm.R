fitGlm <- function(sum_stat, sim_data, ...) UseMethod("fitGlm")
fitGlm.default <- function(sum_stat, sim_data...) stop("Unkown Summary Statistic")

#' Fits a Generalized Linear Model within a block.  
#'
#' The model describes how the summary statistics depend on the parameters.
#' 
#' @param sim_data Simulation data in this block.
#' @param jaatha A Jaatha Object.
#' @return A list containing a list of fitted GLMs for each summary
#' statistic.
fitGlm.Jaatha <- function(jaatha, sim_data) { 
  glm_fitted <- list()
  for (i in seq(along = jaatha@sum_stats)) {
    name <- names(jaatha@sum_stats)[i]
    glm_fitted[[name]] <- fitGlm(jaatha@sum_stats[[name]], sim_data)
  }
  glm_fitted
}


#' Fits a GLM for a summary statistics of type "poisson.independent"
#'
#' @param sim_data Results from simulations
#' @param sum_stat Name of the summary statistics
#' @return A list of fitted GLMs, one for each function
fitGlm.Stat_PoiInd <- function(sum_stat, sim_data) { 
  stat_sim <- t(sapply(sim_data, 
                       function(data) c(data$pars.normal, data[[sum_stat$get_name()]])))
  
  par_names <- names(sim_data[[1]]$pars)
  stat_names <- paste("S", 1:(ncol(stat_sim)-length(par_names)), sep="")
  colnames(stat_sim) <- c(par_names, stat_names)
  
  formulas <- paste0(stat_names, "~", paste(par_names ,collapse= "+"))
  glms <- lapply(formulas, glm, data=data.frame(stat_sim), family=poisson,
                 model = FALSE, x = FALSE, y = FALSE, control = list(maxit = 200))
  sapply(glms, function(x){if (!x$converged) stop("GLM did not converge")})
  glms
}


#' Fits a GLM for a summary statistics of type "poisson.smoothed"
#'
#' @param sim_data Results from simulations
#' @param sum_stat Name of the summary statistics
#' @return A list with one fitted GLM
fitGlm.Stat_PoiSmooth <- function(sum_stat, sim_data) {
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
