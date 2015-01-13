# --------------------------------------------------------------
# Functions to fit GLMs to learn how the summary statistics depend
# on the parameters.
# --------------------------------------------------------------

fitGlm <- function(sum_stat, sim_data, ...) UseMethod("fitGlm")
fitGlm.default <- function(sum_stat, sim_data...) stop('Unkown Summary Statistic')

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
  for (i in seq(along = jaatha@sum.stats)) {
    name <- names(jaatha@sum.stats)[i]
    glm_fitted[[name]] <- fitGlm(jaatha@sum.stats[[i]], sim_data)
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
                       function(data) c(data$pars.normal, sum_stat$transform(data))))
  
  par_names <- names(sim_data[[1]]$pars)
  stat_names <- paste("S", 1:(ncol(stat_sim)-length(par_names)), sep="")
  colnames(stat_sim) <- c(par_names, stat_names)
  
  formulas <- paste0(stat_names, "~", paste(par_names ,collapse= "+"))
  glms <- lapply(formulas, glm, data=data.frame(stat_sim), family=poisson,
                 model = FALSE, x = FALSE, y = FALSE, control = list(maxit = 200))
  sapply(glms, function(x){if (!x$converged) stop('GLM did not converge')})
  glms
}


#' Fits a GLM for a summary statistics of type "poisson.smoothed"
#'
#' @param sim_data Results from simulations
#' @param sum_stat Name of the summary statistics
#' @param jaatha A Jaatha Object.
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
    data.frame(pars,  sum_stat$transform(sim_result))
  }))

  smooth_glm  <- glm(model, data=sim_data_df, family=poisson("log"), 
                     model = FALSE, x = FALSE, y = FALSE,
                     control = list(maxit = 200))
  if (!smooth_glm$converged) stop('GLM did not converge')
  
  #if (!is.null(jaatha@sum.stats[[sum.stat]]$border.transformation)) {
  #  glms <- list(smooth=smooth.glm,
  #               border=fitGlmPoiTransformed(sim.data, sum.stat,
  #                      jaatha@sum.stats[[sum.stat]]$border.transformation,
  #                      weighting, jaatha))
  #} else { 
    glms <- list(smooth=smooth_glm)
  #}
  glms
}


#' Converts simulation results into a data frame that is usable for fitting a
#' glm.
#'
#' Currently only works with nx2 matix summary statistics.
#' 
#' @param sim.data Results from simulations
#' @param sum.stat Name of the summary statistics which should get converted
#' @param mask Boolean vector of positions to exclude in the data.frame
#' @return The summary statistics as data.frame 
#' @importFrom reshape2 melt
convertSimResultsToDataFrame <- function(sim_data, sum_stat, mask=NULL) {
  #do.call(rbind, lapply(sim_data, function(sim_result) {

  dim_names <- lapply(dim(sim_data), function(x) 1:x)
  names(dim_names) <- paste0('X', 1:length(dim(sim_data)))
  dimnames(sim_data) <- dim_names
  sum_stat_df <- melt(data, value.name = 'sum.stat')

  # Add the corsponding parameters
  #if (!is.null(sim_result$pars.normal)) {
  #  pars <- matrix(sim_result$pars.normal, nrow(sum_stat_df),
  #                 length(sim_result$pars.normal), byrow=TRUE)
  #  colnames(pars) <- names(sim_result$pars.normal)
  #  sum_stat_df <- data.frame(pars, sum_stat_df)
  #}

  # Remove masked values (if any)
  if (!is.null(mask)) sum_stat_df <- sum_stat_df[!mask, ]
  sum_stat_df
  #}))
}
