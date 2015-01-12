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
#' @param sim.data Results from simulations
#' @param sum.stat Name of the summary statistics
#' @param weighting Potentially weights for the simulations.
#' @param jaatha A Jaatha Object.
#' @return A list with one fitted GLM
fitPoiSmoothed <- function(sim.data, sum.stat, weighting, jaatha) {
  model <- paste0("sum.stat ~ ",
                  "(", jaatha@sum.stats[[sum.stat]]$model, ")",  
                  "*(", paste(getParNames(jaatha), collapse="+"), ")") 

  sim.data.df <- convertSimResultsToDataFrame(sim.data, sum.stat,
                                              jaatha@sum.stats[[sum.stat]]$border.mask)

  smooth.glm  <- glm(model, data=sim.data.df, family=poisson("log"), 
                     model = FALSE, x = FALSE, y = FALSE,
                     control = list(maxit = 200))
  if (!smooth.glm$converged) stop('GLM did not converge')
  
  if (!is.null(jaatha@sum.stats[[sum.stat]]$border.transformation)) {
    glms <- list(smooth=smooth.glm,
                 border=fitGlmPoiTransformed(sim.data, sum.stat,
                        jaatha@sum.stats[[sum.stat]]$border.transformation,
                        weighting, jaatha))
  } else { 
    glms <- list(smooth=smooth.glm)
  }
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
convertSimResultsToDataFrame <- function(sim.data, sum.stat, mask=NULL) {
  do.call(rbind, lapply(sim.data, function(sim.result) {
    # Convert array to data.frame
    dim.names <- lapply(dim(sim.result[[sum.stat]]), function(x) 1:x)
    names(dim.names) <- paste0('X', 1:length(dim(sim.result[[sum.stat]])))
    dimnames(sim.result[[sum.stat]]) <- dim.names
    sum.stat.df <- melt(sim.result[[sum.stat]], value.name = 'sum.stat')

    # Add the corsponding parameters
    pars <- matrix(sim.result$pars.normal, nrow(sum.stat.df),
                   length(sim.result$pars.normal), byrow=TRUE)
    colnames(pars) <- names(sim.result$pars.normal)
    sum.stat.df <- data.frame(pars, sum.stat.df)
    
    # Remove masked values (if any)
    if (!is.null(mask)) sum.stat.df <- sum.stat.df[!mask, ]
    sum.stat.df
  }))
}
