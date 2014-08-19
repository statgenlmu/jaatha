# --------------------------------------------------------------
# Functions to fit GLMs to learn how the summary statistics depend
# on the parameters.
# 
# Authors:  Paul R. Staab & Lisha Mathew
# Date:     2013-11-14
# Licence:  GPLv3 or later
# --------------------------------------------------------------

#' Fits a Generalized Linear Model within a block.  
#'
#' The model describes how the summary statistics depend on the parameters.
#' 
#' @param sim.data Simulation data in this block.
#' @param jaatha A Jaatha Object.
#' @param weighting Potentially weights for the simulations.
#' @return A list containing a list of fitted GLMs for each summary
#' statistic.
fitGlm <- function(sim.data, jaatha, weighting=NULL){ 
  glm.fitted <- list()
  for (i in seq(along = jaatha@sum.stats)) {
    name <- names(jaatha@sum.stats)[i] 
    if (jaatha@sum.stats[[i]]$method %in% c("poisson.transformed",
                                            "poisson.independent")) {
      glm.fitted[[name]] <- fitGlmPoiTransformed(sim.data, name, 
                                      jaatha@sum.stats[[i]]$transformation, 
                                      weighting, jaatha)

      stopifnot(length(glm.fitted[[name]]) == 
                length(jaatha@sum.stats[[name]]$value.transformed))
    }
    else if (jaatha@sum.stats[[i]]$method == "poisson.smoothing") {
      glm.fitted[[name]] <- fitPoiSmoothed(sim.data, name, weighting, jaatha)
    } 
    else stop("method not found")
  }
  return(glm.fitted)
}

#' Fits a GLM for a summary statistics of type "poisson.transformed" and
#' "poisson.independent"
#'
#' @param sim.data Results from simulations
#' @param sum.stat Name of the summary statistics
#' @param transformation Transformation function to apply to the data
#' @param weighting Potentially weights for the simulations.
#' @param jaatha A Jaatha Object.
#' @return A list of fitted GLMs, one for each function
fitGlmPoiTransformed <- function(sim.data, sum.stat, transformation, weighting, jaatha) {
  stopifnot(!is.null(sum.stat))

  stats.sim <- t(sapply(sim.data, 
                        function(x) c(x$pars.normal, transformation(x[[sum.stat]])))) 
  stats.names <- paste("S", 1:(ncol(stats.sim)-length(getParNames(jaatha))), sep="")
  colnames(stats.sim) <- c(getParNames(jaatha), stats.names) 

  formulas <- paste0(stats.names, "~", paste(getParNames(jaatha) ,collapse= "+"))
  glms <- lapply(formulas, glm, data=data.frame(stats.sim), family=poisson,
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
