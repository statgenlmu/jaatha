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
  lapply(formulas, glm, data=data.frame(stats.sim), family=poisson,
         control = list(maxit = 200))
}

#' Fits a GLM for a summary statistics of type "poisson.smoothed"
#'
#' @param sim.data Results from simulations
#' @param sum.stat Name of the summary statistics
#' @param weighting Potentially weights for the simulations.
#' @param jaatha A Jaatha Object.
#' @return A list with one fitted GLMs
fitPoiSmoothed <- function(sim.data, sum.stat, weighting, jaatha) {
  model <- paste0("sum.stat ~ ",
                  "(", jaatha@sum.stats[[sum.stat]]$model, ")",  
                  "*(", paste(getParNames(jaatha), collapse="+"), ")") 

  sim.data.df <- convertSimResultsToDataFrame(sim.data, sum.stat)
  list(glm(model, data=sim.data.df, family=poisson("log")))
}


#' Converts simulation results into a data frame that is usable for fitting a
#' glm.
#'
#' Currently only works with nx2 matix summary statistics and vectorizes thoose.
#' 
#' @param sim.data Results from simulations
#' @param sum.stat Name of the summary statistics which should get converted
#' @return The summary statistics as data.frame 
convertSimResultsToDataFrame <- function(sim.data, sum.stat) {
  do.call(rbind, lapply(sim.data, function(sim.result) {
    i <- as.vector(row(sim.result[[sum.stat]]))
    j <- as.vector(col(sim.result[[sum.stat]]))
    value <- as.vector(sim.result[[sum.stat]])
    pars <- matrix(sim.result$pars.normal, length(i),
                   length(sim.result$pars.normal), byrow=TRUE)
    colnames(pars) <- names(sim.result$pars.normal)
    data.frame(pars, i, j, sum.stat=value)
  }))
}
