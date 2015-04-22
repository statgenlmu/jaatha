#' Estimates the likelihood for a parameter combination
#'
#' @param param The parameter combination. Must be a named vector where the
#'              names are the name of the parameters.
#' @param glm_fitted A list of the fitted GLMs, as produced by fitGLM().
#' @param sum_stats The summary statistic description, as in jaatha's sum.stats
#' @param scaling_factor The scaling factor used for the simulations
#' slot.
#' @return The estimated log-likelihood
estimateLogLikelihood <- function(param, glm_fitted, sum_stats, 
                                  scaling_factor = 1) {
  
  sum(sapply(seq(along=sum_stats), function(i) {
    name <- names(sum_stats)[i]
    calcStatLLH(sum_stats[[name]], glm_fitted[[name]], param, scaling_factor)
  }))
}

calcStatLLH <- function(sum_stat, ...) UseMethod("calcStatLLH")
calcStatLLH.default <- function(sum_stat, ...) stop('Unkown Summary Statistic')

calcStatLLH.Stat_PoiInd <- function(sum_stat, glm_fitted, param, 
                                    scaling_factor = 1) {
  
  loglambda <- sapply(glm_fitted, predict, 
                      newdata=data.frame(t(as.matrix(param))))
  
  #if glm did not converge, take sum(SS[s]) or a small number like 0.5 
  loglambda[!sapply(glm_fitted, function(x) x$converged)] <- 0.5 
  
  # Upscale predicted expectation value if we use scaling
  if (scaling_factor != 1) loglambda <- loglambda + log(scaling_factor)
  
  # Calculate the Poission log-likelihood
  sum(sum_stat$get_data() * loglambda - 
      exp(loglambda) - calcLogFactorial(sum_stat$get_data())) 
}

calcStatLLH.Stat_PoiSmooth <- function(sum_stat, glm_fitted, param, 
                                       scaling_factor = 1) {
  
  # Create a with parameters for prediction
  pars <- matrix(param, 1, length(param), byrow=TRUE)
  colnames(pars) <- names(param)
  data <- data.frame(pars,  sum_stat$get_data())
  
  suppressWarnings(loglambda <- predict(glm_fitted[['smooth']], newdata=data))
  
  sum(data$sum.stat * loglambda - 
        exp(loglambda) - calcLogFactorial(data$sum.stat))
#   if (!is.null(glm.fitted[[name]][['border']])) {
#     loglambda <- sapply(glm.fitted[[name]]$border, 
#                         predict, newdata=data.frame(t(as.matrix(param))))
#     loglambda[!sapply(glm.fitted[[name]]$border, function(x) x$converged)] <- 0.5 
#     
#     # Upscale predicted expectation value if we use scaling
#     if (scaling_factor != 1) loglambda <- loglambda + log(scaling_factor)
#     
#     sum.stat.value <- sum.stats[[name]]$border.transformed
#     log.li <- log.li + sum(sum.stat.value * loglambda - 
#                              exp(loglambda) - calcLogFactorial(sum.stat.value))
#   }
#   return(log.li)
}
