# --------------------------------------------------------------
# Authors:  Paul R. Staab & Lisha Mathew
# Date:     2013-11-28
# Licence:  GPLv3 or later
# --------------------------------------------------------------

#' Estimates the likelihood for a parameter combination
#'
#' @param param The parameter combination. Must be a named vector where the
#'              names are the name of the parameters.
#' @param glm.fitted A list of the fitted GLMs, as produced by fitGLM().
#' @param sum.stats The summary statistic description, as in jaatha's sum.stats
#' @param scaling_factor The scaling factor used for the simulations
#' slot.
#' @return The estimated log-likelihood
estimateLogLikelihood <- function(param, glm.fitted, sum.stats, 
                                  scaling_factor = 1) {
  log.likelihood <- 0

  log.likelihood <- sum(sapply(seq(along=sum.stats), function(i) {
    name <- names(sum.stats)[i]
    if (sum.stats[[name]]$method %in% c("poisson.transformed", "poisson.independent")) {
      loglambda <- sapply(glm.fitted[[name]], predict, newdata=data.frame(t(as.matrix(param))))

      #if glm did not converge, take sum(SS[s]) or a small number like 0.5 
      loglambda[!sapply(glm.fitted[[name]], function(x) x$converged)] <- 0.5 

      # Upscale predicted expectation value if we use scaling
      if (scaling_factor != 1) loglambda <- loglambda + log(scaling_factor)
      
      sum.stat.value <- sum.stats[[name]]$value.transformed
      return(sum(sum.stat.value * loglambda - exp(loglambda) - calcLogFactorial(sum.stat.value))) 
    }
    else if(sum.stats[[name]]$method == "poisson.smoothing") {
      fake.sim.data <- list(pars.normal=param)
      fake.sim.data[[name]] <- sum.stats[[name]]$value
      new.data <- convertSimResultsToDataFrame(list(fake.sim.data), name, 
                                               sum.stats[[name]]$border.mask)
      suppressWarnings(loglambda <- predict(glm.fitted[[name]][['smooth']], newdata=new.data))

      sum.stat.value <- new.data$sum.stat
      log.li <- sum(sum.stat.value * loglambda - exp(loglambda) - calcLogFactorial(sum.stat.value))
      if (!is.null(glm.fitted[[name]][['border']])) {
        loglambda <- sapply(glm.fitted[[name]]$border, 
                            predict, newdata=data.frame(t(as.matrix(param))))
        loglambda[!sapply(glm.fitted[[name]]$border, function(x) x$converged)] <- 0.5 

        # Upscale predicted expectation value if we use scaling
        if (scaling_factor != 1) loglambda <- loglambda + log(scaling_factor)
        
        sum.stat.value <- sum.stats[[name]]$border.transformed
        log.li <- log.li + sum(sum.stat.value * loglambda - 
                               exp(loglambda) - calcLogFactorial(sum.stat.value))
      }
      return(log.li)
    }
    else stop("SumStat method not supported")
  }))

  return(log.likelihood)
}
