# --------------------------------------------------------------
# Authors:  Paul R. Staab & Lisha Mathew
# Date:     2013-11-28
# Licence:  GPLv3 or later
# --------------------------------------------------------------

estimateLogLikelihood <- function(param, glm.fitted, sum.stats) {
  log.likelihood <- 0

  log.likelihood <- sum(sapply(seq(along=sum.stats), function(i) {
    name <- names(sum.stats)[i]
    if (sum.stats[[name]]$method %in% c("poisson.transformed", "poisson.independent")) {
      loglambda <- sapply(glm.fitted[[name]], predict, newdata=data.frame(t(as.matrix(param))))

      #if glm did not converge, take sum(SS[s]) or a small number like 0.5 
      loglambda[!sapply(glm.fitted[[name]], function(x) x$converged)] <- 0.5 

      sum.stat.value <- sum.stats[[name]]$value.transformed
      return(sum(sum.stat.value * loglambda - exp(loglambda) - calcLogFactorial(sum.stat.value))) 
    }
    else if(sum.stats[[name]]$method == "poisson.smoothing") {
      fake.sim.data <- list(pars.normal=param)
      fake.sim.data[[name]] <- sum.stats[[name]]$value
      new.data <- convertSimResultsToDataFrame(list(fake.sim.data), name)
      loglambda <- predict(glm.fitted[[name]][[1]], newdata=new.data)
      sum.stat.value <- new.data$sum.stat
      return(sum(sum.stat.value * loglambda - exp(loglambda) - calcLogFactorial(sum.stat.value)))
    }
    else stop("SumStat method not supported")
  }))

  return(log.likelihood)
}
