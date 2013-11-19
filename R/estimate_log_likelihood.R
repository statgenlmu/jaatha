# --------------------------------------------------------------
# Authors:  Paul R. Staab & Lisha Mathew
# Date:     2013-11-14
# Licence:  GPLv3 or later
# --------------------------------------------------------------

estimateLogLikelihood <- function(param, glm.fitted, sum.stats) {
  if(min(param)<0 || max(param)>1){ stop("Optimization outside par space") }
  log.likelihood <- 0

  log.likelihood <- sum(sapply(seq(along=sum.stats), function(i) {
    name <- names(sum.stats)[i]
    if (sum.stats[[name]]$method %in% c("poisson.transformed", "poisson.independent")) {
      loglambda <- sapply(glm.fitted[[name]], predict, newdata=data.frame(t(as.matrix(param))))

      #if glm did not converge, take sum(SS[s]) or a small number like 0.5 
      loglambda[!sapply(glm.fitted[[name]], function(x) x$converged)] <- 0.5 

      sum.stat.value <- sum.stats[[name]]$value.transformed
      sum(sum.stat.value * loglambda - exp(loglambda) - calcLogFactorial(sum.stat.value)) 
    }
  }))

  return(log.likelihood)
}
