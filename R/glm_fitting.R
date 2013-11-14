# --------------------------------------------------------------
# glm_fitting.R
# Functions for the machine learning part of Jaatha. 
# 
# Authors:  Paul R. Staab & Lisha Mathew
# Date:     2013-11-13
# Licence:  GPLv3 or later
# --------------------------------------------------------------

fitGlm <- function(sim.data, jaatha, weighting=NULL){ 
  for (i in seq(along = jaatha@sum.stats)) {
    name <- names(jaatha@sum.stats)[i] 
    if (jaatha@sum.stats[[i]]$method %in% c("poisson.transformed",
                                            "poisson.independent")) {
      glm.fitted <- fitGlmTransformed(sim.data, name, 
                                      jaatha@sum.stats[[i]]$transformation, 
                                      weighting, jaatha)
    }
    else stop("method not found")
  }
  return(glm.fitted)
}

fitGlmTransformed <- function(sim.data, sum.stat, transformation, weighting, jaatha) {
  stopifnot(!is.null(sum.stat))
  stopifnot(is(sim.data[[1]]) != "try-error")

  stats.sim <- t(sapply(sim.data, 
                        function(x) c(x$pars.normal, transformation(x[[sum.stat]])))) 
  stats.names <- paste("S", 1:(ncol(stats.sim)-length(getParNames(jaatha))), sep="")
  colnames(stats.sim) <- c(getParNames(jaatha), stats.names) 

  formulas <- paste0(stats.names, "~", paste(getParNames(jaatha) ,collapse= "+"))
  lapply(formulas, glm, data=data.frame(stats.sim), family=poisson,
         control = list(maxit = 200))
}


## Function to estimate the best parameters within the block and for
## those the highest composite likelihood for the specific model
## (specified in Simulator.simulateWithinBlock()). Input:
## Block-'object' to be optimized within, 'jObject' for getting the
## parameter ranges for the search, 'modFeld' which holds the
## coefficients, convergence and sum of the fitted glm for each
## summary statistic, 'ssData' the summary statistic of the data for
## which parameters to search. 'boarder' determines how much of the
## around the boundaries of the blocks should not be used for the
## optimization procedure. 
estimateMlInBlock <- function(block, glm.fitted, sum.stats) {
  block.size <- block@border[,2] - block@border[,1]
  block.middle <- block.size/2 + block@border[,1]

  ##describes 'boarder'% of values that will be excluded
  ##on either side of the block in optimization
  best.value <- optim(block.middle, calcLogLikelihood, 
                      glm.fitted=glm.fitted, sum.stats=sum.stats,  
                      lower=block@border[ ,1], upper=block@border[ ,2],
                      method="L-BFGS-B", control=list(fnscale=-1))

  return(list(est=best.value$par, score=best.value$value))                   
}


calcLogLikelihood <- function(param, glm.fitted, sum.stats) {
  if(min(param)<0 || max(param)>1){ stop("Optimization outside par space") }
  #print(param)
  log.likelihood <- 0

  log.likelihood <- sum(sapply(sum.stats, function(sum.stat) {
    if (sum.stat$method %in% c("poisson.transformed", "poisson.independent")) {
      loglambda <- sapply(glm.fitted, predict, newdata=data.frame(t(as.matrix(param))))

      #if glm did not converge, take sum(SS[s]) or a small number like 0.5 
      loglambda[!sapply(glm.fitted, function(x) x$converged)] <- 0.5 

      sum.stat.value <- sum.stat$transformation(sum.stat$value)
      sum(sum.stat.value * loglambda - exp(loglambda) - calcLogFactorial(sum.stat.value)) 
    }
  }))

  #print(log.likelihood)
  return(log.likelihood)
}
