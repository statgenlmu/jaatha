# --------------------------------------------------------------
# glm_fitting.R
# Functions for the machine learning part of Jaatha. 
# 
# Authors:  Paul R. Staab & Lisha Mathew
# Date:     2013-09-04
# Licence:  GPLv3 or later
# --------------------------------------------------------------

## Function to fit a glm for each summary statistic with the generated
## parameter combinations in the block. The first 'bObject@nPar'
## columns are assumed to be the columns for the parameters, the
## following 'nTotalSumstat' colums contain the results for the summary
## statistic.
glmFitting <- function(sim.result, jaatha, weighting=NULL){ 
  #cat("Fitting model ... \n")
  ##+3 bc intercept, convergence, sumOfSumstat
  nTotalSumstat <- ncol(sim.result)-jaatha@nPar
  modFeld <- array(-1,dim=c(nTotalSumstat,(jaatha@nPar+3)))
  colnames(modFeld) <- c("Intercept", jaatha@par.names, "conv", "ss.sum")

  explanatory <- paste(jaatha@par.names ,collapse= "+")

  ## suppresses warnings; which ocurr whenever sumstat is 0 "suppressWarnings()"		   
  for(s in 1:nTotalSumstat) {
    if ( sum(sim.result[,s+jaatha@nPar]) == 0 ) {
      modFeld[s, ] <- 0 
      next()
    }

    ## if glm function did not converge, set coefficients & convergence to 0
    tryCatch({
      mod <- suppressWarnings(glm(as.formula(paste0("SS", s," ~ ", explanatory)),
                                  data=sim.result, family=poisson, weights=weighting,
                                  control = list(maxit = 200)))}, 
             error=function(e) {
               print(list("Caught error of GLM function!"))
               print(e)
               cat("sumstat",s," sum=",sum(sim.result[,s+jaatha@nPar]), ": ", "\n")
               mod <- list(coef=rep(0,jaatha@nPar+1),conv=0) })

    modFeld[s,] <- c(mod$coef,mod$conv,sum(sim.result[,s+jaatha@nPar])) 
    if (!mod$conv){
      cat("WARNING: sumstat",s," did not converge, sum = ",
          sum(sim.result[,s+jaatha@nPar]),"\n")
    }
  }  
  return (modFeld)
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
estimate <- function(bObject, jaatha, modFeld, boarder=0.25){
  dimSize <- bObject@border[,2] - bObject@border[,1]
  mitte <- dimSize/2 + bObject@border[,1]

  # function for optimization  
  optfunk <- function(par) {                                 
    if(min(par)<0 || max(par)>1){ return(1e11)} else{}
    score <- .calcScore(param=par, modFeld, jaatha)
    return(-score)  ## '-' enables Maximization
  }

  ##boarder to leave free in the block for optimization
  ##(depends on boundaries)
  puffer <- boarder*dimSize   
  ##calculate limits for the optimization procedure
  untere <- sapply(1:jaatha@nPar,
                   function(p) min(bObject@border[p,1] + puffer,1))
  obere <-  sapply(1:jaatha@nPar,
                   function(p) max(bObject@border[p,2] - puffer,0))

  ##describes 'boarder'% of values that will be excluded
  ##on either side of the block in optimization
  OOO <- optim(mitte, optfunk, lower=untere, upper=obere,
               method="L-BFGS-B")
  mitte <- OOO$par    ##the optimal parameters are contained in the vector mitte
  score <- -OOO$value  

  return(list(est=mitte, score=score))                   
}


.calcScore <- function(param, model.coef, jaatha){
  loglambda <- model.coef[ ,1:(jaatha@nPar+1)] %*% c(1, param) 
  #if glm did not converge, take sum(SS[s]) or a small number like 0.5 
  loglambda[model.coef[ , 'conv'] < 0.5] <- 0.5 
  return(sum( jaatha@sumStats * loglambda - exp(loglambda) )) 
}
