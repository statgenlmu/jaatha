# --------------------------------------------------------------
# Simulator.R
# Functions related to the simulation of DNA sequence data.
# 
# Authors:  Lisha Naduvilezhath & Paul R. Staab
# Date:     2012-10-05
# Licence:  GPLv3 or later
# --------------------------------------------------------------

## This function simulates nSamp random parameter combinations within
## lower and upper bounds for each parameter (boarders values
## needed to be between 0 and 1!!) and with the demographic
## model specified in simulate() with theta=5.  Returns the
## random parameters and the corresponding summary statistics. 
simulateWithinBlock<- function(sim, block, jaatha) {
  # Sample random simulation parameters
  randompar <- aperm(array(runif(jaatha@nPar*sim,
                                 min=block@border[,1],
                                 max=block@border[,2]),
                           dim=c(jaatha@nPar,sim)))
  
  # Add the corners of the block to sim parameters
  randompar <- rbind(randompar, getCorners(block))

  # Create "packages" of parameters combinations for possible parallelization.
  sim.packages <- createSimulationPackages(randompar, jaatha@sim.package.size)
  seeds <- generateSeeds(length(sim.packages)+1)

  i <- NULL # To make R CMD check stop complaining
  # Simulate each package, maybe on different cores
  sumStats <- foreach(i = seq(along = sim.packages), .combine='rbind') %dopar% {
    set.seed(seeds[i])
    sim.pars <- .deNormalize(jaatha, 
                             sim.packages[[i]])
    sumStats <- jaatha@simFunc(jaatha, sim.pars)
    return(sumStats)
  }

  set.seed(seeds[length(seeds)])

  # Create combined output
  sim.result <- data.frame(cbind(randompar, sumStats))
  colnames(sim.result) <- c(jaatha@par.names, paste("SS", 1:ncol(sumStats),
                                                    sep=""))  
  .log2("Finished simulating for this block")
  return(sim.result)
}

createSimulationPackages <- function(random.par, package.size) {
  if (package.size == 0) return(list(random.par))

  sim.packages <- list()
  num.pars <- nrow(random.par)
  i <- 0
  
  while (i < num.pars/package.size) {
    i <- i + 1
    lower <- (i-1)*package.size+1
    upper <-  min(i*package.size, num.pars)
    sim.packages[[i]] <- random.par[lower:upper, , drop=F]
  }

  return(sim.packages)  
}

##Function to map 'value' in oldRange to values between 0 and
##1. Returned will be a value between 0 and 1.
Jaatha.normalize01 <- function(oldRange, value){
  oldRange <- log(oldRange)
  value <- log(value)
  return ((value-min(oldRange))/(max(oldRange)-min(oldRange)))
}


##Function to map value between 0 and 1 to oldRange
## Returns single value (in oldRange).
Jaatha.deNormalize01 <- function(oldRange, value){
  oldRange <- log(oldRange)    
  return (exp(value*(max(oldRange)-min(oldRange))+min(oldRange)))
}

.deNormalize <- function(jObject, values, withoutTheta=F){
  if (!is.jaatha(jObject)) stop("jObject is no Jaatha object")
  if (!is.matrix(values)) stop("values is no matrix!") 
  
  result <- apply(values, 1, .deNormalizeVector,
                  jObject=jObject, withoutTheta=withoutTheta)
  if (!is.matrix(result)) result <- matrix(result,1)
  result <- t(result)
  return(result)
}

.deNormalizeVector <- function(jObject, values, withoutTheta){	
  #.log(jObject,"Called .deNormalizeVector")
  .log3("Denormalizing parameters...")
  .log3("values:",values,"| withoutTheta:",withoutTheta)
  if (!is.jaatha(jObject)) stop("jObject is no Jaatha object!")
  if (!is.numeric(values)) stop("trying to deNomalize non-numeric values")
  if (!is.vector(values)) stop("trying to deNormalize non-vector values")
  nPar <- jObject@nPar
  .log3("expecting",nPar,"parmeters")
  if (length(values) != nPar)
    stop("trying to deNormalize vector of wrong length")

  ret <- rep(0,nPar)
  ranges <- jObject@par.ranges
  for (i in 1:nPar){
    ret[i] <- Jaatha.deNormalize01(ranges[i,],values[i])
  }

  #Add names
  names(ret) <- jObject@par.names 
  #.log(jObject,"Finished .deNormalizeVector. Result:",ret)
  return(ret)
}
