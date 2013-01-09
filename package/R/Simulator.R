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
simulateWithinBlock<- function(bObject, jaathaObject) {
  #time1<-Sys.time()
  ##values in [0-1]
  ##dim=c(lower /upper Boundry, #parameters)
  allBoundry <- array(sapply(1:bObject@nPar,
                             function(p) c(bObject@lowerBound[p],
                                           bObject@upperBound[p])),
                      dim=c(2,bObject@nPar))
  randompar <- aperm(array(runif(bObject@nPar*bObject@nSamp,
                                 min=allBoundry[1,],
                                 max=allBoundry[2,]),
                           dim=c(bObject@nPar,bObject@nSamp)))

  #Include corner points into simulation as well
  ##print(allBoundry)
  nCorners <- 2^bObject@nPar
  ##number of corners, corners will also be simulated
  for (c in 1:nCorners){
    ## converts 'c-1' to binary system,
    ##binary system bc corner is either at lower or upper Bound
    ##of parRange for each parameter
    digitalCorner <- .index2blocks(value=c-1, newBase=2,
                                   expo=bObject@nPar) + 1
    ## +1 bc R indices start at 1 (i.e. 1=lower and 2=upper bound)
    #cat("digital:",digitalCorner,"\n")
    corner <- sapply(1: bObject@nPar,
                     function(p) allBoundry[digitalCorner[p],p])
    #cat("   c",c," parameters:",corner,"\n")
    randompar <- rbind(randompar,corner,deparse.level = 0)
    ##deparse.level=0 makes no labels
  }

  # Create "packages" of parameters combinations for possible parallelization.
  sim.packages <- createSimulationPackages(randompar,
                                           jaathaObject@sim.package.size)

  seeds <- generateSeeds(length(sim.packages)+1)

  i <- NULL # To make R CMD check stop complaining
  # Simulate each package, maybe on different cores
  sumStats <- foreach(i = seq(along = sim.packages), .combine='rbind') %dopar% {
    set.seed(seeds[i])
    sim.pars <- .deNormalize(jaathaObject, 
                             sim.packages[[i]],
                             withoutTheta=jaathaObject@externalTheta)
    sumStats <- dm.simSumStats(jaathaObject@dm, sim.pars, jaathaObject@sum.stats.func)
    return(sumStats)
  }

  set.seed(seeds[length(seeds)])

  # Scale SumStats if we use scaling
  sumStats <- sumStats * jaathaObject@scaling.factor

  # Create combined output
  paraNsumstat <- cbind(randompar, sumStats) 
  .log2("Finished simulating for this block")
  return (paraNsumstat)
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
  if ( all(is.na(oldRange)) ) return(value) #external Theta
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
  nPar <- jObject@nPar+jObject@externalTheta-withoutTheta
  .log3("expecting",nPar,"parmeters")
  if (length(values) != nPar)
    stop("trying to deNormalize vector of wrong length")

  ret <- rep(0,nPar)
  ranges <- dm.getParRanges(jObject@dm)
  for (i in 1:nPar){
    ret[i] <- Jaatha.deNormalize01(ranges[i,],values[i])
  }

  #Add names
  names(ret) <- dm.getParameters(jObject@dm)[1:jObject@nPar]
  #.log(jObject,"Finished .deNormalizeVector. Result:",ret)
  return(ret)
}
