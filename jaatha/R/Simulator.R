## This function simulates nSamp random parameter combinations within
## lower and upper bounds for each parameter (boarders values
## needed to be between 0 and 1!!) and with the demographic
## model specified in simulate() with theta=5.  Returns the
## random parameters and the corresponding summary statistics. 
simulateWithinBlock<- function(bObject,jaathaObject){
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
  ##This are the random samples #aperm transposes matrices

  ## include corner points into simulation as well
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

  ## dim= c(#repetitions, parameters +sumstats)
  ## random values between 0 and 1 will be 'saved' in paraNsumstat

  #used fixed value for theta if we are estimating it externally
  #if (jaathaObject@externalTheta) randompar <- cbind(randompar,5)

  #simulate a line of sum stats for each line in randompars
  sumStats <- dm.simSumStats(jaathaObject@dm,
                             .deNormalize(jaathaObject,randompar,
                                          withoutTheta=jaathaObject@externalTheta))
  paraNsumstat <- cbind(randompar, sumStats) 

  return (paraNsumstat)
}



## Function that calls Hudson's ms (and/or seq-gen) with the given parameters. 
## Output
## is written to file called 'simOutput', if not specified.  
## Each parameter needs to be entered as par["paramIndex"] (each
## parameter has its number; needs to be the same order as in
## parameter ranges) To keep the order for the other functions: last
## parameter should be theta!!!
## If finiteSites version is used: 
## Each demo-model needs to be specified in 2 ways, an infinite sites way (i.e.
## as the ms command) and as a finite sites way (i.e. first call ms that 
## generates treeFile which will be read in by seq-gen.)  
#Jaatha.simulate <- function(dm,par,nLoci,fileName="simOutput",
#			    nSample=NA,finite=FALSE # unused
#			   )
#{
#       if (missing(par)) {
#		cat("No model parametes given. Using random values:\n")
#		par <- runif(dm.getNPar(dm))
#		print( .calcAbsParamValue(dm,par) )  
#       }
#	system( dm.simulationCmd(dm,relParamValue=par,nLoci=nLoci,fileName) )
#} 


##Function to map 'value' in oldRange to values between 0 and
##1. Returned will be a value between 0 and 1.
Jaatha.normalize01 <- function(oldRange, value){
  oldRange <- log(oldRange)
  value <- log(value)
  return ((value-min(oldRange))/(max(oldRange)-min(oldRange)))
}


##Function to map value between 0 and 1 to oldRange
## Returns single value (in oldRange).
Jaatha.deNormalize01 <- function(oldRange,value){
  if ( all(is.na(oldRange)) ) return(value) #external Theta
  oldRange <- log(oldRange)    
  return (exp(value*(max(oldRange)-min(oldRange))+min(oldRange)))
}

.deNormalize <- function(jObject, values, withoutTheta=F){
  if (!is.jaatha(jObject)) stop("jObject is no Jaatha object")
  if (!is.matrix(values)) values <- matrix(values,1)

  return( t(apply(values, 1, .deNormalizeVector,
                  jObject=jObject, withoutTheta=withoutTheta)) )
}

.deNormalizeVector <- function(jObject, values, withoutTheta){	
  #.log(jObject,"Called .deNormalizeVector")
  #.log(jObject,"values:",values,"| withoutTheta:",withoutTheta)
  if (!is.jaatha(jObject)) stop("jObject is no Jaatha object!")
  if (!is.numeric(values)) stop("trying to deNomalize non-numeric values")
  if (!is.vector(values)) stop("trying to deNormalize non-vector values")
  nPar <- jObject@nPar+jObject@externalTheta-withoutTheta
  #.log(jObject,"expecting",nPar,"parmeters")
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
