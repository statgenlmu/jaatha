# --------------------------------------------------------------
# Block.R
# Class for repesenting a rectangular area of the parameter space
# 
# Authors:  Lisha Naduvilezhath & Paul R. Staab
# Date:     2012-10-05
# Licence:  GPLv3 or later
# --------------------------------------------------------------


slots <- representation(nPar="numeric",         # number of parameters = number of block dimensions
                        nSamp="numeric",        # how many samples to take
                        nLoci="numeric",        #number of loci to use for each sampling point
                        lowerBound="numeric",   # vector of lower bounds [in range: 0-1] 
                        upperBound="numeric",   # vector of upper bounds [in range: 0-1] 
                        score="numeric",        # log of the maximum (composite)

                        # likelihood within block 
                        MLest="numeric",        # ML estimations of parameters [in range: 0-1] 
                        parNsumstat ="array",   #contains the random parameters in this block and corresponding summary statistics 
                        weight="numeric")       #each simulation within this block will be weighted with this value in glmFitting


setClass("Block" , representation=slots)


## Shows the content of the slots of the Block object.
showBlock <- function(block) {
  cat("*** Object of class BLOCK ***\n")           
  cat(" nPar =",block@nPar,"\n")            
  cat(" lowerBound =",block@lowerBound,"\n")           
  cat(" upperBound =",block@upperBound,"\n")
  cat(" nSamp =",block@nSamp,"\n")             
  cat(" nLoci =",block@nLoci,"\n")             
  cat(" weight =",block@weight,"\n")                      
  cat(" dimensions of parNsumstat =",
      dim(block@parNsumstat)," ")
  if (dim(block@parNsumstat)[1]!=0) {
    cat(block@parNsumstat[1,1],block@parNsumstat[2,1]," ...\n")
  } else {
    cat("\n")
  }
  if (length(block@MLest)!=0){
    cat(" MLest [0..1]=", round(block@MLest,6),
        "with score",block@score,"\n")
  } else{}
  cat("*** End of Object of class BLOCK ***\n")   
}

setMethod("show", "Block", showBlock)
rm(showBlock)


##Function to determine if a given point is within the parameter range
##of the given block. Bounds and point are rounded to the 5th decimal
##place.
isInBlock <- function(object,point) {
  if(length(point)!=object@nPar) {
    print(list(ERROR="Parameter dimensions unequal to block dimensions! \n")) 
  } else{}
  lower <- round(object@lowerBound,5)
  upper <- round(object@upperBound,5)
  point <- round(point,5)
  return(all(c(lower<=point,point<=upper)))
}


## Function to determine the (starting) blockIndex of a given point in
## the parameter range specified in jObject.
## Returns the index of the starting block which contains the point.
###########################################################################
# Not used anywhere => commented
# Paul 29.05.12
###########################################################################
#.getBlockNumber <- function(jObject, point){
#	## pRange contains the boarders of all starting blocks
#	## dim=c(#parameter,#blocks per dimension, start&end)
#	pRange <- .calcBlockParRanges01(jObject@nPar,jObject@nBlocksPerPar)  
#	nTotalBlocks <- (jObject@nBlocksPerPar^jObject@nPar)
#	#if(jObject@externalTheta){nPar <- jObject@nPar+1} else {nPar <- jObject@nPar}
#	
#	## check if dimensions of point and jObject are same
#	if (length(point)!=jObject@nPar){
#		return(print(list(ERROR=paste("Dimension of point and parRange are not 
#					       the same! dim(point): ",
#					       length(point),"!=",jObject@nPar," \n"))))
#	}
#	
#	## convert point into [0-1]-range
#	point01 <- sapply(1:jObject@nPar, function(d) normalize01(
#						jObject@parRange[[d]],point[d]))
#	i <- 1
#	## go through each block until you find the containing block
#	while (i <= nTotalBlocks){
#		## 'b' determines which block-index for each parameter
#		## that is being considered
#		b <- .index2blocks(value=i-1, newBase=jObject@nBlocksPerPar,
#				expo=jObject@nPar) + 1  ##+1 bc R indices start with 0
#		boundry <- sapply(1:jObject@nPar, function(p) pRange[p,b[p],])  #dim=c(2,jObject@nPar)
#		block <- new("Block", nPar=jObject@nPar,
#				lowerBound= boundry[1,], nLoci=1,
#				upperBound= boundry[2,],
#				nSamp=1)
#		if (isInBlock(object=block,point=as.vector(point01))) {
#			return (i)
#		} else{
#			i <- i+1
#		}
#	}
#	cat("Point is not in the given parameter range!\n")
#}
