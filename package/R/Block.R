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
showBlock <- function(object) {
  cat("*** Object of class BLOCK ***\n")           
  cat(" nPar =",object@nPar,"\n")            
  cat(" lowerBound =",object@lowerBound,"\n")           
  cat(" upperBound =",object@upperBound,"\n")
  cat(" nSamp =",object@nSamp,"\n")             
  cat(" nLoci =",object@nLoci,"\n")             
  cat(" weight =",object@weight,"\n")                      
  cat(" dimensions of parNsumstat =",
      dim(object@parNsumstat)," ")
  if (dim(object@parNsumstat)[1]!=0) {
    cat(object@parNsumstat[1,1],object@parNsumstat[2,1]," ...\n")
  } else {
    cat("\n")
  }
  if (length(object@MLest)!=0){
    cat(" MLest [0..1]=", round(object@MLest,6),
        "with score",object@score,"\n")
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
