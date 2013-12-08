# --------------------------------------------------------------
# Block.R
# Class for repesenting a rectangular area of the parameter space
# 
# Authors:  Lisha Naduvilezhath & Paul R. Staab
# Date:     2013-09-04
# Licence:  GPLv3 or later
# --------------------------------------------------------------


slots <- representation(border="matrix", #[in range: 0-1] 
                        score="numeric", # log of the maximum (composite)

                        # likelihood within block 
                        MLest="numeric",        # ML estimations of parameters [in range: 0-1] 
                        parNsumstat ="data.frame",   #contains the random parameters in this block and corresponding summary statistics 
                        weight="numeric")       #each simulation within this block will be weighted with this value in glmFitting


setClass("Block" , representation=slots)


## Shows the content of the slots of the Block object.
showBlock <- function(object) {
  cat("*** Object of class BLOCK ***\n")           
  print(object@border)
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


isInBlock <- function(block, point) {
  length(point) != nrow(block@border) && stop("Point and block dimensions
                                              mismatch")
  res <- all(block@border[,1]-1e-15<=point & point<=block@border[,2]+1e-15)
  return(res)
}

getCorners <- function(block) {
  corners <- foreach(c=1:2^nrow(block@border), .combine=rbind) %do% {
    ## converts 'c-1' to binary system,
    ##binary system bc corner is either at lower or upper Bound
    ##of parRange for each parameter
    digitalCorner <- .index2blocks(value=c-1, newBase=2,
                                   expo=nrow(block@border)) + 1
    ## +1 bc R indices start at 1 (i.e. 1=lower and 2=upper bound)
    corner <- sapply(1:nrow(block@border), function(p) block@border[p, digitalCorner[p]])
    return(corner)
  }
  rownames(corners) <- NULL
  colnames(corners) <- rownames(block@border)
  return(corners)
}

printBorder <- function(block, jaatha) {
  lower <- denormalize(block@border[ ,1], jaatha)
  upper <- denormalize(block@border[ ,2], jaatha)
  return(paste0(round(lower,3), "-", round(upper,3),
                collapse=" x "))
}
