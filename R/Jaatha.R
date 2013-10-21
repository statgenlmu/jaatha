# --------------------------------------------------------------
# Jaatha.R
# This file contains the Jaatha S4-Class and a few related 
# helper functions 
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Date:     2013-09-04
# Licence:  GPLv3 or later
# --------------------------------------------------------------


#' Fast estimation of demographic parameters
#' 
#' Jaatha is a composite maximum likelihood method to estimate parameters
#' of a speciation model of closely related biological species out of SNP
#' data.
#' 
#' @author
#' Lisha Naduvilezhath \email{lisha (at) biologie.uni-muenchen.de},
#' Paul R. Staab \email{staab (at) biologie.uni-muenchen.de} and 
#' Dirk Metzler \email{metzler (at) biologie.uni-muenchen.de}
#'
#' Maintainer: Paul R. Staab \email{staab (at) biologie.uni-muenchen.de}
#' @name jaatha-package
#' @docType package
#' @title Fast estimation of demographic parameters
#' @keywords package
#' @importFrom phyclust ms
#' @importFrom foreach foreach 
#' @importFrom foreach registerDoSEQ
#' @importFrom foreach %do% 
#' @importFrom foreach %dopar% 
#' @importFrom methods new 
#' @importFrom methods representation 
#' @import Rcpp
#' @useDynLib jaatha
NULL

#' The "Jaatha" S4 class saves the basic parameters for a Jaatha estimation
#' procedure
#'
#' Slots:
#' \describe{
#'    \item{dm}{The demographic model to use}
#'    \item{nPar}{The number of parameters to estimate}
#'    \item{logMLmax}{The maximum log composite likelihood over all block}
#'    \item{nTotalSumstat}{The total number of summary statistics}
#'    \item{nNonJsfsSumstat}{The number of summary statistics that are not calculated from the JSFS}
#'    \item{popSampleSizes=}{Sample sizes of both populations. Vector of length 2.}
#'    \item{nLoci}{The number of loci}
#'    \item{MLest}{ML estimations of parameters (transformed between 0 and 1)}
#'    \item{seeds}{A set of random seeds. First one to generate the other two.
#'                 The next one is for the initial search, the last is for the
#'                 refined search}
#'    \item{externalTheta}{Option removed in Version 2.1. Use a Version 2.0.2 if
#'    you want to use external Theta.}
#'    \item{finiteSites}{If TRUE, we use a finite sites mutation model instead of an infinite sites one}
#'    \item{parNames}{The name of the parameters we estimate}
#'    \item{debugMode}{If TRUE, a debug output will be produced}
#'    \item{logFile}{If set, the debug output will be written into this file}
#'    \item{sumStats}{The summary statistics from the real data}
#'    \item{likelihood.table}{A matrix with the best composite log likelihood values and 
#'                            corresponding parameters}
#'    \item{sum.stats.func}{The function for summarizing the JSFS}
#'    \item{sim.package.size}{Execute this number of simulations on a core in a
#'                            row}
#'    \item{cores}{The number of CPU cores to use for simulations}
#'    \item{scaling.factor}{Only simulate this part of the data and interpolate
#'                          the rest}
#'    \item{route}{Tracks the best estimates of each step.}
#'    \item{use.shm}{Use the shared memory /dev/shm for temporary files (Linux only)}
#' }
#'
#' @name Jaatha-class
#' @rdname Jaatha-class
#' @exportClass Jaatha
setClass("Jaatha",
  representation=representation(
      # Settings
      simFunc="function",
      nPar="numeric",
      par.names = "character",
      par.ranges = "matrix",
      seeds="numeric",
      sumStats = "numeric",
      sim.package.size = "numeric",
      cores = "numeric",
      scaling.factor = "numeric",
      use.shm = "logical",
      opts = "list",

      # Results
      logMLmax = "numeric",
      MLest="numeric", 
      starting.positions = "list",
      likelihood.table = "matrix",
      route = "list"
    ),
)

## constructor method for Jaatha object
.init <- function(.Object, sim.func, par.ranges, 
                  sum.stats, seed, cores, 
                  sim.package.size, use.shm=FALSE) {

  .log3("Starting initialization")
  
  .Object@simFunc <- sim.func
  .Object@sumStats <- sum.stats
  .Object@use.shm <- use.shm
  .Object@opts <- list()

  checkType(par.ranges, c("matrix"))
  dim(par.ranges)[2] == 2 || stop("par.ranges must have two columns")
  .Object@par.ranges <- par.ranges
  .Object@nPar <- dim(par.ranges)[1]
  
  if (is.null(rownames(par.ranges))) 
    rownames(par.ranges) <- as.character(1:.Object@nPar) 
  .Object@par.names <- rownames(par.ranges)
  
  # Parallelization options
  if (missing(sim.package.size)) sim.package.size <- 10
  checkType(sim.package.size, c("num","single"))
  .Object@sim.package.size <- sim.package.size

  if (missing(cores)) cores <- 1
  checkType(cores, c("num","single"))
  .Object@cores <- cores

  .Object@likelihood.table <- matrix()
  .Object@starting.positions <- list()

  # Seeds
  # Jaatha uses three seeds. The first is the "main seed" used to generate the
  # other two seeds if provided, the second is the seed for the initial search
  # and the refined search.
  if (missing(seed) || length(seed) == 0) seed <- generateSeeds(3)
  else if (length(seed) == 1) {
    checkType(seed, "num")
    set.seed(seed)
    seed[2:3] <- generateSeeds(2)
  }
  else stop("Malformated argument: seed")
  .Object@seeds <- seed

  .log3("Finished initialization")

  return (.Object)
}

setMethod(f="initialize", signature ="Jaatha", definition=.init)
rm(.init)


#' Initialization of a Jaatha estimation for population genetics
#'
#' This function sets the basic parameters for an analysis with
#' Jaatha and is the first step for each application of it.
#'
#' @param demographic.model The demographic model to use
#' @param jsfs Your observed Joint Site Frequency Spectrum (JSFS). Jaatha uses
#'        the JSFS as summary statistics.   
#' @param folded If 'TRUE', Jaatha will assume that the JSFS is folded.
#' @param seed An integer used as seed for both Jaatha and the simulation software
#' @param sim.package.size When running Jaatha on multiple cores, a singe core
#'              will always execute a whole "package" of simulations the reduce
#'              the inter thread communication overhead. This gives the number
#'              of simulations in such a package. Choose a number such that the
#'              execution of the simulations takes at least 15 seconds.
#' @param cores The number of cores to use in parallel. If 0, it tries to
#'              guess the number of available cores and use them all.
#' @param scaling.factor You can use this option if you have a large dataset. If
#'              so, Jaatha only simulates only a fraction 1/scaling.factor of the
#'              dataset and interpolates the missing data.
#' @param use.shm Logical. Many modern linux distributions have a shared memory
#'              file system available under /dev/shm. Set this to TRUE to use it for
#'              temporary files. Usually gives a huge performance boost.
#' @return A S4-Object of type jaatha containing the settings
#' @examples
#' dm <- dm.createThetaTauModel(c(20,25), 100) 
#' jsfs <- matrix(rpois(21*26, 5), 21, 26)
#' jaatha <- Jaatha.initialize(dm, jsfs) 
#' 
#' @export
Jaatha.initialize <- function(demographic.model, jsfs,
                              seed, sim.package.size=10,
                              cores=1, scaling.factor=1,
                              use.shm=FALSE, folded=FALSE) {

  if (is.list(jsfs)) jsfs <- jsfs[[1]]$jsfs

  checkType(demographic.model, "dm")
  checkType(jsfs, "num")
  checkType(folded, c("bool", "single"))

  if (!folded) {
    sum.stats <- summarizeJSFS(jsfs)
    sim.func <- simulateDemographicModel
  }
  else {
    sum.stats <- summarizeFoldedJSFS(jsfs)
    sim.func <- simulateDemographicModelFolded
  }

  if (missing(seed)) seed <- numeric()

  jaatha <- new("Jaatha", sim.func=sim.func, 
                par.ranges=as.matrix(dm.getParRanges(demographic.model)),  
                sum.stats=sum.stats,
                seed=seed,
                sim.package.size=sim.package.size,
                cores=cores,
                use.shm=use.shm)

  checkType(scaling.factor, c("num","single"))
  demographic.model <- scaleDemographicModel(demographic.model, scaling.factor)
  jaatha@opts[['scaling.factor']] <- scaling.factor

  demographic.model <- finalizeDM(demographic.model)
  jaatha@opts[['dm']] <- demographic.model
  return(jaatha)
}


## Shows the content of the slots of the Jaatha object.
.show <- function(object) {
  initial.done <- length(object@starting.positions) > 0
  refined.done <- !is.na(object@likelihood.table[1,1])
  cat("This is a container object for everything related to a Jaatha analysis.\n\n")
  cat("Status of this analysis:\n")
  cat("Initialization... done\n")
  cat("Initial Search... ")
  if (initial.done) cat("done\n")
  else              cat("\n")
  cat("Refined Search... ")
  if (refined.done) cat("done\n")
  else              cat("\n")
  cat("\n")
  
  if(initial.done & !refined.done) {
    cat("Possible starting positios:\n")
    print(Jaatha.getStartingPoints(object))
    cat("\n")
  }
  if(refined.done) {
    cat("Maximum composite likelihood estimates:\n")
    print(Jaatha.getLikelihoods(object))
    cat("\n")
  }

  if(!initial.done & !refined.done) {
    cat("Next Step: Run Jaatha.initialSearch()\n")
  }
  if(initial.done & !refined.done) {
    cat("Next Step: Run Jaatha.refinedSearch()\n")
  }
}

setMethod("show","Jaatha",.show)
rm(.show)







Jaatha.pickBestStartPoints <- function(blocks, best){
  returnPoints <- list()
  nBlocks <- length(blocks)
  sortedL <- sort(sapply(1:nBlocks, function(x) blocks[[x]]@score),
      decreasing=TRUE)
  #cat("There used to be",nBlocks,"blocks in the list.\n")
  #print(sortedL)
  #cat("Keeping Block: ")
  if (best>length(blocks)){
    stop("There are only ",length(blocks)," blocks to choose from not ",best,"!")
  }else{
    for (s in 1:best){
      for (p in seq(along = blocks)){
        if (sortedL[s] == blocks[[p]]@score){
          returnPoints <- c(returnPoints,blocks[[p]])
          #cat(p," ")
        }else{}
      }
    }
   # cat("\n")
  }
  
  return(returnPoints)
}


## Function to convert 'value' into a 'newBase'-system.  'expo'
## determines the length of the return vector, i.e. how many positions
## the result has. Each position has value: ('newBase'^('expo'-1)).
## value: [0 .. ('newBase'^'expo')-1]
.index2blocks <- function(value,newBase,expo){
  if(value>= newBase^expo){
    print(list(ERROR="Value is too big, i.e. not convertible! \n"))
  }
  else{
    res <- c()
    expo <- expo - 1
    while (expo>-1){
      pos <- newBase^expo
      whole <- floor(value/pos)
      res <- c(res,whole)
      value <- value- whole*pos
      expo <- expo -1 
    }
    if (value!=0) print(list(ERROR="Value is not convertible! \n"))
    else res
  }
}



is.jaatha <- function(jObject){
  return(class(jObject)[1] == "Jaatha")
}


#' Print Start points
#'
#' Method to print the start Points given by an initial Jaatha
#' search sorted by score.
#'
#' @param jObject The Jaatha options
#' @return a matrix with score and parameters of each start point
#' @export
Jaatha.getStartingPoints <- function(jObject){
  startPoints <- jObject@starting.positions
  width <- jObject@nPar + 1
  mat <- matrix(0,length(startPoints),width)
  col.names <- c("score", jObject@par.names)
  colnames(mat) <- col.names

  for (i in 1:length(startPoints)){
    mat[i,1] <- round(startPoints[[i]]@score,2)
    mat[i,-1] <- round(.deNormalize(jObject,t(startPoints[[i]]@MLest))[1:(width-1)], 3)
  }
  perm <- sort.list(mat[,1],decreasing=T) 
  return(mat[perm,])
}

#' Gives the best estimates after a Jaatha search
#'
#' This method extracts the best estimates with log composite likelihood
#' vales from an Jaatha object.
#'
#' @param jObject The Jaatha options
#' @param max.entries If given, no more than this number of entries will be 
#'                returned.
#' @return A matrix with log composite likelihoods and parameters of The
#' best estimates
#' @export
Jaatha.getLikelihoods <- function(jObject, max.entries=NULL) {
  .log3("Called Jaatha.getLikelihoods")
  lt <- jObject@likelihood.table
  lt[,-(1:2)] <- .deNormalize(jObject, lt[,-(1:2), drop=F])
  perm <- sort.list(lt[,1],decreasing=T)  
  lt <- lt[perm, , drop=F]
  .log3("Finished Jaatha.getLikelihoods")
  return(lt[1:min(max.entries, nrow(lt)), , drop=F])
}


printBestPar <- function(jObject, block) {
  .print("Best parameters: ", 
           round(.deNormalize(jObject, t(block@MLest)), 3), 
           "| Score:",  block@score)
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
