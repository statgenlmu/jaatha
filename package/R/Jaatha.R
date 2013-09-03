# --------------------------------------------------------------
# Jaatha.R
# This file contains the main algorithm of Jaatha. 
# 
# Authors:  Lisha Naduvilezhath & Paul R. Staab
# Date:     2012-10-05
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
                  sim.package.size, use.shm) {

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


## Function to save the ten best parameters along the search path with
## their likelihoods.
.saveBestTen <- function (currentTopTen,numSteps,newOptimum){
  if (numSteps<10){  # the first 9 estimates are kept
    currentTopTen[numSteps,] <- c(newOptimum$score,newOptimum$est)
  } else{  # the minimum score in the array is smaller than the new score -> replace
    minScore <- min(currentTopTen[1:9,1])
    #print(minScore)
    if(minScore < newOptimum$score){
      minIndex <- (1:9) [currentTopTen[,1] == minScore]
      #cat("minIndex:",minIndex,minScore,"\n")
      currentTopTen[minIndex,] <- c(newOptimum$score,
                                    newOptimum$est)
    }#else{}
  }
  return(currentTopTen)
}


## Function to concatenate all parNsumstat-fields of the blocks in
## 'blockList' with newParNsumstat, so they can later be used for the
## glm-fitting.
.concatWithPrevSumstats <- function(newParNss,blockList){
  data <- newParNss
  if(length(blockList)!=0){
    for (b in seq(along = blockList)){
      data <- rbind(data,blockList[[b]]@parNsumstat)
    }
  } else{}   
  return (data)
}


## Function to go through the list of 'blockList' (old blocks) and
## determine which blocks to keep. A block will be kept if 'MLpoint'
## falls into the block.  A list of blocks that contain 'MLpoint' will
## be returned.  Each round a block is kept its weight is halved.
## If verbose=FALSE (default) only number of deleted and kept blocks will be 
## written out. If verbose=TRUE details of the blocks that will be kept and 
## deleted will be written into logFile. 
## 
.findReusableBlocks <- function(MLpoint,blockList,weighOld,verbose=FALSE){  
  reusableBlocks <- list()
  listLen <- 1
  nKeep <- 0
  nDel <- 0
  if (length(blockList)!=0) {
    for (b in seq(along = blockList)){
      if(isInBlock(blockList[[b]],MLpoint)){
        blockList[[b]]@weight <- blockList[[b]]@weight*weighOld
        reusableBlocks[[listLen]] <- blockList[[b]]
        if (verbose){ 
          cat("Keeping BLOCK",b," (lower:",round(blockList[[b]]@lowerBound,3),
            " upper:",round(blockList[[b]]@upperBound,3),")\n")
        } else{}
        nKeep <- nKeep +1 
        listLen <- listLen+1        
      } else{        
        if (verbose){ 
          cat("Deleting BLOCK",b,"(lower:",round(blockList[[b]]@lowerBound,3),
            " upper:",round(blockList[[b]]@upperBound,3),")\n")
        }else{}
        #cat(MLpoint,blockList[[b]]@lowerBound,blockList[[b]]@upperBound)
        #print(isInBlock(blockList[[b]],MLpoint))
        #cat(blockList[[b]]@lowerBound<=MLpoint,MLpoint<=blockList[[b]]@upperBound)
        nDel <- nDel +1
      }
    }
  } else{}
  
  if (!verbose){ 
    cat("Number of blocks kept:",nKeep," / Number of blocks deleted:",nDel,"\n")
  }
  return(reusableBlocks)
}


## Function to determine the boarders (i.e. 'around' to either side)
## of the new block where to continue the search. The new boarders are
## determined by substracting and adding 'around' to each dimension of
## point.
.defineBoarders <- function(point,radius){
  if (radius<0) {
    print(list(ERROR= radius," is a negativ number! Cannot proceed! \n"))}
  else if (!all(point>= 0)) {
    print(list(ERROR= point," has negativ part(s)! Cannot proceed! \n"))}
  else if (!all(point<= 1)) {
    print(list(ERROR= point," > than 1! Cannot proceed! \n"))}
  else  {
    ##cat("point",point)
    dimensions <- length(point)
    newBoarder <- sapply(1:dimensions, function(x) c(-radius,radius)+
              point[x])
    ##cat("newB",newBoarder)
    
    ## if lowerB is smaller than 0, increase upperB that much and set
    ## lowerB to 0 if upperBound is bigger than 1, decrease lowerB and
    ## set upperB to 1 points between (1-'radius') till 1 produce same
    ## boarders: (1-2*'radius') and 1 (and 0 till 'radius': 0 and
    ## (2*'radius'))
    for (d in 1:dimensions){
      if(newBoarder[1,d]<0){ 
        newBoarder[2,d] <- min(max(newBoarder[2,d] +
                    abs(newBoarder[1,d]),0),1)
        newBoarder[1,d] <- 0        
      }
      else if(newBoarder[2,d]>1){
        newBoarder[1,d] <- min(max(newBoarder[1,d] -
                    (newBoarder[2,d]-1),0),1)   
        newBoarder[2,d] <- 1
      }
    }
    return (newBoarder)
  }
}


## Function that divides the parameter range [0-1] by the number of
## nBlocksPerPar. Returns an array which holds the new parameter
## boarders [0-1]. These boarders are the same for all parameters but
## to give the same array output as .calcBlockParRanges for all
## parameters the boundrys are returned.
.calcBlockParRanges01 <- function(nPar,nBlocksPerPar){
  intervalSize <- 1/nBlocksPerPar
  lowers <- (1:nBlocksPerPar-1)*intervalSize
  ##dim=c(#blocks per dimension, start&end)
  pRange <- array(rep(c(lowers,lowers+intervalSize),each=nPar),dim=c(nPar,nBlocksPerPar,2))
    return (pRange)
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


## Sets the MLest of the Jaatha or Block object
Jaatha.setMLest <- function(object,estimate){
  object@MLest <- estimate
  cat("'MLest' set to *",estimate,"*\n")
  return (object)
}

## Returns the value of the slot MLest of the Jaatha or Block object in the
## [0..1]-range.
Jaatha.getMLest01 <- function(object){
  return (object@MLest)
}

## Returns the value of the slot MLest of the Jaatha object in the
## acutal parameter range by calling convert2parRange.
## Only works for Jaatha object.
Jaatha.getMLest <- function(jObject){
  if (class(jObject)[1]!="Jaatha"){
    stop("Jaatha object as input needed!")
  }else{
    est <- Jaatha.getMLest01(object=jObject)  
    #print(est)
    estLen <- length(est)
    #cat("printLen",estLen,"\n")
    #cat(length(getparRange(jObject)),"\n")
    #if(length(getparRange(jObject))!=estLen){
    # stop("Dimensions of parameter range and MLest are not matching!")
    #} else{
      return (.deNormalize(jObject, est))
    #}
  }
}

## Returns the value of the slot nPar of the Jaatha or Block object
Jaatha.getnPar <- function(object){
  return (object@nPar)
}

## Returns the value of the slot logMLmax of the Jaatha object
Jaatha.getMLmax <- function(jObject){
  return (jObject@logMLmax)
}

## Sets the value of the slot nSampIni of the Jaatha object
Jaatha.setnSampIni <- function(jObject, value){
  if (value>0 && !missing(value)){
    jObject@nSampIni <- value
    message("'nSampIni' set to ",value)
  }
  else{
    stop("Please specify a number greater than 0!")
  }
  return (jObject)
}

## Returns the value of the slot nSampIni of the Jaatha object
Jaatha.getnSampIni <- function(jObject){
  return (jObject@nSampIni)
}

## Sets the value of the slot nSampMain of the Jaatha object
Jaatha.setnSampMain <- function(jObject, value){
  if (value>0 && !missing(value)){
    jObject@nSampMain <- value
    message("'nSampMain' set to ",value)
  }
  else{
    stop("Please specify a number greater than 0!")
  }
  return (jObject)
}

## Returns the value of the slot nSampMain of the Jaatha object
Jaatha.getnSampMain <- function(jObject){
  return (jObject@nSampMain)
}

## Sets the value of the slot finiteSites of the Jaatha object
Jaatha.setfiniteSites <- function(jObject,value){
  if (jObject@finiteSites == value){
    message("'finiteSites' is already ",value)
  }
  else{
    jObject@finiteSites <- value
    message("'finiteSites' is set to ",value)
  }
  return (jObject)
}

## Returns the value of the slot finiteSites of the Jaatha object
Jaatha.getfiniteSites <- function(jObject){
  return (jObject@finiteSites)
}

## Returns the value of the slot nTotalSumstat of the Jaatha object
Jaatha.getnTotalSumstat <- function(jObject){
  return (jObject@nTotalSumstat)
}

## Sets the value of the slot nTotalSumstat of the Jaatha object
Jaatha.setnTotalSumstat <- function(jObject, value){
  if (value>0 && !missing(value)){
    jObject@nTotalSumstat <- value
    message("'nTotalSumstat' set to ",value)
  }
  else{
    stop("Please specify a number greater than 0!")
  }
  return (jObject)
}

## Returns the value of the slot nNonJsfsSumstat of the Jaatha object
Jaatha.getnNonJsfsSumstat <- function(jObject){
  return (jObject@nNonJsfsSumstat)
}

## Sets the value of the slot nNonJsfsSumstat of the Jaatha object
Jaatha.setnNonJsfsSumstat <- function(jObject, value){
  if (value>0 && !missing(value)){
    jObject@nNonJsfsSumstat <- value
    message("'nNonJsfsSumstat' set to ",value)
  }
  else{
    stop("Please specify a number greater than 0!")
  }
  return (jObject)
}

## Returns the value of the slot parNames of the Jaatha object
Jaatha.getparNames <- function(jObject){
  return (jObject@parNames)
}

## Returns the value of the slot ssFunc of the Jaatha object
Jaatha.getssFunc <- function(jObject){
  return (jObject@ssFunc)
}

## Sets the value of the slot ssFunc of the Jaatha object
Jaatha.setssFunc <- function(jObject,newFunc){
  if (class(newFunc)!="function"){
    stop("Error: Function provided is not of class 'function'!")
  }
  else{
    jObject@ssFunc <- newFunc
    message("'ssFunc' is set new. Use getssFunc(jaathaObject) to view it. Please check if 'nTotalSumstat' and/or 'nNonJsfsSumstat' need to be reset also.")
    return (jObject)
  }
}

## Calls the garbage collector of R. Has to be called explicitly to make memory free
## twice because otherwise .Last.value allocates some memory.
.emptyGarbage <- function(){
  gc()   #verbose=TRUE
  gc()
}

#################
## Trial CODE####
#################

## version to find blocks and concatenate parNsumstat in one step -> is slower than above version
#time1<-Sys.time()             
#reusableBlocksNparNsumstat <- .findBlocksNConcatenate2(searchBlock@MLest[1:jObject@nPar],newParNsumstat,currentBlocks)
#currentBlocks <- reusableBlocksNparNsumstat[[1]]
#searchBlock@parNsumstat <- reusableBlocksNparNsumstat[[2]]
#time2<-Sys.time()
#cat ("Zeit zusammen2:",as.numeric(time2-time1,units="mins"),"\n")



## Finds blocks in 'blockList' that are reusable and concatenates
## their parNsumstat slots. A block will be kept if 'MLpoint' falls
## into the block. Returns a list with [[1]] new list of blocks and
## [[2]] concatenated-parNsumstat array.
.findBlocksNConcatenate <- function(MLpoint,newParNss,blockList){  
  data <- newParNss
  reusableBlocks <- list()
  if (length(blockList)!=0) {
    rlistLen <- 0
    for (b in 1:length(blockList)){
      if(isInBlock(blockList[[b]],MLpoint)){
        rlistLen <- rlistLen+1
        reusableBlocks[[rlistLen]] <- blockList[[b]]
        data <- rbind(data,blockList[[b]]@parNsumstat)
        cat("Keeping BLOCK",b," (lower:",blockList[[b]]@lowerBound,
            " upper:",blockList[[b]]@upperBound,")\n")
      }
      else{        
        cat("Deleting BLOCK",b," (lower:",blockList[[b]]@lowerBound,
            " upper:",blockList[[b]]@upperBound,")\n")
        #cat(MLpoint,blockList[[b]]@lowerBound,blockList[[b]]@upperBound)
        #print(isInBlock(blockList[[b]],MLpoint))
        #cat(blockList[[b]]@lowerBound<=MLpoint,MLpoint<=blockList[[b]]@upperBound)
      }
    }
  }
  return(list(reusableBlocks,data))
}


## Finds blocks in 'blockList' that are reusable and concatenates
## their parNsumstat slots. A block will be kept if 'MLpoint' falls
## into the block. Returns a list with [[1]] new list of blocks and
## [[2]] concatenated-parNsumstat array.
.findBlocksNConcatenate2 <- function(MLpoint,newParNss,blockList){  
  reusableBlocks <- list()
  data <- newParNss
  if (length(blockList)!=0){     # if there are blocks to check
    ## vector of logicals which blocks contain MLest
    keepers <- sapply(blockList,function(b) isInBlock(b,MLpoint))
    cat("Keeping which blocks?",keepers,"\n")
    ## indices of those blocks
    keeperIndex <- (1:length(blockList))[keepers]
    #cat("Indices",keeperIndex"\n")
    if (length(keeperIndex)!=0){ # if there are blocks to keep
      for (k in 1:length(keeperIndex)) {
        reusableBlocks[[k]] <- blockList[[keeperIndex[k]]]
        data <- rbind(data,blockList[[keeperIndex[k]]]@parNsumstat)
        #cat("Keeping BLOCK",keeperIndex[k]," (lower:",blockList[[keeperIndex[k]]]@lowerBound,
        #     " upper:",blockList[[keeperIndex[k]]]@upperBound,")\n")
      }       
    }
  }
  
  return(list(reusableBlocks,data))
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

scaleDemographicModel <- function(dm, scaling.factor) {
  dm@nLoci <- round(dm@nLoci / scaling.factor)
  return(dm)
}
