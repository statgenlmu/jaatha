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
#'    \item{externalTheta}{If TRUE theta will estimated using Watersons-estimator}
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
      dm="DemographicModel",
      nPar = "numeric", 
      logMLmax = "numeric",
      nTotalSumstat = "numeric",
      nNonJsfsSumstat = "numeric",
      popSampleSizes= "numeric", 
      nLoci = "numeric",
      MLest="numeric", 
      seeds="numeric",
      externalTheta= "logical", 
      finiteSites= "logical",
      parNames = "character",
      sumStats = "numeric",
      starting.positions = "list",
      likelihood.table = "matrix",
      sum.stats.func = "function",
      sim.package.size = "numeric",
      cores = "numeric",
      scaling.factor = "numeric",
      route = "list",
      use.shm = "logical"
    ),
)

## constructor method for Jaatha object
.init <- function(.Object, demographicModel=NA, jsfs=NA, folded=F, seed=numeric(),
                  summary.statistics=NA, cores, sim.package.size, scaling.factor,
                  use.shm) {

  .log3("Starting initialization")
  # Set demographic model
  .Object@dm <- finalizeDM(demographicModel)
  
  # Model parameters
  # TODO: remove.
  .Object@nLoci <- demographicModel@nLoci
  .Object@popSampleSizes <- demographicModel@sampleSizes
  .Object@externalTheta <- demographicModel@externalTheta
  .Object@finiteSites <- demographicModel@finiteSites
  .Object@nPar <- dm.getNPar(.Object@dm)

  # Use summary statistics for folded JSFS?
  if (folded) .Object@sum.stats.func <- Jaatha.defaultFoldedSumStats 
  else        .Object@sum.stats.func <- Jaatha.defaultSumStats 

  # Calculate summary statistics of real data if needed
  if (any(is.na(summary.statistics)) & any(is.na(jsfs))) 
    stop("Either summary.statistics or jsfs must be given")
  if (any(!is.na(summary.statistics)) & any(!is.na(jsfs))) 
    stop("Only summary.statistics or jsfs can be used, but not both")

  if (!any(is.na(jsfs))) {
    .log2("JSFS given. Calculating summary statistics")
    summary.statistics <- .Object@sum.stats.func(jsfs = jsfs)
  }

  if ( is.matrix(summary.statistics) && dim(summary.statistics)[1] == 1) 
    summary.statistics <- as.vector(summary.statistics)

  .log2("summary.statistics:",summary.statistics)
  .Object@sumStats <- summary.statistics
  .Object@nTotalSumstat <- length(summary.statistics)

  .Object@likelihood.table <- matrix()
  .Object@starting.positions <- list()

  .Object@use.shm <- use.shm

  # Seeds
  # Jaatha uses three seeds. The first is the "main seed" used to generate the
  # other two seeds if provided, the second is the seed for the initial search
  # and the refined search.
  if (length(seed) == 0) seed <- generateSeeds(3)
  else if (length(seed) == 1) {
    set.seed(seed)
    seed[2:3] <- generateSeeds(2)
  }
  else stop("Malformated argument: seed")
  .Object@seeds <- seed


  #--------------------------------------------------------------------------
  # Parallelization options
  #--------------------------------------------------------------------------
  if (missing(sim.package.size)) sim.package.size <- 10
  checkType(sim.package.size, c("num","single"))
  .Object@sim.package.size <- sim.package.size

  if (missing(cores)) cores <- 1
  .Object@cores <- cores
  
  .log3("Finished initialization")


  #--------------------------------------------------------------------------
  # Scaling options
  #--------------------------------------------------------------------------
  if (missing(scaling.factor)) scaling.factor <- 1
  checkType(scaling.factor, c("num","single"))
  .Object@dm <- scaleDemographicModel(.Object@dm, scaling.factor)
  if (.Object@dm@nLoci < 3) 
    stop("Error: Simulating less then 3 Loci. Use more data or less scaling.")
  .Object@scaling.factor <- scaling.factor
  
  return (.Object)
}

setMethod(f="initialize", signature ="Jaatha",definition=.init)
rm(.init)
    
#' Basic initialization of a Jaatha estimation
#'
#' This function sets the basic parameters for an analysis with
#' Jaatha and is the first step for each application of it.
#'
#' @param demographic.model The demographic model to use
#' @param summary.statistics The summary statistics calculated from the real data
#' @param jsfs Instead of summary.statistics, you can also input the joint site
#'             fequency spectrum of your data. Jaatha will then automatically
#'             calulate summary statistics out of it.
#' @param folded If 'TRUE', Jaatha will assume that the JSFS is folded.
#' @param seed An integer used as seed for both Jaatha and the simulation software
#' @param log.level An integer from 0 to 3 indicating Jaatha's verbosity. 0 is
#'              (almost) no output, 1 is normal output, and 2 and 3 are some and heavy debug
#'              output respectively
#' @param log.file If specified, the output will be redirected to this file
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
#' @param use.shm Logical. Many modern linux distributions have a shared memory ramdisk
#'              available under /dev/shm. Set this to TRUE to use it for
#'              temporary files. Usually gives a huge performance boost.
#' @return A S4-Object of type jaatha containing the settings
#' @export
Jaatha.initialize <- function(demographic.model, summary.statistics, jsfs,
                              folded=FALSE, seed=numeric(), 
                              log.level=1, log.file="", 
                              sim.package.size=10,
                              cores=1, scaling.factor=1,
                              use.shm=FALSE) {

  setLogging(log.level, log.file)

  if (missing(summary.statistics)) summary.statistics <- NA
  if (missing(jsfs)) jsfs <- NA

  jaatha <- new("Jaatha",
                demographicModel=demographic.model,
                summary.statistics=summary.statistics,
                jsfs=jsfs,
                folded=folded,
                seed=seed,
                sim.package.size=sim.package.size,
                cores=cores,
                scaling.factor=scaling.factor,
                use.shm=use.shm)
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


#' Search the parameter space for good starting positions
#'
#' This functions devides the parameter space in different parts (blocks).
#' In each block, simulations for different parameter combinations are run
#' to roughly predict the combination with the highest score (which is 
#' equivalent to the highest composite log likelihood).
#' This points can later be used as starting positions for the secound
#' estimation phase of Jaatha (\code{\link{Jaatha.refinedSearch}})
#'
#' @param jaatha The Jaatha settings (create with \code{\link{Jaatha.initialize}})
#' @param sim Numeric. The number of simulations that are performed in each bin
#' @param blocks.per.par Numeric. The number of block per parameter. 
#'          Will result in Par^blocks.per.par blocks
#'
#' @return The jaatha object with starting positions
#'
#' @export
Jaatha.initialSearch <- function(jaatha, sim=200, blocks.per.par=3){
  jObject <- jaatha
  nSim <- sim 
  nBlocksPerPar <- blocks.per.par
  .log2("Called Jaatha.initialSearch()")
  .log2("nSim:",nSim,"| nBlocksPerPar:",nBlocksPerPar)
  set.seed(jObject@seeds[2])
  .log2("Seeting seed to", jObject@seeds[2])
  tmp.dir <- getTempDir(jObject@use.shm)
  ## change slot values of jObject locally, so that only  
  ## searches with externalTheta=T will be run for initial search
  ## if theta is included into parRange, exclude it and decrese nPar
  extThetaPossible <- !jObject@externalTheta & !jObject@finiteSites & jObject@nPar > 2
  .log2( "extThetaPossible:",extThetaPossible)
  jObject.bu <- jObject
  if ( extThetaPossible ){
    #origParRange <- jObject@parRange
    #jObject@parRange <- jObject@parRange[-jObject@nPar]
    jObject@dm <- dm.setExternalTheta(jObject@dm)
    jObject@nPar <- jObject@nPar - 1
    jObject@externalTheta <- TRUE
    .print("externalTheta set to TRUE for initial search.")
  }
  #print(jObject)

  setParallelization(jObject)
      
  firstBlocks <- list() ## list blocks with simulated summary stats
  ## pRange contains the boarders of all starting blocks
  ## dim=c(#parameter,#blocks per dimension, start&end)
  .log2("Calculation block sizes")
  pRange <- .calcBlockParRanges01(jObject@nPar,nBlocksPerPar)  
  nTotalBlocks <- (nBlocksPerPar^jObject@nPar)
  .print("*** Starting position is being determined ***")
  .print("Creating",nTotalBlocks,"initial blocks ... ")
  for (i in 1:nTotalBlocks){
    ## 'b' determines which block-index for each parameter
    ## is being considered; dim=#nPar
    b <- .index2blocks(value=i-1, newBase=nBlocksPerPar,expo=jObject@nPar) + 1  ##+1 bc R indices start with 0
    boundry <- sapply(1:jObject@nPar, function(p) pRange[p,b[p],])  #dim=c(2,jObject@nPar)
    boundry.readable <- round(.deNormalize(jObject, boundry, withoutTheta=jObject@externalTheta), 3)
    .print("*** Block", i, 
           " (lowerB:", boundry.readable[1, ], 
            " upperB:", boundry.readable[2, ], ")")

    .log3("Creating block",i)
    firstBlocks[[i]] <- new("Block", nPar=jObject@nPar,
            lowerBound= boundry[1,], nLoci=70, weight=1,
            upperBound= boundry[2,],
            nSamp=nSim) 

    .log3("Simulating in block",i)
    firstBlocks[[i]]@parNsumstat <- simulateWithinBlock(bObject=firstBlocks[[i]],
                    jaathaObject=jObject)       
    #print(firstBlocks[[i]]@parNsumstat[,1:jObject@nPar])

    .log3("Fitting GLM in block",i)        
    glm <- glmFitting(bObject=firstBlocks[[i]],
          nTotalSumstat=jObject@nTotalSumstat,
          weighting=rep(1, dim(firstBlocks[[i]]@parNsumstat)[1]))

    .log3("Searching optimal values in block",i)
    optimal <- estimate(bObject=firstBlocks[[i]],jObject=jObject,
            modFeld=glm, ssData=jObject@sumStats, boarder=0, opt=3)
    ##boarder of 0 is important! (otherwise not all param
    ##under consideration)
    firstBlocks[[i]]@score <- optimal$score
    
    #Normalize theta to 0-1 range and ensure that it is inside its parameter range
    if (extThetaPossible) {
      optimal$theta <- min(max(Jaatha.normalize01(getThetaRange(jObject@dm),
                                                  optimal$theta),0),1)
    }
    firstBlocks[[i]]@MLest <- c(optimal$est, optimal$theta)
    printBestPar(jObject, firstBlocks[[i]])
    
    ## parNsumstat will not be needed anymore -> can be
    ## deleted to save memory
    firstBlocks[[i]]@parNsumstat <- array(0,dim=c(0,0))
    .emptyGarbage()
    .print()
  }
  
  #if ( extThetaPossible ) jObject@nPar <- jObject@nPar + 1

  #bestBlockIndex <- which((function(x) max(x)==x)
  #     (sapply(1:nTotalBlocks,function(x) firstBlocks[[x]]@score)))
  #.print("=> Best Block is:",bestBlockIndex,"with score:",
  #     firstBlocks[[bestBlockIndex]]@score,"with\n estimates:\n")
  #print( round( .deNormalize(jObject,firstBlocks[[bestBlockIndex]]@MLest), 3 ) )
  jObject.bu@starting.positions <- firstBlocks
  print(Jaatha.getStartingPoints(jObject.bu, extThetaPossible))

  removeTempFiles()

  return(jObject.bu)
}



#' Iterative search for the maximum composite likelihood parameters
#'
#' This function searches for the parameter combination with the highest
#' composite likelihood. Therefore it iteratively searches the area (=block) 
#' around the last found value for a new maximum using simulations and a
#' generalized linear model.
#' 
#' @param jaatha The Jaatha settings (create with \code{\link{Jaatha.initialize}})
#' @param best.start.pos This is the number of best starting positions
#'      found in the inital search that we will use. Jaatha runs a seperate 
#'      search starting from each of this points.
#' @param sim The number of simulations that are performed in each step
#' @param sim.final The number of simulations that are performed after the search to estimate the 
#'        composite log likelihood. If not specified, the value of \code{nSim} will be used
#' @param epsilon The search stops if the improvement of the score is less than this for 5 times in a row. 
#' @param half.block.size The size of the new block that is created around a new maximum.
#' @param weight The weighting factor that will reduce the influence of old block in the estimation procedure
#' @param max.steps The search will stop at this number of steps if not stopped before (see epsilon)
#'              
#' @return An Jaatha object. The found values are written to the slot likelihood.table.
#'
#' @export
Jaatha.refinedSearch <- 
  function(jaatha, best.start.pos, sim,
           sim.final, epsilon=.2, half.block.size=.05,
           weight=.9, max.steps=200) {

  if (missing(sim.final)) sim.final <- sim

  jObject <- jaatha
  nSim <- sim
  nFinalSim <- sim.final
  halfBlockSize <- half.block.size
  nMaxStep <- max.steps

  # Check parameters
  if (!is.jaatha(jObject)) stop("jObject is not of type Jaatha")
  .log2("Called function Jaatha.refinedSearch()")

  checkType(best.start.pos, c("num", "single"))
  checkType(nSim, c("num", "single"))
  checkType(nFinalSim, c("num", "single"), F)
  checkType(epsilon, c("num", "single"))
  checkType(halfBlockSize, c("num", "single"))
  checkType(weight, c("num", "single"))
  checkType(nMaxStep, c("num", "single"))

  
  if (length(jObject@starting.positions) == 0) 
    stop("No starting positions available. Did you run a initial search first?")
  startPoints <- Jaatha.pickBestStartPoints(blocks=jObject@starting.positions,
                                            best=best.start.pos)
  
  jObject@likelihood.table <- matrix(0,0, jObject@nPar + 2)

  # Setup enviroment for the refined search
  set.seed(jObject@seeds[3])
  .log2("Seeting seed to", jObject@seeds[3])
  setParallelization(jObject)
  tmp.dir <- getTempDir(jObject@use.shm)

  # Start a search for every start point
  for (s in 1:length(startPoints)){
    jObject@MLest <- startPoints[[s]]@MLest
      .print("*** Search with starting Point in Block",s,"of",length(startPoints),"****")
      jObject <- .refinedSearchSingleBlock(jObject,nSim=nSim,nFinalSim=nFinalSim,
                epsilon=epsilon,halfBlockSize=halfBlockSize,
                weight=weight,nMaxStep=nMaxStep,blocknr=s)  
  }

  .print()
  .print("Best log-composite-likelihood values are:")
  print(Jaatha.getLikelihoods(jObject, 5))

  removeTempFiles()
  return(jObject)
}


## This is called from Jaatha.refinedSearch for each block. The actual search is done here.
## Parameters are the same as in Jaatha.refinedSearch
.refinedSearchSingleBlock <- function(jObject, nSim, nFinalSim,
                                     epsilon, halfBlockSize, weight=weight,
                                     nMaxStep=nMaxStep, blocknr){
  ## initialize values 
  .log3("Initializing")
  currentBlocks <- list()
  nSteps <- 1
  noLchangeCount <- 0
  lastNoChange <- -1        # has to be !=0 for the start
  nNewSim <- nSim + 2^jObject@nPar   # no. sim + no. corners      

  # Track route through parameter space
  route <- matrix(0, nMaxStep, jObject@nPar+1)
  route[1,  1] <- -1
  route[1, -1] <- jObject@MLest

  ##since the likelihood estimate for the starting point is
  ##only a very rough estimate we don't keep that value 
  searchBlock <- new("Block", nPar=jObject@nPar, score=-1e11,
                     MLest=jObject@MLest)

  if (jObject@externalTheta) nTotalPar <- jObject@nPar +1
  else nTotalPar <- jObject@nPar

  ## best ten parameters with score are kept for end evaluation
  topTen <- array(0,dim=c(10,(1+nTotalPar)), #dim=(10,likelihood+pars)
                  dimnames= list(1:10, c("score",dm.getParameters(jObject@dm))))     

  ##repeat until likelihood improvement gets smaller than epsilon
  ## 5 times in a row or more than 200 steps used
  repeat{
    .print("-----------------")
    .print("Step No",nSteps)

    ## define parameter range for new block and simulate
    ## within that block newBoarder[1,]=lower und [2,]=upper
    ## Bound values between 0 and 1
    newBoarder <- .defineBoarders(point=
                                  searchBlock@MLest[1:jObject@nPar],
                                  radius=halfBlockSize)  #0.05 halfBlockSize
    ## excludes theta in externalTheta    
    #cat("new Boarder around MLest[0-1]:\n ")
    #print(round(newBoarder,3))             
    searchBlock <- new("Block", lowerBound=newBoarder[1,],
                       upperBound=newBoarder[2,],
                       nPar=jObject@nPar, nSamp=nSim,
                       nLoci=70, weight=1,
                       score=searchBlock@score,
                       MLest=searchBlock@MLest)
    newParNsumstat <-  simulateWithinBlock(bObject=searchBlock,
                                           jaathaObject=jObject)
    ##print(searchBlock@parNsumstat)

    ## use previous simulation results if the MLest is
    ## within that block and fit glm
    currentBlocks <-  .findReusableBlocks(MLpoint=
                                          searchBlock@MLest[1:jObject@nPar],  # excludes theta in externalTheta
                                          blockList=currentBlocks, weighOld=weight)
    ## call R's garbage collector 
    .emptyGarbage()
    ## only for glmFitting step all simulation results are stored in searchBlock
    searchBlock@parNsumstat <-
      .concatWithPrevSumstats(newParNss=newParNsumstat,
                              blockList=currentBlocks)

    currentWeights <- rep(1,nNewSim)
    if (length(currentBlocks)!=0){
      currentWeights <- c(currentWeights, 
                          array(sapply(1:length(currentBlocks),
                                       function(x) rep(currentBlocks[[x]]@weight,nNewSim))))
    }

    glm <- glmFitting(bObject=searchBlock,
                      nTotalSumstat=jObject@nTotalSumstat,
                      weighting=currentWeights)

    ## likelihood of old newOptimum parameters based on new simulated data              
    oldParamLikeli <-  .calcLikelihoodWithModelfeld(param=
                                                    searchBlock@MLest[1:jObject@nPar], # excludes theta in externalTheta
                                                    modelCoefficients=glm,
                                                    observedSS=jObject@sumStats, jObject=jObject, 
                                                    bObject=searchBlock)

    newOptimum <- estimate(bObject=searchBlock,jObject=jObject,
                           modFeld=glm,ssData=jObject@sumStats,boarder=0,opt=3)

    ## keep the best 10 parameter combinations with their score
    topTen <- .saveBestTen(currentTopTen=topTen, numSteps=nSteps,
                           externalTheta=jObject@externalTheta, newOptimum=newOptimum)

    route[nSteps, ] <- c(newOptimum$score, newOptimum$est)

    ## prepare for next round
    nSteps <- nSteps+1        
    ## in parNsumstat only the newest simulation results
    ## should be kept, the old ones are kept in currentBlocks
    searchBlock@parNsumstat <- newParNsumstat
    if (jObject@externalTheta){
      searchBlock@MLest <- c(newOptimum$est,newOptimum$theta)
    } else{
      searchBlock@MLest <- newOptimum$est
    } 
    searchBlock@score <- newOptimum$score 

    currentBlocks[[length(currentBlocks)+1]] <- searchBlock

    # Output current best position
    printBestPar(jObject, searchBlock)

    ## stop criterion 1: likelihood difference less than epsilon 
    if ( (abs(newOptimum$score - oldParamLikeli) < epsilon)){
      #if (lastNoChange==(nSteps-1)){  # if last noChange happend just last step
        noLchangeCount <- noLchangeCount +1           
        .print("No sigificant score changes in the last",noLchangeCount,"Step(s)")
        ## if there has been no change 5 times in a row, stop 
        if (noLchangeCount>4){
          .print()
          .print("*** Finished search ***")
          .print("Score has not change much in the last 5 steps.")
          .print("Seems we have converged.")
          .print()
          break
        } else{}
      } else{ # if first noChange, save old estimates
        noLchangeCount <- 0
      }
    #  lastNoChange <- nSteps
    #} else{}

    ## stop criterion 2: more than nMaxStep search steps
    if (nSteps>(nMaxStep-1)) {
      .print()
      .print("Maximimum number of search steps",nMaxStep,"reached.\n")
      .print()
      break
    }

    .print()

  } #repeat loop end

  ##last best estimates will be taken into topTen
  nBest <- min(nSteps, nrow(topTen))

  ## inlcude last optimum into top ten
  if (jObject@externalTheta){
    topTen[nBest,] <- c(newOptimum$score,newOptimum$est,
                        newOptimum$theta)    # in range[0..1]
  } else{
    topTen[nBest,] <- c(newOptimum$score,newOptimum$est)
  }
  #print(cbind(topTen[1:nBest,1],
  #     .calcAbsParamValue(jObject@dm,topTen[1:nBest,-1]) ))

  route <- route[route[,1] != 0,]
  route[,-1] <- .deNormalize(jObject, route[ ,-1,drop=F])
  jObject@route[[length(jObject@route) + 1]] <- route

  likelihoods <- c()
  .log3("Starting final sim.")
  .print("Calulating log-composite-likelihoods for best estimates:")
  for (t in 1:nBest){
    topPar <- topTen[t,2:(nTotalPar+1)]    # in original parameter range
    .print("* Parameter combination",t,"of",nBest)
    likelihoods[t] <- Jaatha.calcLikelihood(jObject, 
                                            nSimulations=nFinalSim, 
                                            par=topPar)
    #cat(t,likelihoods[t],"\n")
  }
  .log3("Finished final sim.")

  likelihood.table <- cbind(log.cl=likelihoods,block=blocknr,topTen[topTen[,1]!=0,-1])
  jObject@likelihood.table <- rbind(jObject@likelihood.table,likelihood.table)

  ## in case 2 have the same likelihood only first is given back
  best <- (1:nBest) [max(likelihoods)==likelihoods][1] 
  #cat("best",best,"\n")
  jObject@logMLmax <- likelihoods[best]
  #cat("log",jObject@logMLmax,"\n")
  jObject@MLest <- topTen[best,-1]   # in [0..1] range
  #cat("ML",jObject@MLest,"\n")
  .emptyGarbage()
  #if (jObject@finiteSites){ 
  # system(paste("rm treeFile simOutput"))
  #} else{ 
  # system(paste("rm simOutput"))       
  #}
  ## print results in original parameter ranges into file
  #cat(Jaatha.getMLmax(jObject), round(Jaatha.getMLest(jObject),6),
  #   "\n", file=jObject@resultFile, append=T, sep="\t")

  .print()
  return (jObject)
}


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
.saveBestTen <- function (currentTopTen,numSteps,externalTheta,newOptimum){
  if (numSteps<10){  # the first 9 estimates are kept
    if (externalTheta){
      currentTopTen[numSteps,] <- c(newOptimum$score,newOptimum$est,
          newOptimum$theta)
    }else{
      currentTopTen[numSteps,] <- c(newOptimum$score,newOptimum$est)
    }                  
  } else{  # the minimum score in the array is smaller than the new score -> replace
    minScore <- min(currentTopTen[1:9,1])
    #print(minScore)
    if(minScore < newOptimum$score){
      minIndex <- (1:9) [currentTopTen[,1] == minScore]
      #cat("minIndex:",minIndex,minScore,"\n")
      if(externalTheta){
        currentTopTen[minIndex,] <- c(newOptimum$score,
            newOptimum$est,newOptimum$theta)
      }else{
        currentTopTen[minIndex,] <- c(newOptimum$score,
            newOptimum$est)
      } 
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
    if (jObject@externalTheta){
      estLen <- length(est) -1
    } else{
      estLen <- length(est)
    } 
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

## Sets the value of the slot parRange and nPar of the Jaatha or block object.
## if newParNames is specified, parNames will also be set.
#Jaatha.setparRange <- function(object,newParRanges,newParNames){
# if (class(newParRanges)!="list"){
#   stop(cat("Error: newParRanges is not a list! Please change it!\n"))
# }
# else if(!missing(newParNames) && (length(newParRanges)!=length(newParNames))){
#   stop(cat("Error: newParRanges is not of the same length as newParNames! Please change it!\n"))
# }
# else {
#   object@parRange <- newParRanges
#   object@nPar <- length(newParRanges)
#   if (!missing(newParNames)){
#     object@parNames <- newParNames
#   }else{        
#     if (object@externalTheta){ 
#       object@parNames <- c(paste("par",1:object@nPar,sep=""),"theta")
#     }else{
#       p <- 
#       object@parNames <- c(paste("par",1:(object@nPar-1),sep=""),"theta")
#     }
#   } 
#   message("'parNames' is set to *",object@parNames,"*,  'nPar' to *",object@nPar,
#     "*, and 'parRange' set to:\n")
#   print(object@parRange)
#   return (object)
# }
#}

## Returns the value of the slot parRange of the Jaatha object
#Jaatha.getparRange <- function(jObject){
# return (jObject@parRange)
#}

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

## Sets the value of the slot externalTheta of the Jaatha object
Jaatha.setexternalTheta <- function(jObject,value){
  if (jObject@externalTheta == value){
    message("'externalTheta' is already ",value)
  }
  else{
    jObject@externalTheta <- value
    message("'externalTheta' is set to ",value)
  }
  return (jObject)
}

## Returns the value of the slot externalTheta of the Jaatha object
Jaatha.getexternalTheta <- function(jObject){
  return (jObject@externalTheta)
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
#' @param extThetaPossible For internal use only.
#' @return a matrix with score and parameters of each start point
#' @export
Jaatha.getStartingPoints <- function(jObject, extThetaPossible=F){
  startPoints <- jObject@starting.positions
  width <- dm.getNPar(jObject@dm) + 1 + jObject@externalTheta
  mat <- matrix(0,length(startPoints),width)
  col.names <- c("score", dm.getParameters(jObject@dm))
  if (jObject@externalTheta)
    col.names <- c(col.names, getThetaName(jObject@dm))
  colnames(mat) <- col.names

  for (i in 1:length(startPoints)){
    if (jObject@externalTheta & !extThetaPossible) theta <- startPoints[[i]]@MLest[width-1]
    mat[i,1] <- round(startPoints[[i]]@score,2)
    mat[i,-1] <- round(.deNormalize(jObject,t(startPoints[[i]]@MLest))[1:(width-1)], 3)
    if (jObject@externalTheta & !extThetaPossible) mat[i, width] <- round(theta,3)
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
