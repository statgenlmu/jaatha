# --------------------------------------------------------------
# refined_search.R
# This file contains refinedSearch() and helper functions used only 
# by it.
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Date:     2013-09-04
# Licence:  GPLv3 or later
# --------------------------------------------------------------

#' Iterative search for the maximum composite likelihood parameters
#'
#' This function searches for the parameter combination with the highest
#' composite likelihood. Therefore it iteratively searches the area (=block) 
#' around the last found value for a new maximum using simulations and a
#' generalized linear model.
#' 
#' @param jaatha The Jaatha settings (create with \code{\link{Jaatha.initialize}})
#' @param best.start.pos This is the number of best starting positions
#'      found in the initial search that we will use. Jaatha runs a separate 
#'      search starting from each of this points.
#' @param sim The number of simulations that are performed in each step
#' @param sim.final The number of simulations that are performed after the search to estimate the 
#'        composite log likelihood. If not specified, the value of \code{sim} will be used
#' @param epsilon The search stops if the improvement of the score is less than this for 5 times in a row. 
#' @param half.block.size The size of the new block that is created around a new maximum.
#' @param weight The weighting factor that will reduce the influence of old block in the estimation procedure
#' @param max.steps The search will stop at this number of steps if not stopped
#' before (see \code{epsilon})
#' @param rerun You can repeat a previously done refined search in Jaatha.
#'        Do do so, just call the refined search function with the jaatha 
#'        object result of the first refined search and set rerun to 'TRUE'.
#'              
#' @return An Jaatha object. The found values are written to the slot likelihood.table.
#'
#' @export
Jaatha.refinedSearch <- function(jaatha, best.start.pos=2,
                                 sim=length(getParNames(jaatha))*25,
                                 sim.final=min(sim, 100), epsilon=2, half.block.size=.025,
                                 weight=.9, max.steps=200, rerun=FALSE) {

  if (rerun) {
    if( is.null(jaatha@calls[['refined.search']]) ) 
      stop("No arguments found. Did you run the refined search before?")

    best.start.pos <- jaatha@calls[['refined.search']]$best.start.pos
    sim <- jaatha@calls[['refined.search']]$sim
    sim.final <- jaatha@calls[['refined.search']]$sim.final
    epsilon <- jaatha@calls[['refined.search']]$epsilon
    half.block.size <- jaatha@calls[['refined.search']]$half.block.size
    weight <- jaatha@calls[['refined.search']]$weight
    max.steps <- jaatha@calls[['refined.search']]$max.steps

  } else {
    arguments <- list(best.start.pos=best.start.pos,
                      sim=sim,
                      sim.final=sim.final,
                      epsilon=epsilon,
                      half.block.size=half.block.size,
                      weight=weight,
                      max.steps=max.steps)
    jaatha@calls[['refined.search']] <- arguments 
  }

  # Check parameters
  .log2("Called function Jaatha.refinedSearch()")

  checkType(jaatha, c("jaatha", "single"))
  checkType(best.start.pos, c("num", "single"))
  checkType(sim, c("num", "single"))
  checkType(sim.final, c("num", "single"), F)
  checkType(epsilon, c("num", "single"))
  checkType(half.block.size, c("num", "single"))
  checkType(weight, c("num", "single"))
  checkType(max.steps, c("num", "single"))


  if (length(jaatha@starting.positions) == 0) 
    stop("No starting positions available. Did you run a initial search first?")
  startPoints <- Jaatha.pickBestStartPoints(blocks=jaatha@starting.positions,
                                            best=best.start.pos)

  jaatha@likelihood.table <- matrix(0,0, getParNumber(jaatha) + 2)

  # Setup environment for the refined search
  set.seed(jaatha@seeds[3])
  .log2("Seeting seed to", jaatha@seeds[3])
  setParallelization(jaatha@cores)
  tmp.dir <- getTempDir(jaatha@use.shm)

  # Start a search for every start point
  for (s in 1:length(startPoints)){
    .print("*** Search with starting Point in Block",s,"of",length(startPoints),"****")
    jaatha <- refinedSearchSingleBlock(jaatha, startPoints[[s]]@MLest, sim=sim,sim.final=sim.final,
                                        epsilon=epsilon,half.block.size=half.block.size,
                                        weight=weight,max.steps=max.steps,block.nr=s)  
  }

  .print()
  .print("Best log-composite-likelihood values are:")
  print(Jaatha.getLikelihoods(jaatha, 5))

  removeTempFiles()
  return(jaatha)
}


## This is called from Jaatha.refinedSearch for each block. The actual search is done here.
## Parameters are the same as in Jaatha.refinedSearch
refinedSearchSingleBlock <- function(jaatha, start.point, sim, sim.final,
                                     epsilon, half.block.size, weight=weight,
                                     max.steps=max.steps, block.nr){
  ## initialize values 
  .log3("Initializing")
  nSteps <- 1
  noLchangeCount <- 0
  lastNoChange <- -1        # has to be !=0 for the start

  # Track route through parameter space
  route <- matrix(0, max.steps, getParNumber(jaatha)+1)
  route[1,  1] <- -1
  route[1, -1] <- denormalize(start.point, jaatha) 

  ##since the likelihood estimate for the starting point is
  ##only a very rough estimate we don't keep that value 
  searchBlock <- new("Block", score=-1e11, MLest=start.point)
  sim.saved <- list()

  ## best ten parameters with score are kept for end evaluation
  topTen <- array(0, dim=c(10,(1+getParNumber(jaatha))), 
                  dimnames= list(1:10, c("score", getParNames(jaatha))))     

  repeat{
    .print("-----------------")
    .print("Step No",nSteps)

    border <- calcBorders(searchBlock@MLest, radius=half.block.size)

    searchBlock <- new("Block", border=border,
                       weight=1,
                       score=searchBlock@score,
                       MLest=searchBlock@MLest)

    # Simulate
    sim.data <- simulateWithinBlock(sim, searchBlock, jaatha)
    sim.saved <- getReusableSimulations(searchBlock, jaatha, sim.saved,
                                        sim.data, nSteps)

    .print("Using", length(sim.saved), "Simulations")
    glm.fitted <- fitGlm(sim.data, jaatha)

    ## likelihood of old newOptimum parameters based on new simulated data              
    oldParamLikeli <- estimateLogLikelihood(searchBlock@MLest[1:getParNumber(jaatha)], 
                                            glm.fitted, jaatha@sum.stats)

    optimal <- findBestParInBlock(searchBlock, glm.fitted, jaatha@sum.stats) 

    ## keep the best 10 parameter combinations with their score
    topTen <- .saveBestTen(currentTopTen=topTen, numSteps=nSteps,
                           newOptimum=optimal)

    route[nSteps, ] <- c(optimal$score, denormalize(optimal$est, jaatha))

    ## prepare for next round
    nSteps <- nSteps+1        
    ## in parNsumstat only the newest simulation results
    searchBlock@MLest <- optimal$est
    searchBlock@score <- optimal$score 

    # Output current best position
    printBestPar(jaatha, searchBlock)

    ## stop criterion 1: likelihood difference less than epsilon 
    if ( (abs(optimal$score - oldParamLikeli) < epsilon)){
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

    ## stop criterion 2: more than max.steps search steps
    if (nSteps>(max.steps-1)) {
      .print()
      .print("Maximimum number of search steps",max.steps,"reached.\n")
      .print()
      break
    }

    .print()
  }

  #last best estimates will be taken into topTen
  nBest <- min(nSteps, nrow(topTen))

  ## inlcude last optimum into top ten
  topTen[nBest,] <- c(optimal$score, optimal$est)

  route <- route[route[,1] != 0,]
  jaatha@route[[length(jaatha@route) + 1]] <- route

  likelihoods <- c()
  .log3("Starting final sim.")
  .print("Calculating log-composite-likelihoods for best estimates:")
  for (t in 1:nBest){
    topPar <- topTen[t,-1]
    .print("* Parameter combination",t,"of",nBest)
    likelihoods[t] <- simLikelihood(jaatha, sim.final, topPar)
  }
  .log3("Finished final sim.")

  likelihood.table <- cbind(log.cl=likelihoods, block=block.nr,topTen[topTen[,1]!=0,-1])
  jaatha@likelihood.table <- rbind(jaatha@likelihood.table,likelihood.table)

  .emptyGarbage()

  .print()
  return(jaatha)
}


## Function to determine the boarders (i.e. 'around' to either side)
## of the new block where to continue the search. The new boarders are
## determined by subtracting and adding 'around' to each dimension of
## point.
calcBorders <- function(point, radius) {
  radius <= 0 && stop("Radius is non-positive")
  any(point < 0 || point > 1) && stop("Point coordinates outside [0,1]")
  
  border <- t(sapply(point, function(x) x+c(-radius,radius)))
  colnames(border) <- c("lower", "upper")

  ## if lowerB is smaller than 0, increase upperB that much and set
  ## lowerB to 0 
  ## if upperBound is bigger than 1, decrease lowerB and
  ## set upperB to 1 
  b.offset <- rep(0, nrow(border))
  b.offset[border[,1]<0] <- -border[,1][border[,1]<0]
  b.offset[border[,2]>1] <- -border[,2][border[,2]>1] + 1
  return (border + b.offset)
}


getReusableSimulations <- function(block, jaatha, sim.saved, sim.data, step.current) {
  if (length(sim.saved) == 0) return(lapply(sim.data, function(x) { x$step <- step.current; x }))
  c(lapply(sim.data, function(x) { x$step <- step.current; x }),
    sim.saved[sapply(sim.saved, 
                     function(x) isInBlock(block, normalize(x$pars, jaatha)))]) 
}



## Function to go through the list of 'blockList' (old blocks) and
## determine which blocks to keep. A block will be kept if 'MLpoint'
## falls into the block.  A list of blocks that contain 'MLpoint' will
## be returned.  Each round a block is kept its weight is halved.
## If verbose=FALSE (default) only number of deleted and kept blocks will be 
## written out. If verbose=TRUE details of the blocks that will be kept and 
## deleted will be written into logFile. 
## 
findReusableBlocks <- function(MLpoint, blockList, weighOld, verbose=FALSE){  
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
        } 
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
  }
  
  if (!verbose){ 
    cat("Number of blocks kept:",nKeep," / Number of blocks deleted:",nDel,"\n")
  }
  return(reusableBlocks)
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


## Calls the garbage collector of R. Has to be called explicitly to make memory free
## twice because otherwise .Last.value allocates some memory.
.emptyGarbage <- function(){
  gc()   #verbose=TRUE
  gc()
}
