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
#' @param epsilon Obsolete. Has no effect anymore and will be remove on next
#'        major release.
#' @param half.block.size The size of the new block that is created around a new maximum.
#' @param weight Obsolete. Has no effect anymore and will be remove on next
#'        major release.
#' @param max.steps The search will stop at this number of steps if not stopped
#'        before (see \code{epsilon}).
#' @param rerun You can repeat a previously done refined search in Jaatha.
#'        Do do so, just call the refined search function with the jaatha 
#'        object result of the first refined search and set rerun to 'TRUE'.
#'              
#' @return An Jaatha object. The found values are written to the slot likelihood.table.
#'
#' @export
Jaatha.refinedSearch <- function(jaatha, best.start.pos=2,
                                 sim=length(getParNames(jaatha))*25,
                                 sim.final=min(sim, 100), epsilon=NULL, half.block.size=.025,
                                 weight=NULL, max.steps=200, rerun=FALSE) {

  if (!is.null(epsilon)) warning('Parameter "epsilon" is obsolete and will be removed soon.')
  if (!is.null(weight)) warning('Parameter "weight" is obsolete and will be removed soon.')

  if (rerun) {
    if( is.null(jaatha@calls[['refined.search']]) ) 
      stop("No arguments found. Did you run the refined search before?")

    best.start.pos <- jaatha@calls[['refined.search']]$best.start.pos
    sim <- jaatha@calls[['refined.search']]$sim
    sim.final <- jaatha@calls[['refined.search']]$sim.final
    half.block.size <- jaatha@calls[['refined.search']]$half.block.size
    max.steps <- jaatha@calls[['refined.search']]$max.steps

  } else {
    arguments <- list(best.start.pos=best.start.pos,
                      sim=sim,
                      sim.final=sim.final,
                      half.block.size=half.block.size,
                      max.steps=max.steps)
    jaatha@calls[['refined.search']] <- arguments 
  }

  # Check parameters

  checkType(jaatha, c("jaatha", "single"))
  checkType(best.start.pos, c("num", "single"))
  checkType(sim, c("num", "single"))
  checkType(sim.final, c("num", "single"), F)
  checkType(half.block.size, c("num", "single"))
  checkType(max.steps, c("num", "single"))


  if (length(jaatha@starting.positions) == 0) 
    stop("No starting positions available. Did you run a initial search first?")
  startPoints <- Jaatha.pickBestStartPoints(blocks=jaatha@starting.positions,
                                            best=best.start.pos)

  jaatha@likelihood.table <- matrix(0,0, getParNumber(jaatha) + 2)

  # Setup environment for the refined search
  set.seed(jaatha@seeds[3])
  setParallelization(jaatha@cores)
  tmp.dir <- getTempDir(jaatha@use.shm)

  # Start a search for every start point
  for (s in 1:length(startPoints)){
    .print("*** Search with starting Point in Block",s,"of",length(startPoints),"****")
    jaatha <- refinedSearchSingleBlock(jaatha, startPoints[[s]]@MLest, sim=sim,sim.final=sim.final,
                                        half.block.size=half.block.size,
                                        max.steps=max.steps,block.nr=s)  
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
                                     half.block.size, 
                                     max.steps=max.steps, block.nr){
  ## initialize values 
  step.current <- 0

  ##since the likelihood estimate for the starting point is
  ##only a very rough estimate we don't keep that value 
  search.block <- new("Block", score=-Inf, MLest=start.point)
  sim.saved <- list()

  ## best ten parameters with score are kept for end evaluation
  topTen <- array(0, dim=c(10,(1+getParNumber(jaatha))), 
                  dimnames= list(1:10, c("score", getParNames(jaatha))))     

  optimum.li <- -Inf
  optimum.step <- 0

  repeat{
    .print("-----------------")
    step.current <- step.current + 1
    .print("Step No", step.current)

    # Update the Search Blocks Border
    search.block@border <- calcBorders(search.block@MLest, radius=half.block.size)

    # Simulate
    sim.data <- simulateWithinBlock(sim, search.block, jaatha)
    sim.saved <- getReusableSimulations(search.block, jaatha, sim.saved,
                                        sim.data, step.current)

    # Fit the GLM
    glm.fitted <- fitGlm(sim.data, jaatha)

    # Update likelihood of last steps estimate, based on new simulated data.
    # Should be a bit more accurate as previous estimate of the likelihood,
    # as it is in the center of the block now.
    search.block@score <- estimateLogLikelihood(search.block@MLest, glm.fitted, 
                                                jaatha@sum.stats)
    
    # Keep track of the optimal likelihood, and in which step it has be reached.
    if (search.block@score > optimum.li) {
      optimum.li <- search.block@score
      optimum.step <- step.current - 1
    }

    # Keep the best 10 parameter combinations with their score
    topTen <- .saveBestTen(topTen, step.current-1, search.block)

    # Output current best position
    printBestPar(jaatha, search.block)

    # Stop criterion 1: maximal likelihood did not increase in the last 10 steps
    if ( step.current >= optimum.step + 10 ) {
      .print()
      .print("*** Finished search ***")
      .print("Seems we can not improve the likelihood anymore.")
      break
    } 

    # Stop criterion 2: more than max.steps search steps
    if ( step.current > (max.steps-1) ) {
      .print()
      .print("Maximimum number of search steps",max.steps,"reached.\n")
      .print()
      break
    }

    # Estimate the best parameters in the current block.
    search.block@MLest <- findBestParInBlock(search.block, glm.fitted, jaatha@sum.stats)$est 
    .print()
  }

  topTen <- topTen[topTen[,1] != 0, ]

  likelihoods <- c()
  .print("Calculating log-composite-likelihoods for best estimates:")
  for (t in 1:nrow(topTen)){
    topPar <- topTen[t,-1]
    .print("* Parameter combination",t,"of", nrow(topTen))
    likelihoods[t] <- simLikelihood(jaatha, sim.final, topPar)
  }

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
.saveBestTen <- function (currentTopTen, numSteps, search.block){
  if (numSteps <= 10){  # the first 9 estimates are kept
    currentTopTen[numSteps, ] <- c(search.block@score, search.block@MLest)
  } else{
    minScore <- min(currentTopTen[ ,1])
    if(minScore < search.block@score){
      minIndex <- which.min(currentTopTen[ ,1])
      currentTopTen[minIndex,] <- c(search.block@score,
                                    search.block@MLest)
    }
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
