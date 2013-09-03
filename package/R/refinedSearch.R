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

  ## best ten parameters with score are kept for end evaluation
  topTen <- array(0, dim=c(10,(1+jObject@nPar)), 
                  dimnames= list(1:10, c("score", jObject@par.names)))     

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
    #cat("new Boarder around MLest[0-1]:\n ")
    #print(round(newBoarder,3))             
    searchBlock <- new("Block", lowerBound=newBoarder[1,],
                       upperBound=newBoarder[2,],
                       nPar=jObject@nPar, nSamp=nSim,
                       nLoci=70, weight=1,
                       score=searchBlock@score,
                       MLest=searchBlock@MLest)

    # Simulate
    newParNsumstat <-  simulateWithinBlock(searchBlock, jObject)

    ## use previous simulation results if the MLest is
    ## within that block and fit glm
    currentBlocks <-  .findReusableBlocks(MLpoint=
                                          searchBlock@MLest[1:jObject@nPar],  
                                          blockList=currentBlocks, weighOld=weight)
    ## call R's garbage collector 
    .emptyGarbage()
    ## only for glmFitting step all simulation results are stored in searchBlock
    searchBlock@parNsumstat <-
      .concatWithPrevSumstats(newParNss=newParNsumstat,
                              blockList=currentBlocks)

    currentWeights <- rep(1, nNewSim)
    if (length(currentBlocks)!=0){
      currentWeights <- c(currentWeights, 
                          array(sapply(1:length(currentBlocks),
                                       function(x) rep(currentBlocks[[x]]@weight,nNewSim))))
    }

    glm <- glmFitting(searchBlock@parNsumstat, jObject, currentWeights)

    ## likelihood of old newOptimum parameters based on new simulated data              
    oldParamLikeli <-  .calcScore(param=searchBlock@MLest[1:jObject@nPar], glm,
                                  jObject)

    newOptimum <- estimate(searchBlock, jObject,
                           modFeld=glm, boarder=0)

    ## keep the best 10 parameter combinations with their score
    topTen <- .saveBestTen(currentTopTen=topTen, numSteps=nSteps,
                           newOptimum=newOptimum)

    route[nSteps, ] <- c(newOptimum$score, newOptimum$est)

    ## prepare for next round
    nSteps <- nSteps+1        
    ## in parNsumstat only the newest simulation results
    ## should be kept, the old ones are kept in currentBlocks
    searchBlock@parNsumstat <- newParNsumstat
    searchBlock@MLest <- newOptimum$est
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

  }

  #last best estimates will be taken into topTen
  nBest <- min(nSteps, nrow(topTen))

  ## inlcude last optimum into top ten
  topTen[nBest,] <- c(newOptimum$score,newOptimum$est)

  route <- route[route[,1] != 0,]
  route[,-1] <- .deNormalize(jObject, route[ ,-1,drop=F])
  jObject@route[[length(jObject@route) + 1]] <- route

  likelihoods <- c()
  .log3("Starting final sim.")
  .print("Calulating log-composite-likelihoods for best estimates:")
  for (t in 1:nBest){
    topPar <- topTen[t,-1]    # in original parameter range
    .print("* Parameter combination",t,"of",nBest)
    likelihoods[t] <- Jaatha.calcLikelihood(jObject, nFinalSim, topPar)
    #cat(t,likelihoods[t],"\n")
  }
  .log3("Finished final sim.")

  likelihood.table <- cbind(log.cl=likelihoods,block=blocknr,topTen[topTen[,1]!=0,-1])
  jObject@likelihood.table <- rbind(jObject@likelihood.table,likelihood.table)

  ## in case 2 have the same likelihood only first is given back
  best <- (1:nBest) [max(likelihoods)==likelihoods][1] 
  jObject@logMLmax <- likelihoods[best]
  jObject@MLest <- topTen[best,-1]   # in [0..1] range
  .emptyGarbage()

  .print()
  return (jObject)
}
