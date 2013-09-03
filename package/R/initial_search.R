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
  .log2("Called Jaatha.initialSearch()")
  .log2("sim:", sim, "| blocks.per.par:", blocks.per.par)
  set.seed(jaatha@seeds[2])
  .log2("Seeting seed to", jaatha@seeds[2])
  tmp.dir <- getTempDir(jaatha@use.shm)

  setParallelization(jaatha)

  firstBlocks <- list() ## list blocks with simulated summary stats
  ## pRange contains the boarders of all starting blocks
  ## dim=c(#parameter,#blocks per dimension, start&end)

  .log2("Calculating block sizes")
  pRange <- .calcBlockParRanges01(jaatha@nPar,blocks.per.par)  
  nTotalBlocks <- (blocks.per.par^jaatha@nPar)

  .print("*** Determining Starting positions ***")
  .print("Creating", nTotalBlocks, "initial blocks ... ")
  for (i in 1:nTotalBlocks){
    ## 'b' determines which block-index for each parameter
    ## is being considered; dim=#nPar
    b <- .index2blocks(value=i-1, newBase=blocks.per.par, expo=jaatha@nPar) + 1  
    ##+1 bc R indices start with 0

    boundry <- sapply(1:jaatha@nPar, function(p) pRange[p,b[p],])  #dim=c(2,jaatha@nPar)
    boundry.readable <- round(.deNormalize(jaatha, boundry), 3)
    .print("*** Block", i, 
           " (lowerB:", boundry.readable[1, ], 
           " upperB:", boundry.readable[2, ], ")")

    .log3("Creating block",i)
    firstBlocks[[i]] <- new("Block", nPar=jaatha@nPar,
                            lowerBound= boundry[1,], nLoci=70, weight=1,
                            upperBound= boundry[2,],
                            nSamp=sim) 

    .log3("Simulating in block",i)
    parNsumstat <- simulateWithinBlock(firstBlocks[[i]], jaatha)       

    .log3("Fitting GLM in block",i)        
    glm <- glmFitting(parNsumstat, jaatha)

    .log3("Searching optimal values in block",i)
    optimal <- estimate(bObject=firstBlocks[[i]], jaatha,
                        modFeld=glm, boarder=0)

    ##boarder of 0 is important! (otherwise not all param
    ##under consideration)
    firstBlocks[[i]]@score <- optimal$score

    firstBlocks[[i]]@MLest <- c(optimal$est, optimal$theta)
    printBestPar(jaatha, firstBlocks[[i]])

    .print()
  }


  #bestBlockIndex <- which((function(x) max(x)==x)
  #     (sapply(1:nTotalBlocks,function(x) firstBlocks[[x]]@score)))
  #.print("=> Best Block is:",bestBlockIndex,"with score:",
  #     firstBlocks[[bestBlockIndex]]@score,"with\n estimates:\n")
  #print( round( .deNormalize(jObject,firstBlocks[[bestBlockIndex]]@MLest), 3 ) )
  jaatha@starting.positions <- firstBlocks
  print(Jaatha.getStartingPoints(jaatha))

  removeTempFiles()

  return(jaatha)
}
