# --------------------------------------------------------------
# Contains the initialSearch() and helper functions used nowhere else 
# 
# Authors:  Lisha Naduvilezhath & Paul R. Staab
# Date:     2013-09-04
# Licence:  GPLv3 or later
# --------------------------------------------------------------


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
#' @param rerun You can repeat a previously done initial search in Jaatha.
#'        Do do so, just call the initial search function with the jaatha 
#'        object result of the first initial search and set rerun to 'TRUE'.
#'
#' @return The jaatha object with starting positions
#'
#' @export
Jaatha.initialSearch <- function(jaatha, sim=200, blocks.per.par=2, rerun=FALSE){

  if (rerun) {
    if( is.null(jaatha@calls[['initial.search']]) ) 
      stop("No arguments found. Did you run the initial search before?")
    sim <- jaatha@calls[['initial.search']]$sim
    blocks.per.par <- jaatha@calls[['initial.search']]$blocks.per.par
  } else {
    jaatha@calls[['initial.search']] <- list(sim=sim,
                                             blocks.per.par=blocks.per.par) 
  }

  .log2("Called Jaatha.initialSearch()")
  .log2("sim:", sim, "| blocks.per.par:", blocks.per.par)
  set.seed(jaatha@seeds[2])
  .log2("Seeting seed to", jaatha@seeds[2])
  tmp.dir <- getTempDir(jaatha@use.shm)

  setParallelization(jaatha@cores)

  .print("*** Searching starting positions ***")
  .print("Creating initial blocks ... ")

  firstBlocks <- createInitialBlocks(jaatha@par.ranges, blocks.per.par)

  .log2("Calculating block sizes")
  for (i in seq(along=firstBlocks)){
    .print("*** Block", i, ":", printBorder(firstBlocks[[i]], jaatha))

    .log3("Simulating in block", i)
    sim.data <- simulateWithinBlock(sim, firstBlocks[[i]], jaatha)       

    .log3("Fitting GLM in block", i)        
    glms.fitted <- fitGlm(sim.data, jaatha)

    .log3("Searching optimal values in block",i)
    optimal <- findBestParInBlock(firstBlocks[[i]], glms.fitted, jaatha@sum.stats) 
    
    firstBlocks[[i]]@score <- optimal$score

    firstBlocks[[i]]@MLest <- c(optimal$est, optimal$theta)
    printBestPar(jaatha, firstBlocks[[i]])

    .print()
  }

  jaatha@starting.positions <- firstBlocks
  print(Jaatha.getStartingPoints(jaatha))

  removeTempFiles()

  return(jaatha)
}

createInitialBlocks <- function(par.ranges, blocks.per.par) {
  basic.block <- par.ranges # Just to get dimensions & names
  basic.block[,1] <- 0
  basic.block[,2] <- 1/blocks.per.par

  getIthBlock <- function(i) {
    b <- (.index2blocks(value=i-1, newBase=blocks.per.par, expo=nrow(par.ranges))) / blocks.per.par  
    return(new("Block", border=(basic.block + b)))
  }
  return(lapply(1:blocks.per.par^nrow(par.ranges), getIthBlock)) 
}
