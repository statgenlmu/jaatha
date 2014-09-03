# --------------------------------------------------------------
# Contains the initialSearch() and helper functions used nowhere else 
# 
# Authors:  Lisha Naduvilezhath & Paul R. Staab
# Date:     2013-09-04
# Licence:  GPLv3 or later
# --------------------------------------------------------------


#' Search the parameter space for good starting positions
#'
#' This functions divides the parameter space in different parts (blocks).
#' In each block, simulations for different parameter combinations are run
#' to roughly predict the combination with the highest score (which is 
#' equivalent to the highest composite log likelihood).
#' This points can later be used as starting positions for the second
#' estimation phase of Jaatha (\code{\link{Jaatha.refinedSearch}})
#'
#' @param jaatha The Jaatha settings (create with \code{\link{Jaatha.initialize}})
#' @param sim An integer stating the number of simulations that are performed for each
#'        block.
#' @param blocks.per.par Integer. For each parameter axis in turn, we divide the
#'        parameter space into \code{blocks.per.par} equally sized blocks by
#'        restricting them to only a fraction of this axis.
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

  set.seed(jaatha@seeds[2])

  firstBlocks <- createInitialBlocks(jaatha@par.ranges, blocks.per.par)

  for (i in seq(along=firstBlocks)){
    .print("*** Block", i, "of", length(firstBlocks), ":", 
           printBorder(firstBlocks[[i]], jaatha))

    sim.data <- list()     
    
    # Simulate Data and Fit Model
    # If Glm does not converge, try using more simulations
    for (j in 1:5) {
      sim.data = c(sim.data, simulateWithinBlock(sim, firstBlocks[[i]], jaatha))
      tryCatch({
        suppressWarnings( glms.fitted <- fitGlm(sim.data, jaatha) )
        break
      }, error = function(e) {
        if (j < 5) .print("Failed to fit the GLM. Retrying with more simulations...")
        else stop('Failed to fit the GLM. Try disabeling smoothing or using more simulations')
      })
    }

    optimal <- findBestParInBlock(firstBlocks[[i]], glms.fitted, jaatha@sum.stats) 
    
    firstBlocks[[i]]@score <- optimal$score

    firstBlocks[[i]]@MLest <- c(optimal$est, optimal$theta)
    printBestPar(jaatha, firstBlocks[[i]])

    .print()
  }

  jaatha@starting.positions <- firstBlocks
  print(Jaatha.getStartingPoints(jaatha))

  return(jaatha)
}

createInitialBlocks <- function(par.ranges, blocks.per.par) {
  basic.block <- par.ranges # Just to get dimensions & names
  basic.block[,1] <- 0
  basic.block[,2] <- 1
  
  if (blocks.per.par == 1) {
    return(list(new("Block", border=basic.block)))
  }

  blocks <- vector("list", nrow(par.ranges)*blocks.per.par)
  for (i in 1:nrow(par.ranges)) {
    for (j in 1:blocks.per.par) {
      new.block <- basic.block
      new.block[i,2] <- 1/blocks.per.par
      new.block[i, ] <- new.block[i, ] + (j-1)/blocks.per.par
      blocks[(j-1)*nrow(par.ranges)+i] <- new("Block", border=new.block)
    }
  }

  blocks
}
