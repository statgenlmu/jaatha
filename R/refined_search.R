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
#' @param half.block.size The size of the new block that is created around a new maximum.
#' @param max.steps The search will stop at this number of steps if not stopped
#'        before (see \code{epsilon}).
#' @param rerun You can repeat a previously done refined search in Jaatha.
#'        Do do so, just call the refined search function with the jaatha 
#'        object result of the first refined search and set rerun to "TRUE".
#'              
#' @return An Jaatha object. The found values are written to the slot likelihood.table.
#' 
#' @include block.R
#' @export
Jaatha.refinedSearch <- function(jaatha, best.start.pos=2,
                                 sim=length(getParNames(jaatha))*25,
                                 sim.final=min(sim, 100), half.block.size=.025,
                                 max.steps=200, rerun=FALSE) {

  if (rerun) {
    if( is.null(jaatha@calls[["refined.search"]]) ) 
      stop("No arguments found. Did you run the refined search before?")

    best.start.pos <- jaatha@calls[["refined.search"]]$best.start.pos
    sim <- jaatha@calls[["refined.search"]]$sim
    sim.final <- jaatha@calls[["refined.search"]]$sim.final
    half.block.size <- jaatha@calls[["refined.search"]]$half.block.size
    max.steps <- jaatha@calls[["refined.search"]]$max.steps

  } else {
    arguments <- list(best.start.pos=best.start.pos,
                      sim=sim,
                      sim.final=sim.final,
                      half.block.size=half.block.size,
                      max.steps=max.steps)
    jaatha@calls[["refined.search"]] <- arguments 
  }

  # Check parameters

  assert_that(is_jaatha(jaatha))
  assert_that(is_single_numeric(best.start.pos))
  assert_that(is_single_numeric(sim))
  assert_that(is_single_numeric(sim.final))
  assert_that(is_single_numeric(half.block.size))
  assert_that(is_single_numeric(max.steps))

  if (nrow(jaatha@likelihoods_is) == 0) 
    stop("No starting positions available. Did you run a initial search first?")

  jaatha@likelihoods_rs <- create_likelihood_table(jaatha, 0)

  # Setup environment for the refined search
  set.seed(jaatha@seeds[3])

  # Start a search for every start point
  for (s in 1:best.start.pos) {
    jaatha <- refinedSearchSingleBlock(jaatha, start_pos=s, 
                                       sim=sim, sim.final=sim.final,
                                       half.block.size=half.block.size,
                                       max.steps=max.steps)  
  }

  .print()
  .print("Best log-composite-likelihood values are:")
  print(Jaatha.getLikelihoods(jaatha, 5))

  return(jaatha)
}


## This is called from Jaatha.refinedSearch for each block. The actual search is done here.
## Parameters are the same as in Jaatha.refinedSearch
refinedSearchSingleBlock <- function(jaatha, start_pos, sim, sim.final,
                                     half.block.size, 
                                     max.steps=max.steps){
  
  .print("*** Search with starting point ", start_pos, " ***")
  
  ## initialize values 
  step.current <- 0
  estimate <- sort_likelihood_table(jaatha@likelihoods_is)[start_pos, -(1:2)]
  score = -Inf
  likelihoods <- create_likelihood_table(jaatha, max.steps)

  sim.saved <- list()

  optimum.li <- -Inf
  optimum.step <- 0

  repeat{
    .print("-----------------")
    step.current <- step.current + 1
    .print("Step No ", step.current)

    # Update the Search Block
    search.block <- block_class$new(calcBorders(estimate, half.block.size))

    # Simulate Data and Fit Model
    # If Glm does not converge, try using more simulations
    for (j in 1:5) {
      sim.data <- simulateWithinBlock(sim, search.block, jaatha)
      sim.saved <- getReusableSimulations(search.block, jaatha, sim.saved,
                                          sim.data, step.current)
      tryCatch({
        # Fit the GLM
        glm.fitted <- fit_glm(jaatha, sim.saved)
        break
      }, error = function(e) {
        if (j < 5) .print("Failed to fit the GLM. Retrying with more simulations...")
        else stop("Failed to fit the GLM. Try disabeling smoothing")
      })
    }
    stopifnot(exists("glm.fitted"))

    # Update likelihood of last steps estimate, based on new simulated data.
    # Should be a bit more accurate as previous estimate of the likelihood,
    # as it is in the center of the block now.
    score <- estimateLogLikelihood(estimate, glm.fitted, jaatha@sum_stats, 
                                   getScalingFactor(jaatha))
    
    # Keep track of the optimal likelihood, and in which step it has be reached.
    if (score > optimum.li) {
      optimum.li <- score
      optimum.step <- step.current - 1
    }

    # Record the parameter combination
    likelihoods[step.current, ] <- c(score, start_pos, estimate)

    # Output current best position
    printBestPar(estimate, score, jaatha)

    # Stop criterion 1: maximal likelihood did not increase in the last 10 steps
    if ( step.current >= optimum.step + 10 ) {
      .print()
      .print("*** Finished search ***")
      .print("Seems we can not improve the likelihood anymore.")
      break
    } 

    # Stop criterion 2: more than max.steps search steps
    if ( step.current >= (max.steps) ) {
      .print()
      .print("Maximimum number of search steps",max.steps,"reached.\n")
      .print()
      break
    }

    # Estimate the best parameters in the current block.
    estimate <- search_best_par(search.block, glm.fitted, 
                                   jaatha@sum_stats, 
                                   getScalingFactor(jaatha))$est 
    .print()
  }

  likelihoods <- sort_likelihood_table(likelihoods, 10)
  
  .print("Calculating log-likelihoods for best estimates:")

  for (t in 1:nrow(likelihoods)) {
    .print("* Parameter combination", t, "of", nrow(likelihoods))
    likelihoods[t, 1] <- simLikelihood(jaatha, sim.final, likelihoods[t, -(1:2)])
  }
  
  jaatha@likelihoods_rs <- rbind(jaatha@likelihoods_rs, likelihoods)

  .print()
  return(jaatha)
}


## Function to determine the boarders (i.e. "around" to either side)
## of the new block where to continue the search. The new boarders are
## determined by subtracting and adding "around" to each dimension of
## point.
calcBorders <- function(point, radius) {
  if (any(point < 0 || point > 1)) {
    # Due to foalting point precisition, the point may be slight outside of the 
    # parmeter space. Move it to the border if it is.
    point[point < 0-1e-10] <- 0
    point[point > 1+1e-10] <- 1
    
    # Something is wrong if it is more than 1e-10 ouside the parameter space
    if (any(point < 0 || point > 1)) {
      stop("Point coordinates outside [0,1]: ", point[1], "x", point[2])
    }
  }
  
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
                     function(x) block$includes(normalize(x$pars, jaatha)))]) 
}


