# --------------------------------------------------------------
# confidence_intervals.R 
# Contains a function for calculating bias corrected bootstrap 
# confidence intervals 
# 
# Authors:  Paul R. Staab
# Date:     2013-10-21
# Licence:  GPLv3 or later
# --------------------------------------------------------------

#' Function for calculating bootstrap confidence intervals for 
#' jaatha estimates.
#' 
#' This functions calculates bias-corrected and accelerated (BCa)
#' bootstrap confidence intervals as described in Section 14.3 of
#' "An introduction to the bootstrap" by Efron and Tibshirani. 
#' Basically, we simulate many datasets under the model using the estimated
#' parameters and do a complete Jaatha estimation procedure on each dataset.
#' We can then use the thereby estimated values
#' to calculate approximate confidence intervals.
#'
#' Warning: This requires a large number of Jaatha runs and should
#' be best executed on a small cluster rather than on a desktop computer.
#'
#' @param jaatha A jaatha object that was returned by Jaatha.refinedSearch.
#' @param conf.level The intended confidence level of the interval. The actual
#' level can vary slightly.
#' @param replicas The number of Jaatha runs that we perform for calculating 
#'  the confidence interval. Should be reasonable large.
#' @param cores The number of CPU cores that will be used.
#' @param log.folder A folder were log-files for the indiviual runs are placed.
#' @return The Jaatha Object with confidence intervals included.
#' @export
Jaatha.confidenceIntervals <- function(jaatha, conf.level=0.95, 
                                      replicas=100, cores = 1, 
                                      log.folder=paste(tempdir(), "jaatha-logs", sep="/")) {
  
  # Get a seed for each replica plus for simulating data 
  set.seed(jaatha@seeds[1])
  seeds <- generateSeeds(replicas+3)[-(1:2)]

  #if (is.null(log.folder)) log.folder <- paste(tempdir(), "jaatha-logs", sep="/") 
  dir.create(log.folder, showWarnings=FALSE)

  # Simulate data under the fitted model
  .print("Simulating data...\n")
  set.seed(seeds[length(seeds)])
  est.pars <- Jaatha.getLikelihoods(jaatha, 1)[-(1:2)]
  sim.pars <- matrix(est.pars, replicas, getParNumber(jaatha), byrow=TRUE)
  sim.data <- jaatha@simFunc(jaatha, sim.pars)

  setParallelization(cores)

  bs.results <- foreach(i=1:replicas, .combine=rbind) %dopar% {
    .print("Starting run", i, "...")
    sink(paste0(log.folder, "/run_", i, ".log"))

    # Initialize a copy of the jaatha object
    set.seed(seeds[i])
    jaatha.rep <- jaatha
    jaatha.rep@seeds <- c(seeds[i], generateSeeds(2))
    jaatha.rep@sumStats <- sim.data[i, ]
    jaatha.rep@cores <- 1

    jaatha.rep <- Jaatha.initialSearch(jaatha.rep, rerun=TRUE) 
    jaatha.rep <- Jaatha.refinedSearch(jaatha.rep, rerun=TRUE)

    save(jaatha.rep, file=paste0(log.folder, "/run_", i, ".Rda"))
    sink()

    return(Jaatha.getLikelihoods(jaatha.rep, 1)[-(1:2)])
  }

  jaatha@conf.ints <- foreach(i=1:ncol(bs.results), .combine=rbind) %do% {
    par.name <- getParNames(jaatha)[i]
    return( calcBCaConfInt(conf.level, bs.results[,i], est.pars[i], replicas) )
  }
  rownames(jaatha@conf.ints) <- getParNames(jaatha)

  .print("\nConfidence Intervals are:")
  print(jaatha@conf.ints)
  return(jaatha)
}

calcBCaConfInt <- function(conf.level, bs.values, estimates, replicas) {
  z.hat.null <- calcBiasCorrection(log(bs.values), log(estimates), replicas)
  a.hat <- calcAcceleration(log(bs.values))
  z.alpha <- qnorm(p=c((1-conf.level)/2, 1-(1-conf.level)/2)) 
  quantiles.corrected <- pnorm(z.hat.null + (z.hat.null + z.alpha) / (1-a.hat*(z.hat.null + z.alpha)))  
  conf.int <- quantile(log(bs.values), probs=quantiles.corrected) 
  names(conf.int) <- c('lower', 'upper')
  return(exp(conf.int))
}

calcBiasCorrection <- function(bs.values, estimates, replicas){
  bias <- sum(bs.values < estimates)
  # Maybe 0 if the estimate hit the lower boundary
  if (bias == 0) bias <- sum(bs.values <= estimates) 
  if (bias == 0) stop("Error: All bootstrap estimate are larger 
                       than the estimate. Please try more replicas")
  if (bias == 1) stop("Error: All bootstrap estimate are smaller 
                       than the estimate. Please try more replicas")
  bc <- qnorm(bias/replicas)
  return(bc)
}

calcAcceleration <- function(bs.values) {
  m <- mean(bs.values)
  nominator <- sum((m-bs.values)^3)
  denom <- 6*(sum((m-bs.values)^2))^(3/2)
  return(nominator/denom)
}
