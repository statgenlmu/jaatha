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
#' @param conf.level The intended confidence level for the interval.
#' @param replicas The number of Jaatha runs that we perform for calculating 
#'  the confidence interval. Should be reasonable large.
#' @param cores The number of CPU cores that will be used.
#' @param log.folder A folder were log-files for the indiviual runs are placed.
#'  The default is to use a temporary folder, even though it is highly 
#'  recommended to specify a folder and save the logs.
#' @param subset This setting allows you to distribute the CI calculation on
#'  multiple machines. You can specify which replicas should be
#'  run on this machine. Each run is identified by an integer between one and
#'  replicas. All runs which integers are passed as a vector in this arguments
#'  are accutally executed. After all runs have finished, you need to manally
#'  copy all logs into a single folder and call 
#'  \code{\link{Jaatha.getCIsFromLogs}} on this folder.
#' @return The Jaatha Object with confidence intervals included if 'subset' was 
#'  not used. Nothing otherwise.
#' @export
Jaatha.confidenceIntervals <- function(jaatha, conf.level=0.95, 
                                      replicas=100, cores = 1, 
                                      log.folder=tempfile('jaatha-logs'),
                                      subset=1:replicas) {
  
  # Get a seed for each replica plus for simulating data 
  set.seed(jaatha@seeds[1])
  # First two seeds are already used for initial & refined search. 
  # Last seed will be used for generating the bootstrap data.
  seeds <- generateSeeds(replicas+3)[-(1:2)] 

  dir.create(log.folder, showWarnings=FALSE)
  message("Storing logs in ", log.folder)

  est.pars.unscaled <- Jaatha.getLikelihoods(jaatha, 1)[, -(1:2)]
  est.pars <- normalize(est.pars.unscaled, jaatha)
  .print("ML estimates are:", round(est.pars.unscaled, 3), "\n")

  # Simulate data under the fitted model
  .print("Simulating data...")
  set.seed(seeds[length(seeds)])
  sim.pars <- matrix(est.pars, replicas, getParNumber(jaatha), byrow=TRUE)
  sim.data <- runSimulations(sim.pars, cores, jaatha) 
  sum.stats <- lapply(sim.data, convertSimDataToSumStats, 
                      sum.stats=jaatha@sum.stats)

  .print("Conducting Bootstrap Runs...")
  bs.results <- mclapply(subset, rerunAnalysis, 
                         seeds=seeds, jaatha=jaatha, sum.stats=sum.stats,
                         log.folder=log.folder, mc.cores=cores)

  if (all(subset == 1:replicas)) {
    return(invisible(Jaatha.getCIsFromLogs(jaatha, conf.level, log.folder)))
  } 
}

#' Function for calculating bootstrap confidence intervals from logs of a 
#' previous run of Jaatha.confidenceIntervals.
#'
#' @param jaatha The Jaatha object that was used when calling 
#'               Jaatha.confidenceIntervals.
#' @param conf_level The intended confidence level of the interval. The actual
#'               level can vary slightly.
#' @param log_folder The folder with logs from the previous run. Just use the
#'               folder that was given as 'log.folder' argument there.
#' @return The Jaatha Object with confidence intervals included.
#' @export
Jaatha.getCIsFromLogs <- function(jaatha, conf_level=0.95, log_folder) {
  results <- list.files(log_folder, 'run_[0-9]+.Rda$', full.names = TRUE)
  message("Using ", length(results), " completed runs.")
  
  est_pars <- Jaatha.getLikelihoods(jaatha, 1)[, -(1:2)]
  
  bs_estimates <- t(sapply(results, function(result) {
    load(result)
    pars <- Jaatha.getLikelihoods(jaatha, max.entries=1)[,-(1:2)]
    pars_scaled <- normalize(pars, jaatha)
    if (any(pars == 0 | pars == 1))
      warning("Bootstrap estimate hit boundary of parameter space. Confidence Interval might be inaccurate.")
    pars
  }))
  
  jaatha@conf.ints <- t(sapply(1:ncol(bs_estimates), function(i) {
    par.name <- getParNames(jaatha)[i]
    return(calcBCaConfInt(conf_level, bs_estimates[,i], 
                          est_pars[i], length(results)) )
  }))
  rownames(jaatha@conf.ints) <- getParNames(jaatha)
  
  cat("Confidence Intervals are:\n")
  print(jaatha@conf.ints)
  return(invisible(jaatha))
}

rerunAnalysis <- function(idx, jaatha, seeds, sum.stats=NULL, log.folder) {
  message("Starting run ", idx, " ...")

  # Initialize a copy of the jaatha object
  set.seed(seeds[idx])
  jaatha@seeds <- c(seeds[idx], generateSeeds(2))
  sink(paste0(log.folder, "/run_", idx, ".log"))
  if( !is.null(sum.stats) ) jaatha@sum.stats <- sum.stats[[idx]]
  jaatha@cores <- 1

  jaatha <- Jaatha.initialSearch(jaatha, rerun=TRUE) 
  jaatha <- Jaatha.refinedSearch(jaatha, rerun=TRUE)

  save(jaatha, file=paste0(log.folder, "/run_", idx, ".Rda"))
  sink(NULL)

  return(jaatha)
}


convertSimDataToSumStats <- function(sim.data, sum.stats) {
  for (sum.stat in names(sum.stats)) {
    sum.stats[[sum.stat]]$value <- sim.data[[sum.stat]]
    if(!is.null(sum.stats[[sum.stat]]$transformation)) {
      sum.stats[[sum.stat]]$value.transformed <-
        sum.stats[[sum.stat]]$transformation(sum.stats[[sum.stat]]$value) 
    }
  }
  return(sum.stats)
}


calcBCaConfInt <- function(conf.level, bs.values, estimates, replicas) {
  z.hat.null <- calcBiasCorrection(bs.values, estimates, replicas)
  a.hat <- calcAcceleration(bs.values)
  z.alpha <- qnorm(p=c((1-conf.level)/2, 1-(1-conf.level)/2))
  quantiles.corrected <- pnorm(z.hat.null + (z.hat.null + z.alpha) / 
                                 (1-a.hat*(z.hat.null + z.alpha)))  
  conf.int <- quantile(bs.values, probs=quantiles.corrected) 
  names(conf.int) <- c('lower', 'upper')
  return(conf.int)
}


calcBiasCorrection <- function(bs.values, estimates, replicas){
  bias <- sum(bs.values < estimates) / replicas
  # Maybe 0 if the estimate hit the lower boundary
  if (bias == 0) bias <- sum(bs.values <= estimates) / replicas
  if (bias == 0) stop("Error: All bootstrap estimate are larger 
                       than the estimate. Please try more replicas")
  if (bias == 1) stop("Error: All bootstrap estimate are smaller 
                       than the estimate. Please try more replicas")
  bc <- qnorm(bias)
  return(bc)
}


calcAcceleration <- function(bs.values) {
  m <- mean(bs.values)
  nominator <- sum((m-bs.values)^3)
  denom <- 6*(sum((m-bs.values)^2))^(3/2)
  return(nominator/denom)
}
