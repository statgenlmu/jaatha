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
#' @keywords package
#' @importFrom parallel mclapply
NULL


#' The "Jaatha" S4 class saves the basic parameters for a Jaatha estimation
#' procedure
#'
#' Slots:
#' \describe{
#'    \item{simFunc}{Function used for simulating}
#'    \item{par.ranges}{A nx2 matrix stating the ranges for the n model
#'    parameters we want to estimate. The first row gives the lower range for
#'    the parameters, the second row the upper ranges. Each row stands for one
#'    parameter and the row-names will be used as names for the parameters.} 
#'    \item{sumstats}{The observed summary statistics.}
#'    \item{seeds}{A set of random seeds. First one to generate the other two.
#'                 The next one is for the initial search, the last is for the
#'                 refined search}
#'    \item{cores}{The number of CPU cores to use for simulations}
#'    \item{use.shm}{Use the shared memory /dev/shm for temporary files (Linux only)}
#'    \item{opts}{Placeholder for additional arguments, for instance ones that
#'    needs to be passed to simFunc.}
#'    \item{calls}{The function calls for the initial & refined search. These
#'    are used to rerun the searches with the exactly same settings for
#'    generating bootstrap confidence intervals and the likelihood-ratio
#'    statistic.} 
#'    \item{likelihoods_rs}{A matrix with the best composite log likelihood values and 
#'                            corresponding parameters}
#'    \item{conf.ints}{Confidence Intervals for parameter estimates produced by
#'    Jaatha.confidenceIntervals} 
#'    \item{route}{Tracks the best estimates of each step.}
#' }
#'
#' @name Jaatha-class
#' @rdname Jaatha-class
#' @exportClass Jaatha
setClass("Jaatha",
  representation=representation(
      # Settings
      simFunc="function",
      par.ranges = "matrix",
      sum_stats = "list",
      seeds="numeric",
      cores = "numeric",
      use.shm = "logical",

      opts = "list",
      calls = "list",
      likelihoods_rs = "matrix",
      likelihoods_is = "matrix",
      conf.ints = "matrix",
      route = "list",
      scaling_factor = "numeric"
    ),
)

## constructor method for Jaatha object
#' @importFrom methods setMethod
init <- function(.Object, sim_func, par_ranges, sum_stats, 
                 cores = 1, options = list(), scaling_factor = 1,
                 sim_test = TRUE) {
  # Check sim.func
  assert_that(is.function(sim_func))
  .Object@simFunc <- sim_func

  # Check par.ranges
  assert_that(is.matrix(par_ranges))
  nrow(par_ranges) > 0 || stop("No parameters specified")
  ncol(par_ranges) == 2 || stop("par.ranges must have two columns")
  colnames(par_ranges) <- c("min", "max")
  if (is.null(rownames(par_ranges))) 
    rownames(par_ranges) <- as.character(1:nrow(par_ranges)) 
  .Object@par.ranges <- par_ranges

  # Add sumstats
  is.list(sum_stats) || stop("sumstats needs to be a list")
  .Object@sum_stats <- list()
  for (sum_stat in sum_stats) {
    if (!any(class(sum_stat) %in% c("Stat_PoiInd", "Stat_PoiSmooth")))
      stop("Unknown summary statistic of type ", class(sum_stat))
    
    if (sum_stat$get_name() %in% names(.Object@sum_stats)) {
      stop("There is already a summary statistic with name ", 
           sum_stat$get_name())
    }

    .Object@sum_stats[[sum_stat$get_name()]] <- sum_stat
  }

  # Sample seeds
  # Jaatha uses three seeds. The first is the "main seed" used to generate the
  # other two seeds if provided, the second is the seed for the initial search
  # and the refined search.
  .Object@seeds <- sample_seed(3)

  # Check cores 
  assert_that(is.numeric(cores))
  assert_that(length(cores) == 1)
  assert_that(cores %% 1 == 0)
  .Object <- setCores(.Object, cores)

  # Placeholders
  .Object@opts <- options
  .Object@calls <- list()
  .Object@conf.ints <- matrix()
  .Object@likelihoods_rs <- matrix()
  
  .Object@scaling_factor <- scaling_factor
  
  if (sim_test) test_simulation(.Object)

  return(.Object)
}
setMethod(f="initialize", signature ="Jaatha", definition=init)
rm(init)

is_jaatha <- function(jaatha) any("Jaatha" == class(jaatha))




## Shows the content of the slots of the Jaatha object.
.show <- function(object) {
  initial.done <- !is.na(object@likelihoods_is[1,1])
  refined.done <- !is.na(object@likelihoods_rs[1,1])
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
    print(Jaatha.getLikelihoods(object, initial_search = TRUE))
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

printBestPar <- function(estimate, likelihood, jaatha) {
  .print("Best parameters", 
         round(denormalize(estimate, jaatha), 3),
         "with estimated log-likelihood ", round(likelihood, 3))
}


getParNumber <- function(jaatha) nrow(jaatha@par.ranges) 
getParNames <- function(jaatha) rownames(jaatha@par.ranges) 

getScalingFactor <- function(jaatha) {
  jaatha@scaling_factor
}

getStatName <- function(stat, group, pop) {
  if (!missing(pop)) stat <- paste0(stat, "_pop", pop)
  if (group > 0) stat <- paste0(stat, ".", group)
  stat
}
