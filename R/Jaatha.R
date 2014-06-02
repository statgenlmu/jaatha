# --------------------------------------------------------------
# Jaatha.R
# This file contains the Jaatha S4-Class and a few related 
# helper functions 
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Date:     2013-09-04
# Licence:  GPLv3 or later
# --------------------------------------------------------------


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
#' @title Fast estimation of demographic parameters
#' @keywords package
#' @importFrom phyclust ms
#' @importFrom foreach foreach 
#' @importFrom foreach %do% 
#' @importFrom methods new 
#' @importFrom methods representation 
#' @importFrom parallel mclapply
#' @importFrom Rcpp evalCpp
#' @importFrom dplyr as.tbl_cube
#' @useDynLib jaatha
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
#'    \item{sum.stats}{The observed summary statistics.}
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
#'    \item{starting.positions}{A list of the starting positions, returned by
#'    the initial search}
#'    \item{likelihood.table}{A matrix with the best composite log likelihood values and 
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
      sum.stats = "list",
      seeds="numeric",
      cores = "numeric",
      use.shm = "logical",

      opts = "list",
      calls = "list",
      starting.positions = "list",
      likelihood.table = "matrix",
      conf.ints = "matrix",
      route = "list"
    ),
)

## constructor method for Jaatha object
init <- function(.Object, sim.func, par.ranges, 
                  sum.stats, seed, cores, use.shm = FALSE) {

  # Check sim.func
  checkType(sim.func, c("fun", "s"))
  .Object@simFunc <- sim.func

  # Check par.ranges
  checkType(par.ranges, c("mat", "num"))
  dim(par.ranges)[2] == 2 || stop("par.ranges must have two columns")
  colnames(par.ranges) <- c("min", "max")
  if (is.null(rownames(par.ranges))) 
    rownames(par.ranges) <- as.character(1:nrow(par.ranges)) 
  .Object@par.ranges <- par.ranges

  # Check sum.stats
  is.list(sum.stats) || stop("sum.stats needs to be a list")
  for (i in names(sum.stats)) {
    checkType(sum.stats[[i]]$value, c("num"))
    checkType(sum.stats[[i]]$method, c("char", "s"))

    if (sum.stats[[i]]$method == "poisson.independent") {
      sum.stats[[i]]$transformation <- as.vector 
      sum.stats[[i]]$value.transformed <- as.vector(sum.stats[[i]]$value)
    }
    else if (sum.stats[[i]]$method == "poisson.transformed") {
      checkType(sum.stats[[i]]$transformation, c("fun", "s"))
      sum.stats[[i]]$value.transformed <- sum.stats[[i]]$transformation(sum.stats[[i]]$value)
    }
    else if (sum.stats[[i]]$method == "poisson.smoothing") {
      checkType(sum.stats[[i]]$model, c("char", "s"))      
      stopifnot(length(dim(sum.stats[[i]]$value)) == 2)
      if (!is.null(sum.stats[[i]]$border.transformation)) {
        stopifnot(!is.null(sum.stats[[i]]$border.mask))
        sum.stats[[i]]$border.transformed <- 
          sum.stats[[i]]$border.transformation(sum.stats[[i]]$value)
      }
    }
    else {
      stop("Unknown summary statistic type: ", sum.stats[[i]]$method)
    }
  }
  .Object@sum.stats <- sum.stats

  # Check seed
  # Jaatha uses three seeds. The first is the "main seed" used to generate the
  # other two seeds if provided, the second is the seed for the initial search
  # and the refined search.
  if (missing(seed) || length(seed) == 0) seed <- generateSeeds(3)
  else if (length(seed) == 1) {
    checkType(seed, "num")
    set.seed(seed)
    seed[2:3] <- generateSeeds(2)
  }
  else stop("Malformated argument: seed")
  .Object@seeds <- seed

  # Check use.shm
  checkType(use.shm, c("bool", "s"))
  .Object@use.shm <- use.shm

  # Check cores 
  if (missing(cores)) cores <- 1
  checkType(cores, c("num","single"))
  if (cores > 1) setParallelization(cores)
  .Object@cores <- cores

  # Placeholders
  .Object@opts <- list()
  .Object@calls <- list()
  .Object@conf.ints <- matrix()
  .Object@likelihood.table <- matrix()
  .Object@starting.positions <- list()

  return (.Object)
}
setMethod(f="initialize", signature ="Jaatha", definition=init)
rm(init)


#' Initialization of a Jaatha estimation for population genetics
#'
#' This function sets the basic parameters for an analysis with
#' Jaatha and is the first step for each application of it.
#'
#' @param demographic.model The demographic model to use
#' @param jsfs Your observed Joint Site Frequency Spectrum (JSFS). Jaatha uses
#'        the JSFS as summary statistics.   
#' @param folded If 'TRUE', Jaatha will assume that the JSFS is folded.
#' @param seed An integer used as seed for both Jaatha and the simulation software
#' @param cores The number of cores to use in parallel. If 0, it tries to
#'              guess the number of available cores and use them all.
#' @param scaling.factor You can use this option if you have a large dataset. If
#'              so, Jaatha only simulates only a fraction 1/scaling.factor of the
#'              dataset and interpolates the missing data.
#' @param use.shm Logical. Many modern linux distributions have a shared memory
#'              file system available under /dev/shm. Set this to TRUE to use it for
#'              temporary files. Usually gives a huge performance boost.
#'              Warning: This option will be removed in a future version of
#'              jaatha. The cleaner way to achieve this is to move your complete
#'              R-temp directory to the memory drive. This is explained on 
#'              http://www.paulstaab.de/2013/11/r-shm .
#' @param smoothing If set to true, Jaatha uses a different way to summaries the
#'              JSFS. Instead of binning certain areas, and fitting a glm per
#'              area, only one glm is fitted for the complete JSFS, and the
#'              position of the different entries is treated as a model
#'              parameter. This feature is still experimental and not
#'              recommended for productive use at the moment.  
#' @return A S4-Object of type jaatha containing the settings
#' @examples
#' dm <- dm.createThetaTauModel(c(20,25), 100) 
#' jsfs <- matrix(rpois(21*26, 5), 21, 26)
#' jaatha <- Jaatha.initialize(dm, jsfs) 
#' 
#' @export
Jaatha.initialize <- function(demographic.model, jsfs,
                              seed, cores=1, scaling.factor=1,
                              use.shm=FALSE, folded=FALSE, 
                              smoothing=FALSE) {

  checkType(demographic.model, c("dm", "s"))
  checkType(folded, c("bool", "single"))
  checkType(smoothing, c("bool", "single"))
  checkType(scaling.factor, c("num","single"))
  if (smoothing && folded) 
    stop("You can't use smoothing together with a folded JSFS")

  if (missing(seed)) seed <- numeric()

  if (use.shm) warning("'use.shm' will be removed in a future version of Jaatha.
                       Manually move your complete R-tmp to your memory disk
                       instead. See http://www.paulstaab.de/2013/11/r-shm")

  sum.stats <- list()
  groups <- dm.getGroups(demographic.model)

  for (group in groups) {
    if (group == 1 & all(demographic.model@features$group == 0)) {
      name <- 'jsfs'
      group <- 0
      if (is.list(jsfs)) { 
        jsfs.cur <- jsfs$jsfs
      } else {
        jsfs.cur <- jsfs
      }
    } else {
      name <- paste('jsfs', group, sep='.')
      stopifnot(is.list(jsfs))
      stopifnot(is.matrix(jsfs[[name]]))
      jsfs.cur <- jsfs[[name]]
    }
    stopifnot(is.matrix(jsfs.cur))

    if (!smoothing) {
      sum.stats[[name]] <- list(method="poisson.transformed",
                                  transformation=summarizeJSFS,
                                  value=jsfs.cur)

      if (folded) sum.stats$jsfs$transformation <- summarizeFoldedJSFS
    } else {
      sample.size <- dm.getSampleSize(demographic.model, group)
      warning("Smoothing is still very experimental")
      model <- paste0("( X1 + I(X1^2) + X2 + I(X2^2) + log(X1) + log(",
                      sample.size[1]+2,
                      "-X1) + log(X2) + log(",
                      sample.size[2]+2,
                      "-X2) )^2")

      border.mask <- jsfs.cur
      border.mask[, ] <- 0
      border.mask[c(1, nrow(jsfs.cur)), ] <- 1
      border.mask[ ,c(1, ncol(jsfs.cur))] <- 1
      border.mask <- as.logical(border.mask)

      sum.stats[[name]] <- list(method="poisson.smoothing",
                                        model=model,
                                        value=jsfs.cur,
                                        border.transformation=summarizeJsfsBorder,
                                        border.mask=border.mask)
    }
  }

  jaatha <- new("Jaatha", 
                sim.func=function(sim.pars, jaatha)
                  dm.simSumStats(jaatha@opts[['dm']], sim.pars, names(jaatha@sum.stats)), 
                par.ranges=as.matrix(dm.getParRanges(demographic.model)),  
                sum.stats=sum.stats,
                seed=seed,
                cores=cores,
                use.shm=use.shm)

  if (scaling.factor != 1) {
    demographic.model <- scaleDemographicModel(demographic.model, scaling.factor)
    jaatha@opts[['scaling.factor']] <- scaling.factor
  }

  jaatha@opts[['dm']] <- finalizeDM(demographic.model)
  jaatha@opts[['jsfs.folded']] <- folded

  return(jaatha)
}


## Shows the content of the slots of the Jaatha object.
.show <- function(object) {
  initial.done <- length(object@starting.positions) > 0
  refined.done <- !is.na(object@likelihood.table[1,1])
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
    print(Jaatha.getStartingPoints(object))
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

Jaatha.pickBestStartPoints <- function(blocks, best){
  returnPoints <- list()
  nBlocks <- length(blocks)
  sortedL <- sort(sapply(1:nBlocks, function(x) blocks[[x]]@score),
      decreasing=TRUE)
  #cat("There used to be",nBlocks,"blocks in the list.\n")
  #print(sortedL)
  #cat("Keeping Block: ")
  if (best>length(blocks)){
    stop("There are only ",length(blocks)," blocks to choose from not ",best,"!")
  }else{
    for (s in 1:best){
      for (p in seq(along = blocks)){
        if (sortedL[s] == blocks[[p]]@score){
          returnPoints <- c(returnPoints,blocks[[p]])
          #cat(p," ")
        }else{}
      }
    }
   # cat("\n")
  }
  
  return(returnPoints)
}


## Function to convert 'value' into a 'newBase'-system.  'expo'
## determines the length of the return vector, i.e. how many positions
## the result has. Each position has value: ('newBase'^('expo'-1)).
## value: [0 .. ('newBase'^'expo')-1]
.index2blocks <- function(value,newBase,expo){
  if(value>= newBase^expo){
    print(list(ERROR="Value is too big, i.e. not convertible! \n"))
  }
  else{
    res <- c()
    expo <- expo - 1
    while (expo>-1){
      pos <- newBase^expo
      whole <- floor(value/pos)
      res <- c(res,whole)
      value <- value- whole*pos
      expo <- expo -1 
    }
    if (value!=0) print(list(ERROR="Value is not convertible! \n"))
    else res
  }
}

#' Print Start points
#'
#' Method to print the start Points given by an initial Jaatha
#' search sorted by score.
#'
#' @param jaatha The Jaatha options
#' @return a matrix with score and parameters of each start point
#' @export
Jaatha.getStartingPoints <- function(jaatha){
  checkType(jaatha, "jaatha")
  mat <- t(sapply(jaatha@starting.positions, 
                  function(x) round(c(log.likelihood=x@score,
                                      denormalize(x@MLest, jaatha)), 3)) )

  perm <- sort.list(mat[,1],decreasing=T) 
  return(mat[perm,])
}

#' Gives the best estimates after a Jaatha search
#'
#' This method extracts the best estimates with log composite likelihood
#' vales from an Jaatha object.
#'
#' @param jaatha The Jaatha options
#' @param max.entries If given, no more than this number of entries will be 
#'                returned.
#' @return A matrix with log composite likelihoods and parameters of The
#' best estimates
#' @export
Jaatha.getLikelihoods <- function(jaatha, max.entries=NULL) {
  checkType(jaatha, "jaatha")
  lt <- jaatha@likelihood.table
  lt[,-(1:2)] <- t(sapply(1:nrow(lt), function(n) denormalize(lt[n,-(1:2), drop=F], jaatha)))
  perm <- sort.list(lt[,1],decreasing=T)  
  lt <- lt[perm, , drop=F]
  return(lt[1:min(max.entries, nrow(lt)), , drop=F])
}


printBestPar <- function(jaatha, block) {
  .print("Best parameters", 
         round(denormalize(block@MLest, jaatha), 3),
         "with estimated log-likelihood", round(block@score, 3))
}


getParNumber <- function(jaatha) nrow(jaatha@par.ranges) 
getParNames <- function(jaatha) rownames(jaatha@par.ranges) 
