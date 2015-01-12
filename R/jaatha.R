# --------------------------------------------------------------
# Jaatha.R
# This file contains the Jaatha S4-Class and a few related 
# helper functions 
# 
# Authors:  Lisha Mathew & Paul R. Staab
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
#' @keywords package
#' @importFrom phyclust ms
#' @importFrom scrm scrm
#' @importFrom methods new representation 
#' @importFrom parallel mclapply
#' @importFrom Rcpp evalCpp
#' @importFrom reshape2 melt
#' @importFrom checkmate qassert qtest assertClass
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
      route = "list",
      scaling.factor = "numeric"
    ),
)

## constructor method for Jaatha object
init <- function(.Object, sim_func, par_ranges, sum_stats, cores = 1) {
  # Check sim.func
  qassert(sim_func, "f1")
  .Object@simFunc <- sim_func

  # Check par.ranges
  qassert(par_ranges, "M+")
  dim(par_ranges)[2] == 2 || stop("par.ranges must have two columns")
  colnames(par_ranges) <- c("min", "max")
  if (is.null(rownames(par_ranges))) 
    rownames(par_ranges) <- as.character(1:nrow(par_ranges)) 
  .Object@par.ranges <- par_ranges

  # Check sum.stats
  is.list(sum_stats) || stop("sum.stats needs to be a list")
  sapply(sum_stats, function(x) {
    if (!any(class(x) %in% c("Stat_PoiInd", "Stat_PoiSmooth"))) {
      stop("Unknown summary statistic of type ", class(x))
    }
  })
  .Object@sum.stats <- sum_stats

  # Sample seeds
  # Jaatha uses three seeds. The first is the "main seed" used to generate the
  # other two seeds if provided, the second is the seed for the initial search
  # and the refined search.
  .Object@seeds <- sampleSeed(3)

  # Check cores 
  qassert(cores, "R1")
  .Object <- setCores(.Object, cores)

  # Placeholders
  .Object@opts <- list()
  .Object@calls <- list()
  .Object@conf.ints <- matrix()
  .Object@likelihood.table <- matrix()
  .Object@starting.positions <- list()
  
  .Object@scaling.factor <- 1

  return (.Object)
}
setMethod(f="initialize", signature ="Jaatha", definition=init)
rm(init)


#' Initialization of a Jaatha estimation for population genetics
#'
#' This function sets the basic parameters for an analysis with
#' Jaatha and is the first step for each application of it.
#'
#' @param model The demographic model to use
#' @param data  The observed data. Jaatha can use data imported with package 
#'              \pkg{PopGenome}. Please refer the to the vignette 
#'              "The Jaatha HowTo" for more information.
#' @param folded If 'TRUE', Jaatha will assume that the JSFS is folded.
#' @param cores The number of cores to use in parallel. If 0, it tries to
#'              guess the number of available cores and use them all.
#' @param scaling.factor You can use this option if you have a large dataset. If
#'              so, Jaatha only simulates only a fraction 1/scaling.factor of the
#'              dataset and interpolates the missing data.
#' @param smoothing If set to true, Jaatha uses a different way to summaries the
#'              JSFS. Instead of binning certain areas, and fitting a glm per
#'              area, only one glm is fitted for the complete JSFS, and the
#'              position of the different entries is treated as a model
#'              parameter. This feature is still experimental and not
#'              recommended for productive use at the moment.  
#' @param only_synonymous Only use synonymous SNP if set to \code{TRUE}. Requires
#'              to provided \code{data} as a PopGenome "GENOME" object.
#' @return A S4-Object of type jaatha containing the settings
#' @export
Jaatha.initialize <- function(data, model, cores=1, scaling.factor=1,
                              folded=FALSE, smoothing=FALSE, 
                              only_synonymous=FALSE, use_fpc=FALSE,
                              fpc_populations=1:2) {
  
  checkType(model, c("dm", "s"))
  checkType(folded, c("bool", "single"))
  checkType(smoothing, c("bool", "single"))
  checkType(scaling.factor, c("num","single"))
  if (smoothing && folded) 
    stop("You can't use smoothing together with a folded JSFS")
  
  dm <- model
  if (!'jsfs' %in% dm.getSummaryStatistics(dm)) {
    warning("JSFS is not a summary statistic of the model. Adding it.")
    dm <- dm.addSummaryStatistic(dm, 'jsfs')
  }
  
  # Convert the data into a list containing the seg.sites of the different groups
  if ('GENOME' %in% is(data)) {
    checkModelDataConsistency(data, dm)
    data <- convPopGenomeToSegSites(data, only_synonymous)
  }
  if (!is.list(data)) stop('`data` has an unexpected format.')
  
  # ------------------------------------------------------------
  # Create Summary Statistics for each group
  # ------------------------------------------------------------
  sum.stats <- list()
  seg.sites <- NULL
  groups <- dm.getGroups(dm)
  
  for (group in groups) {
    if (length(groups) == 1) {
      grp_name_ext <- ''
      group <- 0
    } else {
      grp_name_ext <- paste0('.', group)
    }
    
    seg.sites <- data[[paste0('seg.sites', grp_name_ext)]]
    if (is.null(seg.sites)) stop('No seg.sites in `data` for group ', group)
    jsfs.value <- calcJsfs(seg.sites, dm.getSampleSize(dm))

    # ------------------------------------------------------------
    # JSFS Summary Statistic
    # ------------------------------------------------------------
    if (!smoothing) {
      sum.stats[[paste0('jsfs', grp_name_ext)]] <- 
        list(method = "poisson.transformed",
             transformation = summarizeJSFS,
             value = jsfs.value,
             data = paste0('jsfs', grp_name_ext))
    
      if (folded) sum.stats$jsfs$transformation <- summarizeFoldedJSFS
    } else {
      sample.size <- dm.getSampleSize(dm)
      warning("Smoothing is experimental")
      model <- paste0("( X1 + I(X1^2) + X2 + I(X2^2) + log(X1) + log(",
                      sample.size[1]+2,
                      "-X1) + log(X2) + log(",
                      sample.size[2]+2,
                      "-X2) )^2")
  
      border.mask <- jsfs.value
      border.mask[, ] <- 0
      border.mask[c(1, nrow(jsfs.value)), ] <- 1
      border.mask[ ,c(1, ncol(jsfs.value))] <- 1
      border.mask <- as.logical(border.mask)
  
      sum.stats[[paste0('jsfs', grp_name_ext)]] <- 
        list(method="poisson.smoothing",
             model=model,
             value=jsfs.value,
             border.mask=border.mask,
             border.transformation=summarizeJsfsBorder,
             data=paste0('jsfs', grp_name_ext))
    }

    # ------------------------------------------------------------
    # FPC Summary Statistic
    # ------------------------------------------------------------
    if (use_fpc) {
      # TODO: Assert that dm contains 'seg.sites' statistic
      for (pop in 1:2) {
        if (pop %in% fpc_populations) {
          sum.stats[[paste0('fpc_pop', pop, grp_name_ext)]] <- 
            Stat_FPC$new(seg.sites, dm, population = pop, group = group)
        }
      }
    }

    # ------------------------------------------------------------
    # PMC Summary Statistic
    # ------------------------------------------------------------
    #if ('pmc' %in% dm.getSummaryStatistics(dm, group)) {
    #  dm <- calcPmcBreaks(dm, seg.sites, group = group)
    #  sum.stats[[paste0('pmc', grp_name_ext)]] <- 
    #    list(method='poisson.transformed', transformation=as.vector,
    #         value=createPolymClasses(seg.sites, dm, group = group),
    #         data = paste0('seg.sites', grp_name_ext))
    #}
  }


  # ------------------------------------------------------------
  # Create the Jaatha object
  # ------------------------------------------------------------
  jaatha <- new("Jaatha", 
                sim_func=function(sim.pars, jaatha) {
                  dm.simSumStats(jaatha@opts[['dm']], sim.pars)
                },
                par_ranges=as.matrix(dm.getParRanges(dm)),  
                sum_stats=sum.stats,
                cores=cores)

  if (scaling.factor != 1) {
    dm <- scaleDemographicModel(dm, scaling.factor)
    jaatha@scaling.factor <- scaling.factor
  }

  jaatha@opts[['dm']] <- dm.finalize(dm)
  jaatha@opts[['jsfs.folded']] <- folded

  invisible(jaatha)
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
  return(mat[perm, , drop=FALSE])
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
         "with estimated log-likelihood ", round(block@score, 3))
}


getParNumber <- function(jaatha) nrow(jaatha@par.ranges) 
getParNames <- function(jaatha) rownames(jaatha@par.ranges) 

getScalingFactor <- function(jaatha) {
  jaatha@scaling.factor
}