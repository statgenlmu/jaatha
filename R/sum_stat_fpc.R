#' @importFrom R6 R6Class
Stat_FPC <- R6Class('Stat_FPC', inherit = Stat_PoiInd,
  private = list(
    stat_name = NA,
    breaks = NA,
    trio_classes = NA
  ),
  public = list(
    initialize = function(seg_sites, model, stat, break_probs = c(.2, .5)) {
      if (any(break_probs > 1)) stop('probs greater then one')
      
      # Setup the cube for the middle locus
      private$stat_name <- stat$get_name()
      fpc_percent <- stat$calculate(seg_sites, NULL, model)
      fake_sim_data <- list()
      fake_sim_data[[private$stat_name]] <- fpc_percent
      
      private$breaks = lapply(1:ncol(fpc_percent), function(i) {
        calcBreaks(fpc_percent[, i], break_probs)
      })
      names(private$breaks) <- colnames(fpc_percent)
      
      super$initialize(fake_sim_data, stat$get_name())
    },
    transform = function(sim_data) {
      fpc <- sim_data[[private$stat_name]]

      central_cube <- generateLociCube(fpc, private$breaks, c(1,2,6))
      outer_cubes <- lapply(private$trio_classes, function(loci) {
        if (length(loci) == 0) return(NULL)
        generateLociCube(fpc, private$breaks, 4:5, loci)
      })
      
      c(central=central_cube, unlist(outer_cubes))
    },
    get_breaks = function() private$breaks
  )
)


calcBreaks <- function(values, props) {
  breaks <- unique(quantile(values, props, na.rm=TRUE))
  breaks[breaks == 0] <- ifelse(length(breaks)==1, 0.01, min(0.01, min(breaks[-1])/2))
  breaks
}


generateLociCube = function(stat, breaks, cols, rows) {
  stopifnot(ncol(stat) == length(breaks))
  stopifnot(all(cols <= ncol(stat)))
  
  # Select rows & remove loci with NaN statistics
  if (missing(rows)) rows <- rep(TRUE, nrow(stat))
  else rows <- 1:nrow(stat) %in% rows
  rows <- rows & apply(!is.nan(stat[ ,cols, drop=FALSE]), 1, all)
  stat <- stat[rows, cols, drop=FALSE]
  breaks <- breaks[cols]
  
  # Classify the loci accordingly to their statistics
  locus_class <- matrix(1, nrow(stat), ncol(stat))
  for (i in 1:ncol(stat)) {
    for (brk in breaks[[i]]) {
      locus_class[,i] <- locus_class[,i] + (stat[,i] > brk)
    }
  }
  
  # Count the classes and return as vector
  dims <- sapply(breaks, length) + 1
  factors <- cumprod(c(1, dims[-length(dims)]))
  classes_int <- apply(locus_class, 1, function(x) sum((x-1)*factors)+1)
  tabulate(classes_int, nbins = prod(dims))
}

#' A function that cassifies locus trios by the distance between the loci
#' 
#' Each pair of loci in a trio can be either 'near' together 
#' (low distance between) or more 'far' appart (somewhat larger distance), 
#' where the distances corresponding to 'near' and 'far' can be defined as
#' function arguments. The function returns trios for with the left to middle 
#' and middle to right comparisons is either 'both near', one near, one far' 
#' or 'both far'. Trio for which at least one distance is to short for 'near' 
#' or to large for 'far' are ignored.
#' 
#' @param llm The Locus Length matrix as produced by 
#'   \code{dm.getLociLengthMatrix}
#' @param near A vector of length two, giving the boundaries for the 'near' 
#'   class.
#' @param far A vector of length two, giving the boundaries for the 'far' 
#'   class.
#' @return A list with entries 'both_near', 'one_one' and 'both_far', which
#'   are vectors of the indexes of the loci that fall into the corresponding
#'   class.
#' @author Paul Staab
classifyTriosByDistance <- function(llm, near=c(5e3, 1e4), far=c(1e4, 2e4)) {
  stopifnot(is.matrix(llm))
  is_near <- llm[ ,c(2,4), drop = FALSE] >= near[1] & 
             llm[ ,c(2,4), drop = FALSE] < near[2]
  is_far <- llm[ ,c(2,4), drop = FALSE] >= far[1] & 
            llm[ ,c(2,4), drop = FALSE] < far[2]
  
  list(both_near=which(is_near[,1, drop = FALSE] & is_near[,2, drop = FALSE]),
       one_one=which((is_near[,1, drop = FALSE] & is_far[,2, drop = FALSE]) | 
                     (is_far[,1, drop = FALSE] & is_near[,2, drop = FALSE])),
       both_far=which(is_far[,1, drop = FALSE] & is_far[,2, drop = FALSE]))
}