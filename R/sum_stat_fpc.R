#' @importFrom R6 R6Class
Stat_FPC <- R6Class('Stat_FPC', inherit = Stat_PoiInd,
  public = list(
    initialize = function(seg_sites, dm, population, group = 0,
                          break_probs = c(.2, .5)) {

      private$individuals = getIndOfPop(dm, population)
      private$llm = dm.getLociLengthMatrix(dm, group)
      if (any(break_probs > 1)) stop('probs greater then one')
      
      fpc.percent <- calcPercentFpcViolation(seg_sites, 
                                             private$individuals, 
                                             private$llm)
      
      private$breaks = list(near=calcBreaks(fpc.percent[, 1], break_probs),
                            far=calcBreaks(fpc.percent[, 2], break_probs),
                            mut=calcBreaks(fpc.percent[, 6], break_probs))
      
      # Calculate observed values
      private$data = self$transform(list(seg.sites = seg_sites))
      if (group > 0) private$seg_sites_name = paste0('seg.sites.', group)
    },
    generate = function(seg_sites, breaks = private$breaks) {
      stopifnot(!is.null(seg_sites))
      percent <- calcPercentFpcViolation(seg_sites,
                                         private$individuals, 
                                         private$llm)
      
      percent <- percent[!(is.nan(percent[,1]) | is.nan(percent[,2]) ), , 
                         drop = FALSE ]
      
      # Classify the loci accordingly
      locus_class <- matrix(1, nrow(percent), 3)
      
      for (brk in breaks$near) {
        locus_class[,1] <- locus_class[,1] + (percent[,1] > brk)
      }
      for (brk in breaks$far) {
        locus_class[,2] <- locus_class[,2] + (percent[,2] > brk)
      }
      for (brk in breaks$mut) {
        locus_class[,3] <- locus_class[,3] + (percent[,6] > brk)
      }
      
      # Count the occurance of each class in a matrix
      dims <- c(length(private$breaks$near) + 1, 
                length(private$breaks$far) + 1, 
                length(private$breaks$mut) + 1)
      
      countClasses(t(locus_class), dims)
    },
    transform = function(sim_data) {
      as.vector(self$generate(sim_data[[private$seg_sites_name]]))
    },
    get_breaks = function() private$breaks
  ),
  private = list(
    seg_sites_name = 'seg.sites',
    individuals = NA,
    llm = NA,
    breaks = NA)
)

calcBreaks <- function(values, props) {
  breaks <- unique(quantile(values, props, na.rm=TRUE))
  breaks[breaks == 0] <- ifelse(length(breaks)==1, 0.01, min(0.01, min(breaks[-1])/2))
  breaks
}

countClasses <- function(classes, dimension) {
  stopifnot(nrow(classes) == length(dimension))
  stat <- array(0, dim = dimension)
  if(ncol(classes) == 0) return(stat)
  
  # Replace NA's with the last value
  for (r in 1:nrow(classes)) {
    classes[r, is.na(classes[r,])] <- dimension[r]
  }
  
  # Count occurences of classes
  if (nrow(classes) == 2) {
    for (r in 1:ncol(classes)) {
      stat[classes[1, r], classes[2, r]] <- 
        stat[classes[1, r], classes[2, r]] + 1
    }
  } else if (nrow(classes) == 3) {
    for (r in 1:ncol(classes)) {
      stat[classes[1, r], classes[2, r], classes[3, r]] <- 
        stat[classes[1, r], classes[2, r], classes[3, r]] + 1
    }
  } else stop("This function can only create 2 and 3-dim. arrays")
  
  stat
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
#' @param dm The demographic model from which the locus trios are taken.
#' @param group The group of loci that is classified
#' @param near A vector of length two, giving the boundaries for the 'near' 
#'   class.
#' @param far A vector of length two, giving the boundaries for the 'far' 
#'   class.
#' @return A list with entries 'both_near', 'one_one' and 'both_far', which
#'   are vectors of the indexes of the loci that fall into the corresponding
#'   class.
#' @author Paul Staab
classifyTriosByDistance <- function(dm, group = 0, 
                                    near=c(5e3, 1e4), far=c(1e4, 2e4)) {
  
  llm <- dm.getLociLengthMatrix(dm, group)
  is_near <- llm[ ,c(2,4)] >= near[1] &  llm[ ,c(2,4)] < near[2]
  is_far <- llm[ ,c(2,4)] >= far[1] &  llm[ ,c(2,4)] < far[2]
  
  list(both_near=which(is_near[,1] & is_near[,2]),
       one_one=which((is_near[,1] & is_far[,2]) | (is_far[,1] & is_near[,2])),
       both_far=which(is_far[,1] & is_far[,2]))
}