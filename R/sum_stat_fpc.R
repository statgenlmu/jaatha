#' @importFrom R6 R6Class
Stat_FPC <- R6Class('Stat_PoiInd', inherit = Stat_Base,
  public = list(
    initialize = function(seg_sites, dm, population, group = 0,
                          break_props = c(.2, .5)) {
      private$group = group
      private$individuals = getIndOfPop(dm, population)
      private$llm = dm.getLociLengthMatrix(dm, group)
      if (any(break_props > 1)) stop('props greater then one')
      
      fpc.percent <- calcPercentFpcViolation(seg_sites, 
                                             private$individuals, 
                                             private$llm)
      
      private$breaks = list(near=calcBreaks(fpc.percent[, 1], break_props),
                            far=calcBreaks(fpc.percent[, 2], break_props),
                            mut=calcBreaks(fpc.percent[, 6], break_props))
      
      # Calculate observed values
      private$data = self$transform(list(seg.sites = seg_sites))
    },
    generate = function(seg_sites, breaks = private$breaks) {
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
      as.vector(self$generate(sim_data$seg.sites))
    },
    get_breaks = function() private$breaks
  ),
  private = list(
    group = NA,
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