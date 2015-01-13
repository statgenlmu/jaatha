# -------------------------------------------------------------------------
# summary_statistics.R
#
# This file contains the default summary statistic functions for jaatha,
# which take sums over different areas of the jsfs.
# 
# Author:   Lisha Mathew & Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2013-09-04
# Licence:  GPLv3 or later
# -------------------------------------------------------------------------

# These are the default summary statistics for Jaatha
# 
# @param jsfs    The joint site frequency spectrum of two populations
# @return        A vector with sums over different areas of the JSFS
#
# @examples
# jsfs <- matrix(rpois(26*26,5),26,26)
# summarizeJSFS(jsfs = jsfs)
summarizeJSFS <- function(jsfs){
  if (is.list(jsfs)) jsfs <- jsfs$jsfs 
  n <- nrow(jsfs)
  m <- ncol(jsfs)
  c(sum(jsfs[1,2:3]),
    sum(jsfs[2:3,1]),
    sum(jsfs[1,4:(m-3)]),
    sum(jsfs[4:(n-3),1]),
    sum(jsfs[1,(m-2):(m-1)]),
    sum(jsfs[(n-2):(n-1),1]),
    sum(jsfs[2:3,2:3]),
    sum(jsfs[2:3,4:(m-3)]),
    sum(jsfs[4:(n-3),2:3]),
    sum(jsfs[(n-2):(n-1),4:(m-3)]),
    sum(jsfs[4:(n-3),(m-2):(m-1)]),
    sum(jsfs[2:3,(m-2):(m-1)]),
    sum(jsfs[(n-2):(n-1),2:3]),
    sum(jsfs[4:(n-3),4:(m-3)]),
    sum(jsfs[(n-2):(n-1),(m-2):(m-1)]),
    jsfs[1,m],
    jsfs[n,1],
    sum(jsfs[n,2:3]),
    sum(jsfs[2:3,m]),
    sum(jsfs[n,4:(m-3)]),
    sum(jsfs[4:(n-3),m]),
    sum(jsfs[n,(m-2):(m-1)]),
    sum(jsfs[(n-2):(n-1),m]) )
}


summarizeJsfsBorder <- function(sim_data) {
  jsfs <- sim_data$jsfs
  n <- nrow(jsfs)
  m <- ncol(jsfs)
  c(sum(jsfs[1,2:3]),
    sum(jsfs[1,4:(m-3)]),
    sum(jsfs[1,(m-2):m]),
    sum(jsfs[2:3,1]),
    sum(jsfs[4:(n-3),1]),
    sum(jsfs[(n-2):(n),1]),
    sum(jsfs[n,2:3]),
    sum(jsfs[n,4:(m-3)]),
    sum(jsfs[n,(m-2):(m-1)]),
    sum(jsfs[2:3,m]),
    sum(jsfs[4:(n-3),m]),
    sum(jsfs[(n-2):(n-1),m]) )
}


summarizeFoldedJSFS <- function(sim_data) {
  jsfs <- sim_data$jsfs 
  n <- nrow(jsfs)
  m <- ncol(jsfs)

  sumstats <- 
    c(sum(jsfs[1, 2:3], jsfs[n, (m-2)]:(m-1)),
      sum(jsfs[1, 4:(m-3)], jsfs[n, 4:(m-3)]),
      sum(jsfs[1, (m-2):(m-1)], jsfs[n, 2:3]),
      sum(jsfs[1, m], jsfs[n, 1]),
      sum(jsfs[2:3,1], jsfs[(n-2):(n-1),m]),
      sum(jsfs[2:3,2:3], jsfs[(n-2):(n-1),(m-2):(m-1)]),
      sum(jsfs[2:3,4:(m-3)], jsfs[(n-2):(n-1),4:(m-3)]),
      sum(jsfs[2:3,(m-2):(m-1)], jsfs[(n-2):(n-1),2:3]),
      sum(jsfs[n,2:3], jsfs[(n-2):(n-1),1]),
    
      sum(jsfs[4:(n-3),1]),
      sum(jsfs[4:(n-3),2:3]),
      sum(jsfs[4:(n-3),(m-2):(m-1)]),
      sum(jsfs[4:(n-3),4:(m-3)]),
      sum(jsfs[4:(n-3),m])              
      )

  return(sumstats)
}

# Binning
Stat_JSFS <- R6Class('Stat_PoiInd', 
  inherit = Stat_Base,
  public = list(
    initialize = function(seg_sites, dm) {
      private$data = self$transform(list(jsfs=calcJsfs(seg_sites, dm.getSampleSize(dm))))
    },
    transform = summarizeJSFS
  )
)

# Binning + Folded JSFS
Stat_JSFS_folded <- R6Class('Stat_PoiInd', 
  inherit = Stat_Base,
  public = list(
    initialize = function(seg_sites, dm) {
       private$data = self$transform(list(jsfs=calcJsfs(seg_sites, 
                                                        dm.getSampleSize(dm))))
    },
  transform = summarizeFoldedJSFS
  )
)

# Smoothing
Stat_JSFS_smooth <- R6Class('Stat_PoiSmooth',
  inherit = Stat_Base,
  private = list(
      model = NA,
      border_mask = NA,
      border_stat = NA
  ),
  public = list(
    initialize = function(seg_sites, dm) {
      sample_size <- dm.getSampleSize(dm)
      jsfs = calcJsfs(seg_sites, sample_size)
      private$model <- paste0("( X1 + I(X1^2) + X2 + I(X2^2) + log(X1) + log(",
                              sample.size[1]+2,
                              "-X1) + log(X2) + log(",
                              sample.size[2]+2,
                              "-X2) )^2")
      
      private$border_mask <- jsfs.value
      private$border_mask[, ] <- 0
      private$border_mask[c(1, nrow(jsfs.value)), ] <- 1
      private$border_mask[ ,c(1, ncol(jsfs.value))] <- 1
      private$border_mask <- as.logical(border.mask)
      
      private$data = self$transform(list(jsfs=jsfs))
    },
    transformation = function(sim_data) sim_data$jsfs,
    get_border_mask = function() private$border_mask
  )
)

  