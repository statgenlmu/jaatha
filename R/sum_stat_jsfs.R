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
summarizeJSFS <- function(jsfs) {
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


summarizeJsfsBorder <- function(jsfs) {
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


summarizeFoldedJSFS <- function(jsfs) {
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
    initialize = function(seg_sites, dm, group=0) {
      private$data = self$transform(list(jsfs=calcJsfs(seg_sites, dm.getSampleSize(dm))))
      if (group > 0) private$jsfs_name = paste0('jsfs.', group)
    },
    transform = function(sim_data) summarizeJSFS(sim_data[[private$jsfs_name]])
  ),
  private = list(
    jsfs_name = 'jsfs'
  )
)

# Binning + Folded JSFS
Stat_JSFS_folded <- R6Class('Stat_PoiInd', 
  inherit = Stat_JSFS,
  public = list(
    transform = function(sim_data) summarizeFoldedJSFS(sim_data[[private$jsfs_name]])
  )
)

Stat_JSFS_border <- R6Class('Stat_PoiInd', 
  inherit = Stat_JSFS,
  public = list(
    transform = function(sim_data) summarizeFoldedJSFS(sim_data[[private$jsfs_name]])
  )
)

# Smoothing
Stat_JSFS_smooth <- R6Class('Stat_PoiSmooth',
  inherit = Stat_PoiSmooth,
  private = list(
      model = NA,
      rows = NA,
      cols = NA,
      jsfs_name = 'jsfs'
  ),
  public = list(
    initialize = function(seg_sites, dm, group=0) {
      sample_size <- dm.getSampleSize(dm)
      jsfs = calcJsfs(seg_sites, sample_size)
      private$model <- paste0("( X1 + I(X1^2) + X2 + I(X2^2) + log(X1) + log(",
                              sample_size[1]+1,
                              "-X1) + log(X2) + log(",
                              sample_size[2]+1,
                              "-X2) )^2")
      
      private$rows <- 2:sample_size[1]
      private$cols <- 2:sample_size[2]
      
      self$set_data(list(jsfs=jsfs))
      if (all(self$get_data()$sum.stat == 0)) 
        stop("Inner JSFS used for smoothing only consists of 0s")
      
      if (group > 0) private$jsfs_name = paste0('jsfs.', group)
    },
    transform = function(sim_data) {
      private$to_data_frame(sim_data[[private$jsfs_name]][private$rows, 
                                                          private$cols])
    }
  )
)


