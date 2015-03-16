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
    sum(jsfs[4:(n-3),m]))
}

# Binning
#' @importFrom coalsimr calc_jsfs get_population_indiviuals
Stat_JSFS <- R6Class('Stat_JSFS', 
  inherit = Stat_PoiInd,
  public = list(
    initialize = function(seg_sites, model, group=0) {
      assert_that(is.list(seg_sites))
      name <- getStatName('jsfs', group)
      fake_sim_data <- list()
      fake_sim_data[[name]] <- calc_jsfs(seg_sites, 
                                         get_population_indiviuals(model, 1),
                                         get_population_indiviuals(model, 2))
      super$initialize(fake_sim_data, name)
    },
    transform = function(sim_data) {
      summarizeJSFS(sim_data[[private$name]])
    }
  )
)

Stat_JSFS_border <- R6Class('Stat_JSFS_border', 
  inherit = Stat_JSFS,
  public = list(
    transform = function(sim_data) summarizeFoldedJSFS(sim_data[[private$name]]),
    get_name = function() paste0('border_', private$name)
  )
)

# Smoothing
#' @importFrom coalsimr calc_jsfs get_population_indiviuals get_sample_size
Stat_JSFS_smooth <- R6Class('Stat_JSFS_smooth',
  inherit = Stat_PoiSmooth,
  private = list(
      model = NA,
      rows = NA,
      cols = NA
  ),
  public = list(
    initialize = function(seg_sites, dm, group=0) {
      sample_size <- get_sample_size(dm)
      model <- paste0("( X1 + I(X1^2) + X2 + I(X2^2) + log(X1) + log(",
                      sample_size[1]+1, "-X1) + log(X2) + log(",
                      sample_size[2]+1, "-X2) )^2")
      
      private$rows <- 2:sample_size[1]
      private$cols <- 2:sample_size[2]
      
      name <- getStatName('jsfs', group)
      fake_sim_data <- list()
      fake_sim_data[[name]] <- calc_jsfs(seg_sites, 
                                         get_population_indiviuals(dm, 1),
                                         get_population_indiviuals(dm, 2))
      super$initialize(fake_sim_data, name, model)
      
      if (all(self$get_data()$sum.stat == 0)) 
        stop("Inner JSFS used for smoothing only consists of 0s")
    },
    transform = function(sim_data) {
      private$to_data_frame(sim_data[[private$name]][private$rows, private$cols])
    }
  )
)


