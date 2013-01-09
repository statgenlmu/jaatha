# -------------------------------------------------------------------------
# summary_statistics.R
#
# This file contains the default summary statistic functions for jaatha,
# which take sums over different areas of the jsfs.
# 
# Author:   Lisha Naduvilezhath & Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-09-14
# Licence:  GPLv3 or later
# -------------------------------------------------------------------------

#' These are the default summary statistics for Jaatha
#' 
#' @param jsfs    The joint site frequency spectrum of two populations
#' @param dm      The corresponding demographic model. Not needed at the moment.
#' @return        A vector with sums over different areas of the JSFS
#' @export
#'
#' @examples
#' jsfs <- matrix(rpois(26*26,5),26,26)
#' Jaatha.defaultSumStats(jsfs = jsfs)
Jaatha.defaultSumStats <- function(dm = NULL, jsfs) {
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


Jaatha.defaultFoldedSumStats <- function(dm = NULL, jsfs) {
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
      sum(jsfs[4:(n-3),m])              )

  return(sumstats)
}
