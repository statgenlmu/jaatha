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
Stat_JSFS <- R6Class('Stat_JSFS', 
  inherit = Stat_PoiInd,
  public = list(
    initialize = function(seg_sites, model, stat) {
      name <- stat$get_name()
      jsfs <- stat$calculate(seg_sites, NULL, model)
      fake_sim_data <- list()
      fake_sim_data[[name]] <- jsfs
      super$initialize(fake_sim_data, stat$get_name())
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
Stat_JSFS_smooth <- R6Class('Stat_JSFS_smooth',
  inherit = Stat_PoiSmooth,
  private = list(
      model = NA,
      rows = NA,
      cols = NA
  ),
  public = list(
    initialize = function(seg_sites, dm, stat) {
      name <- stat$get_name()
      jsfs <- stat$calculate(seg_sites, NULL, dm)
      fake_sim_data <- list()
      fake_sim_data[[name]] <- jsfs
      
      model <- paste0("( X1 + I(X1^2) + X2 + I(X2^2) + log(X1) + log(",
                      nrow(jsfs), "-X1) + log(X2) + log(",
                      ncol(jsfs), "-X2) )^2")
      
      private$rows <- 2:(nrow(jsfs)-1)
      private$cols <- 2:(ncol(jsfs)-1)
      
      super$initialize(fake_sim_data, name, model)
      
      if (all(self$get_data()$sum.stat == 0)) 
        stop("Inner JSFS used for smoothing only consists of 0s")
    },
    transform = function(sim_data) {
      private$to_data_frame(sim_data[[private$name]][private$rows, private$cols])
    }
  )
)
