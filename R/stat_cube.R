#' @importFrom R6 R6Class
stat_cube_class <- R6Class("stat_cube", inherit = stat_basic_class,
  lock_objects = FALSE, lock_class = TRUE,
  private = list(
    breaks = numeric(),
    calc_breaks = function(values, props) {
      "calculate the actual values for breaks"
      breaks <- unique(quantile(values, props, na.rm = TRUE))
      breaks[breaks == 0] <- ifelse(length(breaks)==1, 
                                    0.01, min(0.01, min(breaks[-1])/2))
      breaks
    },
    generate_cube = function(stat, breaks, cols) {
      stopifnot(ncol(stat) == length(breaks))
      stopifnot(all(cols <= ncol(stat)))
      
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
  ),
  public = list(
    initialize = function(name, calc_func, breaks) {
      assert_that(is.function(calc_func))
      private$calculate_matrix <- calc_func
      
      super$initialize(name, function(data, opts) {
        stat <- private$calculate_matrix(data)                  
        private$generate_cube(stat, opts$breaks, 1:ncol(stat))
      })
      
      assert_that(is.numeric(breaks))
      assert_that(length(breaks) > 0)
      if (any(breaks < 0 | breaks > 1)) stop("probs greater then one")
      private$breaks <- breaks
    },
    generate_data_opts = function(data) {
      data_matrix <- private$calculate_matrix(data)
      list(breaks = lapply(1:ncol(data_matrix), function(i) {
        private$calc_breaks(data_matrix[, i], private$breaks)
      }))
    }
  )
)
