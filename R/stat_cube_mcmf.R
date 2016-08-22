stat_cube_class_mcmf <- R6::R6Class("stat_cube_mcmf", inherit = stat_basic_class,
  lock_objects = FALSE, lock_class = TRUE,
  private = list(
    break_values = numeric(),
    bal_breaks = numeric(),
    poly_breaks = numeric(),
    mcmf_col = 1,
    calc_break_values = function(values, props) {
      "calculate the actual values for break_values"
      #browser()
      
      if(private$name == "mcmf") {
        if(private$mcmf_col == 2) props <- private$bal_breaks
        if(private$mcmf_col == 3) props <- private$poly_breaks
      }
      
      #if(private$mcmf_col != 2 || private$name != "mcmf") {
        break_values <- unique(quantile(values, props, na.rm = TRUE))
      #}else{
      #  break_values <- props 
      #}
      
      break_values[break_values == 0] <- ifelse(length(break_values) == 1, 
                                    0.01, min(0.01, min(break_values[-1]) / 2))
      
      if(private$name == "mcmf") private$mcmf_col <- private$mcmf_col + 1 #; browser()}
      break_values
    },
    generate_cube = function(stat, break_values, cols) {
      assert_that(is.matrix(stat))
      assert_that(is.list(break_values))
      assert_that(ncol(stat) == length(break_values))
      assert_that(all(cols <= ncol(stat)))
      #browser()
      # Remove rows that contain NAs or NaNs
      stat <- stat[apply(stat, 1, function(x) all(is.finite(x))), , #nolint
                   drop = FALSE]
      #browser()
      # Classify the loci accordingly to their statistics
      locus_class <- matrix(1, nrow(stat), ncol(stat))
      for (i in 1:ncol(stat)) {
        for (brk in break_values[[i]]) {
          locus_class[, i] <- locus_class[, i] + (stat[, i] > brk)
        }
      }
      #browser()
      # Count the classes and return as vector
      dims <- vapply(break_values, length, numeric(1)) + 1
      factors <- cumprod(c(1, dims[-length(dims)]))
      classes_int <- apply(locus_class, 1, 
                           function(x) sum((x - 1) * factors) + 1) #nolint
      #browser()
      tabulate(classes_int, nbins = prod(dims))
    }
  ),
  public = list(
    initialize = function(name, calc_func, break_values, bal_breaks, poly_breaks) {
      assert_that(is.function(calc_func))
      private$calculate_matrix <- calc_func
      #browser()
      super$initialize(name, function(data, opts) {
        stat <- private$calculate_matrix(data)
        assert_that(is.numeric(stat))
        if (!is.matrix(stat)) stat <- matrix(stat, ncol = 1)
        private$generate_cube(stat, opts$break_values, 1:ncol(stat))
      })
      
      assert_that(is.numeric(break_values))
      assert_that(length(break_values) > 0)
      if (any(break_values < 0 | break_values > 1)) {
        stop("probs greater then one")
      }
      private$break_values <- break_values
      private$bal_breaks <- bal_breaks
      private$poly_breaks <- poly_breaks
    },
    generate_data_opts = function(data) {
      data_matrix <- private$calculate_matrix(data)
      assert_that(is.numeric(data_matrix))
      if (!is.matrix(data_matrix)) data_matrix <- matrix(data_matrix, ncol = 1)
      list(break_values = lapply(seq_len(ncol(data_matrix)), function(i) {
        #browser()
        private$calc_break_values(data_matrix[, i], private$break_values)
      }))
    }
  )
)
