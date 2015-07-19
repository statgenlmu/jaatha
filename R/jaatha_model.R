#' @importFrom R6 R6Class
jaatha_model_class <- R6Class("jaatha_model", 
  lock_objects = FALSE, lock_class = TRUE,
  private = list(
    par_ranges = NA,
    sum_stats = list(),
    add_statistic = function(stat) {
      name <- stat$get_name()
      if (!is.null(private$sum_stats[[name]])) {
        stop("There already is a summary statistic with name '", name, 
             "' in the model")
      }
      private$sum_stats[[name]] <- stat
    }
  ),
  public = list(
    initialize = function(sim_func, par_ranges, sum_stats, test) {
      assert_that(is.function(sim_func))
      private$sim_func <- sim_func
      private$par_ranges <- par_ranges_class$new(par_ranges)
      assert_that(is.list(sum_stats))
      lapply(sum_stats, private$add_statistic)
      if (test) self$test()
    },
    simulate = function(pars, seed, data) {
      assert_that(is_jaatha_data(data))
      
      # Simulate
      set.seed(seed)
      sim_pars <- private$par_ranges$denormalize(pars)
      sim_result <- private$sim_func(sim_pars)
      
      # Calculate Summary Statistics
      sum_stats <- lapply(private$sum_stats, function(sum_stat) {
        sum_stat$calculate(sim_result, data$get_opts(sum_stat))
      })
      
      # Add the parameter values
      sum_stats$pars <- sim_pars
      sum_stats$pars_normal <- pars
      
      sum_stats
    },
    get_par_ranges = function() private$par_ranges,
    get_sum_stats = function() private$sum_stats,
    test = function(quiet = FALSE) {
      time <- system.time(
        sim_data <- private$sim_func(private$par_ranges$get_middle())
      )["elapsed"]
      
      if (!quiet) {
        if (time > 30) warning("Each simulation takes about ", round(time),
                               "s, Jaatha might run for a long time.")
        
        if (time < 1) message("A simulation takes less than a second")
        else message("A simulation takes about ", round(time), "s")
      }
      
      invisible(sim_data)
    }
  )
)


create_jaatha_model <- function(x, ..., test = TRUE) {
  UseMethod("create_jaatha_model")
}


create_jaatha_model.function <- function(x, par_ranges, sum_stats,
                                         test = TRUE) {
  jaatha_model_class$new(x, par_ranges, sum_stats, test = test)
}


is_jaatha_model <- function(x) inherits(x, "jaatha_model")


create_test_model <- function() {
  create_jaatha_model(function(x) rpois(10, x),
                      par_ranges = matrix(c(0.1, 0.1, 10, 10), 2, 2),
                      sum_stats = list(stat_identity(), stat_sum()),
                      test = FALSE)
}
