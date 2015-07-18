#' @importFrom R6 R6Class
jaatha_model_class <- R6Class("jaatha_model", 
  lock_objects = FALSE, lock_class = TRUE,
  private = list(
    par_ranges = NA,
    sum_stats = list(),
    add_statistic = function(stat) {
      private$sum_stats[[stat$get_name()]] <- stat
    },
    opts = list(),
    test = function(quiet = FALSE) {
      time <- system.time(
        a <- self$simulate(self$get_par_ranges()$get_middle(), 1)
      )["elapsed"]
      
      if (time > 30) warning("Each simulation takes about ", round(time),
                             "s, Jaatha might run for a long time.")
      if (!quiet) {
        if (time < 1) message("A simulation takes less than a second")
        else message("A simulation takes about ", round(time), "s")
      }
      
      invisible(NULL)
    }
  ),
  public = list(
    initialize = function(sim_func, par_ranges, sum_stats, ..., test) {
      assert_that(is.function(sim_func))
      private$sim_func <- sim_func
      private$par_ranges <- par_range_class$new(par_ranges)
      assert_that(is.list(sum_stats))
      lapply(sum_stats, private$add_statistic)
      private$sum_stats <- sum_stats
      private$opts <- list(...)
      if (test) private$test()
    },
    simulate = function(pars, seed) {
      # Simulate
      set.seed(seed)
      sim_pars <- private$par_ranges$denormalize(pars)
      sim_result <- private$sim_func(sim_pars)
      
      # Calculate Summary Statistics
      sim_stats <- lapply(private$sum_stats, function(sum_stat) {
        sum_stat$calculate(sim_result)
      })
      
      # Add the parameter values
      sim_stats$pars <- sim_pars
      sim_stats$pars_normal <- pars
      
      sim_stats
    },
    get_par_ranges = function() private$par_ranges,
    get_opts = function(name) private$opts[[name]]
  )
)


create_jaatha_model <- function(x, ..., test = TRUE) {
  UseMethod("create_jaatha_model")
}


create_jaatha_model.function <- function(x, par_ranges, sum_stats, ..., 
                                         test = TRUE) {
  jaatha_model_class$new(x, par_ranges, sum_stats, ..., test)
}


create_test_model <- function() {
  create_jaatha_model(function(x) rpois(10, x),
                      par_ranges = matrix(c(0.1, 0.1, 10, 10), 2, 2),
                      sum_stats = list(stat_identity))
}


create_jaatha_data <- function(...) {}

create_popgen_model <- function(...) {}
create_popgen_data <- function(...) {}