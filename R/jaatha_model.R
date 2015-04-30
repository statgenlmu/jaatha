#' @importFrom R6 R6Class
jaatha_model_class <- R6Class("jaatha_model", lock = FALSE, 
  private = list(
    par_ranges = NA,
    sum_stats = list()
  ),
  public = list(
    initialize = function(sim_func, par_ranges, sum_stats) {
      private$sim_func = sim_func
      private$par_ranges = par_ranges
      private$sum_stats = sum_stats
    },
    simulate = function(pars, seed) {
      # Simulate
      set.seed(seed)
      sim_pars <- self$denormalize_pars(pars)
      sim_result <- private$sim_func(sim_pars)
      
      # Calculate Summary Statistics
      sim_stats <- lapply(private$sum_stats, function(sum_stat) {
        sum_stat$transform(sim_result)
      })
      
      # Add the parameter values
      sim_stats$pars <- sim_pars
      sim_stats$pars_normal <- pars
      
      sim_stats
    }
  )
)


create_jaatha_model <- function(sim_func, par_ranges, sum_stats, test = TRUE) {
  model <- jaatha_model_class$new(sim_func, par_ranges, sum_stats)
  
}





create_jaatha_data <- function(...) {}

create_popgen_model <- function(...) {}
create_popgen_data <- function(...) {}