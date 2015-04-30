#' @importFrom R6 R6Class
jaatha_model_class <- R6Class("jaatha_model", lock = FALSE, 
  private = list(
    par_ranges = NA,
    sum_stats = list(),
    add_statistic = function(stat) {
      private$sum_stats[[stat$get_name()]] <- stat
    }
  ),
  public = list(
    initialize = function(sim_func, par_ranges, sum_stats) {
      assert_that(is.function(sim_func))
      private$sim_func = sim_func
      private$par_ranges = par_range_class$new(par_ranges)
      assert_that(is.list(sum_stats))
      lapply(sum_stats, private$add_statistic)
      private$sum_stats = sum_stats
    },
    simulate = function(pars, seed) {
      # Simulate
      set.seed(seed)
      sim_pars <- private$par_ranges$denormalize(pars)
      sim_result <- private$sim_func(sim_pars)
      
      # Calculate Summary Statistics
      sim_stats <- lapply(private$sum_stats, function(sum_stat) {
        sum_stat$transform(sim_result)
      })
      
      # Add the parameter values
      sim_stats$pars <- sim_pars
      sim_stats$pars_normal <- pars
      
      sim_stats
    },
    get_par_range = function() private$par_ranges
  )
)


create_jaatha_model <- function(sim_func, par_ranges, sum_stats, test = TRUE) {
  jaatha_model_class$new(sim_func, par_ranges, sum_stats)
}


create_test_model <- function() {
  sim_func <- function(x, jaatha) rpois(20, x)
  csi.obs <- csi.sim.func(c(3,5))
  sum_stat <- R6::R6Class("Stat_PoiInd", inherit = jaatha:::Stat_Base, 
                              private = list(mask=rep(c(TRUE,FALSE), 10)),
                              public = list(transform = function(data) {
                                c(sum(data[private$mask]), sum(data[!private$mask]))
                              })
  )$new(csi.obs, "csi")
  par_ranges <- matrix(c(0.1, 0.1, 10, 10), 2, 2)
  rownames(par_ranges) <- c("x", "y")
  create_jaatha_model(sim_func, par_ranges, list(sum_stat))
}



create_jaatha_data <- function(...) {}

create_popgen_model <- function(...) {}
create_popgen_data <- function(...) {}