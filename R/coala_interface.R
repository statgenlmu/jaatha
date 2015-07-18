#' @importFrom coala get_parameter_table get_summary_statistics 
#' @importFrom coala get_locus_number scale_model
create_jaatha_model.coalmodel <- function(x, ..., test = TRUE) {
  sim_func <- function(pars) {
    simulate(self$get_opts("coalmodel"), pars = pars)
  }
  
  # Get parameter ranges
  par_ranges <- as.matrix(get_parameter_table(model)[,-1])
  rownames(par_ranges) <- get_parameter_table(model)$name

  # Get summary statisics
  sum_stats <- list()
  
  create_jaatha_model.function(sim_func, par_ranges, sum_stats, 
                               coalmodel = x, ..., test = test)
}