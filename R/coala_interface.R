#' @importFrom coala get_parameter_table get_summary_statistics 
#' @importFrom coala get_locus_number scale_model
create_jaatha_model.coalmodel <- function(x, ..., test = TRUE) {
  sim_func <- function(pars) {
    simulate(x, pars = pars)
  }
  
  # create parameter ranges
  par_table <- get_parameter_table(x)
  par_ranges <- as.matrix(par_table[,-1])
  rownames(par_ranges) <- par_table$name

  # create summary statisics
  sum_stats <- convert_coala_sumstats(x)
  
  create_jaatha_model.function(sim_func, par_ranges, sum_stats, 
                               ..., test = test)
}


#' @importFrom coala get_summary_statistics 
convert_coala_sumstats <- function(coala_model) {
  coala_sumstats <- get_summary_statistics(coala_model)
  sumstats <- list()
  
  for (stat in coala_sumstats) {
    name <- stat$get_name()
    
    # --- JSFS Summary Statistic ------------------------------------
    if (inherits(stat, "stat_jsfs")) {
      if (!smoothing) {
        sumstats[[name]] <- Stat_JSFS$new(seg_sites, model, stat)
      } else {
        sumstats[[name]] <- 
          Stat_JSFS_smooth$new(seg_sites, model, stat)
        sumstats[[paste0("border_", name)]] <- 
          Stat_JSFS_border$new(seg_sites, model, stat)
      }
    }
    
    # --- JSFS Summary Statistic ------------------------------------
    else if (inherits(stat, "stat_sfs")) {
      sumstats[[name]] <- create_jaatha_stat(name, function(x) x[[name]])
    }
    
    # --- Four Gamete Summary Statistic -----------------------------
    else if (inherits(stat, "stat_four_gamete")) {
      sumstats[[name]] <- Stat_FPC$new(seg_sites, model, stat)
    }
    
    # --- iHH Summary Statistic -------------------------------------
    else if (inherits(stat, "stat_ihh")) {
      sumstats[[name]] <- create_jaatha_stat(name, function(x) {
        do.call(rbind, lapply(x[[name]], function(x) max(x[ , 3])))
      }, poisson = FALSE, breaks = c(.25, .5, .75, .95))
    }
    
    # --- OmegaPrime Summary Statistic ----------------------------------
    else if (inherits(stat, "stat_omega_prime")) {
      sumstats[[name]] <- create_jaatha_stat(name, function(x) {
        matrix(x[[name]], ncol = 1)
      }, poisson = FALSE, breaks = c(.5, .75, .95))
    }
    
    else {
      warning("Summary statistic '", name, "' is not supported. Ignoring it.")
    }
  }
  
  sumstats
}
