#' Use a coala model in Jaatha
#' 
#' This creates a Jaatha model from a coala model. Simulation for this model
#' model are conducted via the \code{simulate} function for the coala model.
#' The parameters that are
#' estimated must be specified via \code{\link[coala]{par_range}} and the
#' model must not have any other named parameters. Summary statistics present 
#' in the coala model are used in Jaatha.
#' 
#' @param x The coala model
#' @param jsfs_summary The way the Joint Site Frquency Spectrum (JSFS) 
#'   is further summarized. Can be \code{sums} (default), \code{none} or 
#'   \code{"smoothing"}. For \code{sums}, 23 different areas of the JSFS
#'   are summed up, and the sums are used as indepented Poission statistcs, 
#'   for \code{none}, all entries are used as indepented Possion statistics.
#'   The value \code{smooth} is experimental so far and should not be used.
#'   This option has no effect if the JSFS used as summary statistic in the
#'   coala model.
#' @param ihs_breaks Quantiles of the real data that will be used as breaks
#'   for binning the iHS statistic if present in the model.
#' @param four_gamete_breaks Quantiles of the real data that will be used as 
#'   breaks for binning the Four Gamete test based statistic if present in the 
#'   model.
#' @param mcmf_breaks Quantiles of the real data that will be used as breaks
#'   for binning the MCMF statistic if present in the model.
#' @inheritParams create_jaatha_model
#' @export
create_jaatha_model.coalmodel <- function(x, 
                                          jsfs_summary = c("sums",
                                                           "none",
                                                           "smooth"),
                                          ihs_breaks = c(.5, .7, .9),
                                          four_gamete_breaks = c(.2, .5),
                                          mcmf_breaks = c(.5, .7, .9),
                                          ...,
                                          scaling_factor = 1,
                                          test = TRUE) {
  
  if (!requireNamespace("coala", quietly = TRUE)) {
    stop("Please install coala to use this function", call. = FALSE)
  }
  
  if (length(jsfs_summary) > 1) jsfs_summary <- jsfs_summary[1]
  
  sim_func <- function(pars, opts = NULL) simulate(x, pars = pars)
  
  # create parameter ranges
  par_table <- coala::get_parameter_table(x)
  par_ranges <- as.matrix(par_table[,-1])
  rownames(par_ranges) <- par_table$name

  # create summary statisics
  sum_stats <- convert_coala_sumstats(x, jsfs_summary, ihs_breaks,
                                      four_gamete_breaks, mcmf_breaks)
  
  create_jaatha_model.function(sim_func, par_ranges, sum_stats, 
                               test = test)
}


convert_coala_sumstats <- function(coala_model, jsfs_summary = "sums",
                                   ihs_breaks, four_gamete_breaks, 
                                   mcmf_breaks) {
  
  if (!requireNamespace("coala", quietly = TRUE)) {
    stop("Please install coala to use this function", call. = FALSE)
  }
  
  assert_that(is.string(jsfs_summary))
  
  lapply(coala::get_summary_statistics(coala_model), function(stat) {
    name <- stat$get_name()
    
    # --- JSFS Summary Statistic ------------------------------------
    if (inherits(stat, "stat_jsfs")) {
      if (jsfs_summary == "sums") {
        return(create_jaatha_stat(name, function(x, opts) {
          sum_jsfs(x[[name]])
        }))
      } else if (jsfs_summary == "none") {
        return(create_jaatha_stat(name, function(x, opts) {
          as.vector(x[[name]])[-c(1, prod(dim(x[[name]])))]
        }))
      } else if (jsfs_summary == "smooth") {
        stop("Smoothing is not suppored right now")
      }
    }
    
    # --- JSFS Summary Statistic ------------------------------------
    if (inherits(stat, "stat_sfs")) {
      return(create_jaatha_stat(name, function(x, opts) x[[name]]))
    }
    
    # --- Four Gamete Summary Statistic -----------------------------
    if (inherits(stat, "stat_four_gamete")) {
      return(create_jaatha_stat(name, function(x, opts) {
        x[[name]][ , c(1, 2, 6), drop = FALSE]
      }, poisson = FALSE, breaks = four_gamete_breaks))
    }
    
    # --- iHH Summary Statistic -------------------------------------
    if (inherits(stat, "stat_ihh")) {
      return(create_jaatha_stat(name, function(x, opts) {
        vapply(x[[name]], function(x) max(x[ , 3]), numeric(1))
      }, poisson = FALSE, breaks = ihs_breaks))
    }
    
    # --- OmegaPrime Summary Statistic ----------------------------------
    if (inherits(stat, "stat_omega_prime") || inherits(stat, "stat_mcmf")) {
      return(create_jaatha_stat(name, function(x, opts) x[[name]],
                                poisson = FALSE, 
                                breaks = mcmf_breaks))
    }
    
    warning("Summary statistic '", name, "' is not supported. Ignoring it.")
    NULL
  })
}


sum_jsfs <- function(jsfs) {
  n <- nrow(jsfs)
  m <- ncol(jsfs)
  c(sum(jsfs[1, 2:3]),
    sum(jsfs[2:3, 1]),
    sum(jsfs[1, 4:(m - 3)]),
    sum(jsfs[4:(n - 3), 1]),
    sum(jsfs[1, (m - 2):(m - 1)]),
    sum(jsfs[(n - 2):(n - 1), 1]),
    sum(jsfs[2:3, 2:3]),
    sum(jsfs[2:3, 4:(m - 3)]),
    sum(jsfs[4:(n - 3), 2:3]),
    sum(jsfs[(n - 2):(n - 1), 4:(m - 3)]),
    sum(jsfs[4:(n - 3), (m - 2):(m - 1)]),
    sum(jsfs[2:3, (m - 2):(m - 1)]),
    sum(jsfs[(n - 2):(n - 1), 2:3]),
    sum(jsfs[4:(n - 3), 4:(m - 3)]),
    sum(jsfs[(n - 2):(n - 1), (m - 2):(m - 1)]),
    jsfs[1, m],
    jsfs[n, 1],
    sum(jsfs[n, 2:3]),
    sum(jsfs[2:3, m]),
    sum(jsfs[n, 4:(m - 3)]),
    sum(jsfs[4:(n - 3), m]),
    sum(jsfs[n, (m - 2):(m - 1)]),
    sum(jsfs[(n - 2):(n - 1), m]) )
}
