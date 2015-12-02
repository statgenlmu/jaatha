#' Parametric Bootstrapping of Jaatha Estimates
#' 
#' This function is a helper function for using the \code{\link[boot]{boot}}
#' function to bootstrap Jaatha estimates. Each bootstap replication requires
#' a complete jaatha estimation on data simulated with the original parameter
#' estimates. Therefore, bootstrapping is normally computationally demanding and
#' should be executed on a computing cluster. The jaatha analyses are
#' resticted to a single CPU cores, so that as many replicas as possible can 
#' be executed in parallel using the corresponding options of 
#' \code{\link[boot]{boot}}.
#' 
#' @param model The jaatha model
#' @param data The jaatha data
#' @param results The results of an \code{\link{jaatha}} analysis performed with 
#'   the same \code{model} and \code{data} as passed to this function. 
#' @param R The number of bootstrapping replicates that are performed.
#' @param ... Additional arguments that are passed on \code{\link[boot]{boot}}.
#'   It is highly recommended to use its \code{parallel} and \code{ncpus} 
#'   options to parallelize the bootstrap replicates.
#' @return The result of \code{\link[boot]{boot}}. This object can be used to
#'   estimate standard errors or confidence intervals of the estimates using
#'   the functions available in package \pkg{boot}.
#' 
#' @importFrom utils capture.output
#' @export
boot_jaatha <- function(model, data, results, R, ...) {
  require_package("boot")
  if (R.Version()$major == 3 && R.Version()$minor < 2.2) {
    stop("This function requires at least R Version 3.2.2")
  }
  
  args <- results$args
  sim_func <- model$get_sim_func()
  
  log_folder <- tempfile("logs_")
  message("Logfiles in folder: ", tempdir())
  message("This might take a while...")
  dir.create(log_folder)
    
  jaatha_stat <- function(data) {
    capture.output({
      results <- jaatha(model, data,
                        repetitions = args$repetition,
                        sim = args$sim,
                        max_steps = args$max_steps,
                        init_method = args$init_method,
                        cores = 1)
    }, 
    file = tempfile(paste0("boot_log_", Sys.getpid(), "_")), 
    type = "message")
    
    results$estimate
  }
  
  simulate_data <- function(data, param) {
    sim_data <- sim_func(param)
    create_jaatha_data.default(sim_data, model)
  }
  
  boot::boot(data, jaatha_stat, R = R, 
             sim = "parametric", 
             ran.gen = simulate_data, 
             mle = results$estimate, ...)
}
