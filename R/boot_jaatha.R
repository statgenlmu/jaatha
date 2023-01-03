#' Parametric Bootstrapping of Jaatha Estimates
#' 
#' This function is a helper function for using the \code{\link[boot]{boot}}
#' function to bootstrap Jaatha estimates. Each bootstap replication requires
#' a complete jaatha estimation on data simulated with the original parameter
#' estimates. Therefore, bootstrapping is normally computationally demanding and
#' should be executed on a computing cluster.
#' 
#' @param results The results of an \code{\link{jaatha}} analysis.
#' @param R The number of bootstrapping replicates that are performed.
#' @param cores_per_run The number of cores that are used for each replicate.
#'        This corresponds to the \code{cores} option of \code{\link{jaatha}}.
#'        Different replicates can be executed in parallel using the 
#'        \code{parallel}, \code{ncpus} and \code{cl} options of 
#'        \code{\link[boot]{boot}}.  The total number of CPU cores
#'        used is \code{ncpus} * \code{cores_per_run}.
#' @param verbose If TRUE (default), each bootstrap estimation is written
#'        as a message.
#' @param ... Additional arguments that are passed on \code{\link[boot]{boot}}.
#'   It is highly recommended to use its \code{parallel} and \code{ncpus} 
#'   options to parallelize the bootstrap replicates.
#' @return The result of \code{\link[boot]{boot}}. This object can be used to
#'   estimate standard errors or confidence intervals of the estimates using
#'   the functions available in package \pkg{boot}.
#'   Note that the returned object contains a vector of parameter values \code{t0}
#'   that is the result of an additional jaatha run for the original data, whereas
#'   the parametric bootstrap simulations used parameter values that are in the vector
#'   \code{mle} in the returned \code{boot} object.
#'   By default, the function \code{boot.ci} of the \code{\link[boot]{boot}} package
#'   uses the parameter values \code{t0} as a reference point.
#'   To use the values in \code{mle} instead, overwrite \code{t0} with \code{mle} before
#'   applying the function \code{boot.ci}.
#' @seealso 
#' \code{\link[boot]{boot}}, \code{\link{jaatha}}
#' 
#' @examples 
#' \donttest{
#' # The original Jaatha anaylsis:
#' model <- create_jaatha_model(function(x) rpois(10, x),
#'                     par_ranges = matrix(c(0.1, 0.1, 10, 10), 2, 2),
#'                     sum_stats = list(create_jaatha_stat("sum", sum)))
#' data <- create_jaatha_data(rpois(10, 9), model)
#' jaatha_result <- jaatha(model, data, cores = 2)
#' 
#' # Bootstrapping the results using 2 CPU cores on localhost:
#' library(boot)
#' library(snow)
#' cl <- makeSOCKcluster(c("localhost", "localhost"))
#' 
#' jaatha_boot_results <- boot_jaatha(jaatha_result, 50, 
#'                                    cores_per_run = 2,
#'                                    parallel = "snow",
#'                                    cl = cl)
#' 
#' stopCluster(cl)
#' boot.ci(jaatha_boot_results, type = "norm")
#' jaatha_boot_results$t0 <- jaatha_boot_results$mle
#' boot.ci(jaatha_boot_results, type = "norm")
#' }
#' 
#' @export
boot_jaatha <- function(results, R, cores_per_run = 1, verbose=TRUE, ...) {
  require_package("boot")
  if (R.Version()$major == 3 && R.Version()$minor < 2.2) {
    stop("This function requires at least R Version 3.2.2")
  }
  
  args <- results$args
  args$cores <- cores_per_run
  model <- args$model
  sim_func <- model$get_sim_func()
  
  log_folder <- tempfile("logs_")
  message("Logfiles in folder: ", tempdir())
  message("This might take a while...")
  dir.create(log_folder)
    
  jaatha_stat <- function(data) {
    args$data <- data  
    utils::capture.output(results <- do.call(jaatha, args),
                          file = tempfile(paste0("boot_log_", 
                                                 Sys.getpid(), "_")), 
                          type = "message")
    if(verbose) {
        message("Bootstrap estimation: ",
                paste(results$estimate, collapse=" "))
    }
   
    results$estimate
  }
  
  simulate_data <- function(data, param) {
    sim_data <- sim_func(param)
    create_jaatha_data.default(sim_data, model)
  }
  
  boot::boot(args$data, jaatha_stat, R = R, 
             sim = "parametric",
             ran.gen = simulate_data, 
             mle = results$estimate, ...)
}
