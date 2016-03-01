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
#' @param ... Additional arguments that are passed on \code{\link[boot]{boot}}.
#'   It is highly recommended to use its \code{parallel} and \code{ncpus} 
#'   options to parallelize the bootstrap replicates.
#' @return The result of \code{\link[boot]{boot}}. This object can be used to
#'   estimate standard errors or confidence intervals of the estimates using
#'   the functions available in package \pkg{boot}.
#' @seealso 
#' \code{\link[boot]{boot}}, \code{\link{jaatha}}
#' 
#' @examples 
#' \dontrun{
#' # The original Jaatha anaylsis:
#' jaatha_result <- jaatha(model, data, cores = 4)
#' 
#' # Bootstrapping the results using 4 CPU cores on host1 and 2 on host2/host3:
#' library(boot)
#' library(snow)
#' cl <- makeSOCKcluster(c("host1", "host1", "host2", "host3"))
#' 
#' jaatha_boot_results <- boot_jaatha(jaatha_result, 200, 
#'                                    cores_per_run = 2,
#'                                    parallel = "snow",
#'                                    cl = cl)
#' 
#' stopCluster(cl)
#' boot.ci(jaatha_boot_results, type = "norm")
#' }
#' 
#' @export
boot_jaatha <- function(results, R, cores_per_run = 1, ...) {
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
    utils::capture.output(results <- do.call(jaatha, args),
                          file = tempfile(paste0("boot_log_", 
                                                 Sys.getpid(), "_")), 
                          type = "message")
    
    results$estimate
  }
  
  simulate_data <- function(data, param) {
    sim_data <- sim_func(param)
    create_jaatha_data.default(sim_data, model)
  }
  
  boot::boot(0, jaatha_stat, R = R, 
             sim = "parametric",
             ran.gen = simulate_data, 
             mle = results$estimate, ...)
}
