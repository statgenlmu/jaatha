#' Function to set the number of cores that Jaatha uses for simulations
#' 
#' Jaatha can distribute the simulation on multiple cores if the operating
#' system supports p-threads (parrallel::mclapply is used internally). This
#' function sets the number of cores that should be used.
#' 
#' @param jaatha The jaatha object for which we set the number of cores
#' @param cores The number of cores to use. On Windows, only one core can be 
#'  used.
#' @return The jaatha object with set number of cores
setCores <- function(jaatha, cores=1) {
  checkType(jaatha, c("jaatha", "single"))
  checkType(cores, c("num", "single"))
  
  if (cores > 1 && .Platform$OS.type == "windows") {
    warning("Parallelization is not supported on Windows. The 'cores' option will be ignored") 
    jaatha@cores <- 1
  } else {
    jaatha@cores <- cores
  }
  
  invisible(jaatha)
}