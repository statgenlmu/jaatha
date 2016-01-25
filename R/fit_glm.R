fit_glm <- function(x, sim_data, ...) UseMethod("fit_glm")
fit_glm.default <- function(x, sim_data, ...) {
  stop("Unknown Summary Statistic")
}


#' @export
fit_glm.jaatha_model <- function(x, sim_data, ...) { #nolint
  "Fits a GLM to the simulation results"
  lapply(x$get_sum_stats(), fit_glm, sim_data, ...)
}


#' @export
#' @importFrom stats glm.fit poisson
fit_glm.jaatha_stat_basic <- function(x, sim_data, ...) {
  "Fits a GLM for each entry of the simulation results"
  Y <- do.call(rbind, lapply(sim_data, function(data) data[[x$get_name()]]))
  X <- cbind(1, 
             do.call(rbind, lapply(sim_data, function(data) data$pars_normal)))
  
  glms <- lapply(1:ncol(Y), function(i) {
    suppressWarnings(
      glm.fit(X, Y[, i], family = poisson("log"), 
              control = list(maxit = 100))[c("coefficients", "converged")]
    )
  })
  
  sapply(glms, function(x) {
    if (!x$converged) stop("GLM did not converge", call. = FALSE)
    if (all(abs(x$coefficients[-1]) < 1e-10)) {
      stop("GLM coefficients are all numerically 0", call. = FALSE)
    }
    NULL
  })

  glms
}


# fit_glm.Stat_PoiSmooth <- function(sum_stat, sim_data) {
#   par_names <- names(sim_data[[1]]$pars)
#   model <- paste0("sum.stat ~ ",
#                   "(", sum_stat$get_model(), ")",  
#                   "*(", paste(par_names, collapse="+"), ")") 
# 
#   sim_data_df <- do.call(rbind, lapply(sim_data, function(sim_result) {
#     pars <- matrix(sim_result$pars.normal, 1,
#                    length(sim_result$pars.normal), byrow=TRUE)
#     colnames(pars) <- names(sim_result$pars.normal)
#     data.frame(pars, sim_result[[sum_stat$get_name()]])
#   }))
# 
#   smooth_glm  <- glm(model, data=sim_data_df, family=poisson("log"), 
#                      model = FALSE, x = FALSE, y = FALSE,
#                      control = list(maxit = 200))
#   if (!smooth_glm$converged) stop("GLM did not converge")
#   
#   list(smooth=smooth_glm)
# }
