#' Function that uses the fitted GLMs to estimate maximum likelihood 
#' parameters within a block.
#'
#' @param block The block we are in.
#' @param glm_fitted The fitted GlMs.
#' @param sum_stats The observed summary statistics
#' @param scaling_factor The scaling factor used for the simulations
#' @return A list with maximum likelihood parameter (\code{est}) 
#'   and log-likelihood (\code{score})
#' @author Paul Staab & Lisha Mathew
findBestParInBlock <- function(block, glm_fitted, sum_stats, scaling_factor=1) {
  border <- block$get_border()
  
  best.value <- optim(block$get_middle(), 
                      estimateLogLikelihood, 
                      glm_fitted = glm_fitted, 
                      sum_stats = sum_stats,
                      scaling_factor = scaling_factor,
                      lower = border[ , 1, drop = FALSE], 
                      upper = border[ , 2, drop = FALSE],
                      method = "L-BFGS-B", 
                      control = list(fnscale = -1))

  assert_that(block$includes(best.value$par))
  list(est = best.value$par, score = best.value$value)
}
