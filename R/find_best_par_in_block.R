# --------------------------------------------------------------
# Authors:  Paul R. Staab & Lisha Mathew
# Date:     2013-11-14
# Licence:  GPLv3 or later
# --------------------------------------------------------------

#' Functions that uses the fitted GLMs to estimate the parameters that have the
#' highest likelihood.
#'
#' @param block The block we are in.
#' @param glm.fitted The fitted GlMs.
#' @param sum.stats The observed summary statistics
#' @return A list with maximum likelihood parameter (est) and log-likelihood
#' (score)
findBestParInBlock <- function(block, glm.fitted, sum.stats) {
  block.size <- block@border[,2] - block@border[,1]
  block.middle <- block.size/2 + block@border[,1]

  ##describes 'boarder'% of values that will be excluded
  ##on either side of the block in optimization
  best.value <- optim(block.middle, estimateLogLikelihood, 
                      glm.fitted=glm.fitted, sum.stats=sum.stats,  
                      lower=block@border[ ,1], upper=block@border[ ,2],
                      method="L-BFGS-B", control=list(fnscale=-1))

  return(list(est=best.value$par, score=best.value$value))                   
}
