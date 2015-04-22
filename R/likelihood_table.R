#' Gives the best estimates after a Jaatha search
#'
#' This method extracts the best estimates with log composite likelihood
#' vales from an Jaatha object.
#'
#' @param jaatha The Jaatha options
#' @param max_rows If given, no more than this number of entries will be 
#'                returned.
#' @param initial_search If \code{TRUE}, then the estimates for the blocks in 
#'   the initial search will be reported.
#' @return A matrix with log composite likelihoods and parameters of The
#' best estimates
#' @export
Jaatha.getLikelihoods <- function(jaatha, max_rows = NULL, 
                                  initial_search = FALSE) {
  
  assert_that(is_jaatha(jaatha))
  
  if (!initial_search) lt <- jaatha@likelihoods_rs
  else lt <- jaatha@likelihoods_is
  
  lt <- sort_likelihood_table(lt, max_rows)
  lt[ ,-(1:2)] <- t(sapply(1:nrow(lt), function(n) {
    denormalize(lt[n,-(1:2), drop=F], jaatha)
  }))
  lt
}


sort_likelihood_table <- function(likelihoods, max_rows = NULL) {
  perm <- order(likelihoods[ , 1], decreasing = TRUE)
  perm <- perm[is.finite(likelihoods[perm , 1])]
  if (!is.null(max_rows)) perm <- perm[1:min(max_rows, length(perm))]
  likelihoods[perm, , drop = FALSE]
}


create_likelihood_table <- function(jaatha, n_rows) {
  lt <- matrix(-Inf, n_rows, nrow(jaatha@par.ranges) + 2)
  colnames(lt) <- c("log-likelihood", "block", getParNames(jaatha))
  lt
}
