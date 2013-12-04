# --------------------------------------------------------------
# Authors:  Paul R. Staab & Lisha Mathew
# Date:     2013-11-13
# Licence:  GPLv3 or later
# --------------------------------------------------------------

#' Convert parameters from the natural scale into Jaatha's internal scale
#'
#' @param value A parameter combination in this natural scale
#' @param jaatha The current jaatha object
#' @return The parameter from value, converted into Jaatha's interval 0-1 scale
normalize <- function(value, jaatha) {
  log.range <- log(jaatha@par.ranges) 
  log.value <- log(value)
  (log.value - log.range[,'min']) / (log.range[,'max'] - log.range[,'min'] ) 
}

#' Convert parameters from Jaatha's internal scale into their natural scale
#'
#' @param value A parameter combination in Jaatha's internal scale
#' @param jaatha The current jaatha object
#' @return The parameter from value, converted into their natural scale
denormalize <- function(value, jaatha) {
  log.range <- log(jaatha@par.ranges) 
  exp(value*(log.range[,'max']-log.range[,'min'])+log.range[,'min'])
}
