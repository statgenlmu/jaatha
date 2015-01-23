Parameter <- R6Class('Parameter', 
  private = list(
    expr = NA,    
    name = NA
  ),
  public = list(
    initialize = function(expr, name=NA) {
      if (!is.expression(expr)) stop("No expression provided: ", expr,
                                     " is of type ", is(expr))
      private$expr <- expr
      
      if (!(is.na(name) | is.character(name)))
        stop('The parameter name must be a character')
      private$name <- name
    },
    eval = function(envir = parent.frame()) {
      eval(private$expr, envir = envir)
    },
    get_expression = function() private$expr,
    get_name = function() private$name
  )
)

#' Define Model Parameters
#' 
#' bla bla bla
#' 
#' @param expr An R command. The command will not be evaluted until a simulation
#'  is performed. I can contain other named parameters, but not parameters 
#'  created with \code{par_expr}. Make sure that the expression always evaluates 
#'  to a suitable parameter value (a single numeric in almost all cases).
#' @describeIn par_expr A parameter whichs value is determined by evalutating an
#'  expression.
#' @export
#' @aliases ModelParameters
#' @author Paul Staab
#' @examples
#' par_expr(5)          # The parameters value is always 5.
#' par_expr(runif(1))   # Creates an parameter which takes a uniformly
#'                      # distributed value in each simulation.    
#' par_range('x', 1, 5) # Creates an parameter with name x with possible values 
#'                      # between 1 and 5.
#' par_expr(2*x)        # The parameters value is always equal two twice the 
#'                      # value of a different model parameter named 'x'.
par_expr <- function(expr) {
  Parameter$new(as.expression(substitute(expr)))
}

#' @describeIn par_expr Creates an parameter that can take a range of possible
#'  values. Used for creating model parameters for \code{Jaatha}.
#' @export
#' @param name A string. The name of the parameter.
#' @param lower A numeric. The lower boundary of the parameter's range.
#' @param upper A numeric. The upper boundary of the parameter's range.
par_range <- function(name, lower, upper) {
  
}