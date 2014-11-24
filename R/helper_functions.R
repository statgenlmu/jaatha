# ------------------------------------------------------------
# helper_functions.R
# A collection of litte helper functions
# 
# Author:   Paul R. Staab
# Email:    staab (at) bio.lmu.de
# Date:     2013-09-04
# Licence:  GPLv3 or later
# ------------------------------------------------------------


#-----------------------------------------------------------------------
# Create and Manage an enviroment for run-time variables
#-----------------------------------------------------------------------

# Create a new enviroment for local variables that won't be looked after package
# loading like the default package enviroment is.
if (!exists(".jaatha")) .jaatha <- new.env()

setJaathaVariable <- function(name, value) {
  .jaatha[[name]] <- value
}

getJaathaVariable <- function(name) {
  return(get(name, envir=.jaatha))
}

isJaathaVariable <- function(name) {
  exists(name, envir=.jaatha)
}



# Function for creating normal user output 
.print <- function(...) { cat(..., "\n") }


#-----------------------------------------------------------------------
# Function to conveniently check the type of user inputs
#-----------------------------------------------------------------------

#' Checks if a variable is of a given type and calls stop() on type mismatch
#'
#' See heading.
#' Side-effect warning: Calls stop() on type mismatch.
#'
#' @param variable the variable to check
#' @param type the name of the type the variable should have. Can be num/numeric, vec/vector
#'             mat/matrix or fun/function. Can also be s/single, in that case it
#'             must be a vector of length one . 
#'             If a vector of type names is given, the variable must be of all types.
#' @param required A boolean that indicates whether the variable must be
#'             specified or can be missing. Value of NULL also counts as missing.
#' @param allow.na If FALSE, an error is returned if the variable contains NA
#'             values. Only works for vectors.
#'
#' @return nothing
checkType <- function(variable, type, required=T, allow.na=T) {
  if (missing(variable)) {
    if (!required) return()
    fun.name <- as.character(sys.call(-1)[[1]])
    var.name <- deparse(substitute(variable))
    stop(fun.name,": Required parameter \"",var.name,"\" is missing.",call.=F)
  } 
      
  if (is.null(variable)) {
    if (!required) return()
    fun.name <- as.character(sys.call(-1)[[1]])
    var.name <- deparse(substitute(variable))
    stop(fun.name,": Required parameter \"",var.name,"\" is NULL.",call.=F)
  } 
  
  if (is.vector(variable)) {
    if (allow.na) {
      variable <- variable[!is.na(variable)]
      if (length(variable) == 0) return()
    } else if (any(sapply(variable, is.na))) {
      fun.name <- as.character(sys.call(-1)[[1]])
      var.name <- deparse(substitute(variable))
      stop(fun.name,": Required parameter \"",var.name,"\" has NA value(s).",call.=F)
    }
  }

  for (i in seq(along=type)){
    if (type[i] == "char" || type[i] == "character") {
      func <- is.character
      error <- "has to be of type character"
    } else if (type[i] == "bool" || type[i] == "boolean") {
      func <- is.logical
      error <- "has to be of type boolean"
    } else if (type[i] == "int" || type[i] == "integer") {
      func <- is.integer
      error <- "has to be of type integer"
    } else if (type[i] == "num" || type[i] == "numeric") {
      func <- is.numeric
      error <- "has to be of type numeric"
    } else if (type[i]  == "vec" || type[i]  == "vector") {
      func <- is.vector
      error <- "has to be a vector"
    } else if (type[i] == "mat" || type[i] == "matrix") {
      func <- is.matrix
      error <- "has to be a matrix"
    } else if (type[i] == "fun" || type[i] == "function") {
      func <- is.function
      error <- "has to be a function"
    } else if (type[i] == "s" || type[i] == "single") {
      func <- function(x) {length(x) == 1}
      error <- "must have length one."
    } else if (type[i] == "dm" || type[i] == "demographicModel") {
      func <- function(x) {class(x)[1] == "DemographicModel"}
      error <- "is no demographic Model"
    } else if (type[i] == "jat" || type[i] == "jaatha") {
      func <- function(x) {class(x)[1] == "Jaatha"}
      error <- "is no Jaatha object"
    } else if (type[i] == "ar" || type[i] == "array") {
      func <- is.array 
      error <- "is not an array"
    } else {
      stop("Unknown type: ",type[i])
    }
    if (!func(variable)) {
      fun.name <- as.character(sys.call(-1)[[1]])
      var.name <- deparse(substitute(variable))
      stop(fun.name,": ",var.name," ",error,call.=F)
    }
  }
}


#-----------------------------------------------------------------------
# Functions to deal with temporary files
#-----------------------------------------------------------------------
getTempFile <- function(file.name="file") {
  tempfile(paste0('jaatha', '_', Sys.getpid(), "_", file.name, "_"))
}


#-----------------------------------------------------------------------
# Function to generate seeds
#-----------------------------------------------------------------------
generateSeeds <- function(n=1) {
  return(sample.int(2^20,n))
}
