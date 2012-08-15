#!/usr/bin/Rscript --vanilla
#
# awesome_helper_functions
# A collection of nice litte helper functions
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-07-30
# Licence:  GPLv3 or later
#


#-----------------------------------------------------------------------
# Functions for easy log creation
#-----------------------------------------------------------------------

# Create a new enviroment for local variables that won't be looked after package
# loading like jaathas enviroment is.
if (!exists(".local")) .local <- new.env()
if (!exists('.local$log.level')) .local$log.level <- 1
if (!exists('.local$log.file'))  .local$log.file  <- ""

# A helper function for easy creation of logging output
#
# @param level An integer indicating the level of the log output.
#              If the level is smaller or equal to the log level stored in
#              'log.level', than the output is printed, otherwise it is
#              discarded. 0 is the default level.
# @param ...   One or more strings/variables to be written to the log stream
# @return      nothing
.log <- function(level, ...) {
  if (level > .local$log.level) return()

  if ( .local$log.level == 1 )  {
    # Normal output without dates
    cat(...,"\n",sep=" ")
    return()
  }

  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n", sep=" ",
      file=.local$log.file, append=T)

  if (.local$log.file != "" && level == 1) {
    # Also print normal output when logging to file
    cat(...,"\n",sep=" ")
  }
}

# Function for creating normal user output 
.print <- function(...){
    .log(0,...)
}
.log1 <- .print

# Creates lvl-2 logging output
.log2 <- function(...) {
  .log(2, ...)
}

# Creates lvl-3 logging output
.log3 <- function(...) {
  .log(3, ...)
}

#' Function to activate/change logging
#'
#' @param log.level An integer indicating the level of the log output.
#'              If the level is smaller or equal to the log level stored in
#'              'log.level', than the output is printed, otherwise it is
#'              discarded. 0 is the default level.
#' @param log.file If given, the logging output will be written into this file.
#' @return      nothing
setLogging <- function(log.level, log.file) {
  checkType(log.level, c("num", "s"), F)
  checkType(log.file, c("char", "s"), F)

  if (!missing(log.level)) .local$log.level <- log.level
  if (!missing(log.file)) .local$log.file <- log.file
}

#' Gets the current logging level
getLogLevel <- function() return(.local$log.level)
#' Gets the current logging file
getLogFile  <- function() return(.local$log.file)



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
  
  if (is.vector(variable) & !allow.na){  
    if (any(sapply(variable, is.na))) {
      fun.name <- as.character(sys.call(-1)[[1]])
      var.name <- deparse(substitute(variable))
      stop(fun.name,": Required parameter \"",var.name,"\" has NA value(s).",call.=F)
    }
  }

  for (i in seq(along=type)){
    if (type[i] == "char" || type[i] == "character") {
      func <- is.character
      error <- "has to be of type character"
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
      func <- function(var) {length(var) == 1}
      error <- "must have length one."
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


generateSeeds <- function(n=1) {
  return(sample.int(2^20,n))
}
