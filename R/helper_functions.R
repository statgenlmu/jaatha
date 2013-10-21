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



#-----------------------------------------------------------------------
# Functions for easy log creation
#-----------------------------------------------------------------------

if (!exists('log.level', envir=.jaatha)) .jaatha$log.level <- 1
if (!exists('log.file', envir=.jaatha))  .jaatha$log.file  <- ""

# A helper function for easy creation of logging output
#
# @param level An integer indicating the level of the log output.
#              If the level is smaller or equal to the log level stored in
#              'log.level', than the output is printed, otherwise it is
#              discarded. 0 is the default level.
# @param ...   One or more strings/variables to be written to the log stream
# @return      nothing
.log <- function(level, ...) {
  if (level > .jaatha$log.level) return()

  if ( .jaatha$log.level == 1 )  {
    # Normal output without dates
    cat(...,"\n",sep=" ")
    return()
  }

  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n", sep=" ",
      file=.jaatha$log.file, append=T)

  if (.jaatha$log.file != "" & level == 1) {
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

  if (!missing(log.level)) .jaatha$log.level <- log.level
  if (!missing(log.file)) .jaatha$log.file <- log.file
}

#' Gets the current logging level
getLogLevel <- function() return(.jaatha$log.level)
#' Gets the current logging file
getLogFile  <- function() return(.jaatha$log.file)



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
      func <- function(x) { is.character(x) }
      error <- "has to be of type character"
    } else if (type[i] == "bool" || type[i] == "boolean") {
      func <- function(x) { is.logical(x) }
      error <- "has to be of type boolean"
    } else if (type[i] == "num" || type[i] == "numeric") {
      func <- function(x) { is.numeric(x) }
      error <- "has to be of type numeric"
    } else if (type[i]  == "vec" || type[i]  == "vector") {
      func <- function(x) { is.vector(x) }
      error <- "has to be a vector"
    } else if (type[i] == "mat" || type[i] == "matrix") {
      func <- function(x) { is.matrix(x) }
      error <- "has to be a matrix"
    } else if (type[i] == "fun" || type[i] == "function") {
      func <- function(x) { is.function(x) }
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
getTempDir <- function(use.shm = FALSE) {
  if (exists("temp.dir", envir=.jaatha) & !use.shm) {
    return(.jaatha$temp.dir)
  }

  if (use.shm) {
    if (!file.exists("/dev/shm")) stop("/dev/shm/ does not exists") 
    tmp.dir <- "/dev/shm/jaatha"
  }
  else tmp.dir <- paste(tempdir(), "/jaatha", sep="")
  
  i <- 1
  while (file.exists(paste(tmp.dir, "-", i, sep="")))
    i <- i + 1

  .jaatha$temp.dir <- paste(tmp.dir, "-", i, sep="")
  dir.create(.jaatha$temp.dir)
  return(.jaatha$temp.dir)
}

getTempFile <- function(file.name="file"){
  if (!exists("temp.file.count", envir=.jaatha))
      .jaatha$temp.file.count <- 0

  .jaatha$temp.file.count <- .jaatha$temp.file.count + 1 %% 1000000
  return(paste(getTempDir(), "/", file.name, "_", Sys.getpid(), "_", .jaatha$temp.file.count, sep=""))
}

removeTempFiles <- function() {
  temp.dir <- NULL
  if (exists("temp.dir", envir=.jaatha)) {
    unlink(.jaatha$temp.dir, recursive=T)
    rm(temp.dir, envir=.jaatha)
  }
}



#-----------------------------------------------------------------------
# Function to generate seeds
#-----------------------------------------------------------------------
generateSeeds <- function(n=1) {
  return(sample.int(2^20,n))
}
