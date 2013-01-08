#---------------------------------------------------------------
# aaa_SimProgram.R
# Abstract class used for connecting a sequence simulation program 
# to the Demographic Model class. Needs to be loaded as first source
# file in the package, hence the "aaa".
# 
# Author:   Paul R. Staab & Lisha Naduvilezhath 
# Email:    staab (at) bio.lmu.de
# Date:     2012-10-05
# Licence:  GPLv3 or later
#----------------------------------------------------------------

# Keep a user modifiable list of availible simulation programs in a private
# enviroment (jaatha's own env is read-only after package load)
if (!exists(".jaatha")) .jaatha <- new.env()
if (!exists("simProgs", envir=.jaatha)) .jaatha$simProgs <- list()

slots <- representation(name="character",
                        executable="character",
                        possible.features="character",
                        possible.sum.stats="character",
                        simFunc="function",
                        singleSimFunc="function",
                        useSingleSimFunc="logical",
                        finalizationFunc="function"
                       )

setClass("SimProgram",slots)
rm(slots)


#-----------------------------------------------------------------------
# Initialization
#-----------------------------------------------------------------------
emptyFunc <- function(){}

.init <- function(.Object, 
                  name=NULL, 
                  executable=NULL, 
                  possible.features=NULL,
                  possible.sum.stats=NULL, 
                  simFunc=NULL,
                  singleSimFunc=NULL,
                  finalizationFunc=NULL) {

  if (is.null(simFunc) & is.null(singleSimFunc))
      stop("One of simFunc and singleSimFunc must be given.")
  if (!is.null(simFunc) & !is.null(singleSimFunc))
      stop("Only one of simFunc and singleSimFunc can be used.")
  if (!is.null(simFunc)) 
    .Object <- sp.setSimFunc(.Object, simFunc)
  if (!is.null(singleSimFunc)) 
    .Object <- sp.setSingleSimFunc(.Object, singleSimFunc)

  .Object <- sp.setFinalizationFunc(.Object, finalizationFunc)

  .Object <- sp.setName(.Object, name)
  .Object <- sp.setExecutable(.Object, executable)
  .Object <- sp.setPossibleFeatures(.Object, possible.features)
  .Object <- sp.setPossibleSumStats(.Object, possible.sum.stats)
  return(.Object)
}

setMethod("initialize","SimProgram",.init)
rm(.init)



#-----------------------------------------------------------------------
# Setters
#-----------------------------------------------------------------------
sp.setName <- function(simProg, name){
  if(!is.character(name)) stop("name must be of type character")
  if(length(name) != 1) stop("name must not be a vector")
  simProg@name <- name
  return(simProg)
}

sp.setExecutable <- function(simProg, executable){
  checkType(executable, "char")
  
  if (executable == "") {
    simProg@executable <- ""
    return(simProg)
  }

  exe.exist <- file.exists(executable)
  
  if(!any(exe.exist)) {
    warning(paste("No executable for simulation program",
                  simProg@name,"found."))
  } 
  
  simProg@executable <- (executable[exe.exist])[1]
  return(simProg)
}

sp.setPossibleFeatures <- function(simProg, features){
  checkType(features, "char")
  simProg@possible.features <- features
  return(simProg)
}

sp.setPossibleSumStats <- function(simProg, sum.stats){
  checkType(sum.stats, "char")
  simProg@possible.sum.stats <- sum.stats
  return(simProg)
}

sp.setSimFunc <- function(simProg, simFunc){
  checkType(simFunc, "fun", F)
  simProg@simFunc <- simFunc
  simProg@useSingleSimFunc <- F
  return(simProg)
}

sp.setSingleSimFunc <- function(simProg, simFunc){
  checkType(simFunc, "fun", F)
  simProg@singleSimFunc <- simFunc
  simProg@useSingleSimFunc <- T
  return(simProg)
}

sp.setFinalizationFunc <- function(simProg, finalFunc){
  if (is.null(finalFunc)) finalFunc <- function(dm){return(dm)}
  checkType(finalFunc, "fun", F)
  simProg@finalizationFunc <- finalFunc
  return(simProg)
}

createSimProgram <- function(name, executable,
                             possible.features,
                             possible.sum.stats,
                             simFunc=NULL,
                             singleSimFunc=NULL,
                             finalizationFunc=NULL) {
  
  simProg <- new("SimProgram",
                 name = name, 
                 executable = executable,
                 possible.features = possible.features,
                 possible.sum.stats = possible.sum.stats,
                 simFunc = simFunc,
                 singleSimFunc = singleSimFunc,
                 finalizationFunc = finalizationFunc) 

  .jaatha$simProgs[[name]] <- simProg
}

getSimProgram <- function(name) {
  return(.jaatha$simProgs[[name]])
}
