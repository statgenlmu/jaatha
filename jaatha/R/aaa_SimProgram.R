slots <- representation(name="character",
                        executable="character",
                        features="character",
                        simParFunc="function",
                        sumStatFunc="function",
                        calcJSFSFunc="function",
                        initialSeedFunc="function")

setClass("SimProgram",slots)
rm(slots)



#-----------------------------------------------------------------------
# Initialization
#-----------------------------------------------------------------------
emptyFunc <- function(){}

.init <- function(.Object, name, executable, features,
                  simParFunc, sumStatFunc, calcJSFSFunc,
                  initialSeedFunc=emptyFunc ){

  .Object <- sp.setName(.Object,name)
  .Object <- sp.setExecutable(.Object,executable)
  .Object <- sp.setFeatures(.Object,features)
  .Object <- sp.setSimParFunc(.Object,simParFunc)
  .Object <- sp.setSumStatFunc(.Object,sumStatFunc)
  .Object <- sp.setCalcJSFSFunc(.Object,calcJSFSFunc)
  .Object <- sp.setInitialSeedFunc(.Object,initialSeedFunc)
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
  if(!is.character(executable)) stop("executable must be of type character")
  exe.exist <- file.exists(executable)

  if(!any(exe.exist)) {
    warning(paste("No executable for simulation program",
                  simProg@name,"found."))
  } else {
    simProg@executable <- (executable[exe.exist])[1]
  }
  return(simProg)
}

sp.setFeatures <- function(simProg, features){
  if(!is.character(features)) stop("features must be of type character")
  simProg@features <- features
  return(simProg)
}

sp.setSimParFunc <- function(simProg, simParFunc){
  if(!is.function(simParFunc)) stop("simParFunc must be of type function")
  simProg@simParFunc <- simParFunc
  return(simProg)
}

sp.setSumStatFunc <- function(simProg, sumStatFunc){
  if(!is.function(sumStatFunc)) stop("sumStatFunc must be of type function")
  simProg@sumStatFunc <- sumStatFunc
  return(simProg)
}

sp.setCalcJSFSFunc <- function(simProg, calcJSFSFunc){
  if(!is.function(calcJSFSFunc)) stop("calcJSFSFunc must be of type function")
  simProg@calcJSFSFunc <- calcJSFSFunc
  return(simProg)
}

sp.setInitialSeedFunc <- function(simProg, initialSeedFunc){
  if(!is.function(initialSeedFunc)) stop("initialSeedFunc must be of type function")
  simProg@initialSeedFunc <- initialSeedFunc
  return(simProg)
}

#-----------------------------------------------------------------------
# Test
#-----------------------------------------------------------------------
#new("SimProgram","test","blub",c("mutation","split"),sum,sin,cos)
