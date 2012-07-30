setClass("DemographicModel" ,
     representation(features="data.frame",
            parameters="character",
            sampleSizes="numeric",
            nLoci="numeric",
            seqLength="numeric",
            tsTvRatio="numeric",
            externalTheta="logical",
            finiteSites="logical",
            profiled="logical",
            simProgs="list",
            currentSimProg="SimProgram")
    )


if(!exists("dm.defaultSimProgs")) {
  dm.defaultSimProgs <- list()
}

#-----------------------------------------------------------------------
# Initialization
#-----------------------------------------------------------------------

.init <- function(.Object, sampleSizes, nLoci, seqLength,
                  finiteSites, tsTvRatio){
    .Object@features <- data.frame( type=character(),
                        parameter=numeric(),
                        population=numeric(),
                        lowerRange=numeric(),
                        upperRange=numeric(),
                        startsAtTime=numeric(),
                        timeLine=numeric()
                      )

    .Object@parameters      <- character()
    .Object@externalTheta   <- F
    .Object@profiled        <- F
    .Object@finiteSites     <- finiteSites
    .Object@sampleSizes     <- sampleSizes
    .Object@seqLength       <- seqLength
    .Object@tsTvRatio       <- tsTvRatio
    .Object@nLoci           <- nLoci

    .Object@simProgs        <- dm.defaultSimProgs
    
    return(.Object)
}

setMethod("initialize","DemographicModel",.init)
rm(.init)

#-----------------------------------------------------------------------
# Print
#-----------------------------------------------------------------------

.show <- function(object){
  dm <- object

  features <- dm@features[!is.na(dm@features$parameter),]

  if (nrow(features) == 0) {
    cat("Your demographic model has no features so far.\n")
    return()
  }

  cat("Your demographic model has the following parameters:\n")
  for (rw in 1:nrow(features)) {
    type <- features[rw,'type'] 
    lR <- as.character(features[rw,'lowerRange']) 
    uR <- as.character(features[rw,'upperRange'])
    parname <- dm@parameters[features[rw,'parameter']]
    
    if (type == "split")
      cat("-",parname,": A split into two populations between",lR,"and",uR,
            "Ne generations ago.\n")
    else if (type == "mutation")
      cat("-",parname,": A scaled mutation rate between",lR,"and",uR,"\n")
    else if (type == "recombination")
      cat("-",parname,": A scaled recombination rate between",lR,"and",uR,"\n")
    else if (type == "migration")
      cat("-",parname,":  A scaled migration rate between",lR,"and",uR,"\n")
    else if (type == "presentSize")
      cat("-",parname,":  At sample time, the second population two was between",
          lR,"and",uR,"times as the first one\n")
    else if (type == "splitSize") 
      cat("-",parname,":  At time of the population split, the second",
          "population two was between",lR,"and",uR,"times as the first one\n")
    else if (type == "migration")
      cat("-",parname,":  A scaled migration rate between",lR,"and",uR,"\n")
    else {
      cat("- Other:\n")
      print(dm@features[rw,])
    }
  }


  features <- dm@features[is.na(dm@features$parameter),]
  if (dim(features)[1] > 0){
    cat("\n")
    cat("Additionally, it has the following features:\n")
     for (rw in 1:(dim(features)[1])) {
       type <- features[rw,'type']
       lR <- as.character(features[rw,'lowerRange']) 
       
       if (type == "split")
         cat("- A split into two populations",lR,"Ne generations ago.\n")
      else if (type == "mutation")
         cat("- A fixed scaled mutation rate of",lR,"\n")
       else if (type == "recombination")
         cat("- A fixed scaled recombination rate of",lR,"\n")
       else if (type == "migration")
         cat("- A fixed scaled migration rate of",lR,"\n")
      else if (type == "presentSize")
        cat("- At sample time, the second population two was",
            lR,"times as the first one\n")
      else if (type == "splitSize") 
        cat("- At time of the population split, the second population
            two was",lR,"times as the first one\n")
       else {
         cat("- Other:\n")
         print(dm@features[rw,])
       }
    }
  }
} 
setMethod("show","DemographicModel",.show)
rm(.show)


  
#-----------------------------------------------------------------------
# "Private" functions
#-----------------------------------------------------------------------

.appendFeature <- function(dm,type,parameter,population=NA,
            lowerRange,upperRange,
            startsAtTime,timeLine){
    
    if (dim(.getFeature(dm,type,population))[1] > 0){
        stop("model feature already exists")
    }

    dm@features <- rbind(dm@features, data.frame(type=type,
                           parameter=parameter,
                           population=population,
                           lowerRange=lowerRange, 
                           upperRange=upperRange, 
                           startsAtTime=startsAtTime,
                           timeLine=timeLine)
                )
    return(dm)
}

.getFeature <- function(dm,type,population=NA){
    .check.dm(dm)
    if (is.na(population))
        return(dm@features[dm@features$type == type 
                       & is.na(dm@features$population),])
    else
        return(dm@features[dm@features$type == type 
                       & dm@features$population == population,])
}
    

.addParameter <- function(dm,parameterName){
    .check.dm(dm)
    if (any(dm@parameters == parameterName)) 
        stop("parameter already exists")
    dm@parameters <- c(dm@parameters,parameterName)
    dm <- .dm.selectSimProg(dm)
    return(dm)
}

.addFeature <- function(dm,type,parName,lowerRange,upperRange,fixedValue,
               population=NA,startsAtTime=NA,timeLine=NA){

    .check.dm(dm)

    #Set range
    if ( ( !missing(fixedValue) & !missing(lowerRange) & !missing(upperRange) ) |
         ( missing(fixedValue) & ( missing(lowerRange) | missing(upperRange) ) ) 
       )
        stop("Exactly either fixedValue or lowerRange and upperRange must be specified.")

    if ( !missing(fixedValue) ) {
        lowerRange <- fixedValue
        upperRange <- fixedValue
    }

    parameter <- NA

    if ( missing(fixedValue) ) {
        dm <- .addParameter(dm,parName)
        parameter <- length(dm@parameters)
    }

    dm <- .appendFeature(dm,type,parameter,population,
                     lowerRange,upperRange,
                 startsAtTime,timeLine)
    
    dm <- .makeThetaLast(dm)

    return(dm)
}

# 
# .setSimProg <- function(dm){
#     .check.dm(dm)
#     if (dm@finiteSites) simProg <- "fsc"
#     else            simProg <- "ms"
#     return(simProg)
# }

.makeThetaLast <- function(dm){
    .check.dm(dm)
    mut <- .getFeature(dm,"mutation",NA)
    if ( dim(mut)[1] == 0 ) return(dm)  #Yet no mutation added
    if ( dim(mut)[1] > 1  ) stop("Are there multiple mutation entries in your model?")
    mutPar <- mut$parameter
    if ( is.na(mutPar) ) return(dm)     #Looks like a fixed mutation rate
    nPar <- length(dm@parameters)
    if ( mutPar == nPar ) return(dm)    #Already last parameter
    #If we are here, theta should not be the last parameter
    mutLine <-  dm@features$type == "mutation"
    otherLine <- dm@features$parameter == nPar
    dm@features$parameter[mutLine] <- nPar
    dm@features$parameter[otherLine] <- mutPar
    dm@parameters[c(mutPar,nPar)] <- dm@parameters[c(nPar,mutPar)]
    return(dm)
}

#----------------------------------------------------------------
# .getParameter
#----------------------------------------------------------------
# returns the parameter that belong to a feature specified via 
# type and population
#
# dm:       the demographic model
# type:     the type of the feature
# population:   the population if the feature or NA
# 
# returns:  0           if feature is not in the model,
#       the fixed Value     if the feature has a "fixedValue" or
#       the parameter name  otherwise
#
.getParameter <- function(dm,type,population=NA){
    .check.dm(dm)
    feature <- .getFeature(dm,type,population)
    if ( dim(feature)[1] == 0 ) return(0)
    if ( feature$lowerRange[1] == feature$upperRange[1] )
        return(feature$lowerRange[1])
    else
        return(dm@parameters[feature$parameter[1]])
}

.getNumOfHistoricalEvent <- function(dm){
    .check.dm(dm)
    return(sum(dm@features$type != "effPopSize" & dm@features$type != "mutation"))
}


#.isFixedFeature <- function(dm,type,population=NA){
#   feature <- .getFeature(dm,type,population)
#   return(feature@lowerRange == feature@upperRange)
#}

.checkFeatureParameters <- function(lowerRange,upperRange,fixedValue,startsAtTime,timeLine){
    if (! .is.NumericOrFalse(lowerRange) ) cat("Wrong value for 'lowerRange'\n")
    if (! .is.NumericOrFalse(upperRange) ) cat("Wrong value for 'upperRange'\n")
    #if (! .is.NumericOrFalse(startsAtTime) ) cat("Wrong value for 'startsAtTime'\n")
    #if (! .is.NumericOrFalse(timeLine) ) cat("Wrong value for 'timeLine'\n")
    if (! .is.NumericOrFalse(fixedValue) ) cat("Wrong value for 'fixedValue'\n")


    if ( fixedValue == F & ( lowerRange == F | upperRange == F ) )
        cat("Wrong parameters. Either both lowerRange and upperRange 
             or fixedValue must be specified\n")

    if ( fixedValue != F && lowerRange != F && upperRange != F ) 
        cat("Wrong parameters. Either enter a fixed value or a range for this parameter \n")
}


.parseFeatureParametes <- function(lowerRange,upperRange,fixedValue,startsAtTime,timeLine){
    if ( startsAtTime == F ) startsAtTime <<- Inf
    if ( timeLine == F ) timeLine <<- Inf

    if ( is.numeric(fixedValue) ) {
        lowerRange <<- fixedValue
        upperRange <<- fixedValue
    }
    return(c(lowerRange,upperRange,fixedValue,startsAtTime,timeLine))
}


.is.NumericOrFalse <- function(x){
    return ( (is.numeric(x) || x == F) && length(x) == 1 )
}

.calcSizeChange <- function(dm,param,percentages){
    .check.dm(dm)
    pop <- dm@features[param,"population"]
    q <- .calcAbsParamValue(dm,percentages[param],param)
    
    # get value of "size" parameter s for the population
    sizeLine <- (1:dim(dm@features)[1])[dm@features$type=="size" & dm@features$population==pop]
    if (sum(sizeLine) == 1) { 
        s <- .calcAbsParamValue(dm,percentages[sizeLine],sizeLine) 
    } else { 
        s <- 1 
    }
    
    # get tau (XXX GLOBAL TAU ONLY ATM XXX)
    tauLine <- (1:dim(dm@features)[1])[dm@features$type=="mutation" & is.na(dm@features$population)]
    tau <- .calcAbsParamValue(dm,percentages[tauLine],tauLine) 

    return(log(q/s)/tau)
}

#.addNeIfMissing <- function(dm){
#   .check.dm(dm)
#   if (.getParameter(dm,"effPopSize",NA) == 0)
#       dm <- dm.addAncPopSize(dm,NA,fixed=10000)
#   return(dm)
#}


.is.dm  <- function(dm){
        return(class(dm)[1] == "DemographicModel")
}

.check.dm <- function(dm){
    if (!.is.dm(dm)) stop("dm is no demographic model")
}

.calcAbsParamValue <- function(dm,percentage,param){

    if (!is.matrix(percentage)) 
        percentage <- matrix(percentage,1)

    return( t(apply(percentage,1,.calcAbsParamValueLine,dm=dm)) )
}

.calcAbsParamValueLine <- function(dm,percentage,param){
    features <- dm@features[!is.na(dm@features$parameter),]

    #If not specified which param to convert, use the first length(percentage) ones
    if ( missing(param) ) param <- 1:length(percentage)
    
    #Calculate scaledValues
    scaledValues <- exp ( percentage * 
                  ( log(features[param,"upperRange"]) 
                    - log(features[param,"lowerRange"]) )) * 
                    features[param,"lowerRange"]

    #Add names
    names(scaledValues) <- dm@parameters[param]

    #Do not scale theta if externalTheta = T
    if ( dm@externalTheta ) {
        #TRUE iff param[.] links to theta
        theta.mask <- param == (1:length(dm@parameters))[dm@parameters == "theta"]
        if (any(theta.mask)) scaledValues[theta.mask] <- percentage[theta.mask]
    }
    
    return( scaledValues )
}

.checkParInRange <- function(dm, param) {
    ranges <- dm.getParRanges(dm,inklExtTheta=F)
    #Seems there can be rounding errors during scalation
    lower <- matrix(ranges[,1],dim(param)[1],dim(param)[2],byrow=T) - 1e-15
    upper <- matrix(ranges[,2],dim(param)[1],dim(param)[2],byrow=T) + 1e-15
    inRange <- lower <= param & param <= upper
    return(all(inRange))
}


.dm.selectSimProg <- function(dm) {
  for (i in seq(along = dm@simProgs)){
    if (all(dm@features$type %in% dm@simProgs[[i]]@features)) {
      dm@currentSimProg <- dm@simProgs[[i]]
      .log1("Using",dm@simProgs[[i]]@name,"for simulations")
      return(dm)
    }
  }
  warning("No suitable simulation software found!")
  return(dm)
}

#-----------------------------------------------------------------------
# Public functions
#-----------------------------------------------------------------------

dm.getParRanges <- function(dm,inklExtTheta=T){
    .check.dm(dm)
    parMask <- !is.na(dm@features$parameter)
    parRanges <- cbind(lower=dm@features$lowerRange[parMask],upper=dm@features$upperRange[parMask])
    parRanges <- parRanges[sort.list(dm@features$parameter[parMask]),]
    if (!inklExtTheta) parRanges <- parRanges[1:dm.getNPar(dm),]
    return( parRanges )
}

#' Create a basic demographic model
#' 
#' This function creates a basic empty demographic model, which
#' is returned. Features like mutation, population splits and 
#' migration can be added afterwards.
#'
#' @param sampleSizes Number of individual that are sampled. If your model 
#'            has multiple populations, this needs to be a vector
#'            containing the sample sizes from each population.
#' @param nLoci       Number of loci that will be simulated
#' @param seqLength   Number of bases for each locus
#' @param finiteSites If 'TRUE', a finite sites mutation model is assumed
#'            instead of an infinite sites one.
#' @param tsTvRatio   Transition transversion ratio
#' @param log.level   An integer specifing the verbosity of the object.
#'                    0 = no output, 1 = normal verbosity, ..., 
#'                    3 = maximal verbosity.
#' @param log.file    If set, the debug output will be written into the given file
#' @return            The demographic model
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.createDemographicModel <- 
  function(sampleSizes, nLoci, seqLength=1000, finiteSites=F,
           tsTvRatio=.33, log.level, log.file){
    setLogging(log.level, log.file)
    dm <- new("DemographicModel",sampleSizes,nLoci,seqLength,
          finiteSites,tsTvRatio)
    return(dm)
}

#' Creates a demographic model by giving the command line parameters 
#' of a simulation program.
#' 
#' In case your model can not be described in terms of the demographic model
#' class or you want to use a simulation program that is not yet supported, you
#' can use this function to create a so called custom model. To do so, you have
#' to state number (par.number) and ranges (par.ranges) of the parameters you
#' want to estimate as well that the simulation program you want to use (sim.exe).
#' Then, you have to create a function that generates the command-line options
#' for the program (sim.options). Finally, you need a function that reads the output
#' of your program and calculates the summary statistics (sumstats).
#' 
#' @param par.number   The number of parameters you want to estimate
#' @param par.names    A character vector with names for your parameters
#' @param par.ranges   The upper and lower boundaries for the values if your parameters.
#'                     This should be a matrix with two columns for the lower and
#'                     upper boundary respectively and a row for each parameter.
#' @param sim.exe      The path of the executable of your simulation program. 
#'                     E.g. "/usr/local/bin/ms" or "C:/msdir/ms.exe"
#' @param sim.options  
#' @param sumstats
#' @return             A demographic model you can use for jaatha
#' @export
dm.createCustomModel <- function(par.number, par.names, par.ranges, sim.exe, 
                                 sim.options, sumstats) {

  checkType(par.number, "num")
  checkType(par.names, c("char","vec"))
  checkType(par.ranges, c("num","mat"))
  checkType(sim.exe, "fun")
  checkType(sim.options, "fun")
  checkType(sumstats, "fun")

  if (length(par.names) != par.number) stop("par.names has wrong length")
  if (all(dim(par.ranges) != c(par.number,2)))
    stop("par.ranges need to a ",par.number,"x2 matrix")

  dm <- dm.createDemographicModel(0,0)
  
}



dm.setExternalTheta <- function(dm){
     .check.dm(dm)
     #newDm <- new("DemographicModel",dm@seqLength,dm@sampleSizes,dm@finiteSites)
     #newDm@features <- dm@features[dm@features$type != "mutation",]
     #newDm@parameters <- dm@parameters[dm@parameters != "theta"]
     #newDm <- dm.addMutation(newDm)
     dm@externalTheta = T
     return(dm)
}


dm.getNPar <- function(dm){
    .check.dm(dm)
    return(length(dm@parameters)-dm@externalTheta)
}



#' Adds mutations to a demographic model
#' 
#' @param dm          The demographic model to which mutations should be added
#' @param lowerRange  If you want to estimate the mutation rate, this will be used 
#'            as the smallest possible value.
#' @param upperRange  If you want to estimate the mutation rate, this will be used 
#'            as the largest possible value.
#' @param fixedValue  If specified, the mutation rate will not be estimated, but assumend
#'            to be fixed at the given value.
#' @return        The demographic model with mutation
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.addMutation <- function(dm,lowerRange,upperRange,fixedValue){
    if ( missing(lowerRange) & missing(upperRange) & missing(fixedValue) ){
        dm@externalTheta <- T
        upperRange <- NA
        lowerRange <- NA
    }

    return(.addFeature(dm,"mutation","theta",lowerRange,upperRange,fixedValue))
}

#' Adds recombination events to a demographic model
#' 
#' @param dm          The demographic model to which recombination events should be added
#' @param lowerRange  If you want to estimate the recombination rate, this will be used 
#'            as the smallest possible value.
#' @param upperRange  If you want to estimate the recombination rate, this will be used 
#'            as the largest possible value.
#' @param fixedValue  If specified, the mutation rate will not be estimated, but assumend
#'            to be fixed at the given value.
#' @return        The demographic model with recombination
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addRecombination(dm,fixed=5)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.addRecombination <- function(dm,lowerRange,upperRange,fixedValue){
    return(.addFeature(dm,"recombination","rho",lowerRange,upperRange,fixedValue))
}

#' Adds symmetric migration to a demographic model
#' 
#' @param dm          The demographic model to which the migration should be added
#' @param lowerRange  If you want to estimate the mutation rate, this will be used 
#'            as the smallest possible value.
#' @param upperRange  If you want to estimate the mutation rate, this will be used 
#'            as the largest possible value.
#' @param fixedValue  If specified, the mutation rate will not be estimated, but assumend
#'            to be fixed at the given value.
#' @return        The demographic model with symmetric migration
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addSymmetricMigration(dm,0.1,5)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.addSymmetricMigration <- function(dm,lowerRange,upperRange,fixedValue){
    return(.addFeature(dm,"migration","m",lowerRange,upperRange,fixedValue))
}

dm.addEffectivePopSize <- function(dm,population,lowerRange,upperRange,fixedValue){
    return(.addFeature(dm,"effPopSize",paste("Ne",population,sep=""),lowerRange,upperRange,fixedValue))
}


#' Adds a speciation event to a demographic model
#' 
#' @param dm          The demographic model to which the speciation event should be added
#' @param lowerRange  If you want to estimate the mutation rate, this will be used 
#'            as the smallest possible value.
#' @param upperRange  If you want to estimate the mutation rate, this will be used 
#'            as the largest possible value.
#' @param fixedValue  If specified, the mutation rate will not be estimated, but assumend
#'            to be fixed at the given value.
#' @return        The extended demographic model
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.addSpeciationEvent <- function(dm,lowerRange,upperRange,fixedValue){
    return(.addFeature(dm,"split","tau",lowerRange,upperRange,fixedValue))
}


#' Adds a different present day population size of one population to a demographic model
#' 
#' @param dm          The demographic model to which the speciation event should be added
#' @param lowerRange  If you want to estimate the size parameter, this will be used 
#'            as the smallest possible value.
#' @param upperRange  If you want to estimate the size parameter, this will be used 
#'            as the largest possible value.
#' @param fixedValue  If specified, the size parameter will not be estimated, but assumend
#'            to be fixed at the given value.
#' @param population  The population of which the size should change. Can be 1 or 2 at the moment.
#' @return        The extended demographic model
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addPresentSize(dm,2,3,population=2)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.addPresentSize <- function(dm,lowerRange,upperRange,fixedValue,population){
    if (missing(population)) stop("No population given!")
    return(.addFeature(dm,"presentSize","q",lowerRange,upperRange,fixedValue,population=population))
}


#' Adds a variable population size of one population at split time to a demographic model
#' 
#' @param dm          The demographic model to which the speciation event should be added
#' @param lowerRange  If you want to estimate the size parameter, this will be used 
#'            as the smallest possible value.
#' @param upperRange  If you want to estimate the size parameter, this will be used 
#'            as the largest possible value.
#' @param fixedValue  If specified, the size parameter will not be estimated, but assumend
#'            to be fixed at the given value.
#' @param population  The population of which the size should change. Can be 1 or 2 at the moment.
#' @return        The extended demographic model
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addSplitSize(dm,.2,5,population=2)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.addSplitSize <- function(dm,lowerRange,upperRange,fixedValue,population){
    if (missing(population)) stop("No population given!")
    return(.addFeature(dm,"splitSize","s",lowerRange,upperRange,fixedValue,population=population))
}


dm.simulationCmd <- function(dm,parameters){
    .log3(3,"Called dm.simulationCmd()")
    .log3("Using simulation program",dm@simProg)
    if  (dm@simProg == "fsc") .fsc.generateCmd(dm,parameters)
    else if (dm@simProg == "ms" ) .ms.generateCmd(dm,parameters)
    else    message("ERROR: unkown simulation programm")
    .log3("Finished dm.simulationCmd()")
}

#' Simulates data according to a demographic model and calculates summary statistics form it
#' 
#' @param dm          The demographic model according to which the simulations should be done
#' @param parameters  A vector of parameters which should be used for the simulations. 
#'            If a matrix is given, a simulation for each row of the matrix will be performed
#' @param sumStatFunc The function to calculate summary statistics of the JSFS. It must return
#'            a numeric vector.
#' @return        A matrix where each row is the vector of summary statistics for 
#'            the parameters in the same row of the "parameter" matrix
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
#' dm.simSumStats(dm,c(1,10))
dm.simSumStats <- function(dm, parameters, sumStatFunc){
    .log3("Called dm.simSumStats()")

    if (!is.matrix(parameters)) parameters <- matrix(parameters,1,length(parameters))

    .check.dm(dm)
    if (dim(parameters)[2] != dm.getNPar(dm)) stop("Wrong number of parameters")
    if ( !.checkParInRange(dm,parameters) ) stop("Parameters out of range")

    simProg   <- dm@currentSimProg
    if (missing(sumStatFunc)) sumStatFunc <- simProg@sumStatFunc

	nSumStats <- length(sumStatFunc(dm,jsfs=matrix(0,dm@sampleSizes[1],dm@sampleSizes[2])))
    nSims	  <- max(dim(parameters)[1],1)
	sumStats  <- matrix(0,nSims,nSumStats)
   
	.log2("Simulating",nSumStats,"summary statistics for",nSims,"parameter combination(s)")
	
	wd <- getwd()
	setwd(tempdir())
    simProg@initialSeedFunc()

	for (n in 1:nSims) {
		suppressWarnings(
          simOutput <- system2(simProg@executable,
                        simProg@simParFunc(dm,parameters[n,]),
                        stdout=T))
        jsfs <- simProg@calcJSFSFunc(dm,simOutput)
		sumStats[n,] <- sumStatFunc(dm,jsfs,simOutput)
		.log2("SumStats:",sumStats[n,])
	}

	setwd(wd)

    .log3(dm,"Finished dm.simSumStats()")
    return(sumStats)
}
