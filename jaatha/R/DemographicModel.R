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
			simProg="character",
			debugMode="logical",
			logFile="character")
	)


#-----------------------------------------------------------------------
# Initialization
#-----------------------------------------------------------------------

.init <- function(.Object,sampleSizes,nLoci,seqLength,finiteSites,tsTvRatio,debugMode,logFile){
	.Object@features <- data.frame( type=character(),
				        parameter=numeric(),
					population=numeric(),
		      		    	lowerRange=numeric(),
		      		    	upperRange=numeric(),
		      		    	startsAtTime=numeric(),
				     	timeLine=numeric()
				      )

	.Object@parameters 	<- character()
	.Object@externalTheta 	<- F
	.Object@profiled 	<- F
	.Object@finiteSites 	<- finiteSites
	.Object@sampleSizes 	<- sampleSizes
	.Object@seqLength 	<- seqLength
	.Object@tsTvRatio 	<- tsTvRatio
	.Object@nLoci		<- nLoci
	.Object@simProg		<- .setSimProg(.Object) 

	.Object <- dm.setDebugMode(.Object,debugMode,logFile)

	return(.Object)
}

setMethod("initialize","DemographicModel",.init)
rm(.init)


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


.setSimProg <- function(dm){
	.check.dm(dm)
	if (dm@finiteSites) simProg <- "fsc"
	else 		    simProg <- "ms"
	return(simProg)
}

.makeThetaLast <- function(dm){
	.check.dm(dm)
	mut <- .getFeature(dm,"mutation",NA)
	if ( dim(mut)[1] == 0 ) return(dm)	#Yet no mutation added
	if ( dim(mut)[1] > 1  ) stop("Are there multiple mutation entries in your model?")
	mutPar <- mut$parameter
	if ( is.na(mutPar) ) return(dm)		#Looks like a fixed mutation rate
	nPar <- length(dm@parameters)
	if ( mutPar == nPar ) return(dm) 	#Already last parameter
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
# dm: 		the demographic model
# type: 	the type of the feature
# population: 	the population if the feature or NA
# 
# returns: 	0 			if feature is not in the model,
#		the fixed Value 	if the feature has a "fixedValue" or
#		the parameter name 	otherwise
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
#	feature <- .getFeature(dm,type,population)
#	return(feature@lowerRange == feature@upperRange)
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
#	.check.dm(dm)
#	if (.getParameter(dm,"effPopSize",NA) == 0)
#		dm <- dm.addAncPopSize(dm,NA,fixed=10000)
#	return(dm)
#}


.dm.log <- function(dm,...){
	.check.dm(dm)
	if (!dm@debugMode) return()
	cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),...,"\n",sep=" ",file=dm@logFile,append=T)
}

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

#-----------------------------------------------------------------------
# Public functions
#-----------------------------------------------------------------------
dm.setDebugMode <- function(dm,debugMode=T,logFile=""){
	.check.dm(dm)
	if (is.logical(debugMode)) dm@debugMode <- debugMode
	if (logFile != "") {
		dm@logFile <- paste(getwd(),"/",logFile,sep="")
		dm@debugMode <- T
	}
	return(dm)
}

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
#' 		      has multiple populations, this needs to be a vector
#'		      containing the sample sizes from each population.
#' @param nLoci	      Number of loci that will be simulated
#' @param seqLength   Number of bases for each locus
#' @param finiteSites If 'TRUE', a finite sites mutation model is assumed
#'		      instead of an infinite sites one.
#' @param tsTvRatio   Transition transversion ratio
#' @param debugMode   If 'TRUE', a debug output will be produced
#' @param logFile     If set, the debug output will be written into the given file
#' @return 	      The demographic model
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.createDemographicModel <- function(sampleSizes,nLoci,seqLength=1000,
				      finiteSites=F,tsTvRatio=.33,
				      debugMode=F,logFile=""){
	dm <- new("DemographicModel",sampleSizes,nLoci,seqLength,
		  finiteSites,tsTvRatio,debugMode,logFile)
	return(dm)
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
#' @param dm	      The demographic model to which mutations should be added
#' @param lowerRange  If you want to estimate the mutation rate, this will be used 
#'		      as the smallest possible value.
#' @param upperRange  If you want to estimate the mutation rate, this will be used 
#'		      as the largest possible value.
#' @param fixedValue  If specified, the mutation rate will not be estimated, but assumend
#'		      to be fixed at the given value.
#' @return 	      The demographic model with mutation
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
#' @param dm	      The demographic model to which recombination events should be added
#' @param lowerRange  If you want to estimate the recombination rate, this will be used 
#'		      as the smallest possible value.
#' @param upperRange  If you want to estimate the recombination rate, this will be used 
#'		      as the largest possible value.
#' @param fixedValue  If specified, the mutation rate will not be estimated, but assumend
#'		      to be fixed at the given value.
#' @return 	      The demographic model with recombination
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
#' @param dm	      The demographic model to which the migration should be added
#' @param lowerRange  If you want to estimate the mutation rate, this will be used 
#'		      as the smallest possible value.
#' @param upperRange  If you want to estimate the mutation rate, this will be used 
#'		      as the largest possible value.
#' @param fixedValue  If specified, the mutation rate will not be estimated, but assumend
#'		      to be fixed at the given value.
#' @return 	      The demographic model with symmetric migration
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
#' @param dm	      The demographic model to which the speciation event should be added
#' @param lowerRange  If you want to estimate the mutation rate, this will be used 
#'		      as the smallest possible value.
#' @param upperRange  If you want to estimate the mutation rate, this will be used 
#'		      as the largest possible value.
#' @param fixedValue  If specified, the mutation rate will not be estimated, but assumend
#'		      to be fixed at the given value.
#' @return 	      The extended demographic model
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
#' @param dm	      The demographic model to which the speciation event should be added
#' @param lowerRange  If you want to estimate the size parameter, this will be used 
#'		      as the smallest possible value.
#' @param upperRange  If you want to estimate the size parameter, this will be used 
#'		      as the largest possible value.
#' @param fixedValue  If specified, the size parameter will not be estimated, but assumend
#'		      to be fixed at the given value.
#' @param population  The population of which the size should change. Can be 1 or 2 at the moment.
#' @return 	      The extended demographic model
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
#' @param dm	      The demographic model to which the speciation event should be added
#' @param lowerRange  If you want to estimate the size parameter, this will be used 
#'		      as the smallest possible value.
#' @param upperRange  If you want to estimate the size parameter, this will be used 
#'		      as the largest possible value.
#' @param fixedValue  If specified, the size parameter will not be estimated, but assumend
#'		      to be fixed at the given value.
#' @param population  The population of which the size should change. Can be 1 or 2 at the moment.
#' @return 	      The extended demographic model
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
	.dm.log(dm,"Called dm.simulationCmd()")
	.dm.log(dm,"Using simulation program",dm@simProg)
	if 	(dm@simProg == "fsc") .fsc.generateCmd(dm,parameters)
	else if (dm@simProg == "ms" ) .ms.generateCmd(dm,parameters)
	else	message("ERROR: unkown simulation programm")
	.dm.log(dm,"Finished dm.simulationCmd()")
}

#' Simulates data according to a demographic model and calculates summary statistics form it
#' 
#' @param dm	      The demographic model according to which the simulations should be done
#' @param parameters  A vector of parameters which should be used for the simulations. 
#'		      If a matrix is given, a simulation for each row of the matrix will be performed
#' @param sumStatFunc The function to calculate summary statistics of the JSFS. It must return
#'		      a numeric vector.
#' @return 	      A matrix where each row is the vector of summary statistics for 
#'		      the parameters in the same row of the "parameter" matrix
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
#' dm.simSumStats(dm,c(1,10))
dm.simSumStats <- function(dm,parameters,sumStatFunc=dm.defaultSumStats){
	.dm.log(dm,"Called dm.simSumStats()")

	if (!is.matrix(parameters)) parameters <- matrix(parameters,1,length(parameters))

	.check.dm(dm)
	if (dim(parameters)[2] != dm.getNPar(dm)) stop("Wrong number of parameters")
	if ( !.checkParInRange(dm,parameters) ) stop("Parameters out of range")

	if 	(dm@simProg == "fsc") {
		sumStats <- .fsc.simSumStats(dm,parameters,sumStatFunc)
	}
	else if (dm@simProg == "ms" ) {
	       	sumStats <- .ms.simSumStats(dm,parameters,sumStatFunc)
	}
	else	message("ERROR: unkown simulation programm")
	.dm.log(dm,"Finished dm.simSumStats()")
	return(sumStats)
}

#' These are the default summary statistics for Jaatha
#' 
#' @param jsfs	      The joint site frequency spectrum of two populations
#' @return 	      A vector with sums over different areas of the JSFS
#' @export
#'
#' @examples
#' jsfs <- matrix(rpois(26*26,5),26,26)
#' dm.defaultSumStats(jsfs)
dm.defaultSumStats <- function(jsfs) {
  n <- nrow(jsfs)
  m <- ncol(jsfs)
  c(sum(jsfs[1,2:3]),
    sum(jsfs[2:3,1]),
    sum(jsfs[1,4:(m-3)]),
    sum(jsfs[4:(n-3),1]),
    sum(jsfs[1,(m-2):(m-1)]),
    sum(jsfs[(n-2):(n-1),1]),
    sum(jsfs[2:3,2:3]),
    sum(jsfs[2:3,4:(m-3)]),
    sum(jsfs[4:(n-3),2:3]),
    sum(jsfs[(n-2):(n-1),4:(m-3)]),
    sum(jsfs[4:(n-3),(m-2):(m-1)]),
    sum(jsfs[2:3,(m-2):(m-1)]),
    sum(jsfs[(n-2):(n-1),2:3]),
    sum(jsfs[4:(n-3),4:(m-3)]),
    sum(jsfs[(n-2):(n-1),(m-2):(m-1)]),
    jsfs[1,m],
    jsfs[n,1],
    sum(jsfs[n,2:3]),
    sum(jsfs[2:3,m]),
    sum(jsfs[n,4:(m-3)]),
    sum(jsfs[4:(n-3),m]),
    sum(jsfs[n,(m-2):(m-1)]),
    sum(jsfs[(n-2):(n-1),m]) )
}
