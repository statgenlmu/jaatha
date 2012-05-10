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

.addFeature <- function(dm,type,parameter,population,
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
	if (is.na(population))
		return(dm@features[dm@features$type == type 
		       		   & is.na(dm@features$population),])
	else
		return(dm@features[dm@features$type == type 
		       		   & dm@features$population == population,])
}
	

.addParameter <- function(dm,parameterName){
	if (any(dm@parameters == parameterName)) 
		stop("parameter already exists")
	dm@parameters <- c(dm@parameters,parameterName)
	return(dm)
}


.setSimProg <- function(dm){
	if (dm@finiteSites) simProg <- "fsc"
	else 		    simProg <- "ms"
	return(simProg)
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
	feature <- .getFeature(dm,type,population)
	if ( dim(feature)[1] == 0 ) return(0)
	if ( feature$lowerRange[1] == feature$upperRange[1] )
		return(feature$lowerRange[1])
	else
		return(dm@parameters[feature$parameter[1]])
}

.getNumOfHistoricalEvent <- function(dm){
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

.calcSizeChange <- function(dm,param,percentages){
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

.addNeIfMissing <- function(dm){
	if (.getParameter(dm,"effPopSize",NA) == 0)
		dm <- dm.addAncPopSize(dm,NA,fixed=10000)
	return(dm)
}


.dm.log <- function(dm,...){
	if (!dm@debugMode) return()
	cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),...,"\n",sep=" ",file=dm@logFile,append=T)
}

#-----------------------------------------------------------------------
# Public functions
#-----------------------------------------------------------------------
dm.setDebugMode <- function(dm,debugMode=T,logFile=""){
	if (is.logical(debugMode)) dm@debugMode <- debugMode
	if (is.character(logFile)) dm@logFile <- logFile
	return(dm)
}

dm.getParRanges <- function(dm){
	parMask <- !is.na(dm@features$parameter)
	return( cbind(lower=dm@features$lowerRange[parMask],upper=dm@features$upperRange[parMask]) )
}

dm.createDemographicModel <- function(sampleSizes,nLoci,seqLength,
				      finiteSites=F,tsTvRatio=.33,
				      debugMode=F,logFile=""){
	dm <- new("DemographicModel",sampleSizes,nLoci,seqLength,finiteSites,tsTvRatio,debugMode,logFile)
	return(dm)
}

dm.setExternalTheta <- function(dm){
	 #newDm <- new("DemographicModel",dm@seqLength,dm@sampleSizes,dm@finiteSites)
	 #newDm@features <- dm@features[dm@features$type != "mutation",]
	 #newDm@parameters <- dm@parameters[dm@parameters != "theta"]
	 #newDm <- dm.addMutation(newDm)
	 dm@externalTheta = T
	 return(dm)
}


dm.getNPar <- function(dm){
	return(length(dm@parameters)-dm@externalTheta)
}

dm.addMutation <- function(dm,lowerRange=F,upperRange=F,fixedValue=F,
			   population=NA,startsAtTime=NA,timeLine=NA){

	#.checkFeatureParameters(lowerRange,upperRange,fixedValue,startsAtTime,timeLine)

	if ( (lowerRange==F) & (upperRange==F) & (fixedValue==F) ){
		dm@externalTheta <- T
		upperRange <- 5
		lowerRange <- 5
	}

	if ( is.numeric(fixedValue) ) {
		lowerRange <- fixedValue
		upperRange <- fixedValue
	}

	parameter <- NA

	if ( fixedValue == F) {
		dm <- .addParameter(dm,"theta")
		parameter <- length(dm@parameters)
	}

	dm <- .addFeature(dm,"mutation",parameter,population,
		             lowerRange,upperRange,startsAtTime,timeLine)
	
	return(dm)
}

dm.addRecombination <- function(dm,lowerRange=F,upperRange=F,fixedValue=F,
			   population=NA,startsAtTime=NA,timeLine=NA){

	#.checkFeatureParameters(lowerRange,upperRange,fixedValue,startsAtTime,timeLine)

	if ( (lowerRange==F) & (upperRange==F) & (fixedValue==F) ){
		dm@externalTheta <- T
		upperRange <- 5
		lowerRange <- 5
	}

	if ( is.numeric(fixedValue) ) {
		lowerRange <- fixedValue
		upperRange <- fixedValue
	}

	parameter <- NA

	if ( fixedValue == F) {
		dm <- .addParameter(dm,"rho")
		parameter <- length(dm@parameters)
	}

	dm <- .addFeature(dm,"recombination",parameter,population,
		             lowerRange,upperRange,startsAtTime,timeLine)
	
	return(dm)
}

dm.addSymmetricMigration <- function(dm,lowerRange=F,upperRange=F,
			   		startsAtTime=NA,timeLine=NA,fixedValue=F){

	if ( is.numeric(fixedValue) ) {
		lowerRange <- fixedValue
		upperRange <- fixedValue
	}

	parameter <- NA
	population <- NA

	if ( fixedValue == F) {
		dm <- .addParameter(dm,"mig")
		parameter <- length(dm@parameters)
	}

	dm <- .addFeature(dm,"migration",parameter,population,
		             lowerRange,upperRange,startsAtTime,timeLine)
	
	return(dm)
}


dm.addEffectivePopSize <- function(dm,population,lowerRange=F,upperRange=F,
			           startsAtTime=NA,timeLine=NA,fixedValue=F){

	if ( is.numeric(fixedValue) ) {
		lowerRange <- fixedValue
		upperRange <- fixedValue
	}

	parameter  <- NA

	if ( fixedValue == F) {
		dm <- .addParameter(dm,paste("Ne",population,sep=""))
		parameter <- length(dm@parameters)
	}

	dm <- .addFeature(dm,"effPopSize",parameter,population,
		             lowerRange,upperRange,startsAtTime,timeLine)

	return(dm)
}


dm.addSpeciationEvent <- function(dm,lowerRange=F,upperRange=F,
				  population=1,timeLine=NA,fixedValue=F){

	#.checkFeatureParameters(lowerRange,upperRange,fixedValue,startsAtTime=F,timeLine)

	if ( is.numeric(fixedValue) ) {
		lowerRange <- fixedValue
		upperRange <- fixedValue
	}

	parameter <- NA 
	startsAtTime <- NA

	if ( fixedValue == F) {
		dm <- .addParameter(dm,"tau")
		parameter <- length(dm@parameters)
	}

	dm <- .addFeature(dm,"split",parameter,population,
		             lowerRange,upperRange,startsAtTime,timeLine)

	return(dm)
}


dm.setPopSize <- function(dm,population=1,lowerRange=F,upperRange=F,fixedValue=F,
			  startsAtTime=NA,timeLine=NA){

	if ( is.numeric(fixedValue) ) {
		lowerRange <- fixedValue
		upperRange <- fixedValue
	}

	parameter <- NA 

	if ( fixedValue == F) {
		dm <- .addParameter(dm,paste("pop",population,"_size",sep=""))
		parameter <- length(dm@parameters)
	}

	dm <- .addFeature(dm,"popSize",parameter,population,
		             lowerRange,upperRange,startsAtTime,timeLine)

	return(dm)
}


dm.addPopSizeChange <- function(dm,population=1,lowerRange=F,upperRange=F,fixedValue=F,
			  startsAtTime=NA,timeLine=NA){

if ( is.numeric(fixedValue) ) {
		lowerRange <- fixedValue
		upperRange <- fixedValue
	}

	parameter <- NA 

	if ( fixedValue == F) {
		dm <- .addParameter(dm,paste("pop",population,"_sizeChange",sep=""))
		parameter <- length(dm@parameters)
	}

	dm <- .addFeature(dm,"sizeChange",parameter,population,
		             lowerRange,upperRange,startsAtTime,timeLine)

	return(dm)
}


dm.simulationCmd <- function(dm,parameters){
	.dm.log(dm,"Called dm.simulationCmd()")
	.dm.log(dm,"Using simulation program",dm@simProg)
	if 	(dm@simProg == "fsc") .fsc.generateCmd(dm,parameters)
	else if (dm@simProg == "ms" ) .ms.generateCmd(dm,parameters)
	else	message("ERROR: unkown simulation programm")
	.dm.log(dm,"Finished dm.simulationCmd()")
}


dm.getJSFS <- function(dm){
	if (dm@finiteSites) .fsc.getJSFS(dm)
}

dm.simSumStats <- function(dm,parameters){
	.dm.log(dm,"Called dm.simSumStats()")
	if (!is.matrix(parameters)) parameters <- matrix(parameters,1,length(parameters))

	if 	(dm@simProg == "fsc") {
		sumStats <- .fsc.simSumStats(dm,parameters)
	}
	else if (dm@simProg == "ms" ) {
	       	sumStats <- .ms.simSumStats(dm,parameters)
	}
	else	message("ERROR: unkown simulation programm")
	.dm.log(dm,"Finished dm.simSumStats()")
	return(sumStats)
}

#source("DemographicModel-ms.R")
#source("DemographicModel-fsc.R")

#-----------------------------------------------------------------------
# CleanUp
#-----------------------------------------------------------------------
