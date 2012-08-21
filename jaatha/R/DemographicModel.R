setClass("DemographicModel" ,
         representation(features="data.frame",
                        parameters="data.frame",
                        sampleSizes="numeric",
                        nLoci="numeric",
                        seqLength="numeric",
                        tsTvRatio="numeric",
                        externalTheta="logical",
                        finiteSites="logical",
                        profiled="logical",
                        currentSimProg="SimProgram")
         )


#-----------------------------------------------------------------------
# Initialization
#-----------------------------------------------------------------------

.init <- function(.Object, sampleSizes, nLoci, seqLength,
                  finiteSites, tsTvRatio){

  .Object@features <- data.frame(type=character(),
                                 parameter=character(),
                                 pop.source=numeric(),
                                 pop.sink=numeric(),
                                 time.point=character(),
                                 group=numeric(),
                                 stringsAsFactors=F )

  .Object@parameters <- data.frame(parameter=character(),
                                   fixed=logical(),
                                   lower.range=numeric(),
                                   upper.range=numeric(),
                                   stringsAsFactors=F )


  .Object@externalTheta   <- F
  .Object@profiled        <- F
  .Object@finiteSites     <- finiteSites
  .Object@sampleSizes     <- sampleSizes
  .Object@seqLength       <- seqLength
  .Object@tsTvRatio       <- tsTvRatio
  .Object@nLoci           <- nLoci

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
    lR <- as.character(features[rw,'lower.range']) 
    uR <- as.character(features[rw,'upper.range'])
    parname <- features[rw,'parameter']

    if (type == "split")
      cat("-",parname,": A split into two pop.sources between",lR,"and",uR,
          "Ne generations ago.\n")
    else if (type == "mutation")
      cat("-",parname,": A scaled mutation rate between",lR,"and",uR,"\n")
    else if (type == "recombination")
      cat("-",parname,": A scaled recombination rate between",lR,"and",uR,"\n")
    else if (type == "migration")
      cat("-",parname,":  A scaled migration rate between",lR,"and",uR,"\n")
    else if (type == "presentSize")
      cat("-",parname,":  At sample time, the second pop.source two was between",
          lR,"and",uR,"times as the first one\n")
    else if (type == "splitSize") 
      cat("-",parname,":  At time of the pop.source split, the second",
          "pop.source two was between",lR,"and",uR,"times as the first one\n")
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
      lR <- as.character(features[rw,'lower.range']) 

      if (type == "split")
        cat("- A split into two pop.sources",lR,"Ne generations ago.\n")
      else if (type == "mutation")
        cat("- A fixed scaled mutation rate of",lR,"\n")
      else if (type == "recombination")
        cat("- A fixed scaled recombination rate of",lR,"\n")
      else if (type == "migration")
        cat("- A fixed scaled migration rate of",lR,"\n")
      else if (type == "presentSize")
        cat("- At sample time, the second pop.source two was",
            lR,"times as the first one\n")
      else if (type == "splitSize") 
        cat("- At time of the pop.source split, the second pop.source
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
# Private functions
#-----------------------------------------------------------------------

dm.addParameter <- function(dm, par.name, lower.range, upper.range, fixed.value) {

  if (missing(lower.range)) lower.range <- NA
  if (missing(upper.range)) upper.range <- NA
  if (missing(fixed.value)) fixed.value <- NA

  checkType(dm,          c("dm",   "s"), T, F)
  checkType(par.name,    c("char", "s"), T, F)
  checkType(lower.range, c("num",  "s"), T, T)
  checkType(upper.range, c("num",  "s"), T, T)
  checkType(fixed.value, c("num",  "s"), T, T)

  if (par.name %in% dm.getParameters(dm, T))
    stop("There is already a parameter with name ",par.name)

  if ( !is.na(fixed.value) ) {
    lower.range <- fixed.value
    upper.range <- fixed.value
  }

  fixed <- F
  if (is.na(lower.range) | is.na(upper.range)) {
    fixed <- T
  } else {
    if (lower.range == upper.range) {
      fixed <- T
    }
  }
 
  dm <- appendToParameters(dm, par.name, fixed, lower.range, upper.range)

  #dm <- makeThetaLast(dm)
  return(dm)
}

#' Low level function for adding a new feature to a demographic Model
#'
#' Use this function to add a feature to a dm if there is no higher level
#' "dm.add*"-function availible.
#'
#' @param dm       The demographic model to which the feature will be added 
#' @param type     The type of the feature coded as a character
#' @param par.name The name of the correspondig parameter (also a character)
#' @param lower.range The lower boundary for the value of the parameter
#' @param upper.range The upper boundary for the value of the parameter
#' @param fixed.value If given, the parameter will be set to a fixed value. This
#'                 is equivalent to seting lower.range equal to upper.range.
#' @param pop.source The source population if availible (think e.g. of migration)
#' @param pop.sink   The target or "sink" population if availible (think e.g. of migration)
#' @param time.point Normally the point in backwards time where the feature
#'                   starts. 
#' @param group     For genomic features, different groups can be created.
#' @return          The extended demographic model.
addFeature <- function(dm, type, parameter=NA,
                       lower.range=NA, upper.range=NA, fixed.value=NA, par.new=T,
                       pop.source=NA, pop.sink=NA,
                       time.point=NA, group=NA) {
  
  if (missing(par.new))     par.new <- T
  if (missing(parameter))   parameter <- NA
  if (missing(pop.source))  pop.source <- NA
  if (missing(pop.sink))    pop.sink <- NA
  if (missing(time.point))  time.point <- NA
  if (missing(group))       group <- NA

  # Check inputs
  checkType(dm,          c("dm",   "s"), T, F)
  checkType(type,        c("char", "s"), T, F)
  checkType(par.new,     c("bool", "s"), T, F)
  checkType(parameter,   c("char", "s"), T, T)
  checkType(pop.source,  c("num",  "s"), T, T)
  checkType(pop.sink,    c("num",  "s"), T, T)
  checkType(time.point,  c("char", "s"), T, T)
  checkType(group,       c("num",  "s"), T, T)

  if (par.new) dm <- dm.addParameter(dm, parameter, lower.range, upper.range, fixed.value)

  # Append the feature
  dm <- appendToFeatures(dm = dm,
                         type = type,
                         parameter = parameter,
                         pop.source = pop.source,
                         pop.sink = pop.sink,
                         time.point = time.point,
                         group = group)

  # Update some technical properties of the model
  dm <- .dm.selectSimProg(dm)

  return(dm)
}


# Helper function that appends a feature to the "feature" dataframe
# Does not check the feature for consistency
# This should only be used by addFeature().
appendToFeatures <- function(dm, type, parameter, pop.source, 
                             pop.sink, time.point, group    ) {
  
  new.feature <- data.frame(type=type,
                            parameter=parameter,
                            pop.source=pop.source,
                            pop.sink=pop.sink,
                            time.point=time.point,
                            group=group,
                            stringsAsFactors=F)

  dm@features <- rbind(dm@features, new.feature)
  return(dm)
}

# Helper function that appends a parameter to the "parameters" dataframe
# Does not check the feature for consistency
# This should only be used by addFeature().
appendToParameters <- function(dm, name, fixed, lower.range, upper.range) {

  new.parameter <- data.frame(name=name,
                              fixed=fixed,
                              lower.range=lower.range,
                              upper.range=upper.range,
                              stringsAsFactors=F)
  
  dm@parameters <- rbind(dm@parameters, new.parameter)
  return(dm)
}


# Function that returns all features with the given criteria in the model.
# Also returns a data.frame.
getFeature <- function(dm, type, pop.source=NA, 
                       pop.sink=NA, time.point=NA, group=NA ){

  checkType(dm, c("dm", "s"), T, F)
  
  feat <- dm@features
  good.rows <- feat$type == type

  if (is.na(pop.source)) good.rows <- good.rows & is.na(feat$pop.source)
  else good.rows <- good.rows & (feat$pop.source == pop.source)

  if (is.na(pop.sink)) good.rows <- good.rows & is.na(feat$pop.sink)
  else good.rows <- good.rows & (feat$pop.sink == pop.sink)

  if (is.na(time.point)) good.rows <- good.rows & is.na(feat$time.point)
  else good.rows <- good.rows & (feat$time.point == time.point)

  if (is.na(group)) good.rows <- good.rows & is.na(feat$group)
  else good.rows <- good.rows & (feat$group == group)

  return(feat[good.rows, ])
}

# Gets the availible populations
getPopulations <- function(dm){
  # Every population other than the ancestral (=0) should be listed
  # the pop.sink field of its speciation event.
  populations <- c(1, dm@features$pop.sink)
  populations <- unique(populations[!is.na(populations)])
}

# Adds an parameter to the parameter list
# This should only be used by addFeature().
addParameter <- function(dm, parameter.name){
  checkType(dm, c("dm", "s"), T, F)
  checkType(parameter.name, c("char", "s"), T, F)
  if (any(dm.getParameters(dm) == parameter.name)) 
    stop("parameter name already exists")
  dm.getParameters <- c(dm.getParameters, parameter.name)
  return(dm)
}

# To use Watersons estimator, Jaatha requires that the mutation parameter
# is the last parameter. This swishes it to the end.
makeThetaLast <- function(dm) {
  checkType(dm, c("dm", "s"), T, F)
  feat <- dm@features
  mut.line <- (1:nrow(feat))[feat$type == "mutation"]
  last.line <- nrow(feat)
  if (length(mut.line) != 1) return(dm)
  if (mut.line == last.line) return(dm)
  new.order <- 1:last.line
  new.order[c(mut.line, last.line)] <- c(last.line, mut.line)
  dm@features <- feat[new.order, ]
  return(dm)
}

# returns the parameter that belong to a feature specified via 
# type and pop.source
#
# @param dm           the demographic model
# @param type         the type of the feature
# @param pop.source   the pop.source if the feature or NA
# 
# @returns  0                   if feature is not in the model,
#           the fixed Value     if the feature has a "fixed.value" or
#           the parameter name  otherwise
getParameter <- function(dm, type, pop.source=NA, pop.sink=NA,
                         time.point=NA, group=NA) {

  checkType(dm, c("dm", "s"), T, F)
  feature <- getFeature(dm, type, pop.source, pop.sink,
                        time.point, group)

  if ( dim(feature)[1] == 0 ) return(0)
  if ( feature$lower.range[1] == feature$upper.range[1] )
    return(feature$lower.range[1])
  else
    return(dm.getParameters[feature$parameter[1]])
}

# .getNumOfHistoricalEvent <- function(dm){
#   .check.dm(dm)
#   return(sum(dm@features$type != "effPopSize" & dm@features$type != "mutation"))
# }


#.isFixedFeature <- function(dm,type,pop.source=NA){
#   feature <- getFeature(dm,type,pop.source)
#   return(feature@lower.range == feature@upper.range)
#}

# .checkFeatureParameters <- function(lower.range,upper.range,fixed.value,time.point,timeLine){
#   if (! .is.NumericOrFalse(lower.range) ) cat("Wrong value for 'lower.range'\n")
#   if (! .is.NumericOrFalse(upper.range) ) cat("Wrong value for 'upper.range'\n")
# if (! .is.NumericOrFalse(time.point) ) cat("Wrong value for 'time.point'\n")
# if (! .is.NumericOrFalse(timeLine) ) cat("Wrong value for 'timeLine'\n")
#   if (! .is.NumericOrFalse(fixed.value) ) cat("Wrong value for 'fixed.value'\n")
# 
# 
#   if ( fixed.value == F & ( lower.range == F | upper.range == F ) )
#     cat("Wrong parameters. Either both lower.range and upper.range 
#         or fixed.value must be specified\n")
# 
#         if ( fixed.value != F && lower.range != F && upper.range != F ) 
#           cat("Wrong parameters. Either enter a fixed value or a range for this parameter \n")
# }


# .parseFeatureParametes <- function(lower.range,upper.range,fixed.value,time.point,timeLine){
#   if ( time.point == F ) time.point <<- Inf
#   if ( timeLine == F ) timeLine <<- Inf
# 
#   if ( is.numeric(fixed.value) ) {
#     lower.range <<- fixed.value
#     upper.range <<- fixed.value
#   }
#   return(c(lower.range,upper.range,fixed.value,time.point,timeLine))
# }


# .is.NumericOrFalse <- function(x){
#   return ( (is.numeric(x) || x == F) && length(x) == 1 )
# }

# .calcSizeChange <- function(dm, param, percentages){
#   checkType(dm, "dm")
# 
#   pop <- dm@features[param,"pop.source"]
#   q <- .calcAbsParamValue(dm,percentages[param],param)
# 
# get value of "size" parameter s for the pop.source
#   sizeLine <- (1:dim(dm@features)[1])[dm@features$type=="size" & dm@features$pop.source==pop]
#   if (sum(sizeLine) == 1) { 
#     s <- .calcAbsParamValue(dm,percentages[sizeLine],sizeLine) 
#   } else { 
#     s <- 1 
#   }
# 
# get tau (XXX GLOBAL TAU ONLY ATM XXX)
#   tauLine <- (1:dim(dm@features)[1])[dm@features$type=="mutation" & is.na(dm@features$pop.source)]
#   tau <- .calcAbsParamValue(dm,percentages[tauLine],tauLine) 
# 
#   return(log(q/s)/tau)
# }


# .calcAbsParamValue <- function(dm,percentage,param){
# 
#   if (!is.matrix(percentage)) 
#     percentage <- matrix(percentage,1)
# 
#   return( t(apply(percentage,1,.calcAbsParamValueLine,dm=dm)) )
# }
# 
# .calcAbsParamValueLine <- function(dm,percentage,param){
#   features <- dm@features[!is.na(dm@features$parameter),]
# 
# If not specified which param to convert, use the first length(percentage) ones
#   if ( missing(param) ) param <- 1:length(percentage)
# 
# Calculate scaledValues
#   scaledValues <- exp ( percentage * 
#                        ( log(features[param,"upper.range"]) 
#                         - log(features[param,"lower.range"]) )) * 
#         features[param,"lower.range"]
# 
# Add names
#         names(scaledValues) <- dm.getParameters[param]
# 
# Do not scale theta if externalTheta = T
#          if ( dm@externalTheta ) {
# TRUE iff param[.] links to theta
#           theta.mask <- param == (1:length(dm.getParameters))[dm.getParameters == "theta"]
#           if (any(theta.mask)) scaledValues[theta.mask] <- percentage[theta.mask]
#         }
# 
#         return( scaledValues )
# }

.checkParInRange <- function(dm, param) {
  ranges <- dm.getParRanges(dm,inklExtTheta=F)
  #Seems there can be rounding errors during scalation
  lower <- matrix(ranges[,1],dim(param)[1],dim(param)[2],byrow=T) - 1e-15
  upper <- matrix(ranges[,2],dim(param)[1],dim(param)[2],byrow=T) + 1e-15
  inRange <- lower <= param & param <= upper
  return(all(inRange))
}


.dm.selectSimProg <- function(dm) {
  for (i in seq(along = .local$simProgs)){
    if (all(dm@features$type %in% .local$simProgs[[i]]@possible.features)) {
      dm@currentSimProg <- .local$simProgs[[i]]
      .log2("Using",dm@currentSimProg@name,"for simulations")
      return(dm)
    }
  }
  warning("No suitable simulation software found!")
  return(dm)
}

#-----------------------------------------------------------------------
# Public functions
#-----------------------------------------------------------------------
dm.getParameters <- function(dm, include.fixed=F) {
  param <- dm@features$parameter
  if (!include.fixed) param <- param[!dm@features$fixed]
  return(param[!is.na(param)])
}

dm.getParRanges <- function(dm,inklExtTheta=T){
  checkType(dm,"dm")
  parMask <- !is.na(dm@features$parameter)
  parRanges <- cbind(lower=dm@features$lower.range[parMask],upper=dm@features$upper.range[parMask])
  parRanges <- parRanges[sort.list(dm@features$parameter[parMask]),]
  if (!inklExtTheta) parRanges <- parRanges[1:dm.getNPar(dm),]
  return( parRanges )
}

#' Create a basic demographic model
#' 
#' This function creates a basic empty demographic model, which
#' is returned. Features like mutation, pop.source splits and 
#' migration can be added afterwards.
#'
#' @param sampleSizes Number of individual that are sampled. If your model 
#'            has multiple pop.sources, this needs to be a vector
#'            containing the sample sizes from each pop.source.
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
##' @export
dm.createCustomModel <- function(par.number, par.names, par.ranges, sim.exe, 
                                 sim.options, sumstats) {

  stop("Not implemented right now!")

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
  checkType(dm, "dm")
  #newDm <- new("DemographicModel",dm@seqLength,dm@sampleSizes,dm@finiteSites)
  #newDm@features <- dm@features[dm@features$type != "mutation",]
  #newDm@parameters <- dm.getParameters[dm.getParameters != "theta"]
  #newDm <- dm.addMutation(newDm)
  dm@externalTheta = T
  return(dm)
}


dm.getNPar <- function(dm){
  checkType(dm, "dm")
  return(length(dm.getParameters(dm))-dm@externalTheta)
}



#' Adds mutations to a demographic model
#' 
#' @param dm          The demographic model to which mutations should be added
#' @param lower.range  If you want to estimate the mutation rate, this will be used 
#'            as the smallest possible value.
#' @param upper.range  If you want to estimate the mutation rate, this will be used 
#'            as the largest possible value.
#' @param fixed.value  If specified, the mutation rate will not be estimated, but assumend
#'            to be fixed at the given value.
#' @return        The demographic model with mutation
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.addMutation <- function(dm,lower.range,upper.range,fixed.value){
  if ( missing(lower.range) & missing(upper.range) & missing(fixed.value) ){
    dm@externalTheta <- T
    upper.range <- 0
    lower.range <- 0
  }

  return(addFeature(dm,"mutation","theta",lower.range,upper.range,fixed.value))
}

#' Adds recombination events to a demographic model
#' 
#' @param dm          The demographic model to which recombination events should be added
#' @param lower.range  If you want to estimate the recombination rate, this will be used 
#'            as the smallest possible value.
#' @param upper.range  If you want to estimate the recombination rate, this will be used 
#'            as the largest possible value.
#' @param fixed.value  If specified, the mutation rate will not be estimated, but assumend
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
dm.addRecombination <- function(dm,lower.range,upper.range,fixed.value){
  return(addFeature(dm,"recombination","rho",lower.range,upper.range,fixed.value))
}

#' Adds symmetric migration to a demographic model
#' 
#' @param dm          The demographic model to which the migration should be added
#' @param lower.range  If you want to estimate the mutation rate, this will be used 
#'            as the smallest possible value.
#' @param upper.range  If you want to estimate the mutation rate, this will be used 
#'            as the largest possible value.
#' @param fixed.value  If specified, the mutation rate will not be estimated, but assumend
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
dm.addSymmetricMigration <- function(dm,lower.range,upper.range,fixed.value){
  return(addFeature(dm,"migration","m",lower.range,upper.range,fixed.value))
}

dm.addEffectivePopSize <- function(dm,pop.source,lower.range,upper.range,fixed.value){
  return(addFeature(dm,"effPopSize",paste("Ne",pop.source,sep=""),lower.range,upper.range,fixed.value))
}

#' Adds a speciation event to a demographic model
#' 
#' @param dm          The demographic model to which the speciation event should be added
#' @param lower.range  If you want to estimate the mutation rate, this will be used 
#'            as the smallest possible value.
#' @param upper.range  If you want to estimate the mutation rate, this will be used 
#'            as the largest possible value.
#' @param fixed.value  If specified, the mutation rate will not be estimated, but assumend
#'            to be fixed at the given value.
#' @return        The extended demographic model
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.addSpeciationEvent <- function(dm, min.time, max.time, fixed.time, 
                                  in.population=1, new.time.point.name=NA,
                                  time.point=NA){

  checkType(dm,                  c("dm",  "s"),  T, F)
  checkType(in.population,       c("num", "s"),  T, F)
  checkType(min.time,            c("num", "s"),  F, F)
  checkType(max.time,            c("num", "s"),  F, F)
  checkType(fixed.time,          c("num", "s"),  F, F)
  checkType(new.time.point.name, c("char","s"),  F, T)
  checkType(time.point,          c("num", "s"),  F, T)
  
  # in.population valid?
  populations <- getPopulations(dm)
  if (!in.population %in% populations) 
    stop("Population", in.population, "not known")

  # time.points
  if (is.na(time.point)) {
    if(is.na(new.time.point.name)) 
      new.time.point.name <- paste("t_split_", in.population, sep="")
    dm <- dm.addParameter(dm, new.time.point.name, min.time, 
                          max.time, fixed.time)
    time.point <- new.time.point.name
  } else {
    # Check if valid
  }

  new.pop <- max(populations) + 1
  .print("The new population is population", new.pop)

  return(addFeature(dm, "split", par.new = F, 
                    pop.source=in.population,
                    pop.sink=new.pop,
                    time.point=time.point))
}


#' Adds a different present day pop.source size of one pop.source to a demographic model
#' 
#' @param dm          The demographic model to which the speciation event should be added
#' @param lower.range  If you want to estimate the size parameter, this will be used 
#'            as the smallest possible value.
#' @param upper.range  If you want to estimate the size parameter, this will be used 
#'            as the largest possible value.
#' @param fixed.value  If specified, the size parameter will not be estimated, but assumend
#'            to be fixed at the given value.
#' @param pop.source  The pop.source of which the size should change. Can be 1 or 2 at the moment.
#' @return        The extended demographic model
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addPresentSize(dm,2,3,pop.source=2)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.addPresentSize <- function(dm,lower.range,upper.range,fixed.value,pop.source){
  if (missing(pop.source)) stop("No pop.source given!")
  return(addFeature(dm,"presentSize","q",lower.range,upper.range,fixed.value,pop.source=pop.source))
}


#' Adds a variable pop.source size of one pop.source at split time to a demographic model
#' 
#' @param dm          The demographic model to which the speciation event should be added
#' @param lower.range  If you want to estimate the size parameter, this will be used 
#'            as the smallest possible value.
#' @param upper.range  If you want to estimate the size parameter, this will be used 
#'            as the largest possible value.
#' @param fixed.value  If specified, the size parameter will not be estimated, but assumend
#'            to be fixed at the given value.
#' @param pop.source  The pop.source of which the size should change. Can be 1 or 2 at the moment.
#' @return        The extended demographic model
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addSplitSize(dm,.2,5,pop.source=2)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.addSplitSize <- function(dm,lower.range,upper.range,fixed.value,pop.source){
  if (missing(pop.source)) stop("No pop.source given!")
  return(addFeature(dm,"splitSize","s",lower.range,upper.range,fixed.value,pop.source=pop.source))
}

#' Creates a standard "Theta/Tau" demopraphic model.
#' 
#' @param sampleSizes   A numeric vector with the sampleSizes of the two pop.sources.
#' @param nLoci         The number of Loci to simulate.
#' @return              A Theta/Tau Model
#' @export
#'
#' @examples
#' dm.createThetaTauModel(c(20,25), 100)
dm.createThetaTauModel <- function(sampleSizes, nLoci) {
  dm <- dm.createDemographicModel(sampleSizes=sampleSizes, nLoci=nLoci, seqLength=1000)
  dm <- dm.addSpeciationEvent(dm,0.01,5)
  dm <- dm.addMutation(dm,1,20)
  return(dm)
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

  checkType(dm, "dm")
  if (dim(parameters)[2] != dm.getNPar(dm)) stop("Wrong number of parameters")
  if ( !.checkParInRange(dm,parameters) ) stop("Parameters out of range")

  simProg   <- dm@currentSimProg
  if (missing(sumStatFunc)) sumStatFunc <- simProg@defaultSumStatFunc

  # nSumStats <- length(sumStatFunc(dm,jsfs=matrix(0,dm@sampleSizes[1],dm@sampleSizes[2])))
  # nSims	  <- max(dim(parameters)[1],1)
  # sumStats  <- matrix(0,nSims,nSumStats)

  #.log2("Simulating",nSumStats,"summary statistics for",nSims,"parameter combination(s)")

  if (!simProg@useSingleSimFunc) {
    # Use the simFunc to simulate all parameter combinations at once.
    # Still needs some thought about how to combine with the
    # defaultSumStatFunc.
    sumStats <- 0

  } else {
    # Use the singleSimFunc to simulate each of the parameter combinations
    sumStats <- 
      apply(parameters, 1,
            function(parameter.combination){
              .log3("Simulating for pars:",parameter.combination)
              sim.out <- simProg@singleSimFunc(dm, parameter.combination)
              sumStats <- sumStatFunc(dm, sim.out)
              .log3("SumStats:", sumStats)
              return(sumStats)
            }) 

    sumStats <- t(sumStats)
  }

  .log3("Finished dm.simSumStats()")
  return(sumStats)
}
