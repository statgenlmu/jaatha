setClass("DemographicModel" ,
         representation(features="data.frame",
                        parameters="data.frame",
                        sampleSizes="numeric",
                        nLoci="numeric",
                        seqLength="numeric",
                        tsTvRatio="numeric",
                        externalTheta="logical",
                        finiteSites="logical",
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



#------------------------------------------------------------------------------
# Private functions
#------------------------------------------------------------------------------

#---------------------------------------------------------------------
# dm.addParameter()
#---------------------------------------------------------------------
#' Create a parameter that can be used for one or more features
#'
#' The function creates a new model parameter. It can either have a fixed value
#' or you can enter a range of possible values if you want to estimate it.
#'
#' @param dm    The demographic model to which the parameter will be added
#' @param par.name The name of the parameter. You can use this name later to
#'                 access the parameter from a model feature
#' @param lower.boundary The lower boundary of the range within which the
#'                       parameter value will be estimated. Don't specify
#'                       'fixed.value' if you want to do so.
#' @param upper.boundary Like 'lower.boundary', but the upper end of the
#'                       parameter range.
#' @param fixed.value    If this argument is given, than rather than being estimated 
#'                       a fixed value will be used.
#' @return The original model extended with the new parameter.
#' @export
#' @examples
#' dm <- dm.createThetaTauModel(c(15,23), 100)
#' dm <- dm.addParameter(dm, "mig", 0.1, 5)
#' dm <- dm.addMigration(dm, par.new=FALSE, parameter="mig", pop.from=1, pop.to=2)
#' dm <- dm.addMigration(dm, par.new=FALSE, parameter="2*mig", pop.from=2, pop.to=1)
dm.addParameter <- function(dm, par.name, lower.boundary, upper.boundary, fixed.value) {

  if (missing(lower.boundary)) lower.boundary <- NA
  if (missing(upper.boundary)) upper.boundary <- NA
  if (missing(fixed.value)) fixed.value <- NA

  checkType(dm,          c("dm",   "s"), T, F)
  checkType(par.name,    c("char", "s"), T, F)
  checkType(lower.boundary, c("num",  "s"), T, T)
  checkType(upper.boundary, c("num",  "s"), T, T)
  checkType(fixed.value, c("num",  "s"), T, T)

  if (par.name %in% dm.getParameters(dm, T))
    stop("There is already a parameter with name ",par.name)

  if ( !is.na(fixed.value) ) {
    lower.boundary <- fixed.value
    upper.boundary <- fixed.value
  }

  fixed <- F
  if (is.na(lower.boundary) | is.na(upper.boundary)) {
    fixed <- T
  } else {
    if (lower.boundary == upper.boundary) {
      fixed <- T
    }
  }
 
  dm <- appendToParameters(dm, par.name, fixed, lower.boundary, upper.boundary)

  dm <- makeThetaLast(dm)
  return(dm)
}

#' Low level function for adding a new feature to a demographic Model
#'
#' Use this function to add a feature to a dm if there is no higher level
#' "dm.add*"-function availible.
#'
#' @param dm       The demographic model to which the feature will be added 
#' @param type     The type of the feature coded as a character
#' @param lower.range The lower boundary for the value of the parameter
#' @param upper.range The upper boundary for the value of the parameter
#' @param fixed.value If given, the parameter will be set to a fixed value. This
#'                 is equivalent to seting lower.range equal to upper.range.
#' @param pop.source The source population if availible (think e.g. of migration)
#' @param pop.sink   The target or "sink" population if availible (think e.g. of migration)
#' @param time.point Normally the point in backwards time where the feature
#'                   starts.
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'lower.range' and 'upper.range' or 
#'            'fixed.value'. It will be named 'parameter'.
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param parameter Either the name of the parameter (par.new=TRUE), or an R expression
#'            possibly containing one or more previously created parameter
#'            names.
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


# Gets the availible populations
getPopulations <- function(dm){
  # Every population other than the ancestral (=0) should be listed
  # the pop.sink field of its speciation event.
  populations <- c(1, dm@features$pop.sink)
  populations <- unique(populations[!is.na(populations)])
}


# To use Watersons estimator, Jaatha requires that the mutation parameter
# is the last parameter. This swishes it to the end.
makeThetaLast <- function(dm) {
  checkType(dm, c("dm", "s"), T, F)
  
  par.name <- dm@features$parameter[dm@features$type == "mutation"]
  if (length(par.name) == 0) return(dm)
  if (length(par.name) > 1) stop("Multiple Mutations features found")

  pars <- dm@parameters
  mut.line <- (1:nrow(pars))[pars$name == par.name]
  last.line <- nrow(pars)
  if (mut.line == last.line) return(dm)

  new.order <- 1:last.line
  new.order[c(mut.line, last.line)] <- c(last.line, mut.line)

  dm@parameters <- pars[new.order, ]
  return(dm)
}


# Checks if a vector of parameters is within the ranges of the model
.checkParInRange <- function(dm, param) {
  ranges <- dm.getParRanges(dm,inklExtTheta=F)
  #Seems there can be rounding errors during scaling
  lower <- t(ranges[, rep("lower.range", nrow(param))]) - 1e-5
  upper <- t(ranges[, rep("upper.range", nrow(param))]) + 1e-5
  inRange <- lower <= param & param <= upper
  all.in.range <- all(inRange)
  if (!all.in.range) {
    inRange <- inRange[, 1] & inRange[, 2]
    .print("The following parameter combinations are out of range:")
    out.of.range <- param[!inRange, ]
    for (i in 1:nrow(out.of.range)) {
      .print(out.of.range[i, ])
    }
  }
  return(all.in.range)
}

# Selects a program for simulation that is capable of all current features
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





#------------------------------------------------------------------------------
# Getters & Setters for Jaatha
#------------------------------------------------------------------------------
dm.getParameters <- function(dm, fixed=F) {
  checkType(dm,"dm")
  param <- dm@parameters$name
  param <- param[dm@parameters$fixed == fixed]
  return(param)
}

dm.getNPar <- function(dm){
  checkType(dm, "dm")
  return(length(dm.getParameters(dm)))
}

dm.getParRanges <- function(dm, inklExtTheta=T){
  #checkType(dm,"dm")
  par.ranges <- dm@parameters[!dm@parameters$fixed, c("lower.range","upper.range")]
  rownames(par.ranges) <- dm@parameters[!dm@parameters$fixed,"name"]
  if (!inklExtTheta) par.ranges <- par.ranges[1:dm.getNPar(dm),]
  return(par.ranges)
}

dm.setExternalTheta <- function(dm){
  checkType(dm, "dm")
  dm@externalTheta = T
  dm <- fixTheta(dm)
  return(dm)
}

fixTheta <- function(dm) {
  param <- dm@parameters
  param[nrow(param), "fixed"] <- T
  dm@parameters <- param
  return(dm)
}

getThetaName <- function(dm){
  return(dm@parameters[nrow(dm@parameters), 'name'])
}

getThetaRange <- function(dm){
  return(dm@parameters[nrow(dm@parameters),c('lower.range', 'upper.range')])
}





#------------------------------------------------------------------------------
# Creation new models
#------------------------------------------------------------------------------

#---------------------------------------------------------------------
# dm.createDemographicModel()
#---------------------------------------------------------------------
#' Create a basic demographic model
#' 
#' This function creates a basic empty demographic model, which
#' is returned. Features like mutation, pop.source splits and 
#' migration can be added afterwards.
#'
#' @param sample.sizes Number of individuals that are sampled. If your model 
#'            consists of multiple populations, this needs to be a vector
#'            containing the sample sizes from each population.
#' @param loci.num     Number of loci that will be simulated
#' @param seq.length   (Average) number of bases for each locus
# @param finiteSites If 'TRUE', a finite sites mutation model is assumed
#            instead of an infinite sites one.
# @param tsTvRatio   Transition transversion ratio
#' @param log.level   An integer specifing the amount of user readable output 
#'                    that will be produced.
#'                    0 = no output, 1 = normal verbosity, ..., 
#'                    3 = maximal verbosity.
#' @param log.file    If set, the debug output will be written into the given file
#' @return            The demographic model
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(sample.sizes=c(25,25), loci.num=100)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
#' dm
dm.createDemographicModel <- function(sample.sizes, loci.num, seq.length=1000, 
                                      #finiteSites=F, tsTvRatio=.33, 
                                      log.level, log.file) {
  setLogging(log.level, log.file)
  dm <- new("DemographicModel", sample.sizes, loci.num, seq.length, F, .33)
  return(dm)
}


#---------------------------------------------------------------------
# dm.createCustomModel()
#---------------------------------------------------------------------
# Creates a demographic model by giving the command line parameters 
# of a simulation program.
# 
# In case your model can not be described in terms of the demographic model
# class or you want to use a simulation program that is not yet supported, you
# can use this function to create a so called custom model. To do so, you have
# to state number (par.number) and ranges (par.ranges) of the parameters you
# want to estimate as well that the simulation program you want to use (sim.exe).
# Then, you have to create a function that generates the command-line options
# for the program (sim.options). Finally, you need a function that reads the output
# of your program and calculates the summary statistics (sumstats).
# 
# @param par.number   The number of parameters you want to estimate
# @param par.names    A character vector with names for your parameters
# @param par.ranges   The upper and lower boundaries for the values if your parameters.
#                     This should be a matrix with two columns for the lower and
#                     upper boundary respectively and a row for each parameter.
# @param sim.exe      The path of the executable of your simulation program. 
#                     E.g. "/usr/local/bin/ms" or "C:/msdir/ms.exe"
# @param sim.options  
# @param sumstats
# @return             A demographic model you can use for jaatha
# @export
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




#---------------------------------------------------------------------------------
# Front end functions for adding features
#---------------------------------------------------------------------------------

#--------------------------------------------------------------------
# dm.addMutation()
#--------------------------------------------------------------------
#' Adds mutations to a demographic model
#' 
#' This functions adds the assumption to the model that neutral mutations
#' occur in the genomes at a constant rate. The rate is quantified through
#' a parameter usually named theta in population genetics. It equals 4*Ne*mu,
#' where Ne is the (effective) number of diploid individuals in the ancestral
#' population and mu is the neutral mutation rate for an entire locus.
#'
#' @param dm  The demographic model to which mutations should be added
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'lower.range' and 'upper.range' or 
#'            'fixed.value'. It will be named 'new.par.name'.
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param lower.range  If you want to estimate the mutation rate, this 
#'            will be used as the smallest possible value.
#' @param upper.range  If you want to estimate the mutation rate, this 
#'            will be used as the largest possible value.
#' @param fixed.value  If specified, the mutation rate will not be 
#'            estimated, but assumed to be fixed at the given value.
#' @param new.par.name The name for the new parameter.
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "theta" will use
#'            an parameter with name theta that you have previously 
#'            created. You can also use R expression here, so "2*theta"
#'            or "5*theta+2*tau" (if tau is another parameter) will also
#'            work (also it does not make much sense).
#' @return    The demographic model with mutation.
#' @export
#'
#' @examples
#' # Create a new parameter
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addSpeciationEvent(dm, 0.01, 5)
#' dm <- dm.addMutation(dm, 1, 20)
#'
#' # Create a new fixed parameter
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addSpeciationEvent(dm, 0.01,5)
#' dm <- dm.addMutation(dm, fixed.value=7)
#'
#' # Use an existing parameter
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addParameter(dm, "theta", 0.01, 5)
#' dm <- dm.addMutation(dm, par.new=FALSE, parameter="2*log(theta)+1")
dm.addMutation <- function(dm, lower.range, upper.range, fixed.value,
                           par.new=T, new.par.name="theta", parameter) {

  if ( missing(lower.range) & missing(upper.range) & missing(fixed.value) ){
    dm@externalTheta <- T
    upper.range <- 0
    lower.range <- 0
  }

  if (par.new) parameter <- new.par.name
  dm <- addFeature(dm, "mutation", parameter, lower.range, upper.range,
                   fixed.value, par.new=par.new, time.point=NA)
  return(dm)
}


#-------------------------------------------------------------------
#  dm.addRecombination()
#-------------------------------------------------------------------
#' Adds recombination events to a demographic model
#'
#' This function add the assumption to the model that recombination
#' events may occur within each locus. The corresponding parameter
#' - usually name rho - equals 4*Ne*r, where Ne is the number of
#' diploid individuals in the ancestral population and r is the 
#' probability that a recombination event within the locus will
#' occur in one generation. Even when using an infinite sites
#' mutation model, this assumes an finite locus length which is given
#' by the 'seq.length' parameter of the demographic model.
#'
#' Please note that it does not make sense to estimate recombination
#' rates with Jaatha because it assumes unlinked loci.
#' 
#' @param dm  The demographic model to which recombination events should be added.
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'lower.range' and 'upper.range' or 
#'            'fixed.value'. It will be named 'new.par.name'.
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param lower.range  If you want to estimate the recombination rate (see note
#             above), this will be used as the smallest possible value.
#' @param upper.range  Same as lower.range, but the largest possible value.
#' @param fixed.value  If specified, the mutation rate will not be estimated,
#'                     but assumed to have the given value.
#' @param new.par.name The name for the new parameter.
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "rho" will use
#'            an parameter with name theta that you have previously 
#'            created. You can also use R expression here, so "2*rho"
#'            or "5*rho+2*tau" (if tau is another parameter) will also
#'            work (also it does not make much sense).
#' @return    The demographic model with recombination
#' @export
#'
#' @examples
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addSpeciationEvent(dm, 0.01, 5)
#' dm <- dm.addRecombination(dm, fixed=20)
#' dm <- dm.addMutation(dm, 1, 20)
dm.addRecombination <- function(dm, lower.range, upper.range, fixed.value,
                                par.new=T, new.par.name="rho", parameter) {

  if (par.new) parameter <- new.par.name
  dm <- addFeature(dm, "recombination", parameter, lower.range, upper.range,
                   fixed.value, par.new=par.new, time.point="0")
  return(dm)
}


#-------------------------------------------------------------------
#  dm.addMigration
#-------------------------------------------------------------------
#' Add migration/gene flow between two populations to a demographic model
#'
#' This function adds the assumption to the model that some individuals
#' 'migrate' from one sub-population to another, i.e. they leave the one
#' and become a member of the other. This is usually used to model ongoing
#' gene flow through hybridisation after the populations separated. 
#'
#' You can enter a time ('time.start') at which the migration is 
#' assumed to start (looking backwards in time). From that time on, a 
#' fixed number of migrants move from population 'pop.from' to
#' population 'pop.to' each generation. This number is given via this 
#' feature's parameter, which equals 4*Ne*m,  where Ne is the number of
#' diploid individuals in the ancestral population and m is the 
#' fraction of 'pop.to' that is replaced with migrants each generation. 
#' If 'pop.to' has also size Ne, than this is just the
#' expected number of individuals that migrate each generation.
#' 
#' You can add different mutation rates at different times to your model.
#' Then each rate will be used for the period from its time point to
#' the next. Migration from and to an population always ends with the 
#' speciation event in which the population is created.
#'
#' @param dm  The demographic model to which recombination events should be added.
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'lower.range' and 'upper.range' or 
#'            'fixed.value'. It will be named 'new.par.name'.
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param lower.range  If you want to estimate the migration parameter (see note
#             above), this will be used as the smallest possible value.
#' @param upper.range  Same as lower.range, but the largest possible value.
#' @param fixed.value  If specified, the migration rate will not be estimated,
#'                     but assumed to have the given value.
#' @param new.par.name The name for the new parameter.
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "M" will use
#'            an parameter with name M that you have previously 
#'            created. You can also use R expression here, so "2*M"
#'            or "5*M+2*tau" (if tau is another parameter) will also
#'            work (also this does not make much sense).
#' @param pop.from The population from which the individuals leave.
#' @param pop.to The population to which the individuals move.
#' @param time.start The time point at which the migration with this rate
#'            starts.
#' @return    The demographic model with migration
#' @export
#'
#' @examples
#' # Constant asymmetric migration
#' dm <- dm.createThetaTauModel(c(25,25), 100)
#' dm <- dm.addMigration(dm, 0.01, 5, pop.from=1, pop.to=2, time.start="0")
#'
#' # Stepwise decreasing mutation
#' dm <- dm.createThetaTauModel(c(25,25), 100)
#' dm <- dm.addMigration(dm, 0.01, 5, pop.from=1, pop.to=2, new.par.name="M", 
#'                       time.start="0")
#' dm <- dm.addMigration(dm, pop.from=1, pop.to=2, par.new=FALSE, 
#'                       parameter="0.5*M", time.start="0.5*tau")
dm.addMigration <- function(dm, lower.range, upper.range, fixed.value,
                            pop.from, pop.to, time.start="0",
                            par.new=T, new.par.name="M",
                            parameter) {

  checkType(pop.from, "num", T, F)
  checkType(pop.to,   "num", T, F)

  if (par.new) parameter <- new.par.name
  dm <- addFeature(dm, "migration", parameter, lower.range, upper.range,
                   fixed.value, par.new, pop.from, pop.to, time.start)
  return(dm)
}


#-------------------------------------------------------------------
#  dm.addSymmetricMigration
#-------------------------------------------------------------------
#' Adds symmetric migration between all populations
#'
#' This adds migration between all subpopulation, all with the same
#' rate. Please look at the documentation for \link{dm.addMigration} for 
#' detailed information about migration.
#' 
#' @param dm  The demographic model to which migration events should be added.
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'lower.range' and 'upper.range' or 
#'            'fixed.value'. It will be named 'new.par.name'.
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param lower.range  If you want to estimate the migration parameter (see note
#             above), this will be used as the smallest possible value.
#' @param upper.range  Same as lower.range, but the largest possible value.
#' @param fixed.value  If specified, the migration rate will not be estimated,
#'                     but assumed to have the given value.
#' @param new.par.name The name for the new parameter.
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "M" will use
#'            an parameter with name M that you have previously 
#'            created. You can also use R expression here, so "2*M"
#'            or "5*M+2*tau" (if tau is another parameter) will also
#'            work (also this does not make much sense).
#' @param time.start The time point at which the migration with this rate
#'            starts.
#' @return    The demographic model with migration
#' @export
#'
#' @examples
#' dm <- dm.createThetaTauModel(c(25,25), 100)
#' dm <- dm.addSymmetricMigration(dm, 0.01, 5, time.start="0.5*tau")
dm.addSymmetricMigration <- function(dm, lower.range, upper.range, fixed.value,
                                     par.new=T, new.par.name="M", parameter, 
                                     time.start="0") {

  if (par.new) dm <- dm.addParameter(dm, new.par.name, lower.range, upper.range, fixed.value)
  if (par.new) parameter <- new.par.name

  for (i in getPopulations(dm)) {
    for (j in getPopulations(dm)) {
      if (i==j) next
      dm <- dm.addMigration(dm, par.new=F, parameter=parameter, time.start=time.start,
                            pop.from=i, pop.to=j)
    }
  }

  return(dm)
}


#-------------------------------------------------------------------
# dm.addSpeciationEvent
#-------------------------------------------------------------------
#' Adds a speciation event to a demographic model
#'
#' You can use this function the create a new population that splits 
#' of from an existing population at a given time in the past. The time
#' can be given as parameter or as an expression based on previously
#' generated time points.
#'
#' As always, time in measured in Number of 4Ne generations in the past,
#' where Ne is the (effective) size of the ancestral population.
#'
#' The command will print the number of the new population, which will
#' be the number of previously existing populations plus one. The ancestral
#' population has number "1".
#'
#' @param dm  The demographic model to which the split should be added.
#' @param new.time.point  If 'TRUE' a new parameter will be created using the
#'            arguments 'min.time' and 'max.time' or 
#'            'fixed.time'. It will be named 'new.time.point.name'.
#'            If 'FALSE' the argument 'time.point' 
#'            will be evaluated instead.
#' @param min.time  If you want to estimate the time point, this will be 
#'            used as the smallest possible value.
#' @param max.time  Same as min.time, but the largest possible value.
#' @param fixed.time  If specified, the time.point will not be estimated,
#'            but assumed to have the given value.
#' @param new.time.point.name The name for the new time point.
#' @param in.population The number of the population in which the spilt
#'            occurs. See above for more information.
#' @param time.point  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "tau" will use
#'            an parameter with name tau that you have previously 
#'            created. You can also use R expression here, i.e. "2*tau"
#'            or "5*M+2*tau" (if M is another parameter) will also
#'            work (also this does not make much sense).
#' @return    The demographic model with a split.
#' @export
#' @examples
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
dm.addSpeciationEvent <- function(dm, min.time, max.time, fixed.time, 
                                  in.population=1, 
                                  new.time.point=T, new.time.point.name=NA,
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
  if (new.time.point) {
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


#-------------------------------------------------------------------
# dm.addSizeChange
#-------------------------------------------------------------------
#' Adds an instantaneous change of the population size of one 
#' population to a model.
#'
#' This function changes the effective population size of one 
#' population. The change is performed at a given time point
#' ('at.time') and applies to the time interval farther into 
#' the past from this point. The population size is set to a
#' factor of the size of the ancestral population Ne.
#'
#' If you want to add a slow, continuous change over some time,
#' then use the \link{dm.addGrowth} function.
#' 
#' @param dm  The demographic model to which the size change should be added.
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'min.size.factor' and 'max.size.factor' or 
#'            'fixed.size.factor'. It will be named 'new.time.point.name'
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param min.size.factor  If you want to estimate the size factor, this will be 
#'            used as the smallest possible value.
#' @param max.size.factor  Same as min.size.factor, but the largest possible value.
#' @param fixed.size.factor  If specified, the size factor not be estimated,
#'            but assumed to have the given value.
#' @param new.par.name Name for the new parameter.
#' @param population The number of the population in which the spilt
#'            occurs. See \link{dm.addSpeciationEvent} for more information.
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "tau" will use
#'            an parameter with name tau that you have previously 
#'            created. You can also use R expression here, i.e. "2*tau"
#'            or "5*M+2*tau" (if M is another parameter) will also
#'            work (also this does not make much sense).
#' @param at.time The time point at which the size changes.
#' @return    The demographic model with a size change.
#' @export
#' @examples
#' # A model with one smaller population
#' dm <- dm.createThetaTauModel(c(20,37), 88)
#' dm <- dm.addSizeChange(dm, 0.1, 1, population=2, at.time="0")
dm.addSizeChange <- function(dm, min.size.factor, max.size.factor,
                             fixed.size.factor, par.new=T, new.par.name="q",
                             parameter, 
                             population, at.time="0") {

  checkType(population, c("num",  "s"), T, F)
  checkType(at.time,    c("char", "s"), T, F)

  if (par.new) parameter <- new.par.name

  dm <- addFeature(dm, "size.change", parameter, min.size.factor, max.size.factor,
                   fixed.size.factor, par.new, population, NA, at.time)

  return(dm)
}


#-------------------------------------------------------------------
# dm.addGrowth
#-------------------------------------------------------------------
#' Adds an growth or decline of the population size of one 
#' population to a model.
#'
#' This function changes the growth factor of a population at given 
#' point in time ('at.time'). This factor than applies to the time 
#' interval farther into the past from this point. 
#'
#' The population size changes by factor exp(-alpha*t), where alpha
#' is the growth parameter and t is the time since the growth has 
#' started. Hence, for positive alpha, the population will 'decline 
#' backwards in time' or grow forwards in time. Similar, will decline
#' in forwards time for a negative value of alpha.
#'
#' If you want to add an instantaneous change of the population size,
#' then use the \link{dm.addSizeChange} function.
#' 
#' @param dm  The demographic model to which the size change should be added.
#' @param par.new  If 'TRUE' a new parameter will be created using the
#'            arguments 'min.growth.rate' and 'max.growth.rate' or 
#'            'fixed.growth.rate'. It will be named 'new.par.name'
#'            If 'FALSE' the argument 'parameter' 
#'            will be evaluated instead.
#' @param min.growth.rate  If you want to estimate the growth rate, this will be 
#'            used as the smallest possible value.
#' @param max.growth.rate  Same as min.growth.rate, but the largest possible value.
#' @param fixed.growth.rate  If specified, the growth rate will not be estimated,
#'            but assumed to have the given value.
#' @param new.par.name Name for the new parameter.
#' @param population The number of the population in which the spilt
#'            occurs. See \link{dm.addSpeciationEvent} for more information.
#' @param parameter  Instead of creating a new parameter, you can also
#'            set the mutation rate to an expression based on existing
#'            parameters. For example setting this to "alpha" will use
#'            an parameter with name tau that you have previously 
#'            created. You can also use R expression here, i.e. "2*alpha"
#'            or "5*M+2*alpha" (if M is another parameter) will also
#'            work (also the latter does not make much sense).
#' @param at.time The time point at which the size changes.
#' @return    The demographic model with a size change.
#' @export
#' @examples
#' # A model with one smaller population
#' dm <- dm.createThetaTauModel(c(20,37), 88)
#' dm <- dm.addGrowth(dm, 0.1, 2, population=2, at.time="0")
dm.addGrowth <- function(dm, min.growth.rate, max.growth.rate, fixed.growth.rate, 
                         par.new=T, new.par.name="alpha", parameter, 
                         population, at.time="0") {

  checkType(population, c("num",  "s"), T, F)
  checkType(at.time,    c("char", "s"), T, F)

  if (par.new) parameter <- new.par.name

  dm <- addFeature(dm, "growth", parameter, min.growth.rate, max.growth.rate,
                   fixed.growth.rate, par.new, population, NA, at.time)

  return(dm)
}


#-------------------------------------------------------------------
# dm.createThetaTauModel
#-------------------------------------------------------------------
#' Creates a standard "Theta/Tau" demopraphic model.
#' 
#' @param sample.sizes   A numeric vector with the sample sizes of the two pop.sources.
#' @param loci.num       The number of Loci to simulate.
#' @param seq.length     For recombination, each locus is assumed to be of this
#'                       length
#' @return               A Theta/Tau Model
#' @export
#'
#' @examples
#' dm <- dm.createThetaTauModel(c(20,25), 100)
#' dm
dm.createThetaTauModel <- function(sample.sizes, loci.num, seq.length=1000) {
  dm <- dm.createDemographicModel(sample.sizes, loci.num, seq.length)
  dm <- dm.addSpeciationEvent(dm, 0.01, 5)
  dm <- dm.addRecombination(dm, fixed.value=20)
  dm <- dm.addMutation(dm, 1, 20)
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
#' dm <- dm.createDemographicModel(c(25,25), 100)
#' dm <- dm.addSpeciationEvent(dm,0.01,5)
#' dm <- dm.addMutation(dm,1,20)
#' dm.simSumStats(dm,c(1,10))
dm.simSumStats <- function(dm, parameters, sumStatFunc){
  .log3("Called dm.simSumStats()")

  jsfs <- F
  if (!is.matrix(parameters)) parameters <- matrix(parameters,1,length(parameters))
  if (missing(sumStatFunc)) {
    if (nrow(parameters) > 1) stop("Only one parameter combination is allow,
                                   when not providing sumStatFunc.")
    jsfs <- T
    sumStatFunc  <- function(dm, jsfs){ return(as.vector(jsfs)) }
  }

  checkType(dm, "dm")
  if (dim(parameters)[2] != dm.getNPar(dm)) stop("Wrong number of parameters")
  if ( !.checkParInRange(dm,parameters) ) stop("Parameters out of range")

  simProg   <- dm@currentSimProg

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
              .log3("Simulating for pars:", parameter.combination)
              sim.out <- simProg@singleSimFunc(dm, parameter.combination)
              sumStats <- sumStatFunc(dm, sim.out)
              .log3("SumStats:", sumStats)
              return(sumStats)
            })

    
    sumStats <- t(sumStats)
    if (jsfs) sumStats <- matrix(sumStats, dm@sampleSizes[1]+1,
                                 dm@sampleSizes[2]+1)
  }
  
  .log3("Finished dm.simSumStats()")
  return(sumStats)
}
