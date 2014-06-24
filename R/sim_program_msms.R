# --------------------------------------------------------------
# sim_prog_msms.R
# Adaptor to calling ms from a demographic model.
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Date:     2012-10-05
# Licence:  GPLv3 or later
# --------------------------------------------------------------

callMsms <- function(jar.path, ms.args, msms.args) {
  out.file = getTempFile("msms")
  seed <- generateSeeds(1)

  # Create the command
  command = paste("java -jar", jar.path, as.character(msms.args), 
                  "-ms", as.character(ms.args), "-seed", seed,
                  ">", out.file)

  # Execute the command
  output <- system(command)

  if(!file.exists(out.file)) stop("msms simulation failed!")
  if(file.info(out.file)$size == 0) stop("msms output is empty!")
  
  return(out.file)
}

checkForMsms <- function(throw.error = TRUE) {
  if (isJaathaVariable('msms.jar')) {
    if (file.exists(getJaathaVariable('msms.jar'))) {  
      return(TRUE)
    }
  }

  # Works on Linux only maybe
  run.path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
  executables <- paste(run.path, "/msms.jar", sep="")
  for (exe in executables) {
    if (file.exists(exe)) {
      .print("Using", exe, "as msms implementation\n")
      setJaathaVariable('msms.jar', exe)     
      return(TRUE)
    }
  }

  if (throw.error) stop("No msms executable found.")
  return(FALSE)
}

possible.features  <- c(getSimProgram('ms')@possible.features, "pos.selection")
possible.sum.stats <- getSimProgram('ms')@possible.sum.stats

# This function generates an string that contains an R command for generating
# an msms call to the current model.
generateMsmsOptionsCommand <- function(dm) {
  cmd <- c('c(')

  for (i in 1:dim(dm@features)[1] ) {
    type <- as.character(dm@features[i,"type"])
    feat <- unlist(dm@features[i, ])

    if (type == "pos.selection") {
      cmd <- c(cmd, '"-SI"', ',', feat['time.point'], ',', length(dm.getSampleSize(dm)), ',')
      if (feat['pop.source'] == 1) { 
        cmd <- c(cmd, 0.0005, ',', 0, ',') 
        cmd <- c(cmd, '"-Sc"', ',', 0, ',', 1, ',', 0, ',', 0, ',', 0, ',') 
      }
      else {
        cmd <- c(cmd, 0, ',', 0.0005, ',') 
        cmd <- c(cmd, '"-Sc"', ',', 0, ',', 0, ',', 0, ',', 0, ',', 0, ',') 
      }
      cmd <- c(cmd, '"-N 100000"', ',') 
      cmd <- c(cmd, '"-SAA"', ',', paste0("2*", feat['parameter']), ',',  '"-SAa"', ',',
               feat['parameter'], ',') 
      cmd <- c(cmd, '"-Sp 0.5"', ',', '"-SForceKeep"', ',')
    }
        }

  cmd <- c(cmd, '" ")')
  cmd
}

generateMsmsOptions <- function(dm, parameters) {
  msms.tmp <- new.env()

  par.names <- dm.getParameters(dm)
  for (i in seq(along = par.names)){
    msms.tmp[[ par.names[i] ]] <- parameters[i]
  }

  fixed.pars <- dm@parameters[dm@parameters$fixed, ]
  if (nrow(fixed.pars) > 0) {
    for (i in 1:nrow(fixed.pars)){
      msms.tmp[[ fixed.pars$name[i] ]] <- fixed.pars$lower.range[i]
    }
  }

  if ( !is.null( dm@options[['msms.cmd']] ) )
    cmd <- dm@options[['msms.cmd']]
  else
    cmd <- generateMsmsOptionsCommand(dm)
  cmd <- eval(parse(text=cmd), envir=msms.tmp)

  return(cmd)
}

msmsSingleSimFunc <- function(dm, parameters) {
  checkType(dm, "dm")
  checkType(parameters, "num")
  if (length(parameters) != dm.getNPar(dm)) 
    stop("Wrong number of parameters!")

  # Generate Options
  ms.options <- paste(sum(dm.getSampleSize(dm)), dm.getLociNumber(dm), 
                      paste(generateMsOptions(dm, parameters), collapse=" "))
  msms.options <- paste(generateMsmsOptions(dm, parameters), collapse= " ") 

  # Simulate
  out.file <- callMsms(getJaathaVariable('msms.jar'), ms.options, msms.options)

  # Parse Output
  sum.stats <- parseMsOutput(out.file, parameters, dm)

  return(sum.stats)
}

finalizeMsms <- function(dm) {
  checkForMsms()
  dm@options[['ms.cmd']] <- generateMsOptionsCommand(dm)
  dm@options[['msms.cmd']] <- generateMsmsOptionsCommand(dm)
  return(dm)
}

createSimProgram("msms", "",
                 possible.features,
                 possible.sum.stats,
                 singleSimFunc=msmsSingleSimFunc,
                 finalizationFunc=finalizeMsms)
