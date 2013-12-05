# --------------------------------------------------------------
# sim_prog_msms.R
# Adaptor to calling ms from a demographic model.
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Date:     2012-10-05
# Licence:  GPLv3 or later
# --------------------------------------------------------------

msms <- function(jar.path, ms.args, msms.args, out.file=NULL, jsfs=FALSE) {
  seed <- generateSeeds(1)
  if (jsfs) msms.args <- paste(msms.args, "-oAFS jAFS")
  command = paste("java -jar", jar.path, as.character(msms.args), 
                  "-ms", as.character(ms.args), "-seed", seed)
  print(command)
  if (!is.null(out.file)) {
    stopifnot(!file.exists(out.file))
    command <- paste(command, ">", out.file)
  }
  output <- system(command, intern=TRUE)

  if (!is.null(out.file)) {
    return(out.file)
  }

  if (jsfs) {
    jsfs.begin <- (1:length(output))[output == "Summary jAFS"]+2
    jsfs.end <- length(output)-1 
    jsfs <- t(sapply(jsfs.begin:jsfs.end, function(x)
                     as.integer(unlist(strsplit(output[x], " ")))))
    return(jsfs)
  }

  return(output)
}

possible.features  <- c("mutation","migration","split",
                        "recombination","size.change","growth", "pos.selection")
possible.sum.stats <- c("jsfs")

# This function generates an string that contains an R command for generating
# an ms call to the current model.
generateMsmsOptionsCommand <- function(dm) {
  nSample <- dm@sampleSizes
  cmd <- c('c(')
  cmd <- c(cmd,'"-I"', ",", length(nSample), ',', 
           paste(nSample, collapse=","), ',')

  for (i in 1:dim(dm@features)[1] ) {
    type <- as.character(dm@features[i,"type"])
    feat <- unlist(dm@features[i, ])

    if (type == "pos.selection") {
      cmd <- c(cmd,'"-t"', ',', feat["parameter"], ',')
    }
  }
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

  ms.options <- paste(sum(dm@sampleSizes), dm@nLoci, generateMsOptions(dm, parameters))
  msms.options <- generateMsmsOptions(dm, parameters) 
  sim.time <- system.time(jsfs <- msms("/home/paul/bin/msms.jar", 
                                          ms.options, msms.optios,
                                          false, true))
  .log3("finished after", sum(sim.time[-3]), "seconds")
  return(list(jsfs=jsfs, pars=parameters))
}

finalizeMs <- function(dm) {
  dm@options[['ms.cmd']] <- generateMsOptionsCommand(dm)
  return(dm)
}

createSimProgram("ms", "",
                 possible.features,
                 possible.sum.stats,
                 singleSimFunc=msmsSingleSimFunc,
                 finalizationFunc=finalizeMs)
