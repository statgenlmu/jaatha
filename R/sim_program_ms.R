# --------------------------------------------------------------
# sim_prog_ms.R
# Adaptor to calling ms from a demographic model.
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Date:     2012-10-05
# Licence:  GPLv3 or later
# --------------------------------------------------------------

possible.features  <- c("mutation","migration","split",
                        "recombination","size.change","growth")
possible.sum.stats <- c("jsfs")

#' Function to perform simulation using ms 
#' 
#' @param opts The options to pass to ms. Must either be a character or character
#' vector.
#' @param dm The demographic model we are using
#' @return The file containing the output of ms
callMs <- function(opts, dm){
  if (missing(opts)) stop("No options given!")
  opts <- unlist(strsplit(opts, " "))
  .log3("Called callMs")
  .log3("Options:", opts)

  ms.file <- getTempFile("ms")

  .log3("Calling ms...")
  ms(sum(dm@sampleSizes), dm@nLoci, opts, ms.file)
  .log3("ms finished. Finished callMs()")
  return(ms.file)
}

# This function generates an string that contains an R command for generating
# an ms call to the current model.
generateMsOptionsCommand <- function(dm) {
  nSample <- dm@sampleSizes
  cmd <- c('c(')
  cmd <- c(cmd,'"-I"', ",", length(nSample), ',', 
           paste(nSample, collapse=","), ',')

  for (i in 1:dim(dm@features)[1] ) {
    type <- as.character(dm@features[i,"type"])
    feat <- unlist(dm@features[i, ])

    if (type == "mutation") {
      cmd <- c(cmd,'"-t"', ',', feat["parameter"], ',')
    }

    if (type == "split") {
      cmd <- c(cmd, '"-ej"', ',', feat["time.point"], ',',
               feat["pop.sink"], ',', feat["pop.source"], ',')
    }

    if (type == "migration")
      cmd <- c(cmd, '"-em"', ',', feat['time.point'], ',',
               feat['pop.sink'], ',', feat['pop.source']  , ',',
               feat['parameter'], ',')

    if (type == "recombination") 
      cmd <- c(cmd, '"-r"', ',', feat['parameter'], ',', dm@seqLength, ',')

    if (type == "size.change"){
      cmd <- c(cmd, '"-en"', ',', feat['time.point'], ',',
               feat["pop.source"], ',', feat['parameter'], ',')
    }

    if (type == "growth"){
      cmd <- c(cmd, '"-eg"', ',' , feat["time.point"], ',',
               feat["pop.source"], ',', feat["parameter"], ',')
    }
  }

  cmd <- c(cmd, '"-T")')
}

generateMsOptions <- function(dm, parameters) {
  .log3("Called .ms.generateCmd()")
  ms.tmp <- new.env()

  par.names <- dm.getParameters(dm)
  for (i in seq(along = par.names)){
    ms.tmp[[ par.names[i] ]] <- parameters[i]
  }

  fixed.pars <- dm@parameters[dm@parameters$fixed, ]
  if (nrow(fixed.pars) > 0) {
    for (i in 1:nrow(fixed.pars)){
      ms.tmp[[ fixed.pars$name[i] ]] <- fixed.pars$lower.range[i]
    }
  }

  if ( !is.null( dm@options[['ms.cmd']] ) )
    cmd <- dm@options[['ms.cmd']]
  else
    cmd <- generateMsOptionsCommand(dm)
  cmd <- eval(parse(text=cmd), envir=ms.tmp)

  .log3("Finished .ms.generateCmd()")

  return(cmd)
}

printMsCommand <- function(dm) {
  cmd <- generateMsOptionsCommand(dm)

  cmd <- cmd[cmd != ","]
  cmd <- cmd[-c(1, length(cmd))]

  cmd <- paste(cmd, collapse=" ")

  cmd <- gsub(",", " ", cmd)
  cmd <- gsub('\"', "", cmd)
  cmd <- gsub('"', " ", cmd)

  return(cmd)
}

msOut2Jsfs <- function(dm, ms.out) {
  .log3("Called .ms.getJSFS()")
  jsfs <- matrix(.Call("msFile2jsfs", ms.out, dm@sampleSizes[1], 
                       dm@sampleSizes[2]),
                 dm@sampleSizes[1] + 1 ,
                 dm@sampleSizes[2] + 1,
                 byrow=T)
  .log3("Finished .ms.getJSFS()")
  return(jsfs)
}

msSingleSimFunc <- function(dm, parameters) {
  .log3("Called msSingleSimFunc()")
  .log3("parameter:",parameters)
  checkType(dm, "dm")
  checkType(parameters, "num")
  if (length(parameters) != dm.getNPar(dm)) 
    stop("Wrong number of parameters!")

  .log2("running ms")
  .log3("executing: \'", printMsCommand(dm), "'")
  ms.options <- generateMsOptions(dm, parameters)
  sim.time <- system.time(ms.out  <- callMs(ms.options, dm))
  .log3("finished after", sum(sim.time[-3]), "seconds")
  .log3("Simulation output in file", ms.out)

  .log2("calculating jsfs")
  jsfs  <- msOut2Jsfs(dm, ms.out)
  .log3("done. Removing tmp files...")
  unlink(ms.out)
  return(list(jsfs=jsfs, pars=parameters))
}


finalizeMs <- function(dm) {
  dm@options[['ms.cmd']] <- generateMsOptionsCommand(dm)
  return(dm)
}

createSimProgram("ms", "",
                 possible.features,
                 possible.sum.stats,
                 singleSimFunc=msSingleSimFunc,
                 finalizationFunc=finalizeMs)

