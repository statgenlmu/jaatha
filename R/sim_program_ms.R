# --------------------------------------------------------------
# Translates an demographic model to an ms command and 
# executes the simulation.
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Licence:  GPLv3 or later
# --------------------------------------------------------------

possible.features  <- c("sample", "loci.number", "loci.length",
                        "mutation", "migration", "split",
                        "recombination", "size.change", "growth",
                        "subgroups")
possible.sum.stats <- c("jsfs", "fpc", "trees", "seg.sites", "pmc", "file")


# This function generates an string that contains an R command for generating
# an ms call to the current model.
generateMsOptionsCommand <- function(dm) {
  nSample <- dm.getSampleSize(dm)
  cmd <- c('c(')
  cmd <- c(cmd,'"-I"', ",", length(nSample), ',', 
           paste(nSample, collapse=","), ',')

  for (i in 1:dim(dm@features)[1] ) {
    type <- as.character(dm@features[i,"type"])
    feat <- unlist(dm@features[i, ])

    if ( type == "mutation" ) {
      cmd <- c(cmd,'"-t"', ',', feat["parameter"], ',')
    }

    else if (type == "split") {
      cmd <- c(cmd, '"-ej"', ',', feat["time.point"], ',',
               feat["pop.sink"], ',', feat["pop.source"], ',')
    }

    else if (type == "migration")
      cmd <- c(cmd, '"-em"', ',', feat['time.point'], ',',
               feat['pop.sink'], ',', feat['pop.source']  , ',',
               feat['parameter'], ',')

    else if (type == "recombination") 
      cmd <- c(cmd, '"-r"', ',', feat['parameter'], ',', dm.getLociLength(dm), ',')

    else if (type == "size.change"){
      cmd <- c(cmd, '"-en"', ',', feat['time.point'], ',',
               feat["pop.source"], ',', feat['parameter'], ',')
    }

    else if (type == "growth"){
      cmd <- c(cmd, '"-eg"', ',' , feat["time.point"], ',',
               feat["pop.source"], ',', feat["parameter"], ',')
      }

    else if (type %in% c("sample", "loci.number", "loci.length", 
                         "pos.selection", "subgroups")) {}
    else stop("Unknown feature:", type)
  }

  if ('trees' %in% dm.getSummaryStatistics(dm)) cmd <- c(cmd, '"-T",')
  cmd <- c(cmd, '" ")')
}

generateMsOptions <- function(dm, parameters, subgroup) {
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

  cmd <- paste("ms", sum(dm.getSampleSize(dm)), dm.getLociNumber(dm), cmd)
  .print(cmd)
}

msSingleSimFunc <- function(dm, parameters) {
  checkType(dm, "dm")
  checkType(parameters, "num")
  if (length(parameters) != dm.getNPar(dm)) stop("Wrong number of parameters!")

  # Run a simulation for each subgroup
  subgroup_sizes <- sampleSubgroupSizes(dm)
  ms.files <- sapply(1:dm.getSubgroupNumber(dm), function(subgroup) {
    if (subgroup_sizes[subgroup] == 0) return("")
    ms.options <- generateMsOptions(dm, parameters, subgroup)
    ms.file <- getTempFile('ms')
    ms(sum(dm.getSampleSize(dm)), 
       subgroup_sizes[subgroup], 
       unlist(strsplit(ms.options, " ")), 
       ms.file)
    ms.file
  })
  
  # Parse the simulation output
  generateSumStats(ms.files, 0, parameters, dm)
}

finalizeMs <- function(dm) {
  dm@options[['ms.cmd']] <- generateMsOptionsCommand(dm)
  return(dm)
}

sampleSubgroupSizes <- function(dm) {
  subgroup_num <- dm.getSubgroupNumber(dm)
  if (subgroup_num == 1) return(dm.getLociNumber(dm))
  rmultinom(1, dm.getLociNumber(dm), rep(1,subgroup_num))[,1]
}

createSimProgram("ms", possible.features, possible.sum.stats, 
                 msSingleSimFunc, finalizeMs, printMsCommand, 100)