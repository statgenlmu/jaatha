# --------------------------------------------------------------
# Translates an demographic model to an ms command and 
# executes the simulation.
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Licence:  GPLv3 or later
# --------------------------------------------------------------

possible.features  <- c("sample", "loci.number", "loci.length",
                        "mutation", "migration", "split",
                        "recombination", "size.change", "growth")
possible.sum.stats <- c("jsfs", "fpc", "tree", "seg.sites", "file")

#' Function to perform simulation using ms 
#' 
#' @param opts The options to pass to ms. Must either be a character or character
#' vector.
#' @param dm The demographic model we are using
#' @return The file containing the output of ms
callMs <- function(opts, dm){
  if (missing(opts)) stop("No options given!")
  opts <- unlist(strsplit(opts, " "))

  ms.file <- getTempFile("ms")

  ms(sum(dm.getSampleSize(dm)), dm.getLociNumber(dm), opts, ms.file)
  return(ms.file)
}

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

    if (type == "mutation") {
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

    else if (type %in% c("sample", "loci.number", "loci.length")) {}
    else stop("Unknown feature:", type)
  }

  if ('trees' %in% dm@sum.stats) cmd <- c(cmd, '"-T",')
  cmd <- c(cmd, '" ")')
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

msSingleSimFunc <- function(dm, parameters) {
  checkType(dm, "dm")
  checkType(parameters, "num")

  if (length(parameters) != dm.getNPar(dm)) stop("Wrong number of parameters!")

  ms.options <- generateMsOptions(dm, parameters)
  sim.time <- system.time(ms.out <- callMs(ms.options, dm))

  breaks.near <- dm@options[['fpc.breaks.near']]
  breaks.far <- dm@options[['fpc.breaks.near']]
  if (is.null(breaks.near))  breaks.near <- c(.25, .5, .75)
  if (is.null(breaks.far))  breaks.far <- c(.25, .5, .75)

  sum.stats <- parseOutput(ms.out, dm.getSampleSize(dm), dm.getLociNumber(dm), 0, 
                           'jsfs' %in% dm@sum.stats, 'seg.sites' %in% dm@sum.stats,
                           'fpc' %in% dm@sum.stats, breaks.near, breaks.far)

  sum.stats[['pars']] <- parameters
  if ("file" %in% dm@sum.stats) {
    sum.stats[['file']] <- ms.out
  }

  if (!'file' %in% dm@sum.stats) unlink(ms.out)
  return(sum.stats)
}

finalizeMs <- function(dm) {
  dm@options[['ms.cmd']] <- generateMsOptionsCommand(dm)
  return(dm)
}


# calcPercentFpcViolations <- function(snp.matrix) {
#   snp.matrix <- snp.matrix[, colSums(snp.matrix)>1, drop=FALSE]
#   if (ncol(snp.matrix) <= 1) return(c(near=NaN, far=NaN, theta=0))
#   snp.state <- apply(combn(1:ncol(snp.matrix), 2), 2, violatesFpc, snp.matrix)
#   return(c(near=sum(snp.state[2, snp.state[1, ]])/sum(snp.state[1, ]),
#            far=sum(snp.state[2, !snp.state[1, ]])/sum(!snp.state[1, ]),
#            theta=ncol(snp.matrix)/sum(1/1:(nrow(snp.matrix)-1)) ))
# }
# 
# violatesFpc <- function(sites, snp.matrix, near=.1) {
#   is.near <- diff(as.numeric(colnames(snp.matrix)[sites])) < near
#   status <- snp.matrix[ ,sites[1]] * 2 + snp.matrix[ ,sites[2]] 
#   if (all(0:3 %in% status)) return(c(near=is.near, violates=TRUE))
#   return(c(near=is.near, violates=FALSE))
# }
# 

createSimProgram("ms", "",
                 possible.features,
                 possible.sum.stats,
                 singleSimFunc=msSingleSimFunc,
                 finalizationFunc=finalizeMs)
