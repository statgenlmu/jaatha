# --------------------------------------------------------------
# sim_prog_seqgen.R
# Adaptor to calling ms from a demographic model.
# 
# Authors:  Paul R. Staab
# Date:     2013-11-21
# Licence:  GPLv3 or later
# --------------------------------------------------------------

# list ms's features + FS related features
sg.features    <- unique(c(getSimProgram('ms')$possible_features,
                           getSimProgram('msms')$possible_features,
                           'mutation.model', 'tstv.ratio', 
                           'base.freq.A', 'base.freq.C', 'base.freq.G', 
                           'base.freq.T',
                           'gtr.rate.1', 'gtr.rate.2', 'gtr.rate.3',
                           'gtr.rate.4','gtr.rate.5','gtr.rate.6',
                           'gamma.categories', 'gamma.rate'))

sg.sum.stats <- c('jsfs', 'file', 'seg.sites', 'fpc')
sg.mutation.models <- c('HKY', 'F84', 'GTR')

checkForSeqgen <- function(throw.error = TRUE) {
  if ( isJaathaVariable('seqgen.exe') ) {
    if (file.exists(getJaathaVariable('seqgen.exe'))) {  
      return(TRUE)
    }
  }

  # Works on Linux only maybe
  run.path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
  executables <- c(paste(run.path, "/seq-gen", sep=""), 
                   paste(run.path, "/seqgen", sep=""))
  for (exe in executables) {
    if (file.exists(exe)) {
      .print("Using", exe, "as seqgen executable\n")
      setJaathaVariable('seqgen.exe', exe)     
      return(TRUE)
    }
  }

  if (throw.error) {
    stop("No seqgen executable found. Please provide one using
          Jaatha.setSeqgenExecutable()")
  }
  return(FALSE)
}

generateTreeModel <- function(dm) {
  stopifnot(all(dm.getGroups(dm) == 1))
  if (any(msms.features %in% dm@features$type)) {
    tree.prog <- getSimProgram('msms')
  } else {
    tree.prog <- getSimProgram('ms')
  }
  
  dm@features <- dm@features[dm@features$type %in% tree.prog$possible_features, ]
  dm@sum.stats <- data.frame(name=c(), group=c())
  dm <- dm.addSummaryStatistic(dm, "seqgen.trees")
  dm <- dm.addSummaryStatistic(dm, "file")
  return(dm)
}


#' Set the path to the executable for seqgen
#'
#' @param seqgen.exe Path to seqgen's executable.
#' @export
Jaatha.setSeqgenExecutable <- function(seqgen.exe) {
  if (file.exists(seqgen.exe)) {
    setJaathaVariable('seqgen.exe', seqgen.exe)     
    .print("Using", seqgen.exe, "as seqgen executable\n")
  } else {
    stop("File", seqgen.exe, "does not exist")
  }
}

# Function to perform simulation using seqgen
# 
# @param opts The options to pass to ms. Must either be a character or character
# vector.
callSeqgen <- function(opts, ms.file) {
  if (missing(opts)) stop("No options given!")
  opts[length(opts) + 1:2] <- c("<", ms.file)
  opts <- paste(opts, collapse=" ")

  if( !file.exists(ms.file) ) stop("ms file not found")
  if( file.info(ms.file)$size == 0 ) stop("ms output is empty")

  seqgen.file <- getTempFile("seqgen")
  cmd <- paste(opts, ">", seqgen.file, sep=" ", collapse=" ")
  
  # Do the acctual simulation
  if (system(cmd, intern = FALSE) != 0) stop("seq-gen simulation failed")
  if( !file.exists(seqgen.file) ) stop("seq-gen simulation failed!")
  if( file.info(seqgen.file)$size == 0 ) stop("seq-gen output is empty!")
  
  seqgen.file
}

generateSeqgenOptions <- function(dm, parameters) {
  seqgen.tmp <- new.env()

  par.names <- dm.getParameters(dm)

  for (i in seq(along = par.names)){
    seqgen.tmp[[ par.names[i] ]] <- parameters[i]
  }

  fixed.pars <- dm@parameters[dm@parameters$fixed, ]
  if (nrow(fixed.pars) > 0) {
    for (i in 1:nrow(fixed.pars)){
      seqgen.tmp[[ fixed.pars$name[i] ]] <- fixed.pars$lower.range[i]
    }
  }

  seqgen.tmp[['seed']] <- generateSeeds(1)
  
  if ( !is.null( dm@options[['seqgen.cmd']] ) )
    cmd <- dm@options[['seqgen.cmd']]
  else
    cmd <- generateSeqgenOptionsCmd(dm)
  cmd <- paste(eval(parse(text=cmd), envir=seqgen.tmp), collapse=" ")

  return(cmd)
}

generateSeqgenOptionsCmd <- function(dm) {  
  base.freqs <- F
  gtr.rates <- F
  includes.model <- F

  opts <- c('c(', paste('"', getJaathaVariable('seqgen.exe'), '"', sep=""), ",")

  for (i in 1:dim(dm@features)[1] ) {
    type <- as.character(dm@features[i,"type"])
    feat <- unlist(dm@features[i, ])

    if (type == "mutation.model") {
      includes.model <- T
      model <- sg.mutation.models[dm@parameters[dm@parameters$name == "mutation.model", 
                                                'lower.range']]
      opts <- c(opts, paste('"-m', model, '"', sep=""), ",")
    }

    else if ( type %in% c('base.freq.A', 'base.freq.C', 
                          'base.freq.G', 'base.freq.T') )
      base.freqs <- T

    else if ( type %in% c('gtr.rate.1', 'gtr.rate.2', 'gtr.rate.3',
                          'gtr.rate.4', 'gtr.rate.5', 'gtr.rate.6') )
      gtr.rates <- T

    else if (type == "tstv.ratio")
      opts <- c(opts, '"-t"', ',', 'tstv.ratio', ',')

    else if (type == "gamma.rate")
      opts <- c(opts, '"-a"', ',', feat['parameter'], ',')

    else if (type == "gamma.categories")
      opts <- c(opts, '"-g"', ',', feat['parameter'], ',')
  }

  if (base.freqs) {
    opts <- c(opts, '"-f"', ',', 'base.freq.A',
              ',', 'base.freq.C',
              ',', 'base.freq.G',  
              ',', 'base.freq.T', ',')
  }

  if (gtr.rates) {
    opts <- c(opts, '"-r"', ',', 'gtr.rate.1',
              ',', 'gtr.rate.2',
              ',', 'gtr.rate.3',  
              ',', 'gtr.rate.4',  
              ',', 'gtr.rate.5',  
              ',', 'gtr.rate.6', ',')
  }

  if (!includes.model) {
    stop("You must specify a finite sites mutation model for this demographic model")
  }

  loci.length <- dm.getLociLength(dm)
  opts <- c(opts, '"-l"', ',', loci.length, ',')
  opts <- c(opts, '"-s"', ',', paste(getThetaName(dm), "/", loci.length), ',')
  opts <- c(opts, '"-p"', ',', loci.length + 1, ',')
  opts <- c(opts, '"-z"', ',', 'seed', ',')
  opts <- c(opts, '"-q"', ')')
  return(opts)
}

printSeqgenCommand <- function(dm) {
  if (is.null(dm@options[['tree.model']])) dm <- dm.finalize(dm)
  tree.model <- dm@options[['tree.model']]
  getSimProgram(tree.model@currentSimProg)$print_cmd_func(tree.model)
  cmd <- generateSeqgenOptionsCmd(dm)

  cmd <- cmd[cmd != ","]
  cmd <- cmd[-c(1, length(cmd))]

  cmd <- paste(cmd, collapse=" ")

  cmd <- gsub(",", " ", cmd)
  cmd <- gsub('\"', "", cmd)
  cmd <- gsub('"', " ", cmd)

  .print(cmd)
}

seqgenSingleSimFunc <- function(dm, parameters) {
  checkType(dm, "dm")
  checkType(parameters, "num")
  checkForSeqgen()
  if (length(parameters) != dm.getNPar(dm)) 
    stop("Wrong number of parameters!")

  # Use ms to simulate the ARG
  tree.model <- dm@options[['tree.model']]
  if (is.null(tree.model)) tree.model <- generateTreeModel(dm)
  sum.stats <- dm.simSumStats(tree.model, parameters)

  # Call seq-gen to distribute mutations
  seqgen.options <- generateSeqgenOptions(dm, parameters)
  seqgen.file <- callSeqgen(seqgen.options, sum.stats[['file']])

  sum.stats[['pars']] <- parameters
  
  # Parse the output & generate additional summary statistics
  if ('fpc' %in% dm.getSummaryStatistics(dm)) {
    breaks.near <- dm@options[['fpc.breaks.near']]
    breaks.far <- dm@options[['fpc.breaks.far']]
    stopifnot(!is.null(breaks.near))
    stopifnot(!is.null(breaks.far))
    
    sum.stats <- parseOutput(seqgen.file, dm.getSampleSize(dm), dm.getLociNumber(dm), 1, 
                             'jsfs' %in% dm.getSummaryStatistics(dm), 
                             'seg.sites' %in% dm.getSummaryStatistics(dm),
                             TRUE, breaks.near, breaks.far)
  } else {
    sum.stats <- parseOutput(seqgen.file, dm.getSampleSize(dm), dm.getLociNumber(dm), 1, 
                             'jsfs' %in% dm.getSummaryStatistics(dm), 
                             'seg.sites' %in% dm.getSummaryStatistics(dm),
                             FALSE)
  }
  
  sum.stats[['pars']] <- parameters
  if ('file' %in% dm.getSummaryStatistics(dm)) {
    sum.stats[['file']] <- c(ms=sum.stats[['file']],
                             seqgen=seqgen.file)
  } else {
    unlink(sum.stats[['file']])
    unlink(seqgen.file)
    sum.stats[['file']] <- NULL
  }

  return(sum.stats)
}

finalizeSeqgen <- function(dm) {
  checkForSeqgen()
  dm@options[['tree.model']] <- dm.finalize(generateTreeModel(dm))
  dm@options[['seqgen.cmd']] <- generateSeqgenOptionsCmd(dm)
  stopifnot(!is.null(dm@options[['tree.model']]))
  return(dm)
}

createSimProgram("seq-gen", sg.features, sg.sum.stats,
                 seqgenSingleSimFunc, finalizeSeqgen, printSeqgenCommand,
                 priority=10)

rm(sg.features, sg.sum.stats, seqgenSingleSimFunc, 
   finalizeSeqgen, printSeqgenCommand)