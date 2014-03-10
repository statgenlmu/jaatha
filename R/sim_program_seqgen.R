# --------------------------------------------------------------
# sim_prog_seqgen.R
# Adaptor to calling ms from a demographic model.
# 
# Authors:  Paul R. Staab
# Date:     2013-11-21
# Licence:  GPLv3 or later
# --------------------------------------------------------------

# list ms's features + FS related features
seqgen.features    <- c('mutation.model', 'tstv.ratio', 
                        'base.freq.A', 'base.freq.C', 'base.freq.G',
                        'base.freq.T',
                        'gtr.rate.1', 'gtr.rate.2', 'gtr.rate.3',
                        'gtr.rate.4','gtr.rate.5','gtr.rate.6',
                        'gamma.categories', 'gamma.rate')

possible.sum.stats <- unique(c("jsfs", "file"), getSimProgram('ms')@possible.sum.stats)
mutation.models    <- c('HKY', 'F84', 'GTR')

possible.features  <- c(getSimProgram('ms')@possible.features, seqgen.features)

checkForSeqgen <- function() {
  if ( isJaathaVariable('seqgen.exe') ) return()

  # Works on Linux only maybe
  run.path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
  executables <- c(paste(run.path, "/seq-gen", sep=""), 
                   paste(run.path, "/seqgen", sep=""))
  for (exe in executables) {
    if (file.exists(exe)) {
      .print("Using", exe, "as seqgen executable\n")
      setJaathaVariable('seqgen.exe', exe)     
      return()
    }
  }

  stop("No seqgen executable found. Please provide one using
       Jaatha.setSeqgenExecutable()")
}

generateMsModel <- function(dm) {
  ms <- getSimProgram('ms')
  dm@features <- dm@features[dm@features$type %in% ms@possible.features, ]
  dm@sum.stats <- dm@sum.stats[dm@sum.stats %in% ms@possible.sum.stats]

  if (!"trees" %in% dm@sum.stats) dm@sum.stats <- append(dm@sum.stats, "trees") 
  if (!"file" %in% dm@sum.stats) dm@sum.stats <- append(dm@sum.stats, "file") 
  if ("jsfs" %in% dm@sum.stats) dm@sum.stats <- dm@sum.stats[dm@sum.stats != 'jsfs']
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
  .log3("executing: '", cmd, "'")
  system(cmd)

  if( !file.exists(seqgen.file) ) stop("seq-gen simulation failed!")
  if( file.info(seqgen.file)$size == 0 ) stop("seq-gen output is empty!")
  return(seqgen.file)
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
      model <- mutation.models[dm@parameters[dm@parameters$name == "mutation.model", 
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
  cmd <- generateSeqgenOptionsCmd(dm)

  cmd <- cmd[cmd != ","]
  cmd <- cmd[-c(1, length(cmd))]

  cmd <- paste(cmd, collapse=" ")

  cmd <- gsub(",", " ", cmd)
  cmd <- gsub('\"', "", cmd)
  cmd <- gsub('"', " ", cmd)

  return(cmd)
}

seqgenOut2Jsfs <- function(dm, seqgen.file) {
  if( ! file.exists(seqgen.file) ) stop("seq-gen simulation failed!")
  if (file.info(seqgen.file)$size == 0) stop("seq-gen output is empty!")

  sample.size <- dm.getSampleSize(dm)
  jsfs <- matrix(.Call("seqgen2jsfs", seqgen.file, sample.size[1], 
                       sample.size[2], dm.getLociNumber(dm)),
                 sample.size[1] + 1 ,
                 sample.size[2] + 1,
                 byrow=T)

  return(jsfs)
}

seqgenSingleSimFunc <- function(dm, parameters) {
  checkType(dm, "dm")
  checkType(parameters, "num")
  checkForSeqgen()
  if (length(parameters) != dm.getNPar(dm)) 
    stop("Wrong number of parameters!")

  sum.stats <- msSingleSimFunc(dm@options[['ms.model']], parameters)

  seqgen.options <- generateSeqgenOptions(dm, parameters)
  seqgen.file <- callSeqgen(seqgen.options, sum.stats[['file']])

  sum.stats[['pars']] <- parameters

  if ('jsfs' %in% dm@sum.stats) {
    sum.stats[['jsfs']] <- seqgenOut2Jsfs(dm, seqgen.file)
  }
  
  if ('file' %in% dm@sum.stats) {
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
  dm@options[['ms.model']] <- finalizeMs(generateMsModel(dm))
  dm@options[['seqgen.cmd']] <- generateSeqgenOptionsCmd(dm)
  stopifnot(!is.null(dm@options[['ms.model']]))
  return(dm)
}

createSimProgram("seq-gen", "",
                 possible.features,
                 possible.sum.stats,
                 singleSimFunc=seqgenSingleSimFunc,
                 finalizationFunc=finalizeSeqgen)
