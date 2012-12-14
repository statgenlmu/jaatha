# --------------------------------------------------------------
# sim_prog_seqgen.R
# Adaptor to calling ms from a demographic model.
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Date:     2012-10-25
# Licence:  GPLv3 or later
# --------------------------------------------------------------

#test code:
#source("./aaa_SimProgram.R")
#source("./helper_functions.R")
#source("./sim_program_ms.R")
#source("./DemographicModel.R")

# list ms's features + FS related features
seqgen.features    <- c('mutation.model', 'tstv.ratio', 
                        'base.freq.A', 'base.freq.C', 'base.freq.G',
                        'base.freq.T',
                        'gtr.rate.1', 'gtr.rate.2', 'gtr.rate.3',
                        'gtr.rate.4','gtr.rate.5','gtr.rate.6')

possible.sum.stats <- c("jsfs")
mutation.models    <- c('HKY', 'F84', 'GTR')

possible.features  <- c(getSimProgram('ms')@possible.features, seqgen.features)

checkForSeqgen <- function() {
  if ( isJaathaVariable('seqgen.exe') ) return()

  # Works on linux
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
  opts[length(opts) + 1:2] <- c(" ", ms.file)
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

  cmd <- generateSeqgenOptionsCmd(dm)
  cmd <- paste(eval(parse(text=cmd), envir=seqgen.tmp), collapse=" ")

  return(cmd)
}

generateSeqgenOptionsCmd <- function(dm, parameters) {  
  base.freqs <- F
  gtr.rates <- F

  opts <- c('c(', paste('"', getJaathaVariable('seqgen.exe'), '"', sep=""), ",")

  for (i in 1:dim(dm@features)[1] ) {
	type <- as.character(dm@features[i,"type"])
    feat <- unlist(dm@features[i, ])
  
    if (type == "mutation.model") {
      model <- mutation.models[dm@parameters[dm@parameters$name == "mutation.model", 
                                             'lower.range']]
      opts <- c(opts, '"-m"', ',', 
                paste('"', model, '"', sep=""), ",")
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
  }

  if (base.freqs) {
    opts <- c(opts, '"-f"', ',', 'base.freq.A',
                            ',', 'base.freq.C',
                            ',', 'base.freq.G',  
                            ',', 'base.freq.T', ',')
  }

  if (gtr.rates) {
    opts <- c(opts, '"-t"', ',', 'base.freq.1',
                            ',', 'base.freq.2',
                            ',', 'base.freq.3',  
                            ',', 'base.freq.4',  
                            ',', 'base.freq.5',  
                            ',', 'base.freq.6', ',')
  }

  opts <- c(opts, '"-l"', ',', dm@seqLength, ',')
  opts <- c(opts, '"-p"', ',', dm@seqLength + 1, ',')
  opts <- c(opts, '"-z"', ',', 'seed', ',')
  opts <- c(opts, '"-q"', ')')
  return(opts)
}

printSeqgenCommand <- function(dm) {
  cmd <- generateSeqgenOptions(dm)

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

  jsfs.size <- (dm@sampleSizes[1]+1)*(dm@sampleSizes[2]+1)
  jsfs <- matrix( .C("seqFile2jsfs",
                     as.character(seqgen.file),
                     as.integer(dm@sampleSizes[1]),
                     as.integer(dm@sampleSizes[2]),
                     as.integer(dm@nLoci),
                     as.integer(jsfs.size),
                     res=integer(jsfs.size),
                     PACKAGE="jaatha")$res,
                 dm@sampleSizes[1] + 1 ,
                 dm@sampleSizes[2] + 1,
                 byrow=T)

  return(jsfs)
}

seqgenSingleSimFunc <- function(dm, parameters) {
  .log3("called msSingleSimFunc()")
  .log3("parameter:",parameters)
  checkType(dm, "dm")
  checkType(parameters, "num")
  checkForSeqgen()
  if (length(parameters) != dm.getNPar(dm)) 
    stop("Wrong number of parameters!")

  .log2("calling ms to generate tree...")
  ms.options <- generateMsOptions(dm, parameters)
  ms.file <- callMs(ms.options)

  .log2("running seq-gen")
  seqgen.options <- generateSeqgenOptions(dm, parameters)
  .log3("options generated")
  #sim.time <- system.time(print(1))
  seqgen.file  <- callSeqgen(seqgen.options, ms.file)
  #.log3("finished after", sum(sim.time[-3]), "seconds")
  .log3("simulation output in file", seqgen.file)

  .log2("calculating jsfs")
  jsfs <- seqgenOut2Jsfs(dm, seqgen.file)
  #jsfs <- matrix(0,  dm@sampleSizes[1] + 1, dm@sampleSizes[2] + 1)

  .log3("done. Removing tmp files...")
  unlink(seqgen.file)
  unlink(ms.file)
  return(jsfs)
}

createSimProgram("seq-gen", "",
                 possible.features,
                 possible.sum.stats,
                 singleSimFunc=seqgenSingleSimFunc)

#test code:
#parameters <- c(1,5)
#dm <- dm.createThetaTauModel(24:25, 100, 1000)
#dm <- jaatha:::dm.addOutgroup(dm, 1, "2*tau")
#jsfs <- jaatha:::seqgenSingleSimFunc(dm, c(1,5))
