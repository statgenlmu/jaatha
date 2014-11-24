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
                           'gamma.categories', 'gamma.rate',
                           'trio.1', 'trio.2', 'trio.3', 
                           'trio.4', 'trio.5'))

sg.sum.stats <- c('jsfs', 'file', 'seg.sites', 'fpc', 'pmc')
sg.mutation.models <- c('HKY', 'F84', 'GTR')

checkForSeqgen <- function(throw.error = TRUE, silent = FALSE) {
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
      if (!silent) message(paste("Using", exe, "as seqgen executable\n"))
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
  dm <- dm.addSummaryStatistic(dm, "trees")
  dm <- dm.addSummaryStatistic(dm, "file")
  dm <- dm.setLociNumber(dm, 1)
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
callSeqgen <- function(opts, ms_files) {
  stopifnot(!missing(opts))
  stopifnot(length(opts) == length(ms_files))

  sapply(seq(along = opts), function(i) {
    if(!file.exists(ms_files[i])) stop("ms file not found")
    if(file.info(ms_files[i])$size == 0 ) stop("ms output is empty")
    
    seqgen_file <- getTempFile("seqgen")
    cmd <- paste(opts[i], "<", ms_files[i], ">", seqgen_file)
    
    # Do the acctual simulation
    if (system(cmd, intern = FALSE) != 0) stop("seq-gen simulation failed")
    
    if( !file.exists(seqgen_file) ) stop("seq-gen simulation failed!")
    if( file.info(seqgen_file)$size == 0 ) stop("seq-gen output is empty!")
    
    seqgen_file
  })
}

generateSeqgenOptions <- function(dm, parameters, locus, trio_opt = NA) {
  # Generate the command template to execute or use the buffered one
  if ( !is.null( dm@options[['seqgen.cmd']] ) ) {
    cmd <- dm@options[['seqgen.cmd']]
  } else {
    cmd <- generateSeqgenOptionsCmd(dm)
  }
  
  # Get the length of the loci we simulate
  if (is.na(trio_opt) || length(trio_opt) == 0) {
    locus_lengths <- dm.getLociLength(dm)
  } else if (length(trio_opt) == 5) {
    locus_lengths <- trio_opt[c(1,3,5)]
  } else {
    print(trio_opt)
    stop('failed to parse trio options')
  }
  
  # Fill the parameters in the template
  sapply(locus_lengths, function(locus_length) {
    par_envir <- createParameterEnv(dm, parameters, locus = locus, 
                                    locus_length = locus_length,
                                    seed = generateSeeds(1))
    paste(eval(parse(text=cmd), envir=par_envir), collapse=" ")
  })
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

  opts <- c(opts, '"-l"', ',', 'locus_length', ',')
  opts <- c(opts, '"-s"', ',', paste(getThetaName(dm), ' / locus_length'), ',')
  opts <- c(opts, '"-p"', ',', 'locus_length + 1', ',')
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
  
  seqgen.files <- lapply(1:dm.getLociNumber(dm), function(locus) {
    # Generate options for seqgen
    seqgen.options <- generateSeqgenOptions(dm, parameters, locus,
                                            dm.getLociTrioOptions(dm))
    
    # Simulate the trees
    sum_stats_ms <- dm.simSumStats(tree.model, parameters)
    #print(sum_stats_ms[['file']])
    tree_files <- parseTrees(sum_stats_ms[['file']][[1]],
                             dm.getLociTrioOptions(dm),
                             getTempFile)
    #print(tree_files)
    
    # Call seq-gen to distribute mutations
    seqgen.file <- callSeqgen(seqgen.options, tree_files)
    
    # Delete tree files
    unlink(c(tree_files, sum_stats_ms[['file']]))
    seqgen.file
  })
  
  # Generate the summary statistics
  generateSumStats(seqgen.files, 'seqgen', parameters, dm)
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