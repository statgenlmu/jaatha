# --------------------------------------------------------------
# sim_prog_msms.R
# Adaptor to calling ms from a demographic model.
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Date:     2012-10-05
# Licence:  GPLv3 or later
# --------------------------------------------------------------

msms <- function(jar.path, ms.args, msms.args, pop.sizes, sum.stats=c('raw'),
                 out.file=NULL) {

  # Check arguments
  if (!all(sum.stats %in% c('raw', 'jsfs', 'seg.sites', 'four.point.condition')))
    stop("Unknown summary statistic")

  # Complement/Modify arguments
  seed <- generateSeeds(1)
  if ('jsfs' %in% sum.stats) msms.args <- paste(msms.args, "-oAFS jAFS")


  # Create the command
  command = paste("java -jar", jar.path, as.character(msms.args), 
                  "-ms", as.character(ms.args), "-seed", seed)
  if (!is.null(out.file)) {
    stopifnot(!file.exists(out.file))
    command <- paste(command, ">", out.file)
  }

  # Execute the command
  output <- system(command, intern=TRUE)

  # Parse output and return summary statistics
  sample.sum.stats <- list()

  if ('jsfs' %in% sum.stats) {
    sample.sum.stats[['jsfs']] <- readJsfsFromOutput(output)
  }

  if ('seg.sites' %in% sum.stats) {
    sample.sum.stats[['seg.sites']] <- readSegSitesFromOutput(output, pop.sizes)
  }

  if ('four.point.condition' %in% sum.stats) {
    if (!is.null(sample.sum.stats[['seg.sites']])) 
      seg.sites <- sample.sum.stats[['seg.sites']] 
    else seg.sites <- readSegSitesFromOutput(output, pop.sizes)
    sample.sum.stats[['four.point.condition']] <- 
      generateFourPointStat(seg.sites, pop.sizes)
  }

  if ('raw' %in% sum.stats) {
    sample.sum.stats[['raw']] <- output
  }

  return(sample.sum.stats)
}

readJsfsFromOutput <- function(output) {
  jsfs.begin <- which(output == "Summary jAFS")+2
  if (output[jsfs.begin-1] != "jAFS 0 vrs 1") 
    stop ("Error parsing msms output")

  jsfs.end <- length(output)-1 
  t(sapply(jsfs.begin:jsfs.end, 
           function(x) as.integer(unlist(strsplit(output[x], " ")))))
}

readSegSitesFromOutput <- function(output, pop.sizes) {
  seg.sites.begin <- which(grepl('^segsites: [0-9]+$', output))

  lapply(seg.sites.begin, function(begin) { 
    positions <- as.numeric(strsplit(output[begin+1], ' ')[[1]][-1])
    stopifnot( length(positions) != 0 )
    seg.sites.char <- output[1:sum(pop.sizes)+begin+1]
    seg.sites <- matrix(as.integer(unlist(strsplit(seg.sites.char, split= ''))), 
                        length(seg.sites.char), byrow=TRUE) 
    colnames(seg.sites) <- positions
    seg.sites
  })
}

generateFourPointStat <- function(seg.sites, pop.sizes) {
  four.point.cond <- rep(0,6)
  pop.one <- c(rep(TRUE, pop.sizes[1]), rep(FALSE, pop.sizes[2]))
  for(locus.snps in seg.sites) {

    four.point.cond <- four.point.cond + 
    c(countViolationsNextSnp(locus.snps[pop.one, ]),
      countViolationsNextSnp(locus.snps[!pop.one, ]),
      countViolationsNextSnp(locus.snps),
      countViolationsAllSnps(locus.snps[pop.one, ]),
      countViolationsAllSnps(locus.snps[!pop.one, ]),
      countViolationsAllSnps(locus.snps))
  }

  four.point.cond
}

countViolationsAllSnps <- function(snp.matrix) {
  count <- 0
  for (i in 1:ncol(snp.matrix)) {
    for (j in 1:ncol(snp.matrix)) {
      if (i >= j) next()
      count = count + violatesFourPointCondidtion(snp.matrix[,i], snp.matrix[,j])
    }
  } 
  count
}

countViolationsNextSnp <- function(snp.matrix) {
  sum(sapply(2:ncol(snp.matrix), 
             function(i) violatesFourPointCondidtion(snp.matrix[,i-1],
                                                     snp.matrix[,i])))
}

violatesFourPointCondidtion <- function(site.one, site.two) {
  status <- site.one * 2 + site.two
  if (all(0:3 %in% status)) return(TRUE)
  return(FALSE) 
}


checkForMsms <- function() {
  if ( isJaathaVariable('msms.jar') ) return()

  # Works on Linux only maybe
  run.path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
  executables <- paste(run.path, "/msms.jar", sep="")
  for (exe in executables) {
    if (file.exists(exe)) {
      .print("Using", exe, "as msms implementation\n")
      setJaathaVariable('msms.jar', exe)     
      return()
    }
  }

  stop("No msms executable found.")
}

possible.features  <- c(getSimProgram('ms')@possible.features,"pos.selection")
possible.sum.stats <- c("jsfs")

# This function generates an string that contains an R command for generating
# an ms call to the current model.
generateMsmsOptionsCommand <- function(dm) {
  cmd <- c('c(')

  for (i in 1:dim(dm@features)[1] ) {
    type <- as.character(dm@features[i,"type"])
    feat <- unlist(dm@features[i, ])

    if (type == "pos.selection") {
      cmd <- c(cmd, '"-SI"', ',', feat['time.point'], ',', length(dm@sampleSizes), ',')
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

  ms.options <- paste(sum(dm@sampleSizes), dm@nLoci, 
                      paste(generateMsOptions(dm, parameters), collapse=" "))
  msms.options <- paste(generateMsmsOptions(dm, parameters), collapse= " ") 

  sim.time <- system.time(sum.stats <- msms(getJaathaVariable('msms.jar'), 
                                       ms.options, msms.options, dm@sampleSizes, 
                                       c('jsfs', 'four.point.condition')) )

  .log3("finished after", sum(sim.time[-3]), "seconds")
  sum.stats[['pars']] <- parameters
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
