# --------------------------------------------------------------
# sim_prog_seqgen.R
# Adaptor to calling ms from a demographic model.
# 
# Authors:  Lisha Mathew & Paul R. Staab
# Date:     2012-10-25
# Licence:  GPLv3 or later
# --------------------------------------------------------------

# list ms's features + FS related features
simulband.features    <- c('band.number', 'band.freq')
possible.sum.stats <- c("jsfs")

possible.features  <- c(getSimProgram('ms')@possible.features, simulband.features)

checkForSimulBands <- function() {
  if ( isJaathaVariable('simulbands.exe') ) return()

  # Works on Linux only maybe
  run.path <- strsplit(Sys.getenv("PATH"), ":")[[1]]
  executables <- c(paste(run.path, "/SimulBands", sep=""), 
                   paste(run.path, "/simulbands", sep=""))
  for (exe in executables) {
    if (file.exists(exe)) {
      .print("Using", exe, "as SimulBands executable\n")
      setJaathaVariable('simulbands.exe', exe)     
      return()
    }
  }

  stop("No SimulBands executable found. Please provide one using
       Jaatha.setSimulBandsExecutable()")
}

#' Set the path to the executable for seqgen
#'
#' @param seqgen.exe Path to seqgen's executable.
#' @export
Jaatha.setSimulBandsExecutable <- function(seqgen.exe) {
  if (file.exists(seqgen.exe)) {
    setJaathaVariable('simulbands.exe', seqgen.exe)     
    .print("Using", seqgen.exe, "as SimulBands executable\n")
  } else {
    stop("File ", seqgen.exe, " does not exist")
  }
}

simulbandsSingleSimFunc <- function(dm, parameters) {
  .log3("called msSingleSimFunc()")
  .log3("parameter:",parameters)
  checkType(dm, "dm")
  checkType(parameters, "num")

  if (length(parameters) != dm.getNPar(dm)) 
    stop("Wrong number of parameters!")

  opts <- c() 

  # Get the executable
  checkForSimulBands()
  opts[1] <- getJaathaVariable('simulbands.exe')
  opts[2] <- "tree"

  # Generate a tree with ms
  .log2("calling ms to generate tree...")
  dm.ms <- dm
  dm.ms@sampleSizes <- dm.ms@sampleSizes * 2
  ms.options <- generateMsOptions(dm.ms, parameters)
  ms.file <- callMs(ms.options, dm.ms)
  tr <- readLines(ms.file)

  # Get remaining options
  band.number.row <- dm@parameters$name == "band.number"
  len <- dm@parameters[band.number.row, 'lower.range']
  opts[3] <- len 
  band.freq.row <- dm@parameters$name == "band.freq"
  opts[4] <- dm@parameters[band.freq.row, 'lower.range']
  
  indiv1 <- dm@sampleSizes[1]
  indiv2 <- dm@sampleSizes[2]
  M <- matrix(rep(0,len*(indiv1+indiv2)*2),nrow=2*(indiv1+indiv2))

  # Get the mutation rate
  rate <- parameters[dm.getParameters(dm) == "theta"]
  opts[5] <- "rate"
  opts[6] <- "seed"
  opts[7] <- " > "
  opts[8] <- "outfile"

  for (i in 3:length(tr)) {
    if (tr[i] == "") stop()
    if (tr[i] == "//") stop()

    trlen <- as.numeric(strsplit(tr[i],"[][]")[[1]][2])
    tree <- paste("\"",strsplit(tr[i],"]")[[1]][2],"\"",sep="")
    opts[2] <- tree
    opts[5] <- rate*trlen/100
    opts[6] <- generateSeeds(1)
  
    simulbands.file <- getTempFile("simulbands") 
    opts[8] <- simulbands.file 

    # Generate the cmd
    simulbands.cmd <- paste(opts, collapse=" ")

    .log2("running SimulBands")
    sim.time <- system.time(system(simulbands.cmd))
    .log3("finished after", sum(sim.time[-3]), "seconds")

    # Parse the output
    dat <- readLines(simulbands.file)
    unlink(simulbands.file)

    for(j in 1:length(dat)) {
      s <- dat[j]
      ss <- strsplit(s," +")[[1]]
      sl <- as.numeric(substring(ss[1], 2))
      s <- as.numeric(strsplit(ss[2],"")[[1]])
      M[sl,] <- M[sl,] | s
    }
  }

  r1 <- sample(1:(2*indiv1))
  r2 <- sample(1:(2*indiv2))
  M <- apply( M[c(r1[1:indiv1],2*indiv1+r2[1:indiv2]),]|
         M[c(r1[(indiv1+1):(2*indiv1)],2*indiv1+r2[(indiv2+1):(2*indiv2)]),] , 
         2 , as.numeric)

  jsfs <- convertBandsToJSFS(M, indiv1, indiv2)
  unlink(ms.file)
  .log3("Simubands simulation succesfully finished")
  
  return(jsfs)
}


finalizeSimulbands <- function(dm) {
  checkForSimulBands()
  dm <- finalizeMs(dm)
  return(dm)
}

convertBandsToJSFS <- function(B,indiv1,indiv2) {
  if(dim(B)[1]!=indiv1+indiv2) {
    print("Error: Total sample size does not match number of lines in B.\n")
  } else {
    n <- dim(B)[2]
    x <- apply(B[1:indiv1,],2,sum)
    y <- apply(B[(indiv1+1):(indiv1+indiv2),],2,sum)
    M <- matrix(rep(0,(indiv1+1)*(indiv2+1)),nrow=indiv1+1)
    dimnames(M) <- list(0:indiv1,0:indiv2)
    for(i in 1:n) {
      M[x[i]+1,y[i]+1] <- M[x[i]+1,y[i]+1]+1
    }
  }
  return(M)
}

createSimProgram("simulbands", "",
                 possible.features,
                 possible.sum.stats,
                 singleSimFunc=simulbandsSingleSimFunc,
                 finalizationFunc=finalizeSimulbands)
