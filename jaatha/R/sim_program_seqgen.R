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
seqgen.features    <- c()
fossible.sum.stats <- c("jsfs")

# possible.features  <- c(getSimProgram('ms')@possible.features, seqgen.features)

# Function to perform simulation using seqgen
# 
# @param opts The options to pass to ms. Must either be a character or character
# vector.
callSeqgen <- function(opts, ms.file){
  if (missing(opts)) stop("No options given!")
  .log3(1)
  opts[length(opts) + 1:2] <- c(" ", ms.file)
  .log3(2)
  opts <- paste(opts, collapse=" ")
  .log3(3)

  if( !file.exists(ms.file) ) stop("ms file not found")
  .log3(4)
  if( file.info(ms.file)$size == 0 ) stop("ms output is empty")
  .log3(5)

  seqgen.file <- getTempFile("seqgen")
  .log3(6)
  cmd <- paste(opts, ">", seqgen.file, sep=" ", collapse=" ")
  .log3(7)
  .log3("executing: '", cmd, "'")
  system(cmd)
  
  if( !file.exists(seqgen.file) ) stop("seq-gen simulation failed!")
  if( file.info(seqgen.file)$size == 0 ) stop("seq-gen output is empty!")
  return(seqgen.file)
}

generateSeqgenOptions <- function(dm, parameters) {
  #return(c("-mHKY", "-l", dm@seqLength, "-p", dm@seqLength + 1))
  return(c("seq-gen", "-mHKY", "-l", dm@seqLength, "-p", dm@seqLength + 1, 
           "-s 0.0002", "-z ", generateSeeds(1), "-q"))
}

printSeqgenCommand <- function(dm) {
  cmd <- generateSeqgenOptions(dm)
  
  cmd <- cmd[cmd != ","]
  #cmd <- cmd[-c(1, length(cmd))]

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

# createSimProgram("seq-gen", "",
#                  possible.features,
#                  possible.sum.stats,
#                  singleSimFunc=seqgenSingleSimFunc)

#test code:
#parameters <- c(1,5)
#dm <- dm.createThetaTauModel(24:25, 100, 1000)
#dm <- jaatha:::dm.addOutgroup(dm, 1, "2*tau")
#jsfs <- jaatha:::seqgenSingleSimFunc(dm, c(1,5))
