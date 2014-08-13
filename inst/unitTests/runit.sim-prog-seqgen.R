test.generateSeqgenOptions <- function() {
  jaatha:::setJaathaVariable('seqgen.exe', 'seq-gen') 
  dm.hky@options$seqgen.cmd <- NULL
  opts <- jaatha:::generateSeqgenOptions(dm.hky, c(1,10))
  opts <- strsplit(opts, " ")[[1]]
  checkTrue(opts[1] == "seq-gen")  
  checkTrue("-l" %in% opts)  
  checkTrue("-p" %in% opts)  
  checkTrue("-z" %in% opts)  
  checkTrue("-q" %in% opts)  
  checkTrue("-mHKY" %in% opts)  
  checkTrue("-t" %in% opts) 
  checkTrue("-f" %in% opts) 
  checkTrue("-s" %in% opts) 
}

test.callSeqgen <- function(opts, ms.file){
  if (!checkForSeqgen(FALSE)) {
    warning('Can not test seqgen: jar not found')
    return()
  }
  
  opts <- c("seq-gen", " -mHKY", " -l", dm.getLociLength(dm.tt), " -p", dm.getLociLength(dm.tt) + 1, " -q")
  
  dm.tt <- dm.addSummaryStatistic(dm.tt, 'trees')
  ms.options <- jaatha:::generateMsOptions(dm.tt, c(1,10))
  ms.file <- jaatha:::callMs(ms.options, dm.tt)
  
  seqgen.file <- callSeqgen(opts, ms.file)
  
  checkTrue( file.exists(seqgen.file) )
  checkTrue( file.info(seqgen.file)$size != 0 )
  unlink(ms.file)
  unlink(seqgen.file)
}

test.seqgenSingleSimFunc <-function() {
  if (!checkForSeqgen(FALSE)) {
    warning('Can not test seqgen: jar not found')
    return()
  }
  
  seqgenSingleSimFunc = getSimProgram('seq-gen')$sim_func
  checkException(seqgenSingleSimFunc(dm.tt, c(1,10)));
  set.seed(100)
  sum.stats <- seqgenSingleSimFunc(dm.hky, c(1,10));
  checkTrue(is.list(sum.stats))
  checkTrue(is.array(sum.stats$jsfs))
  checkTrue(sum(sum.stats$jsfs) > 0)
  
  set.seed(100)
  sum.stats2 <- seqgenSingleSimFunc(dm.hky, c(1,10));
  checkEqualsNumeric(sum.stats$jsfs, sum.stats2$jsfs)
}

test.HkyModel <- function() {
  if (!checkForSeqgen(FALSE)) {
    warning('Can not test seqgen: jar not found')
    return()
  }
  set.seed(12)
  jsfs <- dm.simSumStats(dm.hky, c(1,10));
  checkTrue(sum(jsfs$jsfs) > 0)
}

test.F81Model <- function() {
  if (!checkForSeqgen(FALSE)) {
    warning('Can not test seqgen: jar not found')
    return()
  }
  set.seed(12)
  jsfs <- dm.simSumStats(dm.f81, c(1,10));
  checkTrue(sum(jsfs$jsfs) > 0)
}

test.GtrModel <- function() {
  if (!checkForSeqgen(FALSE)) {
    warning('Can not test seqgen: jar not found')
    return()
  }
  set.seed(12)
  jsfs <- dm.simSumStats(dm.gtr, c(1,10));
  checkTrue(sum(jsfs$jsfs) > 0)
}

test.RateHeterogenity <- function() {
  if (!checkForSeqgen(FALSE)) {
    warning('Can not test seqgen: jar not found')
    return()
  }
  set.seed(12)
  dm.rh <- dm.addMutationRateHeterogenity(dm.hky, 0.1, 5,
                                          categories.number=5)
  
  jsfs <- dm.simSumStats(dm.rh, c(1, 10, 1));
  checkTrue(sum(jsfs$jsfs) > 0)
}

test.finalizeSeqgen <- function() {
  if (!checkForSeqgen(FALSE)) {
    warning('Can not test seqgen: jar not found')
    return()
  }
  
  finalizeSeqgen = getSimProgram('seq-gen')$finalization_func
  dm.hky <- finalizeSeqgen(dm.hky)
  dm.f81 <- finalizeSeqgen(dm.f81)
  dm.gtr <- finalizeSeqgen(dm.gtr)
  
  checkTrue(!is.null(dm.hky@options[['seqgen.cmd']]))
  checkTrue(!is.null(dm.hky@options[['tree.model']]))
  checkTrue(!is.null(dm.f81@options[['seqgen.cmd']]))
  checkTrue(!is.null(dm.f81@options[['tree.model']]))
  checkTrue(!is.null(dm.gtr@options[['seqgen.cmd']]))
  checkTrue(!is.null(dm.gtr@options[['tree.model']]))
}

test.generateTreeModel <- function() {
  for (dm in c(dm.hky, dm.f81, dm.gtr)) { 
    dm.ms <- dm.finalize(generateTreeModel(dm))
    sum.stats <- dm.simSumStats(dm.ms, c(1,5))
    checkTrue( !is.null(sum.stats$file) )
    checkTrue( file.exists(sum.stats$file) )
    unlink(sum.stats$file)
  }
}

test.simulateFpcWithSeqgen <- function() {
  if (!checkForSeqgen(FALSE)) {
    warning('Can not test seqgen: binary not found')
    return()
  }
  
  sum.stats <- dm.simSumStats(dm.sgfpc, c(1,5))
  checkTrue( !is.null(sum.stats$fpc) )
  checkEquals( 5, sum(sum.stats$fpc) )
}

test.seqgenMutationParameterNotLast <- function() {
  dm.test <- dm.hky
  dm.test@parameters <- dm.test@parameters[3:1, ]
  cmd <- paste(generateSeqgenOptionsCmd(dm.test), collapse='')
  checkEquals( 1, grep("theta", cmd) )
  checkEquals( integer(0), grep("tau", cmd) )
  checkEquals( integer(0), grep("rho", cmd) )
}

test.seqgenWithMsms <- function() {
  if ((!checkForSeqgen(FALSE)) | !checkForMsms(FALSE)) {
    warning('Can not test seqgen with msms: not found')
    return()
  }
  dm.selsq <- dm.addPositiveSelection(dm.f81, 100, 500, 
                                      population=1, at.time="0.1")
  dm.selsq <- dm.finalize(dm.selsq)
  checkTrue(!is.null(dm.selsq@options[['seqgen.cmd']]))
  checkTrue(!is.null(dm.selsq@options[['tree.model']]))
  
  set.seed(4444)
  sum.stats <- dm.simSumStats(dm.selsq, c(1, 5, 250))
  checkTrue( !is.null(sum.stats$jsfs) )
  
  set.seed(4444)
  sum.stats2 <- dm.simSumStats(dm.selsq, c(1, 5, 250))
  checkTrue( !is.null(sum.stats2$jsfs) )
  checkEquals(sum.stats, sum.stats2)
}