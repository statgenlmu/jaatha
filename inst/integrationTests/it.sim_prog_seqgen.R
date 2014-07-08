test.callSeqgen <- function(opts, ms.file){
  checkForSeqgen()
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
  set.seed(12)
  jsfs <- dm.simSumStats(dm.hky, c(1,10));
  checkTrue(sum(jsfs$jsfs) > 0)
}

test.F81Model <- function() {
  set.seed(12)
  jsfs <- dm.simSumStats(dm.f81, c(1,10));
  checkTrue(sum(jsfs$jsfs) > 0)
}

test.GtrModel <- function() {
  set.seed(12)
  jsfs <- dm.simSumStats(dm.gtr, c(1,10));
  checkTrue(sum(jsfs$jsfs) > 0)
}

test.RateHeterogenity <- function() {
  set.seed(12)
  dm.rh <- dm.addMutationRateHeterogenity(dm.hky, 0.1, 5,
                                          categories.number=5)

  jsfs <- dm.simSumStats(dm.rh, c(1, 10, 1));
  checkTrue(sum(jsfs$jsfs) > 0)
}

test.finalizeSeqgen <- function() {
  finalizeSeqgen = getSimProgram('seq-gen')$finalization_func
  dm.hky <- finalizeSeqgen(dm.hky)
  dm.f81 <- finalizeSeqgen(dm.f81)
  dm.gtr <- finalizeSeqgen(dm.gtr)

  checkTrue(!is.null(dm.hky@options[['seqgen.cmd']]))
  checkTrue(!is.null(dm.hky@options[['ms.model']]))
  checkTrue(!is.null(dm.f81@options[['seqgen.cmd']]))
  checkTrue(!is.null(dm.f81@options[['ms.model']]))
  checkTrue(!is.null(dm.gtr@options[['seqgen.cmd']]))
  checkTrue(!is.null(dm.gtr@options[['ms.model']]))
}

test.generateMsModel <- function() {
  for (dm in c(dm.hky, dm.f81, dm.gtr)) { 
    dm.ms <- dm.finalize(generateMsModel(dm))
    sum.stats <- dm.simSumStats(dm.ms, c(1,5))
    checkTrue( !is.null(sum.stats$file) )
    checkTrue( file.exists(sum.stats$file) )
    unlink(sum.stats$file)
  }
}

test.simualteFpcWithSeqgen <- function() {
  sum.stats <- dm.simSumStats(dm.sgfpc, c(1,5))
  checkTrue( !is.null(sum.stats$fpc) )
  checkEquals( 5, sum(sum.stats$fpc) )
}