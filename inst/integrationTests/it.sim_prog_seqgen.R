test.callSeqgen <- function(opts, ms.file){
  opts <- c("seq-gen", " -mHKY", " -l", dm.tt@seqLength, " -p", dm.tt@seqLength + 1, " -q")
  
  ms.options <- jaatha:::generateMsOptions(dm.tt, c(1,10))
  ms.file <- jaatha:::callMs(ms.options, dm.tt)

  seqgen.file <- callSeqgen(opts, ms.file)
  checkTrue( file.exists(seqgen.file) )
  checkTrue( file.info(seqgen.file)$size != 0 )
}

test.seqgenSingleSimFunc <-function() {
  checkException(seqgenSingleSimFunc(dm.tt, c(1,10)));
  set.seed(100)
  sum.stats <- jaatha:::seqgenSingleSimFunc(dm.hky, c(1,10));
  checkTrue(is.list(sum.stats))
  checkTrue(is.array(sum.stats$jsfs))
  checkTrue(sum(sum.stats$jsfs) > 0)
  
  set.seed(100)
  sum.stats2 <- jaatha:::seqgenSingleSimFunc(dm.hky, c(1,10));
  checkEqualsNumeric(sum.stats$jsfs, sum.stats2$jsfs)
}

test.HkyModel <- function() {
  set.seed(12)
  jsfs <- jaatha:::seqgenSingleSimFunc(dm.hky, c(1,10));
  checkTrue(sum(jsfs$jsfs) > 0)
}

test.F81Model <- function() {
  set.seed(12)
  jsfs <- jaatha:::seqgenSingleSimFunc(dm.f81, c(1,10));
  checkTrue(sum(jsfs$jsfs) > 0)
}

test.GtrModel <- function() {
  set.seed(12)
  jsfs <- jaatha:::seqgenSingleSimFunc(dm.gtr, c(1,10));
  checkTrue(sum(jsfs$jsfs) > 0)
}

test.RateHeterogenity <- function() {
  set.seed(12)
  dm.rh <- dm.addMutationRateHeterogenity(dm.hky, 0.1, 5, categories.number=5)

  opts <- jaatha:::generateSeqgenOptions(dm.rh, c(1, 1, 10))
  opts <- strsplit(opts, " ")[[1]]
  checkTrue("-a" %in% opts)
  checkTrue("-g" %in% opts)

  jsfs <- jaatha:::seqgenSingleSimFunc(dm.rh, c(1, 1, 10));
  checkTrue(sum(jsfs$jsfs) > 0)
}
