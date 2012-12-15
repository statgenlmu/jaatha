dm <- dm.createThetaTauModel(10:11, 10, 100)
dm <- dm.addOutgroup(dm, "2*tau")
dm.hky <- dm.setMutationModel(dm, "HKY", c(0.2, 0.2, 0.3, 0.3), 2)
dm.f81 <- dm.setMutationModel(dm, "F84", c(0.3, 0.2, 0.3, 0.2), 2)
dm.gtr <- dm.setMutationModel(dm, "GTR", gtr.rates=c(0.2, 0.2, 0.1, 0.1, 0.1, 0.2))
pars <- c(1, 10)

test.generateSeqgenOptions <- function() {
  jaatha:::setJaathaVariable('seqgen.exe', 'seq-gen') 
  opts <- jaatha:::generateSeqgenOptions(dm.hky, pars)
  opts <- strsplit(opts, " ")[[1]]
  checkTrue(opts[1] == "seq-gen")  
  checkTrue("-l" %in% opts)  
  checkTrue("-p" %in% opts)  
  checkTrue("-z" %in% opts)  
  checkTrue("-q" %in% opts)  
  checkTrue("-mHKY" %in% opts)  
  checkTrue("-t" %in% opts) 
  checkTrue("-f" %in% opts) 
}

test.HkyModel <- function() {
  set.seed(12)
  jsfs <- jaatha:::seqgenSingleSimFunc(dm.hky, pars);
  checkTrue(sum(jsfs) > 0)
}

test.F81Model <- function() {
  set.seed(12)
  jsfs <- jaatha:::seqgenSingleSimFunc(dm.f81, pars);
  checkTrue(sum(jsfs) > 0)
}

test.GtrModel <- function() {
  set.seed(12)
  jsfs <- jaatha:::seqgenSingleSimFunc(dm.gtr, pars);
  checkTrue(sum(jsfs) > 0)
}

test.RateHeterogenity <- function() {
  set.seed(12)
  dm.rh <- dm.addMutationRateHeterogenity(dm.hky, 0.1, 5, categories.number=5)

  opts <- jaatha:::generateSeqgenOptions(dm.rh, pars)
  opts <- strsplit(opts, " ")[[1]]
  checkTrue("-a" %in% opts)
  checkTrue("-g" %in% opts)

  jsfs <- jaatha:::seqgenSingleSimFunc(dm.rh, c(1,2,10))
  checkTrue(sum(jsfs) > 0)
}

test.callSeqgen <- function(opts, ms.file){
  opts <- c("seq-gen", " -mHKY", " -l", dm@seqLength, " -p", dm@seqLength + 1, " -q")
  
  ms.options <- jaatha:::generateMsOptions(dm, pars)
  ms.file <- jaatha:::callMs(ms.options)

  seqgen.file <- callSeqgen(opts, ms.file)
  checkTrue( file.exists(seqgen.file) )
  checkTrue( file.info(seqgen.file)$size != 0 )
}

test.seqgenSingleSimFunc <-function() {
  checkException(seqgenSingleSimFunc(dm, pars));
  set.seed(100)
  jsfs <- jaatha:::seqgenSingleSimFunc(dm.hky, pars);
  checkTrue(sum(jsfs) > 0)
  
  set.seed(100)
  jsfs2 <- jaatha:::seqgenSingleSimFunc(dm.hky, pars);
  checkEqualsNumeric(jsfs, jsfs2)
}
