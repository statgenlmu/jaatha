dm <- dm.createThetaTauModel(10:11, 10, 100)
dm <- dm.setMutationModel(dm, "HKY", c(0.1, 0.2, 0.3, 0.4), 2)
dm <- jaatha:::dm.addOutgroup(dm, 1, "2*tau")
pars <- c(1, 10)

test.generateSeqgenOptions <- function() {
  jaatha:::setJaathaVariable('seqgen.exe', 'seq-gen') 
  opts <- jaatha:::generateSeqgenOptions(dm, pars)
  opts <- strsplit(opts, " ")[[1]]
  checkTrue(opts[1] == "seq-gen")  
  checkTrue("-l" %in% opts)  
  checkTrue("-p" %in% opts)  
  checkTrue("-z" %in% opts)  
  checkTrue("-q" %in% opts)  
  checkTrue("-m" %in% opts)  
  checkTrue("HKY" %in% opts) 
  checkTrue("-t" %in% opts) 
  checkTrue("-f" %in% opts) 
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
#  jaatha:::setLogging(3)
#  set.seed(100)
#  jsfs <- jaatha:::seqgenSingleSimFunc(dm, pars);
#  checkTrue(sum(jsfs) > 0)
  
#  set.seed(100)
#  jsfs2 <- jaatha:::seqgenSingleSimFunc(dm, pars);
#  checkEqualsNumeric(jsfs, jsfs2)
}
