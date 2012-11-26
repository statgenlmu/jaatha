dm <- dm.createThetaTauModel(24:25, 100, 1000)
dm <- jaatha:::dm.addOutgroup(dm, 1, "2*tau")
pars <- c(1, 10)

test.generateSeqgenOptions <- function() {
  opts <- generateSeqgenOptions(dm, pars)
  checkTrue(opts[1] == "seq-gen")  
}

test.callSeqgen <- function(opts, ms.file){
  opts <- c("seq-gen", " -mHKY", " -l", dm@seqLength, " -p", dm@seqLength + 1, " -s 0.0002", " -q")
  
  ms.options <- generateMsOptions(dm, pars)
  ms.file <- callMs(ms.options)
 
  seqgen.file <- callSeqgen(opts, ms.file)
  checkTrue( file.exists(seqgen.file) )
  checkTrue( file.info(seqgen.file)$size != 0 )
}

test.seqgenSingleSimFunc <-function() {
  jsfs <- seqgenSingleSimFunc(dm, pars);
  checkTrue(sum(jsfs) > 0)
}
