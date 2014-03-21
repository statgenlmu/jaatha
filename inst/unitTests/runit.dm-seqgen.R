  # Most unit tests for seq-gen are moved to integration tests, as seq-gen is not
  # available on CRAN. 

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
