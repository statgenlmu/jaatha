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

test.finalizeSeqgen <- function() {
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
    dm.ms <- finalizeMs(generateMsModel(dm))
    sum.stats <- msSingleSimFunc(dm.ms, c(1,5))
    checkTrue( !is.null(sum.stats$file) )
    checkTrue( file.exists(sum.stats$file) )
    unlink(sum.stats$file)
  }
}
