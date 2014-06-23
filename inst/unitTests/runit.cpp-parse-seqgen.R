test.parseOutputSeqgen <- function() {
  sg.file <- tempfile("seqgen_")
  write(sg.example, sg.file)
  
  checkException( parseOutput(tempfile("seqgen2_"), dm.getSampleSize(dm.hky), dm.getLociNumber(dm.hky), 1) )
  checkException( parseOutput(sg.file, dm.getSampleSize(dm.hky), dm.getLociNumber(dm.hky)-1, 1) )
  checkException( parseOutput(sg.file, dm.getSampleSize(dm.hky), dm.getLociNumber(dm.hky)+1, 1) )
  
  sum.stats = parseOutput(sg.file, dm.getSampleSize(dm.hky), dm.getLociNumber(dm.hky), 1)
  checkTrue( is.matrix(sum.stats$jsfs) )
  checkTrue( sum(sum.stats$jsfs) > 0 )
  checkTrue( sum.stats$jsfs[1,1] == 0 )
  checkTrue( sum.stats$jsfs[12,13] == 0 )
  
  sum.stats = parseOutput(sg.file, dm.getSampleSize(dm.hky), dm.getLociNumber(dm.hky), 1, generate_seg_sites = TRUE)
  checkTrue( is.matrix(sum.stats$jsfs) )
  checkTrue( is.list(sum.stats$seg.sites) )
  checkEquals( dm.getLociNumber(dm.hky), length(sum.stats$seg.sites) )
  checkEquals( (c(1, 2, 4, 5)-1)/14, as.numeric(colnames(sum.stats$seg.sites[[1]])[1:4]) )
  checkEquals( (c(2, 6, 7)-1)/14, as.numeric(colnames(sum.stats$seg.sites[[2]])[1:3]) )
  
  sum.stats = parseOutput(sg.file, dm.getSampleSize(dm.hky), dm.getLociNumber(dm.hky), 1, generate_fpc = TRUE,
                          fpc_breaks_near = 1:3/4, fpc_breaks_far = 1:3/4)
  checkTrue( is.matrix(sum.stats$fpc) )
  checkEquals( 5, sum(sum.stats$fpc))
}
