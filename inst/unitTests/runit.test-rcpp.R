test.parseMsOutput <- function() {
  dm.tt@sum.stats <- 'file'
  sum.stats <- dm.simSumStats(dm.tt, c(1,5))
  checkException(parseMsOutput('bulb.txt', dm.getSampleSizes(dm.tt)), s=TRUE)
  #checkException(parseMsOutput(sum.stats$file, dm.getSampleSizes(dm.tt)+1), s=TRUE)
  print(parseMsOutput(sum.stats$file, dm.getSampleSizes(dm.tt)))
}
