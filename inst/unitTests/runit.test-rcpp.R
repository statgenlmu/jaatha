test.parseOutput <- function() {
  dm.tt@sum.stats <- 'file'
  ms.file <- dm.simSumStats(dm.tt, c(1,5))$file
  checkException(parseOutput('bulb.txt', dm.getSampleSize(dm.tt), dm.getLociNumber(dm.tt)), s=TRUE)
  #checkException(parseOutput(sum.stats$file, dm.getSampleSize(dm.tt)+1), s=TRUE)
  sum.stats <- parseOutput(ms.file, dm.getSampleSize(dm.tt), dm.getLociNumber(dm.tt), 
                           generate_jsfs = FALSE, generate_seg_sites = TRUE)
  checkEquals( 1, length(sum.stats) )
  checkTrue( is.list(sum.stats$seg.sites) )
  checkEquals( dm.getLociNumber(dm.tt), length(sum.stats$seg.sites) )

  sum.stats <- parseOutput(ms.file, dm.getSampleSize(dm.tt), dm.getLociNumber(dm.tt), 
                           generate_jsfs = TRUE, generate_seg_sites = TRUE)
  checkEquals( 2, length(sum.stats) )
  checkTrue( is.list(sum.stats$seg.sites) )
  checkTrue( is.matrix(sum.stats$jsfs) )
  checkTrue( sum(sum.stats$jsfs) > 0 )

  sum.stats <- parseOutput(ms.file, dm.getSampleSize(dm.tt), dm.getLociNumber(dm.tt), 
                           generate_jsfs = TRUE, generate_seg_sites = FALSE)
  checkEquals( 1, length(sum.stats) )
  checkTrue( is.matrix(sum.stats$jsfs) )
}

test.parseMsPositions <- function() {
  positions <- rep(0, 10)
  positions <- parseMsPositions('positions: 0.0010 0.0474 0.3171')
  checkEquals( '0.0010', positions[1] )
  checkEquals( '0.0474', positions[2] )
  checkEquals( '0.3171', positions[3] )
  checkEquals( 3, length(positions) )
  checkEquals( 5, length(parseMsPositions('positions: 0.1 0.2 0.3 0.4 0.5')) )
  checkEquals( 1, length(parseMsPositions('positions: 0.1')) )

  checkException(parseMsPositions('0.1 0.2 0.3'), s = TRUE)
  checkException(parseMsPositions(' '), s = TRUE)
  checkException(parseMsPositions('segsites: 0'), s = TRUE)
}
