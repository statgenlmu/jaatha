test.parseOutput.ms <- function() {
  set.seed(25)
  dm.tt@sum.stats <- data.frame()
  dm.tt <- dm.addSummaryStatistic(dm.tt, 'file')
  dm.tt <- dm.addSummaryStatistic(dm.tt, 'seg.sites')
  ss <- dm.getSampleSize(dm.tt)
  ln <- dm.getLociNumber(dm.tt)
  breaks <- c(.25, .5, .75)
  ms.file <- dm.simSumStats(dm.tt, c(1,5))$file
  checkException(parseOutput('bulb.txt', ss, ln), s=TRUE)
  sum.stats <- parseOutput(ms.file, ss, ln, 0, FALSE, TRUE, FALSE)
  checkEquals( 1, length(sum.stats) )
  checkTrue( is.list(sum.stats$seg.sites) )
  checkEquals( dm.getLociNumber(dm.tt), length(sum.stats$seg.sites) )

  sum.stats <- parseOutput(ms.file, ss, ln, 0, TRUE, TRUE, FALSE)
  checkEquals( 2, length(sum.stats) )
  checkTrue( is.list(sum.stats$seg.sites) )
  checkTrue( is.matrix(sum.stats$jsfs) )
  checkTrue( sum(sum.stats$jsfs) > 0 )

  checkException( parseOutput(ms.file, ss, ln, 0, TRUE, TRUE, TRUE), s=TRUE )
  sum.stats <- parseOutput(ms.file, ss, ln, 0, TRUE, TRUE, TRUE, breaks, breaks)
  checkEquals( 3, length(sum.stats) )
  checkTrue( is.list(sum.stats$seg.sites) )
  checkTrue( is.matrix(sum.stats$jsfs) )
  checkTrue( sum(sum.stats$jsfs) > 0 )
  checkTrue( is.matrix(sum.stats$fpc) )
  checkEquals( c(5,5), dim(sum.stats$fpc) )
  checkTrue( sum(sum.stats$fpc) > 0 )

  sum.stats <- parseOutput(ms.file, ss, ln, 0, FALSE, FALSE, TRUE, c(breaks, .9), breaks)
  checkEquals( 1, length(sum.stats) )
  checkTrue( is.matrix(sum.stats$fpc) )
  checkEquals( c(6,5), dim(sum.stats$fpc) )
  checkTrue( sum(sum.stats$fpc) > 0 )

  sum.stats <- parseOutput(ms.file, ss, ln, 0, FALSE, FALSE, FALSE)
  checkEquals( 0, length(sum.stats) )
}

test.parseMsPositions <- function() {
  positions <- rep(0, 10)
  positions <- parseMsPositions('positions: 0.0010 0.0474 0.3171')
  checkEquals( 0.0010, positions[1] )
  checkEquals( 0.0474, positions[2] )
  checkEquals( 0.3171, positions[3] )
  checkEquals( 3, length(positions) )
  checkEquals( 5, length(parseMsPositions('positions: 0.1 0.2 0.3 0.4 0.5')) )
  checkEquals( 1, length(parseMsPositions('positions: 0.1')) )

  checkException(parseMsPositions('0.1 0.2 0.3'), s = TRUE)
  checkException(parseMsPositions(' '), s = TRUE)
  checkException(parseMsPositions('segsites: 0'), s = TRUE)
}

