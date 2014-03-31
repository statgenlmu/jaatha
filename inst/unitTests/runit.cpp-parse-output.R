test.parseOutput <- function() {
  set.seed(25)
  dm.tt@sum.stats <- c('file', 'seg.sites')
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

test.addSegSitesToJsfs <- function() {
  seg.sites <- matrix(c(1,0,0,0, 1,1,0,1, 1,0,0,1), 4, 3)  
  jsfs <- matrix(0, 3, 3)
  jsfs.new <- addSegSitesToJsfs(seg.sites, c(2,2), jsfs)
  checkTrue( is.matrix(jsfs.new) )
  checkEquals( c(3,3), dim(jsfs.new) )
  checkEquals( 0, sum(jsfs) )
  checkEquals( 3, sum(jsfs.new) )
  checkEquals( 1, jsfs.new[2,1] )
  checkEquals( 1, jsfs.new[3,2] )
  checkEquals( 1, jsfs.new[2,2] )
}

test.addSegSitesToFpc <- function() {
  seg.sites <- matrix(c(1,1,0,0,0, 1,0,1,0,1, 1,1,1,1,0, 0,1,1,0,1, 0,0,0,0,1, 1,0,0,0,0), 5) 
  positions <- c(0.1, 0.12, 0.5, 0.51, 0.61)
  fpc <- matrix(0, 5, 5)
  fpc.new <- addSegSitesToFpc(seg.sites, positions, c(.25, .5, .75), c(.25, .5, .75), fpc)
  checkEquals(matrix(0, 5, 5), fpc)
  checkEquals(1, sum(fpc.new))
  checkEquals(1, fpc.new[2,2])
  fpc.new <- addSegSitesToFpc(seg.sites, positions, c(.25, .5, .75), c(.25, .5, .75), fpc.new)
  checkEquals(2, sum(fpc.new))
  checkEquals(2, fpc.new[2,2])

  fpc.new <- addSegSitesToFpc(seg.sites, positions, c(.25, .4, .75), c(.25, .5, .75), fpc.new)
  checkEquals(3, sum(fpc.new))
  checkEquals(2, fpc.new[2,2])
  checkEquals(1, fpc.new[3,2])

  fpc.new <- addSegSitesToFpc(seg.sites, positions, c(.25, .4, .75), c(.25, .4, .75), fpc.new)
  checkEquals(4, sum(fpc.new))
  checkEquals(2, fpc.new[2,2])
  checkEquals(1, fpc.new[3,2])
  checkEquals(1, fpc.new[3,3])

  singletons <- matrix(0,5,4)
  singletons[1,] <- 1
  fpc.new <- addSegSitesToFpc(singletons, positions, c(.25, .4, .75), c(.25, .4, .75), fpc)
  checkEquals(1, sum(fpc.new))
  checkEquals(1, fpc.new[5,5])

  positions <- 1:5/6
  fpc.new <- addSegSitesToFpc(seg.sites, positions, c(.25, .4, .75), c(.25, .4, .75), fpc)
  checkEquals(1, sum(fpc.new))
  checkEquals(1, fpc.new[5,3])
}

test.calcPercentFpcViolation <- function() {
  snp.matrix <- matrix(c(1,1,0,0,0,1,0,1,0,1,1,1,1,1,0,0,1,1,0,1,0,0,0,0,1), 5) 
  positions <- c(0.1, 0.12, 0.5, 0.51, 0.61)

  checkEquals(c(0.5, 0.5), calcPercentFpcViolation(snp.matrix, positions))
  positions[4:5] <- c(0.7, 0.75) 
  checkEquals(c(1, 0.4), calcPercentFpcViolation(snp.matrix, positions))

  # No near snps
  positions <- 1:5/5
  checkEquals(c(NA, 0.5), calcPercentFpcViolation(snp.matrix, positions))

  # Only singletons
  snp.matrix[1, ] <- 1
  snp.matrix[-1, ] <- 0
  checkTrue(all( is.na(calcPercentFpcViolation(snp.matrix, positions)) ))

  # No SNPs at all
  checkTrue(all( is.na(calcPercentFpcViolation(matrix(0, 5, 0), numeric(0))) ))
} 
