dm <- dm.createThetaTauModel(c(10, 11), 7)

snp.matrix <- matrix(c(1,1,0,0,0,1,0,1,0,1,1,1,1,1,0,0,1,1,0,1,0,0,0,0,1), 5) 
colnames(snp.matrix) <- c('0.1', '0.12', '0.5', '0.51', '0.61')

test.callMS <- function() {
  checkException(callMs())
  checkTrue(file.exists(callMs("-t 5", dm)))
}

test.msSingleSimFunc <- function() {
  sum.stats <- msSingleSimFunc(dm, c(1,10))
  checkTrue(is.list(sum.stats))
  checkEquals(2, length(sum.stats))

  checkTrue(!is.null(sum.stats$pars))
  checkEquals(c(1,10), sum.stats$pars)

  checkTrue(!is.null(sum.stats$jsfs))
  checkTrue(is.array(sum.stats$jsfs))
  checkEquals(c(11,12), dim(sum.stats$jsfs))
  checkTrue(sum(sum.stats$jsfs) > 0)

  dm@sum.stats <- 'seg.sites'
  sum.stats <- msSingleSimFunc(dm, c(1, .5))
  checkEquals(2, length(sum.stats))
  checkTrue(!is.null(sum.stats$pars))

  checkTrue(!is.null(sum.stats$seg.sites))
  checkTrue( is.list(sum.stats$seg.sites))
  checkEquals(7, length(sum.stats$seg.sites))
  for (i in 1:7) {
    checkTrue( is.matrix(sum.stats$seg.sites[[i]]) )
    checkEquals(21, nrow(sum.stats$seg.sites[[i]]))
  }

  dm@sum.stats <- 'file'
  sum.stats <- msSingleSimFunc(dm, c(1, .5))
  checkEquals(2, length(sum.stats))
  checkTrue(!is.null(sum.stats$pars))
  checkTrue(!is.null(sum.stats$file))
  checkTrue(file.exists(sum.stats$file))
  checkTrue((file.info(sum.stats$file)$size > 0))

  set.seed(123)
  dm@sum.stats <- 'seg.sites'
  dm <- dm.addSummaryStatistic(dm, '4pc')
  sum.stats <- msSingleSimFunc(dm, c(1, 10))
  checkEquals(3, length(sum.stats))
  checkTrue(!is.null(sum.stats$pars))
  checkTrue(!is.null(sum.stats$seg.sites))
  checkTrue(is.matrix(sum.stats[['4pc']]))
  checkEquals(dm.getLociNumber(dm), sum(sum.stats[['4pc']]))

  dm <- dm.addSummaryStatistic(dm, '4pc')
  sum.stats <- msSingleSimFunc(dm, c(1, 0.1))
  checkEquals(3, length(sum.stats))
  checkTrue(!is.null(sum.stats$pars))
  checkTrue(!is.null(sum.stats$seg.sites))
  checkTrue(is.matrix(sum.stats[['4pc']]))
  checkEquals(dm.getLociNumber(dm), sum(sum.stats[['4pc']]))
}

test.calcFpcSumStat <- function() {
  snp.matrix2 <- snp.matrix
  colnames(snp.matrix2) <- 1:5/5

  dm <- dm.addSummaryStatistic(dm, '4pc')
  mat <- calcFpcSumStat(list(snp.matrix, snp.matrix, snp.matrix2, matrix(0, 5, 0)), dm)
  checkTrue(is.matrix(mat))
  checkEquals(4, sum(mat))
  checkEquals(2, mat['3', '3'])
  checkEquals(1, mat['0', '3'])
  checkEquals(1, mat['0', '0'])
}

test.violatesFpc <- function() {
  checkEquals( c(near=TRUE, violates=TRUE), violatesFpc(c(1,2), snp.matrix) )
  checkEquals( c(near=FALSE, violates=FALSE), violatesFpc(c(1,3), snp.matrix) )
  checkEquals( c(near=FALSE, violates=TRUE), violatesFpc(c(1,4), snp.matrix) )
  checkEquals( c(near=FALSE, violates=FALSE), violatesFpc(c(2,3), snp.matrix) )
  checkEquals( c(near=FALSE, violates=TRUE), violatesFpc(c(2,4), snp.matrix) )
  checkEquals( c(near=TRUE, violates=FALSE), violatesFpc(c(3,4), snp.matrix) )
  
  checkEquals( c(near=FALSE, violates=FALSE), violatesFpc(c(1,5), snp.matrix) )
  checkEquals( c(near=FALSE, violates=FALSE), violatesFpc(c(2,5), snp.matrix) )
  checkEquals( c(near=FALSE, violates=FALSE), violatesFpc(c(3,5), snp.matrix) )
  checkEquals( c(near=TRUE, violates=FALSE), violatesFpc(c(4,5), snp.matrix) )
}

test.calcPercentFpcViolations <- function() {
  checkEquals(c(near=0.5, far=0.5), calcPercentFpcViolations(snp.matrix))
  colnames(snp.matrix)[4:5] <- c('0.7', '0.75') 
  checkEquals(c(near=1, far=0.4), calcPercentFpcViolations(snp.matrix))

  # No near snps
  colnames(snp.matrix) <- 1:5/5
  checkEquals(c(near=NaN, far=0.5), calcPercentFpcViolations(snp.matrix))

  # Only singletons
  snp.matrix[1, ] <- 1
  snp.matrix[-1, ] <- 0
  checkEquals(c(near=NaN, far=NaN), calcPercentFpcViolations(snp.matrix))

  # No SNPs at all
  checkEquals(c(near=NaN, far=NaN), calcPercentFpcViolations(matrix(0, 5, 0)))
} 

test.simProg <- function(){
  checkTrue(!is.null(.jaatha$simProgs[["ms"]]))
}
