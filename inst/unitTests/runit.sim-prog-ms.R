dm <- dm.createThetaTauModel(c(10, 11), 7)

snp.matrix <- matrix(c(1,1,0,0,0,1,0,1,0,1,1,1,1,1,0,0,1,1,0,1,0,0,0,0,1), 5) 
colnames(snp.matrix) <- c('0.1', '0.12', '0.5', '0.51', '0.61')

test.callMS <- function() {
  checkException(callMs())
  checkTrue(file.exists(callMs("-t 5", dm)))
}

test.msSingleSimFunc <- function() {
  set.seed(789)
  sum.stats <- msSingleSimFunc(dm.tt, c(1,10))
  checkTrue(is.list(sum.stats))
  checkEquals(2, length(sum.stats))

  checkTrue(!is.null(sum.stats$pars))
  checkEquals(c(1,10), sum.stats$pars)

  checkTrue(!is.null(sum.stats$jsfs))
  checkTrue(is.array(sum.stats$jsfs))
  checkEquals(c(12,13), dim(sum.stats$jsfs))
  checkTrue(sum(sum.stats$jsfs) > 0)

  set.seed(456)
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

  dm@sum.stats <- c('file', 'seg.sites')
  sum.stats <- msSingleSimFunc(dm, c(1, .5))
  checkEquals(3, length(sum.stats))
  checkTrue(!is.null(sum.stats$pars))
  checkTrue(!is.null(sum.stats$seg.sites))
  checkTrue(!is.null(sum.stats$file))
  checkTrue(file.exists(sum.stats$file))
  checkTrue((file.info(sum.stats$file)$size > 0))

  set.seed(123)
  dm@sum.stats <- 'fpc'
  dm <- calcFpcBreaks(dm, sum.stats$seg.sites) 
  sum.stats <- msSingleSimFunc(dm, c(1, 10))
  checkEquals(2, length(sum.stats))
  checkTrue(!is.null(sum.stats$pars))
  checkTrue(is.array(sum.stats[['fpc']]))
  checkEquals(dm.getLociNumber(dm), sum(sum.stats[['fpc']]))

  dm <- dm.addSummaryStatistic(dm, 'seg.sites')
  sum.stats <- msSingleSimFunc(dm, c(1, 0.1))
  checkEquals(3, length(sum.stats))
  checkTrue(!is.null(sum.stats$pars))
  checkTrue(!is.null(sum.stats$seg.sites))
  checkTrue(is.array(sum.stats[['fpc']]))
  checkEquals(dm.getLociNumber(dm), sum(sum.stats[['fpc']]))
}

test.simProg <- function(){
  checkTrue(!is.null(.jaatha$sim_progs[["ms"]]))
}
