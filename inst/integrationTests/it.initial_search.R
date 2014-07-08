test.initialSearch.normal <- function() {
  jaatha <- Jaatha.initialSearch(jaatha.tt, sim=10, blocks.per.par=2)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  checkEquals(4, nrow(pStartPoints))

  # Test reproducibility
  jaatha <- Jaatha.initialSearch(jaatha.tt, sim=10, blocks.per.par=2)
  pStartPoints2 <- Jaatha.getStartingPoints(jaatha)
  checkEquals(pStartPoints, pStartPoints2)

  # Parallel version
  jaatha.tt@cores <- 2 
  jaatha <- Jaatha.initialSearch(jaatha.tt, sim=10, blocks.per.par=2)
  pStartPoints3 <- Jaatha.getStartingPoints(jaatha)
  checkEquals(pStartPoints, pStartPoints3)
}

test.initialSearch.seqgen <- function() {
  dm.sq <- dm.setMutationModel(dm.tt, "HKY", c(.2, .2, .2, .4), 0.5)
  jaatha <- Jaatha.initialize(dm.sq, sum.stats.tt, seed=24)
  jaatha <- Jaatha.initialSearch(jaatha, sim=11, blocks.per.par=2)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  checkEquals(4, nrow(pStartPoints))
}

test.initialSearch.folded <- function() {
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, folded=T, seed=30)
  jaatha <- Jaatha.initialSearch(jaatha, sim=10, blocks.per.par=2)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  checkEquals(4, nrow(pStartPoints))
}

test.initialSearch.smoothing <- function() {
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 30, smoothing=TRUE)
  jaatha <- Jaatha.initialSearch(jaatha, sim=50, blocks.per.par=1)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  #checkEquals(1, nrow(pStartPoints))
}

test.initialSearch.fpc <- function() {
  jaatha <- Jaatha.initialize(dm.fpc, sum.stats.fpc, 1234, cores=2)
  jaatha <- Jaatha.initialSearch(jaatha, sim=20, blocks.per.par=2)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  checkEquals(6, nrow(pStartPoints))
}

test.initialSearch.fpc_groups <- function() {
  set.seed(1234)
  dm <- dm.addSampleSize(dm.fpc, 11:12, group = 2)
  dm <- dm.addSampleSize(dm, 5:6, group = 3)
  sum.stats <- dm.simSumStats(dm.addSummaryStatistic(dm, 'seg.sites'), c(1, 2, 5))
  jaatha.fpc <- Jaatha.initialize(dm, sum.stats, 123)
  jaatha.fpc <- Jaatha.initialSearch(jaatha.fpc, sim=20, blocks.per.par=2)
}