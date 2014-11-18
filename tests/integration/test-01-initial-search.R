context("Initial Search")

test_that("CSI complete run", {
  jaatha <- Jaatha.initialSearch(jaatha.csi, 20, 2)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  expect_equal(nrow(Jaatha.getStartingPoints(jaatha)), 4)
  jaatha@cores <- 1
  jaatha <- Jaatha.initialSearch(jaatha.csi, 20, 2)
  expect_equal(pStartPoints, Jaatha.getStartingPoints(jaatha))
  
  jaatha@cores <- 2
  jaatha <- Jaatha.initialSearch(jaatha.csi, 40, 1)
  expect_equal(nrow(Jaatha.getStartingPoints(jaatha)), 1)
  
  jaatha <- Jaatha.initialSearch(jaatha.csi, 20, 3)
  expect_equal(nrow(Jaatha.getStartingPoints(jaatha)), 6)
  
  jaatha <- Jaatha.refinedSearch(jaatha, 2, sim = 10, sim.final = 10)
  expect_true(is.matrix(jaatha@likelihood.table))
  expect_true(ncol(jaatha@likelihood.table) == 4)
  expect_true(nrow(jaatha@likelihood.table) >= 10)
  jaatha@cores <- 1
  jaatha2 <- Jaatha.refinedSearch(jaatha, 2, sim = 10, sim.final = 10)
  expect_equal(jaatha2@likelihood.table, jaatha@likelihood.table)
})

test_that("test.initialSearch.folded", {
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, folded = T, seed = 30)
  jaatha <- Jaatha.initialSearch(jaatha, sim = 50, blocks.per.par = 1)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  expect_equal(nrow(pStartPoints), 1)
})

test_that("test.initialSearch.fpc", {
  jaatha <- Jaatha.initialize(dm.fpc, sum.stats.fpc, 1234, cores = 2)
  jaatha <- Jaatha.initialSearch(jaatha, sim = 20, blocks.per.par = 1)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  expect_equal(nrow(pStartPoints), 1)
})

test_that("test.initialSearch.fpc_groups", {
  set.seed(1234)
  dm <- dm.addSampleSize(dm.fpc, 11:12, group = 2)
  dm <- dm.addSampleSize(dm, 5:6, group = 3)
  sum.stats <- dm.simSumStats(dm.addSummaryStatistic(dm, "seg.sites"), 
                              c(1, 2, 5))
  jaatha.fpc <- Jaatha.initialize(dm, sum.stats, 123)
  jaatha.fpc <- Jaatha.initialSearch(jaatha.fpc, sim = 10, blocks.per.par = 1)
})

test_that("test.initialSearch.normal", {
  jaatha <- Jaatha.initialSearch(jaatha.tt, sim = 10, blocks.per.par = 1)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  expect_equal(nrow(pStartPoints), 1)
})

test_that("test.initialSearch.seqgen", {
  dm.sq <- dm.setMutationModel(dm.tt, "HKY", c(0.2, 0.2, 0.2, 0.4), 0.5)
  jaatha <- Jaatha.initialize(dm.sq, sum.stats.tt, seed = 24, cores = 2)
  jaatha <- Jaatha.initialSearch(jaatha, sim = 10, blocks.per.par = 1)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  expect_equal(nrow(pStartPoints), 1)
})

test_that("test.initialSearch.seqgenAndMsms", {
  set.seed(333)
  dm.selsq <- jaatha:::dm.addPositiveSelection(dm.hky, 100, 500, population = 1, 
                                               at.time = "0.1")
  dm.selsq <- jaatha:::dm.finalize(dm.selsq)
  sum.stats <- dm.simSumStats(dm.selsq, c(1, 2, 255))
  jaatha.selsq <- Jaatha.initialize(dm.selsq, sum.stats, 123, cores = 2)
  jaatha.selsq <- Jaatha.initialSearch(jaatha.selsq, sim=10, blocks.per.par=1)
})

test_that("test.initialSearch.smoothing", {
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 30, smoothing=TRUE, cores=2)
  jaatha <- Jaatha.initialSearch(jaatha, sim = 50, blocks.per.par = 1)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
})