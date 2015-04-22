context("Initial Search")

test_that("CSI complete run", {
  jaatha <- Jaatha.initialSearch(jaatha.csi, 20, 2)
  pStartPoints <- Jaatha.getLikelihoods(jaatha, initial_search = TRUE)
  expect_equal(nrow(pStartPoints), 4)
  jaatha@cores <- 1
  jaatha <- Jaatha.initialSearch(jaatha.csi, 20, 2)
  expect_equal(pStartPoints, 
               Jaatha.getLikelihoods(jaatha, initial_search = TRUE))
  
  jaatha@cores <- 2
  jaatha <- Jaatha.initialSearch(jaatha.csi, 40, 1)
  expect_equal(nrow(Jaatha.getLikelihoods(jaatha, initial_search = TRUE)), 1)
  
  jaatha <- Jaatha.initialSearch(jaatha.csi, 20, 3)
  expect_equal(nrow(Jaatha.getLikelihoods(jaatha, initial_search = TRUE)), 6)
  
  jaatha <- Jaatha.refinedSearch(jaatha, 2, sim = 10, sim.final = 10)
  expect_true(is.matrix(jaatha@likelihood.table))
  expect_true(ncol(jaatha@likelihood.table) == 4)
  expect_true(nrow(jaatha@likelihood.table) >= 10)
  jaatha@cores <- 1
  jaatha2 <- Jaatha.refinedSearch(jaatha, 2, sim = 10, sim.final = 10)
  expect_equal(jaatha2@likelihood.table, jaatha@likelihood.table)
})

test_that("test.initialSearch.normal", {
  jaatha <- Jaatha.initialSearch(jaatha.tt, sim = 10, blocks.per.par = 1)
  pStartPoints <- Jaatha.getLikelihoods(jaatha, initial_search = TRUE)
  expect_equal(nrow(pStartPoints), 1)
})

test_that("test.initialSearch.folded", {
  jaatha <- Jaatha.initialize(sum.stats.tt, dm.tt, folded = TRUE)
  jaatha <- Jaatha.initialSearch(jaatha, sim = 20, blocks.per.par = 1)
  pStartPoints <- Jaatha.getLikelihoods(jaatha, initial_search = TRUE)
  expect_equal(nrow(pStartPoints), 1)
})

test_that("test.initialSearch.fpc", {
  jaatha <- Jaatha.initialize(sum.stats.mig, dm.mig, cores=2, use_fpc=TRUE)
  jaatha <- Jaatha.initialSearch(jaatha, sim = 20, blocks.per.par = 1)
  pStartPoints <- Jaatha.getLikelihoods(jaatha, initial_search = TRUE)
  expect_equal(nrow(pStartPoints), 1)
})

test_that("test.initialSearch.fpc_groups", {
  jaatha <- Jaatha.initialize(sum.stats.grp, dm.grp, cores=2, use_fpc=TRUE)
  jaatha <- Jaatha.initialSearch(jaatha, sim = 10, blocks.per.par = 1)
  pStartPoints <- Jaatha.getLikelihoods(jaatha, initial_search = TRUE)
  expect_equal(nrow(pStartPoints), 1)
})
