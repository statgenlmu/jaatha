context("Initial Search")

test_that("creation of initial blocks works", {
    par.ranges <- matrix(1:4, 2, 2)
    expect_equal(1, length(createInitialBlocks(par.ranges, 1)))
    expect_equal(4, length(createInitialBlocks(par.ranges, 2)))
    expect_equal(6, length(createInitialBlocks(par.ranges, 3)))
    expect_equal(8, length(createInitialBlocks(par.ranges, 4)))
    par.ranges <- matrix(1:6, 3, 2)
    expect_equal(1, length(createInitialBlocks(par.ranges, 1)))
    expect_equal(6, length(createInitialBlocks(par.ranges, 2)))
    expect_equal(9, length(createInitialBlocks(par.ranges, 3)))
    expect_equal(12, length(createInitialBlocks(par.ranges, 4)))
})


test_that("initial search works", {
  tmp <- tempfile()
  
  jaatha.csi@cores <- 2
  sink(tmp)
  jaatha <- Jaatha.initialSearch(jaatha.csi, 10, 2)
  sink(NULL)
  pStartPoints <- Jaatha.getLikelihoods(jaatha, initial_search = TRUE)
  expect_equal(nrow(pStartPoints), 4)
  
  sink(tmp)
  jaatha.csi@cores <- 1
  jaatha <- Jaatha.initialSearch(jaatha.csi, 10, 2)
  sink(NULL)
  expect_equal(pStartPoints, 
               Jaatha.getLikelihoods(jaatha, initial_search = TRUE))

  sink(tmp)
  jaatha <- Jaatha.initialSearch(jaatha.csi, 10, 1)
  sink(NULL)
  expect_equal(nrow(Jaatha.getLikelihoods(jaatha, initial_search = TRUE)), 1)
  
  sink(tmp)
  jaatha <- Jaatha.initialSearch(jaatha.csi, 10, 3)
  sink(NULL)
  expect_equal(nrow(Jaatha.getLikelihoods(jaatha, initial_search = TRUE)), 6)
  
  unlink(tmp)
})
