context('Refined Search')

test_that('refined search works', {
  tmp <- tempfile()

  sink(tmp)
  jaatha <- Jaatha.initialSearch(jaatha.csi, 10, 2)
  jaatha <- Jaatha.refinedSearch(jaatha, 2, sim = 20, sim.final = 10)
  sink(NULL)

  expect_true(is.matrix(jaatha@likelihood.table))
  expect_true(ncol(jaatha@likelihood.table) == 4)
  expect_true(nrow(jaatha@likelihood.table) >= 10)
  
  unlink(tmp)
})


test_that('refined search works with smoothing', {
  tmp <- tempfile()
  sink(tmp)
  smooth_jaatha <- Jaatha.initialSearch(smooth_jaatha, 25, 2)
  pStartPoints <- Jaatha.getStartingPoints(smooth_jaatha)
  expect_equal(nrow(pStartPoints), 4)
  smooth_jaatha <- Jaatha.refinedSearch(smooth_jaatha, 1, 25, 25, max.steps = 5)
  expect_true(ncol(smooth_jaatha@likelihood.table) == 4)
  expect_true(nrow(smooth_jaatha@likelihood.table) == 4)
  sink(NULL)
  unlink(tmp)
})
