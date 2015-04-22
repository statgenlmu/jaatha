context('Refined Search')

test_that('refined search works', {
  tmp <- tempfile()

  sink(tmp)
  jaatha <- Jaatha.initialSearch(jaatha.csi, 10, 2)
  jaatha <- Jaatha.refinedSearch(jaatha, 2, sim = 20, 
                                 sim.final = 10, max.steps = 10)
  sink(NULL)

  expect_true(is.matrix(jaatha@likelihoods_rs))
  expect_true(ncol(jaatha@likelihoods_rs) == 4)
  expect_true(nrow(jaatha@likelihoods_rs) >= 10)
  
  unlink(tmp)
})


test_that('refined search works with smoothing', {
  tmp <- tempfile()
  sink(tmp)
  smooth_jaatha <- Jaatha.initialSearch(smooth_jaatha, 25, 2)
  sink(NULL)
  pStartPoints <- Jaatha.getLikelihoods(smooth_jaatha, initial_search = TRUE)
  expect_equal(nrow(pStartPoints), 4)
  sink(tmp)
  smooth_jaatha <- Jaatha.refinedSearch(smooth_jaatha, 1, 25, 25, max.steps = 5)
  sink(NULL)
  expect_true(ncol(smooth_jaatha@likelihoods_rs) == 4)
  expect_true(nrow(smooth_jaatha@likelihoods_rs) == 5)

  unlink(tmp)
})
