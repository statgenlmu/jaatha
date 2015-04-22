context("Likelihood Table")

test_that("Sorting a likelihood table works correctly", {
  lt <- matrix(c(-Inf, 3, 2, 1, 2, 3), 3, 2)
  expect_equal(sort_likelihood_table(lt),
               matrix(c(3, 2, 2, 3), 2, 2))
  expect_equal(sort_likelihood_table(lt, 2),
               matrix(c(3, 2, 2, 3), 2, 2))
  expect_equal(sort_likelihood_table(lt, 1),
               matrix(c(3, 2), 1, 2))
  expect_equal(sort_likelihood_table(lt, 5),
               matrix(c(3, 2, 2, 3), 2, 2))
  
  lt <- matrix(c(-Inf, -Inf, -Inf, 1, 2, 3), 3, 2)
  expect_equal(sort_likelihood_table(lt),
               matrix(0, 0, 2))
})
