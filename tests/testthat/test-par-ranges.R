context("Parameter Ranges")

test_that("Ranges can be initialized", {
  par_range_class$new(matrix(1:4, 2, 2))
  par_range_class$new(matrix(1:6, 3, 2))
  par_range_class$new(matrix(-3:2, 3, 2))
  expect_error(par_range_class$new(1:6))
  expect_error(par_range_class$new(matrix(1:6, 2, 3)))
  expect_error(par_range_class$new(matrix(6:1, 3, 2)))
})

test_that("Parameter normalization works", {
  par_range <- par_range_class$new(matrix(-3:2, 3, 2))
  expect_true(all(par_range$normalize(-3:-1) == c(0, 0, 0)))
  expect_true(all(par_range$normalize(0:2) == c(1, 1, 1)))
  expect_true(all(par_range$normalize(-2:0) > 0))
  expect_true(all(par_range$normalize(-2:0) < 1))
})

test_that("Parameter denormalization works", {
  par_range <- par_range_class$new(matrix(-3:2, 3, 2))
  expect_true(all(par_range$denormalize(c(0, 0, 0)) == c(0, 0, 0)))
  expect_true(all(par_range$denormalize(c(1, 1, 1)) == c(1, 1, 1)))
  expect_true(sum(abs(denormalize(c(0, 0), jaatha.csi) - c(0.1, 0.1))) < 1e-11)
  expect_true(sum(abs(denormalize(c(1, 1), jaatha.csi) - c(10, 10))) < 1e-11)
})


test_that("Normalization and Denormalization are inverse", {
    for (x in 0:10/10) {
        pars <- rep(x, 2)
        expect_true(sum(abs(normalize(denormalize(x, jaatha.csi), 
            jaatha.csi) - pars)) < 1e-11)
    }
})

