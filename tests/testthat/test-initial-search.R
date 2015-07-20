context("Initial Search")

test_that("determining block per par works", {
  expect_equal(determine_bpp(2, 1), 1)
  expect_equal(determine_bpp(2, 2), 2)
  expect_equal(determine_bpp(2, 3), 2)
  expect_equal(determine_bpp(2, 4), 2)
  expect_equal(determine_bpp(2, 5), 3)
  expect_equal(determine_bpp(2, 9), 5)
  
  expect_equal(determine_bpp(3, 1), 1)
  expect_equal(determine_bpp(3, 2), 2)
  expect_equal(determine_bpp(3, 3), 2)
  expect_equal(determine_bpp(3, 4), 2)
  expect_equal(determine_bpp(3, 5), 2)
  expect_equal(determine_bpp(3, 9), 3)
})


test_that("creation of initial blocks works", {
  par_ranges <- par_ranges_class$new(matrix(1:4, 2, 2))
  expect_equal(1, length(create_initial_blocks(par_ranges, 1)))
  expect_equal(4, length(create_initial_blocks(par_ranges, 2)))
  expect_equal(6, length(create_initial_blocks(par_ranges, 3)))
  expect_equal(8, length(create_initial_blocks(par_ranges, 4)))
  par_ranges <- par_ranges_class$new(matrix(1:6, 3, 2))
  expect_equal(1, length(create_initial_blocks(par_ranges, 1)))
  expect_equal(6, length(create_initial_blocks(par_ranges, 2)))
  expect_equal(9, length(create_initial_blocks(par_ranges, 3)))
  expect_equal(12, length(create_initial_blocks(par_ranges, 4)))
})


