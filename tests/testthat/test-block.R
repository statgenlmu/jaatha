context("Block")

border <- matrix(c(0.4, 0.4, 0.6, 0.6), 2, 2)
border_1dim <- matrix(c(0, 1), 1, 2)

test_that("blocks can be initalized", {
  block <- block_class$new(border)
  expect_equivalent(block$get_border(), border)
  
  expect_error(block_class$new(t(border)))
  expect_error(block_class$new(border[ , 2:1]))
})


test_that("test if a block includes a point works", {
  block <- block_class$new(border)
  expect_equal(block$includes(c(0.5, 0.5)), TRUE)
  expect_equal(block$includes(c(0.4, 0.5)), TRUE)
  expect_equal(block$includes(c(0.6, 0.5)), TRUE)
  expect_equal(block$includes(c(0.6, 0.6)), TRUE)
  expect_equal(block$includes(c(0.3, 0.6)), FALSE)
  expect_equal(block$includes(c(0.3, 0.3)), FALSE)
  expect_equal(block$includes(c(0.4, 0.61)), FALSE)
  
  expect_error(block$includes(0.5))
})


test_that("block returns its middle", {
  block <- block_class$new(border)
  expect_equivalent(block$get_middle(), c(.5, .5))
  
  block <- block_class$new(border_1dim)
  expect_equivalent(block$get_middle(), .5)
})


test_that("block returns its corners", {
  block <- block_class$new(border)
  expect_equivalent(block$get_corners(), matrix(c(.4, .6, .4, .6,
                                                  .4, .4, .6, .6), 4, 2))
  
  block <- block_class$new(border_1dim)
  expect_equivalent(block$get_corners(), matrix(c(0, 1), 1, 2))
})
