context("Initial Search")

test_that("test.createInitialBlocks", {
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

