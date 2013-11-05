test.createInitialBlocks <- function() {
  par.ranges <- matrix(1:4, 2, 2)
  checkEquals(length(createInitialBlocks(par.ranges, 1)), 1)
  checkEquals(length(createInitialBlocks(par.ranges, 2)), 4)
  checkEquals(length(createInitialBlocks(par.ranges, 3)), 9)
  checkEquals(length(createInitialBlocks(par.ranges, 4)), 16)

  par.ranges <- matrix(1:6, 3, 2)
  checkEquals(length(createInitialBlocks(par.ranges, 1)), 1)
  checkEquals(length(createInitialBlocks(par.ranges, 2)), 8)
  checkEquals(length(createInitialBlocks(par.ranges, 3)), 27)
  checkEquals(length(createInitialBlocks(par.ranges, 4)), 64)
}
