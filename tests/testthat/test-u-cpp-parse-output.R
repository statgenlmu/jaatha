context("Rcpp output parsing")

test_that("parseMsPositions works", {
  positions <- rep(0, 10)
  positions <- parseMsPositions("positions: 0.0010 0.0474 0.3171")
  expect_equal(positions[1], 0.001)
  expect_equal(positions[2], 0.0474)
  expect_equal(positions[3], 0.3171)
  expect_equal(length(positions), 3)
  expect_equal(length(parseMsPositions("positions: 0.1 0.2 0.3 0.4 0.5")), 5)
  expect_equal(length(parseMsPositions("positions: 0.1")), 1)
  expect_error(parseMsPositions("0.1 0.2 0.3"))
  expect_error(parseMsPositions(" "))
  expect_error(parseMsPositions("segsites: 0"))
})

test_that("parseOutput works for ms", {
  set.seed(25)
  dm.tt@sum.stats <- data.frame()
  dm.tt <- dm.addSummaryStatistic(dm.tt, "file")
  ss <- dm.getSampleSize(dm.tt)
  ln <- dm.getLociNumber(dm.tt)

  ms.file <- dm.simSumStats(dm.tt, c(1, 5))$file
  expect_error(parseOutput("bulb.txt", ss, ln))
  
  seg_sites <- parseOutput(ms.file, ss, ln, 0)
  expect_true(is.list(seg_sites))
  expect_equal(length(seg_sites), dm.getLociNumber(dm.tt))
  
  ms.file <- c(ms.file, dm.simSumStats(dm.tt, c(1, 5))$file)
  seg_sites <- parseOutput(ms.file, ss, 2*ln, 0)
  expect_true(is.list(seg_sites))
  expect_equal(length(seg_sites), 2*dm.getLociNumber(dm.tt))
  
  unlink(ms.file)
})

test_that("parseOutput work for seq-gen", {
  if (!test_seqgen) return()
  invisible(files <- dm.simSumStats(dm.hky, c(1, 5))$file)
  seg_sites <- parseOutput(files[['seqgen']], dm.getSampleSize(dm.hky), 
                           dm.getLociNumber(dm.hky), 1)
  
  expect_true(is.list(seg_sites))
  expect_equal(length(seg_sites), dm.getLociNumber(dm.hky))
  
  seg_sites <- parseOutput(c(files[['seqgen']], files[['seqgen']]), 
                           dm.getSampleSize(dm.hky), 
                           2*dm.getLociNumber(dm.hky), 1)
  expect_true(is.list(seg_sites))
  expect_equal(length(seg_sites), 2*dm.getLociNumber(dm.hky))

  unlink(files)
})