context("FGC Summary Statistic")

test_that("test.calcFpcBreaks", {
    dm = calcFpcBreaks(dm.fpc, seg.sites)
    expect_false(is.null(dm@options[["fpc.breaks.near"]]))
    expect_false(is.null(dm@options[["fpc.breaks.far"]]))
    dm = calcFpcBreaks(dm.fpc, seg.sites, group = 1)
    expect_false(is.null(dm@options[["group.1"]][["fpc.breaks.near"]]))
    expect_false(is.null(dm@options[["group.1"]][["fpc.breaks.far"]]))
    dm = calcFpcBreaks(dm, seg.sites, group = 2)
    expect_false(is.null(dm@options[["group.1"]][["fpc.breaks.near"]]))
    expect_false(is.null(dm@options[["group.1"]][["fpc.breaks.far"]]))
    expect_false(is.null(dm@options[["group.2"]][["fpc.breaks.near"]]))
    expect_false(is.null(dm@options[["group.2"]][["fpc.breaks.far"]]))
})

test_that("generateFpcStat works", {
  seg.sites <- list(matrix(c(1, 1, 0, 0, 0, 
                             1, 0, 1, 0, 1, 
                             1, 1, 1, 1, 0, 
                             0, 1, 1, 0, 1, 
                             0, 0, 0, 0, 1, 
                             1, 0, 0, 0, 0), 6, byrow=TRUE))
  colnames(seg.sites[[1]]) <- c(0.1, 0.12, 0.5, 0.51, 0.61)
  dm <- dm.createDemographicModel(c(3,3), 1)
  dm@options[["fpc.breaks.near"]] <- .5
  dm@options[["fpc.breaks.far"]] <- .5
    
  fpc <- generateFpcStat(seg.sites, dm)
  expect_equal(sum(abs(fpc)), 1)
  expect_equal(fpc[2, 2], 1)  
  
  dm@options[["fpc.breaks.far"]] <- .9
  fpc <- generateFpcStat(seg.sites, dm)
  expect_equal(sum(abs(fpc)), 1)
  expect_equal(fpc[2, 1], 1) 
  
  dm@options[["fpc.breaks.far"]] <- c(.25, .85)
  seg.sites[[2]] <- matrix(c(1, 1, 0, 0, 0, 
                             1, 0, 1, 0, 1, 
                             1, 1, 1, 1, 0, 
                             0, 1, 1, 0, 1, 
                             0, 0, 0, 0, 1, 
                             1, 0, 0, 0, 0), 6, byrow=TRUE)
  colnames(seg.sites[[2]]) <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  fpc <- generateFpcStat(seg.sites, dm)
  expect_equal(sum(abs(fpc)), 2)
  expect_equal(fpc[2, 2], 1)
  expect_equal(fpc[3, 2], 1)

  colnames(seg.sites[[2]]) <- c(0.1, 0.11, 0.12, 0.13, 0.14)
  fpc <- generateFpcStat(seg.sites, dm)
  expect_equal(sum(abs(fpc)), 2)
  expect_equal(fpc[2, 2], 1)
  expect_equal(fpc[2, 4], 1)
  
  seg.sites[[3]] <- matrix(0, 6, 0)
  fpc <- generateFpcStat(seg.sites, dm)
  expect_equal(sum(abs(fpc)), 3)
  expect_equal(fpc[3, 4], 1)
  
  seg.sites[[4]] <- seg.sites[[1]]
  fpc <- generateFpcStat(seg.sites, dm)
  expect_equal(sum(abs(fpc)), 4)
  expect_equal(fpc[2, 2], 2)
})

test_that("calcPercentFpcViolation works", {
  snp.matrix <- matrix(c(1, 1, 0, 0, 0, 
                         1, 0, 1, 0, 1, 
                         1, 1, 1, 1, 0, 
                         0, 1, 1, 0, 1, 
                         0, 0, 0, 0, 1), 5)
  
  positions <- c(0.1, 0.12, 0.5, 0.51, 0.61)
  expect_equal(calcPercentFpcViolation(snp.matrix, positions), 
               c(0.5, 0.5))
  positions[4:5] <- c(0.7, 0.75)
  expect_equal(calcPercentFpcViolation(snp.matrix, positions), 
               c(1, 0.4))
  positions <- 1:5/5
  expect_equal(calcPercentFpcViolation(snp.matrix, positions), 
               c(NA, 0.5))
  snp.matrix[1, ] <- 1
  snp.matrix[-1, ] <- 0
  expect_true(all(is.na(calcPercentFpcViolation(snp.matrix, 
                                                positions))))
  expect_true(all(is.na(calcPercentFpcViolation(matrix(0, 5, 
                                                       0), numeric(0)))))
})

test_that("countClasses works", {
  classes <- rbind(c(1, 1, 2, 2, NA, NA),
                   c(1, 2, 1, NA, NA, NA))
  dimension <- c(3, 3)
  expect_equal(countClasses(classes, dimension), 
               array(c(1, 1, 0, 1, 0, 0, 0, 1, 2), c(3,3)))
  
  classes <- rbind(c(1, 1, 1, 1, 2, 2),
                   c(1, 1, 2, 2, 1, 1),
                   c(1, 1, 1, 2, 1, 2))
  expect_equal(countClasses(classes, c(2,2,2)), 
               array(c(2, 1, 1, 0, 0, 1, 1, 0), c(2,2,2)))
})
