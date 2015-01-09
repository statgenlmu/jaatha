context("FGC Summary Statistic")

test_that("Fpc Breaks calculation works", {
    dm = calcFpcBreaks(dm.fpc, seg.sites, population = 1)
    expect_false(is.null(dm@options[["fpc_breaks_pop1"]]))
    expect_false(is.null(dm@options[["fpc_breaks_pop1"]]$near))
    expect_false(is.null(dm@options[["fpc_breaks_pop1"]]$far))
    expect_false(is.null(dm@options[["fpc_breaks_pop1"]]$mut))
    
    dm = calcFpcBreaks(dm.fpc, seg.sites, group = 1, population = 1)
    expect_false(is.null(dm@options[["group.1"]][["fpc_breaks_pop1"]]))
    expect_false(is.null(dm@options[["group.1"]][["fpc_breaks_pop1"]]$near))
    expect_false(is.null(dm@options[["group.1"]][["fpc_breaks_pop1"]]$far))
    expect_false(is.null(dm@options[["group.1"]][["fpc_breaks_pop1"]]$mut))
    
    dm = calcFpcBreaks(dm, seg.sites, group = 2, population = 1)
    expect_false(is.null(dm@options[["group.1"]][["fpc_breaks_pop1"]]))
    expect_false(is.null(dm@options[["group.2"]][["fpc_breaks_pop1"]]))
    expect_false(is.null(dm@options[["group.2"]][["fpc_breaks_pop1"]]$near))
    expect_false(is.null(dm@options[["group.2"]][["fpc_breaks_pop1"]]$far))
    expect_false(is.null(dm@options[["group.2"]][["fpc_breaks_pop1"]]$mut))
    
    dm = calcFpcBreaks(dm.fpc, seg.sites, group = 2, population = 2)
    expect_true(is.null(dm@options[["group.1"]][["fpc_breaks_pop2"]]))
    expect_false(is.null(dm@options[["group.2"]][["fpc_breaks_pop2"]]))
    expect_false(is.null(dm@options[["group.2"]][["fpc_breaks_pop2"]]$near))
    expect_false(is.null(dm@options[["group.2"]][["fpc_breaks_pop2"]]$far))
    expect_false(is.null(dm@options[["group.2"]][["fpc_breaks_pop2"]]$mut))    
    
    #dm.lt <- dm.addLocusTrio(dm.fpc, locus_length = c(200, 400, 200),
    #                         distance = c(100, 100), group = 2)
    #dm = calcFpcBreaks(dm.lt, seg.sites)
    #expect_false(is.null(dm@options[["fpc.breaks.near"]]))
    #expect_false(is.null(dm@options[["fpc.breaks.far"]]))
    #expect_false(is.null(dm@options[["fpc.breaks.between"]]))
    
    #dm = calcFpcBreaks(dm.lt, list(seg.sites[[1]]), group = 2)
    #expect_false(is.null(dm@options[['group.2']][["fpc.breaks.between"]]))
})

test_that("generateFpcStat works", {
  seg.sites <- list(matrix(c(1, 1, 0, 0, 0, 
                             1, 0, 1, 0, 1, 
                             1, 1, 1, 1, 0, 
                             0, 1, 1, 0, 1, 
                             0, 0, 0, 0, 1, 
                             1, 0, 0, 0, 0), 6, byrow=TRUE))
  attr(seg.sites[[1]], "positions") <- c(0.1, 0.12, 0.5, 0.51, 0.61)
  dm <- dm.createDemographicModel(c(6,0), 1)
  dm@options[["fpc_breaks_pop1"]] <- list(near=.5, far=.5, between=.5)
    
  fpc <- generateFpcStat(seg.sites, dm, population = 1)
  expect_equal(sum(abs(fpc)), 1)
  expect_equal(fpc[2, 2, 1], 1)  
  
  dm@options[["fpc_breaks_pop1"]]$far <- .9
  fpc <- generateFpcStat(seg.sites, dm, population = 1)
  expect_equal(sum(abs(fpc)), 1)
  expect_equal(fpc[2, 1, 1], 1) 
  
  dm@options[["fpc_breaks_pop1"]]$far <- c(.25, .85)
  seg.sites[[2]] <- matrix(c(1, 1, 0, 0, 0, 
                             1, 0, 1, 0, 1, 
                             1, 1, 1, 1, 0, 
                             0, 1, 1, 0, 1, 
                             0, 0, 0, 0, 1, 
                             1, 0, 0, 0, 0), 6, byrow=TRUE)
  attr(seg.sites[[2]], "positions") <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  dm <- dm.setLociNumber(dm, 2)
  fpc <- generateFpcStat(seg.sites, dm, population = 1)
  expect_equal(sum(abs(fpc)), 2)
  expect_equal(fpc[2, 2, 1], 2)

  attr(seg.sites[[2]], "positions") <- c(0.1, 0.11, 0.12, 0.13, 0.14)
  fpc <- generateFpcStat(seg.sites, dm, population = 1)
  expect_equal(sum(abs(fpc)), 2)
  expect_equal(fpc[2, 2, 1], 1)
  expect_equal(fpc[2, 3, 1], 1)
  
  seg.sites[[3]] <- matrix(0, 6, 0)
  attr(seg.sites[[3]], "positions") <- numeric()
  dm <- dm.setLociNumber(dm, 3)  
  fpc <- generateFpcStat(seg.sites, dm, population = 1)
  expect_equal(sum(abs(fpc)), 3)
})

test_that("calcPercentFpcViolation works", {
  seg_sites <- list(matrix(c(1, 1, 0, 0, 0, 
                             1, 0, 1, 0, 1, 
                             1, 1, 1, 1, 0, 
                             0, 1, 1, 0, 1, 
                             0, 0, 0, 0, 1), 5))
  
  attr(seg_sites[[1]], 'positions') <- c(0.1, 0.12, 0.5, 0.51, 0.61)
  locus_length <- matrix(c(0, 0, 100, 0, 0), 1, 5)
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length)  
  expect_equal(fpc_violations[1, ], c(mid_near=0.5, mid_far=0.5, outer=NaN, 
                                      between=NaN, mid=0.5, perc_polym=0.05))
  
  seg_sites[[2]] <- seg_sites[[1]]
  locus_length <- rbind(locus_length, matrix(c(0, 0, 50, 0, 0), 1, 5))
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[1, ], c(mid_near=0.5, mid_far=0.5, outer=NaN, 
                                      between=NaN, mid=0.5, perc_polym=0.05))
  expect_equal(fpc_violations[2, ], c(mid_near=0.5, mid_far=0.5, outer=NaN, 
                                      between=NaN, mid=0.5, perc_polym=0.10))
  
  attr(seg_sites[[2]], 'positions')[4:5] <- c(0.7, 0.75)
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length)  
  expect_equal(fpc_violations[2, ], c(mid_near=1, mid_far=0.4, outer=NaN, 
                                      between=NaN, mid=0.5, perc_polym=0.10))
  
  attr(seg_sites[[1]], 'positions') <- 1:5/5
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length) 
  expect_equal(fpc_violations[1, ], c(mid_near=NaN, mid_far=0.5, outer=NaN, 
                                      between=NaN, mid=0.5, perc_polym=0.05))
  
  seg_sites[[2]][1, ] <- 1
  seg_sites[[2]][-1, ] <- 0
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[2, ], c(mid_near=NA, mid_far=NA, outer=NA, 
                                      between=NA, mid=NA, perc_polym=0.10))

  
  seg_sites[[3]] <- matrix(0, 5, 0)
  attr(seg_sites[[3]], 'positions') <- numeric(0)
  locus_length <- rbind(locus_length, c(0, 0, 50, 0, 0))
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length) 
  expect_equal(fpc_violations[3, ], c(mid_near=NaN, mid_far=NaN, outer=NaN, 
                                      between=NaN, mid=NaN, perc_polym=0))
  

  # With locus-trios
  seg_sites[[4]] <- matrix(c(1, 1, 0, 0, 0, 
                             1, 0, 1, 0, 1, 
                             1, 1, 1, 1, 0, 
                             0, 1, 1, 0, 1, 
                             0, 0, 0, 0, 1), 5)
  attr(seg_sites[[4]], 'positions') <- c(0.15, 0.55, 0.05, 0.08, 0.30)
  attr(seg_sites[[4]], 'locus') <- c(-1, -1, 0, 0, 0)
  locus_length <- rbind(locus_length, c(10, 5, 6, 5, 10))
  
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length) 
  expect_equal(fpc_violations[4, ], c(mid_near=0, mid_far=NaN, outer=1, 
                                      between=0.5, mid=0, perc_polym=0.5))
  
  attr(seg_sites[[4]], 'locus') <- c(0, 0, 1, 1, 1)
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length) 
  expect_equal(fpc_violations[4, ], c(mid_near=NaN, mid_far=1, outer=0, 
                                      between=0.5, mid=1, perc_polym=1/3))
  
  attr(seg_sites[[4]], 'locus') <- c(-1, -1, 0, 0, 1)
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length) 
  expect_equal(fpc_violations[4, ], c(mid_near=0, mid_far=NaN, outer=1, 
                                      between=0.5, mid=0, perc_polym=1/3))
  
  attr(seg_sites[[4]], 'locus') <- c(-1, -1, 0, 1, 1)
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length) 
  expect_equal(fpc_violations[4, ], c(mid_near=NaN, mid_far=NaN, outer=1, 
                                      between=0, mid=NaN, perc_polym=1/6)) 
  
  
  attr(seg_sites[[4]], 'locus') <- c(-1, -1, -1, 1, 1)
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length) 
  expect_equal(fpc_violations[4, ], c(mid_near=NaN, mid_far=NaN, outer=1/3, 
                                      between=NaN, mid=NaN, perc_polym=0)) 
  
  # Not all individuals
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:2, locus_length) 
  expect_equal(fpc_violations[4, ], c(mid_near=NaN, mid_far=NaN, outer=0, 
                                      between=NaN, mid=NaN, perc_polym=0)) 
  
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:4, locus_length) 
  expect_equal(fpc_violations[4, ], c(mid_near=NaN, mid_far=NaN, outer=0.25, 
                                      between=NaN, mid=NaN, perc_polym=0)) 
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
