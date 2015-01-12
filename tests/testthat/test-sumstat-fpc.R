context("FGC Summary Statistic")

test_that("Fpc Breaks calculation works", {
  fpc = Stat_FPC$new(sum.stats.tt$seg.sites, dm.tt, population = 1, group = 0)
  expect_false(is.null(fpc$get_breaks()))
  expect_false(is.null(fpc$get_breaks()$near))
  expect_false(is.null(fpc$get_breaks()$far))
  expect_false(is.null(fpc$get_breaks()$mut))
  rm(fpc)
  
  fpc = Stat_FPC$new(sum.stats.tt$seg.sites, dm.tt, population = 2, group = 0)
  expect_false(is.null(fpc$get_breaks()))
  expect_false(is.null(fpc$get_breaks()$near))
  expect_false(is.null(fpc$get_breaks()$far))
  expect_false(is.null(fpc$get_breaks()$mut))
  rm(fpc)
  
  fpc = Stat_FPC$new(sum.stats.grp$seg.sites.1, dm.grp, 
                     population = 1, group = 1)
  expect_false(is.null(fpc$get_breaks()))
  expect_false(is.null(fpc$get_breaks()$near))
  expect_false(is.null(fpc$get_breaks()$far))
  expect_false(is.null(fpc$get_breaks()$mut))
  rm(fpc)  
  
  fpc = Stat_FPC$new(sum.stats.grp$seg.sites.2, dm.grp, 
                     population = 2, group = 2)
  expect_false(is.null(fpc$get_breaks()))
  expect_false(is.null(fpc$get_breaks()$near))
  expect_false(is.null(fpc$get_breaks()$far))
  expect_false(is.null(fpc$get_breaks()$mut))
  rm(fpc)
})

test_that("generateFpcStat works", {
  # One Locus
  seg.sites <- list(matrix(c(1, 1, 0, 0, 0, 
                             1, 0, 1, 0, 1, 
                             1, 1, 1, 1, 0, 
                             0, 1, 1, 0, 1, 
                             0, 0, 0, 0, 1, 
                             1, 0, 0, 0, 0), 6, byrow=TRUE))
  attr(seg.sites[[1]], "positions") <- c(0.1, 0.12, 0.5, 0.51, 0.61)
  dm <- dm.createDemographicModel(c(6,0), 1)
  fpc = Stat_FPC$new(seg.sites, dm, population = 1) 
  
  breaks = fpc$generate(seg.sites, breaks = list(near=.5, far=.5, mut=.5))  
  expect_equal(sum(abs(breaks)), 1)
  expect_equal(breaks[2, 2, 1], 1)  
  
  breaks = fpc$generate(seg.sites, breaks = list(near=.5, far=.9, mut=.5))  
  expect_equal(sum(abs(breaks)), 1)
  expect_equal(breaks[2, 1, 1], 1)
  
  # Check that NaNs are ignored
  attr(seg.sites[[1]], "positions") <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  stat = fpc$generate(seg.sites, breaks = list(near=.5, far=.9, mut=.5))  
  expect_equal(sum(abs(stat)), 0) 
  
  # Two Loci
  attr(seg.sites[[1]], "positions") <- c(0.1, 0.12, 0.5, 0.51, 0.61)
  seg.sites[[2]] <- matrix(c(1, 1, 0, 0, 0, 
                             1, 0, 1, 0, 1, 
                             1, 1, 1, 1, 0, 
                             0, 1, 1, 0, 1, 
                             0, 0, 0, 0, 1, 
                             1, 0, 0, 0, 0), 6, byrow=TRUE)
  attr(seg.sites[[2]], "positions") <- c(0.1, 0.3, 0.5, 0.7, 0.9)
  dm <- dm.setLociNumber(dm, 2)
  fpc = Stat_FPC$new(seg.sites, dm, population = 1) 
  
  stat = fpc$generate(seg.sites, breaks = list(near=.5, far=c(.25, .85), mut=.5))  
  expect_equal(sum(abs(stat)), 1)
  expect_equal(stat[2, 2, 1], 1)

  attr(seg.sites[[2]], "positions") <- c(0.1, 0.11, 0.12, 0.13, 0.14)
  stat = fpc$generate(seg.sites, breaks = list(near=.5, far=c(.25, .85), mut=.5))  
  expect_equal(sum(abs(stat)), 1)
  expect_equal(stat[2, 2, 1], 1)

  attr(seg.sites[[2]], "positions") <- c(0.1, 0.11, 0.12, 0.6, 0.7)
  stat = fpc$generate(seg.sites, breaks = list(near=.5, far=c(.25, .85), mut=.5))  
  expect_equal(sum(abs(stat)), 2)
  expect_equal(stat[2, 2, 1], 2)
  
  attr(seg.sites[[2]], "positions") <- c(0.1, 0.11, 0.12, 0.6, 0.7)
  stat = fpc$generate(seg.sites, breaks = list(near=.5, far=c(.25, .7), mut=.001))  
  expect_equal(sum(abs(stat)), 2)
  expect_equal(stat[2, 2, 2], 1) 
  expect_equal(stat[2, 3, 2], 1)   
  
  seg.sites[[3]] <- matrix(0, 6, 0)
  attr(seg.sites[[3]], "positions") <- numeric()
  dm <- dm.setLociNumber(dm, 3)  
  fpc = Stat_FPC$new(seg.sites, dm, population = 1) 
  stat = fpc$generate(seg.sites, breaks = list(near=.5, far=c(.25, .7), mut=.001))
  expect_equal(sum(abs(stat)), 2)
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
  attr(seg_sites[[1]], 'positions') <- 1:5/50
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length) 
  expect_equal(fpc_violations[1, ], c(mid_near=0.5, mid_far=NaN, outer=NaN, 
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
