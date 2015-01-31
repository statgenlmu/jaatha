context("SumStat FPC")


test_that("Fpc Breaks calculation works", {
  fpc = Stat_FPC$new(sumstat_tt$seg.sites, dm_tt, population = 1, group = 0)
  expect_false(is.null(fpc$get_breaks()))
  expect_false(is.null(fpc$get_breaks()$mid_near))
  expect_false(is.null(fpc$get_breaks()$mid_far))
  expect_false(is.null(fpc$get_breaks()$perc_polym))
  rm(fpc)
  
  fpc = Stat_FPC$new(sumstat_tt$seg.sites, dm_tt, population = 2, group = 0)
  expect_false(is.null(fpc$get_breaks()))
  expect_false(is.null(fpc$get_breaks()$mid_near))
  expect_false(is.null(fpc$get_breaks()$mid_far))
  expect_false(is.null(fpc$get_breaks()$perc_polym))
  rm(fpc)
  
  fpc = Stat_FPC$new(sum_stat_grps$seg.sites.1, dm_grps, 
                     population = 1, group = 1)
  expect_false(is.null(fpc$get_breaks()))
  expect_false(is.null(fpc$get_breaks()$mid_near))
  expect_false(is.null(fpc$get_breaks()$mid_far))
  expect_false(is.null(fpc$get_breaks()$perc_polym))
  rm(fpc)  
  
  fpc = Stat_FPC$new(sum_stat_grps$seg.sites.2, dm_grps, 
                     population = 2, group = 2)
  expect_false(is.null(fpc$get_breaks()))
  expect_false(is.null(fpc$get_breaks()$mid_near))
  expect_false(is.null(fpc$get_breaks()$mid_far))
  expect_false(is.null(fpc$get_breaks()$perc_polym))
  rm(fpc)
})

test_that("generateLociCube works", {
  stat = cbind(1:6, 1:6, 1:6)
  breaks = list(1:2+.5, 3.5, c(1,3,5)+.5)
  cube = array(generateLociCube(stat, breaks, 1:3), c(3,2,4))
  expect_equal(sum(cube), 6)
  expect_equal(cube[1, 1, 1], 1)
  expect_equal(cube[2, 1, 2], 1)
  expect_equal(cube[3, 1, 2], 1)
  expect_equal(cube[3, 2, 3], 2)
  expect_equal(cube[3, 2, 4], 1)
  
  # Check that NaNs are ignored
  stat[3,1] = NaN
  stat[4,2] = NaN 
  cube = array(generateLociCube(stat, breaks, 1:3), c(3,2,4))
  expect_equal(sum(cube), 4)
  expect_equal(cube[1, 1, 1], 1)
  expect_equal(cube[2, 1, 2], 1)
  expect_equal(cube[3, 2, 3], 1)
  expect_equal(cube[3, 2, 4], 1)
  
  # No valid loci
  stat[,1] = NaN
  cube = array(generateLociCube(stat, breaks, 1:3), c(3,2,4))
  expect_equal(sum(cube), 0)
  
  # Rows work
  stat = cbind(1:6, 1:6, 1:6)
  cube = array(generateLociCube(stat, breaks, 1:3, rows=1:4), c(3,2,4))
  expect_equal(sum(cube), 4)
  expect_equal(cube[1, 1, 1], 1)
  expect_equal(cube[2, 1, 2], 1)
  expect_equal(cube[3, 1, 2], 1)
  expect_equal(cube[3, 2, 3], 1)
  
  # 2D
  cube = array(generateLociCube(stat, breaks, 1:2, rows=1:4), c(3,2))
  expect_equal(sum(cube), 4)
  expect_equal(cube[1, 1], 1)
  expect_equal(cube[2, 1], 1)
  expect_equal(cube[3, 1], 1)
  expect_equal(cube[3, 2], 1)
})

test_that("calcPercentFpcViolation works", {
  seg_sites <- list(matrix(c(1, 1, 0, 0, 0, 
                             1, 0, 1, 0, 1, 
                             1, 1, 0, 1, 0, 
                             0, 1, 1, 0, 1, 
                             0, 0, 0, 0, 1), 5))
  
  attr(seg_sites[[1]], 'positions') <- c(0.1, 0.12, 0.5, 0.51, 0.61)
  locus_length <- matrix(c(0, 0, 100, 0, 0), 1, 5)
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length)  
  expect_equal(fpc_violations[1, ], c(mid_near=.5, mid_far=.5, outer=NaN, 
                                      between=NaN, mid=.5, perc_polym=0.05))
  
  seg_sites[[2]] <- seg_sites[[1]]
  locus_length <- rbind(locus_length, matrix(c(0, 0, 50, 0, 0), 1, 5))
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length)
  expect_equal(fpc_violations[1, ], c(mid_near=.5, mid_far=.5, outer=NaN, 
                                      between=NaN, mid=.5, perc_polym=0.05))
  expect_equal(fpc_violations[2, ], c(mid_near=.5, mid_far=.5, outer=NaN, 
                                      between=NaN, mid=.5, perc_polym=0.10))
  
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
  expect_equal(fpc_violations[4, ], c(mid_near=NaN, mid_far=NaN, outer=1, 
                                      between=1, mid=NaN, perc_polym=0.5))
  
  attr(seg_sites[[4]], 'locus') <- c(0, 0, 1, 1, 1)
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length) 
  expect_equal(fpc_violations[4, ], c(mid_near=NaN, mid_far=1, outer=NaN, 
                                      between=1, mid=1, perc_polym=1/3))
  
  attr(seg_sites[[4]], 'locus') <- c(-1, -1, 0, 0, 1)
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length) 
  expect_equal(fpc_violations[4, ], c(mid_near=NaN, mid_far=NaN, outer=1, 
                                      between=1, mid=NaN, perc_polym=1/3))
  
  attr(seg_sites[[4]], 'locus') <- c(-1, -1, 0, 1, 1)
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length) 
  expect_equal(fpc_violations[4, ], c(mid_near=NaN, mid_far=NaN, outer=1, 
                                      between=NaN, mid=NaN, perc_polym=1/6)) 
  
  attr(seg_sites[[4]], 'locus') <- c(-1, -1, -1, 1, 1)
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:5, locus_length) 
  expect_equal(fpc_violations[4, ], c(mid_near=NaN, mid_far=NaN, outer=1, 
                                      between=NaN, mid=NaN, perc_polym=0)) 
  
  # Not all individuals
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:2, locus_length) 
  expect_equal(fpc_violations[4, ], c(mid_near=NaN, mid_far=NaN, outer=NaN, 
                                      between=NaN, mid=NaN, perc_polym=0)) 
  
  fpc_violations <- calcPercentFpcViolation(seg_sites, 1:4, locus_length) 
  expect_equal(fpc_violations[4, ], c(mid_near=NaN, mid_far=NaN, outer=1, 
                                      between=NaN, mid=NaN, perc_polym=0)) 
})


test_that('Distance based classification of trios works', {
  llm <- matrix(c(0,0,1000,0,0), 5, 5, byrow = TRUE)
  expect_equal(classifyTriosByDistance(llm),
               list(both_near=numeric(), one_one=numeric(), both_far=numeric()))
  
  llm <- matrix(c(1000, 7502, 1000, 9050, 1000,
                  1000, 17502, 1000, 10050, 1000,
                  1000, 7502, 1000, 15050, 1000,
                  1000, 502, 1000, 15050, 1000,
                  1000, 6502, 1000, 35050, 1000), 5, 5, byrow = TRUE)
                  
  expect_equal(classifyTriosByDistance(llm), 
               list(both_near=1, one_one=3, both_far=2))
  expect_equal(classifyTriosByDistance(llm, near=c(5e2, 1e3), far=c(1e3, 1e5)),
               list(both_near=numeric(), one_one=4, both_far=(1:5)[-4]))
  
  # Test with only one locus
  llm <- matrix(c(0,0,1000,0,0), 1, 5, byrow = TRUE)
  expect_equal(classifyTriosByDistance(llm),
               list(both_near=numeric(), one_one=numeric(), both_far=numeric()))
})


test_that('Stat_FPC works with groups', {
  fpc = Stat_FPC$new(sumstat_tt$seg.sites, dm_tt, 1)
  expect_that(sum(fpc$get_data()), is_more_than(0))
  expect_that(sum(fpc$get_data()), 
              is_less_than(coalsimr::get_locus_number(dm_tt)+1))
  expect_equal(fpc$transform(sumstat_tt), fpc$get_data())
  
  # With groups
  fpc = Stat_FPC$new(sum_stat_grps$seg.sites.2, dm_grps, 1, group = 2)
  expect_that(sum(fpc$get_data()), is_more_than(0))
  expect_that(sum(fpc$get_data()), 
              is_less_than(coalsimr::get_locus_number(dm_grps, 2)+1))
  expect_equal(fpc$transform(sum_stat_grps), fpc$get_data())
  
  # With trios
  skip_on_cran()
  dm_trios = dm_tt + coalsimr::locus_trio(locus_length = c(10, 30, 50), 
                                          distance = c(20, 40),
                                          number = 5,
                                          group = 2)
  sumstat_trios <- simulate(dm_trios, pars=c(1,5))
  fpc = Stat_FPC$new(sumstat_trios$seg.sites.2, dm_trios, 1, group = 2) 
  expect_that(sum(fpc$get_data()), is_more_than(0))
  expect_equal(fpc$transform(sumstat_trios), fpc$get_data())
})

