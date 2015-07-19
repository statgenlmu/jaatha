context("SumStat FPC")

test_that("fpc initialization works", {
  stat <- coala::sumstat_four_gamete()
  fpc = Stat_FPC$new(sumstat_tt$seg_sites, dm_tt, stat)
  
  expect_that(fpc$get_data(), is_a("integer"))
  expect_that(sum(fpc$get_data()), is_more_than(0))
  expect_that(fpc$get_breaks(), is_a("list"))
  expect_equal(length(fpc$get_breaks()), 6)
})


test_that("fpc break calculation works", {
  stat <- coala::sumstat_four_gamete()
  fpc = Stat_FPC$new(sumstat_tt$seg_sites, dm_tt, stat)
  expect_false(is.null(fpc$get_breaks()))
  expect_false(is.null(fpc$get_breaks()$mid_near))
  expect_false(is.null(fpc$get_breaks()$mid_far))
  expect_false(is.null(fpc$get_breaks()$perc_polym))
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


test_that("Distance based classification of trios works", {
  skip("Awaiting removal")
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


test_that("Stat_FPC works with groups", {
  skip("Awaiting removal")
  fpc = Stat_FPC$new(sumstat_tt$seg_sites, dm_tt, 1)
  expect_that(sum(fpc$get_data()), is_more_than(0))
  expect_that(sum(fpc$get_data()), 
              is_less_than(coala::get_locus_number(dm_tt)+1))
  expect_equal(fpc$transform(sumstat_tt), fpc$get_data())
  
  # With groups
  fpc = Stat_FPC$new(sum_stat_grps$seg_sites.2, dm_grps, 1, group = 2)
  expect_that(sum(fpc$get_data()), is_more_than(0))
  expect_that(sum(fpc$get_data()), 
              is_less_than(coala::get_locus_number(dm_grps, 2)+1))
  expect_equal(fpc$transform(sum_stat_grps), fpc$get_data())
  
  # With trios
  skip_on_cran()
  dm_trios = dm_tt + coala::locus_trio(locus_length = c(10, 30, 50), 
                                          distance = c(20, 40),
                                          number = 5,
                                          group = 2)
  sumstat_trios <- simulate(dm_trios, pars=c(1,5))
  fpc = Stat_FPC$new(sumstat_trios$seg_sites.2, dm_trios, 1, group = 2) 
  expect_that(sum(fpc$get_data()), is_more_than(0))
  expect_equal(fpc$transform(sumstat_trios), fpc$get_data())
})

