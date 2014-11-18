context("ms simulation interface")

test_that("msSimFunc is working", {
  set.seed(789)
  sum_stats <- msSingleSimFunc(dm.tt, c(1, 10))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
  
  set.seed(789)
  sum_stats2 <- msSingleSimFunc(dm.tt, c(1, 10))
  expect_equal(sum_stats, sum_stats2)
})

test_that("msSimFunc works with inter-locus variation", {
  dm_tmp <- dm.addInterLocusVariation(dm.tt)
  set.seed(117)
  sum_stats <- msSingleSimFunc(dm_tmp, c(1, 5))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
  
  set.seed(117)
  sum_stats2 <- msSingleSimFunc(dm_tmp, c(1, 5))
  expect_equal(sum_stats$jsfs, sum_stats2$jsfs)
})

test_that("generation of PMC statistic works for ms", {
  set.seed(941)
  dm.tt <- dm.addSummaryStatistic(dm.tt, "pmc")
  dm.tt@options[['pmc_breaks_private']] <- .5
  dm.tt@options[['pmc_breaks_fixed']] <- .5
  sum.stats <- msSingleSimFunc(dm.tt, c(1, 0.1))
  expect_equal(length(sum.stats), 3)
  expect_false(is.null(sum.stats$pars))
  expect_false(is.null(sum.stats$pmc))
  expect_true(is.array(sum.stats[["pmc"]]))
  expect_equal(sum(sum.stats[["pmc"]]), dm.getLociNumber(dm.tt))
})

test_that("the ms sim program exists", {
  expect_false(is.null(.jaatha$sim_progs[["ms"]]))
})

test_that("ms can simulate subgroups", {
  dm_tmp <- dm.addSubgroups(dm.tt, 3)
  dm_tmp <- dm.addSummaryStatistic(dm_tmp, 'seg.sites')
  sum_stats <- dm.simSumStats(dm_tmp, c(1, 5))
  expect_equal(length(sum_stats$seg.sites), 10)
  for (seg_sites in sum_stats$seg.sites) {
    expect_true(is.matrix(seg_sites))
  }
  expect_false(any(is.na(sum_stats$jsfs)))
  expect_true(sum(sum_stats$jsfs) > 0)
  
  dm_tmp <- dm.addSubgroups(dm.tt, 20)
  dm_tmp <- dm.addSummaryStatistic(dm_tmp, 'seg.sites')
  sum_stats <- dm.simSumStats(dm_tmp, c(1, 5))
  expect_equal(length(sum_stats$seg.sites), 10)
})

test_that("ms can simulate variable rates", {
  dm_tmp <- dm.createDemographicModel(5:6, 10, 1000)
  dm_tmp <- dm.addSpeciationEvent(dm_tmp, 0.01, 5, new.time.point.name="tau")  
  dm_tmp <- dm.addRecombination(dm_tmp, fixed.value=2)
  dm_tmp <- dm.addMutation(dm_tmp, 1, 20, variance = 20)
  dm_tmp <- dm.addSummaryStatistic(dm_tmp, 'seg.sites')
  sum_stats <- dm.simSumStats(dm_tmp, c(1, 5))
  expect_equal(length(sum_stats$seg.sites), 10)
  for (seg_sites in sum_stats$seg.sites) {
    expect_true(is.matrix(seg_sites))
  }
  expect_false(any(is.na(sum_stats$jsfs)))
  expect_true(sum(sum_stats$jsfs) > 0)
})

test_that("sample subgroup sizes works", {
  expect_equal(sampleSubgroupSizes(dm.tt), 10)
  expect_equal(sampleSubgroupSizes(dm.mig), 10)
  expect_equal(sampleSubgroupSizes(dm.hky), 5)
  
  dm_tmp <- dm.addSubgroups(dm.tt, 5)
  for (i in 1:10) {
    subg_sizes <- sampleSubgroupSizes(dm_tmp, c(1, 5))
    expect_equal(length(subg_sizes), 5)
    expect_equal(sum(subg_sizes), 10)
  }
  
  # Zero inflation
  dm_tmp <- dm.addSubgroups(dm.tt, 2, zero_inflation = 0.5)
  subg_sizes <- sampleSubgroupSizes(dm_tmp, c(1, 5))
  expect_equal(length(subg_sizes), 3)
  expect_equal(sum(subg_sizes), 10)
  expect_equal(subg_sizes[1], 5)
  
  dm_tmp <- dm.addParameter(dm.tt, 0.01, 0.99, par.name = 'zi')
  dm_tmp <- dm.addSubgroups(dm_tmp, 2, zero_inflation = 'zi')
  expect_equal(sampleSubgroupSizes(dm_tmp, c(1, 5, 0.11))[1], 1)
  expect_equal(sampleSubgroupSizes(dm_tmp, c(1, 5, 0.13))[1], 1)
  expect_equal(sampleSubgroupSizes(dm_tmp, c(1, 5, 0.17))[1], 2)
  expect_equal(sampleSubgroupSizes(dm_tmp, c(1, 5, 0.251))[1], 3)
  expect_equal(sampleSubgroupSizes(dm_tmp, c(1, 5, 0.71))[1], 7)
  expect_equal(sampleSubgroupSizes(dm_tmp, c(1, 5, 0.97))[1], 10)
})