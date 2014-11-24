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
