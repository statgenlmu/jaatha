context("msms simulation interface")

test_that("test.callMsms", {
  if (!test_msms) return()
  jar.path = getJaathaVariable("msms.jar")
  ms.args <- "5 1 -r 10 100 -t 5 -I 2 3 2 1"
  msms.args <- ""
  set.seed(17)
  out.file <- callMsms(jar.path, ms.args, msms.args)
  set.seed(17)
  out.file.2 <- callMsms(jar.path, ms.args, msms.args)
  set.seed(20)
  out.file.3 <- callMsms(jar.path, ms.args, msms.args)
  expect_equal(file.info(out.file.2)$size, file.info(out.file)$size)
  expect_true(file.info(out.file)$size != file.info(out.file.3)$size)
  unlink(c(out.file, out.file.2, out.file.3))
})

test_that("test.generateMsmsOptionsCommand", {
  dm <- dm.addPositiveSelection(dm.tt, 100, 500, population = 1, 
                                at.time = "2")
  opts <- generateMsmsOptionsCommand(dm)
  s <- 5
  expect_true("-SI" %in% eval(parse(text = opts)))
  expect_true("-SAA" %in% eval(parse(text = opts)))
  expect_true("-SAa" %in% eval(parse(text = opts)))
})

test_that("test.msmsPrint", {
  if (!test_msms) return()
  tmp_file <- tempfile()
  sink(tmp_file)
  print(dm.sel)
  sink(NULL)
  unlink(tmp_file)
})

test_that("test.msmsSimFunc", {
  if (!test_msms) return()
  sum.stats <- msmsSimFunc(dm.sel, c(1, 1.5, 1500, 5))
  expect_true(is.list(sum.stats))
  expect_false(is.null(sum.stats$pars))
  expect_true(is.matrix(sum.stats$jsfs))
  expect_equal(length(sum.stats), 2)
  dm.sel@sum.stats <- data.frame()
  dm.sel <- dm.addSummaryStatistic(dm.sel, "seg.sites")
  dm.sel <- dm.addSummaryStatistic(dm.sel, "file")
  sum.stats <- msmsSimFunc(dm.sel, c(1, 1.5, 1500, 5))
  expect_true(is.list(sum.stats))
  expect_false(is.null(sum.stats$pars))
  expect_true(is.list(sum.stats$seg.sites))
  expect_equal(length(sum.stats$seg.sites), 3)
  expect_true(file.exists(sum.stats$file))
  expect_equal(length(sum.stats), 3)
  unlink(sum.stats$file)
  dm.sel@sum.stats <- data.frame()
  dm.sel <- dm.addSummaryStatistic(dm.sel, "fpc")
  dm.sel <- calcFpcBreaks(dm.sel, sum.stats$seg.sites)
  sum.stats <- msmsSimFunc(dm.sel, c(1, 1.5, 1500, 5))
  expect_true(is.list(sum.stats))
  expect_false(is.null(sum.stats$pars))
  expect_true(is.matrix(sum.stats$fpc))
  expect_equal(sum(sum.stats$fpc), 3)
  expect_equal(length(sum.stats), 2)
})

test_that("Generation of PMC statistic works", {
  if (!test_msms) return()
  set.seed(941)
  dm.sel <- dm.addSummaryStatistic(dm.sel, "pmc")
  dm.sel@options[['pmc_breaks_private']] <- .5
  dm.sel@options[['pmc_breaks_fixed']] <- .5
  sum.stats <- msmsSimFunc(dm.sel, c(0.1, 2, 2, 500))
  expect_equal(length(sum.stats), 3)
  expect_false(is.null(sum.stats$pars))
  expect_false(is.null(sum.stats$pmc))
  expect_true(is.array(sum.stats[["pmc"]]))
  expect_equal(sum(sum.stats[["pmc"]]), dm.getLociNumber(dm.sel))
})

test_that("msms can simulate subgroups", {
  if (!test_msms) return()
  dm_tmp <- dm.addSubgroups(dm.sel, 3)
  dm_tmp <- dm.addSummaryStatistic(dm_tmp, 'seg.sites')
  sum_stats <- dm.simSumStats(dm_tmp, c(1, 5, 2, 1500))
  expect_equal(length(sum_stats$seg.sites), 3)
  for (seg_sites in sum_stats$seg.sites) {
    expect_true(is.matrix(seg_sites))
  }
  expect_false(any(is.na(sum_stats$jsfs)))
  expect_true(sum(sum_stats$jsfs) > 0)
})

test_that("Gamma Distributed Rates work", {
  if (!test_msms) return()
  dm_tmp <- dm.addPositiveSelection(dm.mig, 1000, 2000, population=2,
                                    at.time='tau/2', variance='s^2',
                                    var.classes = 3)
  dm_tmp <- dm.addSummaryStatistic(dm_tmp, 'seg.sites')
  
  sum_stats <- dm.simSumStats(dm_tmp, c(1, 5, 2, 1500))
  for (seg_sites in sum_stats$seg.sites) {
    expect_true(is.matrix(seg_sites))
  }
  expect_false(any(is.na(sum_stats$jsfs)))
  expect_true(sum(sum_stats$jsfs) > 0)
  
  dm_tmp <- dm.addPositiveSelection(dm.mig, 1000, 2000, population=2,
                                    at.time='tau/2', variance='s^2',
                                    var.classes = 6)
  dm_tmp <- dm.addSummaryStatistic(dm_tmp, 'seg.sites')
  
  sum_stats <- dm.simSumStats(dm_tmp, c(1, 5, 2, 1500))
  for (seg_sites in sum_stats$seg.sites) {
    expect_true(is.matrix(seg_sites))
  }
  expect_false(any(is.na(sum_stats$jsfs)))
  expect_true(sum(sum_stats$jsfs) > 0)
})

