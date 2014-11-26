context("seqgen simulation interface")

test_that("test.F81Model", {
  if (!test_seqgen) skip('seq-gen not installed')
  set.seed(12)
  jsfs <- dm.simSumStats(dm.f81, c(1, 10))
  expect_true(sum(jsfs$jsfs) > 0)
})

test_that("test.GtrModel", {
  if (!test_seqgen) skip('seq-gen not installed')
  set.seed(12)
  jsfs <- dm.simSumStats(dm.gtr, c(1, 10))
  expect_true(sum(jsfs$jsfs) > 0)
})

test_that("test.HkyModel", {
  if (!test_seqgen) skip('seq-gen not installed')
  set.seed(12)
  jsfs <- dm.simSumStats(dm.hky, c(1, 10))
  expect_true(sum(jsfs$jsfs) > 0)
})

test_that("test.RateHeterogenity", {
  if (!test_seqgen) skip('seq-gen not installed')
  set.seed(12)
  dm.rh <- dm.addMutationRateHeterogenity(dm.hky, 0.1, 5, categories.number = 5)
  jsfs <- dm.simSumStats(dm.rh, c(1, 10, 1))
  expect_true(sum(jsfs$jsfs) > 0)
})

test_that("test.finalizeSeqgen", {
  if (!test_seqgen) skip('seq-gen not installed')
  finalizeSeqgen = getSimProgram("seq-gen")$finalization_func
  dm.hky <- finalizeSeqgen(dm.hky)
  dm.f81 <- finalizeSeqgen(dm.f81)
  dm.gtr <- finalizeSeqgen(dm.gtr)
  expect_false(is.null(dm.hky@options[["seqgen.cmd"]]))
  expect_false(is.null(dm.hky@options[["tree.model"]]))
  expect_false(is.null(dm.f81@options[["seqgen.cmd"]]))
  expect_false(is.null(dm.f81@options[["tree.model"]]))
  expect_false(is.null(dm.gtr@options[["seqgen.cmd"]]))
  expect_false(is.null(dm.gtr@options[["tree.model"]]))
})

test_that("test.generateSeqgenOptions", {
  if (!test_seqgen) skip('seq-gen not installed')
  dm.hky@options$seqgen.cmd <- NULL
  opts <- generateSeqgenOptions(dm.hky, c(1, 10), 1)
  opts <- strsplit(opts, " ")[[1]]
  expect_true("-l" %in% opts)
  expect_true("-p" %in% opts)
  expect_true("-z" %in% opts)
  expect_true("-q" %in% opts)
  expect_true("-mHKY" %in% opts)
  expect_true("-t" %in% opts)
  expect_true("-f" %in% opts)
  expect_true("-s" %in% opts)
})

test_that("test.generateTreeModel", {
  if (!test_seqgen) skip('seq-gen not installed')
  for (dm in c(dm.hky, dm.f81, dm.gtr)) {
    dm.ms <- dm.finalize(generateTreeModel(dm))
    sum.stats <- dm.simSumStats(dm.ms, c(1, 5))
    expect_false(is.null(sum.stats$file))
    expect_true(file.exists(sum.stats$file[[1]]))
    unlink(sum.stats$file)
  }
})

test_that("test.seqgenMutationParameterNotLast", {
  if (!test_seqgen) skip('seq-gen not installed')
  dm.test <- dm.hky
  dm.test@parameters <- dm.test@parameters[3:1, ]
  cmd <- paste(generateSeqgenOptionsCmd(dm.test), collapse = "")
  expect_equal(grep("theta", cmd), 1)
  expect_equal(grep("tau", cmd), integer(0))
  expect_equal(grep("rho", cmd), integer(0))
})

test_that("test.seqgenSingleSimFunc", {
  if (!test_seqgen) skip('seq-gen not installed')
  seqgenSingleSimFunc = getSimProgram("seq-gen")$sim_func
  invisible(expect_error(seqgenSingleSimFunc(dm.tt, c(1, 10))))
  
  set.seed(100)
  sum.stats <- seqgenSingleSimFunc(dm.hky, c(1, 10))
  expect_true(is.list(sum.stats))
  expect_true(is.array(sum.stats$jsfs))
  expect_true(sum(sum.stats$jsfs) > 0)
  
  set.seed(100)
  sum.stats2 <- seqgenSingleSimFunc(dm.hky, c(1, 10))
  expect_equal(sum.stats2$jsfs, sum.stats$jsfs)
})

test_that("test.seqgenWithMsms", {
  if (!test_seqgen) skip('seq-gen not installed')
  if (!test_msms) skip('seq-gen not installed')
  dm.selsq <- dm.addPositiveSelection(dm.f81, 100, 500, population = 1, 
                                      at.time = "0.1")
  dm.selsq <- dm.finalize(dm.selsq)
  expect_false(is.null(dm.selsq@options[["seqgen.cmd"]]))
  expect_false(is.null(dm.selsq@options[["tree.model"]]))
  set.seed(4444)
  sum.stats <- dm.simSumStats(dm.selsq, c(1, 5, 250))
  expect_false(is.null(sum.stats$jsfs))
  set.seed(4444)
  sum.stats2 <- dm.simSumStats(dm.selsq, c(1, 5, 250))
  expect_false(is.null(sum.stats2$jsfs))
  expect_equal(sum.stats2, sum.stats)
})

test_that("test.simulateFpcWithSeqgen", {
  if (!test_seqgen) skip('seq-gen not installed')
  seg.sites <- dm.simSumStats(dm.addSummaryStatistic(dm.hky, 'seg.sites'),
                              c(1, 5))$seg.sites
  dm.sgfpc <- dm.addSummaryStatistic(dm.hky, 'fpc')
  dm.sgfpc <- jaatha:::calcFpcBreaks(dm.sgfpc, seg.sites, 3)
  sum.stats <- dm.simSumStats(dm.sgfpc, c(1, 5))
  expect_false(is.null(sum.stats$fpc))
  expect_equal(sum(sum.stats$fpc), 5)
})

test_that("seq-gen can simulate trios", {
  if (!test_seqgen) skip('seq-gen not installed')
  dm.lt <- dm.useLociTrios(dm.setLociLength(dm.f81, 50), c(10, 5, 20, 5, 10))
  dm.lt <- dm.addSummaryStatistic(dm.lt, 'seg.sites')
  
  sum.stats <- dm.simSumStats(dm.lt, c(1, 10))
  expect_that(sum(sum.stats$jsfs), is_less_than(sum(sapply(sum.stats$seg.sites, ncol))))
})

test_that("Generation of PMC statistic works", {
  if (!test_seqgen) skip('seq-gen not installed')
  set.seed(941)
  dm.f81 <- dm.addSummaryStatistic(dm.f81, "pmc")
  dm.f81@options[['pmc_breaks_private']] <- .5
  dm.f81@options[['pmc_breaks_fixed']] <- .5
  sum.stats <- dm.simSumStats(dm.f81, c(1, 10))
  expect_equal(length(sum.stats), 3)
  expect_false(is.null(sum.stats$pars))
  expect_false(is.null(sum.stats$pmc))
  expect_true(is.array(sum.stats[["pmc"]]))
  expect_equal(sum(sum.stats[["pmc"]]), dm.getLociNumber(dm.f81))
})
