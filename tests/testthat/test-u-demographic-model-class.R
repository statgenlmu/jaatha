context("Demographic Model")

test_that("test.addFeature", {
  dm <- dm.createDemographicModel(11:12, 100)
  n.feat <- nrow(dm@features)
  dm <- addFeature(dm, type = "split", parameter = "tau", lower.range = 1, 
                   upper.range = 10)
  expect_equal(nrow(dm@features), n.feat + 1)
  dm <- addFeature(dm, "mutation", "theta", fixed.value = 5, 
                   pop.source = 1, pop.sink = 2, time.point = "t2", group = 3)
  expect_equal(nrow(dm@features), n.feat + 2)
})

test_that("test.addParameter", {
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- dm.addParameter(dm, "theta", 1, 5)
  dm <- dm.addParameter(dm, "rho", fixed = 20)
  expect_equal(c("theta", "rho"), dm@parameters$name)
  expect_equal(c(F, T), dm@parameters$fixed)
  expect_equal(c(1, 20), dm@parameters$lower.range)
  expect_equal(c(5, 20), dm@parameters$upper.range)
})

test_that("test.addPositiveSelection", {
  dm <- dm.addPositiveSelection(dm.tt, 1, 2, population = 1, 
                                at.time = "2")
  expect_true("pos.selection" %in% dm@features$type)
  dm <- dm.addPositiveSelection(dm.tt, fixed.strength = 1, 
                                population = 1, at.time = "2")
  expect_true("pos.selection" %in% dm@features$type)
  expect_error(dm.addPositiveSelection(dm.tt, 1, 2, at.time = "2"))
  expect_error(dm.addPositiveSelection(dm.tt, 1, 2, population = 1))
})

test_that("test.addSampleSize", {
  dm <- dm.tt
  expect_equal(sum(dm@features$type == "sample"), 2)
  expect_equal(subset(dm@features, type == "sample" & pop.source == 
                        1)$parameter, "11")
  expect_equal(subset(dm@features, type == "sample" & pop.source == 
                        2)$parameter, "12")
  expect_equal(nrow(subset(dm@features, type == "sample" & 
                             group == 0)), 2)
  expect_equal(nrow(subset(dm@features, type == "sample" & 
                             time.point == "0")), 2)
  dm <- dm.addSampleSize(dm, c(2, 5), group = 1)
  expect_equal(sum(dm@features$type == "sample"), 4)
  expect_equal(nrow(subset(dm@features, type == "sample" & 
                             group == 0)), 2)
  expect_equal(nrow(subset(dm@features, type == "sample" & 
                             group == 1)), 2)
})

test_that("test.dm.addSummaryStatistics", {
  dm <- dm.createDemographicModel(11:12, 100)
  expect_equal(length(dm.getSummaryStatistics(dm)), 1)
  expect_true(dm.getSummaryStatistics(dm) == "jsfs")
  dm <- dm.addSummaryStatistic(dm, "seg.sites")
  expect_equal(length(dm.getSummaryStatistics(dm)), 2)
  expect_true(all(dm.getSummaryStatistics(dm) == c("jsfs", 
                                                   "seg.sites")))
  dm <- dm.addSummaryStatistic(dm, "seg.sites")
  expect_equal(length(dm.getSummaryStatistics(dm)), 2)
  expect_true(all(dm.getSummaryStatistics(dm) == c("jsfs", 
                                                   "seg.sites")))
  dm <- dm.addSummaryStatistic(dm, "file", group = 2)
  expect_equal(length(dm.getSummaryStatistics(dm, group = 1)), 
               2)
  expect_equal(length(dm.getSummaryStatistics(dm, group = 2)), 
               3)
  dm <- dm.addSummaryStatistic(dm, "fpc", group = 1)
  expect_equal(length(dm.getSummaryStatistics(dm, group = 1)), 
               3)
  expect_equal(length(dm.getSummaryStatistics(dm, group = 2)), 
               3)
  expect_error(dm.addSummaryStatistic(dm, "no.existing.sumstat"))
  expect_error(dm.addSummaryStatistic(dm, 1:10))
})

test_that("test.dm.finalize", {
  dm <- dm.finalize(dm.grp)
  dm.1 <- dm@options$grp.models[["1"]]
  dm.2 <- dm@options$grp.models[["2"]]
  dm.3 <- dm@options$grp.models[["3"]]
  expect_false(is.null(dm.1))
  expect_false(is.null(dm.2))
  expect_false(is.null(dm.3))
  expect_false(is.null(dm.1@options$ms.cmd))
  expect_false(is.null(dm.2@options$ms.cmd))
  expect_false(is.null(dm.3@options$ms.cmd))
})

test_that("test.dm.getSummaryStatistic", {
  expect_equal(dm.getSummaryStatistics(dm.tt),  as.factor("jsfs"))
  expect_equal(length(dm.getSummaryStatistics(dm.tt)), 1)
  expect_true(dm.getSummaryStatistics(dm.grp, 1) == "jsfs")
  expect_equal(length(dm.getSummaryStatistics(dm.grp, 1)),  1)
  expect_true(dm.getSummaryStatistics(dm.grp, 2) == "jsfs")
  expect_equal(length(dm.getSummaryStatistics(dm.grp, 2)), 1)
  expect_true(all(c("jsfs", "fpc") %in% dm.getSummaryStatistics(dm.fpc)))
  expect_equal(length(dm.getSummaryStatistics(dm.fpc)), 2)
})

test_that("test.generateGroupModel", {
  dm <- generateGroupModel(dm.tt, 1)
  expect_equal(dm@features, dm.tt@features)
  expect_equal(dm@sum.stats, dm.tt@sum.stats)
  expect_equal(dm@options, dm.tt@options)
  dm <- dm.setLociLength(dm.tt, 23, group = 1)
  dm <- generateGroupModel(dm, 1)
  expect_equal(nrow(dm@features), nrow(dm.tt@features))
  expect_true(all(dm@features$group == 0))
  expect_equal(dm.getLociLength(dm), 23)
  dm.3 <- dm.setLociLength(dm.tt, 23, group = 1)
  dm.3 <- dm.setLociLength(dm.3, 30, group = 2)
  dm.3 <- dm.setLociNumber(dm.3, 31, group = 2)
  dm <- generateGroupModel(dm.3, 1)
  expect_equal(nrow(dm@features), nrow(dm.tt@features))
  expect_true(all(dm@features$group == 0))
  expect_equal(dm.getLociLength(dm), 23)
  dm <- generateGroupModel(dm.3, 2)
  expect_equal(nrow(dm@features), nrow(dm.tt@features))
  expect_true(all(dm@features$group == 0))
  expect_equal(dm.getLociLength(dm), 30)
  expect_equal(dm.getLociNumber(dm), 31)
  sum.stats <- dm.tt@sum.stats
  dm <- dm.addSummaryStatistic(dm.tt, "file", group = 2)
  dm.1 <- generateGroupModel(dm, 1)
  expect_true(sum.stats$name %in% dm.1@sum.stats$name)
  expect_true(dm.1@sum.stats$name %in% sum.stats$name)
  dm.2 <- generateGroupModel(dm, 2)
  expect_equal(nrow(dm.2@sum.stats), nrow(sum.stats) + 1)
  expect_true("file" %in% dm.getSummaryStatistics(dm.2))
  dm@options[["bli"]] <- 0
  dm@options[["bla"]] <- 0
  dm@options[["blub"]] <- 0
  dm@options[["group.1"]] <- list(blub = 1, bli = 1)
  dm@options[["group.2"]] <- list(blub = 2, bla = 2)
  dm.1 <- generateGroupModel(dm, 1)
  dm.2 <- generateGroupModel(dm, 2)
  expect_equal(dm.1@options$bli, 1)
  expect_equal(dm.1@options$bla, 0)
  expect_equal(dm.1@options$blub, 1)
  expect_equal(dm.2@options$bli, 0)
  expect_equal(dm.2@options$bla, 2)
  expect_equal(dm.2@options$blub, 2)
})

test_that("test.getGroups", {
  expect_equal(dm.getGroups(dm.tt), 1)
  dm <- dm.setLociLength(dm.tt, 23, group = 1)
  expect_equal(dm.getGroups(dm), 1)
  dm <- dm.setLociLength(dm.tt, 32, group = 2)
  expect_equal(dm.getGroups(dm), 1:2)
  dm <- dm.setLociNumber(dm, 32, group = 4)
  expect_equal(dm.getGroups(dm), c(1:2, 4))
})

test_that("test.getSampleSize", {
  dm.test <- dm.tt
  expect_equal(dm.getSampleSize(dm.test), c(11L, 12L))
  dm.test <- dm.addSampleSize(dm.test, c(2, 5), group = 2)
  expect_equal(dm.getSampleSize(dm.test, 0), c(11L, 12L))
  expect_equal(dm.getSampleSize(dm.test, 1), c(11L, 12L))
  expect_equal(dm.getSampleSize(dm.test, 2), c(2L, 5L))
  expect_error(dm.getSampleSize(dm.test))
  dm.test <- dm.setLociLength(dm.test, 15, group = 3)
  expect_equal(dm.getSampleSize(dm.test, 3), c(11L, 12L))
})

test_that("test.getThetaName", {
  expect_equal(getThetaName(dm.tt), "theta")
  expect_equal(getThetaName(dm.fpc), "theta")
  if (test_seqgen) {
    expect_equal(getThetaName(dm.hky), "theta")
    expect_equal(getThetaName(dm.f81), "theta")
  }
  dm.test <- dm.createDemographicModel(11:12, 100)
  dm.test <- dm.addMutation(dm.test, 1, 5, new.par.name = "abcd")
  dm.test <- dm.addRecombination(dm.test, 1, 5)
  expect_equal(getThetaName(dm.test), "abcd")
})

test_that("test.parInRange", {
  checkParInRange(dm.tt, c(1, 5))
  checkParInRange(dm.tt, c(2, 7))
  checkParInRange(dm.tt, c(2.1, 7.7))
  expect_error(checkParInRange(dm.tt, c(0, 5)))
  expect_error(checkParInRange(dm.tt, c(0, -1)))
  expect_error(checkParInRange(dm.tt, c(10, 1)))
  expect_error(checkParInRange(dm.tt, c(100, 100)))
  expect_error(checkParInRange(dm.tt, 1))
  expect_error(checkParInRange(dm.tt, matrix(1, 2, 2)))
  expect_error(checkParInRange(dm.tt, NULL))
})

test_that("test.scaleDemographicModel", {
  dm <- dm.setLociNumber(dm.tt, 25, group = 1)
  dm <- dm.setLociNumber(dm, 27, group = 2)
  dm <- scaleDemographicModel(dm, 5)
  expect_equal(dm.getLociNumber(dm, 0), 2L)
  expect_equal(dm.getLociNumber(dm, 1), 5L)
  expect_equal(dm.getLociNumber(dm, 2), 5L)
})

test_that("searchFeature", {
  expect_equal(nrow(searchFeature(dm.tt, type = "sample")), 2)
  expect_equal(nrow(searchFeature(dm.tt, type = c("split", "sample"))), 3)
  expect_equal(nrow(searchFeature(dm.tt, type = c("split", "sample"), 
                                  pop.sink = NA)), 2)
  expect_equal(nrow(searchFeature(dm.tt, time.point = "tau")), 1)
  expect_equal(nrow(searchFeature(dm.tt, time.point = NA)), 3)
  
  # With groups
  expect_equal(searchFeature(dm.grp, 'loci.length', group=0)$parameter, '1000')
  expect_equal(searchFeature(dm.grp, 'loci.length', group=1)$parameter, '100')
  expect_equal(searchFeature(dm.grp, 'loci.length', group=2)$parameter, '200')
  expect_equal(searchFeature(dm.grp, 'loci.length', group=3)$parameter, '1000')
  
  expect_equal(nrow(searchFeature(dm.grp, group = 0)), 7)
  expect_equal(nrow(searchFeature(dm.grp, group = 1)), 7)
  expect_equal(nrow(searchFeature(dm.grp, group = 2)), 7)
  expect_equal(nrow(searchFeature(dm.grp, group = 3)), 7)
})

test_that("set and get loci length works", {
  dm_tmp <- dm.setLociLength(dm.tt, 17)
  expect_equal(sum(dm_tmp@features$type == "loci.length"), 1)
  expect_equal(searchFeature(dm_tmp, "loci.length")$parameter, "17")
  expect_equal(dm.getLociLength(dm_tmp), 17)
  
  dm_tmp <- dm.setLociLength(dm_tmp, 23, group = 2)
  expect_equal(sum(dm_tmp@features$type == "loci.length"), 2)
  expect_equal(dm.getLociLength(dm_tmp), 17)
  expect_equal(dm.getLociLength(dm_tmp, 1), 17)
  expect_equal(dm.getLociLength(dm_tmp, 2), 23)
  
  dm_tmp <- dm.setLociLength(dm_tmp, 32, group = 1)
  expect_equal(dm.getLociLength(dm_tmp), 32)
  expect_equal(dm.getLociLength(dm_tmp, 1), 32)
  expect_equal(dm.getLociLength(dm_tmp, 2), 23)
})

test_that("set and get loci number works", {
  dm <- dm.setLociNumber(dm.tt, 17)
  expect_equal(sum(dm@features$type == "loci.number"), 1)
  expect_equal(dm@features$parameter[dm@features$type == "loci.number"], "17")
  expect_equal(dm.getLociNumber(dm), 17)
  
  dm <- dm.setLociNumber(dm, 23, group = 2)
  expect_equal(sum(dm@features$type == "loci.number"), 2)
  expect_equal(dm.getLociNumber(dm), 17)
  expect_equal(dm.getLociNumber(dm, 1), 17)
  expect_equal(dm.getLociNumber(dm, 2), 23)
  
  dm <- dm.setLociNumber(dm, 32, group = 1)
  expect_equal(sum(dm@features$type == "loci.number"), 3)
  expect_equal(dm.getLociNumber(dm), 32)
  expect_equal(dm.getLociNumber(dm, 1), 32)
  expect_equal(dm.getLociNumber(dm, 2), 23)
})

test_that("test.setMutationModel", {
  dm <- dm.createThetaTauModel(11:12, 100)
  dm <- dm.setMutationModel(dm, "HKY")
  expect_true("mutation.model" %in% dm@parameters$name)
  expect_error(dm <- dm.setMutationModel(dm, "bla"))
})

test_that("test.simSumStats", {
  expect_error(dm.simSumStats(1))
  expect_error(dm.simSumStats(dm.tt, 1))
  expect_error(dm.simSumStats(dm.tt, 1:3))
  expect_error(dm.simSumStats(dm.tt, c(2, 50)))
  sum.stats <- dm.simSumStats(dm.tt, c(1, 5), "jsfs")
  expect_true(is.list(sum.stats))
  expect_false(is.null(sum.stats$jsfs))
  expect_true(sum(sum.stats$jsfs) > 0)
  sum.stats <- dm.simSumStats(dm.grp, c(1, 5), "jsfs")
  expect_false(is.null(sum.stats))
  expect_false(is.null(sum.stats$jsfs.1))
  expect_true(sum(sum.stats$jsfs.1) > 0)
  expect_false(is.null(sum.stats$jsfs.2))
  expect_true(sum(sum.stats$jsfs.2) > 0)
  expect_false(is.null(sum.stats$jsfs.3))
  expect_true(sum(sum.stats$jsfs.3) > 0)
})

test_that("Loci trios are added to model", {
  if (!test_seqgen) return()
  dm.lt <- dm.useLociTrios(dm.hky, c(3, 2, 5, 1, 4))
  expect_true("trio.1" %in% dm.lt@features$type)
  expect_true("trio.2" %in% dm.lt@features$type)
  expect_true("trio.3" %in% dm.lt@features$type)
  expect_true("trio.4" %in% dm.lt@features$type)
  expect_true("trio.5" %in% dm.lt@features$type)
  expect_error(dm.useLociTrios(dm.hky, c(5, 5, 5, 1, 4)))
  expect_error(dm.useLociTrios(dm.hky, c(5, 5, 5)))
})

test_that("getTrioOptions works", {
  if (!test_seqgen) return()
  dm.lt <- dm.useLociTrios(dm.setLociLength(dm.f81, 50), c(10, 5, 20, 5, 10))
  
  expect_equal(dm.getLociTrioOptions(dm.lt), c(10, 5, 20, 5, 10))
  expect_true(is.na(dm.getLociTrioOptions(dm.tt)))
  
  dm.lt <- dm.useLociTrios(dm.setLociLength(dm.f81, 50), c(10, 5, 20, 5, 10), 
                           group=2)
  expect_true(is.na(dm.getLociTrioOptions(dm.lt)))
  expect_equal(dm.getLociTrioOptions(dm.lt, group=2), c(10, 5, 20, 5, 10))
})