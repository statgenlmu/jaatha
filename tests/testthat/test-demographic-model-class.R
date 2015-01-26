context("Demographic Model")

test_that("creating models works", {
  dm <- dm.createDemographicModel(11:12, 111, 1234) 
  expect_equal(dm.getSampleSize(dm), 11:12)
  expect_equal(dm.getLociNumber(dm), 111)
  expect_equal(dm.getLociLength(dm), 1234)
})


test_that("adding parameters works", {
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- dm.addParameter(dm, "theta", 1, 5)
  expect_equal("theta", dm@parameters$name)
  expect_equal(1, dm@parameters$lower.range)
  expect_equal(5, dm@parameters$upper.range)
})


test_that("adding features works", {
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- dm + Feature$new('blub', 5)
  expect_equal(nrow(searchFeature(dm, 'blub')), 1)
  
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- dm + Feature$new('bli', par_range('bla', 1, 5))
  expect_equal(nrow(searchFeature(dm, 'bli')), 1)
  expect_true('bla' %in% dm@parameters$name)
  
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- dm + Feature$new('bli', par_range('bla', 1, 5), variance='15')
  expect_true(hasInterLocusVariation(dm))
  
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- dm + Feature$new('bli', par_range('bla', 1, 5), group=1, variance='15')
  expect_false(hasInterLocusVariation(dm, 0))
  expect_true(hasInterLocusVariation(dm, 1))
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
  expect_equal(length(dm.getSummaryStatistics(dm, group = 1)), 2)
  expect_equal(length(dm.getSummaryStatistics(dm, group = 2)), 3)
  dm <- dm.addSummaryStatistic(dm, "trees", group = 1)
  expect_equal(length(dm.getSummaryStatistics(dm, group = 1)), 3)
  expect_equal(length(dm.getSummaryStatistics(dm, group = 2)), 3)
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
})

test_that("test.generateGroupModel", {
  dm <- generateGroupModel(dm.tt, 1)
  expect_equal(dm@features, dm.tt@features)
  expect_equal(dm@sum.stats, dm.tt@sum.stats)
  expect_equal(dm@options, dm.tt@options)
  
  dm <- dm.addLocus(dm.tt, 23, 10, group = 1)
  dm <- generateGroupModel(dm, 1)
  expect_equal(nrow(dm@features), nrow(dm.tt@features))
  expect_true(all(dm@features$group == 0))
  expect_equal(dm.getLociLength(dm), 23)
  
  dm.3 <- dm.addLocus(dm.tt, 23, 5, group = 1)
  dm.3 <- dm.addLocus(dm.3, 30, 31, group = 2)
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
  dm <- dm.addLocus(dm.tt, 23, 10, group = 1)
  expect_equal(dm.getGroups(dm), 1)
  dm <- dm.addLocus(dm.tt, 32, 10, group = 2)
  expect_equal(dm.getGroups(dm), 1:2)
  dm <- dm.addLocus(dm, 32, group = 4)
  expect_equal(dm.getGroups(dm), c(1:2, 4))
})

test_that("test.getSampleSize", {
  dm.test <- dm.tt
  expect_equal(dm.getSampleSize(dm.test), c(11L, 12L))
})

test_that("test.getThetaName", {
  expect_equal(getThetaName(dm.tt), "theta")
  if (test_seqgen) {
    expect_equal(getThetaName(dm.hky), "theta")
    expect_equal(getThetaName(dm.f81), "theta")
  }
  dm.test <- dm.createDemographicModel(11:12, 100)
  dm.test <- dm.addMutation(dm.test, 1, 5, parameter = "abcd")
  dm.test <- dm.addRecombination(dm.test, 1, 5)
  expect_equal(getThetaName(dm.test), "abcd")
})

test_that("test.parInRange", {
  checkParInRange(dm.tt, c(1, 5))
  checkParInRange(dm.tt, c(2, 7))
  checkParInRange(dm.tt, c(0.5, 7.7))
  expect_error(checkParInRange(dm.tt, c(0, 5)))
  expect_error(checkParInRange(dm.tt, c(0, -1)))
  expect_error(checkParInRange(dm.tt, c(10, 1)))
  expect_error(checkParInRange(dm.tt, c(100, 100)))
  expect_error(checkParInRange(dm.tt, 1))
  expect_error(checkParInRange(dm.tt, matrix(1, 2, 2)))
  expect_error(checkParInRange(dm.tt, NULL))
})

test_that("test.scaleDemographicModel", {
  dm <- dm.createDemographicModel(11:12, 10)
  dm <- dm.addLocus(dm, 10, 25, group = 1)
  dm <- dm.addLocus(dm, 15, 25, group = 2)
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
  expect_equal(nrow(searchFeature(dm.tt, time.point = NA)), 1)
})

test_that("set and get loci length works", {
  dm_tmp <- dm.setLociLength(dm.tt, 17)
  expect_equal(dm.getLociLength(dm_tmp), 17)
  
  dm_tmp <- dm.addLocus(dm_tmp, 22, 10, group = 2)
  dm_tmp <- dm.setLociLength(dm_tmp, 23, group = 2)
  expect_equal(dm.getLociLength(dm_tmp), 17)
  expect_equal(dm.getLociLength(dm_tmp, 1), 17)
  expect_equal(dm.getLociLength(dm_tmp, 2), 23)

  dm_tmp <- dm.addLocus(dm_tmp, 32, 10, group = 1)
  dm_tmp <- dm.setLociLength(dm_tmp, 32, group = 1)
  expect_equal(dm.getLociLength(dm_tmp), 32)
  expect_equal(dm.getLociLength(dm_tmp, 1), 32)
  expect_equal(dm.getLociLength(dm_tmp, 2), 23)
})

test_that("set and get loci number works", {
  dm <- dm.setLociNumber(dm.tt, 17)
  expect_equal(dm.getLociNumber(dm), 17)
  
  dm <- dm.addLocus(dm, 22, 10, group = 2)
  dm <- dm.setLociNumber(dm, 23, group = 2)
  expect_equal(dm.getLociNumber(dm), 17)
  expect_equal(dm.getLociNumber(dm, 1), 17)
  expect_equal(dm.getLociNumber(dm, 2), 23)

  dm <- dm.addLocus(dm, 32, 10, group = 1)
  dm <- dm.setLociNumber(dm, 32, group = 1)
  expect_equal(dm.getLociNumber(dm), 32)
  expect_equal(dm.getLociNumber(dm, 1), 32)
  expect_equal(dm.getLociNumber(dm, 2), 23)
})

test_that('locus length matrix generations works', {
  # Multiple loci with equal length
  dimnames <- list(NULL,  c('length_l', 'length_il', 'length_m', 
                            'length_ir', 'length_r') )
  
  expect_equal(dm.getLociLengthMatrix(dm.tt), 
               matrix(c(0, 0, 100, 0, 0), 5, 5, TRUE, dimnames))
  
  # Multiple loci with differnt length 
  dm <- dm.addLocus(dm.tt, 21, 1, group = 2)
  dm <- dm.addLocus(dm, 22, 1, group = 2)
  dm <- dm.addLocus(dm, 23, 1, group = 2)
  expect_equal(dm.getLociLengthMatrix(dm, group = 2),
               matrix(c(0, 0, 21, 0, 0,
                        0, 0, 22, 0, 0,
                        0, 0, 23, 0, 0), 3, 5, TRUE, dimnames))
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
  
  sum.stats <- sum.stats.grp
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
  dm.lt <- dm.addLocusTrio(dm.hky, locus_length =  c(3, 5, 7),
                           distance = c(4, 6), group = 2)
  
  expect_true(all(dm.getLociLengthMatrix(dm.lt, 2) == 3:7))
  expect_true(nrow(searchFeature(dm.lt, 'locus_trios', group = 2)) > 0)
})

test_that("Adding and Getting inter locus variation works", {
  expect_false(dm.hasInterLocusVariation(dm.tt))
  
  dm_tmp <- dm.addInterLocusVariation(dm.tt, 0)
  expect_true(dm.hasInterLocusVariation(dm_tmp))
  dm_tmp <- dm.addInterLocusVariation(dm_tmp, 0)
  expect_true(dm.hasInterLocusVariation(dm_tmp))
  expect_equal(nrow(searchFeature(dm_tmp, 'inter_locus_variation')), 1)
  
  dm_tmp <- dm.addInterLocusVariation(dm.tt, 2)
  expect_false(dm.hasInterLocusVariation(dm_tmp))  
  expect_false(dm.hasInterLocusVariation(dm_tmp, 1))
  expect_true(dm.hasInterLocusVariation(dm_tmp, 2))  
})

test_that("test.printEmptyDM", {
  tmp_file <- tempfile()
  sink(tmp_file)
  dm <- dm.createDemographicModel(25:26, 100)
  print(dm)
  sink(NULL)
  unlink(tmp_file)
})

test_that("test.printGroupDM", {
  tmp_file <- tempfile()
  sink(tmp_file)
  print(dm.grp)
  sink(NULL)
  unlink(tmp_file)
})


test_that('setTrioMutationsRates works', {
  dm <- dm.setTrioMutationRates(dm_trios, '17', 'theta', group=2)
  expect_equal(nrow(searchFeature(dm, 'mutation', group=2)), 1)
  expect_equal(searchFeature(dm, 'mutation', group=2)$parameter, "17")
  expect_equal(nrow(searchFeature(dm, 'mutation_outer', group=2)), 1)
  expect_equal(searchFeature(dm, 'mutation_outer', group=2)$parameter, "theta")
})

test_that('adding Migration works', {
  dm <- dm.createDemographicModel(10:11, 100)
  dm <- dm.addSymmetricMigration(dm, 1, 5)
  dm <- dm.finalize(dm)
  expect_that(nrow(searchFeature(dm, 'migration')), is_equivalent_to(2))
  expect_equal(grep('M', paste(dm@options$ms.cmd, collapse=' ')), 1)
  
  dm <- dm.createDemographicModel(10:11, 100)
  dm <- dm.addMigration(dm, 1, 5, pop.from = 1, pop.to = 2)
  dm <- dm.finalize(dm)
  expect_that(nrow(searchFeature(dm, 'migration')), is_equivalent_to(1))
  expect_equal(grep('M', paste(dm@options$ms.cmd, collapse=' ')), 1)
})

test_that('getting the available Populations works', {
  dm <- dm.createDemographicModel(10:11, 100)
  expect_equal(getPopulations(dm), 1:2)
  expect_equal(getPopulations(dm.tt), 1:2)
  expect_equal(getPopulations(dm.hky), 1:3)
})