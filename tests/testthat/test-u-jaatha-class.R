context("Jaatha class")

test_that("JSFS is added when not in model", {
  dm <- resetSumStats(dm.tt)
  expect_warning(j <- Jaatha.initialize(dm, sum.stats.tt$jsfs))
  expect_equal(length(j@sum.stats), 1)
  expect_true(is.list(j@sum.stats[["jsfs"]]))
  expect_equal(dm.getSummaryStatistics(j@opts[['dm']]), as.factor(c('jsfs')))
})


test_that("Initialization of FPC statistic", {
  expect_error(Jaatha.initialize(dm.fpc, sum.stats.fpc$jsfs))
       
  # Without groups
  jaatha.fpc <- Jaatha.initialize(dm.fpc, sum.stats.fpc)
  expect_true(sum(jaatha.fpc@sum.stats[["fpc"]]$value) > 0)
  expect_false(is.null(jaatha.fpc@opts$dm@options[["fpc.breaks.near"]]))
  expect_false(is.null(jaatha.fpc@opts$dm@options[["fpc.breaks.far"]]))
  
  # With groups
  dm.fpc <- dm.addSummaryStatistic(dm.grp, 'fpc')
  jaatha.fpc <- Jaatha.initialize(dm.fpc, sum.stats.grp, 123)
  expect_equal(length(jaatha.fpc@sum.stats), 6)
  expect_false(is.null(jaatha.fpc@sum.stats$fpc.1))
  expect_true(sum(jaatha.fpc@sum.stats[["fpc.1"]]$value) > 0)
  expect_false(is.null(jaatha.fpc@sum.stats$fpc.2))
  expect_true(sum(jaatha.fpc@sum.stats[["fpc.2"]]$value) > 0)
  expect_false(is.null(jaatha.fpc@sum.stats$fpc.3))
  expect_true(sum(jaatha.fpc@sum.stats[["fpc.3"]]$value) > 0)
})


test_that("Initialization of PMC statistic", {
  dm.pmc <- dm.addSummaryStatistic(dm.fpc, 'pmc')
  expect_error(Jaatha.initialize(dm.pmc, sum.stats.fpc$jsfs))
  
  # Without groups
  jaatha.pmc <- Jaatha.initialize(dm.pmc, sum.stats.fpc)
  expect_true(sum(jaatha.pmc@sum.stats[["fpc"]]$value) > 0)
  expect_true(sum(jaatha.pmc@sum.stats[["pmc"]]$value) > 0)
  expect_false(is.null(jaatha.pmc@opts$dm@options[["pmc_breaks_private"]]))
  expect_false(is.null(jaatha.pmc@opts$dm@options[["pmc_breaks_fixed"]]))  
  
  # With groups
  dm.pmc <- dm.addSummaryStatistic(dm.grp, 'pmc', 1)
  dm.pmc <- dm.addSummaryStatistic(dm.pmc, 'pmc', 3)
  
  jaatha.pmc <- Jaatha.initialize(dm.pmc, sum.stats.grp, 123)
  expect_equal(length(jaatha.pmc@sum.stats), 5)
  expect_false(is.null(jaatha.pmc@sum.stats$pmc.1))
  expect_true(sum(jaatha.pmc@sum.stats[["pmc.1"]]$value) > 0)
  expect_true(is.null(jaatha.pmc@sum.stats$pmc.2))
  expect_false(is.null(jaatha.pmc@sum.stats$pmc.3))
  expect_true(sum(jaatha.pmc@sum.stats[["pmc.3"]]$value) > 0)
})


test_that("Jaatha.initialization.groups", {
  sum.stats.grp <- dm.simSumStats(dm.grp, c(1, 5))
  jaatha.grp <- Jaatha.initialize(dm.grp, sum.stats.grp, 123, 
                                  folded = FALSE, smoothing = FALSE)
  expect_true(is.list(jaatha.grp@sum.stats[["jsfs.1"]]))
  expect_true(jaatha.grp@sum.stats[["jsfs.1"]]$method == "poisson.transformed")
  expect_true(sum(jaatha.grp@sum.stats[["jsfs.1"]]$value) > 
                0)
  expect_true(is.list(jaatha.grp@sum.stats[["jsfs.2"]]))
  expect_true(jaatha.grp@sum.stats[["jsfs.2"]]$method == "poisson.transformed")
  expect_true(sum(jaatha.grp@sum.stats[["jsfs.2"]]$value) > 
                0)
  expect_true(is.list(jaatha.grp@sum.stats[["jsfs.3"]]))
  expect_true(jaatha.grp@sum.stats[["jsfs.3"]]$method == "poisson.transformed")
  expect_true(sum(jaatha.grp@sum.stats[["jsfs.3"]]$value) > 
                0)
  suppressWarnings(
    jaatha.grp <- Jaatha.initialize(dm.grp, sum.stats.grp, 123, 
                                    folded = FALSE, smoothing = TRUE)
  )
  expect_true(is.list(jaatha.grp@sum.stats[["jsfs.1"]]))
  expect_true(jaatha.grp@sum.stats[["jsfs.1"]]$method == "poisson.smoothing")
  expect_true(sum(jaatha.grp@sum.stats[["jsfs.1"]]$value) > 
                0)
  expect_true(is.list(jaatha.grp@sum.stats[["jsfs.2"]]))
  expect_true(jaatha.grp@sum.stats[["jsfs.2"]]$method == "poisson.smoothing")
  expect_true(sum(jaatha.grp@sum.stats[["jsfs.2"]]$value) > 
                0)
  expect_true(is.list(jaatha.grp@sum.stats[["jsfs.3"]]))
  expect_true(jaatha.grp@sum.stats[["jsfs.3"]]$method == "poisson.smoothing")
  expect_true(sum(jaatha.grp@sum.stats[["jsfs.3"]]$value) > 
                0)
})

test_that("JaathaInitialize", {
  jsfs <- sum.stats.tt$jsfs
  
  # dm + jsfs
  checkType(Jaatha.initialize(dm.tt, jsfs), "jaatha")
  checkType(Jaatha.initialize(demographic.model = dm.tt, jsfs = jsfs), 
            "jaatha")
  checkType(Jaatha.initialize(dm.tt, sum.stats.tt), "jaatha")
  expect_error(Jaatha.initialize(dm.tt))
  expect_error(Jaatha.initialize(NULL))
  expect_error(Jaatha.initialize(jsfs = NULL))
  
  # seed
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 123)
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, seed = 123)
  expect_true(jaatha@seeds[1] == 123)
  
  # cores
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, cores = 1)
  expect_true(jaatha@cores == 1)
  
  # scaling
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, scaling.factor = 17)
  expect_true(jaatha@opts[["scaling.factor"]] == 17)
  
  # use.shm
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, use.shm = FALSE)
  suppressWarnings(
    jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, use.shm = TRUE)
  )
  
  # folded
  # jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, folded = TRUE)
  
  
  # smoothing
  suppressWarnings(
    jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, smoothing = TRUE)
  )
  expect_true(jaatha@sum.stats$jsfs$method == "poisson.smoothing")
  expect_false(is.null(jaatha@sum.stats$jsfs$model))
  checkType(jaatha@sum.stats$jsfs$model, "char")    
  expect_error(Jaatha.initialize(dm.tt, sum.stats.tt, folded = TRUE, 
                                 smoothing = TRUE))
})

test_that("getParNames", {
  expect_equal(dm.getParameters(dm.tt), getParNames(jaatha.tt))
  expect_equal(dm.getParameters(dm.mig), getParNames(jaatha.mig))
})

test_that("getParNumber", {
  expect_equal(dm.getNPar(dm.tt), getParNumber(jaatha.tt))
  expect_equal(dm.getNPar(dm.mig), getParNumber(jaatha.mig))
})

test_that("init", {
  sum.stat1 <- list(method = "poisson.independent", value = 1:10)
  sum.stat2 <- list(method = "poisson.transformed", 
                    value = matrix(1, 3, 3), transformation = diag)
  sum.stat3 <- list(method = "poisson.smoothing", model = "a+b", 
                    value = matrix(1, 3, 3))
  sum.stats <- list(ss1 = sum.stat1, ss2 = sum.stat2, ss3 = sum.stat3)
  jaatha <- new("Jaatha", csi.sim.func, csi.par.ranges, csi.sum.stats)
  jaatha <- new("Jaatha", csi.sim.func, csi.par.ranges, sum.stats)
  jaatha <- new("Jaatha", csi.sim.func, csi.par.ranges, sum.stats, 
                123, 1, FALSE)
  expect_true(all(jaatha@par.ranges == csi.par.ranges))
  expect_equal(3, length(jaatha@seeds))
  expect_equal(123, jaatha@seeds[1])
  expect_equal(1, jaatha@cores)
  expect_equal(FALSE, jaatha@use.shm)
  expect_true(length(jaatha@sum.stats) == 3)
  expect_false(is.null(jaatha@sum.stats$ss1))
  expect_false(is.null(jaatha@sum.stats$ss2))
  expect_false(is.null(jaatha@sum.stats$ss3))
  expect_true(all(jaatha@sum.stats$ss1$value == 1:10))
  expect_false(is.null(jaatha@sum.stats$ss1$transformation))
  expect_true(all(jaatha@sum.stats$ss2$value == 1))
  expect_false(is.null(jaatha@sum.stats$ss2$transformation))
  expect_true(all(jaatha@sum.stats$ss2$value.transformed == c(1, 1, 1)))
  expect_true(all(jaatha@sum.stats$ss3$value == 1))
})