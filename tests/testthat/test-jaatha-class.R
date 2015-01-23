context("Jaatha class")

test_that("Custom initialization works", {
  jaatha <- new("Jaatha", csi.sim.func, csi.par.ranges, list(csi=csi.sum.stat), 2)
  expect_equal(length(jaatha@sum.stats), 1)
  expect_equal(jaatha@cores, 2)
  
  expect_true(all(jaatha@par.ranges == csi.par.ranges))
  expect_equal(3, length(jaatha@seeds))
})


test_that("PG initialization works", {
  # dm + jsfs
  checkType(Jaatha.initialize(sum.stats.tt, dm.tt), "jaatha")
  checkType(Jaatha.initialize(model = dm.tt, data = sum.stats.tt), "jaatha")

  expect_error(Jaatha.initialize(dm.tt))
  expect_error(Jaatha.initialize(NULL))
  
  # cores
  jaatha <- Jaatha.initialize(sum.stats.tt, dm.tt)
  expect_equal(jaatha@cores, 1)
  jaatha <- Jaatha.initialize(sum.stats.tt, dm.tt, cores = 2)
  expect_equal(jaatha@cores, 2)
  
  # scaling
  jaatha <- Jaatha.initialize(sum.stats.tt, dm.tt, scaling.factor = 5)
  expect_equal(jaatha@scaling.factor, 5)
  expect_equal(dm.getLociNumber(jaatha@opts[['dm']]), 1L)
    
  # folded
  jaatha <- Jaatha.initialize(sum.stats.tt, dm.tt, folded = TRUE)
  
  # smoothing
  jaatha <- Jaatha.initialize(sum.stats.mig, dm.mig, smoothing = TRUE)
  expect_equal(length(jaatha@sum.stats), 2)
  expect_that(jaatha@sum.stats$jsfs$get_data(), is_a('data.frame'))
  expect_equal(colnames(jaatha@sum.stats$jsfs$get_data()), c('X1', 'X2', 'sum.stat'))
  expect_that(sum(jaatha@sum.stats$jsfs$get_data()), is_more_than(0))
  expect_that(jaatha@sum.stats$jsfs$get_model(), is_a('character'))
  
  expect_that(jaatha@sum.stats$border_jsfs, is_a('Stat_PoiInd'))
  expect_that(sum(jaatha@sum.stats$border_jsfs$get_data()), is_more_than(0))
  
  expect_error(Jaatha.initialize(dm.tt, sum.stats.tt, 
                                 folded = TRUE,  smoothing = TRUE))
})


test_that("JSFS is added when not in model", {
  dm <- resetSumStats(dm.tt)
  expect_warning(j <- Jaatha.initialize(sum.stats.tt, dm))
  expect_equal(length(j@sum.stats), 1)
  expect_false(is.null(j@sum.stats[["jsfs"]]))
  expect_equal(dm.getSummaryStatistics(j@opts[['dm']]), as.factor(c('jsfs')))
})

test_that("PG initialization with groups works", {
  jaatha.grp <- Jaatha.initialize(sum.stats.grp, dm.grp)
  for (i in 1:3) {
    name <- getStatName('jsfs', i)
    expect_that(jaatha.grp@sum.stats[[name]], is_a('Stat_PoiInd'))
    expect_that(sum(jaatha.grp@sum.stats[[name]]$get_data()), is_more_than(0))
  }
    
  jaatha.grp <- Jaatha.initialize(sum.stats.grp, dm.grp, smoothing = TRUE)
  for (i in 1:3) {
    name <- getStatName('jsfs', i)
    expect_that(jaatha.grp@sum.stats[[name]], is_a('Stat_PoiSmooth'))
    expect_that(sum(jaatha.grp@sum.stats[[name]]$get_data()), is_more_than(0))
  }
})


test_that("PG initialization with FPC statistic", {
  # Without groups
  jaatha.fpc <- Jaatha.initialize(sum.stats.tt, dm.tt, use_fpc = TRUE)
  for (pop in 1:2) {
    fpc_stat <- jaatha.fpc@sum.stats[[getStatName('fpc',0,pop)]]
    expect_false(is.null(fpc_stat))
    expect_that(sum(fpc_stat$get_data()), is_more_than(0))
    expect_that(sum(fpc_stat$get_data()), is_less_than(dm.getLociNumber(dm.tt)+1))
  }
  
  jaatha.fpc <- Jaatha.initialize(sum.stats.tt, dm.tt, 
                                  use_fpc = TRUE, fpc_populations = 1)
  expect_false(is.null(jaatha.fpc@sum.stats[[getStatName('fpc',0,1)]]))
  expect_true(is.null(jaatha.fpc@sum.stats[[getStatName('fpc',0,2)]]))  
  
  jaatha.fpc <- Jaatha.initialize(sum.stats.tt, dm.tt,
                                  use_fpc = TRUE, fpc_populations = 2)
  expect_true(is.null(jaatha.fpc@sum.stats[["fpc_pop1"]]))
  expect_false(is.null(jaatha.fpc@sum.stats[["fpc_pop2"]]))
  
  # With groups
  jaatha.fpc <- Jaatha.initialize(sum.stats.grp, dm.grp, use_fpc = TRUE)
  expect_false(is.null(jaatha.fpc@sum.stats[[getStatName('fpc',1,1)]]))
  expect_false(is.null(jaatha.fpc@sum.stats[[getStatName('fpc',1,2)]]))
  expect_false(is.null(jaatha.fpc@sum.stats[[getStatName('fpc',2,1)]]))
  expect_false(is.null(jaatha.fpc@sum.stats[[getStatName('fpc',2,2)]]))
  expect_false(is.null(jaatha.fpc@sum.stats[[getStatName('fpc',3,1)]]))
  expect_false(is.null(jaatha.fpc@sum.stats[[getStatName('fpc',3,2)]]))
})

test_that("seg.sites are added when not in model", {
  j <- Jaatha.initialize(sum.stats.tt, dm.tt, use_fpc = TRUE)
  expect_true('seg.sites' %in% dm.getSummaryStatistics(j@opts[['dm']]))
})


# test_that("Initialization of PMC statistic", {
#   dm.pmc <- dm.addSummaryStatistic(dm.fpc, 'pmc')
#   
#   # Without groups
#   jaatha.pmc <- Jaatha.initialize(sum.stats.fpc, dm.pmc)
#   expect_true(sum(jaatha.pmc@sum.stats[["pmc"]]$value) > 0)
#   expect_false(is.null(jaatha.pmc@opts$dm@options[["pmc_breaks_private"]]))
#   expect_false(is.null(jaatha.pmc@opts$dm@options[["pmc_breaks_fixed"]]))  
#   
#   # With groups
#   dm.pmc <- dm.addSummaryStatistic(dm.grp, 'pmc', group=1)
#   dm.pmc <- dm.addSummaryStatistic(dm.pmc, 'pmc', group=3)
#   
#   jaatha.pmc <- Jaatha.initialize(sum.stats.grp, dm.pmc)
#   expect_equal(length(jaatha.pmc@sum.stats), 5)
#   expect_false(is.null(jaatha.pmc@sum.stats$pmc.1))
#   expect_true(sum(jaatha.pmc@sum.stats[["pmc.1"]]$value) > 0)
#   expect_true(is.null(jaatha.pmc@sum.stats$pmc.2))
#   expect_false(is.null(jaatha.pmc@sum.stats$pmc.3))
#   expect_true(sum(jaatha.pmc@sum.stats[["pmc.3"]]$value) > 0)
# })






test_that("getParNames", {
  expect_equal(dm.getParameters(dm.tt), getParNames(jaatha.tt))
  expect_equal(dm.getParameters(dm.mig), getParNames(jaatha.mig))
})

test_that("getParNumber", {
  expect_equal(dm.getNPar(dm.tt), getParNumber(jaatha.tt))
  expect_equal(dm.getNPar(dm.mig), getParNumber(jaatha.mig))
})

