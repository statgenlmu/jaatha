context("PMC Summary Statistic")

test_that("generatePolyClasses workds", {
  dm <- dm.createDemographicModel(c(2,2), 1)
  
  seg_sites <- list(matrix(c(1, 0, 0, 0, 
                             1, 1, 0, 1, 
                             1, 1, 0, 0,
                             0, 0, 1, 1), 4, 4))
  
  expect_error(createPolymClasses(seg_sites, dm))
  
  dm@options['pmc_breaks_private'] <- 0.5
  dm@options['pmc_breaks_fixed'] <- 0.3
  pmc <- createPolymClasses(seg_sites, dm)
  expect_equal(sum(abs(pmc)), 1)
  expect_equal(pmc[1, 2], 1)
  
  dm@options$pmc_breaks_private <- 0.1
  dm@options$pmc_breaks_fixed <- c(0.1, 0.6)
  pmc <- createPolymClasses(seg_sites, dm)
  expect_equal(sum(abs(pmc)), 1)
  expect_equal(pmc[2, 2], 1)
  
  
  seg_sites <- list(seg_sites[[1]],
                    matrix(c(1, 0, 0, 0, 
                             1, 0, 0, 1, 
                             1, 0, 0, 0,
                             0, 0, 0, 1), 4, 4))
  pmc <- createPolymClasses(seg_sites, dm)
  expect_equal(sum(abs(pmc)), 2)
  expect_equal(pmc[2, 2], 1)
  expect_equal(pmc[2, 1], 1)
  
  
  seg_sites <- list(seg_sites[[1]], seg_sites[[1]],
                    matrix(0, 4, 0))
  pmc <- createPolymClasses(seg_sites, dm)
  expect_equal(sum(abs(pmc)), 3)
  expect_equal(pmc[2, 2], 2)
  expect_equal(pmc[3, 4], 1)
})