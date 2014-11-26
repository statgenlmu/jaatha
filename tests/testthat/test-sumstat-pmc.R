context("PMC Summary Statistic")

test_that("Poly Classes are correct", {
  seg.sites <- matrix(c(1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1), 
                      4, 3)
  expect_equal(classifyPolym(seg.sites, c(2, 2)), 
               c(private=1/3, fixed=0))
  expect_equal(classifyPolym(seg.sites, c(1, 3)), 
               c(private=0, fixed=1/3))
  expect_equal(classifyPolym(seg.sites, c(3, 1)), 
               c(private=1/3, fixed=0))
  
  seg.sites <- matrix(c(1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0), 4, 3)
  expect_equal(classifyPolym(seg.sites, c(2, 2)), 
               c(private=0, fixed=2/3))
  expect_equal(classifyPolym(seg.sites, c(1, 3)), 
               c(private=1/3, fixed=0))
  expect_equal(classifyPolym(seg.sites, c(3, 1)), 
               c(private=1/3, fixed=0))
})

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

test_that('PMC tranfomration works', {
  pmc <- matrix(1:12, 3, 4)
  expect_equal(transformPmc(pmc), c(1, 2, 4, 5, 7, 8, 12))
  
  pmc <- matrix(1:15 , 3, 5)
  expect_equal(transformPmc(pmc), c(1, 2, 4, 5, 7, 8, 10, 11, 15))

  pmc <- matrix(1:12 , 4, 3)
  expect_equal(transformPmc(pmc), c(1, 2, 3, 5, 6, 7, 12))
})
