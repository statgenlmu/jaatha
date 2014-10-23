context("JSFS Summary Statistic")

test_that("calcJsfs works", {
  seg.sites <- list(matrix(c(1, 1, 0, 0, 0, 
                             1, 0, 1, 0, 1, 
                             1, 1, 1, 1, 0, 
                             0, 1, 1, 0, 1, 
                             0, 0, 0, 0, 1, 
                             1, 0, 0, 0, 0), 6, byrow=TRUE))
  
  dm <- dm.createDemographicModel(c(3,3), 1)
  jsfs <- calcJsfs(seg.sites, dm.getSampleSize(dm))
  expect_equal(dim(jsfs), c(4, 4))
  expect_equal(sum(jsfs), 5)
  expect_true(all(jsfs >= 0))
  expect_equal(jsfs[4, 2], 1)
  expect_equal(jsfs[3, 2], 2)
  expect_equal(jsfs[2, 1], 1)
  expect_equal(jsfs[2, 3], 1)
  
  seg.sites[[2]] <- seg.sites[[1]]
  dm <- dm.createDemographicModel(c(3,3), 2)
  jsfs <- calcJsfs(seg.sites, dm.getSampleSize(dm))
  expect_equal(dim(jsfs), c(4, 4))
  expect_equal(sum(jsfs), 10)
  expect_true(all(jsfs >= 0))
  expect_equal(jsfs[4, 2], 2)
  expect_equal(jsfs[3, 2], 4)
  expect_equal(jsfs[2, 1], 2)
  expect_equal(jsfs[2, 3], 2)
})