context("PopGenome import")

test_that("PopGenome data import works", {  
  seg_sites <- convPopGenomeToSegSites(data_pg)
  expect_is(seg_sites, "list")
  expect_equal(length(seg_sites), 1)
  expect_is(seg_sites$seg.sites, "list")
  expect_equal(length(seg_sites$seg.sites), 1)  
  expect_is(seg_sites$seg.sites[[1]], "matrix")
  expect_equal(nrow(seg_sites$seg.sites[[1]]), 12)
  expect_equal(grep("Individual_1", row.names(seg_sites$seg.sites[[1]])), 1:5)
  expect_equal(grep("Individual_2", row.names(seg_sites$seg.sites[[1]])), 6:10)  
  expect_equal(grep("Out", row.names(seg_sites$seg.sites[[1]])), 11:12)
  
  expect_false(is.null(attr(seg_sites$seg.sites[[1]], "positions")))
  expect_true(all(attr(seg_sites$seg.sites[[1]], "positions") >= 0)) 
  expect_true(all(attr(seg_sites$seg.sites[[1]], "positions") <= 1))
})
