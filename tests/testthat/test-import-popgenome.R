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

test_that("PopGenome Model creation works", {
  dm_pg <- dm.createModelFromPopGenome(data_pg)
  expect_equal(dm.getSampleSize(dm_pg), c(5,5))
  expect_equal(dm.getLociNumber(dm_pg), 1)
  expect_equal(dm.getLociLength(dm_pg), 16)
})

test_that("Initialization with PopGenome-Data works", {
  if (!test_seqgen) skip('seq-gen not installed')
  dm_pg <- dm.createModelFromPopGenome(data_pg)
  dm_pg <- dm.addMutation(dm_pg, 1, 5)
  dm_pg <- dm.addSpeciationEvent(dm_pg, .1, .5, 'tau', 1, 2)
  
  # Outgroup is missing
  expect_error(jaatha <- Jaatha.initialize(data_pg, dm_pg))
  
  # Wrong outgroup size
  expect_error(Jaatha.initialize(data_pg, dm.addOutgroup(dm_pg, "2*tau", 3)))
  
  dm_pg <- dm.addOutgroup(dm_pg, "2*tau", sample_size = 2)
  jaatha <- Jaatha.initialize(data_pg, dm_pg)
})