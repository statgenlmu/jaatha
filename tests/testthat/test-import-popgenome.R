context("PopGenome import")

# Create Test Data
output <- tempfile("output")
sink(output)
data_pg <- PopGenome::readData(system.file('example_fasta_files',  package='jaatha'), 
                               progress_bar_switch = FALSE)
data_pg <- PopGenome::set.outgroup(data_pg, c("Individual_Out-1", "Individual_Out-2"))
data_pg <- PopGenome::set.populations(data_pg, list(paste0("Individual_1-", 1:5), 
                                                    paste0("Individual_2-", 1:5)))
sink(NULL)
unlink(output)


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
  dm_pg <- createModelFromPopGenome(data_pg, quiet=TRUE)
  suppressMessages(dm_pg <- createModelFromPopGenome(data_pg))
  expect_equal(coalsimr::get_sample_size(dm_pg), c(5,5,2))
  expect_equal(coalsimr::get_outgroup(dm_pg), 3)
  expect_equal(coalsimr::get_outgroup_size(dm_pg), 2)
  expect_equal(coalsimr::get_locus_number(dm_pg), 1)
  expect_equal(coalsimr::get_locus_length(dm_pg), 16)
})


test_that("Initialization with PopGenome-Data works", {
  dm_pg <- createModelFromPopGenome(data_pg, quiet = TRUE) +
    coalsimr::feat_mutation(coalsimr::par_range('theta', 1, 5), model = 'HKY') +
    coalsimr::feat_pop_merge(coalsimr::par_range('tau', .1, .5), 2, 1) +
    coalsimr::feat_migration(coalsimr::par_const(.5), symmetric = TRUE) +
    coalsimr::feat_recombination(coalsimr::par_const(.05)) +
    coalsimr::sumstat_jsfs()
  
  jaatha <- Jaatha.initialize(data_pg, dm_pg)
})