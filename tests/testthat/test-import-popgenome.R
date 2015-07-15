context("PopGenome import")

create_popgenome_test_data <- function() {
  if (!requireNamespace("PopGenome", quietly = TRUE)) {
    stop("Please install package PopGenome")
  }
  
  # Create Test Data
  output <- tempfile("output")
  capture.output({
    fasta <- system.file("example_fasta_files", package = "jaatha")
    data_pg <- PopGenome::readData(fasta, progress_bar_switch = FALSE)
    data_pg <- PopGenome::set.outgroup(data_pg, c("Individual_Out-1", 
                                                  "Individual_Out-2"))
    data_pg <- PopGenome::set.populations(data_pg, 
                                          list(paste0("Individual_1-", 1:5), 
                                               paste0("Individual_2-", 1:5)))
  })
  data_pg
}


test_that("PopGenome data import works", {
  skip_if_not_installed("PopGenome")
  data_pg <- create_popgenome_test_data()
  
  seg_sites <- convPopGenomeToSegSites(data_pg)
  expect_is(seg_sites, "list")
  expect_equal(length(seg_sites), 1)
  expect_is(seg_sites$seg_sites, "list")
  expect_equal(length(seg_sites$seg_sites), 1)  
  expect_is(seg_sites$seg_sites[[1]], "matrix")
  expect_equal(nrow(seg_sites$seg_sites[[1]]), 12)
  expect_equal(grep("Individual_1", row.names(seg_sites$seg_sites[[1]])), 1:5)
  expect_equal(grep("Individual_2", row.names(seg_sites$seg_sites[[1]])), 6:10)  
  expect_equal(grep("Out", row.names(seg_sites$seg_sites[[1]])), 11:12)
  
  expect_false(is.null(attr(seg_sites$seg_sites[[1]], "positions")))
  expect_true(all(attr(seg_sites$seg_sites[[1]], "positions") >= 0)) 
  expect_true(all(attr(seg_sites$seg_sites[[1]], "positions") <= 1))
})


test_that("Trios can be imported from PopGenome", {
  skip_if_not_installed("PopGenome")
  data_pg <- create_popgenome_test_data()
  
  seg_sites_list <- convPopGenomeToSegSites(data_pg, trios = list(rep(1, 3), 
                                                                  rep(1, 3)))
  expect_is(seg_sites_list, "list")
  expect_is(seg_sites_list$seg_sites, "list")
  expect_equal(length(seg_sites_list$seg_sites), 2)
  expect_equal(seg_sites_list$seg_sites[[1]], seg_sites_list$seg_sites[[2]])
  expect_equal(ncol(seg_sites_list$seg_sites[[1]]), 48)
  
  expect_false(is.null(attr(seg_sites_list$seg_sites[[1]], "positions")))
  expect_true(all(attr(seg_sites_list$seg_sites[[1]], "positions") >= 0)) 
  expect_true(all(attr(seg_sites_list$seg_sites[[1]], "positions") <= 1))
  expect_equal(attr(seg_sites_list$seg_sites[[1]], "positions")[1:16],
               attr(seg_sites_list$seg_sites[[1]], "positions")[17:32])
})


test_that("PopGenome Model creation works", {
  skip_if_not_installed("PopGenome")
  data_pg <- create_popgenome_test_data()
  
  dm_pg <- createModelFromPopGenome(data_pg, quiet = TRUE)
  suppressMessages(dm_pg <- createModelFromPopGenome(data_pg))
  expect_equal(coala::get_sample_size(dm_pg), c(5,5,2))
  expect_equal(coala::get_outgroup(dm_pg), 3)
  expect_equal(coala::get_outgroup_size(dm_pg), 2)
  expect_equal(coala::get_locus_number(dm_pg), 1)
  expect_equal(coala:::get_locus_length(dm_pg, 1), 16)
})


test_that("Initialization with PopGenome-Data works", {
#   skip_on_cran()
#   skip_if_not_installed("PopGenome")
#   data_pg <- create_popgenome_test_data()
#   
#   if (!coala:::sg_find_exe(FALSE, TRUE)) skip("seqgen not installed")
# 
#   dm_pg <- createModelFromPopGenome(data_pg, quiet = TRUE) +
#     coala::feat_mutation(coala::par_range("theta", 1, 5), model = "HKY",
#                          tstv_ratio = 1.6, base_frequencies = rep(.25, 4)) +
#     coala::feat_pop_merge(coala::par_range("tau", .1, .5), 2, 1) +
#     coala::feat_migration(coala::par_const(.5), symmetric = TRUE) +
#     coala::feat_recombination(coala::par_const(.05)) +
#     coala::sumstat_jsfs()
#   
#   jaatha <- Jaatha.initialize(data_pg, dm_pg)
})
