context("PopGenome Interface")

test_that("getting segsites from PopGenome works for individual loci", {
  skip_if_not_installed("PopGenome")
  data_pg <- create_popgenome_test_data()
  
  seg_sites <- get_popgenome_segsites(data_pg, FALSE)
  expect_equal(length(seg_sites), 1)
  expect_true(coala::is_segsites(seg_sites[[1]]))
  expect_equal(nrow(seg_sites[[1]]), 12)
  
  expect_false(is.null(coala::get_positions(seg_sites[[1]])))
  expect_true(all(coala::get_positions(seg_sites[[1]]) >= 0)) 
  expect_true(all(coala::get_positions(seg_sites[[1]]) <= 1))
})


test_that("it imports PopGenome data", {
  skip_if_not_installed("PopGenome")
  skip_if_not_installed("coala")
  
  coala_model <- coala::coal_model(c(5, 5, 2), 1) +
    coala::feat_outgroup(3) +
    coala::feat_migration(coala::par_range("m", .1, 5), symmetric = TRUE) +
    coala::feat_mutation(coala::par_range("theta", .1, 5), model = "HKY", 
                         tstv_ratio = 1.6, 
                         base_frequencies = rep(.25, 4)) +
    coala::sumstat_sfs("sfs", population = 1) +
    coala::sumstat_jsfs("jsfs")
  jaatha_model <- create_jaatha_model(coala_model, test = FALSE, 
                                      jsfs_part = 1, jsfs_part_hi = 1)
  data_pg <- create_popgenome_test_data()
  
  jaatha_data <- create_jaatha_data(data_pg, jaatha_model, coala_model)
  expect_true(is_jaatha_data(jaatha_data))
  
  # Check sample size consistency
  expect_error(create_jaatha_data(data_pg, jaatha_model, 
                                  coala::coal_model(c(6, 4, 2), 1) +
                                    coala::feat_outgroup(3)))
})
