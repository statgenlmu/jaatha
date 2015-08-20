context("PopGenome Interface")

test_that("getting segsites from PopGenome works for individual loci", {
  skip_if_not_installed("PopGenome")
  data_pg <- create_popgenome_test_data()
  
  seg_sites <- get_popgenome_segsites(data_pg, FALSE, NULL)
  expect_is(seg_sites, "list")
  expect_equal(length(seg_sites), 1)
  expect_is(seg_sites[[1]], "matrix")
  expect_equal(nrow(seg_sites[[1]]), 12)
  expect_equal(grep("Individual_1", row.names(seg_sites[[1]])), 1:5)
  expect_equal(grep("Individual_2", row.names(seg_sites[[1]])), 6:10)  
  expect_equal(grep("Out", row.names(seg_sites[[1]])), 11:12)
  
  expect_false(is.null(attr(seg_sites[[1]], "positions")))
  expect_true(all(attr(seg_sites[[1]], "positions") >= 0)) 
  expect_true(all(attr(seg_sites[[1]], "positions") <= 1))
})


test_that("getting segsites from PopGenome works for trios", {
  skip_if_not_installed("PopGenome")
  data_pg <- create_popgenome_test_data()
  
  seg_sites <- get_popgenome_segsites(data_pg, FALSE,
                                      trios = list(rep(1, 3), rep(1, 3)))
  expect_is(seg_sites, "list")
  expect_equal(length(seg_sites), 2)
  expect_equal(seg_sites[[1]], seg_sites[[2]])
  expect_equal(ncol(seg_sites[[1]]), 48)
  
  expect_false(is.null(attr(seg_sites[[1]], "positions")))
  expect_true(all(attr(seg_sites[[1]], "positions") >= 0)) 
  expect_true(all(attr(seg_sites[[1]], "positions") <= 1))
  expect_equal(attr(seg_sites[[1]], "positions")[1:16],
               attr(seg_sites[[1]], "positions")[17:32])
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
  jaatha_model <- create_jaatha_model(coala_model, test = FALSE)
  data_pg <- create_popgenome_test_data()
  
  jaatha_data <- create_jaatha_data(data_pg, jaatha_model, coala_model)
  expect_true(is_jaatha_data(jaatha_data))
  
  expect_error(create_jaatha_data(data_pg, jaatha_model, 
                                  coala::coal_model(c(6, 4, 2), 1) +
                                    coala::feat_outgroup(3)))
  expect_error(create_jaatha_data(data_pg, jaatha_model, 
                                  coala::coal_model(c(5, 5, 3), 1) +
                                    coala::feat_outgroup(3)))
  expect_error(create_jaatha_data(data_pg, jaatha_model, 
                                  coala::coal_model(c(5, 5, 3), 1)))
})
