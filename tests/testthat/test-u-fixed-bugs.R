context("Fixed bugs")

test_that("test.ApeCodeInVignette", {
    require("ape") || stop("Package \"ape\" is required for unit tests")
    sample.file <- system.file("example_fasta_files/sample.fasta", 
        package = "jaatha")
    sample.data <- read.dna(sample.file, format = "fasta", as.character = TRUE)
    sample.jsfs <- calculateJsfs(sample.data, pop1.rows = 3:7, 
        pop2.rows = 8:12, outgroup.row = 1:2)
    expect_true(sum(sample.jsfs) > 0)
})

test_that("test.printEmptyDM", {
    tmp_file <- tempfile()
    sink(tmp_file)
    dm <- dm.createDemographicModel(25:26, 100)
    print(dm)
    sink(NULL)
    unlink(tmp_file)
})

test_that("test.printGroupDM", {
  tmp_file <- tempfile()
  sink(tmp_file)
  print(dm.grp)
  sink(NULL)
  unlink(tmp_file)
})