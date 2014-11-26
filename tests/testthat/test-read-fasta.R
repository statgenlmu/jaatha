context("Fasta input")

test_that("reading Fasta using ape works", {
  if (!test_ape) return()
  sample.file <- system.file("example_fasta_files/sample.fasta", 
                             package = "jaatha")
  sample.data <- read.dna(sample.file, format = "fasta", as.character = TRUE)
  sample.jsfs <- calculateJsfs(sample.data, pop1.rows = 3:7, 
                               pop2.rows = 8:12, outgroup.row = 1:2)
  expect_true(is.matrix(sample.jsfs))
  expect_equal(17, sum(sample.jsfs))
})