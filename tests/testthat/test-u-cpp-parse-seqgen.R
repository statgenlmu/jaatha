context("Rcpp output parsing - seqgen")

test_that("test.parseOutputSeqgen", {
  if (!test_seqgen) return()
  sg.file <- dm.simSumStats(dm.hky, c(1, 5))$file['seqgen']
  expect_error(parseOutput(tempfile("seqgen2_"), dm.getSampleSize(dm.hky), 
                           dm.getLociNumber(dm.hky), 1))
  expect_error(parseOutput(sg.file, dm.getSampleSize(dm.hky), 
                           dm.getLociNumber(dm.hky) - 1, 1))
  expect_error(parseOutput(sg.file, dm.getSampleSize(dm.hky), 
                           dm.getLociNumber(dm.hky) + 1, 1))
  sum.stats = parseOutput(sg.file, dm.getSampleSize(dm.hky), 
                          dm.getLociNumber(dm.hky), 1)
  expect_true(is.matrix(sum.stats$jsfs))
  expect_true(sum(sum.stats$jsfs) > 0)
  expect_true(sum.stats$jsfs[1, 1] == 0)
  expect_true(sum.stats$jsfs[12, 13] == 0)
  sum.stats = parseOutput(sg.file, dm.getSampleSize(dm.hky), 
                          dm.getLociNumber(dm.hky), 1, generate_seg_sites = TRUE)
  expect_true(is.matrix(sum.stats$jsfs))
  expect_true(is.list(sum.stats$seg.sites))
  expect_equal(length(sum.stats$seg.sites), dm.getLociNumber(dm.hky))
  sum.stats = parseOutput(sg.file, dm.getSampleSize(dm.hky), 
                          dm.getLociNumber(dm.hky), 1, generate_fpc = TRUE, fpc_breaks_near = 1:3/4, 
                          fpc_breaks_far = 1:3/4)
  expect_true(is.matrix(sum.stats$fpc))
  expect_equal(sum(sum.stats$fpc), 5)
  unlink(sg.file)
})