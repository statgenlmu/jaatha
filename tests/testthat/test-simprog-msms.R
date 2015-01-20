context("msms simulation interface")

test_that("test.callMsms", {
  if (!test_msms) skip('msms not installed')
  jar.path = getJaathaVariable("msms.jar")
  ms.args <- "5 1 -r 10 100 -t 5 -I 2 3 2 1"
  msms.args <- ""
  set.seed(17)
  out.file <- callMsms(jar.path, ms.args, msms.args)
  set.seed(17)
  out.file.2 <- callMsms(jar.path, ms.args, msms.args)
  set.seed(20)
  out.file.3 <- callMsms(jar.path, ms.args, msms.args)
  expect_equal(file.info(out.file.2)$size, file.info(out.file)$size)
  expect_true(file.info(out.file)$size != file.info(out.file.3)$size)
  unlink(c(out.file, out.file.2, out.file.3))
})

test_that("test.generateMsmsOptionsCommand", {
  dm <- dm.addPositiveSelection(dm.tt, 100, 500, population = 1, 
                                at.time = "2")
  s <- 5
  opts <- paste(eval(parse(text = generateMsmsOptionsCommand(dm))), 
                collapse = " ")
  expect_equal(grep("-SI", opts), 1)
  expect_equal(grep("-SA 5", opts), 1)
})

test_that("generateMsmsOptionsCommand works with balancing selection", {
  dm <- dm.addBalancingSelection(dm.tt, 100, 500, population = 1, at.time = "2")
  s <- 5
  opts <- paste(eval(parse(text = generateMsmsOptionsCommand(dm))), 
                collapse = " ")

  expect_equal(grep("-SI", opts), 1)
  expect_equal(grep("-SAA 0", opts), 1)
  expect_equal(grep("-SAa 5", opts), 1)
})

test_that("test.msmsPrint", {
  if (!test_msms) skip('msms not installed')
  tmp_file <- tempfile()
  sink(tmp_file)
  print(dm.sel)
  sink(NULL)
  unlink(tmp_file)
})

test_that("msmsSimFunc works", {
  if (!test_msms) skip('msms not installed')
  set.seed(6688)
  sum_stats <- msmsSimFunc(dm.sel, c(1, 1.5, 1500, 5))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
  
  set.seed(6688)
  sum_stats2 <- msmsSimFunc(dm.sel, c(1, 1.5, 1500, 5))
  expect_equal(sum_stats, sum_stats2)
})

test_that("msmsSimFunc works for balacing selection", {
  if (!test_msms) skip('msms not installed')
  dm <- dm.addBalancingSelection(dm.tt, 100, 500, population = 1, at.time = "2")
  set.seed(6688)
  sum_stats <- msmsSimFunc(dm, c(1, 5, 300))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
  
  set.seed(6688)
  sum_stats2 <- msmsSimFunc(dm, c(1, 5, 300))
  expect_equal(sum_stats, sum_stats2)
})

test_that("msmsSimFunc works with inter-locus variation", {
  if (!test_msms) skip('msms not installed')
  dm_tmp <- dm.addInterLocusVariation(dm.sel)
  
  set.seed(1100)
  sum_stats <- msmsSimFunc(dm_tmp, c(1, 1.5, 1500, 5))
  expect_true(is.matrix(sum_stats$jsfs))
  expect_true(sum(sum_stats$jsfs) > 0)
  
  set.seed(1100)
  sum_stats2 <- msmsSimFunc(dm_tmp, c(1, 1.5, 1500, 5))
  expect_equal(sum_stats$jsfs, sum_stats2$jsfs)
})

test_that("Creation of parameter enviroment works", {
  par_envir <- createParameterEnv(dm.tt, c(1,5))
  expect_equal(par_envir[['tau']], 1)
  expect_equal(par_envir[['theta']], 5)
  
  par_envir <- createParameterEnv(dm.tt, c(1,5), locus = 17)
  expect_equal(par_envir[['locus']], 17)
  
  par_envir <- createParameterEnv(dm.tt, c(1,5), locus = 23, seed = 115)
  expect_equal(par_envir[['locus']], 23)
  expect_equal(par_envir[['seed']], 115)
})