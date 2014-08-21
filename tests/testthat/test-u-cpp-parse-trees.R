context("Rcpp tree parsing")

test_that("Trees are extracted from simulation output", {
  sim_output <- tempfile("sim_output")
  tree_file <- tempfile("tree_file")
  
  cat("ms 3 1 -t 1 -r 1 20 -T 
30461 15911 34727

//
[2](2:0.865,(1:0.015,3:0.015):0.850);
[3](2:0.865,(1:0.015,3:0.015):0.850);
[4](2:1.261,(1:0.015,3:0.015):1.246);
[11](2:1.261,(1:0.015,3:0.015):1.246);
segsites: 5
positions: 0.2046 0.2234 0.2904 0.6209 0.9527 
01100
10011
01100

//
[2](3:0.613,(1:0.076,2:0.076):0.537);
[18](3:0.460,(1:0.076,2:0.076):0.384);
segsites: 2
positions: 0.3718 0.8443 
01
01
10
", file=sim_output);
  
  expect_equal(parseTrees(sim_output, tree_file, NA), tree_file)
  expect_true(file.exists(tree_file))
  trees <- scan(tree_file, "character", quiet=TRUE)
  expect_equal(length(trees), 6)
  expect_equal(trees[1], '[2](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[2], '[3](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[3], '[4](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[4], '[11](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[5], '[2](3:0.613,(1:0.076,2:0.076):0.537);')
  expect_equal(trees[6], '[18](3:0.460,(1:0.076,2:0.076):0.384);')
  unlink(tree_file)
  
  parseTrees(sim_output, tree_file, c(2, 4, 8, 2, 4))
  expect_true(file.exists(tree_file))
  trees <- scan(tree_file, "character", quiet=TRUE)
  expect_equal(length(trees), 7)
  expect_equal(trees[1], '[2](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[2], '[3](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[3], '[5](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[4], '[4](2:1.261,(1:0.015,3:0.015):1.246);')
  
  expect_equal(trees[5], '[2](3:0.613,(1:0.076,2:0.076):0.537);')
  expect_equal(trees[6], '[8](3:0.460,(1:0.076,2:0.076):0.384);')
  expect_equal(trees[7], '[4](3:0.460,(1:0.076,2:0.076):0.384);')
  
  parseTrees(sim_output, tree_file, c(9, 2, 2, 5, 2))
  expect_true(file.exists(tree_file))
  trees <- scan(tree_file, "character", quiet=TRUE)
  expect_equal(length(trees), 9)
  expect_equal(trees[1], '[2](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[2], '[3](2:0.865,(1:0.015,3:0.015):0.850);')
  expect_equal(trees[3], '[4](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[4], '[2](2:1.261,(1:0.015,3:0.015):1.246);')
  expect_equal(trees[5], '[2](2:1.261,(1:0.015,3:0.015):1.246);')
  
  expect_equal(trees[6], '[2](3:0.613,(1:0.076,2:0.076):0.537);')
  expect_equal(trees[7], '[7](3:0.460,(1:0.076,2:0.076):0.384);')  
  expect_equal(trees[8], '[2](3:0.460,(1:0.076,2:0.076):0.384);')
  expect_equal(trees[9], '[2](3:0.460,(1:0.076,2:0.076):0.384);')
  
  expect_error(parseTrees(sim_output, tree_file, c(9, 2, 2, 5, 1)))
  expect_error(parseTrees(sim_output, tree_file, c(19, 2, 2, 5, 1)))
})