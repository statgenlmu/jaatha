context("Simulation Cache")

test_that("caching of simulations works", {
  sc <- create_sim_cache()
  expect_equal(sc$get_size(), 0)
  
  sim_data <- lapply(1:10 / 10, function(x) list(pars_normal = rep(x, 3)))
  sc$add(sim_data)
  expect_equal(sc$get_size(), 10)
  block_all <- create_block(matrix(c(0, 1), 3, 2, byrow = TRUE))
  expect_equal(sc$get_sim_data(block_all), sim_data)
  block_small <- create_block(matrix(c(0, .25), 3, 2, byrow = TRUE))
  expect_equal(sc$get_sim_data(block_small), sim_data[1:2])
  
  sc$add(sim_data)
  expect_equal(sc$get_size(), 20)
  expect_equal(sc$get_sim_data(block_small), sim_data[c(1:2, 1:2)])
  
  sc$add(sim_data)
  expect_equal(sc$get_size(), 30)
  expect_equal(sc$get_sim_data(block_small), sim_data[c(1:2, 1:2, 1:2)])
  
  expect_error(sc$add(5))
  capture.output(expect_error(sc$add(list(5))))
  capture.output(expect_error(sc$add(list(list(5)))))
  
  # One parameter
  sc <- create_sim_cache()
  sim_data <- lapply(1:10 / 10, function(x) list(pars_normal = x))
  sc$add(sim_data)
  expect_equal(sc$get_sim_data(create_block(matrix(0:1, 1, 2))), sim_data)
})
