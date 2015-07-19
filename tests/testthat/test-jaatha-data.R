context("Jaatha Data")


test_that("default creation of data works", {
  model <- create_test_model()
  real_data <- rep(c(17, 5), 5)
  
  jaatha_data <- create_jaatha_data(real_data, model)
  expect_true(is_jaatha_data(jaatha_data))
  
  expect_equal(jaatha_data$get_values(), 
               list(id = real_data, sum = sum(real_data)))
  expect_equal(jaatha_data$get_values("id"), real_data)
  expect_equal(jaatha_data$get_values(stat_identity()), real_data)
  
  expect_equal(jaatha_data$get_options(), list())
  expect_equal(jaatha_data$get_options("id"), NULL)
  expect_equal(jaatha_data$get_options(stat_identity()), NULL)
})
