context('Parameter Class')

test_that('Getting and Setting Expressions works', {
  expect_error(Parameter$new(2*x))
  expect_error(Parameter$new(2))
  expect_error(Parameter$new('2'))
  basic_par <- Parameter$new(expression(2*x))
  x <- 5
  expect_equal(basic_par$eval(), 10)
  
  test_env <- new.env()
  test_env[['x']] <- 6
  expect_equal(basic_par$eval(envir = test_env), 12)
  expect_equal(basic_par$eval(), 10)
  
  expr <- basic_par$get_expression()
  expect_is(expr, 'expression')
  expect_equal(eval(expr), 10)
})


test_that('Getting and Setting Names work', {
  basic_par <- Parameter$new(expression(2*x), 'blub')
  expect_equal(basic_par$get_name(), 'blub')
  
  expect_error(Parameter$new(expression(2*x), 1))
  suppressWarnings(expect_error(Parameter$new(expression(2*x), mean)))
})


test_that('par_expr works', {
  basic_par <- par_expr(2*x)
  expect_is(basic_par, 'Parameter')
  x <- 5
  expect_equal(basic_par$eval(), 10)
  x <- 6
  expect_equal(basic_par$eval(), 12)
  
  basic_par <- par_expr(sqrt(y))
  y <- 4
  expect_equal(basic_par$eval(), 2)
})