context("Jaatha Model")

test_that("jaatha model can be initialized", {
  model <- create_test_model()
})


test_that("simulation works", {
  model <- create_test_model()
  res <- model$simulate(pars = c(1, 1), seed = 1)
  expect_that(res, is_a("list"))
  expect_equivalent(res$pars, c(10, 10))
  expect_equivalent(res$pars_normal, c(1, 1))
})
