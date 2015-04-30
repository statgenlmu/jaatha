context("Jaatha Model")

test_that("jaatha model can be initialized", {
  stat <- R6::R6Class("Stat_PoiInd", inherit = jaatha:::Stat_Base, 
                      public = list(transform = function(data) sum(data)))
  
  create_jaatha_model(function(x, jaatha) rpois(20, x),
                      matrix(c(1, 2), 1, 2),
                      list(stat))
})

test_that("simulation works", {
  model <- create_test_model()
  res <- model$simulate(pars = c(1, 1), seed = 1)
  expect_that(res, is_a("list"))
  expect_equivalent(res$pars, c(10, 10))
  expect_equivalent(res$pars_normal, c(1, 1))
})
