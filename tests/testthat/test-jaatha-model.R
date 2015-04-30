context("Jaatha Model")

test_that("jaatha model can be initialized", {
  stat <- R6::R6Class("Stat_PoiInd", inherit = jaatha:::Stat_Base, 
                      public = list(transform = function(data) sum(data)))
  
  create_jaatha_model(function(x, jaatha) rpois(20, x),
                      matrix(c(1, 2), 1, 2),
                      list(stat))
})