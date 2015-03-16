context("Jaatha class")

test_that("Custom initialization works", {
  jaatha <- new("Jaatha", 
                csi.sim.func, 
                csi.par.ranges, 
                list(csi=csi.sum.stat), 2)
  expect_equal(length(jaatha@sum_stats), 1)
  expect_equal(jaatha@cores, 2)
  
  expect_true(all(jaatha@par.ranges == csi.par.ranges))
  expect_equal(3, length(jaatha@seeds))
})
