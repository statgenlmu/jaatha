context("Refined Search")

test_that("test.refinedSearch.csi", {
    jaatha <- Jaatha.refinedSearch(jaatha.csi, 2, sim = 10, sim.final = 10)
    expect_true(is.matrix(jaatha@likelihood.table))
    expect_true(ncol(jaatha@likelihood.table) == 4)
    expect_true(nrow(jaatha@likelihood.table) >= 10)
    jaatha2 <- Jaatha.refinedSearch(jaatha.csi, 2, sim = 10, 
        sim.final = 10)
    expect_equal(jaatha2@likelihood.table, jaatha@likelihood.table)
})

test_that("test.refinedSearch.dm", {
    jaatha <- Jaatha.initialSearch(jaatha.tt, 20, 1)
    jaatha <- Jaatha.refinedSearch(jaatha, 1, 10, sim.final = 10, 
        max.steps = 10)
    expect_true(is.matrix(jaatha@likelihood.table))
    expect_true(ncol(jaatha@likelihood.table) == 4)
    expect_true(nrow(jaatha@likelihood.table) >= 9)
    jaatha2 <- Jaatha.initialSearch(jaatha.tt, 20, 1)
    jaatha2 <- Jaatha.refinedSearch(jaatha2, 1, 10, sim.final = 10, 
        max.steps = 10)
    expect_equal(jaatha2@likelihood.table, jaatha@likelihood.table)
})

