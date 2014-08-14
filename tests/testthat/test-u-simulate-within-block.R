context("Simulation (Jaatha side)")

test_that("test.runSimulatinos", {
    set.seed(15)
    pars.test <- matrix(0.5, 3, 2)
    sum.stats1 <- runSimulations(pars.test, 1, jaatha.csi)
    expect_equal(3, length(sum.stats1))
    for (i in 1:3) {
        expect_true(all(denormalize(pars.test[i, ], jaatha.csi) == 
            sum.stats1[[i]]$pars))
        expect_true(all(pars.test[i, ] == sum.stats1[[i]]$pars.normal))
        expect_false(is.null(sum.stats1[[i]]$poisson.vector))
        expect_true(sum(sum.stats1[[i]]$poisson.vector) > 0)
    }
})

test_that("test.simulateWithinBlock", {
    checkSumStat <- function(x, block) {
        if (is.null(x$pars) || is.null(x$pars.normal)) 
            return(FALSE)
        isInBlock(block, x$pars.normal)
    }
    sum.stats <- simulateWithinBlock(10, block.test, jaatha.csi)
    expect_true(is.list(sum.stats))
    expect_equal(14, length(sum.stats))
    expect_true(all(sapply(sum.stats, checkSumStat, block = block.test)))
    sum.stats <- simulateWithinBlock(2, block.test, jaatha.tt)
    expect_equal(6, length(sum.stats))
    expect_true(all(sapply(sum.stats, checkSumStat, block = block.test)))
})

