context("Custom simulator")

test_that("customSimulator", {
    expect_true(length(jaatha.csi@starting.positions) == 4)
    estimates <- Jaatha.getLikelihoods(jaatha.csi)
    expect_true(all(estimates[, 1] != 0))
    expect_true(all(0.1 <= estimates[, 3] & estimates[, 4] <= 
        10))
})

