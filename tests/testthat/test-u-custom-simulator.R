context("Custom simulator")

test_that("customSimulator", {
    expect_true(length(jaatha.csi@starting.positions) == 4)
})

