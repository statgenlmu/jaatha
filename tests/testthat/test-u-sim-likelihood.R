context("Final Likelihood estimation")

test_that("test.simLikelihood", {
    lh1 <- simLikelihood(jaatha.csi, 10, c(0.5, 0.5))
    expect_true(is.numeric(lh1))
    expect_true(lh1 != 0)
    lh2 <- simLikelihood(smooth.jaatha, 20, c(0.5, 0.5))
    expect_true(is.numeric(lh2))
    expect_true(lh2 != 0)
    lh3 <- simLikelihood(smooth.border.jaatha, 20, c(0.5, 0.5))
    expect_true(is.numeric(lh3))
    expect_true(lh3 != 0)
    expect_true(lh2 != lh3)
})

