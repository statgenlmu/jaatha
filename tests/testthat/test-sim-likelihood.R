context("Final Likelihood estimation")

test_that("direct simulation of likelihoods works", {
    lh1 <- simLikelihood(jaatha.csi, 10, c(0.5, 0.5))
    expect_true(is.numeric(lh1))
    expect_true(lh1 != 0)
    
    lh2 <- simLikelihood(smooth_jaatha, 20, c(0.5, 0.5))
    expect_true(is.numeric(lh2))
    expect_true(lh2 != 0)
})

