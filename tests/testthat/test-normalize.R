context("Parameter scaling")

test_that("Parameter denormalization works", {
    expect_true(sum(abs(denormalize(c(0, 0), jaatha.csi) - c(0.1, 0.1))) < 1e-11)
    expect_true(sum(abs(denormalize(c(1, 1), jaatha.csi) - c(10, 10))) < 1e-11)
})

test_that("Parameter normalization works", {
    expect_true(all(normalize(c(0.1, 0.1), jaatha.csi) == c(0, 0)))
    expect_true(all(normalize(c(10, 10), jaatha.csi) == c(1, 1)))
})

test_that("Normalization and Denormalization are inverse", {
    for (x in 0:10/10) {
        pars <- rep(x, 2)
        expect_true(sum(abs(normalize(denormalize(x, jaatha.csi), 
            jaatha.csi) - pars)) < 1e-11)
    }
})

