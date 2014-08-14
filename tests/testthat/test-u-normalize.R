context("Parameter scaling")

test_that("test.denormalize", {
    expect_true(sum(abs(denormalize(c(0, 0), jaatha.tt) - c(0.01, 
        1))) < 1e-11)
    expect_true(sum(abs(denormalize(c(1, 1), jaatha.tt) - c(5, 
        20))) < 1e-11)
})

test_that("test.normalize", {
    expect_true(all(normalize(c(0.01, 1), jaatha.tt) == c(0, 
        0)))
    expect_true(all(normalize(c(5, 20), jaatha.tt) == c(1, 1)))
})

test_that("test.normalizeDenormalize", {
    for (x in 0:10/10) {
        pars <- rep(x, 2)
        expect_true(sum(abs(normalize(denormalize(x, jaatha.tt), 
            jaatha.tt) - pars)) < 1e-11)
    }
})

