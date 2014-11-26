context("Rcpp generate Polym Classes statistic")

test_that("Poly Classes are correct", {
    seg.sites <- matrix(c(1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1), 
        4, 3)
    expect_equal(classifyPolym(seg.sites, c(2, 2)), 
                c(private=1/3, fixed=0))
    expect_equal(classifyPolym(seg.sites, c(1, 3)), 
                 c(private=0, fixed=1/3))
    expect_equal(classifyPolym(seg.sites, c(3, 1)), 
                 c(private=1/3, fixed=0))
    
    seg.sites <- matrix(c(1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0), 4, 3)
    expect_equal(classifyPolym(seg.sites, c(2, 2)), 
                 c(private=0, fixed=2/3))
    expect_equal(classifyPolym(seg.sites, c(1, 3)), 
                 c(private=1/3, fixed=0))
    expect_equal(classifyPolym(seg.sites, c(3, 1)), 
                 c(private=1/3, fixed=0))
})