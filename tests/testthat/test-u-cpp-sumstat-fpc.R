context("Rcpp generate FGC statistic")

test_that("test.addSegSitesToFpc", {
    seg.sites <- matrix(c(1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 
        1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0), 
        5)
    positions <- c(0.1, 0.12, 0.5, 0.51, 0.61)
    fpc <- matrix(0, 5, 5)
    fpc.new <- addSegSitesToFpc(seg.sites, positions, c(0.25, 
        0.5, 0.75), c(0.25, 0.5, 0.75), fpc)
    expect_equal(fpc, matrix(0, 5, 5))
    expect_equal(sum(fpc.new), 1)
    expect_equal(fpc.new[2, 2], 1)
    fpc.new <- addSegSitesToFpc(seg.sites, positions, c(0.25, 
        0.5, 0.75), c(0.25, 0.5, 0.75), fpc.new)
    expect_equal(sum(fpc.new), 2)
    expect_equal(fpc.new[2, 2], 2)
    fpc.new <- addSegSitesToFpc(seg.sites, positions, c(0.25, 
        0.4, 0.75), c(0.25, 0.5, 0.75), fpc.new)
    expect_equal(sum(fpc.new), 3)
    expect_equal(fpc.new[2, 2], 2)
    expect_equal(fpc.new[3, 2], 1)
    fpc.new <- addSegSitesToFpc(seg.sites, positions, c(0.25, 
        0.4, 0.75), c(0.25, 0.4, 0.75), fpc.new)
    expect_equal(sum(fpc.new), 4)
    expect_equal(fpc.new[2, 2], 2)
    expect_equal(fpc.new[3, 2], 1)
    expect_equal(fpc.new[3, 3], 1)
    singletons <- matrix(0, 5, 4)
    singletons[1, ] <- 1
    fpc.new <- addSegSitesToFpc(singletons, positions, c(0.25, 
        0.4, 0.75), c(0.25, 0.4, 0.75), fpc)
    expect_equal(sum(fpc.new), 1)
    expect_equal(fpc.new[5, 5], 1)
    positions <- 1:5/6
    fpc.new <- addSegSitesToFpc(seg.sites, positions, c(0.25, 
        0.4, 0.75), c(0.25, 0.4, 0.75), fpc)
    expect_equal(sum(fpc.new), 1)
    expect_equal(fpc.new[5, 3], 1)
})

test_that("test.calcPercentFpcViolation", {
    snp.matrix <- matrix(c(1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 
        1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1), 5)
    positions <- c(0.1, 0.12, 0.5, 0.51, 0.61)
    expect_equal(calcPercentFpcViolation(snp.matrix, positions), 
        c(0.5, 0.5))
    positions[4:5] <- c(0.7, 0.75)
    expect_equal(calcPercentFpcViolation(snp.matrix, positions), 
        c(1, 0.4))
    positions <- 1:5/5
    expect_equal(calcPercentFpcViolation(snp.matrix, positions), 
        c(NA, 0.5))
    snp.matrix[1, ] <- 1
    snp.matrix[-1, ] <- 0
    expect_true(all(is.na(calcPercentFpcViolation(snp.matrix, 
        positions))))
    expect_true(all(is.na(calcPercentFpcViolation(matrix(0, 5, 
        0), numeric(0)))))
})

