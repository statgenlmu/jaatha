test.addSegSitesToFpc <- function() {
  seg.sites <- matrix(c(1,1,0,0,0, 1,0,1,0,1, 1,1,1,1,0, 0,1,1,0,1, 0,0,0,0,1, 1,0,0,0,0), 5) 
  positions <- c(0.1, 0.12, 0.5, 0.51, 0.61)
  fpc <- matrix(0, 5, 5)
  fpc.new <- addSegSitesToFpc(seg.sites, positions, c(.25, .5, .75), c(.25, .5, .75), fpc)
  checkEquals(matrix(0, 5, 5), fpc)
  checkEquals(1, sum(fpc.new))
  checkEquals(1, fpc.new[2,2])
  fpc.new <- addSegSitesToFpc(seg.sites, positions, c(.25, .5, .75), c(.25, .5, .75), fpc.new)
  checkEquals(2, sum(fpc.new))
  checkEquals(2, fpc.new[2,2])

  fpc.new <- addSegSitesToFpc(seg.sites, positions, c(.25, .4, .75), c(.25, .5, .75), fpc.new)
  checkEquals(3, sum(fpc.new))
  checkEquals(2, fpc.new[2,2])
  checkEquals(1, fpc.new[3,2])

  fpc.new <- addSegSitesToFpc(seg.sites, positions, c(.25, .4, .75), c(.25, .4, .75), fpc.new)
  checkEquals(4, sum(fpc.new))
  checkEquals(2, fpc.new[2,2])
  checkEquals(1, fpc.new[3,2])
  checkEquals(1, fpc.new[3,3])

  singletons <- matrix(0,5,4)
  singletons[1,] <- 1
  fpc.new <- addSegSitesToFpc(singletons, positions, c(.25, .4, .75), c(.25, .4, .75), fpc)
  checkEquals(1, sum(fpc.new))
  checkEquals(1, fpc.new[5,5])

  positions <- 1:5/6
  fpc.new <- addSegSitesToFpc(seg.sites, positions, c(.25, .4, .75), c(.25, .4, .75), fpc)
  checkEquals(1, sum(fpc.new))
  checkEquals(1, fpc.new[5,3])
}

test.calcPercentFpcViolation <- function() {
  snp.matrix <- matrix(c(1,1,0,0,0,1,0,1,0,1,1,1,1,1,0,0,1,1,0,1,0,0,0,0,1), 5) 
  positions <- c(0.1, 0.12, 0.5, 0.51, 0.61)

  checkEquals(c(0.5, 0.5), calcPercentFpcViolation(snp.matrix, positions))
  positions[4:5] <- c(0.7, 0.75) 
  checkEquals(c(1, 0.4), calcPercentFpcViolation(snp.matrix, positions))

  # No near snps
  positions <- 1:5/5
  checkEquals(c(NA, 0.5), calcPercentFpcViolation(snp.matrix, positions))

  # Only singletons
  snp.matrix[1, ] <- 1
  snp.matrix[-1, ] <- 0
  checkTrue(all( is.na(calcPercentFpcViolation(snp.matrix, positions)) ))

  # No SNPs at all
  checkTrue(all( is.na(calcPercentFpcViolation(matrix(0, 5, 0), numeric(0))) ))
} 
