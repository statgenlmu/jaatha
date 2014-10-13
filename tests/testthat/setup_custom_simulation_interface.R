# A Block
block.test <- new("Block")
block.test@border <- matrix(c(0, 0, 0.1, 0.1), 2, 2) 


csi.sim.func <- function(x, jaatha) {
  list(poisson.vector=c(rpois(3, x[1]), rpois(3, x[2])))
}
csi.obs <- csi.sim.func(c(2:3))
csi.sum.stats <- list("poisson.vector"=list(method="poisson.independent",
                                            value=csi.obs$poisson.vector))
csi.par.ranges <- matrix(c(0.1, 0.1, 10, 10), 2, 2)
rownames(csi.par.ranges) <- c('x', 'y')
jaatha.csi <- new("Jaatha", csi.sim.func, csi.par.ranges, csi.sum.stats, 
                  123, cores = 2)
sim.data.csi <- jaatha:::simulateWithinBlock(10, block.test, jaatha.csi)

# A Smoothing Model
smooth.func <- function(x) {
  # x[1:2] = Matrix cell 
  # x[3:4] = Model Parameters
  x[1] * x[3]^2 + x[2] * x[4]^1.5 + prod(log(x[1:2]))
}

smooth.simfunc <- function(x, jaatha) {
  stopifnot(length(x) == 2)
  idxs <- cbind(rep(1:10, each=12), rep(1:12, 10), x[1], x[2])
  sampled.values <- sapply(apply(idxs, 1, smooth.func), function(x) rpois(1, x))
  list(mat=matrix(sampled.values, 10, 12, byrow=TRUE))
}

smooth.obs <- smooth.simfunc(c(3, 4))
smooth.sum.stats <- list("mat"=list(method="poisson.smoothing",
                                    model="(X1^2)*(X2^2)+log(X1)*log(X2)",
                                    value=smooth.obs$mat))

border.mask <- smooth.obs$mat
border.mask[, ] <- 0
border.mask[c(1, nrow(smooth.obs$mat)), ] <- 1
border.mask[ ,c(1, ncol(smooth.obs$mat))] <- 1
border.mask <- as.logical(border.mask)

border.transformation <- function(x) {
  c(x[1, 1:12], x[2:9, 1], x[10, 1:12], x[2:9, 12])
}

smooth.border.sum.stats <- list("mat"=list(method="poisson.smoothing",
                                           model="I(X1^2)*I(X2^2)+log(X1)*log(X2)",
                                           value=smooth.obs$mat,
                                           border.transformation=border.transformation,
                                           border.mask=border.mask))
rm(border.transformation, border.mask)

smooth.par.ranges <- matrix(c(2, 1, 7, 7), 2, 2)
rownames(smooth.par.ranges) <- c('x', 'y')

smooth.jaatha <- new("Jaatha", smooth.simfunc, 
                     smooth.par.ranges, smooth.sum.stats, 123)

smooth.border.jaatha <- new("Jaatha", smooth.simfunc, smooth.par.ranges, 
                            smooth.border.sum.stats, 123)
smooth.sim.data <- jaatha:::simulateWithinBlock(10, block.test, smooth.jaatha)