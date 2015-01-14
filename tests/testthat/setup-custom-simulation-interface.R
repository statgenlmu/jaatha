set.seed(111222555)

# A Block
block.test <- new("Block")
block.test@border <- matrix(c(0.4, 0.4, 0.6, 0.6), 2, 2) 

csi.sim.func <- function(x, jaatha) {
  list(data=c(rpois(5, x[1]), rpois(5, x[2])))
}
csi.obs <- csi.sim.func(c(3,5))
csi.sum.stat <- R6::R6Class("Stat_PoiInd", inherit = jaatha:::Stat_Base)$new(csi.obs)
csi.par.ranges <- matrix(c(0.1, 0.1, 10, 10), 2, 2)
rownames(csi.par.ranges) <- c('x', 'y')
jaatha.csi <- new("Jaatha", csi.sim.func, csi.par.ranges, list(csi=csi.sum.stat), 2)
sim.data.csi <- jaatha:::simulateWithinBlock(10, block.test, jaatha.csi)

# A Smoothing Model
smooth_func <- function(x) {
  # x[1:2] = Matrix cell 
  # x[3:4] = Model Parameters
  x[1] * x[3]^2 + x[2] * x[4]^1.5 + prod(log(x[1:2]))
}
 
smooth_simfunc <- function(x, jaatha) {
  stopifnot(length(x) == 2)
  idxs <- cbind(rep(1:10, each=12), rep(1:12, 10), x[1], x[2])
  sampled.values <- sapply(apply(idxs, 1, smooth_func), function(x) rpois(1, x))
  list(data=matrix(sampled.values, 10, 12, byrow=TRUE))
}

smooth_obs <- smooth_simfunc(c(3, 4))
smooth_stat <- jaatha:::Stat_PoiSmooth$new(smooth_obs$data, "(X1^2)*(X2^2)+log(X1)*log(X2)")
  
smooth_par_ranges <- matrix(c(2, 1, 7, 7), 2, 2)
rownames(smooth_par_ranges) <- c('x', 'y')
 
smooth_jaatha <- new("Jaatha", smooth_simfunc, smooth_par_ranges, list(mat=smooth_stat))
smooth_sim_data <- jaatha:::simulateWithinBlock(10, block.test, smooth_jaatha)
