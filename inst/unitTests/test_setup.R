# --------------------------------------------------------------
# Authors:  Paul R. Staab
# Date:     2013-11-13
# Licence:  GPLv3 or later
# --------------------------------------------------------------

# Setup for tests
dm.tt        <- dm.createThetaTauModel(11:12, 10)
sum.stats.tt <- dm.simSumStats(dm.tt, c(1, 5))
jaatha.tt    <- Jaatha.initialize(dm.tt, jsfs=sum.stats.tt) 

dm.mig        <- dm.addSymmetricMigration(dm.tt, 1, 5)
sum.stats.mig <- dm.simSumStats(dm.mig, c(1, 1, 5))
jaatha.mig    <- Jaatha.initialize(dm.mig, jsfs=sum.stats.mig) 

csi.sim.func <- function(x, jaatha) {
  list(poisson.vector=c(rpois(3, x[1]), rpois(3, x[2])))
}
csi.obs <- csi.sim.func(c(2:3))
csi.sum.stats <- list("poisson.vector"=list(method="poisson.independent",
                                            value=csi.obs$poisson.vector))
csi.par.ranges <- matrix(c(0.1, 0.1, 10, 10), 2, 2)
rownames(csi.par.ranges) <- c('x', 'y')

jaatha.csi <- new("Jaatha", csi.sim.func, csi.par.ranges, csi.sum.stats, 123)
jaatha.csi <- Jaatha.initialSearch(jaatha.csi, 20, 2)
jaatha.csi <- Jaatha.refinedSearch(jaatha.csi, 1, 20, max.steps=10)


block.test <- new("Block")
block.test@border <- matrix(c(0, 0, 0.1, 0.1), 2, 2) 
