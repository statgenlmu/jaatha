#!/usr/bin/Rscript
#
# --------------------------------------------------------------
# Authors:  Paul R. Staab
# Date:     2013-11-13
# Licence:  GPLv3 or later
# --------------------------------------------------------------

library(jaatha)

if (is.element("jaatha", loadedNamespaces()))
  attach(loadNamespace("jaatha"), name=paste("namespace", "jaatha", sep=":"), pos=3)

set.seed(13579)

# Setup for tests
dm.tt        <- dm.createThetaTauModel(11:12, 10)
sum.stats.tt <- dm.simSumStats(dm.tt, c(1, 5))
jaatha.tt    <- Jaatha.initialize(dm.tt, jsfs=sum.stats.tt) 

dm.mig        <- dm.addSymmetricMigration(dm.tt, 1, 5)
sum.stats.mig <- dm.simSumStats(dm.mig, c(1, 1, 5))
jaatha.mig    <- Jaatha.initialize(dm.mig, jsfs=sum.stats.mig) 

dm.sg <- dm.addOutgroup(dm.tt, "2*tau")
dm.hky <- dm.setMutationModel(dm.sg, "HKY", c(0.2, 0.2, 0.3, 0.3), 2)
dm.f81 <- dm.setMutationModel(dm.sg, "F84", c(0.3, 0.2, 0.3, 0.2), 2)
dm.gtr <- dm.setMutationModel(dm.sg, "GTR", gtr.rates=c(0.2, 0.2, 0.1, 0.1, 0.1, 0.2))

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

sim.data.csi <- simulateWithinBlock(10, block.test, jaatha.csi)
sim.data.tt  <- simulateWithinBlock(10, block.test, jaatha.tt)

save(list=ls(), file="test_setup.Rda")
