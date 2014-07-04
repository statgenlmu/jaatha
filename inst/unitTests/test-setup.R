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

# A Block
block.test <- new("Block")
block.test@border <- matrix(c(0, 0, 0.1, 0.1), 2, 2) 


# Theta-Tau Model 
dm.tt        <- dm.createThetaTauModel(11:12, 10)
sum.stats.tt <- dm.simSumStats(dm.tt, c(1, 5))
jaatha.tt    <- Jaatha.initialize(dm.tt, jsfs=sum.stats.tt) 
sim.data.tt  <- simulateWithinBlock(10, block.test, jaatha.tt)


# Migration Model
dm.mig        <- dm.addSymmetricMigration(dm.tt, 1, 5)
sum.stats.mig <- dm.simSumStats(dm.mig, c(1, 1, 5))
jaatha.mig    <- Jaatha.initialize(dm.mig, jsfs=sum.stats.mig) 


# Groups
dm.grp <- dm.tt
dm.grp <- dm.setLociLength(dm.grp, 100, 1) 
dm.grp <- dm.setLociNumber(dm.grp, 15, 1) 
dm.grp <- dm.setLociLength(dm.grp, 200, 2) 
dm.grp <- dm.addSampleSize(dm.grp, 5:6, 3)
sum.stats.grp <- dm.simSumStats(dm.grp, c(1, 5))


# Finite Sites Models
dm.sg <-  dm.addOutgroup(dm.tt, "2*tau")
dm.hky <- dm.setMutationModel(dm.sg, "HKY", c(0.2, 0.2, 0.3, 0.3), 2)
dm.hky <- dm.setLociNumber(dm.hky, 5)
dm.hky <- dm.setLociLength(dm.hky, 15)
dm.hky@sum.stats <- data.frame()
dm.hky <- dm.addSummaryStatistic(dm.hky, 'jsfs')
dm.hky <- dm.addSummaryStatistic(dm.hky, 'file')
dm.f81 <- dm.setMutationModel(dm.sg, "F84", c(0.3, 0.2, 0.3, 0.2), 2)
dm.gtr <- dm.setMutationModel(dm.sg, "GTR", gtr.rates=c(0.2, 0.2, 0.1, 0.1, 0.1, 0.2))
set.seed(20)
sg.file <- dm.simSumStats(dm.hky, c(1, 5))$file['seqgen']
sg.example <- scan(sg.file, 'char', sep="\n")
unlink(sg.file)
rm(sg.file)

# fpc model
dm.fpc <- dm.createDemographicModel(c(15,20), 100, 1000)
dm.fpc <- dm.addSpeciationEvent(dm.fpc, .1, 5)
dm.fpc <- dm.addRecombination(dm.fpc, 1, 5)
dm.fpc <- dm.addMutation(dm.fpc, 1, 10)
dm.fpc <- dm.addSymmetricMigration(dm.fpc, fixed=.75)
seg.sites <- dm.simSumStats(dm.addSummaryStatistic(dm.fpc, 'seg.sites'), c(1, 2, 3))$seg.sites
dm.fpc <- dm.addSummaryStatistic(dm.fpc, 'fpc')
dm.fpc <- calcFpcBreaks(dm.fpc, seg.sites)
sum.stats.fpc <- dm.simSumStats(dm.addSummaryStatistic(dm.fpc, 'seg.sites'), c(1, 2, 3))

# fpc + seq-gen
seg.sites <- dm.simSumStats(dm.addSummaryStatistic(dm.hky, 'seg.sites'), c(1, 5))$seg.sites
dm.sgfpc <- dm.addSummaryStatistic(dm.hky, 'fpc')
dm.sgfpc <- calcFpcBreaks(dm.sgfpc, seg.sites, 3)

# Selection Model
dm.sel <- dm.addPositiveSelection(dm.mig, 1000, 2000, population=2,
                                  at.time='tau/2')
dm.sel <- dm.setLociNumber(dm.sel, 3)
dm.sel <- calcFpcBreaks(dm.sel, seg.sites)

# Custom Simulation Interface
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

sim.data.csi <- simulateWithinBlock(10, block.test, jaatha.csi)


# Smoothing
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

smooth.jaatha <- new("Jaatha", smooth.simfunc, smooth.par.ranges, smooth.sum.stats, 123)
smooth.border.jaatha <- new("Jaatha", smooth.simfunc, smooth.par.ranges, 
                            smooth.border.sum.stats, 123)
smooth.sim.data <- simulateWithinBlock(10, block.test, smooth.jaatha)

save(list=ls(), file="test-setup.Rda")
