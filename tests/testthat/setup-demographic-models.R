# Theta-Tau Model 
dm.tt        <- dm.createThetaTauModel(11:12, 10)
sum.stats.tt <- dm.simSumStats(dm.tt, c(1, 5))
jaatha.tt    <- Jaatha.initialize(dm.tt, jsfs=sum.stats.tt, cores = 2) 

# Migration Model
dm.mig        <- dm.addSymmetricMigration(dm.tt, 1, 5)
sum.stats.mig <- dm.simSumStats(dm.mig, c(1, 1, 5))
jaatha.mig    <- Jaatha.initialize(dm.mig, jsfs=sum.stats.mig, cores = 2) 

# Groups
dm.grp <- dm.tt
dm.grp <- dm.setLociLength(dm.grp, 100, 1) 
dm.grp <- dm.setLociNumber(dm.grp, 15, 1) 
dm.grp <- dm.setLociLength(dm.grp, 200, 2) 
dm.grp <- dm.addSampleSize(dm.grp, 5:6, 3)
sum.stats.grp <- dm.simSumStats(dm.addSummaryStatistic(dm.grp, 'seg.sites'), 
                                                       c(1, 3))

# fpc model
dm.fpc <- dm.createDemographicModel(c(15,20), 100, 1000)
dm.fpc <- dm.addSpeciationEvent(dm.fpc, .1, 5)
dm.fpc <- dm.addRecombination(dm.fpc, 1, 5)
dm.fpc <- dm.addMutation(dm.fpc, 1, 10)
dm.fpc <- dm.addSymmetricMigration(dm.fpc, fixed=.75)
seg.sites <- dm.simSumStats(dm.addSummaryStatistic(dm.fpc, 'seg.sites'), 
                            c(1, 2, 3))$seg.sites
dm.fpc <- dm.addSummaryStatistic(dm.fpc, 'fpc')
dm.fpc <- jaatha:::calcFpcBreaks(dm.fpc, seg.sites)
sum.stats.fpc <- dm.simSumStats(dm.addSummaryStatistic(dm.fpc, 'seg.sites'), 
                                c(1, 2, 3))



# Finite Sites Models
if (jaatha:::checkForSeqgen(FALSE, TRUE)) {
  test_seqgen <- TRUE
  dm.sg <-  dm.addOutgroup(dm.tt, "2*tau")
  dm.hky <- dm.setMutationModel(dm.sg, "HKY", c(0.2, 0.2, 0.3, 0.3), 2)
  dm.hky <- dm.setLociNumber(dm.hky, 5)
  dm.hky <- dm.setLociLength(dm.hky, 15)
  dm.hky@sum.stats <- data.frame()
  dm.hky <- dm.addSummaryStatistic(dm.hky, 'jsfs')
  dm.f81 <- dm.setMutationModel(dm.sg, "F84", c(0.3, 0.2, 0.3, 0.2), 2)
  dm.gtr <- dm.setMutationModel(dm.sg, "GTR", 
                                gtr.rates=c(0.2, 0.2, 0.1, 0.1, 0.1, 0.2))
} else {
  test_seqgen <- FALSE
}

# Selection Models
if (jaatha:::checkForMsms(FALSE, TRUE)) {
  test_msms <- TRUE
  # Selection Model
  dm.sel <- jaatha:::dm.addPositiveSelection(dm.mig, 1000, 2000, population=2,
                                             at.time='tau/2')
  dm.sel <- dm.setLociNumber(dm.sel, 3)
} else {
  warning("Msms not found. Skipping tests.")
  test_msms <- FALSE
}

if (require('ape', quietly = TRUE)) {
  test_ape <- TRUE
} else {
  test_ape <- FALSE
  warning("Package ape not available. Skipping some tests.")
}