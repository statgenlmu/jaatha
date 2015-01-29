set.seed(10121416)

# Theta-Tau Model 
dm.tt <- dm.createDemographicModel(11:12, 5, 100) +
  feat_pop_merge(par_range('tau', .1, 2), 2, 1) +
  feat_recombination(par_const(.5)) +
  feat_mutation(par_range('theta', 1, 10))

sum.stats.tt <- dm.simSumStats(dm.addSummaryStatistic(dm.tt, 'seg.sites'), 
                               c(1, 5))
jaatha.tt    <- Jaatha.initialize(sum.stats.tt, dm.tt, cores = 2) 

# Migration Model
dm.mig <- dm.createDemographicModel(11:12, 5) +
  feat_migration(par_range('M', .1, 5), symmetric = TRUE) +
  feat_pop_merge(par_range('tau', .1, 2), 2, 1) +
  feat_recombination(par_const(.5)) +
  feat_mutation(par_range('theta', 1, 10))

sum.stats.mig <- dm.simSumStats(dm.addSummaryStatistic(dm.mig, 'seg.sites'), 
                                c(.3, 1, 5))
jaatha.mig    <- Jaatha.initialize(sum.stats.mig, dm.mig, cores = 2)

# Groups
dm.grp <- dm.mig
dm.grp <- dm.addLocus(dm.grp, 100, 15, 1)
dm.grp <- dm.addLocus(dm.grp, 200, 10, 2)
dm.grp <- dm.addLocus(dm.grp, 50, 5, 3)
sum.stats.grp <- dm.simSumStats(dm.addSummaryStatistic(dm.grp, 'seg.sites'), 
                                c(1, 1, 3))


# Finite Sites Models
if (jaatha:::checkForSeqgen(FALSE, TRUE)) {
  test_seqgen <- TRUE
  dm.sg <-  dm.addOutgroup(dm.tt, "2*tau")
  dm.hky <- dm.setMutationModel(dm.sg, "HKY", c(0.2, 0.2, 0.3, 0.3), 2)
  dm.hky <- jaatha:::dm.setLociNumber(dm.hky, 5)
  dm.hky <- jaatha:::dm.setLociLength(dm.hky, 15)
  dm.hky <- jaatha:::resetSumStats(dm.hky)
  dm.hky <- dm.addSummaryStatistic(dm.hky, 'jsfs')
  dm.f81 <- dm.setMutationModel(dm.sg, "F84", c(0.3, 0.2, 0.3, 0.2), 2)
  dm.gtr <- dm.setMutationModel(dm.sg, "GTR", 
                                gtr.rates=c(0.2, 0.2, 0.1, 0.1, 0.1, 0.2))
  
  dm_trios <- dm.addLocusTrio(dm.hky, locus_length = rep(1000, 3), 
                              distance = c(7502, 9050), group = 2)
  dm_trios <- dm.addLocusTrio(dm_trios, locus_length = rep(950, 3), 
                              distance = c(17502, 10050), group = 2)
  dm_trios <- dm.addLocusTrio(dm_trios, locus_length = rep(1050, 3), 
                              distance = c(7502, 15050), group = 2)
  dm_trios <- dm.addLocusTrio(dm_trios, locus_length = rep(1060, 3), 
                              distance = c(502, 15050), group = 2)
  dm_trios <- dm.addLocusTrio(dm_trios, locus_length = rep(600, 3), 
                              distance = c(6502, 35050), group = 2)
  trios_sum_stats <- dm.simSumStats(dm.addSummaryStatistic(dm_trios, 'seg.sites'), 
                                    c(1, 5))
} else {
  test_seqgen <- FALSE
}

# Selection Models
if (jaatha:::checkForMsms(FALSE, TRUE)) {
  test_msms <- TRUE
  # Selection Model
  dm.sel <- jaatha:::dm.addPositiveSelection(dm.mig, 1000, 2000, population=2,
                                             at.time='tau/2')
  dm.sel <- jaatha:::dm.setLociNumber(dm.sel, 3)
} else {
  test_msms <- FALSE
}

if (require('ape', quietly = TRUE)) {
  test_ape <- TRUE
} else {
  test_ape <- FALSE
}



