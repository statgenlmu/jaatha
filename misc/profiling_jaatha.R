#!/usr/bin/Rscript --vanilla
#
# profiling_jaatha
# %DESCRIPTION%
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-12-07
# Licence:  GPLv3 or later
#

library(jaatha)

dm <- dm.createThetaTauModel(c(20,25), 20000)
sim_data <- dm.simSumStats(dm.addSummaryStatistic(dm, 'seg.sites'), c(.1,5))
jaatha <- Jaatha.initialize(jsfs, dm, use_fpc=TRUE, smoothing=FALSE)

profile.is <- "./jaatha_is2.prof"
gc()
Rprof(profile.is, memory.profiling=TRUE)
jaatha <- Jaatha.initialSearch(jaatha)
Rprof(NULL)
cat("Now run: R CMD Rprof", profile.is, "| less\n")

profile.rs <- "./jaatha_rs2.prof"
gc()
Rprof(profile.rs, memory.profiling=TRUE)
jaatha <- Jaatha.refinedSearch(jaatha, 2)
Rprof(NULL)
cat("Now run: R CMD Rprof", profile.rs, "| less\n")
