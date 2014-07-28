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

dm <- dm.createThetaTauModel(12:13, 100)
jsfs <- dm.simSumStats(dm, c(1,5))
jaatha <- Jaatha.initialize(dm, jsfs=jsfs, smoothing=TRUE)

profile.is <- tempfile("jaatha_is")
gc()
Rprof(profile.is)
jaatha <- Jaatha.initialSearch(jaatha)
Rprof(NULL)
cat("Now run: R CMD Rprof", profile.is, "| less\n")

profile.rs <- tempfile("jaatha_is")
gc()
Rprof(profile.rs)
jaatha <- Jaatha.refinedSearch(jaatha, 2)
Rprof(NULL)
cat("Now run: R CMD Rprof", profile.rs, "| less\n")
