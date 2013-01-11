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

dm <- dm.createThetaTauModel(12:13, 50)
jsfs <- dm.simSumStats(dm, c(1,5))
jaatha <- Jaatha.initialize(dm, jsfs=jsfs, use.shm=TRUE)

profile.is<- "./tmp/jaatha_is"
Rprof(profile.is)
jaatha <- Jaatha.initialSearch(jaatha, sim=50)
Rprof(NULL)

profile.rs <- "./tmp/jaatha_rs"
Rprof(profile.rs)
jaatha <- Jaatha.refinedSearch(jaatha, 2, sim=50)
Rprof(NULL)

cat("Now run\n")
cat("R CMD Rprof", profile.is, "| less\n")
cat("R CMD Rprof", profile.rs, "| less\n")
