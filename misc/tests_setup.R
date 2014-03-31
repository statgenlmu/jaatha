#!/usr/bin/Rscript --vanilla
#
# tests_setup
# %DESCRIPTION%
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2013-09-04
# Licence:  GPLv3 or later
#

dm <- dm.createThetaTauModel(11:12, 20)
jsfs <- dm.simSumStats(dm, c(1,5))
jaatha <- Jaatha.initialize(dm, jsfs=jsfs, cores=2, use.shm=TRUE) 
