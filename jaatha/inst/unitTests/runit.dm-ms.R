#!/usr/bin/Rscript --vanilla
#
# DemographicModel-ms.test
# Unit test for the ms-adapter to the demographic model class
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-07-27
# Licence:  GPLv3 or later
#

# - Test setup -----------------------------------
 
library("RUnit")
library("jaatha")

# ------------------------------------------------

dm <- dm.createThetaTauModel(c(10,25), 100)

test.callMS <- function() {
  checkException(callMs())
  checkTrue(file.exists(callMs("2 1 -t 5")))
}

test.msSingleSimFunc <- function(){
  set.seed(100)
  jsfs <- msSingleSimFunc(dm, c(1,10))
  print(jsfs[1, 2])
  checkEquals(jsfs[1, 2], 1058)
}

test.simProg <- function(){
  checkTrue(!is.null(.local$simProgs[["ms"]]))
}
