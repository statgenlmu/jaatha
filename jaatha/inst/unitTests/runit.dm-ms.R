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

test.callMS <- function() {
  checkException(callMs())
  set.seed(14)
  checkEquals(callMs(2,1,"-t 5")[4], "positions:    0.6319236931 ")
}
