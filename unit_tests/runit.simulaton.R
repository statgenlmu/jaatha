#!/usr/bin/Rscript --vanilla
#
# runit.simulaton.R
# Unit test for simulator.R
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-08-09
# Licence:  GPLv3 or later
#

library(jaatha)
library(RUnit)


test.createSimulationPackages <- function() {
  random.par <- matrix(runif(100*5),100,5)

  l1 <- createSimulationPackages(random.par, 10)
  checkEquals(length(l1), 10)
  checkEquals(l1[[10]], random.par[91:100, ])
  l1 <- createSimulationPackages(random.par,  7)
  checkEquals(length(l1), 15)
  l1 <- createSimulationPackages(random.par,  200)
  checkEquals(length(l1), 1)
}
