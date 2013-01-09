# unit_tests/runit.fixed-bugs
# Unit test ensuring regressions of already fixed bugs
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-11-30
# Licence:  GPLv3 or later
#

## -- Fixed bugs ----------------------------------------

#print() failed for empty demographic models
test.showEmptyModel <- function() {
  dm <- dm.createDemographicModel(25:26, 100)
  print(dm)
}

#Jaatha was not able to run on one parameter models
test.oneParModel <- function(){
  dm.onePar <- dm.createDemographicModel(10:11, 20)
  dm.onePar <- dm.addSpeciationEvent(dm.onePar, fixed.time=0)
  dm.onePar <- dm.addMutation(dm.onePar, 1, 20)
  jsfs <- matrix(rpois(11*12, 10), 11, 12)
  jaatha <- Jaatha.initialize(dm.onePar, jsfs=jsfs, seed=100)
  jaatha <- Jaatha.initialSearch(jaatha, sim=10, blocks.per.par=2)
  checkEquals(nrow(Jaatha.getStartingPoints(jaatha)), 2)
}
