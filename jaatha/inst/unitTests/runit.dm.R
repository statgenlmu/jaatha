library("RUnit")
library("jaatha")

dm <- dm.createDemographicModel(c(20,25),nLoci=100)
dm <- dm.addSpeciationEvent(dm,0.1,5)
dm <- dm.addMutation(dm,5,20)

#------------------------------------------------------------------
# Tests 
#------------------------------------------------------------------

test.getFeature <- function() {
  checkEquals(nrow(getFeature(dm, "split")), 1)  
  checkEquals(nrow(getFeature(dm, "split", NA, NA, NA, NA)), 1)
  checkEquals(nrow(getFeature(dm, "mutation")), 1)
}

test.appendFeature <- function() {
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- appendFeature(dm, "split", 1, 1, 1)
  checkEquals(nrow(dm@features), 1)
  dm <- appendFeature(dm, "mut", 2, 1, 5, 1, 2, "t1", 1)
  checkEquals(nrow(dm@features), 2)
}

test.addParameter <- function() {
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- addParameter(dm, "theta")
  dm <- addParameter(dm, "rho")
  checkEquals(dm@parameters, c("theta","rho"))
  checkException(addParameter(dm, "theta"))
  checkEquals(dm@parameters, c("theta","rho"))
}

test.addFeature <- function() {
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- addFeature(dm, "split", "tau", 1, 10)
  checkEquals(nrow(dm@features), 1)
  dm <- addFeature(dm, "mutation", "theta", fixed.value=5,
                   pop.source=1, pop.sink=2,
                   time.point="t2", group=3)
  checkEquals(nrow(dm@features), 2)
}

test.makeThetaLast <- function() {
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- dm.addMutation(dm,5,20)
  dm <- dm.addSpeciationEvent(dm,0.1,5)
  checkEquals(dm@parameters[length(dm@parameters)], "theta")
}



## -- Fixed bugs ----------------------------------------
#print() failed for empty demographic models
test.showEmptyModel <- function() {
  dm <- dm.createDemographicModel(25:26, 100)
  print(dm)
}
