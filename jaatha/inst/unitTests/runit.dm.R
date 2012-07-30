### --- Test setup ---
 
library("RUnit")
library("jaatha")

# dm.thetaTau <- dm.createDemographicModel(c(20,25),nLoci=100)
# dm.thetaTau <- dm.addSpeciationEvent(dm.thetaTau,0.1,5)
# dm.thetaTau <- dm.addMutation(dm.thetaTau,5,20)
# 
# dm.extTheta <- dm.createDemographicModel(c(20,25),nLoci=100)
# dm.extTheta <- dm.addSpeciationEvent(dm.extTheta,0.1,5)
# dm.extTheta <- dm.addSymmetricMigration(dm.extTheta,1,5)
# dm.extTheta <- dm.addMutation(dm.extTheta)
# 
# dm.eTp <- dm.createDemographicModel(c(20,25),nLoci=100)
# dm.eTp <- dm.addSpeciationEvent(dm.eTp,0.1,5)
# dm.eTp <- dm.addSymmetricMigration(dm.eTp,1,5)
# dm.eTp <- dm.addMutation(dm.eTp,1,5)
# 
# load("samples.save")

### --- Test functions ---


## -- Fixed bugs ----------------------------------------
#print() failed for empty demographic models
test.showEmptyModel <- function() {
  dm <- dm.createDemographicModel(25:26, 100)
  print(dm)
}
