### --- Test setup ---
 
library("RUnit")
library("jaatha")

dm.thetaTau <- dm.createDemographicModel(c(20,25),nLoci=100)
dm.thetaTau <- dm.addSpeciationEvent(dm.thetaTau,0.1,5)
dm.thetaTau <- dm.addMutation(dm.thetaTau,5,20)

dm.extTheta <- dm.createDemographicModel(c(20,25),nLoci=100)
dm.extTheta <- dm.addSpeciationEvent(dm.extTheta,0.1,5)
dm.extTheta <- dm.addSymmetricMigration(dm.extTheta,1,5)
dm.extTheta <- dm.addMutation(dm.extTheta)

dm.eTp <- dm.createDemographicModel(c(20,25),nLoci=100)
dm.eTp <- dm.addSpeciationEvent(dm.eTp,0.1,5)
dm.eTp <- dm.addSymmetricMigration(dm.eTp,1,5)
dm.eTp <- dm.addMutation(dm.eTp,1,5)

load("samples.save")

### --- Test functions ---

test.dm.simSumStats <- function(){
	checkException( dm.simSumStats( 1 ) )
	checkException( dm.simSumStats(dm.thetaTau, 1 ) )
	checkException( dm.simSumStats(dm.thetaTau, 1:3 ) )
	checkException( dm.simSumStats(dm.thetaTau, c(2,50)) )
	set.seed(1)
	sim <- dm.simSumStats(dm.thetaTau, c(2,10))
	#samples[["simSumStats"]] <- sim
	#save(samples,file="samples.save")
	checkEquals( sim , samples[["simSumStats"]] )
}

test.initialSearch.normal <- function(){
	jaatha <- Jaatha.initialize(dm.thetaTau, samples[["simSumStats"]], seed=1)
	startPoints <- Jaatha.initialSearch(jaatha,nSim=10,nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha,startPoints)
	#samples[["initialSearch.normal"]] <- pStartPoints
	#save(samples,file="samples.save")
	checkEquals(pStartPoints, samples[["initialSearch.normal"]])
}

test.initialSearch.extTheta <- function(){
	jaatha <- Jaatha.initialize(dm.extTheta, samples[["sumStats.extTheta"]], seed=1)
	startPoints <- Jaatha.initialSearch(jaatha,nSim=10,nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha,startPoints)
	#samples[["initialSearch.extTheta"]] <- pStartPoints
	#save(samples,file="samples.save")
	checkEquals(pStartPoints, samples[["initialSearch.extTheta"]])
}

test.initialSearch.extThetaPossible <- function(){
	jaatha <- Jaatha.initialize(dm.eTp, samples[["sumStats.extTheta"]], seed=1)
	startPoints <- Jaatha.initialSearch(jaatha,nSim=10,nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha,startPoints)
	#samples[["initialSearch.eTp"]] <- pStartPoints
	#save(samples,file="samples.save")
	checkEquals(pStartPoints, samples[["initialSearch.eTp"]])
}

test.refineSearch <- function(){
	set.seed(1)
	startPoints <- Jaatha.pickBestStartPoints(samples[["refineSearch.startpoints"]],best=2)
	jaatha <- Jaatha.initialize(dm.thetaTau, samples[["simSumStats"]], seed=1)
	jaatha <- Jaatha.refineSearch(jaatha,startPoints,nSim=10,epsilon=.2,
				      halfBlockSize=.05,weight=.9,nMaxStep=10)
	lt <- Jaatha.printLikelihoods(jaatha)
	#samples[["refineSearch.result"]] <- lt; save(samples,file="samples.save")
	checkEquals(lt, samples[["refineSearch.result"]])
}

## -- Fixed bugs ----------------------------------------
#print() failed for empty demographic models
test.showEmptyModel <- function() {
  dm <- dm.createDemographicModel(25:26, 100)
  print(dm)
}
