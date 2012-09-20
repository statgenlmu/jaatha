### --- Test setup ---

rerecord.results = F
 
library("RUnit")

dm.thetaTau <- dm.createDemographicModel(c(20,25), 20)
dm.thetaTau <- dm.addSpeciationEvent(dm.thetaTau,0.1,5)
dm.thetaTau <- dm.addMutation(dm.thetaTau,5,20)

dm.extTheta <- dm.createDemographicModel(c(25,15), 80)
dm.extTheta <- dm.addSpeciationEvent(dm.extTheta,0.1,5)
dm.extTheta <- dm.addSymmetricMigration(dm.extTheta,1,5)
dm.extTheta <- dm.addMutation(dm.extTheta)

dm.eTp <- dm.createDemographicModel(c(20,40), 50)
dm.eTp <- dm.addSpeciationEvent(dm.eTp,0.1,5)
dm.eTp <- dm.addSymmetricMigration(dm.eTp,1,5)
dm.eTp <- dm.addMutation(dm.eTp,1,5)


### --- Test functions ---

test.dm.simSumStats <- function(){
    load("samples.save")
	checkException( dm.simSumStats( 1 ) )
	checkException( dm.simSumStats(dm.thetaTau, 1 ) )
	checkException( dm.simSumStats(dm.thetaTau, 1:3 ) )
	checkException( dm.simSumStats(dm.thetaTau, c(2,50)) )
	set.seed(1)
	sim <- dm.simSumStats(dm.thetaTau, c(2,10), Jaatha.defaultSumStats )
	if (rerecord.results) samples[["simSumStats"]] <- sim
	checkEquals( sim , samples[["simSumStats"]] )
    if (rerecord.results) save(samples,file="samples.save")
}

test.initialSearch.normal <- function(){
    load("samples.save")
	jaatha <- Jaatha.initialize(dm.thetaTau, samples[["simSumStats"]], seed=1)
	startPoints <- Jaatha.initialSearch(jaatha,nSim=10,nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha,startPoints)
	if (rerecord.results) samples[["initialSearch.normal"]] <- pStartPoints
	checkEquals(pStartPoints, samples[["initialSearch.normal"]])
    if (rerecord.results) save(samples,file="samples.save")
}

test.initialSearch.extTheta <- function(){
    load("samples.save")
	jaatha <- Jaatha.initialize(dm.extTheta, samples[["sumStats.extTheta"]], seed=1)
	startPoints <- Jaatha.initialSearch(jaatha,nSim=10,nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha,startPoints)
	if (rerecord.results) samples[["initialSearch.extTheta"]] <- pStartPoints
	checkEquals(pStartPoints, samples[["initialSearch.extTheta"]])
    if (rerecord.results) save(samples,file="samples.save")
}

test.initialSearch.extThetaPossible <- function(){
    load("samples.save")
	jaatha <- Jaatha.initialize(dm.eTp, samples[["sumStats.extTheta"]], seed=1)
	startPoints <- Jaatha.initialSearch(jaatha,nSim=10,nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha, startPoints)
	if (rerecord.results) samples[["initialSearch.eTp"]] <- pStartPoints
	checkEquals(pStartPoints, samples[["initialSearch.eTp"]])
    if (rerecord.results) save(samples,file="samples.save")
}

test.initialSearch.folded <- function() {
    load("samples.save")
    jsfs <- matrix(1, 21, 26)
	jaatha <- Jaatha.initialize(dm.thetaTau, jsfs=jsfs, folded=T, seed=1)
	startPoints <- Jaatha.initialSearch(jaatha, nSim=10, nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha, startPoints)
	if (rerecord.results) samples[["initialSearch.folded"]] <- pStartPoints
	checkEquals(pStartPoints, samples[["initialSearch.folded"]])
    if (rerecord.results) save(samples,file="samples.save")
}

test.refineSearch <- function(){
    load("samples.save")
	set.seed(1)
	startPoints <- Jaatha.pickBestStartPoints(samples[["refineSearch.startpoints"]],best=2)
	jaatha <- Jaatha.initialize(dm.thetaTau, samples[["simSumStats"]], seed=1)
	jaatha <- Jaatha.refineSearch(jaatha,startPoints,nSim=10,epsilon=.2,
				      halfBlockSize=.05,weight=.9,nMaxStep=10)
	lt <- Jaatha.printLikelihoods(jaatha)
	if (rerecord.results) samples[["refineSearch.result"]] <- lt 
	checkEquals(lt, samples[["refineSearch.result"]])
    if (rerecord.results) save(samples,file="samples.save")
}

## -- Fixed bugs ----------------------------------------
#print() failed for empty demographic models
test.showEmptyModel <- function() {
  dm <- dm.createDemographicModel(25:26, 100)
  print(dm)
}

