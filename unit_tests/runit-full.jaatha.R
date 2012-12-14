### --- Test setup ---

rerecord.results = F
 
library("RUnit")

dm.thetaTau <- dm.createDemographicModel(c(10,15), 20, 100)
dm.thetaTau <- dm.addSpeciationEvent(dm.thetaTau,0.1,5)
dm.thetaTau <- dm.addMutation(dm.thetaTau,5,20)

dm.extTheta <- dm.createDemographicModel(c(12,11), 30)
dm.extTheta <- dm.addSpeciationEvent(dm.extTheta,0.1,5)
dm.extTheta <- dm.addSymmetricMigration(dm.extTheta,1,5)
dm.extTheta <- dm.addMutation(dm.extTheta)

dm.eTp <- dm.createDemographicModel(c(12,13), 25)
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
	jaatha <- Jaatha.initialize(dm.thetaTau, samples[["simSumStats"]], seed=1,
                                sim.package.size=3)
	startPoints <- Jaatha.initialSearch(jaatha, nSim=10, nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha,startPoints)
	if (rerecord.results) samples[["initialSearch.normal"]] <- pStartPoints
	checkEquals(pStartPoints, samples[["initialSearch.normal"]])
    if (rerecord.results) save(samples,file="samples.save")

    # Test reproducibility
	jaatha2 <- Jaatha.initialize(dm.thetaTau, samples[["simSumStats"]], seed=1,
                                sim.package.size=3)
	startPoints2 <- Jaatha.initialSearch(jaatha2, nSim=10, nBlocksPerPar=2)
	pStartPoints2 <- Jaatha.printStartPoints(jaatha, startPoints2)
	checkEquals(pStartPoints[1:2,], pStartPoints2[1:2,])
}

test.initialSearch.extThet <- function(){
    load("samples.save")
	jaatha <- Jaatha.initialize(dm.extTheta, samples[["sumStats.extTheta"]],
                                seed=5)
	startPoints <- Jaatha.initialSearch(jaatha,nSim=10,nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha,startPoints)
	if (rerecord.results) samples[["initialSearch.extTheta"]] <- pStartPoints
	checkEquals(pStartPoints, samples[["initialSearch.extTheta"]])
    if (rerecord.results) save(samples,file="samples.save")
}

test.initialSearch.seqgen <- function(){
    load("samples.save")
    dm.sq <- dm.setMutationModel(dm.thetaTau, "HKY", c(.2, .2, .2, .4), 0.5)
	jaatha <- Jaatha.initialize(dm.sq, samples[["sumStats.extTheta"]],
                                seed=10)
	startPoints <- Jaatha.initialSearch(jaatha,nSim=10,nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha,startPoints)
	if (rerecord.results) samples[["initialSearch.seqgen"]] <- pStartPoints
	checkEquals(pStartPoints, samples[["initialSearch.seqgen"]])
    if (rerecord.results) save(samples,file="samples.save")
}

test.initialSearch.extThetaPossible <- function(){
    load("samples.save")
	jaatha <- Jaatha.initialize(dm.eTp, samples[["sumStats.extTheta"]], seed=10)
	startPoints <- Jaatha.initialSearch(jaatha,nSim=10,nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha, startPoints)
	if (rerecord.results) samples[["initialSearch.eTp"]] <- pStartPoints
	checkEquals(pStartPoints, samples[["initialSearch.eTp"]])
    if (rerecord.results) save(samples,file="samples.save")
}

test.initialSearch.folded <- function() {
    load("samples.save")
    jsfs <- matrix(1, 21, 26)
	jaatha <- Jaatha.initialize(dm.thetaTau, jsfs=jsfs, folded=T, seed=20)
	startPoints <- Jaatha.initialSearch(jaatha, nSim=10, nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha, startPoints)
	if (rerecord.results) samples[["initialSearch.folded"]] <- pStartPoints
	checkEquals(pStartPoints, samples[["initialSearch.folded"]])
    if (rerecord.results) save(samples,file="samples.save")
}

test.initialSearch.parallel <- function() {
    load("samples.save")
	jaatha <- Jaatha.initialize(dm.thetaTau, samples[["simSumStats"]], seed=100,
                                sim.package.size=2)
	startPoints <- Jaatha.initialSearch(jaatha, nSim=10, nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha,startPoints)
	if (rerecord.results) samples[["initialSearch.parallel.simple"]] <- pStartPoints
	checkEquals(pStartPoints, samples[["initialSearch.parallel.simple"]])
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


