### --- Test setup ---

library("RUnit")
library("jaatha")

dm.thetaTau <- dm.createDemographicModel(c(15,20),nLoci=20)
dm.thetaTau <- dm.addSpeciationEvent(dm.thetaTau,0.1,5)
dm.thetaTau <- dm.addMutation(dm.thetaTau,5,20)

### --- Test functions ---

test.initialSearch.parallel.simple <- function(){
    load("samples.save")
	jaatha <- Jaatha.initialize(dm.thetaTau, samples[["simSumStats"]], seed=1, 
                                parallelization.model="simple",
                                sim.package.size=3)
	startPoints <- Jaatha.initialSearch(jaatha, nSim=10, nBlocksPerPar=2)
	pStartPoints <- Jaatha.printStartPoints(jaatha,startPoints)
	checkEquals(pStartPoints, samples[["initialSearch.normal"]])
}

# test.refineSearch <- function(){
#     load("samples.save")
#     set.seed(1)
#     startPoints <- Jaatha.pickBestStartPoints(samples[["refineSearch.startpoints"]],best=2)
#     jaatha <- Jaatha.initialize(dm.thetaTau, samples[["simSumStats"]], seed=1)
#     jaatha <- Jaatha.refineSearch(jaatha,startPoints,nSim=10,epsilon=.2,
#                       halfBlockSize=.05,weight=.9,nMaxStep=10)
#     lt <- Jaatha.printLikelihoods(jaatha)
#     if (rerecord.results) samples[["refineSearch.result"]] <- lt 
#     checkEquals(lt, samples[["refineSearch.result"]])
#     if (rerecord.results) save(samples,file="samples.save")
# }
