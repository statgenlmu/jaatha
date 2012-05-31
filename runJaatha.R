library("jaatha")

dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000,debug=T,logFile="debug.log")
dm <- dm.addSpeciationEvent(dm,0.017,5)
dm <- dm.addMutation(dm,5,20)
dm <- dm.addPresentSize(dm,1,10,population=2)
dm <- dm.addSplitSize(dm,1,5,population=2)
dm <- dm.addRecombination(dm,fixed=20)
dm <- dm.addSymmetricMigration(dm,fixed=.5)

sumStats <- dm.simSumStats(dm,c(1,4,2,10))
dm2 <- dm 
rm(dm)

dm.extTheta <- dm.createDemographicModel(c(20,25),nLoci=100,debug=T, logFile="debug.log")
dm.extTheta <- dm.addSpeciationEvent(dm.extTheta,0.1,5)
dm.extTheta <- dm.addSymmetricMigration(dm.extTheta,1,5)
dm.extTheta <- dm.addMutation(dm.extTheta,1,5)
sumStats <- dm.simSumStats(dm.extTheta,c(1,2,10))
dm2 <- dm.extTheta

jaatha <- Jaatha.initialize(dm2, sumStats, debug=T, logFile="debug.log")
startPoints <- Jaatha.initialSearch(jaatha,nSim=50,nBlocksPerPar=2)
startPoints <- Jaatha.pickBestStartPoints(startPoints,best=1)
jaatha <- Jaatha.refineSearch(jaatha,startPoints,nSim=50,epsilon=2,
			             halfBlockSize=.05,weight=.9,nMaxStep=100)
