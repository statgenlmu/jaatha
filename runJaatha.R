library("jaatha")

dm <- dm.createDemographicModel(sampleSizes=c(25,25),nLoci=100,seqLength=1000,debug=T,logFile="debug.log")
dm <- dm.addSpeciationEvent(dm,0.017,5)
dm <- dm.addMutation(dm,5,20)
dm <- dm.addPresentSize(dm,1,10,population=2)
dm <- dm.addSplitSize(dm,1,5,population=2)
dm <- dm.addRecombination(dm,fixed=20)
dm <- dm.addSymmetricMigration(dm,fixed=.5)

sumStats <- dm.simSumStats(dm,c(1,4,2,10))

jaatha <- Jaatha.initialize(dm, sumStats,debug=T,logFile="debug.log")
startPoints <- Jaatha.initialSearch(jaatha,nSim=50,nBlocksPerPar=2,multiple=T)
startPoints <- Jaatha.pickBestStartPoints(startPoints,best=1)
jaatha <- Jaatha.refineSearch(jaatha,startPoints,nSim=50,epsilon=2,
			             halfBlockSize=.05,weight=.9,nMaxStep=100)
