### --- Test setup ---

rerecord.results = FALSE
#rerecord.results = TRUE

library("RUnit")

dm.thetaTau <- dm.createThetaTauModel(c(10,15), 20, 100)


### --- Test functions ---

test.dm.simSumStats <- function(){
  load("samples.save")
  checkException( dm.simSumStats( 1 ) )
  checkException( dm.simSumStats(dm.thetaTau, 1 ) )
  checkException( dm.simSumStats(dm.thetaTau, 1:3 ) )
  checkException( dm.simSumStats(dm.thetaTau, c(2,50)) )
  set.seed(1)
  sim <- dm.simSumStats(dm.thetaTau, c(2,10))[[1]]$jsfs
  if (rerecord.results) samples[["jsfs"]] <- sim
  checkEquals( sim , samples[["jsfs"]] )
  if (rerecord.results) save(samples,file="samples.save")
}

test.initialSearch.normal <- function(){
  load("samples.save")
  jaatha <- Jaatha.initialize(dm.thetaTau, samples[["jsfs"]], seed=1,
                              sim.package.size=3)
  jaatha <- Jaatha.initialSearch(jaatha, sim=10, blocks.per.par=2)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  if (rerecord.results) samples[["refineSearch.jaatha"]] <- jaatha
  if (rerecord.results) samples[["initialSearch.normal"]] <- pStartPoints
  checkEquals(pStartPoints, samples[["initialSearch.normal"]])
  if (rerecord.results) save(samples, file="samples.save")

  # Test reproducibility
  jaatha2 <- Jaatha.initialize(dm.thetaTau, samples[["jsfs"]], seed=1,
                               sim.package.size=3)
  jaatha2 <- Jaatha.initialSearch(jaatha2, sim=10, blocks.per.par=2)
  pStartPoints2 <- Jaatha.getStartingPoints(jaatha2)
  checkEquals(pStartPoints[1:2,], pStartPoints2[1:2,])
}

test.initialSearch.seqgen <- function(){
  load("samples.save")
  set.seed(20)
  dm.sq <- dm.setMutationModel(dm.thetaTau, "HKY", c(.2, .2, .2, .4), 0.5)
  jsfs <- dm.simSumStats(dm.sq, c(1.5, 7))
  jaatha <- Jaatha.initialize(dm.sq, jsfs=jsfs, seed=18)
  jaatha <- Jaatha.initialSearch(jaatha, sim=10,blocks.per.par=2)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  if (rerecord.results) samples[["initialSearch.seqgen"]] <- pStartPoints
  checkEquals(pStartPoints, samples[["initialSearch.seqgen"]])
  if (rerecord.results) save(samples,file="samples.save")
}


test.initialSearch.folded <- function() {
  load("samples.save")
  jsfs <- matrix(1, 21, 26)
  jaatha <- Jaatha.initialize(dm.thetaTau, jsfs=jsfs, folded=T, seed=20)
  jaatha <- Jaatha.initialSearch(jaatha, sim=10, blocks.per.par=2)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  if (rerecord.results) samples[["initialSearch.folded"]] <- pStartPoints
  checkEquals(pStartPoints, samples[["initialSearch.folded"]])
  if (rerecord.results) save(samples,file="samples.save")
}

test.initialSearch.parallel <- function() {
  load("samples.save")
  jaatha <- Jaatha.initialize(dm.thetaTau, samples[["jsfs"]], seed=100,
                              sim.package.size=2)
  jaatha <- Jaatha.initialSearch(jaatha, sim=10, blocks.per.par=2)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  if (rerecord.results) samples[["initialSearch.parallel.simple"]] <- pStartPoints
  checkEquals(pStartPoints, samples[["initialSearch.parallel.simple"]])
  if (rerecord.results) save(samples,file="samples.save")
}

test.refinedSearch <- function(){
  load("samples.save")
  set.seed(1)
  jaatha <- samples[["refineSearch.jaatha"]]
  jaatha@cores <- 2
  jaatha@sim.package.size <- 5
  jaatha <- Jaatha.refinedSearch(jaatha, 2, sim=10,epsilon=.2,
                                 half.block.size=.05,weight=.9,max.steps=50)
  lt <- Jaatha.getLikelihoods(jaatha)
  if (rerecord.results) samples[["refineSearch.result"]] <- lt 
  checkEquals(lt, samples[["refineSearch.result"]])
  if (rerecord.results) save(samples,file="samples.save")
}
