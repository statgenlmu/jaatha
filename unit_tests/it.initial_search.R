### --- Test setup ---

#rerecord.results = FALSE
rerecord.results = TRUE

library("RUnit")

test.dm.simSumStats <- function(){
  load("samples.save")
  checkException( dm.simSumStats( 1 ) )
  checkException( dm.simSumStats(dm.tt, 1 ) )
  checkException( dm.simSumStats(dm.tt, 1:3 ) )
  checkException( dm.simSumStats(dm.tt, c(2,50)) )
  set.seed(1)
  sim <- dm.simSumStats(dm.tt, c(2,10))[[1]]$jsfs
  if (rerecord.results) samples[["jsfs"]] <- sim
  checkEquals( sim , samples[["jsfs"]] )
  if (rerecord.results) save(samples,file="samples.save")
}

test.initialSearch.normal <- function(){
  load("samples.save")
  jaatha <- Jaatha.initialize(dm.tt, samples[["jsfs"]], seed=1)
  jaatha <- Jaatha.initialSearch(jaatha, sim=10, blocks.per.par=2)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  if (rerecord.results) samples[["initialSearch.normal"]] <- pStartPoints
  checkEquals(pStartPoints, samples[["initialSearch.normal"]])
  if (rerecord.results) save(samples, file="samples.save")

  # Test reproducibility
  jaatha2 <- Jaatha.initialize(dm.tt, samples[["jsfs"]], seed=1)
  jaatha2 <- Jaatha.initialSearch(jaatha2, sim=10, blocks.per.par=2)
  pStartPoints2 <- Jaatha.getStartingPoints(jaatha2)
  checkEquals(pStartPoints[1:2,], pStartPoints2[1:2,])

  # Parallel version
  jaatha <- Jaatha.initialize(dm.tt, samples[["jsfs"]], 
                              seed=1, cores=2)
  jaatha <- Jaatha.initialSearch(jaatha, sim=10, blocks.per.par=2)
  pStartPoints3 <- Jaatha.getStartingPoints(jaatha)
  checkEquals(pStartPoints, pStartPoints3)
}

test.initialSearch.seqgen <- function(){
  load("samples.save")
  set.seed(20)
  dm.sq <- dm.setMutationModel(dm.tt, "HKY", c(.2, .2, .2, .4), 0.5)
  jsfs <- dm.simSumStats(dm.sq, c(1.5, 7))
  jaatha <- Jaatha.initialize(dm.sq, jsfs=jsfs, seed=18)
  jaatha <- Jaatha.initialSearch(jaatha, sim=11, blocks.per.par=3)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  if (rerecord.results) samples[["initialSearch.seqgen"]] <- pStartPoints
  checkEquals(pStartPoints, samples[["initialSearch.seqgen"]])
  if (rerecord.results) save(samples,file="samples.save")
}

test.initialSearch.folded <- function() {
  load("samples.save")
  jsfs <- matrix(1, 21, 26)
  jaatha <- Jaatha.initialize(dm.tt, jsfs=jsfs, folded=T, seed=20)
  jaatha <- Jaatha.initialSearch(jaatha, sim=10, blocks.per.par=2)
  pStartPoints <- Jaatha.getStartingPoints(jaatha)
  if (rerecord.results) samples[["initialSearch.folded"]] <- pStartPoints
  checkEquals(pStartPoints, samples[["initialSearch.folded"]])
  if (rerecord.results) save(samples,file="samples.save")
}
