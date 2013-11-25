test.addParameter <- function() {
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- dm.addParameter(dm, "theta", 1, 5)
  dm <- dm.addParameter(dm, "rho", fixed=20)
  checkEquals(dm@parameters$name,  c("theta","rho"))
  checkEquals(dm@parameters$fixed, c(F,T))
  checkEquals(dm@parameters$lower.range, c(1,20))
  checkEquals(dm@parameters$upper.range, c(5,20))
}

test.addFeature <- function() {
  dm <- dm.createDemographicModel(11:12, 100)
  dm <- addFeature(dm, type="split", parameter="tau", 
                   lower.range=1, upper.range=10)
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
  checkEquals(dm.getParameters(dm)[dm.getNPar(dm)], "theta")
}

test.setMutationModel <- function() {
  dm <- dm.createThetaTauModel(11:12, 100)
  dm <- dm.setMutationModel(dm, "HKY")
  checkTrue("mutation.model" %in% dm@parameters$name)
  checkException(dm <- dm.setMutationModel(dm, "bla"))
}

test.simSumStats <- function() {
  checkException( dm.simSumStats( 1 ) )
  checkException( dm.simSumStats(dm.tt, 1 ) )
  checkException( dm.simSumStats(dm.tt, 1:3 ) )
  checkException( dm.simSumStats(dm.tt, c(2,50)) )

  sum.stats <- dm.simSumStats(dm.tt, c(1,5), "jsfs")
  checkTrue( is.list(sum.stats) )
  checkTrue( !is.null(sum.stats$jsfs) )
  checkTrue( sum(sum.stats$jsfs) > 0 )
}

test.parInRange <- function() {
  checkParInRange(dm.tt, c(1, 5))
  checkParInRange(dm.tt, c(2, 7))
  checkParInRange(dm.tt, c(2.1, 7.7))
  checkException( checkParInRange(dm.tt, c(0, 5)) )
  checkException( checkParInRange(dm.tt, c(0, -1)) )
  checkException( checkParInRange(dm.tt, c(10, 1)) )
  checkException( checkParInRange(dm.tt, c(100, 100)) )
  checkException( checkParInRange(dm.tt, 1) )
  checkException( checkParInRange(dm.tt, matrix(1, 2, 2) ))
  checkException( checkParInRange(dm.tt, NULL ))
}

## -- Fixed bugs ----------------------------------------
test.showEmptyModel <- function() {
  dm <- dm.createDemographicModel(25:26, 100)
  invisible(print(dm))
}
