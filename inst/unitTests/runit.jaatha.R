test.jaathaInitialiation <- function() {
}

test.getParNames <- function() {
  checkEquals( getParNames(jaatha.tt), dm.getParameters(dm.tt) )
  checkEquals( getParNames(jaatha.mig), dm.getParameters(dm.mig) )
}

test.getParNumber <- function() {
  checkEquals( getParNumber(jaatha.tt), dm.getNPar(dm.tt) )
  checkEquals( getParNumber(jaatha.mig), dm.getNPar(dm.mig) )
}
