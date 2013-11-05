test.jaathaInitialiation <- function() {
  dm <- dm.createThetaTauModel(10:11, 100)
  jsfs <- matrix(rpois(11*12, 10), 11, 12)
  jaatha <- Jaatha.initialize(dm, jsfs)
  checkEquals( jaatha@nPar, 2 )

  dm <- dm.addSymmetricMigration(dm, 1, 3)
  jaatha <- Jaatha.initialize(dm, jsfs)
  checkEquals( jaatha@nPar, 3 )
}
