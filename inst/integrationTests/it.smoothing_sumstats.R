test.smoothing <- function() {
  smooth.jaatha <- Jaatha.initialSearch(smooth.jaatha, 25, 2)
  pStartPoints <- Jaatha.getStartingPoints(smooth.jaatha)
  checkEquals(4, nrow(pStartPoints))

  smooth.jaatha <- Jaatha.refinedSearch(smooth.jaatha, 1, 25, 25)
  checkTrue( ncol(smooth.jaatha@likelihood.table) == 4 )
  checkTrue( nrow(smooth.jaatha@likelihood.table) >= 5 )
}
