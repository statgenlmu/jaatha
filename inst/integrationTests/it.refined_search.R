test.refinedSearch.csi <- function() {
  jaatha <- Jaatha.refinedSearch(jaatha.csi, 2, sim=10, sim.final=10)
  checkTrue( is.matrix(jaatha@likelihood.table) )
  checkTrue( ncol(jaatha@likelihood.table) == 4 )
  checkTrue( nrow(jaatha@likelihood.table) >= 10 )

  # Test reproducibility
  jaatha2 <- Jaatha.refinedSearch(jaatha.csi, 2, sim=10, sim.final=10)

  checkEquals(jaatha@likelihood.table, jaatha2@likelihood.table)
}

test.refinedSearch.dm <- function() {
  jaatha <- Jaatha.initialSearch(jaatha.tt, 20, 1)
  jaatha <- Jaatha.refinedSearch(jaatha, 1, 10, sim.final=10, max.steps=10)
  checkTrue( is.matrix(jaatha@likelihood.table) )
  checkTrue( ncol(jaatha@likelihood.table) == 4 )
  checkTrue( nrow(jaatha@likelihood.table) >= 9 )
  
  # Test reproducibility
  jaatha2 <- Jaatha.initialSearch(jaatha.tt, 20, 1)
  jaatha2 <- Jaatha.refinedSearch(jaatha2, 1, 10, sim.final=10, max.steps=10)
  checkEquals(jaatha@likelihood.table, jaatha2@likelihood.table)
}