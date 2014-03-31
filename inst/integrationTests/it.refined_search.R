test.refinedSearch.normal <- function() {
  jaatha <- Jaatha.refinedSearch(jaatha.csi, 2, sim=10, sim.final=10,
                                 max.steps=25)
  checkTrue( is.matrix(jaatha@likelihood.table) )
  checkTrue( ncol(jaatha@likelihood.table) == 4 )
  checkTrue( nrow(jaatha@likelihood.table) >= 10 )

  # Test reproducibility
  jaatha2 <- Jaatha.refinedSearch(jaatha.csi, 2, sim=10, sim.final=10,
                                 max.steps=25)

  checkEquals(jaatha@likelihood.table, jaatha2@likelihood.table)
}
