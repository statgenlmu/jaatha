test.refinedSearch.normal <- function() {
  jaatha <- Jaatha.refinedSearch(jaatha.csi, 2, sim=10, sim.final=10,
                                 max.steps=25)

  # Test reproducibility
  jaatha2 <- Jaatha.refinedSearch(jaatha.csi, 2, sim=10, sim.final=10,
                                 max.steps=25)

  checkEquals(jaatha@likelihood.table, jaatha2@likelihood.table)
}
