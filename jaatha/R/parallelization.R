# ----------------------------------------------------------------------
# Parallelization.R
# Functions for parallelizing Jaatha using the 'foreach' mechanism.
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-08-09
# Licence:  GPLv3 or later
# ----------------------------------------------------------------------

setSingleCoreMode <- function(jaatha) {
  registerDoSEQ()
}

setSimpleParallelization <- function(jaatha) {
  require(doMC)
  cores <- jaatha@cores
  if (cores == 0) cores <- NULL
  registerDoMC(cores)
}

setParrallelizationForInitialSearch <- function(jaatha) {
  if (jaatha@parallelization.model == "none")   setSingleCoreMode(jaatha)
  if (jaatha@parallelization.model == "simple") setSimpleParallelization(jaatha)
  if (jaatha@parallelization.model == "nodes")
    stop('Model "nodes" is not yet implemented')
}

setParrallelizationForRefineSearch <- function(jaatha) {
  if (jaatha@parallelization.model == "none")   setSingleCoreMode(jaatha)
  if (jaatha@parallelization.model == "simple") setSimpleParallelization(jaatha)
  if (jaatha@parallelization.model == "nodes")
    stop('Model "nodes" is not yet implemented')
}
