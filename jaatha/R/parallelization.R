# Parallelization.R
# Functions for parallelizing Jaatha using the 'foreach' mechanism.
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-08-09
# Licence:  GPLv3 or later
#

#require(foreach)

setSingleCoreMode <- function() {
  registerDoSEQ()
}

setSimpleParallelization <- function() {
  require(doMC)
  registerDoMC()
}

setParrallelizationForInitialSearch <- function(jaatha) {
  if (jaatha@parallelization.model == "none") setSingleCoreMode()
  if (jaatha@parallelization.model == "simple") setSimpleParallelization()
  if (jaatha@parallelization.model == "nodes")
    stop('Model "nodes" is not yet implemented')
}
