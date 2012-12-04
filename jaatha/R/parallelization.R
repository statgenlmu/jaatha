# ----------------------------------------------------------------------
# Parallelization.R
# Functions for parallelizing Jaatha using the 'foreach' mechanism.
# 
# Author:   Paul R. Staab & Lisha Naduvilezhath 
# Email:    staab (at) bio.lmu.de
# Date:     2012-08-09
# Licence:  GPLv3 or later
# ----------------------------------------------------------------------

setSingleCoreMode <- function(jaatha) {
  registerDoSEQ()
}

setParallelization <- function(jaatha) {
  require(doMC)
  cores <- jaatha@cores
  if (cores == 0) cores <- NULL
  registerDoMC(cores)
}
