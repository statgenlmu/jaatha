# ----------------------------------------------------------------------
# Parallelization.R
# Functions for parallelizing Jaatha using the 'foreach' mechanism.
# 
# Author:   Paul R. Staab & Lisha Naduvilezhath 
# Email:    staab (at) bio.lmu.de
# Date:     2012-08-09
# Licence:  GPLv3 or later
# ----------------------------------------------------------------------

setParallelization <- function(jaatha) {
  cores <- jaatha@cores
  if (cores > 1) {
    require(doMC)
    registerDoMC(cores)
  } else {
    registerDoSEQ()
  }
}
