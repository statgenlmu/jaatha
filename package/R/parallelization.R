# ----------------------------------------------------------------------
# parallelization.R
# Functions for parallelizing Jaatha using 'foreach'
# 
# Author:   Paul R. Staab & Lisha Mathew 
# Email:    staab (at) bio.lmu.de
# Date:     2013-09-04
# Licence:  GPLv3 or later
# ----------------------------------------------------------------------

setParallelization <- function(jaatha) {
  cores <- jaatha@cores
  if (cores > 1) {
    require(doMC)
    .registerDoMC(cores)
  } else {
    registerDoSEQ()
  }
}

# Ugly workaround to prevent note that registerDoMC is not visible at the moment
# setParallelization is loaded
.registerDoMC <- function(cores) {
    registerDoMC(cores)
}
