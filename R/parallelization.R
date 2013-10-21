# ----------------------------------------------------------------------
# parallelization.R
# Functions for parallelizing Jaatha using 'foreach'
# 
# Author:   Paul R. Staab & Lisha Mathew 
# Email:    staab (at) bio.lmu.de
# Date:     2013-10-21
# Licence:  GPLv3 or later
# ----------------------------------------------------------------------

setParallelization <- function(cores=1) {
  checkType(cores, c("num", "single"))
  if (cores > 1) {
    if (!require(doMC)) 
      stop("You need the package 'doMC' to run Jaatha on multiple cores.
            This package is not available on Windows systems.") 
    doMC::registerDoMC(cores)
  } else {
    foreach::registerDoSEQ()
  }
}
