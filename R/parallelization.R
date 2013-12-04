# ----------------------------------------------------------------------
# parallelization.R
# Functions for parallelizing Jaatha using 'foreach'
# 
# Author:   Paul R. Staab & Lisha Mathew 
# Email:    staab (at) bio.lmu.de
# Date:     2013-11-15
# Licence:  GPLv3 or later
# ----------------------------------------------------------------------

setParallelization <- function(cores=1) {
  checkType(cores, c("num", "single"))
  if (cores > 1 && .Platform$OS.type == "windows")
    stop("Parallelization is not supported on Windows.") 
}
