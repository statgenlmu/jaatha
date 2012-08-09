# Parallelization.R
# Functions for parallelizing Jaatha using the automagic 'foreach' mechanism.
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-08-09
# Licence:  GPLv3 or later
#

require(foreach)

setSingleCoreMode <- function() {
  registerDoSEQ()
}

setSimpleParallelization <- function() {
  require(doMC)
  registerDoMC()
}
