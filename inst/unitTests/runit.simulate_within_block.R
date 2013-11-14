# --------------------------------------------------------------
# Authors:  Paul R. Staab
# Date:     2013-11-14
# Licence:  GPLv3 or later
# --------------------------------------------------------------

test.simulateWithinBlock <- function() {
  checkSumStat <- function(x, block) { 
    if (is.null(x$pars) || is.null(x$pars.normal)) return(FALSE)
    isInBlock(block, x$pars.normal)
  }

  sum.stats <- simulateWithinBlock(10, block.test, jaatha.csi)
  checkTrue( is.list(sum.stats) )
  checkEquals( length(sum.stats), 14 )
  checkTrue( all(sapply(sum.stats, checkSumStat, block=block.test) ) ) 

  sum.stats <- simulateWithinBlock(2, block.test, jaatha.tt)
  checkEquals( length(sum.stats), 6 )
  checkTrue( all(sapply(sum.stats, checkSumStat, block=block.test) ) ) 
}

test.runSimulatinos <- function() {
  set.seed(15)
  pars.test <- matrix(0.5, 3, 2)
  sum.stats1 <- runSimulations(pars.test, 1, jaatha.csi)
  checkEquals( length(sum.stats1), 3 )
  for (i in 1:3) {
    checkTrue( all(denormalize(pars.test[i, ], jaatha.csi) == sum.stats1[[i]]$pars) )
    checkTrue( all(pars.test[i, ] == sum.stats1[[i]]$pars.normal) )
    checkTrue( !is.null(sum.stats1[[i]]$poisson.vector) )
    checkTrue( sum(sum.stats1[[i]]$poisson.vector) > 0 )
  }

  if (require(multicore)) {
    set.seed(15)
    sum.stats2 <- runSimulations(pars.test, 2, jaatha.csi)
    checkEquals( length(sum.stats2), 3 )
    for (i in 1:3) {
      checkTrue( all(denormalize(pars.test[i, ], jaatha.csi) == sum.stats1[[i]]$pars) )
      checkTrue( all(pars.test[i, ] == sum.stats1[[i]]$pars.normal) )
      checkTrue( !is.null(sum.stats2[[i]]$poisson.vector) )
      checkTrue( sum(sum.stats2[[i]]$poisson.vector) > 0 )
      checkEquals( sum.stats1[[i]]$poisson.vector,
                   sum.stats2[[i]]$poisson.vector )
    }
  }
}
