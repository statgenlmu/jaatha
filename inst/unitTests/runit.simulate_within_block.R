# --------------------------------------------------------------
# Authors:  Paul R. Staab
# Date:     2013-11-13
# Licence:  GPLv3 or later
# --------------------------------------------------------------

test.simulateWithinBlock <- function() {
  dm <- dm.createThetaTauModel(11:12, 10)
  jaatha <- Jaatha.initialize(dm, jsfs=matrix(1,12,13)) 

  block <- new("Block")
  block@border <- matrix(c(0, 0, 0.1, 0.1), 2, 2) 

  sum.stats <- simulateWithinBlock(3, block, jaatha)
  checkEquals( length(sum.stats), 7 )
}

test.runSimulatinos <- function() {
  dm <- dm.createThetaTauModel(11:12, 10)
  jaatha <- Jaatha.initialize(dm, jsfs=matrix(1,12,13)) 

  set.seed(15)
  pars <- matrix(0.5, 3, 2)
  sum.stats1 <- runSimulations(pars, 1, jaatha)
  checkEquals( length(sum.stats1), 3 )
  for (i in 1:3) {
    checkTrue( all(denormalize(pars[i, ], jaatha) == sum.stats1[[i]]$pars) )
    checkTrue( all(pars[i, ] == sum.stats1[[i]]$pars.normal) )
    checkTrue( !is.null(sum.stats1[[i]]$jsfs) )
    checkTrue( sum(sum.stats1[[i]]$jsfs) > 0 )
  }

  if (require(multicore)) {
    set.seed(15)
    sum.stats2 <- runSimulations(pars, 2, jaatha)
    checkEquals( length(sum.stats2), 3 )
    for (i in 1:3) {
      checkTrue( all(denormalize(pars[i, ], jaatha) == sum.stats1[[i]]$pars) )
      checkTrue( all(pars[i, ] == sum.stats1[[i]]$pars.normal) )
      checkTrue( !is.null(sum.stats2[[i]]$jsfs) )
      checkTrue( sum(sum.stats2[[i]]$jsfs) > 0 )
      checkEquals( sum.stats1[[i]]$jsfs, sum.stats2[[i]]$jsfs )
    }
  }
}
