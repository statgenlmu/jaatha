# --------------------------------------------------------------
# Authors:  Paul R. Staab
# Date:     2013-11-14
# Licence:  GPLv3 or later
# --------------------------------------------------------------

test.init <- function() {
  sum.stat1 <- list(method="poisson.independent", value=1:10)
  sum.stat2 <- list(method="poisson.transformed", value=matrix(1, 3, 3),
                    transformation=diag)
  sum.stats <- list("ss1"=sum.stat1, "ss2"=sum.stat2)

  jaatha <- new("Jaatha", csi.sim.func, csi.par.ranges, csi.sum.stats)
  jaatha <- new("Jaatha", csi.sim.func, csi.par.ranges, sum.stats)

  jaatha <- new("Jaatha", csi.sim.func, csi.par.ranges, 
                sum.stats, 123, 1, FALSE)

  checkTrue( all(jaatha@par.ranges == csi.par.ranges) )
  checkEquals( length(jaatha@seeds), 3 )
  checkEquals( jaatha@seeds[1], 123 )
  checkEquals( jaatha@cores, 1 )
  checkEquals( jaatha@use.shm, FALSE )

  checkTrue( length(jaatha@sum.stats) == 2 )
  checkTrue( !is.null(jaatha@sum.stats$ss1) )
  checkTrue( !is.null(jaatha@sum.stats$ss2) )

  checkTrue( all(jaatha@sum.stats$ss1$value == 1:10) )
  checkTrue( !is.null(jaatha@sum.stats$ss1$transformation) )

  checkTrue( all(jaatha@sum.stats$ss2$value == 1) )
  checkTrue( !is.null(jaatha@sum.stats$ss2$transformation) )
  checkTrue( all(jaatha@sum.stats$ss2$value.transformed == c(1, 1, 1)) )
}

test.getParNames <- function() {
  checkEquals( getParNames(jaatha.tt), dm.getParameters(dm.tt) )
  checkEquals( getParNames(jaatha.mig), dm.getParameters(dm.mig) )
}

test.getParNumber <- function() {
  checkEquals( getParNumber(jaatha.tt), dm.getNPar(dm.tt) )
  checkEquals( getParNumber(jaatha.mig), dm.getNPar(dm.mig) )
}
