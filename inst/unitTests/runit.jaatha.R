# --------------------------------------------------------------
# Authors:  Paul R. Staab
# Date:     2013-11-14
# Licence:  GPLv3 or later
# --------------------------------------------------------------

test.init <- function() {
  sum.stat1 <- list(method="poisson.independent", value=1:10)
  sum.stat2 <- list(method="poisson.transformed", value=matrix(1, 3, 3),
                    transformation=diag)
  sum.stat3 <- list(method="poisson.smoothing", model="a+b", value=matrix(1, 3, 3))
  sum.stats <- list("ss1"=sum.stat1, "ss2"=sum.stat2, "ss3"=sum.stat3)

  jaatha <- new("Jaatha", csi.sim.func, csi.par.ranges, csi.sum.stats)
  jaatha <- new("Jaatha", csi.sim.func, csi.par.ranges, sum.stats)

  jaatha <- new("Jaatha", csi.sim.func, csi.par.ranges, 
                sum.stats, 123, 1, FALSE)

  checkTrue( all(jaatha@par.ranges == csi.par.ranges) )
  checkEquals( length(jaatha@seeds), 3 )
  checkEquals( jaatha@seeds[1], 123 )
  checkEquals( jaatha@cores, 1 )
  checkEquals( jaatha@use.shm, FALSE )

  checkTrue( length(jaatha@sum.stats) == 3 )
  checkTrue( !is.null(jaatha@sum.stats$ss1) )
  checkTrue( !is.null(jaatha@sum.stats$ss2) )
  checkTrue( !is.null(jaatha@sum.stats$ss3) )

  checkTrue( all(jaatha@sum.stats$ss1$value == 1:10) )
  checkTrue( !is.null(jaatha@sum.stats$ss1$transformation) )

  checkTrue( all(jaatha@sum.stats$ss2$value == 1) )
  checkTrue( !is.null(jaatha@sum.stats$ss2$transformation) )
  checkTrue( all(jaatha@sum.stats$ss2$value.transformed == c(1, 1, 1)) )

  checkTrue( all(jaatha@sum.stats$ss3$value == 1) )
}

test.getParNames <- function() {
  checkEquals( getParNames(jaatha.tt), dm.getParameters(dm.tt) )
  checkEquals( getParNames(jaatha.mig), dm.getParameters(dm.mig) )
}

test.getParNumber <- function() {
  checkEquals( getParNumber(jaatha.tt), dm.getNPar(dm.tt) )
  checkEquals( getParNumber(jaatha.mig), dm.getNPar(dm.mig) )
}

test.JaathaInitialize <- function() {
  jsfs <- sum.stats.tt$jsfs
  checkType(Jaatha.initialize(dm.tt, jsfs), "jaatha") 
  checkType(Jaatha.initialize(demographic.model=dm.tt, jsfs=jsfs), "jaatha") 
  checkType(Jaatha.initialize(dm.tt, sum.stats.tt), "jaatha") 
  checkException(Jaatha.initialize(dm.tt))
  checkException(Jaatha.initialize(NULL))
  checkException(Jaatha.initialize(jsfs=NULL))

  # Seed
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 123)
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, seed=123)
  checkType(jaatha, "jaatha")
  checkTrue(jaatha@seeds[1] == 123)

  # Cores
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 123, 1)
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 123, cores=1)
  checkType(jaatha, "jaatha")
  checkTrue(jaatha@cores == 1)

  # Scaling Factor
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 123, 1, scaling.factor=17)
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 123, 1, 17)
  checkType(jaatha, "jaatha")
  checkTrue(jaatha@opts[['scaling.factor']] == 17)

  # Use Shm
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 123, 1, 17, TRUE)
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 123, 1, 17, use.shm=TRUE)
  checkType(jaatha, "jaatha")
  checkTrue(jaatha@use.shm)

  # Folded
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 123, 1, 17, TRUE, TRUE)
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 123, 1, 17, TRUE, folded=TRUE)
  checkType(jaatha, "jaatha")

  # Smoothing  
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 123, 1, 17, TRUE, FALSE, TRUE)
  jaatha <- Jaatha.initialize(dm.tt, sum.stats.tt, 123, 1, 17, TRUE, FALSE,
                              smoothing=TRUE)
  checkType(jaatha, "jaatha")
  checkTrue(jaatha@sum.stats$jsfs$method == "poisson.smoothing")
  checkTrue(!is.null(jaatha@sum.stats$jsfs$model))
  checkType(jaatha@sum.stats$jsfs$model, "char")

  checkException(Jaatha.initialize(dm.tt, sum.stats.tt, folded=TRUE,
                                   smoothing=TRUE))
}
