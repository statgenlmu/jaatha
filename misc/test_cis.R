#!/usr/bin/Rscript --vanilla
#
# test_cis
# %DESCRIPTION%
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2013-10-22
# Licence:  GPLv3 or later
#

args <- commandArgs(TRUE)
set.seed(18)
cores <- 2

library(jaatha)

sampleFromModel <- function(x, y, z){
  return( c(rpois(10, x), rpois(10, y), rpois(10, z)) )
}

sim.func <- function(jaatha, sim.pars) {
  t(apply(sim.pars, 1, 
          function(x) sampleFromModel(x[1], x[2], x[3])))
}

in.interval <- rep(0, 3)

for (i in 1:250) {
  print(i)
  par.ranges <- matrix(c(1, 1, 1, 10, 10, 10), 3, 2)
  rownames(par.ranges) <- c('x', 'y', 'z')
  colnames(par.ranges) <- c('min', 'max')

  pars.true <- runif(3)*(par.ranges[,2] - par.ranges[,1])+par.ranges[,1]

  data.sim <- sampleFromModel(pars.true[1], pars.true[2], pars.true[3])

  jaatha <- new('Jaatha', sim.func, par.ranges, data.sim, 
                cores=cores, sim.package.size=10) 
  jaatha <- Jaatha.initialSearch(jaatha, 100, 2)
  jaatha <- Jaatha.refinedSearch(jaatha, 2, 100)
  jaatha <- Jaatha.confidenceIntervals(jaatha, 0.95, 96, cores)

  conf.ints <- jaatha@conf.ints
  in.interval <- in.interval + (conf.ints[,1] <= pars.true & pars.true <= conf.ints[,2]) 
}

print(in.interval/250)
