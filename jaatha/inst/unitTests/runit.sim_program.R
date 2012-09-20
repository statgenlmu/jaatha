#!/usr/bin/Rscript --vanilla
#
# runit.sim_program.R
# Unit tests for the simProgram class
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-08-07
# Licence:  GPLv3 or later
#

# - Test setup -----------------------------------
 
library("RUnit")
library("jaatha")

# ------------------------------------------------

test.InstanceCreation <- function() {
  # Should both work; only one of simFunc/singleSimFunc given
  s1 <- new("SimProgram",
            name="name",
            executable="ms",
            possible.features="feature",
            possible.sum.stats="sum.stat",
            simFunc=sin) 

  s2 <- new("SimProgram",
            name="name",
            executable="",
            possible.features="feature",
            possible.sum.stats="sum.stat",
            singleSimFunc=sin) 

  #Print working?
  print(s1)
  print(s2)

  #Should not work because of wrong arguments
  checkException(
    new("SimProgram", name="name", executable="ms", possible.features="feature",
        possible.sum.stats="sum.stat",
        simFunc=sin,
        singleSimFunc=cos) 
  )

  checkException(
    new("SimProgram", executable="ms", possible.features="feature",
        simFunc=sin,
        possible.sum.stats="sum.stat") 
  )
  
  checkException(
    new("SimProgram", name="name", executable="ms", possible.features="feature",
        possible.sum.stats="sum.stat") 
  )
  
  checkException(
    new("SimProgram", name="name", executable="ms", possible.features="feature",
        simFunc=sin) 
  )
}

