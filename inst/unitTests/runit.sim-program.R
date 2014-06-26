#!/usr/bin/Rscript --vanilla
#
# Unit tests for Jaatha's sim_program mechanism
# 
# Author:   Paul R. Staab 
# Email:    staab (at) bio.lmu.de
# Date:     2012-08-07
# Licence:  GPLv3 or later
#

test.createSimProgram <- function() {
  # Is the list of sim_programs created
  checkTrue(exists(".jaatha"))
  checkTrue(exists("sim_progs", envir=.jaatha))
  sim_prog_nr = length(.jaatha$sim_progs)

  # Should not work because of wrong arguments
  checkException(createSimProgram("name", "feature", sin, cos))
  checkException(createSimProgram("name"))
  
  # Should work
  createSimProgram("test1", "feature", "sum.stat", sin)
  checkEquals(sim_prog_nr + 1, length(.jaatha$sim_progs))
  
  createSimProgram("test2", "feature", "sum.stat", sin, cos)
  checkEquals(sim_prog_nr + 2, length(.jaatha$sim_progs))
  
  createSimProgram("test2", "feature", "sum.stat", sin, cos, tan)
  checkEquals(sim_prog_nr + 2, length(.jaatha$sim_progs))
}

