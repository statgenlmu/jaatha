# --------------------------------------------------------------
# Authors:  Paul R. Staab
# Date:     2013-11-13
# Licence:  GPLv3 or later
# --------------------------------------------------------------

test.normalize <- function() {
  dm <- dm.createThetaTauModel(11:12, 10)
  jaatha <- Jaatha.initialize(dm, jsfs=matrix(1,12,13)) 

  checkTrue( all( normalize(c(0.01, 1), jaatha) == c(0,0) ) ) 
  checkTrue( all( normalize(c(5, 20), jaatha) == c(1,1) ) ) 
}

test.denormalize <- function() {
  dm <- dm.createThetaTauModel(11:12, 10)
  jaatha <- Jaatha.initialize(dm, jsfs=matrix(1,12,13)) 

  checkTrue( sum(abs(denormalize(c(0, 0), jaatha)-c(0.01,1))) < 1e-11 ) 
  checkTrue( sum(abs(denormalize(c(1, 1), jaatha)-c(5,20))) < 1e-11 ) 
}

test.normalizeDenormalize <- function() {
  dm <- dm.createThetaTauModel(11:12, 10)
  jaatha <- Jaatha.initialize(dm, jsfs=matrix(1,12,13)) 
  
  for (x in 0:10/10) {
    pars <- rep(x, 2) 
    checkTrue( sum(abs(normalize(denormalize(x, jaatha), jaatha) - pars)) < 1e-11 ) 
  }
}
