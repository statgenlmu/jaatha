# --------------------------------------------------------------
# Authors:  Paul R. Staab
# Date:     2013-11-13
# Licence:  GPLv3 or later
# --------------------------------------------------------------

test.normalize <- function() {
  checkTrue( all( normalize(c(0.01, 1), jaatha.tt) == c(0,0) ) ) 
  checkTrue( all( normalize(c(5, 20), jaatha.tt) == c(1,1) ) ) 
}

test.denormalize <- function() {
  checkTrue( sum(abs(denormalize(c(0, 0), jaatha.tt)-c(0.01,1))) < 1e-11 ) 
  checkTrue( sum(abs(denormalize(c(1, 1), jaatha.tt)-c(5,20))) < 1e-11 ) 
}

test.normalizeDenormalize <- function() {
  for (x in 0:10/10) {
    pars <- rep(x, 2) 
    checkTrue( sum(abs(normalize(denormalize(x, jaatha.tt), jaatha.tt) - pars)) < 1e-11 ) 
  }
}
