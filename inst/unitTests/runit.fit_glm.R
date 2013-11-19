# --------------------------------------------------------------
# Author:   Paul R. Staab
# Date:     2013-11-14
# Licence:  GPLv3 or later
# --------------------------------------------------------------

test.fitGlm <- function() {
  # poisson.independent
  glms.fitted.csi <- fitGlm(sim.data.csi, jaatha.csi) 
  checkTrue( is.list(glms.fitted.csi) )
  checkEquals( length(glms.fitted.csi), 1 )

  checkEquals( length(glms.fitted.csi$poisson.vector), 
              length(sim.data.csi[[1]]$poisson.vector) )
  checkTrue( all(sapply(glms.fitted.csi$poisson.vector, is)[1, ] == "glm") )

  # poisson.transformed
  glms.fitted.tt <- fitGlm(sim.data.tt, jaatha.tt) 
  checkTrue( is.list(glms.fitted.tt) )
  checkEquals( length(glms.fitted.tt), 1 )

  checkTrue( is.list(glms.fitted.tt$jsfs) )
  checkEquals( length(glms.fitted.tt$jsfs), 
              length(jaatha.tt@sum.stats$jsfs$value.transformed) )
  checkTrue( all(sapply(glms.fitted.tt$jsfs, is)[1, ] == "glm") )
}
