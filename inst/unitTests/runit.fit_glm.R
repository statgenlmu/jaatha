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

  # poisson.smoothed
  glms.fitted.smooth <- fitGlm(smooth.sim.data, smooth.jaatha)
  checkTrue( is.list(glms.fitted.smooth) )
  checkEquals( length(glms.fitted.smooth), 1 )

  checkTrue( is.list(glms.fitted.smooth$mat) )
  checkEquals( length(glms.fitted.smooth$mat), 1 )
  checkTrue( is(glms.fitted.smooth$mat[[1]]) == "glm")
}

test.convertSimResultsToDataFrame <- function() {
  smooth.df <- convertSimResultsToDataFrame(smooth.sim.data, "mat")
  checkTrue( is.data.frame(smooth.df) )
  checkEquals( ncol(smooth.df), 5 )
  checkEquals( nrow(smooth.df),
              length(as.vector(smooth.sim.data[[1]]$mat))*length(smooth.sim.data))
  checkTrue( !is.null(colnames(smooth.df)) )
  checkTrue( all(colnames(smooth.df) == c("x", "y", "i", "j", "sum.stat")) )
}
