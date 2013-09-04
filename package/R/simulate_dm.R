# --------------------------------------------------------------
# simulate_dm.R
# Default simulation function for demographic models 
# 
# Authors:  Paul R. Staab
# Date:     2013-09-04
# Licence:  GPLv3 or later
# --------------------------------------------------------------

simulateDemographicModel <- function(jaatha, sim.pars) {
  sumStats <- dm.simSumStats(jaatha@opts[['dm']], sim.pars)
  return(t(sapply(sumStats, summarizeJSFS)))
}

simulateDemographicModelFolded <- function(jaatha, sim.pars) {
  sumStats <- dm.simSumStats(jaatha@opts[['dm']], sim.pars)
  return(t(sapply(sumStats, summarizeFoldedJSFS)))
}


  # Scale SumStats if we use scaling
  # sumStats <- sumStats * jaatha@scaling.factor
