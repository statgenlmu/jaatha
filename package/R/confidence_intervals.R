# --------------------------------------------------------------
# confidence_intervals.R 
# Contains a function for calculating bias corrected bootstrap 
# confidence intervals 
# 
# Authors:  Paul R. Staab
# Date:     2013-09-24
# Licence:  GPLv3 or later
# --------------------------------------------------------------

require(foreach)
require(doMC)

Jaatha.confidenceIntervals <- function(jaatha, conf.level=0.95, replicas, cores, log.folder=NULL) {
  # Get a seed for each replica plus for simulating data 
  set.seed(jaatha@seeds[1])
  seeds <- jaatha:::generateSeeds(replicas+3)[-(1:2)]

  if (!is.null(log.folder)) dir.create(log.folder, showWarnings=FALSE)

  # Simulate data under the fitted model
  set.seed(seeds[length(seeds)])
  est.pars <- Jaatha.getLikelihoods(jaatha, 1)[-(1:2)]
  sim.pars <- matrix(est.pars, replicas, jaatha@nPar, byrow=TRUE)
  data.sim <- dm.simSumStats(jaatha@dm, sim.pars, jaatha::Jaatha.defaultSumStats)

  registerDoMC(cores)
  bs.results <- foreach(i=1:replicas, .combine=rbind) %dopar% {
    cat("Starting run", i, "...\n")
    if (!is.null(log.folder)) sink(paste0(log.folder, "/run_", i, ".log"))
    jaatha.rep <- Jaatha.initialize(jaatha@dm, data.sim[i, ], seed=seeds[i]) 
    # Important: Use the same settings as in the original calls here:
    jaatha.rep <- Jaatha.initialSearch(jaatha.rep, 10, 2) 
    jaatha.rep <- Jaatha.refinedSearch(jaatha.rep, 2, 10)

    if (!is.null(log.folder)) save(jaatha.rep, file=paste0(log.folder, "/run_", i, ".Rda"))
    sink()
    return(Jaatha.getLikelihoods(jaatha.rep, 1)[-(1:2)])
  }

  conf.intervals <- list()
  for (i in 1:ncol(bs.results)) {
    par.name <- dm.getParameters(jaatha@dm)[i]
    conf.intervals[[par.name]] <- calcBCaConfInt(conf.level, bs.results[,i], est.pars[i], replicas)
  }
  
  return(conf.intervals)
}


calcBCaConfInt <- function(conf.level, boot, ML, replicas) {
  z.hat.null <- calcBiasCorrection(log(boot), ML=log(ML), replicas)
  a.hat <- calcAcceleration(log(boot))
  z.alpha <- qnorm(p=c((1-conf.level)/2, 1-(1-conf.level)/2)) 
  corr.quantiles <- pnorm(z.hat.null + (z.hat.null + z.alpha) / (1-a.hat*(z.hat.null + z.alpha)))  
  conf.int <- quantile(log(boot), probs=corr.quantiles) 
  return(exp(conf.int))
}


calcBiasCorrection <- function(parE, ML, replicas){
  return(qnorm(sum(parE < ML)/replicas))
}

calcAcceleration <- function(parE) {
  m <- mean(parE)
  nominator <- sum((m-parE)^3)
  denom <- 6*(sum((m-parE)^2))^(3/2)
  return(nominator/denom)
}

getCI <- function(a,b){
  zAlpha <- qnorm(0.975)
  pnorm(b+((zAlpha*c(-1,1)+b)/(1-a*(b+zAlpha*c(-1,1)))))
}

