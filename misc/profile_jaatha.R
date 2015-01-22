# Example script for profiling the initiual and refined search.

dm <- dm.createThetaTauModel(c(20,25), 200)
sim_data <- dm.simSumStats(dm.addSummaryStatistic(dm, 'seg.sites'), c(.1,5))
jaatha <- Jaatha.initialize(sim_data, dm, use_fpc=TRUE, smoothing=FALSE)

profile.is <- tempfile('jaatha_is2.prof')
gc()
Rprof(profile.is, memory.profiling=TRUE)
jaatha <- Jaatha.initialSearch(jaatha)
Rprof(NULL)
summaryRprof(profile.is, memory = 'both')

profile.rs <- "./jaatha_rs2.prof"
gc()
Rprof(profile.rs, memory.profiling=TRUE)
jaatha <- Jaatha.refinedSearch(jaatha, 2)
Rprof(NULL)
summaryRprof(profile.rs, memory = 'both')
