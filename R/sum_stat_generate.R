generateSumStats <- function(file, program, parameters, dm) {
  model_stats <- dm.getSummaryStatistics(dm)
  
  # Detemine if we need to generate seg.sites
  generate_seg_sites <- 'seg.sites' %in% model_stats
  if (any(c('fpc', 'pmc') %in% model_stats)) {
    generate_seg_sites <- TRUE
  }
  
  # Parse the simulation output
  sum_stats <- parseOutput(file, dm.getSampleSize(dm), dm.getLociNumber(dm), 
                           program, 'jsfs' %in% model_stats, generate_seg_sites,
                           dm.getLociTrioOptions(dm))
  
  # Add the parameters of the simulation
  sum_stats[['pars']] <- parameters
  
  # Add pmc statistic if needed
  if ('pmc' %in% model_stats) {
    stopifnot(!is.null(sum_stats$seg.sites))
    sum_stats[['pmc']] <- createPolymClasses(sum_stats$seg.sites, dm)
  }
  
  # Add fpc statistic if needed
  if ('fpc' %in% model_stats) {
    stopifnot(!is.null(sum_stats$seg.sites))
    sum_stats[['fpc']] <- generateFpcStat(sum_stats$seg.sites, dm)
  }
  
  # Remove seg.sites if it was just there to generate pmc or fpc.
  if (!'seg.sites' %in% model_stats) sum_stats[['seg.sites']] <- NULL
  sum_stats  
  
  # Add the simulation file if needed, or delete it otherwise
  if ("file" %in% model_stats) {
    sum_stats[['file']] <- file
  } else {
    unlink(file)
  }
  
  sum_stats
}