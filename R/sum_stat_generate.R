generateSumStats <- function(files, program, parameters, dm, seg_sites) {
  model_stats <- dm.getSummaryStatistics(dm)
  
  calc_seg_sites <- any(c('seg.sites', 'jsfs', 'pmc', 'fpc') %in% model_stats)
  if (missing(seg_sites) & calc_seg_sites) {
    if (program == 'ms') {
      seg_sites <- parseMsOutput(files, 
                                 dm.getSampleSize(dm), 
                                 dm.getLociNumber(dm))
    } else if (program == 'seqgen') {
      seg_sites <- parseSeqgenOutput(files, 
                                     sum(dm.getSampleSize(dm)),
                                     dm.getLociLength(dm),
                                     dm.getLociNumber(dm),
                                     dm.getLociTrioOptions(dm))
    } else {
      stop("Unknown program: ", program)
    }
  }

  # Add the parameters of the simulation
  sum_stats <- list(pars=parameters)

  # Add seg_sites
  if ('seg.sites' %in% model_stats) {
    sum_stats[['seg.sites']] <- seg_sites
  }
  
  # Add JSFS
  if ('jsfs' %in% model_stats) {
    sum_stats[['jsfs']] <- calcJsfs(seg_sites, dm.getSampleSize(dm))
  }
  
  # Add pmc statistic if needed
  if ('pmc' %in% model_stats) {
    sum_stats[['pmc']] <- createPolymClasses(seg_sites, dm)
  }
  
  # Add fpc statistic if needed
  if ('fpc' %in% model_stats) {
    sum_stats[['fpc']] <- generateFpcStat(seg_sites, dm)
  }
  
  # Add files if needed
  if ('file' %in% model_stats) {
    sum_stats[['file']] <- files
  } else {
    unlink(unlist(files))
  }
  
  sum_stats
}