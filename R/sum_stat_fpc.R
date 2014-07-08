calcFpcSumStat <- function(seg.sites, dm, group=0) {
  if (group == 0) {
    breaks.near <- dm@options[['fpc.breaks.near']]
    breaks.far  <- dm@options[['fpc.breaks.far']]
  } else {
    group.name <- paste("group", group, sep='.')
    breaks.near <- dm@options[[group.name]][['fpc.breaks.near']]
    breaks.far  <- dm@options[[group.name]][['fpc.breaks.far']]
  }
  stopifnot( !is.null(breaks.near) )
  stopifnot( !is.null(breaks.far) )

  sum.stat <- matrix(0, length(breaks.near)+2, length(breaks.far)+2)

  for(x in seg.sites) { 
    sum.stat <- addSegSitesToFpc(x, as.numeric(colnames(x)), breaks.near, 
                                 breaks.far, sum.stat)
  }

  sum.stat
}
  
calcFpcBreaks <- function(dm, seg.sites, number=5, group=0) {
  fpc.percent <- t(sapply(seg.sites, function(x) 
                          calcPercentFpcViolation(x, as.numeric(colnames(x))) ))
  props <- seq(0, 1, length.out = number + 2)[-c(1, number+2)]
    
  if (group == 0) {
    # Options for the default group are directly in options
    dm@options[['fpc.breaks.near']] <- calcBreaks(fpc.percent[, 1], props)
    dm@options[['fpc.breaks.far']] <- calcBreaks(fpc.percent[, 2], props)
  } else {
    # Options for other groups go into a group vector and are copied to the 
    # options vector once the model for the individual groups are created.
    group.name <- paste("group", group, sep='.')
    if(is.null(dm@options[[group.name]])) dm@options[[group.name]] <- list()
    dm@options[[group.name]][['fpc.breaks.near']] <- calcBreaks(fpc.percent[, 1], props)
    dm@options[[group.name]][['fpc.breaks.far']] <- calcBreaks(fpc.percent[, 2], props)
  }
  
  dm
}

calcBreaks <- function(values, props) {
  breaks <- unique(quantile(values, props, na.rm=TRUE))
  breaks[breaks == 0] <- ifelse(length(breaks)==1, 0.01, min(0.01, min(breaks[-1])/2))
  breaks
}