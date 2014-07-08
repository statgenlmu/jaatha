calcFpcSumStat <- function(seg.sites, dm) {
  breaks.near <- dm@options[['fpc.breaks.near']]
  breaks.far  <- dm@options[['fpc.breaks.far']]

  sum.stat <- matrix(0, length(breaks.near)+2, length(breaks.far)+2)

  for(x in seg.sites) { 
    sum.stat <- addSegSitesToFpc(x, as.numeric(colnames(x)), breaks.near, 
                                 breaks.far, sum.stat)
  }

  sum.stat
}
  
calcFpcBreaks <- function(dm, seg.sites, number=5) {
  fpc.percent <- t(sapply(seg.sites, function(x) 
                          calcPercentFpcViolation(x, as.numeric(colnames(x))) ))
  props <- seq(0, 1, length.out = number + 2)[-c(1, number+2)]
  
  list(near=calcBreaks(fpc.percent[, 1], props), 
       far=calcBreaks(fpc.percent[, 2], props))
}

calcBreaks <- function(values, props) {
  breaks <- unique(quantile(values, props, na.rm=TRUE))
  breaks[breaks == 0] <- ifelse(length(breaks)==1, 0.01, min(0.01, min(breaks[-1])/2))
  breaks
}