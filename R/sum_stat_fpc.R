calcFpcBreaks <- function(dm, seg.sites, number=5, group=0) {
  trio_opts <- dm.getLociTrioOptions(dm, group, relative=TRUE)
  if (any(is.na(trio_opts))) trio_opts <- numeric(0)
    
  fpc.percent <- t(sapply(seg.sites, function(x) 
                          calcPercentFpcViolation(x, trio_opts) ))
  props <- seq(0, 1, length.out = number + 2)[-c(1, number+2)]
    
  if (group == 0) {
    # Options for the default group are directly in options
    dm@options[['fpc.breaks.near']] <- calcBreaks(fpc.percent[, 1], props)
    dm@options[['fpc.breaks.far']] <- calcBreaks(fpc.percent[, 2], props)
    if (length(trio_opts) > 0) 
      dm@options[['fpc.breaks.between']] <- calcBreaks(fpc.percent[, 3], props)
  } else {
    # Options for other groups go into a group vector and are copied to the 
    # options vector once the model for the individual groups are created.
    group.name <- paste("group", group, sep='.')
    if(is.null(dm@options[[group.name]])) dm@options[[group.name]] <- list()
    dm@options[[group.name]][['fpc.breaks.near']] <- calcBreaks(fpc.percent[, 1], props)
    dm@options[[group.name]][['fpc.breaks.far']] <- calcBreaks(fpc.percent[, 2], props)
    if (length(trio_opts) > 0)
      dm@options[[group.name]][['fpc.breaks.between']] <- calcBreaks(fpc.percent[, 3], props)
  }
  
  dm
}

calcBreaks <- function(values, props) {
  breaks <- unique(quantile(values, props, na.rm=TRUE))
  breaks[breaks == 0] <- ifelse(length(breaks)==1, 0.01, min(0.01, min(breaks[-1])/2))
  breaks
}

generateFpcStat <- function(seg_sites, dm, group = 0) {
  # Get breaks between classes for each dimension
  if (group == 0) {
    breaks_near <- dm@options[['fpc.breaks.near']]
    breaks_far <- dm@options[['fpc.breaks.far']]
    breaks_betw <- dm@options[['fpc.breaks.between']]
  } else {
    group.name <- paste("group", group, sep='.')
    breaks_near <- dm@options[[group.name]][['fpc.breaks.near']]
    breaks_far <- dm@options[[group.name]][['fpc.breaks.far']]
    breaks_betw <- dm@options[[group.name]][['fpc.breaks.between']]
  }
  if (is.null(breaks_near) || is.null(breaks_far)) 
    stop("Missing classes breaks for calculating fpc statistic.")
  
  # Calculate percent of SNPs that a private and fixed polym
  trio_opts <- dm.getLociTrioOptions(dm, group, relative=TRUE)
  stopifnot(is.null(breaks_betw) || all(!is.na(trio_opts)))
            
  if (any(is.na(trio_opts))) trio_opts <- numeric(0)
  percent <- sapply(seg_sites, function(x){
    calcPercentFpcViolation(x, trio_opts)
  })
  
  # Classify the loci accordingly
  locus_class <- matrix(1, nrow(percent), ncol(percent))
  for (brk in breaks_near) {
    locus_class[1,] <- locus_class[1,] + (percent[1,] > brk)
  }
  for (brk in breaks_far) {
    locus_class[2,] <- locus_class[2,] + (percent[2,] > brk)
  }
  if (!is.null(breaks_betw)) {
    for (brk in breaks_betw) {
      locus_class[3,] <- locus_class[3,] + (percent[3,] > brk)
    }
  }
  
  # Count the occurance of each class in a matrix
  if (is.null(breaks_betw)) {
    dims <- c(length(breaks_near) + 2, length(breaks_far) + 2)
  } else {
    dims <- c(length(breaks_near) + 2, 
              length(breaks_far) + 2, 
              length(breaks_betw) + 2)
  }
  countClasses(locus_class, dims)
}

countClasses <- function(classes, dimension) {
  stopifnot(nrow(classes) == length(dimension))
  stat <- array(0, dim = dimension)
  
  # Replace NA's with the last value
  for (r in 1:nrow(classes)) {
    classes[r, is.na(classes[r,])] <- dimension[r]
  }
  
  # Count occurences of classes
  if (nrow(classes) == 2) {
    for (r in 1:ncol(classes)) {
      stat[classes[1, r], classes[2, r]] <- 
        stat[classes[1, r], classes[2, r]] + 1
    }
  } else if (nrow(classes) == 3) {
    for (r in 1:ncol(classes)) {
      stat[classes[1, r], classes[2, r], classes[3, r]] <- 
        stat[classes[1, r], classes[2, r], classes[3, r]] + 1
    }
  } else stop("This function can only create 2 and 3-dim. arrays")
  
  stat
}