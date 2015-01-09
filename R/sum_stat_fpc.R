calcFpcBreaks <- function(dm, seg.sites, population, props=c(.2, .5), group=0) {
  if (any(props > 1)) stop('props greater then one')
    
  fpc.percent <- calcPercentFpcViolation(seg.sites, 
                                         getIndOfPop(dm, population),
                                         dm.getLociLengthMatrix(dm, group))
    
  key <- paste0('fpc_breaks_pop', population)
  breaks <- list(near=calcBreaks(fpc.percent[, 1], props),
                 far=calcBreaks(fpc.percent[, 2], props),
                 mut=calcBreaks(fpc.percent[, 6], props))
  
  if (group == 0) {
    # Options for the default group are directly in options
    dm@options[[key]] <- breaks
  } else {
    # Options for other groups go into a group vector and are copied to the 
    # options vector once the model for the individual groups are created.
    group.name <- paste("group", group, sep='.')
    if(is.null(dm@options[[group.name]])) dm@options[[group.name]] <- list()
    dm@options[[group.name]][[key]] <- breaks
  }
  
  dm
}

calcBreaks <- function(values, props) {
  breaks <- unique(quantile(values, props, na.rm=TRUE))
  breaks[breaks == 0] <- ifelse(length(breaks)==1, 0.01, min(0.01, min(breaks[-1])/2))
  breaks
}

generateFpcStat <- function(seg_sites, dm, population, group = 0) {
  # Get breaks between classes for each dimension
  
  key <- paste0('fpc_breaks_pop', population)
  if (group == 0) {
    breaks_near <- dm@options[[key]]$near
    breaks_far <- dm@options[[key]]$far
    breaks_betw <- dm@options[[key]]$mut
  } else {
    group.name <- paste("group", group, sep='.')
    breaks_near <- dm@options[[group.name]][[key]]$near
    breaks_far <- dm@options[[group.name]][[key]]$far
    breaks_betw <- dm@options[[group.name]][[key]]$mut
  }
  if (is.null(breaks_near) || is.null(breaks_far)) 
    stop("Missing classes breaks for calculating fpc statistic.")
  
  percent <- calcPercentFpcViolation(seg_sites, 
                                     getIndOfPop(dm, population),
                                     dm.getLociLengthMatrix(dm, group))
  
  # Classify the loci accordingly
  locus_class <- matrix(1, nrow(percent), 3)
  for (brk in breaks_near) {
    locus_class[,1] <- locus_class[,1] + (percent[,1] > brk)
  }
  for (brk in breaks_far) {
    locus_class[,2] <- locus_class[,2] + (percent[,2] > brk)
  }
  if (!is.null(breaks_betw)) {
    for (brk in breaks_betw) {
      locus_class[,3] <- locus_class[,3] + (percent[,6] > brk)
    }
  }
  
  # Count the occurance of each class in a matrix
  dims <- c(length(breaks_near) + 1, 
            length(breaks_far) + 1, 
            length(breaks_betw) + 1)
  countClasses(t(locus_class), dims)
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