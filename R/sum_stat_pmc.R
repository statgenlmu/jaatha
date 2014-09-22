createPolymClasses <- function(seg_sites, dm, group = 0) {
  # Get breaks between classes for each dimension
  if (group == 0) {
    breaks_priv <- dm@options[['pmc_breaks_private']]
    breaks_fixed <- dm@options[['pmc_breaks_fixed']]
  } else {
    group.name <- paste("group", group, sep='.')
    breaks_priv <- dm@options[[group.name]][['pmc_breaks_private']]
    breaks_fixed <- dm@options[[group.name]][['pmc_breaks_fixed']]
  }
  if (is.null(breaks_priv) || is.null(breaks_fixed)) 
    stop("Missing classes breaks for calculating Polym. statistic.")
  
  # Calculate percent of SNPs that a private and fixed polym
  percent <- sapply(seg_sites, classifyPolym, dm.getSampleSize(dm, group))
  
  # Classify the loci accordingly
  locus_class <- matrix(1, nrow(percent), ncol(percent))
  for (brk in breaks_priv) {
    locus_class[1,] <- locus_class[1,] + (percent[1,] > brk)
  }
  for (brk in breaks_fixed) {
    locus_class[2,] <- locus_class[2,] + (percent[2,] > brk)
  }
  
  # Count the occurance of each class in a matrix
  countClasses(locus_class, c(length(breaks_priv)+2, length(breaks_fixed)+2))
}

calcPmcBreaks <- function(dm, seg_sites, number=2, group=0) {
  percent <- sapply(seg_sites, classifyPolym, dm.getSampleSize(dm, group))
  probs <- seq(0, 1, length.out = number + 2)[-c(1, number+2)]
  
  if (group == 0) {
    # Options for the default group are directly in options
    dm@options[['pmc_breaks_private']] <- calcBreaks(percent[1, ], probs)
    dm@options[['pmc_breaks_fixed']] <- calcBreaks(percent[2, ], probs)
  } else {
    # Options for other groups go into a group vector and are copied to the 
    # options vector once the model for the individual groups are created.
    group.name <- paste("group", group, sep='.')
    if(is.null(dm@options[[group.name]])) dm@options[[group.name]] <- list()
    dm@options[[group.name]][['pmc_breaks_private']] <- 
      calcBreaks(percent[1, ], probs)
    dm@options[[group.name]][['pmc_breaks_fixed']] <- 
      calcBreaks(percent[2, ], probs)
  }
  
  dm
}

transformPmc <- function(pmc) {
  c(as.vector(pmc[-nrow(pmc), -ncol(pmc)]), pmc[nrow(pmc), ncol(pmc)])
}