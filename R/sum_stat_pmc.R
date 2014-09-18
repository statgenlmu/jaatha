createPolymClasses <- function(seg_sites, dm) {
  # Get breaks between classes for each dimension
  breaks_priv <- dm@options[['pmc_breaks_private']]
  breaks_fixed <- dm@options[['pmc_breaks_fixed']]
  if (is.null(breaks_priv) || is.null(breaks_fixed)) 
    stop("Missing classes breaks for calculating Polym. statistic.")
  
  # Calculate percent of SNPs that a private and fixed polym
  percent <- sapply(seg_sites, classifyPolym, dm.getSampleSize(dm))
  
  # Classify the loci accordingly
  locus_class <- matrix(1, nrow(percent), ncol(percent))
  for (brk in breaks_priv) {
    locus_class[1,] <- locus_class[1,] + (percent[1,] > brk)
  }
  for (brk in breaks_fixed) {
    locus_class[2,] <- locus_class[2,] + (percent[2,] > brk)
  }
  
  # Count the occurance of each class in a matrix
  stat <- matrix(0, length(breaks_priv)+2, length(breaks_fixed)+2)
  for (r in 1:ncol(locus_class)) {
    if (any(is.na(locus_class[, r]))) {
      stat[length(breaks_priv)+2, length(breaks_fixed)+2] <-
        stat[length(breaks_priv)+2, length(breaks_fixed)+2] + 1  
    } else {
      stat[locus_class[1, r], locus_class[2, r]] <- 
        stat[locus_class[1, r], locus_class[2, r]] + 1
    }
  }
  stat
}