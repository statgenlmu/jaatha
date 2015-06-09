#' @importFrom coala get_outgroup_size
checkModelDataConsistency <- function(data, model) {
  if (!requireNamespace("PopGenome", quietly = TRUE)) {
    stop("Please install package 'PopGenome'")
  }
  
  if (get_outgroup_size(model) != length(data@outgroup)) {
    stop("Expecting an outgroup size of ", length(data@outgroup), 
         " from data, but is ", get_outgroup_size(model), " in model.")
  }
}


convPopGenomeToSegSites <- function(data, only_synonymous=FALSE) {
  if (!requireNamespace("PopGenome", quietly = TRUE)) {
    stop("Please install package 'PopGenome'")
  }
  
  seg_sites_list <- lapply(1:length(data@n.valid.sites), function(i) {
    if (data@n.valid.sites[[i]] == 0) return(NULL)
    seg_sites <- PopGenome::get.biallelic.matrix(data, i)
    
    # Sort individuals as Pop1, Pop2, Outgroup
    pop1 <- which(row.names(seg_sites) %in% data@populations[[1]])
    pop2 <- which(row.names(seg_sites) %in% data@populations[[2]])
    outgroup <- which(row.names(seg_sites) %in% data@outgroup)
    
    if (only_synonymous) {
      syn <- data@region.data@synonymous[[i]]
      syn[is.na(syn)] <- FALSE
      seg_sites <- seg_sites[c(pop1, pop2, outgroup), syn, drop=FALSE]
    } else {
      seg_sites <- seg_sites[c(pop1, pop2, outgroup), , drop=FALSE]
    }

    attr(seg_sites, "positions") <- 
      as.numeric(colnames(seg_sites)) / data@n.sites[[i]]
    seg_sites
  })
  
  list(seg_sites = seg_sites_list[!sapply(seg_sites_list, is.null)])
}


#' @importFrom coala coal_model feat_outgroup
createModelFromPopGenome <- function(data, quiet=FALSE) {
  if (!requireNamespace("PopGenome", quietly = TRUE)) {
    stop("Please install package 'PopGenome'")
  }
  
  stopifnot("GENOME" %in% is(data))
  sample_sizes <- sapply(data@populations, length)
  if (!quiet) { 
    message("Sample Sizes: ")
    for (pop in seq(along = sample_sizes)) {
      message(" Population ", pop, ": ", sample_sizes[pop])
    }
  }
  
  outgroup_size <- length(data@outgroup)
  outgroup_number <- length(sample_sizes) + 1
  if (!quiet) message("Outgroup size: ", paste(outgroup_size, collapse=" "), 
                      " (will be population ", outgroup_number, ")")


  loci_mask <- data@n.valid.sites > 0
  loci_length <- round(mean(data@n.valid.sites[loci_mask]))
  loci_number <- sum(loci_mask)
  
  if (!quiet) message("Number of Loci: ", loci_number)
  if (!quiet) message("Average Loci Length: ", loci_length, "bp")
  
  # Calculate TS/TV, but don"t add it to the model
  tstv_ratio <- sum(sapply(data@region.data@transitions[loci_mask], sum)) /
      sum(sapply(data@region.data@transitions[loci_mask], function(x) sum(1-x)))
  if (!quiet) message("Observed TS/TV: ", tstv_ratio, " (Not added to Model)")
  
  coal_model(c(sample_sizes, outgroup_size), loci_number, loci_length) +
    feat_outgroup(outgroup_number)
}
