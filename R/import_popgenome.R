convPopGenomeToSegSites <- function(data, only_synonymous=FALSE) {
  seg_sites_list <- lapply(1:length(data@n.valid.sites), function(i) {
    if (data@n.valid.sites[[i]] == 0) return(NULL)
    seg_sites <- data@region.data@biallelic.matrix[[i]]
    
    # Sort individuals as Pop1, Pop2, Outgroup
    pop1 <- which(row.names(seg_sites) %in% data@populations[[1]])
    pop2 <- which(row.names(seg_sites) %in% data@populations[[2]])
    outgroup <- which(row.names(seg_sites) %in% data@outgroup)
    
    if (only_synonymous) {
      syn <- data@region.data@synonymous[[i]]
      syn[is.na(syn)] <- FALSE
      seg_sites <- seg_sites[c(pop1, pop2, outgroup), syn]
    } else {
      seg_sites <- seg_sites[c(pop1, pop2, outgroup), ]
    }

    attr(seg_sites, "positions") <- 
      as.numeric(colnames(seg_sites)) / data@n.sites[[i]]
    seg_sites
  })
  
  list(seg.sites = seg_sites_list[!sapply(seg_sites, is.null)])
}

dm.createModelFromPopGenome <- function(data, finite_sites = TRUE, 
                                        base_frequencies = c(A=.25, C=.25, G=.25, T=.25)) {
  
  stopifnot("GENOME" %in% is(data))
  sample_sizes <- sapply(data@populations, length)
  
  loci_mask <- data@n.valid.sites > 0
  loci_length <- data@n.valid.sites[loci_mask]
  
  dm <- dm.createDemographicModel(sample_sizes, 
                                  length(loci_length), 
                                  round(mean(loci_length)))
  
  if (finite_sites) {
    tstv_ration <- sum(sapply(data@region.data@transitions[loci_mask], sum)) /
      sum(sapply(data@region.data@transitions[loci_mask], function(x) sum(1-x)))
    
    dm <- dm.setMutationModel(dm, mutation.model = "HKY", 
                              base.frequencies = base_frequencies,
                              tstv.ratio = tstv_ration)
  }
  
  dm
}